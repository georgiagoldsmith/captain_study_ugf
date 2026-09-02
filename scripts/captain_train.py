#!/usr/bin/env python
"""CAPTAIN 3 port of the UGF 'prioritization_9 / seed 8052026 original' run.

Reproduces config_train_empirical_K_growth_seed8052026.txt +
config_predict_prioritization_9_seed8052026.txt as closely as the v3 preview
allows. Points at the existing captain2-main/ugf_data tree -- no data is copied.

Run with the captain3preview environment:
    uv run python captain_train.py

Grid facts (computed from the rasters, for reference):
    raster                235 x 578 = 135,830 cells
    mask (AOH union)      44,004 valid cells      <- v2's n_pus planning graph
    cost extent           46,134 cells            <- v2's total_cells_T
    existing PAs in mask   8,374 cells
    disturbance >= 1.0    25,938 cells (hard-excluded)
    free & selectable     13,337 cells
    area target           floor(0.3 * 46,134) = 13,840 total -> 5,466 new
"""

from __future__ import annotations

import gc
import glob
import os
import resource
from pathlib import Path

import numpy as np
import torch

import captain as cn

# =============================================================================
# Configuration -- mirrors the v2 config files
# =============================================================================

V2_DIR = Path(os.environ.get(
    "V2_DIR",
    "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data",
))
if not V2_DIR.exists():
    raise FileNotFoundError(f"v2 data dir not found -- update V2_DIR: {V2_DIR}")
SDM_DIR = V2_DIR / "present_habitat_suitability"
ENV_DIR = V2_DIR / "environmental_layers"
TRAIT_FILE = V2_DIR / "species_traits.csv"

RESULTS_DIR = Path(
    os.environ.get(
        "RESULTS_DIR",
        Path(__file__).resolve().parent.parent / "outputs" / "captain3_seed8052026",
    )
)
MASK_FILE = RESULTS_DIR / "aoh_union_mask.npy"

# Disturbance drives carrying capacity. Two layers are available and they are
# NOT interchangeable in every role:
#   disturbance.tif      categorical LULC-derived, 4 unique values, 26,569 cells
#                        at exactly 1.0 (urban/cropland/bare/water)
#   disturbance_ghm.tif  continuous Global Human Modification, 883 unique
#                        values, range 0.001-0.974, NOTHING at exactly 1.0
DISTURBANCE_FILE = os.environ.get("DISTURBANCE_FILE", "disturbance.tif")
# The hard-exclusion mask is a LAND-ELIGIBILITY rule (you cannot protect a city),
# not a statement about disturbance intensity -- so it always comes from the
# categorical layer. Deriving it from gHM instead would block zero cells and
# balloon the eligible pool from 13,337 to ~37,600, which would confound any
# gHM-vs-LULC comparison with a completely different feasible set.
EXCLUSION_FILE = "disturbance.tif"
# Cost layer. costs.tif is NOT cocoa suitability: prioritzr.R line 161 normalises
# the real layer into `cocoa_suitability_ugf_norm`, then line 208 exports
# `cocoa_ugf_norm` instead -- a variable never assigned in that script, left over
# in the session from the coarse cocoa_ugf.tif. Confirmed: costs.tif correlates
# 1.000 with cocoa_ugf.tif and -0.015 with the real cocoa suitability.
# costs_cocoa_suitability.tif is line 161's raster, resampled to the 3km grid.
COST_FILE = os.environ.get("COST_FILE", "costs.tif")
# Optional future disturbance -> SpatialData gets a delta_per_step and carrying
# capacity declines over the episode, so populations can fall and loss_impact can
# fire. Without it K is static, nothing declines, and the temporal dimension is
# inert (verified: 0/18 species lose any population).
# MUST pair like with like: disturbance_future_ssp4_rcp34_2075.tif is CATEGORICAL
# (3 values), matching disturbance.tif (4 values), NOT gHM (7,330 values).
# Differencing gHM against it would produce a meaningless per-step delta.
FUTURE_DISTURBANCE_FILE = os.environ.get("FUTURE_DISTURBANCE_FILE", "")

# [general] seed. Controls policy weight init, the top-k tie-break RNG, and the
# ES noise stream -- i.e. exactly the sources of run-to-run variation the v2
# seed sweep was probing.
SEED = int(os.environ.get("SEED", 8052026))
N_TIME_STEPS = 30               # [general] steps
N_EPOCHS = int(os.environ.get("N_EPOCHS", 100))   # [general] epochs
# v2's batch_size = 3 does NOT carry over: v3 estimates both the NES gradient
# and the Jaccard index used to adapt sigma from these samples, and 3 is far
# too few for either (it drove sigma 0.2 -> 1.8 over 100 epochs).
N_PERTURBATIONS = int(os.environ.get("N_PERTURBATIONS", 30))
# Number of EpisodeRunners = parallel worker processes. Affects wall-clock
# only, never the result. ~830 MB each (18 per-species dispersal matrices).
N_RUNNERS = int(os.environ.get("N_RUNNERS", 8))
TOTAL_TARGET_CELLS = int(os.environ.get("TOTAL_TARGET_CELLS", 13_840))
# Protection schedule. Default = one-shot: the whole budget is spent at t=0 and
# the remaining 29 steps only simulate the consequences. Setting CELLS_PER_STEP
# smaller staggers protection across the episode, so the policy re-decides each
# step against the disturbance level *at that step*. That only means anything
# when FUTURE_DISTURBANCE_FILE is set -- with static disturbance nothing changes
# between steps, so staging just re-solves the same problem 30 times.
# GlobalBudgetManager caps each step at (total_target - already_protected), and
# the 8,374 locked-in PAs already count as protected, so the per-step figure
# must divide the NEW cells: (13,201 - 8,374) / 30 = 161.
CELLS_PER_STEP = int(os.environ.get("CELLS_PER_STEP", TOTAL_TARGET_CELLS))
MIN_SUITABILITY = 0.5           # [env_settings] prob_threshold
MIN_COST = 0.1                  # [env_settings] min_cost
DISPERSAL_WINDOW = int(os.environ.get("DISPERSAL_WINDOW", 3))   # in CELLS, not km
# Per-species dispersal builds one sparse matrix per species. At 1km with a 9km
# range that is 154M nonzeros x 18 species = 22 GB per runner, which does not fit
# in 24 GiB. SCALAR_DISPERSAL=1 uses a single shared matrix (1.2 GB) at the cost
# of collapsing the AVONET Hand-Wing spread (lambda 0.30-2.23) to a uniform 1.0.
SCALAR_DISPERSAL = os.environ.get("SCALAR_DISPERSAL", "0") == "1"
MEAN_DISPERSAL_RATE = 1.0       # [species_settings] mean_dispersal_rate
REWARD_WEIGHTS = np.array([0.7, 0.3])   # species_risk, cost (carbon weight was 0)
# [policy] nn_nodes = 3 2. That default is very tight: 13 features squeezed through
# 3 then 2 ReLU units, 53 parameters. Measured consequence -- in a degenerate run
# every hidden unit dies for 78% of cells, so 10,437 distinct cells collapse onto
# ONE score and top-k breaks the tie by array order, painting raster-row stripes.
# Widen (e.g. "16,8") or switch activation (tanh never zeroes) to give the network
# room to discriminate.
HIDDEN_DIM = [int(x) for x in os.environ.get("HIDDEN_DIM", "3,2").split(",")]
ACTIVATION = os.environ.get("ACTIVATION", "relu")
# v3 multiplies reward_weights by a SECOND, auto-fitted calibration vector.
# Leave False to keep 0.7/0.3 meaning what it meant in v2 (the two components
# are already on comparable scales here, so calibration buys nothing and
# silently re-weights cost ~4x above extinction risk). Set True only if you
# intend to re-tune REWARD_WEIGHTS against variance-normalised units.
CALIBRATE_REWARDS = False
# [extinction_risk] risk_weights, reversed: v2 lists CR..LC, v3 indexes LC..CR
THREAT_WEIGHTS = np.array([-1, -8, -16, -32, -64])
# [extinction_risk] relative_protected_range_thresholds, padded to n_classes + 1
PROTECT_THRESHOLDS = np.array([0.0, 0.1, 0.25, 0.50, 0.80, 1.0])

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

os.makedirs(RESULTS_DIR, exist_ok=True)


# =============================================================================
# Mask -- v3 needs ONE shared mask so every layer flattens to the same cells
# =============================================================================


def build_mask() -> np.ndarray:
    """Union of all species' AOH == v2's reference_grid_pu (44,004 cells).

    Must be float: grid_utils.flatten_grid() assigns np.nan into it in place.
    """
    if MASK_FILE.exists():
        return np.load(MASK_FILE)

    acc = None
    for f in sorted(glob.glob(str(SDM_DIR / "*.tif"))):
        arr, _ = cn.data_loader.load_map(f, nan_to_num=True)
        acc = arr if acc is None else acc + arr

    mask = (acc > 0).astype(np.float32)
    np.save(MASK_FILE, mask)
    print(f"mask: {int(mask.sum())} valid cells")
    return mask


# =============================================================================
# Episode setup
# =============================================================================


def create_episode_runner() -> cn.EpisodeRunner:
    mask = build_mask()

    # Habitat suitability: binary AOH (1 inside range, NaN outside).
    # future_dir omitted -- use_future_sdms = False in the v2 config.
    sdm = cn.load_spatial_data_from_dir(
        dir=SDM_DIR,
        mask=mask,
        lower_bound=0,
        upper_bound=1,
        min_threshold=MIN_SUITABILITY,
    )

    disturbance = cn.load_spatial_data(
        file=ENV_DIR / DISTURBANCE_FILE,
        future_file=(ENV_DIR / FUTURE_DISTURBANCE_FILE) if FUTURE_DISTURBANCE_FILE else None,
        n_time_steps=N_TIME_STEPS if FUTURE_DISTURBANCE_FILE else None,
        mask=mask,
        lower_bound=0,
        upper_bound=1,
    )

    # min_cost 0.1 applied as the clip floor
    costs = cn.load_spatial_data(
        file=ENV_DIR / COST_FILE, mask=mask, lower_bound=MIN_COST, upper_bound=1
    )

    # add_to_existing_protected_areas = True -> seed the matrix with real PAs
    # rather than zeros. These cells are then excluded from further action
    # automatically via env.no_action_mask, and count toward TOTAL_TARGET_CELLS.
    protection = cn.load_spatial_data(
        file=ENV_DIR / "protected_areas.tif", mask=mask, lower_bound=0, upper_bound=1
    )

    # Hard-excluded cells (urban/cropland/bare/water, disturbance == 1.0).
    # In v2 you made these unaffordable via graph_cost = budget + 1; v3 has a
    # first-class action_mask, so the cost layer stays undistorted.
    dist_grid, _ = cn.data_loader.load_map(ENV_DIR / EXCLUSION_FILE, nan_to_num=True)
    action_mask = cn.SpatialData(
        data=(dist_grid >= 1.0).astype(np.float32), mask=mask, lower_bound=0, upper_bound=1
    )

    # -------------------------------------------------------------------------
    # Species parameters
    # -------------------------------------------------------------------------
    traits = cn.data_loader.load_trait_table(
        TRAIT_FILE, sdm.names, ref_column="species", fill_gaps=True
    )

    # NOTE: species_traits.csv already stores growth as a MULTIPLIER
    # (1.008 - 1.219). Do NOT add 1.0 the way examples/train_policy.py does --
    # that example's trait table stores the increment instead.
    growth_rates = traits["growth_rate"].to_numpy(copy=True)

    carrying_capacity = traits["K_tetradensity"].to_numpy(copy=True)
    sensitivity = traits["sensitivity_disturbance"].to_numpy(copy=True)[:, np.newaxis]
    conservation_status = traits["conservation_status"].to_numpy(copy=True) - 1

    # AVONET Hand-Wing Index, normalised to mean = MEAN_DISPERSAL_RATE.
    # Per-species rates require cached_dispersal_matrix=None so BioEnv builds
    # D^(1/lambda_s) for each species at init.
    hwi = traits["dispersal_ability"].to_numpy(copy=True)
    dispersal_rates = hwi / np.nanmean(hwi) * MEAN_DISPERSAL_RATE
    if SCALAR_DISPERSAL:
        dispersal_rates = MEAN_DISPERSAL_RATE      # one shared matrix

    ext_risk = cn.ExtinctionRisk(
        init_status=conservation_status,
        n_classes=5,
        protect_thresholds=PROTECT_THRESHOLDS,
        # loss_thresholds left at the default beta**(1/alpha) ladder -- v2's
        # pop_decrease_threshold = 0.01 has no direct v3 equivalent (see notes).
        alpha=0.5,
    )

    env = cn.BioEnv(
        sdms=sdm,
        disturbance=disturbance,
        costs=costs,
        protection_matrix=protection,
        action_mask=action_mask,
        species_k=carrying_capacity,
        growth_rates=growth_rates,
        sensitivity_rates=sensitivity,
        dispersal_rates=dispersal_rates,
        dispersal_cutoff=DISPERSAL_WINDOW,
        ext_risk=ext_risk,
        species_traits=traits,
        device=DEVICE,
    )

    feature_extractor = cn.FeatureExtractor(
        env, feature_set=None, time_rescale=N_TIME_STEPS / 2, device=DEVICE
    )

    # [policy] nn_nodes = 3 2
    model = cn.CellNN(
        input_dim=feature_extractor.n_features,
        hidden_dim=HIDDEN_DIM,
        activation=ACTIVATION,
    )
    policy = cn.PolicyNetwork(model, seed=SEED, device=DEVICE)

    reward_objs = [
        cn.CalcRewardExtRisk(threat_weights=THREAT_WEIGHTS, device=DEVICE),
        cn.CalcRewardPersistentCost(rescaler=float(1.0 / costs.data.sum())),
    ]
    rewards = cn.Rewards(reward_obj_list=reward_objs, reward_weights=REWARD_WEIGHTS)

    budget_manager = cn.GlobalBudgetManager(
        total_target=TOTAL_TARGET_CELLS,
        cells_per_time_step=CELLS_PER_STEP,
        feature_updates_per_time_step=1,
    )

    return cn.EpisodeRunner(
        env=env,
        feature_extractor=feature_extractor,
        policy_network=policy,
        rewards=rewards,
        budget_manager=budget_manager,
        n_steps=N_TIME_STEPS,
    )


# =============================================================================
# Train
# =============================================================================


def main():
    import time

    # One runner per worker process: run_episode() mutates env state, so
    # runners cannot be shared across processes. A single runner puts the
    # trainer in sequential mode (evolution_train.py:148).
    n_runners = 1 if DEVICE == "cuda" else N_RUNNERS
    print(f"building {n_runners} runner(s)...")
    runners = [create_episode_runner() for _ in range(n_runners)]
    episode = runners[0]
    print(f"device: {DEVICE}")
    print(f"cells: {episode.env.n_cells} | species: {episode.env.n_species}")
    print(f"features: {episode.feature_extractor.n_features}")
    print(f"policy params: {len(episode.policy.get_flat_weights())}")
    print(f"already protected: {int(episode.env.protected_cells_mask.sum())}")
    print(f"blocked from action: {int(episode.env.no_action_mask.sum())}")

    print(f"perturbations/epoch: {N_PERTURBATIONS} | runners: {n_runners}")
    print(f"disturbance: {DISTURBANCE_FILE} | future: {FUTURE_DISTURBANCE_FILE or 'none'} | exclusion: {EXCLUSION_FILE}")
    print(f"network: hidden {HIDDEN_DIM} | activation {ACTIVATION}")
    print(f"cost layer: {COST_FILE}")
    print(f"dispersal: cutoff {DISPERSAL_WINDOW} cells | "
          f"{'SCALAR (uniform)' if SCALAR_DISPERSAL else 'per-species HWI'}")

    _rss = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # macOS: bytes
    print(f"parent RSS before trainer: {_rss():.0f} MB")

    trainer = cn.EvolStrategiesTrainer(
        runners,
        initial_coeffs=episode.policy.get_flat_weights(),
        scheduler=cn.LearningScheduler(initial_alpha=0.2, initial_sigma=0.2),
        epsilon_reward=0.5,
        n_perturbations=N_PERTURBATIONS,
        seed=SEED,
    )

    # EvolStrategiesTrainer pickles each runner out to its own worker process and
    # then never stores the list (it only reads list_of_env_params[0].rewards.names).
    # The parent's remaining copies are dead weight -- ~0.5 GB each -- and only
    # runners[0] is still needed, by TrainingLogger. Dropping the rest roughly
    # halves the run's footprint with no effect on results.
    if len(runners) > 1 and os.environ.get("RELEASE_SPARE_RUNNERS", "1") == "1":
        del runners[1:]
        gc.collect()
        print(f"parent RSS after releasing spare runners: {_rss():.0f} MB")

    if CALIBRATE_REWARDS:
        calib = trainer.get_reward_calibrated_weights(n_probes=10, verbose=True)
        trainer.calibrate_reward_scales(calib, verbose=True)
        trainer.save_reward_calibration(calib, RESULTS_DIR / "reward_calibration.json")
    else:
        calib = np.ones(len(episode.rewards.names))

    # Always report what the optimiser is ACTUALLY maximising, since
    # reward_weights and reward_calibration multiply rather than override.
    print("effective reward weights (reward_weights x calibration):")
    for name, w, c in zip(episode.rewards.names, episode.rewards._reward_weights.numpy(), calib):
        print(f"  {name:<18} {w:.3f} x {c:.3f} = {w * c:.3f}")

    logger = cn.algorithms.TrainingLogger(
        trainer=trainer,
        episode=episode,
        results_dir=RESULTS_DIR,
        log_file="training_log.tsv",
        weights_file="trained_weights.npy",
        plot_freq=int(os.environ.get("PLOT_FREQ", 10)),
    )

    for epoch in range(N_EPOCHS):
        t0 = time.time()
        avg_reward, summary = trainer.train_epoch()
        logger.log_epoch(epoch, avg_reward, summary, time.time() - t0)

    trainer.close()
    print(f"done -> {RESULTS_DIR}")


if __name__ == "__main__":
    main()
