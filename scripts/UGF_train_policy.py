#!/usr/bin/env python
"""
train.py — general CAPTAIN-3 training template

Example: Train a conservation policy with Evolution Strategies.
This script demonstrates how to:
1. Load real spatial data (GeoTIFF habitat suitability maps)
2. Set up multiple parallel episode runners
3. Train a policy using Evolution Strategies
4. Log training progress

Requirements:
- Example data in DATA_DIR (see below)
- Species trait CSV file
- A reward calibration file (run calibrate_rewards.py once first)

Every optional capability (future scenarios, cost in the reward/features, pre-existing
protected areas) is behind a toggle so you can turn pieces on incrementally
instead of maintaining separate scripts per feature combination.

Each toggle is labeled with a bracketed tag, e.g. [TOGGLE: FUTURE]. Blocks
marked ### are inactive by default - can be toggled on

Run with:  uv run cpar_train_policy.py
"""

import warnings
warnings.filterwarnings("ignore", message="Sparse CSR tensor support is in beta state")

import logging
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

import captain as cn

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[logging.StreamHandler()],
)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

###############################################################################
## USER MANUAL CHANGES ##

#-------Seed keep the same for train.py and inference.py
SEED = 

# =============================================================================
# CONFIG — DATA PATHS
# =============================================================================
DATA_DIR = Path("/Users/kalena/captain3preview/ugf3_data")  # <-- change this
OUTPUTS_DIR = Path("/Users/kalena/captain3preview/outputs")

# --- present-day layers (always required) -----------------------------------#
PRESENT_SDMS_DIR = "present_habitat_suitability"        # AOH rasters folder, one per species
SPECIES_TRAIT_FILE = "species_traits.csv"               # must have a "species" column matching AOH filenames
DISTURBANCE_FILE = "environmental_layers/disturbance.tif"
COST_FILE = "environmental_layers/cost.tif"
DATA_MASK = "environmental_layers/area_mask.npy"     # build with cn.data_loader.create_mask_from_map()

#-------- future scenario layers (optional) -----------------------------------#
# Set INCLUDE_FUTURE = True below and fill these in to project SDMs/disturbance/cost forward across N_TIME_STEPS
FUTURE_SDMS_DIR = None                      # e.g. "future_sdms"
FUTURE_DISTURBANCE_FILE = None                  # e.g. "environmental_layers/future_area_swept_disturbance.tif"
FUTURE_COST_FILE = None                     # e.g. "environmental_layers/future_cost.tif"

#-------------pre-existing protected areas-------------------------------------#
# [TOGGLE: EXISTING_PAS]
EXISTING_PROTECTED_AREAS_FILE = None      # e.g. "environmental_layers/protected_areas.tif"


# =============================================================================
# CONFIG — EXPERIMENTAL TOGGLES
# =============================================================================

# [TOGGLE: SENSITIVITY] "empirical" uses trait table's sensitivity_disturbance column; "flat" overrides every species to the
# same constant, for testing how much empirical variation matters.
SENSITIVITY_MODE = "empirical"              # "empirical" or "flat"
FLAT_SENSITIVITY_VALUE = 0.5

# [TOGGLE: COST] When True, cost re-enters both the reward (CalcRewardPersistentCost) and what the policy observes
# (the "cost" feature). costs.data is always loaded either way
INCLUDE_COST = True

# [TOGGLE: FUTURE] When True, uses the FUTURE_* paths above so SDMs,
# disturbance, and cost change across N_TIME_STEPS instead of staying static.
INCLUDE_FUTURE = True

# [TOGGLE: EXISTING_PAS] When True, protection matrix initializes from
# EXISTING_PROTECTED_AREAS_FILE instead of all zeros.
INCLUDE_EXISTING_PROTECTED_AREAS = True

########### Regional-agent setup: USE_REGIONAL_AGENTS=True 
# trains one coordinated agent per region (see scripts/UGF_create_country_regions.py). 
# When False, REGION_ID selects a single global agent (None = unrestricted, or a region NAME from
# regions_tbl.csv to restrict that single agent to just that region).
USE_REGIONAL_AGENTS = False
REGION_ID = None             #ignored if true, or can set region name for just one
REGION_FILE = "environmental_layers/regions.tif"
REGION_TABLE = "environmental_layers/regions_tbl.csv"
REWARD_WEIGHTS = np.array([1.0, 1.0]) if INCLUDE_COST else np.array([1.0])

# =============================================================================
# CONFIG — TRAINING / POLICY PARAMETERS
# =============================================================================

N_EPOCHS = 100  # Number of training iterations
N_PERTURBATIONS = 6  # Number of parallel episode evaluations (sequential on GPU)
N_PARALLEL_WORKERS = 4  # Number of CPUs (if not CUDA)
N_TIME_STEPS = 50  # Time duration of each episode
TARGET_PROTECTED_CELLS_FRACTION = 0.30  # Fraction of valid cells to be protected
CELLS_PER_STEP = 1000  # number of new cells protected at each time step

# =============================================================================
# CONFIG — SPECIES / DISPERSAL PARAMETERS
# =============================================================================
# AVG_CARRYING_CAPACITY = 100  # 'individuals' per cells (* empirical relative abundance)
# Turned off since species_traits.csv has species_k called "K_tetradensity"
DISPERSAL_RATE = 0.5  # can be an array (per-species values)
DISPERSAL_WINDOW = 3
MIN_HABITAT_SUITABILITY = 0.5   # any value in (0, 1) works identically on binary AOH  #Can be overwritten below by species-specific thresholds

###############################################################################

# =============================================================================
# OUTPUT
# =============================================================================

if USE_REGIONAL_AGENTS:
    res_name = "regions"
elif REGION_ID is None:
    res_name = "global"
else:
    res_name = REGION_ID

res_name = (
    res_name
    + f"_{SENSITIVITY_MODE}"
    + ("_cost" if INCLUDE_COST else "")
    + ("_future" if INCLUDE_FUTURE else "")
)

RESULTS_DIR = OUTPUTS_DIR / (
        res_name
        + "_w"
        + "_w".join([str(r) for r in REWARD_WEIGHTS])
        + f"_p{TARGET_PROTECTED_CELLS_FRACTION}"
)
LOG_FILE = res_name + "_log.tsv"
MODEL_FILE = "trained_weights.npy"
CALIBRATION_FILE = DATA_DIR / "reward_calibration.json"
PLOT_FEATURES = False
PLOT_TRAIN_FREQ = 1  # plot intermediate protection results during training
PLOT_DATA = False
os.makedirs(RESULTS_DIR, exist_ok=True)
if PLOT_FEATURES:
    RESULTS_DIR_FEATURE_PLOTS = RESULTS_DIR / "feature_plots"
    os.makedirs(RESULTS_DIR_FEATURE_PLOTS, exist_ok=True)

# =============================================================================
# Episode setup
# =============================================================================

def create_episode_runner() -> cn.EpisodeRunner:
     """Create an episode runner with real data.

    This function loads spatial data and creates all components
    needed for one episode runner. Called once per parallel worker.
    """
     # Load present and future species distribution maps

    # --- SDMs -----------------------------------------------------------
    # [TOGGLE: FUTURE] future_dir switches between static (None) and a
    # projected scenario (confirmed optional in spatial_data.py).

     global PLOT_DATA, PLOT_FEATURES
     mask, _ = cn.data_loader.load_map(DATA_DIR / DATA_MASK)

     sdm = cn.load_spatial_data_from_dir(
        dir=DATA_DIR / PRESENT_SDMS_DIR,
        future_dir=(DATA_DIR / FUTURE_SDMS_DIR) if (INCLUDE_FUTURE and FUTURE_SDMS_DIR) else None,
        mask=mask,
        lower_bound=0,
        upper_bound=1,
        n_time_steps=N_TIME_STEPS,
        min_threshold=MIN_HABITAT_SUITABILITY,
    )
     mask, _ = cn.data_loader.load_map(DATA_DIR / DATA_MASK)

    # --- Disturbance ------------------------------------------------------
    # Load disturbance layer with predicted future change
    ### [TOGGLE: FUTURE]
     disturbance = cn.load_spatial_data(
        file=DATA_DIR / DISTURBANCE_FILE,
        future_file=(DATA_DIR / FUTURE_DISTURBANCE_FILE) if (INCLUDE_FUTURE and FUTURE_DISTURBANCE_FILE) else None,
        mask=mask,
        lower_bound=0,
        upper_bound=1,
        n_time_steps=N_TIME_STEPS,
    )

    # --- Protection matrix --------------------------------------------------
    ### [TOGGLE: EXISTING_PAS] Default: zeros = "no protected areas" (starts empty).
     if INCLUDE_EXISTING_PROTECTED_AREAS and EXISTING_PROTECTED_AREAS_FILE:
        protection = cn.load_spatial_data(
            file=DATA_DIR / EXISTING_PROTECTED_AREAS_FILE,
            future_file=None,
            mask=mask,
            lower_bound=0,
            upper_bound=1,
            n_time_steps=N_TIME_STEPS,
        )
        ###########################################
     else:
        protection = cn.SpatialData(
            data=np.zeros(disturbance.shape),
            mask=mask,
            lower_bound=0,
            upper_bound=1,
        )

    # --- Costs --------------------------------------------------------------
    ### [TOGGLE: FUTURE] Loaded regardless of INCLUDE_COST — BioEnv() requires
    # a costs argument structurally (confirmed bioenv.py). INCLUDE_COST only
    # controls whether it affects the reward/features further down.
     costs = cn.load_spatial_data(
        file=DATA_DIR / COST_FILE,
        future_file=(DATA_DIR / FUTURE_COST_FILE) if (INCLUDE_FUTURE and FUTURE_COST_FILE) else None,
        mask=mask,
        lower_bound=0,
        upper_bound=1,
        n_time_steps=N_TIME_STEPS,
    )

    # Regional-agent setup: build per-region masks/targets, or a single
    # action_mask restricting a global agent to one region.
     TARGET_PROTECTED_CELLS = int(TARGET_PROTECTED_CELLS_FRACTION * np.nansum(mask))

     if USE_REGIONAL_AGENTS:
        tmp, _ = cn.data_loader.load_map(DATA_DIR / REGION_FILE)
        region_tbl = pd.read_csv(DATA_DIR / REGION_TABLE)

        regional_totals, per_step_regional_targets, region_masks = {}, {}, {}
        for region_num, name in zip(region_tbl["REGION_ID"], region_tbl["NAME"], strict=True):
            r_mask = cn.SpatialData(
                data=tmp == region_num, mask=mask, lower_bound=0, upper_bound=1
            )
            target = int(TARGET_PROTECTED_CELLS_FRACTION * r_mask.data.sum())
            regional_totals[name] = target
            per_step_regional_targets[name] = int(
                CELLS_PER_STEP * (r_mask.data.sum() / np.nansum(mask))
            )
            region_masks[name] = r_mask._nonzero_cells_mask
            print(
                f"Region {region_num} ({name}): size={r_mask.data.sum()}, "
                f"target={target}, per_step={per_step_regional_targets[name]}"
            )

        region_mask = None  # regional mode does not restrict env.action_mask

     elif REGION_ID is None:
        region_mask = None

     else:
        tmp, _ = cn.data_loader.load_map(DATA_DIR / REGION_FILE)
        region_tbl = pd.read_csv(DATA_DIR / REGION_TABLE)
        match = region_tbl.loc[region_tbl["NAME"] == REGION_ID, "REGION_ID"]
        if match.empty:
            raise ValueError(f"REGION_ID {REGION_ID!r} not found in {REGION_TABLE}")
        region_num = match.iloc[0]

        region_mask = cn.SpatialData(
            data=tmp != region_num, mask=mask, lower_bound=0, upper_bound=1
        )
        TARGET_PROTECTED_CELLS = int(
            TARGET_PROTECTED_CELLS_FRACTION * (1 - region_mask.data).sum()
            + protection.data.sum()
        )

    # --- Species traits -------------------------------------------------
    # simple imputation of missing data (could be replaced e.g. RF imputation)
     traits = cn.data_loader.load_trait_table(
        DATA_DIR / SPECIES_TRAIT_FILE, sdm.names, ref_column="species", fill_gaps=True
    )
    # Per-species minimum habitat suitability (overrides the scalar
    # MIN_HABITAT_SUITABILITY set above; below this threshold a cell doesn't
    # contribute to that species' carrying capacity)
    #sdm.reset_threshold(traits["min_habitat_suitability"].to_numpy())

# extract parameters for simulation
    # [TOGGLE: SENSITIVITY]
     if SENSITIVITY_MODE == "empirical":
        sensitivity = traits["sensitivity_disturbance"].to_numpy(copy=True)[:, np.newaxis]
     elif SENSITIVITY_MODE == "flat":
        sensitivity = np.full((sdm.shape[0], 1), FLAT_SENSITIVITY_VALUE, dtype=np.float32)
     else:
        raise ValueError("SENSITIVITY_MODE must be 'empirical' or 'flat'")

    # Growth rate — confirmed BioEnv warns if growth_rates.min() < 1.0.
    # Use directly if your column is already a multiplicative factor >= 1;
    # otherwise adjust (e.g. "+ 1.0") before passing in.
     growth_rates = traits["growth_rate"].to_numpy(copy=True)

    # Carrying capacity — empirical
     species_k = traits["K_tetradensity"].to_numpy(copy=True)
    # or for estimated proxy: carrying_capacity = AVG_CARRYING_CAPACITY / traits["conservation_status"].to_numpy(copy=True)

    # conservation_status -> Initial extinction risk from conservation status
     conservation_status = traits["conservation_status"].to_numpy(copy=True) - 1
     ext_risk = cn.ExtinctionRisk(init_status=conservation_status, n_classes=5, alpha=0.5)

    # Dispersal — empirical per-species array from trait table. bioenv.py: dispersal_rates
    # accepts scalar or per-species array; passing an array
    # (with no cached_dispersal_matrix) makes BioEnv build species-specific
    # matrices internally.
     dispersal_rates = traits["dispersal_ability"].to_numpy(copy=True)

    ### trait-weighted reward ########## -----------
    # [TOGGLE: SPECIES_VALUE_REWARD] Confirmed from rewards.py:
    # CalcRewardSpecieValue reads env._species_traits (a DataFrame set via
    # BioEnv(species_traits=...)) to add a reward term = dot product of a
    # chosen trait column with species abundance. NOT read by
    # FeatureExtractor (confirmed feature_extractor.py) — this only affects
    # the reward, not what the policy observes. Uncomment both this table
    # and the CalcRewardSpecieValue line in the rewards section below to use.
    ### species_traits_table = traits[["species", "biomass", "rel_abundance"]]
    ###########################################################################

     env = cn.BioEnv(
        sdms=sdm,
        disturbance=disturbance,
        costs=costs,
        protection_matrix=protection,
        species_k=species_k,
        growth_rates=growth_rates,
        sensitivity_rates=sensitivity,
        dispersal_rates=dispersal_rates,
        dispersal_cutoff=DISPERSAL_WINDOW,
        ext_risk=ext_risk,
        action_mask=region_mask,  # None, or restricts actions to one region
        device=DEVICE,
     )

    # Protection budget as a fraction of valid cells (env.n_cells confirmed
    # in bioenv.py: count of non-masked cells).
    # target_protected_cells = int(TARGET_PROTECTED_FRACTION * env.n_cells)

    # --- Feature set ----------------------------------------------------
    # [TOGGLE: COST] Confirmed default_feature_set in feature_extractor.py:
    # ["time", "disturbance", "disturbance_conv", "species_richness",
    #  "total_population", "current_ext_risk", "cost",
    #  "protection_matrix", "protection_matrix_conv"]
    # feature_set=None uses that full list. Species-only excludes "cost" so
    # the policy doesn't observe it at all (not just absent from reward).
     if INCLUDE_COST:
        feature_set = None  # full default set, includes "cost"
     else:
        feature_set = [
            "time", "disturbance", "disturbance_conv", "species_richness",
            "total_population", "current_ext_risk",
            "protection_matrix", "protection_matrix_conv",
        ]
    # Create agent components
     feature_extractor = cn.FeatureExtractor(
        env,
        feature_set=feature_set,  # Use default feature set (cost toggle) (can be customized)
        time_rescale=N_TIME_STEPS / 2,
        device=DEVICE,
     )
     if PLOT_FEATURES:
        feature_extractor.plot_features(
            env, rescale=False, outdir=RESULTS_DIR_FEATURE_PLOTS
        )
        PLOT_FEATURES = False

     env.ext_risk.species_per_class(env.current_ext_risk)

     model = cn.CellNN(input_dim=feature_extractor.n_features, hidden_dim=16)
     if USE_REGIONAL_AGENTS:
        policy = cn.RegionalPolicyNetwork(model, seed=SEED, device=DEVICE)
     else:
        policy = cn.PolicyNetwork(model, seed=SEED, device=DEVICE)

    # --- Rewards ----------------------------------------------------------
    # [TOGGLE: COST] Extinction-risk term always included. Cost penalty
    # term added only when INCLUDE_COST=True.
     reward_obj_list = [
        cn.CalcRewardExtRisk(threat_weights=np.array([1, 0, -8, -16, -32]), device=DEVICE),
     ]

     if INCLUDE_COST:
        reward_obj_list.append(
            cn.CalcRewardPersistentCost(rescaler=float(1.0 / costs.data.sum()))
        )

     rewards = cn.Rewards(
        reward_obj_list=reward_obj_list,
        reward_weights=REWARD_WEIGHTS,
     )

     if PLOT_DATA:
        cn.plots.plot_grid(
            np.sum(env.reconstruct_h_grid > 1, axis=0),
            title="species richness",
            outfile=RESULTS_DIR_FEATURE_PLOTS / "species_richness",
        )
        PLOT_DATA = False

# Create episode runner
     if USE_REGIONAL_AGENTS:
        budget_manager = cn.RegionalBudgetManager(
            masks=region_masks,
            total_targets=regional_totals,
            cells_per_time_step=per_step_regional_targets,
        )
     else:
        budget_manager = cn.GlobalBudgetManager(
            total_target=TARGET_PROTECTED_CELLS,
            cells_per_time_step=CELLS_PER_STEP,
            feature_updates_per_time_step=1,
        )

     ep = cn.EpisodeRunner(
        env=env,
        feature_extractor=feature_extractor,
        policy_network=policy,
        rewards=rewards,
        n_steps=N_TIME_STEPS,
        budget_manager=budget_manager,
    )

    # plot end features
     if PLOT_FEATURES:
        ep.run_episode()
        feature_extractor.plot_features(
            env, rescale=False, outdir=RESULTS_DIR_FEATURE_PLOTS
        )

     return ep


# =============================================================================
# Main Training loop
# =============================================================================

def main():
    print("=" * 60)
    print("CAPTAIN-3 Training")
    print(f"  sensitivity={SENSITIVITY_MODE}  cost={INCLUDE_COST}  "
          f"future={INCLUDE_FUTURE}  existing_PAs={INCLUDE_EXISTING_PROTECTED_AREAS}")
    print("=" * 60)

    if not DATA_DIR.exists():
        raise FileNotFoundError("Update DATA_DIR to point to your data.")

    print(f"Device: {DEVICE}")

    if DEVICE == "cuda":
        # GPU mode: single runner, sequential perturbations (avoids pickling CUDA tensors)
        print(
            f"Creating 1 episode runner on {DEVICE} ({N_PERTURBATIONS} perturbations)..."
        )
        episode_runners = [create_episode_runner()]
    else:
        # CPU mode: multiple runners with multiprocessing pool
        print(f"Creating {N_PARALLEL_WORKERS} parallel episode runners on CPU...")
        episode_runners = [create_episode_runner() for _ in range(N_PARALLEL_WORKERS)]

    # Get reference to first runner for logging
    episode = episode_runners[0]
    print(f"  Grid: {episode.env.n_cells} cells, {episode.env.n_species} species")
    print(f"  Features: {episode.feature_extractor.n_features}")
    print(f"  Policy parameters: {len(episode.policy.get_flat_weights())}")


    # Create trainer
    trainer = cn.EvolStrategiesTrainer(
        episode_runners,
        initial_coeffs=episode.policy.get_flat_weights(),
        scheduler=cn.LearningScheduler(initial_alpha=0.05, initial_sigma=0.3),
        epsilon_reward=0.75,
        n_perturbations=N_PERTURBATIONS,
        seed=SEED,
    )


    # Load reward calibration. This is always computed once against the
    # global agent (see calibrate_rewards.py) and reused here regardless of
    # USE_REGIONAL_AGENTS/REGION_ID, so reward values stay comparable across
    # global/regional/single-region runs.
    if not CALIBRATION_FILE.exists():
        raise FileNotFoundError(
            f"{CALIBRATION_FILE} not found — run calibrate_rewards.py first."
        )
    trainer.load_reward_calibration(CALIBRATION_FILE, verbose=True)

    logger = cn.algorithms.TrainingLogger(
        trainer=trainer,
        episode=episode,
        results_dir=RESULTS_DIR,
        log_file=LOG_FILE,
        weights_file=MODEL_FILE,
        plot_freq=PLOT_TRAIN_FREQ,
    )

    # Training loop
    print(f"\nTraining for {N_EPOCHS} epochs...")
    print("-" * 60)

    t_start = time.time()

    for epoch in range(N_EPOCHS):
        t0 = time.time()
        avg_reward, summary = trainer.train_epoch()
        # Log progress
        logger.log_epoch(epoch, avg_reward, summary, time.time() - t0)

    # Summary
    print("-" * 60)
    print(f"Training complete in {time.time() - t_start:.1f}s")
    print(f"Log saved to: {logger.log_path}")
    print(f"Weights saved to: {logger.weights_path}")

    # Plot reward across epochs
    reward_plot_path = RESULTS_DIR / "reward_over_training.png"
    cn.plots.plot_rl_rewards(
        logger.log_path,
        title=f"Reward over training ({res_name})",
        outfile=reward_plot_path,
        reward_col="avg_reward",
    )
    print(f"Reward plot saved to: {reward_plot_path}")

    # Cleanup
    trainer.close()

if __name__ == "__main__":
    main()
