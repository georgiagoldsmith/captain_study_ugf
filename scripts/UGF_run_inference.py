#!/usr/bin/env python
"""
run_inference.py — general CAPTAIN-3 inference template

Loads a policy trained with the matching train.py (Ra_train_policy.py) toggle
structure and runs it to produce a protection matrix + extinction-risk plots.

Requirements:
- Same DATA_DIR / data layout used for training
- Same species trait CSV
- trained_weights.npy produced by a matching training run

CRITICAL: every CONFIG — TOGGLES / PARAMETERS value below must match the
training run whose weights you're loading. The policy's input feature
count and dispersal/growth/sensitivity setup all depend on these — a
mismatch will either raise a shape error or silently give meaningless
output. This is not enforced automatically by the code; it's on you to
keep these two scripts' configs in sync per experiment.

Each toggle is labeled with a bracketed tag, e.g. [TOGGLE: FUTURE]. Blocks
marked ### are inactive by default — can be toggled on.

Run with:  uv run Ra_run_inference.py
"""

import warnings
warnings.filterwarnings("ignore", message="Sparse CSR tensor support is in beta state")

import logging
import os
import shutil
import subprocess
from datetime import datetime
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
# Policy weights are overwritten by the loaded trained_weights.npy before any
# action is taken, so this seed has no effect on inference output — it's
# kept only for API parity with the training script.
SEED =

# =============================================================================
# CONFIG — DATA PATHS
# (must match the training run's config exactly — see module docstring)
# =============================================================================
DATA_DIR = Path("/Users/kalena/captain3preview/ugf3_data")  # <-- change this
OUTPUTS_DIR = Path("/Users/kalena/captain3preview/outputs")

# --- present-day layers (always required) -----------------------------------#
PRESENT_SDMS_DIR = "present_habitat_suitability"        # SDMS or AOH rasters, one per species
SPECIES_TRAIT_FILE = "species_traits.csv"               # must have a "species" column matching filenames
DISTURBANCE_FILE = "environmental_layers/disturbance.tif"
COST_FILE = "environmental_layers/cost.tif"          # [TOGGLE: COST] — still loaded even if unused
DATA_MASK = "environmental_layers/area_mask.npy"     # build with cn.data_loader.create_mask_from_map()

#-------- future scenario layers (optional) -----------------------------------#
# [TOGGLE: FUTURE] Must match whatever the training run used. If these stay
# None, delta_sdm/delta stays None and the layer doesn't change over time (present-day-only scenario)
FUTURE_SDMS_DIR = None                      # e.g. "future_sdms"
FUTURE_DISTURBANCE_FILE = None             # e.g. "environmental_layers/future_area_swept_disturbance.tif"
FUTURE_COST_FILE = None                     # e.g. "environmental_layers/future_cost.tif"


#-------------pre-existing protected areas-------------------------------------#
# [TOGGLE: EXISTING_PAS] Must match the training run's starting protection state.
EXISTING_PROTECTED_AREAS_FILE = None

# =============================================================================
# CONFIG — EXPERIMENTAL TOGGLES
# (must match the training run — these change env/feature/reward structure,
# not just numeric values, so a mismatch changes n_features / policy shape)
# =============================================================================

# [TOGGLE: SENSITIVITY] "empirical" uses your trait table's
# sensitivity_disturbance column; "flat" overrides every species to the same constant.
SENSITIVITY_MODE = "empirical"              # "empirical" or "flat"
FLAT_SENSITIVITY_VALUE = 0.5

# [TOGGLE: COST] When True, cost re-enters both the reward (irrelevant here,
# inference uses NoRewards) and what the policyobserves (the "cost" feature)
# This must match training, since it changes feature_extractor.n_features and therefore the policy's input shape. 
# costs.data is still loaded either way — BioEnv() requires a costs argument.
INCLUDE_COST = False

# [TOGGLE: FUTURE] When True, uses the FUTURE_* paths above so SDMs,
# disturbance, and cost change across N_TIME_STEPS instead of staying static.
INCLUDE_FUTURE = False

# [TOGGLE: EXISTING_PAS] When True, protection matrix initializes from
# EXISTING_PROTECTED_AREAS_FILE instead of all zeros.
INCLUDE_EXISTING_PROTECTED_AREAS = True

########### Regional-agent setup: must match training. USE_REGIONAL_AGENTS=True
# loads a coordinated per-region policy (see scripts/create_toy_regions.py).
# When False, REGION_ID selects a single global agent (None = unrestricted,
# or a region NAME from regions_tbl.csv to restrict that single agent to
# just that region).
USE_REGIONAL_AGENTS = False
REGION_ID = None             # ignored if True, or can set region name for just one
REGION_FILE = "environmental_layers/regions.tif"
REGION_TABLE = "environmental_layers/regions_tbl.csv"
REWARD_WEIGHTS = np.array([1.0, 1.0]) if INCLUDE_COST else np.array([1.0])

# =============================================================================
# CONFIG — INFERENCE / POLICY PARAMETERS
# (N_TIME_STEPS, TARGET_PROTECTED_CELLS_FRACTION, CELLS_PER_STEP must match
# training — no N_EPOCHS / N_PERTURBATIONS here, inference doesn't train)
# =============================================================================

N_TIME_STEPS = 50  # Time duration of each episode
TARGET_PROTECTED_CELLS_FRACTION = 0.30  # Fraction of valid cells to be protected
CELLS_PER_STEP = 1000  # number of new cells protected at each time step

# =============================================================================
# CONFIG — SPECIES / DISPERSAL PARAMETERS
# =============================================================================
DISPERSAL_RATE = 0.5  # kept for parity with train.py; not directly used below —
                       # dispersal is built from the trait table's
                       # dispersal_ability column instead (see create_episode_runner)
DISPERSAL_WINDOW = 3
MIN_HABITAT_SUITABILITY = 0.5   # any value in (0, 1) works identically on binary AOH

# =============================================================================
# OUTPUT
# (identical naming formula to train.py's RESULTS_DIR, so inference always
# finds the weights produced by the matching training config)
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
MODEL_FILE = "trained_weights.npy"
TRAINED_MODEL = RESULTS_DIR / MODEL_FILE
RES_DIR = RESULTS_DIR / "inference"
os.makedirs(RES_DIR, exist_ok=True)

PLOT_FEATURES = False
PLOT_DATA = True


if PLOT_FEATURES:
    RESULTS_DIR_FEATURE_PLOTS = RES_DIR / "feature_plots"
    os.makedirs(RESULTS_DIR_FEATURE_PLOTS, exist_ok=True)

# =============================================================================
# Episode setup
# =============================================================================

def create_episode_runner() -> cn.EpisodeRunner:
    """Create an episode runner with real data and load trained weights.

    Mirrors train.py's create_episode_runner(), with two differences:
    rewards uses cn.NoRewards() (inference doesn't need a reward signal),
    and the policy's weights are loaded from TRAINED_MODEL before return.
    """
    global PLOT_DATA, PLOT_FEATURES

    mask, _ = cn.data_loader.load_map(DATA_DIR / DATA_MASK)

    # --- SDMs -----------------------------------------------------------
    # [TOGGLE: FUTURE] future_dir switches between static (None) and a
    # projected scenario (confirmed optional in spatial_data.py).
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
    # controls whether it affects the feature set further down.
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
    traits = cn.data_loader.load_trait_table(
        DATA_DIR / SPECIES_TRAIT_FILE, sdm.names, ref_column="species", fill_gaps=True
    )

    # [TOGGLE: SENSITIVITY]
    if SENSITIVITY_MODE == "empirical":
        sensitivity = traits["sensitivity_disturbance"].to_numpy(copy=True)[:, np.newaxis]
    elif SENSITIVITY_MODE == "flat":
        sensitivity = np.full((sdm.shape[0], 1), FLAT_SENSITIVITY_VALUE, dtype=np.float32)
    else:
        raise ValueError("SENSITIVITY_MODE must be 'empirical' or 'flat'")

    # Growth rate — confirmed BioEnv warns if growth_rates.min() < 1.0.
    growth_rates = traits["growth_rate"].to_numpy(copy=True)

    # Carrying capacity — empirical
    species_k = traits["K_tetradensity"].to_numpy(copy=True)

    # conservation_status -> Initial extinction risk from conservation status
    conservation_status = traits["conservation_status"].to_numpy(copy=True) - 1
    ext_risk = cn.ExtinctionRisk(init_status=conservation_status, n_classes=5, alpha=0.5)

    # Dispersal — empirical per-species array from trait table.
    dispersal_rates = traits["dispersal_ability"].to_numpy(copy=True)

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

    # --- Feature set ----------------------------------------------------
    # [TOGGLE: COST] Confirmed default_feature_set in feature_extractor.py:
    # ["time", "disturbance", "disturbance_conv", "species_richness",
    #  "total_population", "current_ext_risk", "cost",
    #  "protection_matrix", "protection_matrix_conv"]
    # feature_set=None uses that full list. Species-only excludes "cost" so
    # the policy doesn't observe it at all (not just absent from reward).
    # This MUST match the training run — it changes feature_extractor.n_features.
    if INCLUDE_COST:
        feature_set = None  # full default set, includes "cost"
    else:
        feature_set = [
            "time", "disturbance", "disturbance_conv", "species_richness",
            "total_population", "current_ext_risk",
            "protection_matrix", "protection_matrix_conv",
        ]

    feature_extractor = cn.FeatureExtractor(
        env,
        feature_set=feature_set,
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

    # Load trained weights — must come from a run with the same toggles above,
    # otherwise get_flat_weights()/set_flat_weights() shapes won't match.
    if not TRAINED_MODEL.exists():
        raise FileNotFoundError(
            f"{TRAINED_MODEL} not found — check TRAINED_MODEL/RESULTS_DIR "
            f"naming matches the training run's config, or run training first."
        )
    policy.set_flat_weights(np.load(TRAINED_MODEL))

    # Inference doesn't need a reward signal — NoRewards() confirmed in
    # reward_aggregator.py.
    rewards = cn.NoRewards()

    if PLOT_DATA:
        cn.plots.plot_grid(
            np.sum(env.reconstruct_h_grid > 1, axis=0),
            title="species richness",
            outfile=RES_DIR / "species_richness",
        )
        PLOT_DATA = False

    # Create budget manager
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
        save_protection_history=True,
    )

    return ep


# =============================================================================
# Main inference run
# =============================================================================

def main():
    print("=" * 60)
    print("CAPTAIN-3 Inference")
    print(f"  sensitivity={SENSITIVITY_MODE}  cost={INCLUDE_COST}  "
          f"future={INCLUDE_FUTURE}  existing_PAs={INCLUDE_EXISTING_PROTECTED_AREAS}")
    print(f"  model: {TRAINED_MODEL}")
    print("=" * 60)

    if not DATA_DIR.exists():
        raise FileNotFoundError("Update DATA_DIR to point to your data.")

    print(f"Device: {DEVICE}")

    ep = create_episode_runner()
    print(f"  Grid: {ep.env.n_cells} cells, {ep.env.n_species} species")
    print(f"  Features: {ep.feature_extractor.n_features}")
    print(f"  Policy parameters: {len(ep.policy.get_flat_weights())}")

    weights = np.load(TRAINED_MODEL)
    res, _ = ep.run_episode(weights)
    protection_grid_snapshot = ep.env.protection_matrix.reconstruct_grid[0].copy()

    cn.plots.plot_grid(
        ep.env.protection_matrix.reconstruct_grid[0],
        title="protection matrix",
        outfile=RES_DIR / "protection_matrix",
        dpi=300,
        figsize=(6, 8),
    )

    history = (res["protection_history"] > 0).int() * (
            1 + res["protection_history"].max() - res["protection_history"]
    )
    mask, _ = cn.data_loader.load_map(DATA_DIR / DATA_MASK)
    protection_res = cn.SpatialData(
        data=np.zeros(ep.env.disturbance.shape),
        mask=mask,
        lower_bound=0,
        upper_bound=1,
    )
    protection_res._data += history

    cn.plots.plot_grid(
        protection_res.reconstruct_grid[0],
        title="protection matrix through time",
        outfile=RES_DIR / "protection_matrix_through_time",
        dpi=300,
        figsize=(6, 8),
        cmap="viridis",
    )

    # present extinction risks
    cn.plots.plot_extinction_risk(
        ep.env.ext_risk.init_status,
        labels=["LC", "NT", "VU", "EN", "CR"],
        outfile=RES_DIR / "Extinction_risk",
        title="Present extinction risk",
        dpi=200,
    )

    # (predicted) future extinction risks, with protection
    cn.plots.plot_extinction_risk(
        ep.env.current_ext_risk,
        labels=["LC", "NT", "VU", "EN", "CR"],
        outfile=RES_DIR / "Extinction_risk_future",
        title="Future extinction risk (protection)",
        dpi=200,
    )

    # run without protection for comparison
    # (NoBudgetManager's step context isn't compatible with
    # RegionalPolicyNetwork, which requires region_masks/region_k rather than
    # n_cells, so skip this comparison run in regional mode.)
    if not USE_REGIONAL_AGENTS:
        ep_noprot = cn.EpisodeRunner(
            env=ep.env,
            feature_extractor=ep.feature_extractor,
            policy_network=ep.policy,
            rewards=cn.NoRewards(),
            n_steps=N_TIME_STEPS,
            budget_manager=cn.NoBudgetManager(),
            save_protection_history=True,
        )

        res_noprot, _ = ep_noprot.run_episode(weights)

        cn.plots.plot_extinction_risk(
            ep_noprot.env.current_ext_risk,
            labels=["LC", "NT", "VU", "EN", "CR"],
            outfile=RES_DIR / "Extinction_risk_future_no_protection",
            title="Future extinction risk (no protection)",
            dpi=200,
        )
    else:
        print(
            "Skipping no-protection comparison run in regional multi-agent mode "
            "(RegionalPolicyNetwork requires region_masks/region_k, not n_cells)."
        )


