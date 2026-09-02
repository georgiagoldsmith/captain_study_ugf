#!/usr/bin/env python
"""CAPTAIN 3 inference for the UGF seed-8052026 port.

The v3 analogue of predict_ugf_model.py. Rebuilds the trained environment,
loads trained_weights.npy, and writes:

    priority_score_grid.csv   235 x 578, continuous 0-1 policy score
                              (rank-normalised) -- the surface top-k ranks on
    selected_grid.csv         235 x 578, binary 1 = protected at end of episode
    *.png                     quick-look maps
    summary.txt               outcome stats + agreement with the v2 run

NOTE on the priority gradient: v2 produced 0/0.2/.../1 by averaging 5 binary
prediction runs whose weights were sampled from the last 5 training iterations.
v3's top-k is deterministic given weights (tie-break noise is 1e-7), and
TrainingLogger overwrites trained_weights.npy each epoch rather than keeping a
history -- so that ensemble cannot be reconstructed post hoc. The continuous
policy score is used instead: it is what the selection is actually ranked on,
and needs no arbitrary perturbation scale.

    python scripts/captain_predict.py
"""

from __future__ import annotations

import numpy as np
import torch

import captain as cn

from captain_train import (  # noqa: E402
    RESULTS_DIR,
    SEED,
    N_TIME_STEPS,
    TOTAL_TARGET_CELLS,
    create_episode_runner,
)

WEIGHTS = RESULTS_DIR / "trained_weights.npy"
# Compare against the v2 run for the SAME seed.
V2_TIF = (
    "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/"
    f"captain2-main/ugf_data/predictions_{SEED}/"
    f"captain_priority_prioritization_9_seed{SEED}.tif"
)


def main():
    if not WEIGHTS.exists():
        raise FileNotFoundError(f"train first -- no weights at {WEIGHTS}")

    ep = create_episode_runner()
    env = ep.env
    weights = np.load(WEIGHTS)

    # -- priority surface -----------------------------------------------------
    # Observe at t=0, BEFORE any new protection is applied: this is the state
    # the one-shot allocation actually scores.
    env.reset()
    ep.policy.set_flat_weights(weights)
    obs = ep.feature_extractor.observe(env)
    scores = ep.policy.get_scores(obs)

    blocked = env.no_action_mask
    eligible = ~blocked

    # Rank-normalise eligible cells to 0-1 so the map is comparable to v2's
    # selection-frequency gradient. Blocked cells (existing PAs + excluded)
    # get 0; the R script overlays PAs separately anyway.
    priority = torch.zeros_like(scores)
    elig_scores = scores[eligible]
    order = torch.argsort(torch.argsort(elig_scores)).float()
    priority[eligible] = order / max(len(order) - 1, 1)

    # -- run the episode for the actual selection + outcome -------------------
    res, _ = ep.run_episode(weights)
    selected = env.protected_cells_mask.cpu().numpy().astype(float)

    coords = env.sdms._coords
    shape2d = env.sdms._data_shape[1:]

    priority_grid = cn.grid_utils.reconstruct_grid(
        priority.cpu().numpy()[np.newaxis, :], coords, shape2d
    )[0]
    selected_grid = cn.grid_utils.reconstruct_grid(
        selected[np.newaxis, :], coords, shape2d
    )[0]

    # reconstruct_grid leaves non-planning cells as NaN; v2's graph_to_grid
    # wrote 0 there and the R script re-masks, so match that.
    priority_grid = np.nan_to_num(priority_grid, nan=0.0)
    selected_grid = np.nan_to_num(selected_grid, nan=0.0)

    np.savetxt(RESULTS_DIR / "priority_score_grid.csv", priority_grid, delimiter=",")
    np.savetxt(RESULTS_DIR / "selected_grid.csv", selected_grid, delimiter=",")

    cn.plots.plot_grid(
        priority_grid,
        title="CAPTAIN 3 priority score (seed 8052026)",
        outfile=RESULTS_DIR / "priority_score",
        cmap="YlGnBu",
        dpi=200,
        figsize=(6, 8),
    )
    cn.plots.plot_grid(
        selected_grid,
        title="CAPTAIN 3 selected cells (seed 8052026)",
        outfile=RESULTS_DIR / "selected_cells",
        cmap="YlGnBu",
        dpi=200,
        figsize=(6, 8),
    )

    # -- outcome stats --------------------------------------------------------
    lines = []
    lines.append(f"protected cells: {int(selected.sum())} (target {TOTAL_TARGET_CELLS})")
    lines.append(f"timesteps: {N_TIME_STEPS}")
    lines.append("")
    lines.append("extinction risk, end of episode:")
    for k, v in res["extinction_risk"].items():
        lines.append(f"  {k:<10} {v}")

    # -- agreement with the v2 run -------------------------------------------
    try:
        import rasterio

        with rasterio.open(V2_TIF) as src:
            v2 = src.read(1).astype(float)
        v2 = np.nan_to_num(v2, nan=0.0)

        # v2 grid is selection frequency; treat any selected cell as chosen
        v2_sel = v2 > 0
        v3_sel = selected_grid > 0
        inter = int((v2_sel & v3_sel).sum())
        union = int((v2_sel | v3_sel).sum())
        lines.append("")
        lines.append("agreement with v2 seed 8052026:")
        lines.append(f"  v2 selected      {int(v2_sel.sum())}")
        lines.append(f"  v3 selected      {int(v3_sel.sum())}")
        lines.append(f"  intersection     {inter}")
        lines.append(f"  Jaccard (IoU)    {inter / max(union, 1):.4f}")
        lines.append(f"  % of v2 matched  {100 * inter / max(int(v2_sel.sum()), 1):.1f}%")
    except Exception as e:  # noqa: BLE001
        lines.append("")
        lines.append(f"v2 comparison skipped: {e}")

    text = "\n".join(lines)
    (RESULTS_DIR / "summary.txt").write_text(text + "\n")
    print(text)
    print(f"\nwrote -> {RESULTS_DIR}")


if __name__ == "__main__":
    main()
