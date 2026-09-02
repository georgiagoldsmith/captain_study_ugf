"""Connectivity reward for CAPTAIN 3.

v3 ships no connectivity term (rewards.py has only CalcReward,
CalcRewardPersistentCost, CalcRewardExtRisk and CalcRewardSpecieValue), so this
port follows CAPTAIN 2's `BioDivEnv.calc_connectivity_reward` -- the same term
Horn et al. 2026 enable via `rewards_variables = species_risk connectivity cost`.

v2's formulation:

    c = convolve(protection_matrix, ones(3,3))   # protected count in each 3x3
    c[c == 1] = 0                                # a lone protected cell scores 0
    reward = sum(c) / n_protected * min(10, n_species)

Reimplemented here as a sparse adjacency matvec over v3's flattened valid-cell
representation, which is equivalent and cheap (~400k nonzeros, one matvec per
step). Two deliberate departures:

  * Normalised by 9 so a single step returns [0, 1]: 1.0 means every protected
    cell sits in a fully protected 3x3 block, 0.0 means every one is isolated.
    v2's `* min(10, n_species)` scaling is dropped -- magnitude belongs in
    `rescaler`, not baked into the metric.
  * Edge cells have fewer than 9 neighbours, so a solution hugging the study-area
    boundary cannot reach exactly 1.0. This matches v2, which used
    mode='constant' (zero padding) for the same reason.

Because protection here is allocated one-shot at t=0 and then static, this
returns the same value every step, so the episode total is n_steps x the
per-step score. Set rescaler to about 1/n_steps to keep the episode total on the
same order as the extinction-risk (~13) and cost (~-9) terms.
"""

from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors

import captain as cn


class CalcRewardConnectivity(cn.CalcReward):
    """Reward spatial compactness of the protected set."""

    def __init__(
        self,
        env,
        name: str = "connectivity",
        rescaler: float = 1.0,
        positive: bool = True,
        radius: int = 1,
    ):
        """
        Args:
            env: BioEnv -- used once, at construction, for the cell coordinates.
            rescaler: scaling applied to the per-step score.
            radius: Chebyshev radius; 1 gives v2's 3x3 block.
        """
        super().__init__(name, rescaler, positive)

        pts = np.column_stack(env.sdms._coords)
        nn = NearestNeighbors(radius=radius, metric="chebyshev").fit(pts)
        # mode="connectivity" -> 1 for each neighbour within radius, self included
        self._adj = nn.radius_neighbors_graph(
            pts, radius=radius, mode="connectivity"
        ).tocsr()
        # a full block is (2r+1)^2 cells; used to normalise to [0, 1]
        self._block = float((2 * radius + 1) ** 2)

    def calc_reward(self, env) -> float:
        p = env.protected_cells_mask.detach().cpu().numpy().astype(np.float32)
        n_prot = float(p.sum())
        if n_prot == 0:
            return 0.0

        # protected-cell count within each cell's neighbourhood
        c = self._adj @ p
        # v2: a neighbourhood containing exactly one protected cell earns nothing,
        # so isolated selections are not rewarded merely for existing
        c[c == 1.0] = 0.0

        return float(c.sum() / n_prot / self._block) * self._rescaler
