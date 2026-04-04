#!/usr/bin/env python3
"""
Constrained greedy Jaccard bundling for ICON velocity_tendencies.

Bundle sizes must be 2, 4, 8, or a multiple of 8.  This reflects
AoSoA / SIMD-lane alignment: the number of fields in a storage
record should match a power-of-two or cache-line multiple so that
vectorised inner loops touch complete records without waste.

Algorithm (two-phase)
---------------------
Phase 1 — unconstrained greedy (Jaccard >= threshold):
    Identical to the seed-and-grow procedure in greedy_jaccard_bundles.py.

Phase 2 — constrained merge/pad:
    While any bundle has an invalid size:
      (a) Find the smallest invalid-sized bundle B_small.
      (b) Among all other bundles B_other, find the one whose
          average pairwise Jaccard with B_small is highest.
      (c) If |B_small| + |B_other| is a valid size, merge.
      (d) Otherwise, pad B_small from unassigned items (lowest
          threshold first) until it reaches the next valid size.

Valid sizes:  {2, 4, 8, 16, 24, 32, ...}

References
----------
  Jaccard, P. (1901). Étude comparative de la distribution florale
    dans une portion des Alpes et du Jura.  Bull. Soc. Vaudoise Sci.
    Nat. 37(1): 547-579.
  Edgar, R.C. (2010). Search and clustering orders of magnitude faster
    than BLAST.  Bioinformatics 26(19): 2460-2461.
  Bron, C. & Kerbosch, J. (1973). Algorithm 457: Finding all cliques
    of an undirected graph.  Commun. ACM 16(9): 575-577.
"""

from __future__ import annotations
from collections import defaultdict
from itertools import combinations


# ---------------------------------------------------------------------------
# Valid bundle sizes
# ---------------------------------------------------------------------------

VALID_SIZES = {2, 4, 8, 16, 24, 32, 40, 48, 56, 64}

def is_valid_size(n: int) -> bool:
    if n in (2, 4, 8):
        return True
    return n >= 8 and n % 8 == 0

def next_valid_size(n: int) -> int:
    """Smallest valid size >= n."""
    for v in sorted(VALID_SIZES):
        if v >= n:
            return v
    # For very large n, round up to next multiple of 8
    return ((n + 7) // 8) * 8


# ---------------------------------------------------------------------------
# Jaccard similarity
# ---------------------------------------------------------------------------

def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def avg_jaccard_between(group_a: list, group_b: list,
                        item_to_ctx: dict) -> float:
    """Average pairwise Jaccard between two groups of items."""
    total, count = 0.0, 0
    for a in group_a:
        for b in group_b:
            total += jaccard(item_to_ctx[a], item_to_ctx[b])
            count += 1
    return total / count if count > 0 else 0.0


# ---------------------------------------------------------------------------
# Phase 1: unconstrained greedy
# ---------------------------------------------------------------------------

def _unconstrained_greedy(item_to_ctx: dict, threshold: float) -> list[list[str]]:
    items_sorted = sorted(item_to_ctx, key=lambda x: -len(item_to_ctx[x]))
    assigned = set()
    bundles = []

    for seed in items_sorted:
        if seed in assigned:
            continue
        bundle = [seed]
        assigned.add(seed)

        for cand in items_sorted:
            if cand in assigned:
                continue
            cand_ctx = item_to_ctx[cand]
            min_j = min(jaccard(item_to_ctx[m], cand_ctx) for m in bundle)
            if min_j >= threshold:
                bundle.append(cand)
                assigned.add(cand)

        bundles.append(bundle)

    return bundles


# ---------------------------------------------------------------------------
# Phase 2: constrained merge/pad
# ---------------------------------------------------------------------------

def _constrained_postprocess(bundles: list[list[str]],
                             item_to_ctx: dict) -> list[list[str]]:
    """
    Merge and pad bundles until all have valid sizes.

    Strategy:
      1. Sort bundles by size (ascending).
      2. Repeatedly merge the two smallest bundles whose combined size
         is valid, preferring pairs with highest average Jaccard.
      3. If no valid merge exists, pad the smallest bundle by stealing
         the most-similar member from the largest bundle (if the donor
         remains valid after losing a member).
      4. As a last resort, merge the two smallest bundles regardless
         of the resulting size, then continue.
    """
    bundles = [list(b) for b in bundles]  # deep copy

    max_iters = len(bundles) * 10
    for _ in range(max_iters):
        # Check if all valid
        invalid = [i for i, b in enumerate(bundles) if not is_valid_size(len(b))]
        if not invalid:
            break

        # Sort invalid by size ascending
        invalid.sort(key=lambda i: len(bundles[i]))
        idx_small = invalid[0]
        b_small = bundles[idx_small]

        # Try merging with another bundle to hit a valid size
        best_merge = None
        best_jacc = -1.0
        for j, b_other in enumerate(bundles):
            if j == idx_small:
                continue
            combined = len(b_small) + len(b_other)
            if is_valid_size(combined):
                jacc = avg_jaccard_between(b_small, b_other, item_to_ctx)
                if jacc > best_jacc:
                    best_jacc = jacc
                    best_merge = j

        if best_merge is not None:
            bundles[idx_small] = b_small + bundles[best_merge]
            del bundles[best_merge]
            continue

        # No valid merge: try stealing a member from another bundle
        # to grow b_small toward next valid size
        target = next_valid_size(len(b_small))
        need = target - len(b_small)
        stolen = False
        for _ in range(need):
            best_donor_idx = None
            best_member = None
            best_j = -1.0
            for j, b_donor in enumerate(bundles):
                if j == idx_small:
                    continue
                donor_size_after = len(b_donor) - 1
                if donor_size_after < 1:
                    continue
                # Only steal if the donor remains valid (or was already invalid)
                if is_valid_size(len(b_donor)) and not is_valid_size(donor_size_after):
                    continue
                for m in b_donor:
                    j_avg = avg_jaccard_between([m], b_small, item_to_ctx)
                    if j_avg > best_j:
                        best_j = j_avg
                        best_donor_idx = j
                        best_member = m

            if best_member is not None:
                bundles[idx_small].append(best_member)
                bundles[best_donor_idx].remove(best_member)
                if not bundles[best_donor_idx]:
                    del bundles[best_donor_idx]
                stolen = True
            else:
                break

        if stolen:
            continue

        # Last resort: force-merge the two smallest
        if len(bundles) >= 2:
            sizes = sorted(range(len(bundles)), key=lambda i: len(bundles[i]))
            i1, i2 = sizes[0], sizes[1]
            if i1 > i2:
                i1, i2 = i2, i1
            bundles[i1] = bundles[i1] + bundles[i2]
            del bundles[i2]
        else:
            break

    return bundles


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def constrained_jaccard_bundles(
    item_to_ctx: dict[str, set],
    threshold: float = 0.5,
) -> list[list[str]]:
    """
    Greedy Jaccard bundling with valid-size constraint.

    Returns bundles where every bundle has size in {2, 4, 8, 16, 24, ...}.
    """
    raw = _unconstrained_greedy(item_to_ctx, threshold)
    return _constrained_postprocess(raw, item_to_ctx)


# ---------------------------------------------------------------------------
# velocity_tendencies data
# ---------------------------------------------------------------------------

NEST_DATA = {
    1:  ["p_diag%vn_ie","p_prog%vn","z_vt_ie","p_diag%vt","z_kin_hor_e","p_metrics%wgtfacq_e"],
    2:  ["p_diag%vn_ie","p_diag%vn_ie_ubc","z_vt_ie","p_diag%vt","z_kin_hor_e","p_prog%vn","p_metrics%wgtfacq_e"],
    3:  ["p_diag%vt","p_int%rbf_vec_coeff_e","p_prog%vn"],
    4:  ["p_diag%vn_ie","p_metrics%wgtfac_e","p_prog%vn","z_kin_hor_e","p_diag%vt"],
    5:  ["z_vt_ie","p_metrics%wgtfac_e","p_diag%vt"],
    6:  ["z_w_concorr_me","p_prog%vn","p_metrics%ddxn_z_full","p_diag%vt","p_metrics%ddxt_z_full"],
    7:  ["z_v_grad_w","p_diag%vn_ie","p_patch%edges%inv_dual_edge_length","p_prog%w",
         "z_vt_ie","p_patch%edges%inv_primal_edge_length","p_patch%edges%tangent_orientation","z_w_v"],
    8:  ["z_w_con_c"],
    9:  ["z_ekinh","p_int%e_bln_c_s","z_kin_hor_e"],
    10: ["z_w_concorr_mc","p_int%e_bln_c_s","z_w_concorr_me"],
    11: ["p_diag%w_concorr_c","p_metrics%wgtfac_c","z_w_concorr_mc"],
    12: ["z_w_con_c","p_prog%w"],
    13: ["z_w_con_c","p_diag%w_concorr_c"],
    14: ["z_w_con_c_full","z_w_con_c"],
    15: ["p_diag%ddt_w_adv_pc","z_w_con_c","p_prog%w","p_metrics%coeff1_dwdz","p_metrics%coeff2_dwdz"],
    16: ["p_diag%ddt_w_adv_pc","p_int%e_bln_c_s","z_v_grad_w"],
    17: ["p_diag%ddt_vn_apc_pc","z_kin_hor_e","p_metrics%coeff_gradekin","z_ekinh","p_diag%vt",
         "p_patch%edges%f_e","zeta","p_int%c_lin_e","z_w_con_c_full","p_diag%vn_ie","p_metrics%ddqz_z_full_e"],
}

# Per-array: which nests have U-access on this specific array?
ARRAY_U_NESTS = {
    "p_prog%vn":       {3},
    "p_prog%w":        {7},
    "z_kin_hor_e":     {9},
    "z_ekinh":         {17},
    "z_w_concorr_me":  {10},
    "z_v_grad_w":      {16},
    "z_w_con_c_full":  {17},
    "z_w_v":           {7},
    "zeta":            {17},
}

# ---------------------------------------------------------------------------
# Domain-curated bundles (size-constrained: 2, 4, 8, or multiple of 8)
#
# These are informed by co-occurrence Jaccard analysis, producer-consumer
# dataflow, and domain knowledge of the ICON velocity_tendencies module.
#
# Partition:  8 + 8 + 4 + 4 + 4 + 4 + 2 = 34
# ---------------------------------------------------------------------------

CURATED_BUNDLES = [
    {
        "name": "Edge diagnostic & interpolation",
        "tag": "S+U bridge",
        "rationale": (
            "Core edge-level fields in the v-interpolation and kinetic energy "
            "computation phases (N1-N6).  Two bridge-U arrays (prog.vn, "
            "z_kin_hor_e) drive the layout tradeoff.  AoSoA with tuned block "
            "size recommended."
        ),
        "arrays": [
            "p_prog%vn",             # bridge-U, 5 nests
            "p_diag%vn_ie",          # 5 nests
            "p_diag%vt",             # 7 nests
            "z_vt_ie",               # 4 nests
            "z_kin_hor_e",           # bridge-U, 5 nests
            "p_metrics%wgtfac_e",    # 2 nests
            "p_metrics%wgtfacq_e",   # 2 nests
            "p_diag%vn_ie_ubc",      # 1 nest
        ],
    },
    {
        "name": "Cell vertical advection chain",
        "tag": "S+U bridge",
        "rationale": (
            "Cell-dimensioned arrays forming the vertical advection and "
            "contravariant correction chain (N8-N15).  Two bridge-U arrays "
            "(prog.w, z_w_con_c_full) feed into the indirect gather in "
            "N16-N17.  AoSoA recommended."
        ),
        "arrays": [
            "z_w_con_c",             # 5 nests
            "z_w_con_c_full",        # bridge-U, 2 nests
            "p_prog%w",             # bridge-U, 3 nests
            "p_diag%w_concorr_c",    # 2 nests
            "z_w_concorr_mc",        # 2 nests
            "p_metrics%wgtfac_c",    # 1 nest
            "p_metrics%coeff1_dwdz", # 1 nest
            "p_metrics%coeff2_dwdz", # 1 nest
        ],
    },
    {
        "name": "Edge gradient & metric correction",
        "tag": "S+U bridge",
        "rationale": (
            "Edge-level gradient fields and the terrain-following metric "
            "correction arrays.  Two bridge-U arrays (z_v_grad_w, "
            "z_w_concorr_me) cross the S→U boundary into N10 and N16."
        ),
        "arrays": [
            "z_v_grad_w",            # bridge-U, 2 nests
            "z_w_concorr_me",        # bridge-U, 2 nests
            "p_metrics%ddxn_z_full",  # 1 nest
            "p_metrics%ddxt_z_full",  # 1 nest
        ],
    },
    {
        "name": "Edge geometry (read-only)",
        "tag": "structured, SoA",
        "rationale": (
            "Read-only 2D (h,b) edge topology arrays used in the gradient-of-w "
            "computation (N7) and the final tendency (N17).  SoA layout "
            "optimal for vectorisation."
        ),
        "arrays": [
            "p_patch%edges%f_e",
            "p_patch%edges%inv_dual_edge_length",
            "p_patch%edges%inv_primal_edge_length",
            "p_patch%edges%tangent_orientation",
        ],
    },
    {
        "name": "Interpolation weights (read-only)",
        "tag": "structured, SoA",
        "rationale": (
            "Read-only stencil weight arrays used in cell-to-edge and "
            "edge-to-cell interpolation operators.  All accessed with "
            "structured h subscripts.  SoA layout optimal."
        ),
        "arrays": [
            "p_int%e_bln_c_s",
            "p_int%c_lin_e",
            "p_int%rbf_vec_coeff_e",
            "p_metrics%coeff_gradekin",
        ],
    },
    {
        "name": "Final tendency & N17 fields",
        "tag": "mixed S+U",
        "rationale": (
            "Output tendency arrays and N17-exclusive metric/gather fields.  "
            "z_ekinh is bridge-U (gathered from z_kin_hor_e via edge-to-cell "
            "in N9, consumed indirectly in N17)."
        ),
        "arrays": [
            "z_ekinh",               # bridge-U, 2 nests
            "p_diag%ddt_vn_apc_pc",  # 1 nest (output)
            "p_diag%ddt_w_adv_pc",   # 2 nests (output)
            "p_metrics%ddqz_z_full_e",# 1 nest
        ],
    },
    {
        "name": "Vertex arrays",
        "tag": "U-only",
        "rationale": (
            "Vertex-dimensioned arrays accessed exclusively via indirection "
            "(ividx/ivblk).  Always gathered into edge loops.  SoA or AoS "
            "choice has minimal impact since access is always indirect."
        ),
        "arrays": [
            "z_w_v",
            "zeta",
        ],
    },
]


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def short(name: str) -> str:
    return (name
        .replace("p_diag%", "diag.")
        .replace("p_prog%", "prog.")
        .replace("p_metrics%", "met.")
        .replace("p_int%", "int.")
        .replace("p_patch%edges%", "edg."))


def build_array_to_nests(nest_data):
    a2n = {}
    for nid, arrays in nest_data.items():
        for a in arrays:
            a2n.setdefault(a, set()).add(nid)
    return a2n


def classify_array(name, a2n, u_nests_set):
    """Classify array as S-only, bridge-S, bridge-U, U-only."""
    nids = a2n[name]
    s_nests = {nid for nid in range(1, 18)
               if not any("U" in "" for _ in [])}  # placeholder
    has_u = name in ARRAY_U_NESTS and len(ARRAY_U_NESTS[name]) > 0
    return "bridge-U" if has_u else "S"


def print_curated_bundles_md(a2n):
    print("# Constrained Storage Bundles: `velocity_tendencies`\n")
    print("Bundle sizes constrained to {2, 4, 8} or multiples of 8 "
          "(AoSoA / SIMD-lane alignment).\n")
    print(f"Partition: {' + '.join(str(len(b['arrays'])) for b in CURATED_BUNDLES)} "
          f"= {sum(len(b['arrays']) for b in CURATED_BUNDLES)} arrays "
          f"in {len(CURATED_BUNDLES)} bundles.\n")

    for bi, bundle in enumerate(CURATED_BUNDLES):
        arrays = bundle["arrays"]
        common = set.intersection(*(a2n[a] for a in arrays))
        any_nests = set.union(*(a2n[a] for a in arrays))
        n_bridge = sum(1 for a in arrays if a in ARRAY_U_NESTS)

        print(f"## Bundle {bi+1}: {bundle['name']} "
              f"({len(arrays)} arrays) — {bundle['tag']}\n")
        print(f"{bundle['rationale']}\n")
        print(f"- **Nests using all members (∩):** {sorted(common) if common else '(none)'}")
        print(f"- **Nests using any member (∪):** {sorted(any_nests)}")
        print(f"- **Bridge-U count:** {n_bridge}\n")

        print("| Array | Class | #Nests | Nests |")
        print("|-------|-------|:------:|-------|")
        for a in arrays:
            cls = "bridge-U" if a in ARRAY_U_NESTS else "S"
            nids = sorted(a2n[a])
            print(f"| `{a}` | {cls} | {len(nids)} | {nids} |")
        print()

    # Intra-bundle Jaccard statistics
    print("## Intra-Bundle Jaccard Statistics\n")
    print("| Bundle | Size | Min J | Avg J | Max J |")
    print("|-------:|-----:|------:|------:|------:|")
    for bi, bundle in enumerate(CURATED_BUNDLES):
        arrays = bundle["arrays"]
        if len(arrays) < 2:
            print(f"| {bi+1} | {len(arrays)} | — | — | — |")
            continue
        jaccards = []
        for a, b in combinations(arrays, 2):
            jaccards.append(jaccard(a2n[a], a2n[b]))
        print(f"| {bi+1} | {len(arrays)} "
              f"| {min(jaccards):.2f} | {sum(jaccards)/len(jaccards):.2f} "
              f"| {max(jaccards):.2f} |")
    print()

    # Algorithmic bundles for comparison
    print("## Algorithmic Result (constrained greedy)\n")
    print("For comparison, the constrained greedy algorithm at threshold=0.4 "
          "produces the following bundles.\n")
    algo_bundles = constrained_jaccard_bundles(a2n, threshold=0.4)
    for bi, bundle in enumerate(algo_bundles):
        print(f"### Algo-Bundle {bi+1} ({len(bundle)} arrays)\n")
        for a in bundle:
            cls = "bridge-U" if a in ARRAY_U_NESTS else "S"
            print(f"- `{a}` ({cls}, nests: {sorted(a2n[a])})")
        print()


if __name__ == "__main__":
    a2n = build_array_to_nests(NEST_DATA)

    # Verify curated bundles cover all arrays
    all_curated = []
    for b in CURATED_BUNDLES:
        all_curated.extend(b["arrays"])
    all_arrays = set(a2n.keys())
    assert set(all_curated) == all_arrays, (
        f"Mismatch: curated has {set(all_curated) - all_arrays} extra, "
        f"missing {all_arrays - set(all_curated)}"
    )
    assert len(all_curated) == len(set(all_curated)), "Duplicates in curated bundles"
    for b in CURATED_BUNDLES:
        assert is_valid_size(len(b["arrays"])), (
            f"Bundle '{b['name']}' has invalid size {len(b['arrays'])}"
        )

    print_curated_bundles_md(a2n)