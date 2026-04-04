#!/usr/bin/env python3
"""
Array co-access analysis for ICON subroutines via Loki.

Parses a Fortran source file, discovers loop nests, extracts array
working sets per nest, then reports co-access statistics as Markdown.

Install:  pip install "loki @ git+https://github.com/ecmwf-ifs/loki.git"
Usage:    python co_access_analysis.py <fortran_file> [subroutine_name]

If subroutine_name is omitted, defaults to "velocity_tendencies".

Sections produced
-----------------
  1. Per-nest working sets       — what each nest touches
  2. Per-array nest membership   — where each array appears
  3. Pairwise co-occurrence      — how often two arrays share a nest
  4. Co-occurrence frequency     — distribution of shared-nest counts
  5. Equivalence classes         — arrays with identical membership
  6. Nest intersection matrix    — array overlap between nests
"""

import sys
import io
import builtins
from pathlib import Path
from dataclasses import dataclass, field
from collections import defaultdict
from itertools import combinations

from loki import Sourcefile, FindNodes, FindVariables, fgen
from loki.ir import Loop, Assignment
from loki.expression.symbols import Array, Scalar, RangeIndex


# ===========================================================================
#  Loki helpers (from analyze_veltend.py)
# ===========================================================================

VAR_TO_ROLE = {"jc": "h", "je": "h", "jk": "v", "jb": "b"}

INDIRECT_ARRAYS = {
    "icidx", "icblk", "ieidx", "ieblk", "ividx", "ivblk",
    "iqidx", "iqblk", "incidx", "incblk",
}


def _role(varname):
    if not varname:
        return None
    return VAR_TO_ROLE.get(varname.lower(), varname.lower())


def _classify_range(bounds_str, var):
    b = bounds_str.lower()
    role = _role(var)
    if role == "v":
        if "1" in b and "nlev" in b and "+" not in b and "-" not in b.split(":")[-1]:
            return "full_vert"
        return "partial_vert"
    if role == "h":
        return "full_horiz"
    if role == "b":
        return "full_block"
    return "unknown"


# ---------------------------------------------------------------------------
# Expression helpers
# ---------------------------------------------------------------------------

def get_full_name(var):
    if hasattr(var, "parent") and var.parent is not None:
        return get_full_name(var.parent) + "%" + var.basename
    return str(var.name) if hasattr(var, "name") else str(var)


def expr_contains_var(expr, varname):
    vl = varname.lower()
    if isinstance(expr, Scalar):
        return expr.name.lower() == vl
    if isinstance(expr, Array):
        if expr.name.lower() == vl:
            return True
        if expr.dimensions:
            return any(expr_contains_var(d, varname) for d in expr.dimensions)
    if hasattr(expr, "children"):
        return any(expr_contains_var(c, varname) for c in expr.children)
    return False


def is_indirect_subscript(expr):
    if isinstance(expr, Array):
        return expr.name.lower() in INDIRECT_ARRAYS
    if hasattr(expr, "children"):
        return any(is_indirect_subscript(c) for c in expr.children)
    return False


def dominant_loop_var(expr, loop_vars):
    for lv in loop_vars:
        if expr_contains_var(expr, lv):
            return lv
    return None


# ---------------------------------------------------------------------------
# ArrayAccessInfo
# ---------------------------------------------------------------------------

@dataclass
class ArrayAccessInfo:
    name: str
    dim_vars: list
    is_indirect: list

    @property
    def structuredness(self):
        parts = []
        for v, ind in zip(self.dim_vars, self.is_indirect):
            if ind:
                parts.append("U")
            elif v and v not in ("const", "?"):
                parts.append("S")
            else:
                parts.append("C")
        return tuple(parts)

    @property
    def role_su_signature(self):
        parts = []
        for v, su in zip(self.dim_vars, self.structuredness):
            role = _role(v) if v and v not in ("const", "?") else "-"
            parts.append(f"{role}:{su}")
        return "  ".join(parts)


def extract_array_accesses(node, loop_vars):
    accesses = []
    for var in FindVariables().visit(node):
        if not isinstance(var, Array):
            continue
        if var.dimensions is None or len(var.dimensions) == 0:
            continue
        if var.name.lower() in INDIRECT_ARRAYS:
            continue
        name = get_full_name(var)
        dim_vars, is_indirect = [], []
        for dim in var.dimensions:
            if isinstance(dim, RangeIndex):
                dim_vars.append(None)
                is_indirect.append(False)
                continue
            indirect = is_indirect_subscript(dim)
            is_indirect.append(indirect)
            lv = dominant_loop_var(dim, loop_vars)
            dim_vars.append(lv if lv else ("?" if indirect else "const"))
        accesses.append(ArrayAccessInfo(
            name=name, dim_vars=dim_vars, is_indirect=is_indirect))
    return accesses


# ---------------------------------------------------------------------------
# LoopNestInfo
# ---------------------------------------------------------------------------

@dataclass
class LoopNestInfo:
    vars: list
    bounds: list
    depth: int = 1
    loop_shape: str = ""
    ranges: list = field(default_factory=list)
    behavior: str = "compute"
    accesses: list = field(default_factory=list)
    ir_nodes: list = field(default_factory=list)

    def classify(self):
        self.depth = len(self.vars)
        self.loop_shape = ".".join(_role(v) or "?" for v in self.vars)
        self.ranges = [_classify_range(b, v)
                       for b, v in zip(self.bounds, self.vars)]

    @property
    def unique_array_names(self):
        return sorted(set(acc.name for acc in self.accesses))

    @property
    def has_unstructured(self):
        return any("U" in acc.structuredness for acc in self.accesses)

    @property
    def collapsed_su(self):
        role_worst = {}
        priority = {"U": 3, "S": 2, "C": 1}
        for acc in self.accesses:
            for v, su in zip(acc.dim_vars, acc.structuredness):
                role = _role(v) if v and v not in ("const", "?") else None
                if role is None:
                    continue
                cur = role_worst.get(role, "C")
                if priority.get(su, 0) > priority.get(cur, 0):
                    role_worst[role] = su
        order = {"h": 0, "v": 1, "b": 2}
        parts = []
        for role in sorted(role_worst.keys(),
                           key=lambda r: order.get(r, 99)):
            parts.append(f"{role}:{role_worst[role]}")
        return "  ".join(parts) if parts else "-:C"

    @property
    def source_code(self):
        if not self.ir_nodes:
            return "! (no IR nodes stored)"
        if len(self.ir_nodes) == 1:
            return fgen(self.ir_nodes[0])
        outer_role = _role(self.vars[0])
        if outer_role == "b":
            return (f"! inside DO {self.vars[0]} = {self.bounds[0]}\n"
                    + fgen(self.ir_nodes[1]))
        else:
            return fgen(self.ir_nodes[0])


# ---------------------------------------------------------------------------
# Behavior detection
# ---------------------------------------------------------------------------

def detect_behavior(loop_body):
    has_accumulate = False
    for assign in FindNodes(Assignment).visit(loop_body):
        lhs_name = get_full_name(assign.lhs).lower()
        rhs_str = str(assign.rhs).lower()
        if "max(" in rhs_str or "min(" in rhs_str:
            if lhs_name in rhs_str:
                return "reduction"
        if lhs_name in rhs_str:
            has_accumulate = True
    return "accumulate" if has_accumulate else "compute"


# ---------------------------------------------------------------------------
# Loop nest discovery
# ---------------------------------------------------------------------------

def _direct_child_loops(node):
    result = []
    for attr in ("body", "else_body"):
        children = getattr(node, attr, None)
        if children is None:
            continue
        if not isinstance(children, (list, tuple)):
            children = [children]
        for child in children:
            if isinstance(child, Loop):
                result.append(child)
            elif hasattr(child, "body") or hasattr(child, "else_body"):
                result.extend(_direct_child_loops(child))
    return result


def analyze_loop_nests(routine):
    """
    Walk the routine body and discover all leaf loop nests.

    Returns a list of LoopNestInfo objects, each representing a 1- or
    2-level loop nest with its array access inventory.
    """
    nests = []
    all_loops = FindNodes(Loop).visit(routine.body)
    captured_as_inner = set()

    for outer_loop in all_loops:
        outer_var = str(outer_loop.variable)
        direct_children = _direct_child_loops(outer_loop)
        if not direct_children:
            continue

        for inner_loop in direct_children:
            inner_var = str(inner_loop.variable)
            if _direct_child_loops(inner_loop):
                continue
            captured_as_inner.add(id(inner_loop))
            loop_vars = list(dict.fromkeys(
                [outer_var.lower(), inner_var.lower()] +
                (["jb"] if "jb" not in
                 [outer_var.lower(), inner_var.lower()] else [])
            ))
            accesses = extract_array_accesses(inner_loop.body, loop_vars)
            behavior = detect_behavior(inner_loop.body)
            if not accesses and behavior == "compute":
                continue
            nest = LoopNestInfo(
                vars=[outer_var, inner_var],
                bounds=[str(outer_loop.bounds), str(inner_loop.bounds)],
                accesses=accesses, behavior=behavior,
                ir_nodes=[outer_loop, inner_loop],
            )
            nest.classify()
            nests.append(nest)

    for loop in all_loops:
        if id(loop) in captured_as_inner:
            continue
        if _direct_child_loops(loop):
            continue
        var = str(loop.variable)
        role = _role(var)
        if role not in ("h", "v"):
            continue
        loop_vars = [var.lower()]
        if "jb" not in loop_vars:
            loop_vars.append("jb")
        accesses = extract_array_accesses(loop.body, loop_vars)
        behavior = detect_behavior(loop.body)
        if not accesses and behavior == "compute":
            continue
        nest = LoopNestInfo(
            vars=[var], bounds=[str(loop.bounds)],
            accesses=accesses, behavior=behavior,
            ir_nodes=[loop],
        )
        nest.classify()
        nests.append(nest)

    return nests


# ===========================================================================
#  Co-access analysis
# ===========================================================================

def build_array_to_nests(nests):
    """Invert nest -> arrays into array -> set of nest IDs (1-based)."""
    a2n = {}
    for idx, nest in enumerate(nests):
        nid = idx + 1
        for arr in nest.unique_array_names:
            a2n.setdefault(arr, set()).add(nid)
    return a2n


def build_co_occurrence(a2n):
    """Count how many nests each pair of arrays shares."""
    cooccur = {}
    all_arrays = sorted(a2n)
    for a, b in combinations(all_arrays, 2):
        shared = len(a2n[a] & a2n[b])
        if shared > 0:
            cooccur[(a, b)] = shared
    return cooccur


def build_equivalence_classes(a2n):
    """Group arrays with identical nest-membership sets."""
    key_to_arrays = defaultdict(list)
    for arr, nids in a2n.items():
        key_to_arrays[frozenset(nids)].append(arr)
    groups = []
    for nids, arrays in sorted(key_to_arrays.items(),
                                key=lambda x: (-len(x[1]), sorted(x[0]))):
        groups.append((nids, sorted(arrays)))
    return groups


# ===========================================================================
#  Markdown report
# ===========================================================================

def print_report(routine_name, nests):
    a2n = build_array_to_nests(nests)
    cooccur = build_co_occurrence(a2n)
    equiv = build_equivalence_classes(a2n)
    all_arrays = sorted(a2n)
    n_pairs = sum(1 for _ in combinations(all_arrays, 2))
    n_arr_nest = sum(len(nest.unique_array_names) for nest in nests)

    print(f"# Array co-access analysis: `{routine_name}`\n")
    print(f"{len(all_arrays)} unique arrays across {len(nests)} loop "
          f"nests, forming {n_arr_nest} array-nest pairs.\n")

    # ── 1. Per-nest working sets ─────────────────────────────────────
    print("## 1. Per-nest working sets\n")
    print("Each row lists the arrays that a single loop nest accesses "
          "together. Arrays within the same row share a working set: "
          "they are loaded into cache during the same inner-loop "
          "execution and therefore benefit from co-located storage.\n")

    print("| Nest | Shape | Collapsed | Behavior | Size | Arrays |")
    print("|-----:|-------|-----------|----------|-----:|--------|")
    for idx, nest in enumerate(nests):
        nid = idx + 1
        beh = nest.behavior if nest.behavior != "compute" else ""
        tag = "**U**" if nest.has_unstructured else "S"
        arr_list = ", ".join(f"`{a}`" for a in nest.unique_array_names)
        print(f"| {nid} | `{nest.loop_shape}` {tag} "
              f"| `{nest.collapsed_su}` | {beh} "
              f"| {len(nest.unique_array_names)} | {arr_list} |")
    print()

    # ── 2. Per-array nest membership ─────────────────────────────────
    print("## 2. Per-array nest membership\n")
    print("Each row shows which nests access a given array. Arrays "
          "appearing in many nests have the most constrained layout "
          "requirements.\n")

    sorted_arrays = sorted(all_arrays,
                           key=lambda a: (-len(a2n[a]), a))

    print("| Array | #Nests | Nests |")
    print("|-------|:------:|-------|")
    for arr in sorted_arrays:
        nids = sorted(a2n[arr])
        print(f"| `{arr}` | {len(nids)} | {nids} |")
    print()

    # ── 3. Pairwise co-occurrence ────────────────────────────────────
    print("## 3. Pairwise co-occurrence\n")
    print(f"Of the {n_pairs} possible array pairs, "
          f"{len(cooccur)} share at least one nest.\n")

    ranked = sorted(cooccur.items(), key=lambda x: (-x[1], x[0]))

    print("| Array A | Array B | Shared nests | Which nests |")
    print("|---------|---------|:------------:|-------------|")
    for (a, b), count in ranked[:50]:
        shared = sorted(a2n[a] & a2n[b])
        print(f"| `{a}` | `{b}` | {count} | {shared} |")
    print()

    if len(ranked) > 50:
        print(f"({len(ranked) - 50} additional pairs omitted.)\n")

    # ── 4. Co-occurrence frequency distribution ──────────────────────
    print("## 4. Co-occurrence frequency distribution\n")
    print("How many array pairs share exactly k nests.\n")

    freq = defaultdict(int)
    for count in cooccur.values():
        freq[count] += 1
    freq[0] = n_pairs - len(cooccur)

    print("| Shared nests (k) | #Pairs | Fraction |")
    print("|-----------------:|-------:|---------:|")
    for k in sorted(freq):
        frac = freq[k] / n_pairs
        print(f"| {k} | {freq[k]} | {frac:.1%} |")
    print()

    # ── 5. Equivalence classes ───────────────────────────────────────
    print("## 5. Equivalence classes (identical nest membership)\n")
    print("Arrays that appear in exactly the same set of nests. "
          "They are indistinguishable from a co-access perspective.\n")

    multi = [(nids, arrs) for nids, arrs in equiv if len(arrs) > 1]
    single = [(nids, arrs) for nids, arrs in equiv if len(arrs) == 1]

    if multi:
        print(f"### Multi-member classes ({len(multi)} classes, "
              f"{sum(len(a) for _, a in multi)} arrays)\n")
        for nids, arrs in multi:
            nids_sorted = sorted(nids)
            print(f"**Nests {{{', '.join(str(n) for n in nids_sorted)}}}** "
                  f"— {len(arrs)} arrays:\n")
            for a in arrs:
                print(f"- `{a}`")
            print()
    else:
        print("No multi-member equivalence classes found.\n")

    print(f"### Singletons ({len(single)} arrays)\n")
    print("| Array | Nests |")
    print("|-------|-------|")
    for nids, arrs in single:
        print(f"| `{arrs[0]}` | {sorted(nids)} |")
    print()

    # ── 6. Nest intersection matrix ──────────────────────────────────
    print("## 6. Nest intersection sizes\n")
    print("How many arrays two nests share.\n")

    nest_ids = list(range(1, len(nests) + 1))
    print("| | " + " | ".join(f"N{n}" for n in nest_ids) + " |")
    print("|---" + "|---:" * len(nest_ids) + "|")
    for ni in nest_ids:
        si = set(nests[ni - 1].unique_array_names)
        row = f"| **N{ni}** |"
        for nj in nest_ids:
            sj = set(nests[nj - 1].unique_array_names)
            inter = len(si & sj)
            row += f" {inter if inter > 0 else '·'} |"
        print(row)
    print()


# ===========================================================================
#  Jaccard similarity
# ===========================================================================

def jaccard(set_a, set_b):
    """Jaccard similarity: |A ∩ B| / |A ∪ B|."""
    if not set_a and not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)


def compute_jaccard_matrix(a2n):
    """Pairwise Jaccard on nest-membership sets for all array pairs."""
    all_arrays = sorted(a2n)
    jac = {}
    for a, b in combinations(all_arrays, 2):
        jac[(a, b)] = jaccard(a2n[a], a2n[b])
    return jac


# ===========================================================================
#  Greedy seed-and-grow bundling (Jaccard >= threshold)
# ===========================================================================

def greedy_jaccard_bundles(a2n, threshold=0.5):
    """
    Group arrays into bundles where every pair has Jaccard >= threshold.

    Algorithm:
      1. Sort arrays by descending nest count (most-referenced first).
      2. Take the first unassigned array as seed of a new bundle.
      3. Walk remaining unassigned arrays in sorted order.
         A candidate joins iff J(candidate, m) >= threshold for ALL
         current members m.
      4. Close the bundle. Repeat from (2) until all assigned.

    Returns list of bundles (each a list of array names).
    """
    items = sorted(a2n, key=lambda a: -len(a2n[a]))
    assigned = set()
    bundles = []

    for seed in items:
        if seed in assigned:
            continue
        bundle = [seed]
        assigned.add(seed)

        for cand in items:
            if cand in assigned:
                continue
            min_j = min(jaccard(a2n[m], a2n[cand]) for m in bundle)
            if min_j >= threshold:
                bundle.append(cand)
                assigned.add(cand)

        bundles.append(bundle)
    return bundles


# ===========================================================================
#  Size-constrained post-processing
# ===========================================================================

VALID_BUNDLE_SIZES = {2, 4, 8, 16, 24, 32, 40, 48, 56, 64}


def is_valid_size(n):
    """Check if n is a valid bundle size: 2, 4, 8, or a multiple of 8."""
    if n in (2, 4, 8):
        return True
    return n >= 8 and n % 8 == 0


def next_valid_size(n):
    """Smallest valid size >= n."""
    for v in sorted(VALID_BUNDLE_SIZES):
        if v >= n:
            return v
    return ((n + 7) // 8) * 8


def constrain_bundles(bundles, a2n):
    """
    Merge and pad bundles until all have valid sizes {2, 4, 8, 16, ...}.

    Strategy:
      1. Repeatedly find the smallest invalid bundle.
      2. Try merging it with the most Jaccard-similar other bundle
         if the combined size is valid.
      3. Otherwise steal the most-similar member from a donor bundle
         (if the donor remains valid after losing a member).
      4. Last resort: force-merge the two smallest bundles.
    """
    bundles = [list(b) for b in bundles]

    def avg_j(ga, gb):
        total = sum(jaccard(a2n[a], a2n[b]) for a in ga for b in gb)
        return total / (len(ga) * len(gb)) if ga and gb else 0.0

    for _ in range(len(bundles) * 20):
        invalid = [i for i, b in enumerate(bundles)
                   if not is_valid_size(len(b))]
        if not invalid:
            break

        invalid.sort(key=lambda i: len(bundles[i]))
        idx = invalid[0]
        b_small = bundles[idx]

        # Try merging with another bundle → valid combined size
        best_merge, best_j = None, -1.0
        for j, b_other in enumerate(bundles):
            if j == idx:
                continue
            if is_valid_size(len(b_small) + len(b_other)):
                jval = avg_j(b_small, b_other)
                if jval > best_j:
                    best_j = jval
                    best_merge = j

        if best_merge is not None:
            bundles[idx] = b_small + bundles[best_merge]
            del bundles[best_merge]
            continue

        # Try stealing members to reach next valid size
        target = next_valid_size(len(b_small))
        need = target - len(b_small)
        stolen = False
        for _ in range(need):
            best_di, best_m, best_mj = None, None, -1.0
            for j, b_donor in enumerate(bundles):
                if j == idx:
                    continue
                after = len(b_donor) - 1
                if after < 1:
                    continue
                if is_valid_size(len(b_donor)) and not is_valid_size(after):
                    continue
                for m in b_donor:
                    mj = avg_j([m], b_small)
                    if mj > best_mj:
                        best_mj = mj
                        best_di = j
                        best_m = m
            if best_m is not None:
                bundles[idx].append(best_m)
                bundles[best_di].remove(best_m)
                if not bundles[best_di]:
                    del bundles[best_di]
                stolen = True
            else:
                break
        if stolen:
            continue

        # Last resort: force-merge two smallest
        if len(bundles) >= 2:
            sizes = sorted(range(len(bundles)),
                           key=lambda i: len(bundles[i]))
            i1, i2 = sizes[0], sizes[1]
            if i1 > i2:
                i1, i2 = i2, i1
            bundles[i1] = bundles[i1] + bundles[i2]
            del bundles[i2]
        else:
            break

    return bundles


# ===========================================================================
#  Extended report: Jaccard + bundles (sections 7–10)
# ===========================================================================

def print_jaccard_and_bundles(nests, a2n, threshold=0.5):
    """Append sections 7–10 to the report: Jaccard, bundles, constrained."""
    all_arrays = sorted(a2n)
    jac_matrix = compute_jaccard_matrix(a2n)
    n_pairs = sum(1 for _ in combinations(all_arrays, 2))

    # ── 7. Jaccard similarity (top pairs) ────────────────────────────
    print("## 7. Jaccard similarity (top pairs)\n")
    print("Jaccard index J(A,B) = |nests(A) ∩ nests(B)| / "
          "|nests(A) ∪ nests(B)| normalises co-occurrence by the "
          "combined footprint of both arrays. J = 1 means identical "
          "nest membership; J = 0 means disjoint.\n")

    ranked = sorted(jac_matrix.items(), key=lambda x: (-x[1], x[0]))

    n_above = sum(1 for _, j in ranked if j >= threshold)
    n_zero = sum(1 for _, j in ranked if j == 0.0)
    print(f"Of {n_pairs} pairs: {n_above} have J >= {threshold}, "
          f"{n_zero} have J = 0 (disjoint).\n")

    print("| Array A | Array B | Jaccard | Shared | |A| | |B| |")
    print("|---------|---------|:-------:|:------:|----:|----:|")
    shown = 0
    for (a, b), j in ranked:
        if j == 0.0:
            break
        shared = len(a2n[a] & a2n[b])
        print(f"| `{a}` | `{b}` | {j:.2f} | {shared} "
              f"| {len(a2n[a])} | {len(a2n[b])} |")
        shown += 1
        if shown >= 50:
            break
    print()

    # Jaccard distribution
    print("### Jaccard distribution\n")
    buckets = defaultdict(int)
    for _, j in jac_matrix.items():
        if j == 0.0:
            buckets["0.00"] += 1
        elif j < 0.25:
            buckets["0.01–0.24"] += 1
        elif j < 0.50:
            buckets["0.25–0.49"] += 1
        elif j < 0.75:
            buckets["0.50–0.74"] += 1
        elif j < 1.00:
            buckets["0.75–0.99"] += 1
        else:
            buckets["1.00"] += 1

    print("| Jaccard range | #Pairs | Fraction |")
    print("|:-------------:|-------:|---------:|")
    for label in ["0.00", "0.01–0.24", "0.25–0.49",
                  "0.50–0.74", "0.75–0.99", "1.00"]:
        cnt = buckets.get(label, 0)
        frac = cnt / n_pairs if n_pairs > 0 else 0
        print(f"| {label} | {cnt} | {frac:.1%} |")
    print()

    # ── 8. Greedy Jaccard bundles (unconstrained) ────────────────────
    print(f"## 8. Greedy Jaccard bundles (threshold = {threshold})\n")
    print("Arrays grouped by the seed-and-grow algorithm: every pair "
          f"within a bundle has J >= {threshold}. Seed order is by "
          "descending nest count.\n")

    raw_bundles = greedy_jaccard_bundles(a2n, threshold)

    sizes_raw = [len(b) for b in raw_bundles]
    print(f"{len(raw_bundles)} bundles: sizes {sizes_raw}, "
          f"total {sum(sizes_raw)} arrays.\n")

    for bi, bundle in enumerate(raw_bundles):
        common = set.intersection(*(a2n[a] for a in bundle))
        any_nests = set.union(*(a2n[a] for a in bundle))

        # Intra-bundle Jaccard stats
        if len(bundle) >= 2:
            intra = [jaccard(a2n[a], a2n[b])
                     for a, b in combinations(bundle, 2)]
            jmin, javg, jmax = min(intra), sum(intra)/len(intra), max(intra)
            jstats = f"J: min={jmin:.2f} avg={javg:.2f} max={jmax:.2f}"
        else:
            jstats = "J: —"

        print(f"### Bundle {bi+1} ({len(bundle)} arrays) — {jstats}\n")
        print(f"- Nests ∩: "
              f"{sorted(common) if common else '(none)'}")
        print(f"- Nests ∪: {sorted(any_nests)}\n")
        print("| Array | #Nests | Nests |")
        print("|-------|:------:|-------|")
        for a in bundle:
            print(f"| `{a}` | {len(a2n[a])} | {sorted(a2n[a])} |")
        print()

    # ── 9. Size-constrained bundles ──────────────────────────────────
    print("## 9. Size-constrained bundles "
          "(valid sizes: 2, 4, 8, multiples of 8)\n")
    print("The unconstrained bundles from §8 are post-processed: "
          "invalid-sized bundles are merged (preferring highest average "
          "Jaccard) or padded (stealing the most-similar member from a "
          "donor) until every bundle reaches a valid size.\n")

    constrained = constrain_bundles(raw_bundles, a2n)

    sizes_con = [len(b) for b in constrained]
    all_valid = all(is_valid_size(s) for s in sizes_con)
    print(f"{len(constrained)} bundles: sizes {sizes_con}, "
          f"total {sum(sizes_con)} arrays. "
          f"{'All sizes valid.' if all_valid else '⚠ Some sizes invalid.'}\n")

    for bi, bundle in enumerate(constrained):
        common = set.intersection(*(a2n[a] for a in bundle))
        any_nests = set.union(*(a2n[a] for a in bundle))

        if len(bundle) >= 2:
            intra = [jaccard(a2n[a], a2n[b])
                     for a, b in combinations(bundle, 2)]
            jmin, javg, jmax = min(intra), sum(intra)/len(intra), max(intra)
            jstats = f"J: min={jmin:.2f} avg={javg:.2f} max={jmax:.2f}"
        else:
            jstats = "J: —"

        print(f"### Bundle {bi+1} ({len(bundle)} arrays) — {jstats}\n")
        print(f"- Nests ∩: "
              f"{sorted(common) if common else '(none)'}")
        print(f"- Nests ∪: {sorted(any_nests)}\n")
        print("| Array | #Nests | Nests |")
        print("|-------|:------:|-------|")
        for a in bundle:
            print(f"| `{a}` | {len(a2n[a])} | {sorted(a2n[a])} |")
        print()

    # ── 10. Bundle comparison summary ────────────────────────────────
    print("## 10. Bundle comparison summary\n")
    print("| | Unconstrained | Constrained |")
    print("|---|---:|---:|")
    print(f"| Bundles | {len(raw_bundles)} | {len(constrained)} |")
    print(f"| Sizes | {sizes_raw} | {sizes_con} |")
    print(f"| Singletons | "
          f"{sum(1 for b in raw_bundles if len(b)==1)} | "
          f"{sum(1 for b in constrained if len(b)==1)} |")
    print(f"| Largest | {max(sizes_raw)} | {max(sizes_con)} |")
    n_valid_raw = sum(1 for s in sizes_raw if is_valid_size(s))
    print(f"| Valid-sized | {n_valid_raw}/{len(raw_bundles)} "
          f"| {sum(1 for s in sizes_con if is_valid_size(s))}"
          f"/{len(constrained)} |")
    print()


# ===========================================================================
#  Entry point
# ===========================================================================

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <fortran_file> [subroutine_name]")
        sys.exit(1)

    filepath = Path(sys.argv[1])
    routine_name = (sys.argv[2] if len(sys.argv) > 2
                    else "velocity_tendencies")

    print(f"Parsing {filepath} ...", file=sys.stderr)
    src = Sourcefile.from_file(
        str(filepath), preprocess=True, defines=["_OPENACC"],
        includes=[str(filepath.parent / "include"),
                  str(filepath.parent)],
    )

    routine = None
    for r in src.all_subroutines:
        if r.name.lower() == routine_name.lower():
            routine = r
            break

    if routine is None:
        available = [r.name for r in src.all_subroutines]
        print(f"ERROR: '{routine_name}' not found. "
              f"Available: {available}", file=sys.stderr)
        sys.exit(1)

    n_loops = len(FindNodes(Loop).visit(routine.body))
    print(f"Found: {routine.name}  ({len(routine.variables)} vars, "
          f"{n_loops} loops)", file=sys.stderr)

    # ── Discover loop nests via Loki ──
    nests = analyze_loop_nests(routine)
    print(f"Discovered {len(nests)} leaf loop nests.",
          file=sys.stderr)

    # ── Generate Markdown report ──
    buf = io.StringIO()
    old_print = builtins.print

    def md_print(*args, **kwargs):
        kwargs["file"] = buf
        old_print(*args, **kwargs)

    builtins.print = md_print
    print_report(routine.name, nests)
    a2n = build_array_to_nests(nests)
    print_jaccard_and_bundles(nests, a2n, threshold=0.5)
    builtins.print = old_print

    report = buf.getvalue()

    # Write to file alongside the source
    outpath = filepath.with_suffix(".co_access.md")
    outpath.write_text(report)
    old_print(f"Report written to {outpath}", file=sys.stderr)

    # Also print to stdout
    old_print(report)


if __name__ == "__main__":
    main()