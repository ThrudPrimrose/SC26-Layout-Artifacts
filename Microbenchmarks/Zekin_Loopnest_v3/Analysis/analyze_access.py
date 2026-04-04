#!/usr/bin/env python3
"""
Analyze loop nests and array access patterns in ICON's velocity_tendencies
using Loki (ECMWF source-to-source framework).

Install:  pip install "loki @ git+https://github.com/ecmwf-ifs/loki.git"
Usage:    python analyze_veltend.py <fortran_file> [subroutine_name]

Terminology
-----------
  h  = horizontal (nproma) dimension.  Loop vars: jc, je  (treated as equivalent)
  v  = vertical (nlev) dimension.      Loop var:  jk
  b  = block dimension.                Loop var:  jb

  S  = structured access:   subscript is a direct loop variable (+/- offset ok)
  U  = unstructured access: subscript wraps an indirection array (icidx, ividx, ...)
  C  = constant access:     subscript is loop-invariant (literal, parameter)

Loop classification
-------------------
  shape:      b.h, v.h, h (single), v (single), b.v, etc.
  range:      full_vert (1:nlev), partial_vert (2:nlev, ...), full_horiz, full_block
  behavior:   compute, accumulate (LHS += ...), reduction (MAX/MIN over loop)
"""

import sys
from pathlib import Path
from dataclasses import dataclass, field
from collections import defaultdict
from itertools import combinations

from loki import Sourcefile, FindNodes, FindVariables, fgen
from loki.ir import Loop, Assignment, Conditional, Section
from loki.expression.symbols import Array, Scalar, RangeIndex

# Try to import InlineCall for reduction detection
try:
    from loki.expression.symbols import InlineCall
    HAS_INLINE_CALL = True
except ImportError:
    HAS_INLINE_CALL = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VAR_TO_ROLE = {"jc": "h", "je": "h", "jk": "v", "jb": "b"}

INDIRECT_ARRAYS = {
    "icidx", "icblk", "ieidx", "ieblk", "ividx", "ivblk",
    "iqidx", "iqblk", "incidx", "incblk",
}

def _role(varname):
    if not varname:
        return None
    return VAR_TO_ROLE.get(varname.lower(), varname.lower())


# ---------------------------------------------------------------------------
# ArrayAccessInfo
# ---------------------------------------------------------------------------

@dataclass
class ArrayAccessInfo:
    """Per-reference: name, per-dim loop var, per-dim structuredness."""
    name: str
    dim_vars: list        # detected loop var per dim (or "const"/"?")
    is_indirect: list     # per-dim: True if uses indirection array

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
        """Canonical form: 'h:S v:S b:S' (je==jc -> h)."""
        parts = []
        for v, su in zip(self.dim_vars, self.structuredness):
            role = _role(v) if v and v not in ("const", "?") else "-"
            parts.append(f"{role}:{su}")
        return "  ".join(parts)


# ---------------------------------------------------------------------------
# LoopNestInfo — now supports 1-level and 2-level nests
# ---------------------------------------------------------------------------

@dataclass
class LoopNestInfo:
    """A loop nest (1 or 2 levels) with access inventory and behavior tag."""
    vars: list                    # loop variable names, outer to inner
    bounds: list                  # bounds strings, outer to inner
    depth: int = 1
    loop_shape: str = ""          # e.g. "b.h", "v.h", "h", "v"
    ranges: list = field(default_factory=list)   # e.g. ["full_block", "full_horiz"]
    behavior: str = "compute"     # "compute", "accumulate", "reduction"
    accesses: list = field(default_factory=list)
    ir_nodes: list = field(default_factory=list)  # Loki Loop IR nodes [outer, inner] or [single]

    def classify(self):
        self.depth = len(self.vars)
        self.loop_shape = ".".join(_role(v) or "?" for v in self.vars)
        self.ranges = [_classify_range(b, v) for b, v in zip(self.bounds, self.vars)]

    @property
    def unique_array_names(self):
        return sorted(set(acc.name for acc in self.accesses))

    @property
    def has_unstructured(self):
        return any("U" in acc.structuredness for acc in self.accesses)

    # --- Per-array (fine-grained) ---

    @property
    def role_su_patterns(self):
        return sorted(set(acc.role_su_signature for acc in self.accesses))

    @property
    def has_mixed_su(self):
        sigs = set(acc.structuredness for acc in self.accesses)
        has_pure = any(all(c in ("S", "C") for c in s) for s in sigs)
        has_u = any("U" in s for s in sigs)
        return has_pure and has_u

    # --- Collapsed (worst-case per role) ---

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
        for role in sorted(role_worst.keys(), key=lambda r: order.get(r, 99)):
            parts.append(f"{role}:{role_worst[role]}")
        return "  ".join(parts) if parts else "-:C"

    @property
    def pattern_key(self):
        return (self.loop_shape, tuple(self.ranges), self.behavior, self.collapsed_su)

    @property
    def shape_key(self):
        return (self.loop_shape, tuple(self.ranges), self.behavior)

    @property
    def source_code(self):
        if not self.ir_nodes:
            return "! (no IR nodes stored)"
        if len(self.ir_nodes) == 1:
            return fgen(self.ir_nodes[0])
        outer_role = _role(self.vars[0])
        if outer_role == "b":
            return f"! inside DO {self.vars[0]} = {self.bounds[0]}\n" + fgen(self.ir_nodes[1])
        else:
            return fgen(self.ir_nodes[0])


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
        accesses.append(ArrayAccessInfo(name=name, dim_vars=dim_vars, is_indirect=is_indirect))
    return accesses


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
                (["jb"] if "jb" not in [outer_var.lower(), inner_var.lower()] else [])
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


# ---------------------------------------------------------------------------
# Grouping (original)
# ---------------------------------------------------------------------------

def group_arrays(nests):
    array_role_sigs = defaultdict(set)
    for nest in nests:
        for acc in nest.accesses:
            array_role_sigs[acc.name].add(acc.role_su_signature)
    role_groups = defaultdict(set)
    for name, sigs in array_role_sigs.items():
        role_groups[frozenset(sigs)].add(name)
    return array_role_sigs, role_groups


# ===========================================================================
#  SET ANALYSIS — new section
# ===========================================================================

def _jaccard(set_a, set_b):
    """Jaccard similarity between two sets."""
    if not set_a and not set_b:
        return 0.0
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union > 0 else 0.0


def compute_set_analysis(nests):
    """
    Build all the derived data structures for the set analysis report.

    Returns a dict with:
      all_arrays        – sorted list of every unique array name
      array_to_nests    – {name: set of nest indices (1-based)}
      nest_to_arrays    – {nest_idx: set of array names}
      cooccurrence      – {(a,b): int count}  (symmetric, a<=b)
      jaccard_arrays    – {(a,b): float}       (symmetric, a<=b)
      array_class       – {name: "S-only"|"U-only"|"bridge-U"|"bridge-S"}
      nest_jaccard      – [(i,j,jac,shared)] sorted desc by jac
      bundles           – list of lists (greedy Jaccard>=0.5 grouping)
      producer_consumer – [(producer_nest, consumer_nest, shared_arrays)]
    """
    n = len(nests)
    all_arrays = sorted(set(a for nest in nests for a in nest.unique_array_names))

    # --- Array <-> Nest membership ---
    array_to_nests = {a: set() for a in all_arrays}
    nest_to_arrays = {}
    for idx, nest in enumerate(nests):
        nid = idx + 1  # 1-based
        nest_to_arrays[nid] = set(nest.unique_array_names)
        for a in nest.unique_array_names:
            array_to_nests[a].add(nid)

    # --- Per-array: which nests have U-access on this specific array? ---
    array_u_nests = {a: set() for a in all_arrays}
    for idx, nest in enumerate(nests):
        nid = idx + 1
        for acc in nest.accesses:
            if "U" in acc.structuredness:
                array_u_nests[acc.name].add(nid)

    # --- S-nest / U-nest classification ---
    s_nests = set()
    u_nests = set()
    for idx, nest in enumerate(nests):
        nid = idx + 1
        if nest.has_unstructured:
            u_nests.add(nid)
        else:
            s_nests.add(nid)

    # --- Array classification ---
    array_class = {}
    for a in all_arrays:
        nids = array_to_nests[a]
        in_s = bool(nids & s_nests)
        in_u = bool(nids & u_nests)
        has_direct_u = len(array_u_nests[a]) > 0
        if has_direct_u:
            array_class[a] = "bridge-U"   # this array itself is U-accessed somewhere
        elif in_s and in_u:
            array_class[a] = "bridge-S"   # appears in U-nest but only with S access
        elif in_u:
            array_class[a] = "U-only"
        else:
            array_class[a] = "S-only"

    # --- Co-occurrence (pair count) ---
    cooccurrence = defaultdict(int)
    for nest in nests:
        names = nest.unique_array_names
        for i in range(len(names)):
            for j in range(i, len(names)):
                pair = tuple(sorted([names[i], names[j]]))
                cooccurrence[pair] += 1

    # --- Jaccard similarity between arrays (on nest membership) ---
    jaccard_arrays = {}
    for a, b in combinations(all_arrays, 2):
        pair = tuple(sorted([a, b]))
        jaccard_arrays[pair] = _jaccard(array_to_nests[a], array_to_nests[b])

    # --- Nest-pair overlap (Jaccard on array sets) ---
    nest_jaccard = []
    for i in range(n):
        for j in range(i + 1, n):
            ni, nj = i + 1, j + 1
            sA = nest_to_arrays[ni]
            sB = nest_to_arrays[nj]
            shared = sorted(sA & sB)
            jac = _jaccard(sA, sB)
            if shared:
                nest_jaccard.append((ni, nj, jac, shared))
    nest_jaccard.sort(key=lambda x: -x[2])

    # --- Storage bundles (greedy Jaccard >= 0.5) ---
    bundles = []
    assigned = set()
    # Seed order: most-referenced arrays first
    seeds = sorted(all_arrays, key=lambda a: -len(array_to_nests[a]))
    for seed in seeds:
        if seed in assigned:
            continue
        bundle = [seed]
        assigned.add(seed)
        for cand in seeds:
            if cand in assigned:
                continue
            fits = all(
                _jaccard(array_to_nests[m], array_to_nests[cand]) >= 0.5
                for m in bundle
            )
            if fits:
                bundle.append(cand)
                assigned.add(cand)
        bundles.append(bundle)

    # --- Producer-consumer chains ---
    # A nest "produces" the arrays on LHS; a later nest that reads them is a consumer.
    nest_lhs = {}
    nest_rhs = {}
    for idx, nest in enumerate(nests):
        nid = idx + 1
        lhs_names = set()
        rhs_names = set()
        for assign in FindNodes(Assignment).visit(
            nest.ir_nodes[-1].body if nest.ir_nodes else []
        ):
            lhs_names.add(get_full_name(assign.lhs))
            for var in FindVariables().visit(assign.rhs):
                if isinstance(var, Array) and var.name.lower() not in INDIRECT_ARRAYS:
                    rhs_names.add(get_full_name(var))
        nest_lhs[nid] = lhs_names
        nest_rhs[nid] = rhs_names

    producer_consumer = []
    for i in range(n):
        for j in range(i + 1, n):
            ni, nj = i + 1, j + 1
            # ni produces -> nj consumes
            flow_ij = sorted(nest_lhs[ni] & nest_rhs[nj])
            if flow_ij:
                producer_consumer.append((ni, nj, flow_ij))
            # nj produces -> ni consumes
            flow_ji = sorted(nest_lhs[nj] & nest_rhs[ni])
            if flow_ji:
                producer_consumer.append((nj, ni, flow_ji))

    return {
        "all_arrays": all_arrays,
        "array_to_nests": array_to_nests,
        "array_u_nests": array_u_nests,
        "nest_to_arrays": nest_to_arrays,
        "s_nests": s_nests,
        "u_nests": u_nests,
        "array_class": array_class,
        "cooccurrence": cooccurrence,
        "jaccard_arrays": jaccard_arrays,
        "nest_jaccard": nest_jaccard,
        "bundles": bundles,
        "producer_consumer": producer_consumer,
        "nest_lhs": nest_lhs,
        "nest_rhs": nest_rhs,
    }


# ===========================================================================
#  SET ANALYSIS — Markdown output
# ===========================================================================

def print_set_analysis(routine_name, nests, sa):
    """
    Emit the full set-analysis report as Markdown.

    sa: dict returned by compute_set_analysis()
    """
    p = print

    p(f"\n---\n")
    p(f"# Set Analysis: `{routine_name}`\n")
    p(f"{len(sa['all_arrays'])} unique arrays across {len(nests)} loop nests.\n")

    # ── 1. Array classification ──────────────────────────────────────
    p("## 1. Array Classification\n")
    p("Each array is classified by its participation across structured (S) and unstructured (U) nests.\n")
    p("| Class | Meaning |")
    p("|-------|---------|")
    p("| `S-only` | Appears exclusively in nests with no indirect access |")
    p("| `U-only` | Appears exclusively in nests that have indirect access (rare) |")
    p("| `bridge-U` | The array itself is accessed via indirection in ≥1 nest **and** appears in S-nests — layout conflict |")
    p("| `bridge-S` | Accessed only with structured subscripts, but co-occurs in a U-nest alongside indirect arrays |")
    p()

    class_order = ["S-only", "bridge-S", "bridge-U", "U-only"]
    groups = defaultdict(list)
    for a in sa["all_arrays"]:
        groups[sa["array_class"][a]].append(a)

    for cls in class_order:
        if cls not in groups:
            continue
        arr_list = groups[cls]
        p(f"### {cls} ({len(arr_list)} arrays)\n")
        p("| Array | Nests (count) | U-access nests |")
        p("|-------|---------------|----------------|")
        for a in arr_list:
            nids = sorted(sa["array_to_nests"][a])
            u_nids = sorted(sa["array_u_nests"][a])
            nids_str = ", ".join(str(x) for x in nids)
            u_str = ", ".join(str(x) for x in u_nids) if u_nids else "—"
            p(f"| `{a}` | {nids_str} ({len(nids)}) | {u_str} |")
        p()

    p(f"**Summary**: "
      f"{len(groups.get('S-only',[]))} S-only, "
      f"{len(groups.get('bridge-S',[]))} bridge-S, "
      f"{len(groups.get('bridge-U',[]))} bridge-U, "
      f"{len(groups.get('U-only',[]))} U-only.\n")

    # ── 2. Per-nest working sets ─────────────────────────────────────
    p("## 2. Per-Nest Working Sets\n")
    p("| Nest | Shape | Collapsed | Behavior | #Arrays | Arrays |")
    p("|-----:|-------|-----------|----------|--------:|--------|")
    for idx, nest in enumerate(nests):
        nid = idx + 1
        beh = nest.behavior if nest.behavior != "compute" else ""
        tag = "**U**" if nest.has_unstructured else "S"
        arrays_short = ", ".join(f"`{a}`" for a in nest.unique_array_names)
        p(f"| {nid} | `{nest.loop_shape}` {tag} | `{nest.collapsed_su}` | {beh} "
          f"| {len(nest.unique_array_names)} | {arrays_short} |")
    p()

    # ── 3. Co-occurrence matrix (top pairs) ──────────────────────────
    p("## 3. Array Co-occurrence (top pairs)\n")
    p("Pairs of arrays that appear together in the most nests. High co-occurrence suggests "
      "they benefit from co-located storage.\n")

    sorted_pairs = sorted(
        ((pair, cnt) for pair, cnt in sa["cooccurrence"].items() if pair[0] != pair[1]),
        key=lambda x: -x[1]
    )

    p("| Array A | Array B | Co-occurrence | Jaccard |")
    p("|---------|---------|:------------:|:-------:|")
    seen = set()
    shown = 0
    for (a, b), cnt in sorted_pairs:
        canon = tuple(sorted([a, b]))
        if canon in seen:
            continue
        seen.add(canon)
        jac = sa["jaccard_arrays"].get(canon, 0.0)
        cls_a = sa["array_class"][a]
        cls_b = sa["array_class"][b]
        marker = ""
        if "bridge" in cls_a or "bridge" in cls_b:
            marker = " ⚠"
        p(f"| `{a}` | `{b}` | {cnt} | {jac:.2f}{marker} |")
        shown += 1
        if shown >= 40:
            break
    p()
    p("> ⚠ = at least one array is a bridge (S+U conflict)\n")

    # ── 4. Jaccard clusters (arrays with Jaccard == 1.0) ────────────
    p("## 4. Identical-Membership Clusters (Jaccard = 1.0)\n")
    p("Arrays that appear in exactly the same set of nests — they are inseparable from "
      "a co-access perspective.\n")

    # Build equivalence classes
    membership_key = {}
    for a in sa["all_arrays"]:
        key = frozenset(sa["array_to_nests"][a])
        membership_key[a] = key

    equiv_classes = defaultdict(list)
    for a in sa["all_arrays"]:
        equiv_classes[membership_key[a]].append(a)

    # Only show groups with >1 member
    multi_groups = [(nids, arrs) for nids, arrs in equiv_classes.items() if len(arrs) > 1]
    multi_groups.sort(key=lambda x: -len(x[1]))

    if multi_groups:
        for nids, arrs in multi_groups:
            nids_sorted = sorted(nids)
            has_u = any(sa["array_class"][a] in ("bridge-U", "U-only") for a in arrs)
            tag = "⚠ HAS BRIDGE-U" if has_u else "fully structured"
            p(f"### Cluster: nests {{{', '.join(str(x) for x in nids_sorted)}}} — {tag}\n")
            for a in sorted(arrs):
                p(f"- `{a}` ({sa['array_class'][a]})")
            p()
    else:
        p("No arrays share identical nest membership.\n")

    # ── 5. Nest overlap ──────────────────────────────────────────────
    p("## 5. Nest Overlap (Jaccard on array working sets)\n")
    p("Nest pairs ranked by array-set similarity. High Jaccard ≈ candidates for loop fusion "
      "or shared tiling.\n")

    p("| Nest A | Nest B | Jaccard | |A| | |B| | |A∩B| | Shared arrays |")
    p("|-------:|-------:|:-------:|----:|----:|------:|---------------|")
    for (ni, nj, jac, shared) in sa["nest_jaccard"][:30]:
        nest_i = nests[ni - 1]
        nest_j = nests[nj - 1]
        tag_i = "U" if nest_i.has_unstructured else "S"
        tag_j = "U" if nest_j.has_unstructured else "S"
        sz_i = len(nest_i.unique_array_names)
        sz_j = len(nest_j.unique_array_names)
        shared_str = ", ".join(f"`{a}`" for a in shared)
        p(f"| {ni} ({nest_i.loop_shape} {tag_i}) | {nj} ({nest_j.loop_shape} {tag_j}) "
          f"| {jac:.2f} | {sz_i} | {sz_j} | {len(shared)} | {shared_str} |")
    p()

    # ── 6. Storage bundles ───────────────────────────────────────────
    p("## 6. Storage Bundles (greedy Jaccard ≥ 0.5)\n")
    p("Arrays grouped by nest-membership similarity. Within a bundle, every pair has "
      "Jaccard ≥ 0.5 — they are accessed together often enough to benefit from co-location "
      "in a single AoS or AoSoA record.\n")

    for bi, bundle in enumerate(sa["bundles"]):
        has_u = any(sa["array_class"][a] in ("bridge-U", "U-only") for a in bundle)
        tag = " ⚠ HAS BRIDGE-U" if has_u else ""

        # Which nests use ALL members?
        all_nids = set(range(1, len(nests) + 1))
        common_nests = all_nids
        for a in bundle:
            common_nests = common_nests & sa["array_to_nests"][a]
        common_nests = sorted(common_nests)

        # Which nests use ANY member?
        any_nests = sorted(set().union(*(sa["array_to_nests"][a] for a in bundle)))

        p(f"### Bundle {bi+1}: {len(bundle)} arrays{tag}\n")
        p(f"- **Nests using all members (∩):** "
          f"{', '.join(str(x) for x in common_nests) if common_nests else '(none)'}")
        p(f"- **Nests using any member (∪):** "
          f"{', '.join(str(x) for x in any_nests)}\n")

        p("| Array | Class | Nest count |")
        p("|-------|-------|:----------:|")
        for a in bundle:
            p(f"| `{a}` | {sa['array_class'][a]} | {len(sa['array_to_nests'][a])} |")
        p()

    # ── 7. Producer → Consumer chains ────────────────────────────────
    p("## 7. Producer → Consumer Dataflow\n")
    p("Array-level data dependencies between nests: nest A writes an array that nest B reads. "
      "Ordered by nest index (execution order).\n")

    p("| Producer | Consumer | Dataflow arrays |")
    p("|:--------:|:--------:|-----------------|")
    # Deduplicate and sort
    seen_flows = set()
    for (prod, cons, arrs) in sorted(sa["producer_consumer"], key=lambda x: (x[0], x[1])):
        key = (prod, cons, tuple(arrs))
        if key in seen_flows:
            continue
        seen_flows.add(key)
        nest_p = nests[prod - 1]
        nest_c = nests[cons - 1]
        tag_p = "U" if nest_p.has_unstructured else "S"
        tag_c = "U" if nest_c.has_unstructured else "S"
        arrs_str = ", ".join(f"`{a}`" for a in arrs)
        p(f"| N{prod} ({nest_p.loop_shape} {tag_p}) | N{cons} ({nest_c.loop_shape} {tag_c}) | {arrs_str} |")
    p()

    # ── 8. Bridge-U focus ────────────────────────────────────────────
    p("## 8. Bridge-U Deep Dive\n")
    p("These arrays are the core layout-conflict set: they are accessed with direct "
      "structured subscripts in some nests and via indirection in others.\n")

    bridge_u = [a for a in sa["all_arrays"] if sa["array_class"][a] == "bridge-U"]
    for a in bridge_u:
        s_nest_ids = sorted(sa["array_to_nests"][a] - sa["array_u_nests"][a])
        u_nest_ids = sorted(sa["array_u_nests"][a])
        p(f"### `{a}`\n")
        p(f"- **Structured nests:** {', '.join(str(x) for x in s_nest_ids)} "
          f"({len(s_nest_ids)} nests)")
        p(f"- **Unstructured nests:** {', '.join(str(x) for x in u_nest_ids)} "
          f"({len(u_nest_ids)} nests)")

        # Co-accessed arrays in U-nests
        co_in_u = set()
        for nid in u_nest_ids:
            co_in_u |= sa["nest_to_arrays"][nid]
        co_in_u.discard(a)
        p(f"- **Co-accessed in U-nests:** {', '.join(f'`{x}`' for x in sorted(co_in_u))}")

        # Co-accessed arrays in S-nests
        co_in_s = set()
        for nid in s_nest_ids:
            co_in_s |= sa["nest_to_arrays"][nid]
        co_in_s.discard(a)
        p(f"- **Co-accessed in S-nests:** {', '.join(f'`{x}`' for x in sorted(co_in_s))}")

        # Overlap
        overlap = sorted(co_in_s & co_in_u)
        if overlap:
            p(f"- **Shared co-access (S∩U):** {', '.join(f'`{x}`' for x in overlap)}")
        p()


# ---------------------------------------------------------------------------
# Pretty printing (original report)
# ---------------------------------------------------------------------------

def print_analysis(routine_name, nests, array_role_sigs, role_groups):
    p = print

    p(f"# Loop and Array Access Analysis: `{routine_name}`\n")
    p("| Symbol | Meaning |")
    p("|--------|---------|")
    p("| `h` | horizontal (nproma) dimension — loop vars `jc`, `je` (treated as equivalent) |")
    p("| `v` | vertical (nlev) dimension — loop var `jk` |")
    p("| `b` | block dimension — loop var `jb` |")
    p("| `S` | **structured** access — subscript is a direct loop variable |")
    p("| `U` | **unstructured** access — subscript wraps an indirection array (`icidx`, `ividx`, ...) |")
    p("| `C` | **constant** access — subscript is loop-invariant |")
    p()

    # --- Array groups ---
    p("## Array Groups\n")
    p("Arrays in the same group have identical access patterns across all loops.\n")

    for sigs, names in sorted(role_groups.items(), key=lambda x: -len(x[1])):
        sig_list = sorted(sigs)
        has_u = any("U" in s for s in sig_list)
        tag = "HAS UNSTRUCTURED" if has_u else "FULLY STRUCTURED"
        p(f"### {tag} ({len(names)} arrays)\n")
        p("Patterns observed:\n")
        for sig in sig_list:
            p(f"- `{sig}`")
        p()
        p("Arrays:\n")
        for n in sorted(names):
            p(f"- `{n}`")
        p()

    # --- Per-array summary ---
    p("## Per-Array Access Summary\n")
    p("> `S+U` marks arrays accessed both structured and unstructured (layout tradeoff)\n")
    p("| Array | Role:S/U Patterns | Conflict |")
    p("|-------|-------------------|----------|")
    for name in sorted(array_role_sigs.keys()):
        role_sigs = sorted(array_role_sigs[name])
        has_u = any("U" in s for s in role_sigs)
        has_s_only = any(all(c in ("S", "C", " ", ":", "-") for c in s) for s in role_sigs)
        marker = "S+U" if (has_u and has_s_only) else ""
        sigs_str = " \\| ".join(f"`{s}`" for s in role_sigs)
        p(f"| `{name}` | {sigs_str} | {marker} |")
    p()

    # --- Loop shape counts ---
    p("## Loop Shape Counts\n")
    p("| Count | Shape | Ranges | Behavior |")
    p("|------:|-------|--------|----------|")
    shape_counts = defaultdict(int)
    for nest in nests:
        shape_counts[nest.shape_key] += 1
    for (shape, ranges, behavior), count in sorted(shape_counts.items(), key=lambda x: -x[1]):
        ranges_str = ", ".join(ranges)
        beh = behavior if behavior != "compute" else ""
        p(f"| {count} | `{shape}` | {ranges_str} | {beh} |")
    p()

    # --- Loop patterns (collapsed) ---
    p("## Loop Patterns (collapsed)\n")
    p("> **Collapsed**: if ANY array in the nest is `U` for a role, the whole nest is `U` for that role.\n")
    p("| Count | Shape | Ranges | Behavior | Collapsed S/U |")
    p("|------:|-------|--------|----------|---------------|")
    pattern_groups = defaultdict(list)
    for nest in nests:
        pattern_groups[nest.pattern_key].append(nest)
    for key, group in sorted(pattern_groups.items(), key=lambda x: -len(x[1])):
        shape, ranges, behavior, collapsed = key
        ranges_str = ", ".join(ranges)
        beh = behavior if behavior != "compute" else ""
        p(f"| {len(group)} | `{shape}` | {ranges_str} | {beh} | `{collapsed}` |")
    p()

    # --- Detailed loop nests ---
    p("## Detailed Loop Nests\n")
    p(f"{len(nests)} nests found.\n")

    for i, nest in enumerate(nests):
        beh = f" `{nest.behavior}`" if nest.behavior != "compute" else ""
        p(f"### Nest {i+1}: `{nest.loop_shape}` ({', '.join(nest.ranges)}){beh}\n")
        p(f"Collapsed: `{nest.collapsed_su}`\n")

        p(f"```fortran\n{nest.source_code}\n```\n")

        p("<details><summary>Array access table</summary>\n")
        p("| Array | Role:S/U |")
        p("|-------|----------|")
        seen = set()
        for acc in nest.accesses:
            key = (acc.name, acc.role_su_signature)
            if key not in seen:
                seen.add(key)
                p(f"| `{acc.name}` | `{acc.role_su_signature}` |")
        p("\n</details>\n")


# ---------------------------------------------------------------------------
# Single-block view
# ---------------------------------------------------------------------------

def _strip_block(role_su_sig):
    parts = [p.strip() for p in role_su_sig.split("  ") if p.strip()]
    filtered = [p for p in parts if not p.startswith("b:")]
    return "  ".join(filtered) if filtered else "-"

def _strip_block_collapsed(collapsed):
    parts = [p.strip() for p in collapsed.split("  ") if p.strip()]
    filtered = [p for p in parts if not p.startswith("b:")]
    return "  ".join(filtered) if filtered else "-"

def _strip_block_shape(shape):
    parts = [p for p in shape.split(".") if p != "b"]
    return ".".join(parts) if parts else "-"


def print_single_block_view(routine_name, nests):
    p = print

    p(f"\n---\n")
    p(f"# Single-Block View (`nblks=1`)\n")
    p(f"Block dimension stripped: `b.h` becomes `h`, `b:S`/`b:U` removed from all signatures.\n")

    pattern_groups_nb = defaultdict(list)
    for nest in nests:
        shape_nb = _strip_block_shape(nest.loop_shape)
        ranges_nb = tuple(r for r, v in zip(nest.ranges, nest.vars) if _role(v) != "b")
        collapsed_nb = _strip_block_collapsed(nest.collapsed_su)
        key = (shape_nb, ranges_nb, nest.behavior, collapsed_nb)
        pattern_groups_nb[key].append(nest)

    p("## Loop Patterns (single-block, collapsed)\n")
    p("| Count | Shape | Ranges | Behavior | Collapsed S/U |")
    p("|------:|-------|--------|----------|---------------|")
    for key, group in sorted(pattern_groups_nb.items(), key=lambda x: -len(x[1])):
        shape_nb, ranges_nb, behavior, collapsed_nb = key
        ranges_str = ", ".join(ranges_nb) if ranges_nb else "-"
        beh = behavior if behavior != "compute" else ""
        p(f"| {len(group)} | `{shape_nb}` | {ranges_str} | {beh} | `{collapsed_nb}` |")
    p()

    p("## Example Nest per Pattern Type\n")
    p("One representative nest from each pattern (the one touching the most arrays).\n")

    for key, group in sorted(pattern_groups_nb.items(), key=lambda x: -len(x[1])):
        shape_nb, ranges_nb, behavior, collapsed_nb = key
        ranges_str = ", ".join(ranges_nb) if ranges_nb else "-"
        beh = f" `{behavior}`" if behavior != "compute" else ""

        example = max(group, key=lambda n: len(set(a.name for a in n.accesses)))
        n_unique = len(set(a.name for a in example.accesses))

        p(f"### `{shape_nb}` {ranges_str}{beh} — collapsed `{collapsed_nb}` ({len(group)}x)\n")
        p(f"{n_unique} unique arrays.\n")
        p(f"```fortran\n{example.source_code}\n```\n")

        p("<details><summary>Array access table</summary>\n")
        p("| Array | Role:S/U (no block) |")
        p("|-------|---------------------|")
        seen = set()
        for acc in example.accesses:
            sig_nb = _strip_block(acc.role_su_signature)
            akey = (acc.name, sig_nb)
            if akey not in seen:
                seen.add(akey)
                p(f"| `{acc.name}` | `{sig_nb}` |")
        p("\n</details>\n")

    # Widest nest
    p("## Widest Nest (most unique array references)\n")
    widest = max(nests, key=lambda n: len(set(a.name for a in n.accesses)))
    n_unique = len(set(a.name for a in widest.accesses))
    beh = f" `{widest.behavior}`" if widest.behavior != "compute" else ""
    shape_nb = _strip_block_shape(widest.loop_shape)
    collapsed_nb = _strip_block_collapsed(widest.collapsed_su)

    p(f"**{n_unique} unique arrays** accessed in a single nest.\n")
    p(f"- Shape: `{widest.loop_shape}` (single-block: `{shape_nb}`)")
    p(f"- Collapsed: `{collapsed_nb}`{beh}\n")
    p(f"```fortran\n{widest.source_code}\n```\n")

    p("| Array | Role:S/U (full) | Role:S/U (no block) |")
    p("|-------|-----------------|---------------------|")
    seen = set()
    for acc in widest.accesses:
        sig_nb = _strip_block(acc.role_su_signature)
        akey = (acc.name, acc.role_su_signature)
        if akey not in seen:
            seen.add(akey)
            p(f"| `{acc.name}` | `{acc.role_su_signature}` | `{sig_nb}` |")
    p()

    # Array groups without block
    p("## Array Groups (single-block)\n")
    array_sigs_nb = defaultdict(set)
    for nest in nests:
        for acc in nest.accesses:
            sig_nb = _strip_block(acc.role_su_signature)
            array_sigs_nb[acc.name].add(sig_nb)

    groups_nb = defaultdict(set)
    for name, sigs in array_sigs_nb.items():
        groups_nb[frozenset(sigs)].add(name)

    for sigs, names in sorted(groups_nb.items(), key=lambda x: -len(x[1])):
        sig_list = sorted(sigs)
        has_u = any("U" in s for s in sig_list)
        tag = "HAS UNSTRUCTURED" if has_u else "FULLY STRUCTURED"
        p(f"### {tag} ({len(names)} arrays)\n")
        p("Patterns: " + ", ".join(f"`{s}`" for s in sig_list) + "\n")
        for n in sorted(names):
            p(f"- `{n}`")
        p()

    # Per-array table
    p("## Per-Array Summary (single-block)\n")
    p("| Array | Role:S/U Patterns | Conflict |")
    p("|-------|-------------------|----------|")
    for name in sorted(array_sigs_nb.keys()):
        sigs = sorted(array_sigs_nb[name])
        has_u = any("U" in s for s in sigs)
        has_s = any(all(c in ("S", "C", " ", ":", "-") for c in s) for s in sigs)
        marker = "S+U" if (has_u and has_s) else ""
        sigs_str = " \\| ".join(f"`{s}`" for s in sigs)
        p(f"| `{name}` | {sigs_str} | {marker} |")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <fortran_file> [subroutine_name]")
        sys.exit(1)

    filepath = Path(sys.argv[1])
    routine_name = sys.argv[2] if len(sys.argv) > 2 else "velocity_tendencies"

    print(f"Parsing {filepath} ...")
    src = Sourcefile.from_file(
        str(filepath), preprocess=True, defines=["_OPENACC"],
        includes=[str(filepath.parent / "include"), str(filepath.parent)],
    )

    routine = None
    for r in src.all_subroutines:
        if r.name.lower() == routine_name.lower():
            routine = r
            break

    if routine is None:
        print(f"ERROR: '{routine_name}' not found. Available: {[r.name for r in src.all_subroutines]}")
        sys.exit(1)

    print(f"Found: {routine.name}  ({len(routine.variables)} vars, "
          f"{len(FindNodes(Loop).visit(routine.body))} loops)")

    nests = analyze_loop_nests(routine)
    array_role_sigs, role_groups = group_arrays(nests)

    # ── Compute set analysis ──
    sa = compute_set_analysis(nests)

    # ── Write reports ──
    import io, builtins

    buf = io.StringIO()
    old_print = builtins.print

    def md_print(*args, **kwargs):
        kwargs["file"] = buf
        old_print(*args, **kwargs)

    builtins.print = md_print

    print_analysis(routine.name, nests, array_role_sigs, role_groups)
    print_single_block_view(routine.name, nests)
    print_set_analysis(routine.name, nests, sa)

    builtins.print = old_print

    outpath = filepath.with_suffix(".analysis.md")
    outpath.write_text(buf.getvalue())
    old_print(f"Report written to {outpath}")
    old_print(buf.getvalue())


if __name__ == "__main__":/*
 * bench_ddt_vn_gpu_sweep.cu -- ddt_vn_apc_pc GPU stencil benchmark
 *
 * Group design (all sizes 2/4/8):
 *   grpA[NA=4, N, nlevp1]   = {vn, vt, vn_ie, z_vt_ie}
 *   grpD[ND=2, N_c, nlev]   = {z_ekinh, z_w_con_c_full}
 *   Standalone: z_kin_hor_e, ddqz_z_full_e, f_e, zeta, out, vn_ie(for Full AoS)
 *   Pair-contiguous (IP(je,n)=n+2*je for AoS/grp/AoSoA):
 *     coeff_gradekin, c_lin_e, cell_idx, vert_idx
 *
 * 5 layouts: SoA, FullAoS, Grouped, AoSoA-32, AoSoA-64
 * Sweeps threadblock, coarsening, connectivity. 100 reps, L2 flush.
 *
 * Compile: nvcc -O3 -std=c++17 -arch=sm_90 bench_ddt_vn_gpu_sweep.cu -o bench
 */
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <random>
#include <vector>
#include "icon_data_loader.h"

#include <cuda_runtime.h>
#define GPU_CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){fprintf(stderr,"CUDA %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e));exit(1);}}while(0)
#define GM cudaMalloc
#define GF cudaFree
#define GC(d,s,n) GPU_CHECK(cudaMemcpy(d,s,n,cudaMemcpyHostToDevice))
#define GCD(d,s,n) GPU_CHECK(cudaMemcpy(d,s,n,cudaMemcpyDeviceToHost))
#define GS cudaDeviceSynchronize
#define GE cudaEvent_t
#define GEC cudaEventCreate
#define GER cudaEventRecord
#define GES cudaEventSynchronize
#define GEE cudaEventElapsedTime
#define GED cudaEventDestroy
#define GSET cudaMemset
#define HD __host__ __device__ __forceinline__

static constexpr int NPROMA=81920,NRUNS=100,WARMUP=10;
static constexpr int A_VN=0,A_VT=1,A_VN_IE=2,A_VT_IE=3,NA=4;
static constexpr int D_EKINH=0,D_WCON=1,ND=2;
static constexpr int E3_VT=0,E3_DDQZ=1,NE3=2;
static constexpr int E2_CG0=0,E2_CG1=1,E2_CL0=2,E2_CL1=3,NE2=4;
static constexpr int CN_CI0=0,CN_CI1=1,CN_VI0=2,CN_VI1=3,NCONN=4;

static inline uint64_t splitmix64(uint64_t x){x+=0x9E3779B97F4A7C15ULL;x=(x^(x>>30))*0xBF58476D1CE4E5B9ULL;x=(x^(x>>27))*0x94D049BB133111EBULL;return x^(x>>31);}
static void fill_rand(double*a,size_t n,unsigned s){for(size_t i=0;i<n;i++){uint64_t h=splitmix64((uint64_t)s*2654435761ULL+i);a[i]=(double)(int64_t)(h&0xFFFFF)/100000.0-5.0;}}

HD int I2(int je,int jk,int N){return je+jk*N;}
HD int IX(int je,int n,int N){return je+n*N;}       /* SoA stride-N pairs */
HD int IP(int je,int n){return n+2*je;}              /* pair-contiguous */
HD int IA(int f,int je,int jk,int N){return f+NA*(je+N*jk);}
HD int ID(int f,int c,int jk,int Nc){return f+ND*(c+Nc*jk);}
HD int IE3(int f,int je,int jk,int N){return f+NE3*(je+N*jk);}
HD int IE2(int f,int je){return f+NE2*je;}
HD int ICN(int f,int je){return f+NCONN*je;}
template<int V> HD int IAao(int f,int je,int jk,int N){int t=(N+V-1)/V;return jk*t*NA*V+(je/V)*NA*V+f*V+(je%V);}
template<int V> HD int IDao(int f,int c,int jk,int Nc){int t=(Nc+V-1)/V;return jk*t*ND*V+(c/V)*ND*V+f*V+(c%V);}
template<int V> size_t szAao(int N,int K){return(size_t)((N+V-1)/V)*NA*V*K;}
template<int V> size_t szDao(int Nc,int K){return(size_t)((Nc+V-1)/V)*ND*V*K;}

/* Stencil bodies — uses STENCIL_EXPR from common locals */
#define STENCIL_EXPR(oi) out[oi]=-(ekin_e*(cg0-cg1)+cg1*eh1-cg0*eh0+vt_e*(fe_e+0.5*(zeta0+zeta1))+(cl0*w0+cl1*w1)*(vk-vk1)/ddqz_e);

#define SOA_BODY() {int c2=I2(je,jk,N);\
    int ci0=cell_idx[IX(je,0,N)],ci1=cell_idx[IX(je,1,N)],vi0=vert_idx[IX(je,0,N)],vi1=vert_idx[IX(je,1,N)];\
    double cg0=coeff_gradekin[IX(je,0,N)],cg1=coeff_gradekin[IX(je,1,N)],cl0=c_lin_e[IX(je,0,N)],cl1=c_lin_e[IX(je,1,N)];\
    double ekin_e=z_kin_hor_e[c2],vt_e=vt[c2],ddqz_e=ddqz_z_full_e[c2],fe_e=f_e[je];\
    double eh0=z_ekinh[I2(ci0,jk,N_c)],eh1=z_ekinh[I2(ci1,jk,N_c)],w0=z_w_con_c_full[I2(ci0,jk,N_c)],w1=z_w_con_c_full[I2(ci1,jk,N_c)];\
    double vk=vn_ie[I2(je,jk,N)],vk1=vn_ie[I2(je,jk+1,N)],zeta0=zeta[I2(vi0,jk,N_v)],zeta1=zeta[I2(vi1,jk,N_v)];\
    STENCIL_EXPR(c2)}

#define AOS_BODY() {\
    int ci0=aos_conn[ICN(CN_CI0,je)],ci1=aos_conn[ICN(CN_CI1,je)],vi0=aos_conn[ICN(CN_VI0,je)],vi1=aos_conn[ICN(CN_VI1,je)];\
    double cg0=aos_e2d[IE2(E2_CG0,je)],cg1=aos_e2d[IE2(E2_CG1,je)],cl0=aos_e2d[IE2(E2_CL0,je)],cl1=aos_e2d[IE2(E2_CL1,je)];\
    int eb=NE3*(je+N*jk);double vt_e=aos_e3d[eb+E3_VT],ddqz_e=aos_e3d[eb+E3_DDQZ];\
    double ekin_e=z_kin_hor_e[I2(je,jk,N)],fe_e=f_e[je];\
    double eh0=aos_cell[ID(D_EKINH,ci0,jk,N_c)],eh1=aos_cell[ID(D_EKINH,ci1,jk,N_c)],w0=aos_cell[ID(D_WCON,ci0,jk,N_c)],w1=aos_cell[ID(D_WCON,ci1,jk,N_c)];\
    double vk=vn_ie[I2(je,jk,N)],vk1=vn_ie[I2(je,jk+1,N)],zeta0=zeta[I2(vi0,jk,N_v)],zeta1=zeta[I2(vi1,jk,N_v)];\
    STENCIL_EXPR(I2(je,jk,N))}

#define GRP_BODY() {\
    int ci0=cell_idx[IP(je,0)],ci1=cell_idx[IP(je,1)],vi0=vert_idx[IP(je,0)],vi1=vert_idx[IP(je,1)];\
    double cg0=coeff_gradekin[IP(je,0)],cg1=coeff_gradekin[IP(je,1)],cl0=c_lin_e[IP(je,0)],cl1=c_lin_e[IP(je,1)];\
    double vt_e=grpA[IA(A_VT,je,jk,N)],vk=grpA[IA(A_VN_IE,je,jk,N)],vk1=grpA[IA(A_VN_IE,je,jk+1,N)];\
    double ekin_e=z_kin_hor_e[I2(je,jk,N)],ddqz_e=ddqz_z_full_e[I2(je,jk,N)],fe_e=f_e[je];\
    double eh0=grpD[ID(D_EKINH,ci0,jk,N_c)],eh1=grpD[ID(D_EKINH,ci1,jk,N_c)],w0=grpD[ID(D_WCON,ci0,jk,N_c)],w1=grpD[ID(D_WCON,ci1,jk,N_c)];\
    double zeta0=zeta[I2(vi0,jk,N_v)],zeta1=zeta[I2(vi1,jk,N_v)];\
    STENCIL_EXPR(I2(je,jk,N))}

#define AOSOA_BODY(V) {\
    int ci0=cell_idx[IP(je,0)],ci1=cell_idx[IP(je,1)],vi0=vert_idx[IP(je,0)],vi1=vert_idx[IP(je,1)];\
    double cg0=coeff_gradekin[IP(je,0)],cg1=coeff_gradekin[IP(je,1)],cl0=c_lin_e[IP(je,0)],cl1=c_lin_e[IP(je,1)];\
    double vt_e=grpA[IAao<V>(A_VT,je,jk,N)],vk=grpA[IAao<V>(A_VN_IE,je,jk,N)],vk1=grpA[IAao<V>(A_VN_IE,je,jk+1,N)];\
    double ekin_e=z_kin_hor_e[I2(je,jk,N)],ddqz_e=ddqz_z_full_e[I2(je,jk,N)],fe_e=f_e[je];\
    double eh0=grpD[IDao<V>(D_EKINH,ci0,jk,N_c)],eh1=grpD[IDao<V>(D_EKINH,ci1,jk,N_c)],w0=grpD[IDao<V>(D_WCON,ci0,jk,N_c)],w1=grpD[IDao<V>(D_WCON,ci1,jk,N_c)];\
    double zeta0=zeta[I2(vi0,jk,N_v)],zeta1=zeta[I2(vi1,jk,N_v)];\
    STENCIL_EXPR(I2(je,jk,N))}

enum CellDist{UNIFORM=0,NORMAL1=1,NORMAL4=2,SEQUENTIAL=3,EXACT=4};
static const char*dist_names[]={"uniform","normal_var1","normal_var4","sequential","exact"};
static void gen_conn(int*L,int N,int Nt,CellDist d,std::mt19937&rng){switch(d){case UNIFORM:{std::uniform_int_distribution<int>u(0,Nt-1);for(int i=0;i<N;i++){L[i*2]=u(rng);L[i*2+1]=u(rng);}break;}case NORMAL1:{std::normal_distribution<double>nd(0,1);for(int i=0;i<N;i++){L[i*2]=((i+1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;L[i*2+1]=((i-1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;}break;}case NORMAL4:{std::normal_distribution<double>nd(0,2);for(int i=0;i<N;i++){L[i*2]=((i+1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;L[i*2+1]=((i-1+(int)std::round(nd(rng)))%Nt+Nt)%Nt;}break;}case SEQUENTIAL:for(int i=0;i<N;i++){L[i*2]=(i+1)%Nt;L[i*2+1]=(i+1)%Nt;}break;}}
static void pack_pairs(double*d,const double*s,int N){for(int i=0;i<N;i++){d[2*i]=s[i];d[2*i+1]=s[i+N];}}
static void pack_pairs_i(int*d,const int*s,int N){for(int i=0;i<N;i++){d[2*i]=s[i];d[2*i+1]=s[i+N];}}

static constexpr size_t FLUSH_N=48ULL*1024*1024;static double*d_fb=nullptr;
__global__ void fk(double*b,size_t n){size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;if(i<n)b[i]=b[i]*1.00001+1e-12;}
static void fi(){GPU_CHECK(GM(&d_fb,FLUSH_N*8));GPU_CHECK(GSET(d_fb,0,FLUSH_N*8));fk<<<(int)((FLUSH_N+255)/256),256>>>(d_fb,FLUSH_N);GPU_CHECK(GS());}
static void fl(){fk<<<(int)((FLUSH_N+255)/256),256>>>(d_fb,FLUSH_N);GPU_CHECK(GS());}
static void fd(){GF(d_fb);}

/* Host data */
struct SoAH{int N,Nc,Nv,nl;size_t se,sc,sv;
    double*vt,*vni,*ek,*dq,*cg,*cl,*fe,*eh,*wc,*zt;int*ci,*vi;double*out;
    void alloc(int N_,int c,int v,int l){N=N_;Nc=c;Nv=v;nl=l;se=(size_t)N*nl;sc=(size_t)c*nl;sv=(size_t)v*nl;
        vt=new double[se];vni=new double[(size_t)N*(nl+1)];ek=new double[se];dq=new double[se];
        cg=new double[N*2];cl=new double[N*2];fe=new double[N];eh=new double[sc];wc=new double[sc];zt=new double[sv];
        ci=new int[N*2];vi=new int[N*2];out=new double[se];}
    void fill(){fill_rand(vt,se,101);fill_rand(vni,(size_t)N*(nl+1),102);fill_rand(ek,se,103);fill_rand(dq,se,104);
        fill_rand(cg,N*2,105);fill_rand(cl,N*2,106);fill_rand(fe,N,107);fill_rand(eh,sc,108);fill_rand(wc,sc,109);fill_rand(zt,sv,110);
        memset(out,0,se*8);for(size_t i=0;i<se;i++)if(std::abs(dq[i])<1e-10)dq[i]=1.0;}
    void free_all(){delete[]vt;delete[]vni;delete[]ek;delete[]dq;delete[]cg;delete[]cl;delete[]fe;delete[]eh;delete[]wc;delete[]zt;delete[]ci;delete[]vi;delete[]out;}};

struct AoSH{int N,Nc,Nv,nl;size_t se,s3,s2,sn,sl;double*e3,*vni,*ek,*fe,*e2;int*cn;double*cd,*zt,*out;
    void from(const SoAH&s){N=s.N;Nc=s.Nc;Nv=s.Nv;nl=s.nl;se=s.se;s3=(size_t)NE3*N*nl;s2=(size_t)NE2*N;sn=(size_t)NCONN*N;sl=(size_t)ND*Nc*nl;
        e3=new double[s3];for(int jk=0;jk<nl;jk++)for(int je=0;je<N;je++){e3[IE3(E3_VT,je,jk,N)]=s.vt[I2(je,jk,N)];e3[IE3(E3_DDQZ,je,jk,N)]=s.dq[I2(je,jk,N)];}
        vni=new double[(size_t)N*(nl+1)];memcpy(vni,s.vni,(size_t)N*(nl+1)*8);ek=new double[se];memcpy(ek,s.ek,se*8);fe=new double[N];memcpy(fe,s.fe,N*8);
        e2=new double[s2];for(int je=0;je<N;je++){e2[IE2(E2_CG0,je)]=s.cg[IX(je,0,N)];e2[IE2(E2_CG1,je)]=s.cg[IX(je,1,N)];e2[IE2(E2_CL0,je)]=s.cl[IX(je,0,N)];e2[IE2(E2_CL1,je)]=s.cl[IX(je,1,N)];}
        cn=new int[sn];for(int je=0;je<N;je++){cn[ICN(CN_CI0,je)]=s.ci[IX(je,0,N)];cn[ICN(CN_CI1,je)]=s.ci[IX(je,1,N)];cn[ICN(CN_VI0,je)]=s.vi[IX(je,0,N)];cn[ICN(CN_VI1,je)]=s.vi[IX(je,1,N)];}
        cd=new double[sl];for(int jk=0;jk<nl;jk++)for(int ic=0;ic<Nc;ic++){cd[ID(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)];cd[ID(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}
        zt=new double[(size_t)Nv*nl];memcpy(zt,s.zt,(size_t)Nv*nl*8);out=new double[se];memset(out,0,se*8);}
    void free_all(){delete[]e3;delete[]vni;delete[]ek;delete[]fe;delete[]e2;delete[]cn;delete[]cd;delete[]zt;delete[]out;}};

struct GrpH{int N,Nc,Nv,nl;size_t se;double*gA,*gD,*ek,*dq,*fe,*cg,*cl,*zt;int*ci,*vi;double*out;
    void from(const SoAH&s){N=s.N;Nc=s.Nc;Nv=s.Nv;nl=s.nl;se=s.se;int np=nl+1;
        gA=new double[(size_t)NA*N*np];memset(gA,0,(size_t)NA*N*np*8);
        for(int jk=0;jk<nl;jk++)for(int je=0;je<N;je++){gA[IA(A_VT,je,jk,N)]=s.vt[I2(je,jk,N)];gA[IA(A_VN_IE,je,jk,N)]=s.vni[I2(je,jk,N)];}
        for(int je=0;je<N;je++)gA[IA(A_VN_IE,je,nl,N)]=s.vni[I2(je,nl,N)];
        gD=new double[(size_t)ND*Nc*nl];for(int jk=0;jk<nl;jk++)for(int ic=0;ic<Nc;ic++){gD[ID(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)];gD[ID(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}
        ek=new double[se];memcpy(ek,s.ek,se*8);dq=new double[se];memcpy(dq,s.dq,se*8);fe=new double[N];memcpy(fe,s.fe,N*8);
        zt=new double[(size_t)Nv*nl];memcpy(zt,s.zt,(size_t)Nv*nl*8);
        cg=new double[N*2];pack_pairs(cg,s.cg,N);cl=new double[N*2];pack_pairs(cl,s.cl,N);
        ci=new int[N*2];pack_pairs_i(ci,s.ci,N);vi=new int[N*2];pack_pairs_i(vi,s.vi,N);
        out=new double[se];memset(out,0,se*8);}
    void free_all(){delete[]gA;delete[]gD;delete[]ek;delete[]dq;delete[]fe;delete[]cg;delete[]cl;delete[]zt;delete[]ci;delete[]vi;delete[]out;}};

template<int V>struct AoSoAH{int N,Nc,Nv,nl;size_t se,sA,sD;double*gA,*gD,*ek,*dq,*fe,*cg,*cl,*zt;int*ci,*vi;double*out;
    void from(const SoAH&s){N=s.N;Nc=s.Nc;Nv=s.Nv;nl=s.nl;se=s.se;int np=nl+1;sA=szAao<V>(N,np);sD=szDao<V>(Nc,nl);
        gA=new double[sA];memset(gA,0,sA*8);for(int jk=0;jk<nl;jk++)for(int je=0;je<N;je++){gA[IAao<V>(A_VT,je,jk,N)]=s.vt[I2(je,jk,N)];gA[IAao<V>(A_VN_IE,je,jk,N)]=s.vni[I2(je,jk,N)];}
        for(int je=0;je<N;je++)gA[IAao<V>(A_VN_IE,je,nl,N)]=s.vni[I2(je,nl,N)];
        gD=new double[sD];memset(gD,0,sD*8);for(int jk=0;jk<nl;jk++)for(int ic=0;ic<Nc;ic++){gD[IDao<V>(D_EKINH,ic,jk,Nc)]=s.eh[I2(ic,jk,Nc)];gD[IDao<V>(D_WCON,ic,jk,Nc)]=s.wc[I2(ic,jk,Nc)];}
        ek=new double[se];memcpy(ek,s.ek,se*8);dq=new double[se];memcpy(dq,s.dq,se*8);fe=new double[N];memcpy(fe,s.fe,N*8);
        zt=new double[(size_t)Nv*nl];memcpy(zt,s.zt,(size_t)Nv*nl*8);
        cg=new double[N*2];pack_pairs(cg,s.cg,N);cl=new double[N*2];pack_pairs(cl,s.cl,N);
        ci=new int[N*2];pack_pairs_i(ci,s.ci,N);vi=new int[N*2];pack_pairs_i(vi,s.vi,N);
        out=new double[se];memset(out,0,se*8);}
    void free_all(){delete[]gA;delete[]gD;delete[]ek;delete[]dq;delete[]fe;delete[]cg;delete[]cl;delete[]zt;delete[]ci;delete[]vi;delete[]out;}};

/* Device structs */
struct DS{double*vt,*vni,*ek,*dq,*cg,*cl,*fe,*eh,*wc,*zt,*out;int*ci,*vi;};
struct DA{double*e3,*vni,*ek,*fe,*e2,*cd,*zt,*out;int*cn;};
struct DG{double*gA,*gD,*ek,*dq,*fe,*cg,*cl,*zt,*out;int*ci,*vi;};
struct DO_{double*gA,*gD,*ek,*dq,*fe,*cg,*cl,*zt,*out;int*ci,*vi;size_t sA,sD;};

/* GPU kernels */
#define KL(TX,TY,B) int jb=((int)blockIdx.x*(int)blockDim.x+(int)threadIdx.x)*TX,kb=((int)blockIdx.y*(int)blockDim.y+(int)threadIdx.y)*TY;\
    _Pragma("unroll")for(int dy=0;dy<TY;dy++){int jk=kb+dy;if(jk>=nlev)break;_Pragma("unroll")for(int dx=0;dx<TX;dx++){int je=jb+dx;if(je>=N)break;B}}

template<int TX,int TY>__global__ void ks(double*__restrict__ out,const double*__restrict__ vt,const double*__restrict__ vn_ie,const double*__restrict__ z_kin_hor_e,const double*__restrict__ ddqz_z_full_e,const double*__restrict__ coeff_gradekin,const double*__restrict__ c_lin_e,const double*__restrict__ f_e,const double*__restrict__ z_ekinh,const double*__restrict__ z_w_con_c_full,const double*__restrict__ zeta,const int*__restrict__ cell_idx,const int*__restrict__ vert_idx,int N,int N_c,int N_v,int nlev){KL(TX,TY,SOA_BODY())}
template<int TX,int TY>__global__ void ka(double*__restrict__ out,const double*__restrict__ aos_e3d,const double*__restrict__ vn_ie,const double*__restrict__ z_kin_hor_e,const double*__restrict__ f_e,const double*__restrict__ aos_e2d,const int*__restrict__ aos_conn,const double*__restrict__ aos_cell,const double*__restrict__ zeta,int N,int N_c,int N_v,int nlev){KL(TX,TY,AOS_BODY())}
template<int TX,int TY>__global__ void kg(double*__restrict__ out,const double*__restrict__ grpA,const double*__restrict__ grpD,const double*__restrict__ z_kin_hor_e,const double*__restrict__ ddqz_z_full_e,const double*__restrict__ f_e,const double*__restrict__ coeff_gradekin,const double*__restrict__ c_lin_e,const double*__restrict__ zeta,const int*__restrict__ cell_idx,const int*__restrict__ vert_idx,int N,int N_c,int N_v,int nlev){KL(TX,TY,GRP_BODY())}
template<int V,int TX,int TY>__global__ void ko(double*__restrict__ out,const double*__restrict__ grpA,const double*__restrict__ grpD,const double*__restrict__ z_kin_hor_e,const double*__restrict__ ddqz_z_full_e,const double*__restrict__ f_e,const double*__restrict__ coeff_gradekin,const double*__restrict__ c_lin_e,const double*__restrict__ zeta,const int*__restrict__ cell_idx,const int*__restrict__ vert_idx,int N,int N_c,int N_v,int nlev){KL(TX,TY,AOSOA_BODY(V))}

static bool verify(const double*g,const double*r,size_t n,int&nf_,double&mr){nf_=0;mr=0;for(size_t i=0;i<n;i++){double d=std::abs(g[i]-r[i]),dn=std::max(std::abs(r[i]),1e-300),rv=d/dn;if(rv>mr)mr=rv;if(d>1e-12+1e-8*std::abs(r[i]))nf_++;}return nf_==0;}

struct KC{int bx,by,tx,ty;};
static const KC cfgs[]={{32,1,1,1},{64,1,1,1},{128,1,1,1},{256,1,1,1},{64,2,1,1},{128,2,1,1},{64,4,1,1},{32,4,1,1},{32,8,1,1},{64,1,2,1},{128,1,2,1},{64,2,2,1},{64,1,4,1},{128,1,4,1},{64,1,1,2},{128,1,1,2},{64,2,1,2},{64,1,1,4},{128,1,1,4},{64,1,2,2},{128,1,2,2},{64,2,2,2},{64,1,4,2},{64,1,2,4}};
static constexpr int NC_=sizeof(cfgs)/sizeof(cfgs[0]);

#define DT(K,tx,ty,g,b,...) do{if(tx==1&&ty==1)K<1,1><<<g,b>>>(__VA_ARGS__);else if(tx==2&&ty==1)K<2,1><<<g,b>>>(__VA_ARGS__);else if(tx==4&&ty==1)K<4,1><<<g,b>>>(__VA_ARGS__);else if(tx==1&&ty==2)K<1,2><<<g,b>>>(__VA_ARGS__);else if(tx==1&&ty==4)K<1,4><<<g,b>>>(__VA_ARGS__);else if(tx==2&&ty==2)K<2,2><<<g,b>>>(__VA_ARGS__);else if(tx==4&&ty==2)K<4,2><<<g,b>>>(__VA_ARGS__);else if(tx==2&&ty==4)K<2,4><<<g,b>>>(__VA_ARGS__);}while(0)
#define DTV(K,V,tx,ty,g,b,...) do{if(tx==1&&ty==1)K<V,1,1><<<g,b>>>(__VA_ARGS__);else if(tx==2&&ty==1)K<V,2,1><<<g,b>>>(__VA_ARGS__);else if(tx==4&&ty==1)K<V,4,1><<<g,b>>>(__VA_ARGS__);else if(tx==1&&ty==2)K<V,1,2><<<g,b>>>(__VA_ARGS__);else if(tx==1&&ty==4)K<V,1,4><<<g,b>>>(__VA_ARGS__);else if(tx==2&&ty==2)K<V,2,2><<<g,b>>>(__VA_ARGS__);else if(tx==4&&ty==2)K<V,4,2><<<g,b>>>(__VA_ARGS__);else if(tx==2&&ty==4)K<V,2,4><<<g,b>>>(__VA_ARGS__);}while(0)

#define AL(p,n) GPU_CHECK(GM(&(p),(n)))
static void ads(DS&d,const SoAH&h){int N=h.N,nl=h.nl;AL(d.vt,h.se*8);AL(d.vni,(size_t)N*(nl+1)*8);AL(d.ek,h.se*8);AL(d.dq,h.se*8);AL(d.cg,N*2*8);AL(d.cl,N*2*8);AL(d.fe,N*8);AL(d.eh,h.sc*8);AL(d.wc,h.sc*8);AL(d.zt,h.sv*8);AL(d.ci,N*2*4);AL(d.vi,N*2*4);AL(d.out,h.se*8);GC(d.vt,h.vt,h.se*8);GC(d.vni,h.vni,(size_t)N*(nl+1)*8);GC(d.ek,h.ek,h.se*8);GC(d.dq,h.dq,h.se*8);GC(d.cg,h.cg,N*2*8);GC(d.cl,h.cl,N*2*8);GC(d.fe,h.fe,N*8);GC(d.eh,h.eh,h.sc*8);GC(d.wc,h.wc,h.sc*8);GC(d.zt,h.zt,h.sv*8);GC(d.ci,h.ci,N*2*4);GC(d.vi,h.vi,N*2*4);}
static void fds(DS&d){GF(d.vt);GF(d.vni);GF(d.ek);GF(d.dq);GF(d.cg);GF(d.cl);GF(d.fe);GF(d.eh);GF(d.wc);GF(d.zt);GF(d.ci);GF(d.vi);GF(d.out);}
static void ada(DA&d,const AoSH&h){int N=h.N,nl=h.nl;AL(d.e3,h.s3*8);AL(d.vni,(size_t)N*(nl+1)*8);AL(d.ek,h.se*8);AL(d.fe,N*8);AL(d.e2,h.s2*8);AL(d.cn,h.sn*4);AL(d.cd,h.sl*8);AL(d.zt,(size_t)h.Nv*nl*8);AL(d.out,h.se*8);GC(d.e3,h.e3,h.s3*8);GC(d.vni,h.vni,(size_t)N*(nl+1)*8);GC(d.ek,h.ek,h.se*8);GC(d.fe,h.fe,N*8);GC(d.e2,h.e2,h.s2*8);GC(d.cn,h.cn,h.sn*4);GC(d.cd,h.cd,h.sl*8);GC(d.zt,h.zt,(size_t)h.Nv*nl*8);}
static void fda(DA&d){GF(d.e3);GF(d.vni);GF(d.ek);GF(d.fe);GF(d.e2);GF(d.cn);GF(d.cd);GF(d.zt);GF(d.out);}
static void adg(DG&d,const GrpH&h){int N=h.N,Nc=h.Nc,Nv=h.Nv,nl=h.nl;AL(d.gA,(size_t)NA*N*(nl+1)*8);AL(d.gD,(size_t)ND*Nc*nl*8);AL(d.ek,h.se*8);AL(d.dq,h.se*8);AL(d.fe,N*8);AL(d.cg,N*2*8);AL(d.cl,N*2*8);AL(d.zt,(size_t)Nv*nl*8);AL(d.ci,N*2*4);AL(d.vi,N*2*4);AL(d.out,h.se*8);GC(d.gA,h.gA,(size_t)NA*N*(nl+1)*8);GC(d.gD,h.gD,(size_t)ND*Nc*nl*8);GC(d.ek,h.ek,h.se*8);GC(d.dq,h.dq,h.se*8);GC(d.fe,h.fe,N*8);GC(d.cg,h.cg,N*2*8);GC(d.cl,h.cl,N*2*8);GC(d.zt,h.zt,(size_t)Nv*nl*8);GC(d.ci,h.ci,N*2*4);GC(d.vi,h.vi,N*2*4);}
static void fdg(DG&d){GF(d.gA);GF(d.gD);GF(d.ek);GF(d.dq);GF(d.fe);GF(d.cg);GF(d.cl);GF(d.zt);GF(d.ci);GF(d.vi);GF(d.out);}
template<int V>static void ado(DO_&d,const AoSoAH<V>&h){int N=h.N,Nv=h.Nv,nl=h.nl;d.sA=h.sA;d.sD=h.sD;AL(d.gA,h.sA*8);AL(d.gD,h.sD*8);AL(d.ek,h.se*8);AL(d.dq,h.se*8);AL(d.fe,N*8);AL(d.cg,N*2*8);AL(d.cl,N*2*8);AL(d.zt,(size_t)Nv*nl*8);AL(d.ci,N*2*4);AL(d.vi,N*2*4);AL(d.out,h.se*8);GC(d.gA,h.gA,h.sA*8);GC(d.gD,h.gD,h.sD*8);GC(d.ek,h.ek,h.se*8);GC(d.dq,h.dq,h.se*8);GC(d.fe,h.fe,N*8);GC(d.cg,h.cg,N*2*8);GC(d.cl,h.cl,N*2*8);GC(d.zt,h.zt,(size_t)Nv*nl*8);GC(d.ci,h.ci,N*2*4);GC(d.vi,h.vi,N*2*4);}
static void fdo(DO_&d){GF(d.gA);GF(d.gD);GF(d.ek);GF(d.dq);GF(d.fe);GF(d.cg);GF(d.cl);GF(d.zt);GF(d.ci);GF(d.vi);GF(d.out);}

#define BENCH(lbl,K_CALL,vw) for(int w=0;w<WARMUP;w++){K_CALL;}GPU_CHECK(GS());\
    for(int r=0;r<NRUNS;r++){fl();GPU_CHECK(GER(e0));K_CALL;GPU_CHECK(GER(e1));GPU_CHECK(GES(e1));float ms;GPU_CHECK(GEE(&ms,e0,e1));fprintf(csv,lbl",%s,%d,%d,%d,%d,"#vw",%d,%.6f\n",dn,c.bx,c.by,c.tx,c.ty,r,(double)ms);}

int main(int argc,char*argv[]){
    int nl=(argc>=2)?atoi(argv[1]):90;int N=NPROMA,Nc=N,Nv=N;
    printf("ddt_vn GPU sweep  N=%d nl=%d NRUNS=%d cfgs=%d NA=%d ND=%d NE3=%d NE2=%d\n",N,nl,NRUNS,NC_,NA,ND,NE3,NE2);
    FILE*csv=fopen("ddt_vn_apc_pc_gpu_sweep.csv","w");
    fprintf(csv,"layout,cell_dist,block_x,block_y,coarsen_x,coarsen_y,vec_width,run_id,time_ms\n");
    std::mt19937 rng(42);fi();GE e0,e1;GPU_CHECK(GEC(&e0));GPU_CHECK(GEC(&e1));

    /* ── Load ICON exact connectivity ─────────────────────────────── */
    int icon_step = 9;
    std::string gp = icon_global_path(icon_step);
    std::string pp = icon_patch_path(icon_step);
    int icon_nproma = icon_read_nproma(gp.c_str());
    if (icon_nproma <= 0)
        fprintf(stderr, "WARNING: could not read nproma from '%s'\n", gp.c_str());
    printf("Loading ICON: %s  (nproma=%d)\n", pp.c_str(), icon_nproma);
    IconEdgeData ied;
    bool have_exact = (icon_nproma > 0) &&
        icon_load_patch(pp.c_str(), icon_nproma, ied);
    if (have_exact)
        printf("ICON: nproma=%d  n_edges=%d (valid=%d)  n_cells=%d  n_verts=%d\n",
               ied.nproma, ied.n_edges, ied.n_edges_valid, ied.n_cells, ied.n_verts);
    else
        printf("ICON data not available — skipping 'exact' distribution\n");

    int *icon_ci = nullptr, *icon_vi = nullptr;
    if (have_exact) {
        icon_ci = new int[N * 2];
        icon_vi = new int[N * 2];
        for (int i = 0; i < N; i++) {
            int src = i % ied.n_edges;
            icon_ci[i * 2 + 0] = ied.cell_idx[src * 2 + 0] % Nc;
            icon_ci[i * 2 + 1] = ied.cell_idx[src * 2 + 1] % Nc;
            icon_vi[i * 2 + 0] = ied.vert_idx[src * 2 + 0] % Nv;
            icon_vi[i * 2 + 1] = ied.vert_idx[src * 2 + 1] % Nv;
        }
        ied.free_all();
    }

    int n_dists = have_exact ? 5 : 4;

    for(int di=0;di<n_dists;di++){printf("\n=== dist=%s ===\n",dist_names[di]);
        SoAH soa;soa.alloc(N,Nc,Nv,nl);soa.fill();
        if (di == EXACT) {
            memcpy(soa.ci, icon_ci, N * 2 * sizeof(int));
            memcpy(soa.vi, icon_vi, N * 2 * sizeof(int));
        } else {
            gen_conn(soa.ci,N,Nc,(CellDist)di,rng);
            gen_conn(soa.vi,N,Nv,(CellDist)di,rng);
        }
        AoSH fao;fao.from(soa);GrpH grp;grp.from(soa);AoSoAH<32>a32;a32.from(soa);AoSoAH<64>a64;a64.from(soa);

        double*hr=new double[soa.se];memset(soa.out,0,soa.se*8);
        for(int jk=0;jk<nl;jk++)for(int je=0;je<N;je++){const double*vt=soa.vt,*vn_ie=soa.vni,*z_kin_hor_e=soa.ek,*ddqz_z_full_e=soa.dq;
            const double*coeff_gradekin=soa.cg,*c_lin_e=soa.cl,*f_e=soa.fe,*z_ekinh=soa.eh,*z_w_con_c_full=soa.wc,*zeta=soa.zt;
            const int*cell_idx=soa.ci,*vert_idx=soa.vi;double*out=soa.out;int N_c=Nc,N_v=Nv;SOA_BODY()}
        memcpy(hr,soa.out,soa.se*8);

        DS ds;ads(ds,soa);DA da;ada(da,fao);DG dg;adg(dg,grp);DO_ d32;ado<32>(d32,a32);DO_ d64;ado<64>(d64,a64);

        {dim3 bk(64,1),gr((N+63)/64,nl);double*ho=new double[soa.se];int nf;double mr;
            ks<1,1><<<gr,bk>>>(ds.out,ds.vt,ds.vni,ds.ek,ds.dq,ds.cg,ds.cl,ds.fe,ds.eh,ds.wc,ds.zt,ds.ci,ds.vi,N,Nc,Nv,nl);GPU_CHECK(GS());GCD(ho,ds.out,soa.se*8);verify(ho,hr,soa.se,nf,mr);printf("  SoA    %s %.2e\n",nf?"FAIL":"OK",mr);
            ka<1,1><<<gr,bk>>>(da.out,da.e3,da.vni,da.ek,da.fe,da.e2,da.cn,da.cd,da.zt,N,Nc,Nv,nl);GPU_CHECK(GS());GCD(ho,da.out,soa.se*8);verify(ho,hr,soa.se,nf,mr);printf("  AoS    %s %.2e\n",nf?"FAIL":"OK",mr);
            kg<1,1><<<gr,bk>>>(dg.out,dg.gA,dg.gD,dg.ek,dg.dq,dg.fe,dg.cg,dg.cl,dg.zt,dg.ci,dg.vi,N,Nc,Nv,nl);GPU_CHECK(GS());GCD(ho,dg.out,soa.se*8);verify(ho,hr,soa.se,nf,mr);printf("  Grp    %s %.2e\n",nf?"FAIL":"OK",mr);
            ko<32,1,1><<<gr,bk>>>(d32.out,d32.gA,d32.gD,d32.ek,d32.dq,d32.fe,d32.cg,d32.cl,d32.zt,d32.ci,d32.vi,N,Nc,Nv,nl);GPU_CHECK(GS());GCD(ho,d32.out,soa.se*8);verify(ho,hr,soa.se,nf,mr);printf("  AoSoA32 %s %.2e\n",nf?"FAIL":"OK",mr);
            ko<64,1,1><<<gr,bk>>>(d64.out,d64.gA,d64.gD,d64.ek,d64.dq,d64.fe,d64.cg,d64.cl,d64.zt,d64.ci,d64.vi,N,Nc,Nv,nl);GPU_CHECK(GS());GCD(ho,d64.out,soa.se*8);verify(ho,hr,soa.se,nf,mr);printf("  AoSoA64 %s %.2e\n",nf?"FAIL":"OK",mr);
            delete[]ho;}

        const char*dn=dist_names[di];
        for(int ci=0;ci<NC_;ci++){auto&c=cfgs[ci];printf("  cfg %2d/%d (%3d,%d) c(%d,%d)...",ci+1,NC_,c.bx,c.by,c.tx,c.ty);fflush(stdout);
            dim3 bk(c.bx,c.by),gr((N+bk.x*c.tx-1)/(bk.x*c.tx),(nl+bk.y*c.ty-1)/(bk.y*c.ty));
            BENCH("soa",DT(ks,c.tx,c.ty,gr,bk,ds.out,ds.vt,ds.vni,ds.ek,ds.dq,ds.cg,ds.cl,ds.fe,ds.eh,ds.wc,ds.zt,ds.ci,ds.vi,N,Nc,Nv,nl),0)
            BENCH("aos",DT(ka,c.tx,c.ty,gr,bk,da.out,da.e3,da.vni,da.ek,da.fe,da.e2,da.cn,da.cd,da.zt,N,Nc,Nv,nl),0)
            BENCH("grp",DT(kg,c.tx,c.ty,gr,bk,dg.out,dg.gA,dg.gD,dg.ek,dg.dq,dg.fe,dg.cg,dg.cl,dg.zt,dg.ci,dg.vi,N,Nc,Nv,nl),0)
            BENCH("aosoa",DTV(ko,32,c.tx,c.ty,gr,bk,d32.out,d32.gA,d32.gD,d32.ek,d32.dq,d32.fe,d32.cg,d32.cl,d32.zt,d32.ci,d32.vi,N,Nc,Nv,nl),32)
            BENCH("aosoa",DTV(ko,64,c.tx,c.ty,gr,bk,d64.out,d64.gA,d64.gD,d64.ek,d64.dq,d64.fe,d64.cg,d64.cl,d64.zt,d64.ci,d64.vi,N,Nc,Nv,nl),64)
            fflush(csv);printf(" done\n");}

        fds(ds);fda(da);fdg(dg);fdo(d32);fdo(d64);soa.free_all();fao.free_all();grp.free_all();a32.free_all();a64.free_all();delete[]hr;}
    if (icon_ci) { delete[] icon_ci; delete[] icon_vi; }
    GPU_CHECK(GED(e0));GPU_CHECK(GED(e1));fd();fclose(csv);
    printf("\nWritten: ddt_vn_apc_pc_gpu_sweep.csv  (%d×%d×5×%d)\n",n_dists,NC_,NRUNS);return 0;}

    main()