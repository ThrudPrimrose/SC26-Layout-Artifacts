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


if __name__ == "__main__":
    main()