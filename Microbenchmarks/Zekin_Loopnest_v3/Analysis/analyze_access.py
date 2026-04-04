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
        """
        Single collapsed signature for the whole nest: one entry per dimension role.
        For each role (h, v, b) present across ALL arrays, take the worst case:
          U > S > C
        If mixed S/U for a role → U wins.

        Example: arrays with patterns {h:S v:S b:S, h:U v:S h:U}
          → h: U (some array is U)  v: S (all structured)  b: S
          → collapsed: "h:U v:S b:S"
        """
        # Collect worst-case structuredness per role
        role_worst = {}  # role -> worst S/U seen
        priority = {"U": 3, "S": 2, "C": 1}

        for acc in self.accesses:
            for v, su in zip(acc.dim_vars, acc.structuredness):
                role = _role(v) if v and v not in ("const", "?") else None
                if role is None:
                    continue  # skip constant dims
                cur = role_worst.get(role, "C")
                if priority.get(su, 0) > priority.get(cur, 0):
                    role_worst[role] = su

        # Build sorted signature: h, v, b order
        order = {"h": 0, "v": 1, "b": 2}
        parts = []
        for role in sorted(role_worst.keys(), key=lambda r: order.get(r, 99)):
            parts.append(f"{role}:{role_worst[role]}")
        return "  ".join(parts) if parts else "-:C"

    @property
    def pattern_key(self):
        """Collapsed identity: shape + ranges + behavior + collapsed S/U."""
        return (self.loop_shape, tuple(self.ranges), self.behavior, self.collapsed_su)

    @property
    def shape_key(self):
        """Shape + ranges + behavior (ignoring S/U details)."""
        return (self.loop_shape, tuple(self.ranges), self.behavior)

    @property
    def source_code(self):
        """
        Fortran source for this nest via Loki's fgen backend.

        For b.h nests (jb outer, jc/je inner): prints only the inner loop,
        since the jb loop body contains many unrelated nests.
        For v.h nests (jk outer, jc/je inner): prints the full outer loop
        (jk wraps only this one inner loop, so it's compact).
        For single loops: prints the loop.
        """
        if not self.ir_nodes:
            return "! (no IR nodes stored)"
        if len(self.ir_nodes) == 1:
            return fgen(self.ir_nodes[0])
        # 2-level nest
        outer_role = _role(self.vars[0])
        if outer_role == "b":
            # jb body is large — show only the inner loop
            return f"! inside DO {self.vars[0]} = {self.bounds[0]}\n" + fgen(self.ir_nodes[1])
        else:
            # v.h or other — show the full outer loop (compact)
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
# Behavior detection: compute / accumulate / reduction
# ---------------------------------------------------------------------------

def detect_behavior(loop_body):
    """
    Detect loop behavior from assignments in the body.

    Returns:
      "reduction"   — scalar LHS accumulated with MAX/MIN (classical reduction)
      "accumulate"  — array LHS appears on RHS too (update-in-place, e.g. A = A + ...)
      "compute"     — pure computation (LHS not on RHS)
    """
    has_accumulate = False

    for assign in FindNodes(Assignment).visit(loop_body):
        lhs_name = get_full_name(assign.lhs).lower()
        rhs_str = str(assign.rhs).lower()

        # Check for MAX/MIN reduction on scalar (classical reduction pattern)
        if "max(" in rhs_str or "min(" in rhs_str:
            if lhs_name in rhs_str:
                return "reduction"

        # Check for accumulation: LHS name appears on RHS
        if lhs_name in rhs_str:
            has_accumulate = True

    return "accumulate" if has_accumulate else "compute"


# ---------------------------------------------------------------------------
# Loop nest discovery — captures 1-level AND 2-level nests
# ---------------------------------------------------------------------------

def _direct_child_loops(node):
    """
    Find Loop nodes that are DIRECT children of 'node' — not nested through
    another Loop.  Recurses through Sections, Conditionals, etc. but STOPS
    at Loop boundaries.

    This prevents jb -> je being captured when je actually lives inside jk.
    """
    result = []
    # node.body / node.else_body are the children to walk
    for attr in ("body", "else_body"):
        children = getattr(node, attr, None)
        if children is None:
            continue
        if not isinstance(children, (list, tuple)):
            children = [children]
        for child in children:
            if isinstance(child, Loop):
                result.append(child)
                # Do NOT recurse into this loop — its children are not direct
            elif hasattr(child, "body") or hasattr(child, "else_body"):
                # Recurse into Sections, Conditionals, PragmaRegions, etc.
                result.extend(_direct_child_loops(child))
    return result


def analyze_loop_nests(routine):
    """
    Walk routine body and find all leaf loop nests using direct-child discovery.

    For each loop, finds its DIRECT child loops (not grandchildren through
    intermediate loops).  Captures:
      - 2-level pairs: parent -> child (where child has no further children)
      - 1-level: leaf loops with no child loops at all
    """
    nests = []
    all_loops = FindNodes(Loop).visit(routine.body)

    # Track loops captured as inner of a pair
    captured_as_inner = set()

    for outer_loop in all_loops:
        outer_var = str(outer_loop.variable)

        # Find DIRECT child loops only (not grandchildren)
        direct_children = _direct_child_loops(outer_loop)
        if not direct_children:
            continue  # leaf — handled in second pass

        for inner_loop in direct_children:
            inner_var = str(inner_loop.variable)

            # Only leaf pairs: inner loop must have no children
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

    # Second pass: capture leaf loops NOT already used as inner of a pair
    for loop in all_loops:
        if id(loop) in captured_as_inner:
            continue
        if _direct_child_loops(loop):
            continue  # not a leaf

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
# Grouping
# ---------------------------------------------------------------------------

def group_arrays(nests):
    """Group arrays by observed role:S/U signatures."""
    array_role_sigs = defaultdict(set)
    for nest in nests:
        for acc in nest.accesses:
            array_role_sigs[acc.name].add(acc.role_su_signature)

    role_groups = defaultdict(set)
    for name, sigs in array_role_sigs.items():
        role_groups[frozenset(sigs)].add(name)

    return array_role_sigs, role_groups


# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

def print_analysis(routine_name, nests, array_role_sigs, role_groups):
    p = print  # shorthand

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
# Single-block view (b=1): strip block dimension entirely
# ---------------------------------------------------------------------------

def _strip_block(role_su_sig):
    """Remove b:* entries from a role:S/U signature string."""
    parts = [p.strip() for p in role_su_sig.split("  ") if p.strip()]
    filtered = [p for p in parts if not p.startswith("b:")]
    return "  ".join(filtered) if filtered else "-"


def _strip_block_collapsed(collapsed):
    """Remove b:* from collapsed signature."""
    parts = [p.strip() for p in collapsed.split("  ") if p.strip()]
    filtered = [p for p in parts if not p.startswith("b:")]
    return "  ".join(filtered) if filtered else "-"


def _strip_block_shape(shape):
    """Remove 'b' from loop shape: 'b.h' -> 'h', 'b.v' -> 'v', 'v.h' stays."""
    parts = [p for p in shape.split(".") if p != "b"]
    return ".".join(parts) if parts else "-"


def print_single_block_view(routine_name, nests):
    """Reprint assuming nblks=1: strip block dimension entirely."""
    p = print

    p(f"\n---\n")
    p(f"# Single-Block View (`nblks=1`)\n")
    p(f"Block dimension stripped: `b.h` becomes `h`, `b:S`/`b:U` removed from all signatures.\n")

    # --- Build single-block pattern groups ---
    pattern_groups_nb = defaultdict(list)
    for nest in nests:
        shape_nb = _strip_block_shape(nest.loop_shape)
        ranges_nb = tuple(r for r, v in zip(nest.ranges, nest.vars) if _role(v) != "b")
        collapsed_nb = _strip_block_collapsed(nest.collapsed_su)
        key = (shape_nb, ranges_nb, nest.behavior, collapsed_nb)
        pattern_groups_nb[key].append(nest)

    # --- Loop patterns table ---
    p("## Loop Patterns (single-block, collapsed)\n")
    p("| Count | Shape | Ranges | Behavior | Collapsed S/U |")
    p("|------:|-------|--------|----------|---------------|")
    for key, group in sorted(pattern_groups_nb.items(), key=lambda x: -len(x[1])):
        shape_nb, ranges_nb, behavior, collapsed_nb = key
        ranges_str = ", ".join(ranges_nb) if ranges_nb else "-"
        beh = behavior if behavior != "compute" else ""
        p(f"| {len(group)} | `{shape_nb}` | {ranges_str} | {beh} | `{collapsed_nb}` |")
    p()

    # --- Example nest per pattern ---
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

    # --- Widest nest ---
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

    # --- Array groups without block ---
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

    # --- Per-array table ---
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

    # Write to markdown file
    outpath = filepath.with_suffix(".analysis.md")
    import io
    buf = io.StringIO()
    _orig_print = __builtins__["print"] if isinstance(__builtins__, dict) else __builtins__.print

    def md_print(*args, **kwargs):
        kwargs["file"] = buf
        _orig_print(*args, **kwargs)

    # Monkey-patch print for the report functions
    import builtins
    old_print = builtins.print
    builtins.print = md_print

    print_analysis(routine.name, nests, array_role_sigs, role_groups)
    print_single_block_view(routine.name, nests)

    builtins.print = old_print

    outpath.write_text(buf.getvalue())
    old_print(f"Report written to {outpath}")

    # Also print to stdout
    old_print(buf.getvalue())


if __name__ == "__main__":
    main()