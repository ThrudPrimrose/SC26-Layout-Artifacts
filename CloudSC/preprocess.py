#!/usr/bin/env python3
"""
Inline ECMWF thermodynamic helper functions in cloudsc_py.

All function bodies are hardcoded below. The inlining pass runs twice
so that nested calls (e.g. foedem→foealfa) are fully resolved.

Usage:
    python inline_thermo.py input.py -o output.py
    python inline_thermo.py input.py              # stdout
"""

import argparse
import ast
import copy
import sys

# ============================================================================
# Hardcoded function table: name -> (param_names, return_expression_source)
#
# These are the closure-style helpers (single arg ptare, capturing ydcst/ydthf).
# The return expression is exactly what appears after "return" in the original.
# ============================================================================

FUNCTIONS: dict[str, tuple[list[str], str]] = {
    "foedelta": (
        ["ptare"],
        "max(0.0, 1.0 * np.sign(ptare - ydcst.rtt))",
    ),
    "foealfa": (
        ["ptare"],
        "min(1.0, ((max(ydthf.rtice, min(ydthf.rtwat, ptare)) - ydthf.rtice) * ydthf.rtwat_rtice_r) ** 2)",
    ),
    "foeewm": (
        ["ptare"],
        "ydthf.r2es * (foealfa(ptare) * np.exp(ydthf.r3les * (ptare - ydcst.rtt) / (ptare - ydthf.r4les))"
        " + (1.0 - foealfa(ptare)) * np.exp(ydthf.r3ies * (ptare - ydcst.rtt) / (ptare - ydthf.r4ies)))",
    ),
    "foedem": (
        ["ptare"],
        "foealfa(ptare) * ydthf.r5alvcp * (1.0 / (ptare - ydthf.r4les) ** 2)"
        " + (1.0 - foealfa(ptare)) * ydthf.r5alscp * (1.0 / (ptare - ydthf.r4ies) ** 2)",
    ),
    "foeldcpm": (
        ["ptare"],
        "foealfa(ptare) * ydthf.ralvdcp + (1.0 - foealfa(ptare)) * ydthf.ralsdcp",
    ),
    "foelhm": (
        ["ptare"],
        "foealfa(ptare) * ydcst.rlvtt + (1.0 - foealfa(ptare)) * ydcst.rlstt",
    ),
    "foeeliq": (
        ["ptare"],
        "ydthf.r2es * np.exp(ydthf.r3les * (ptare - ydcst.rtt) / (ptare - ydthf.r4les))",
    ),
    "foeeice": (
        ["ptare"],
        "ydthf.r2es * np.exp(ydthf.r3ies * (ptare - ydcst.rtt) / (ptare - ydthf.r4ies))",
    ),
    "fokoop": (
        ["ptare"],
        "min(ydthf.rkoop1 - ydthf.rkoop2 * ptare,"
        " foeeliq(ptare) / foeeice(ptare))",
    ),
    "foedelta": (
        ["ptare"],
        "max(0.0, 1.0 * np.sign(ptare - ydcst.rtt))",
    ),
    "foeldcp": (
        ["ptare"],
        "foedelta(ptare) * ydthf.ralvdcp + (1.0 - foedelta(ptare)) * ydthf.ralsdcp",
    ),
    "foelh": (
        ["ptare"],
        "foedelta(ptare) * ydcst.rlvtt + (1.0 - foedelta(ptare)) * ydcst.rlstt",
    ),
    "foede": (
        ["ptare"],
        "(foedelta(ptare) * ydthf.r5alvcp + (1.0 - foedelta(ptare)) * ydthf.r5alscp)"
        " / (ptare - (ydthf.r4les * foedelta(ptare) + ydthf.r4ies * (1.0 - foedelta(ptare)))) ** 2",
    ),
    "foeew": (
        ["ptare"],
        "ydthf.r2es * np.exp((ydthf.r3les * foedelta(ptare) + ydthf.r3ies * (1.0 - foedelta(ptare)))"
        " * (ptare - ydcst.rtt) / (ptare - (ydthf.r4les * foedelta(ptare) + ydthf.r4ies * (1.0 - foedelta(ptare)))))",
    ),
}

# ============================================================================
# Parse each return expression into an AST node (done once at import time)
# ============================================================================

_PARSED: dict[str, tuple[list[str], ast.expr]] = {}
for _name, (_params, _src) in FUNCTIONS.items():
    _expr_node = ast.parse(_src, mode="eval").body
    _PARSED[_name] = (_params, _expr_node)


# ============================================================================
# AST transformer: substitute Name nodes from a mapping
# ============================================================================

class _Substitutor(ast.NodeTransformer):
    def __init__(self, mapping: dict[str, ast.expr]):
        self.mapping = mapping

    def visit_Name(self, node: ast.Name) -> ast.expr:
        if node.id in self.mapping:
            return copy.deepcopy(self.mapping[node.id])
        return node


# ============================================================================
# AST transformer: replace call sites with inlined expressions
# ============================================================================

class _Inliner(ast.NodeTransformer):
    def __init__(self):
        self.count: dict[str, int] = {n: 0 for n in _PARSED}

    def visit_Call(self, node: ast.Call) -> ast.expr:
        # Recurse first (inside-out) so inner calls are inlined before outer
        self.generic_visit(node)

        if not isinstance(node.func, ast.Name):
            return node
        name = node.func.id
        if name not in _PARSED:
            return node
        if node.keywords:
            return node

        params, expr_template = _PARSED[name]
        if len(node.args) != len(params):
            print(
                f"  [warn] {name}: expected {len(params)} args, got {len(node.args)}",
                file=sys.stderr,
            )
            return node

        mapping = dict(zip(params, node.args))
        result = _Substitutor(mapping).visit(copy.deepcopy(expr_template))
        self.count[name] += 1
        return ast.copy_location(result, node)


# ============================================================================
# Also remove original function definitions if requested
# ============================================================================

class _RemoveDefs(ast.NodeTransformer):
    def visit_FunctionDef(self, node: ast.FunctionDef):
        if node.name in _PARSED:
            return None
        self.generic_visit(node)
        return node


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Inline ECMWF thermodynamic helpers in cloudsc_py"
    )
    parser.add_argument("input", help="Input Python source file", default="cloudsc_py.py")
    parser.add_argument("-o", "--output", default="cloudsc_preprocessed.py", help="Output file (default: stdout)")
    parser.add_argument(
        "--remove-defs", action="store_true",
        help="Remove the original function definitions after inlining",
    )
    parser.add_argument(
        "--passes", type=int, default=2,
        help="Number of inlining passes (default: 2)",
    )
    args = parser.parse_args()

    with open(args.input) as f:
        source = f.read()

    tree = ast.parse(source)

    # Run multiple passes so nested calls get fully resolved
    total: dict[str, int] = {n: 0 for n in _PARSED}
    for i in range(args.passes):
        inliner = _Inliner()
        tree = inliner.visit(tree)
        ast.fix_missing_locations(tree)
        any_work = False
        for name, n in inliner.count.items():
            total[name] += n
            if n > 0:
                any_work = True
        if not any_work:
            print(f"  pass {i + 1}: nothing left to inline, stopping early", file=sys.stderr)
            break
        else:
            changed = ", ".join(f"{n}:{c}" for n, c in inliner.count.items() if c > 0)
            print(f"  pass {i + 1}: {changed}", file=sys.stderr)

    # Summary
    for name, n in total.items():
        if n > 0:
            print(f"  {name}: {n} total", file=sys.stderr)

    if args.remove_defs:
        tree = _RemoveDefs().visit(tree)
        ast.fix_missing_locations(tree)

    result = ast.unparse(tree)

    if args.output:
        with open(args.output, "w") as f:
            f.write(result)
        print(f"Written to {args.output}", file=sys.stderr)
    else:
        print(result)


if __name__ == "__main__":
    main()