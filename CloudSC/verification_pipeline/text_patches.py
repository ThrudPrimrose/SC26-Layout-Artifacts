"""text_patches.py — Post-codegen text-level patches for CloudSC.

Ported from velocity/utils/compile_if_propagated_sdfgs.py.
Applied to generated C++ source after DaCe codegen, before compilation.
"""

from __future__ import annotations

import re
from pathlib import Path


def fix_mixed_precision_ambiguity(file_path: Path):
    """Inject fp16_operators.h to resolve double/half operator ambiguity.

    Ported from velocity/utils/compile_if_propagated_sdfgs.py.
    """
    with open(file_path, "r") as f:
        content = f.read()

    include = '#include "fp16_operators.h"'
    if include in content:
        return

    new_content = re.sub(
        r"(#include\s+[<\"][^>\"]+[>\"])",
        rf"\1\n{include}",
        content,
        count=1,
    )

    if new_content != content:
        with open(file_path, "w") as f:
            f.write(new_content)


def floatify_generated_code(file_path: Path):
    """Replace double literals in generated C++ with float equivalents.

    DaCe codegen emits double literals for interstate edge conditions
    and other places that the SDFG-level floatify doesn't reach.
    Ported from compile_if_propagated_sdfgs.py.
    """
    with open(file_path, "r") as f:
        code = f.read()

    new_code = code
    for val in ["0.85", "0.5", "1.0", "0.0", "0.05", "0.65", "1.15"]:
        new_code = re.sub(
            rf"(?<!float\()(?<![0-9]){re.escape(val)}(?!f)(?![0-9])",
            f"{val}f",
            new_code,
        )

    # ipow(x, 2) -> (x * x) — safety net for anything codegen emitted
    new_code = re.sub(
        r"dace::math::ipow\(([^,]+),\s*2\)",
        r"((\1) * (\1))",
        new_code,
    )

    if new_code != code:
        with open(file_path, "w") as f:
            f.write(new_code)


# TODO: BFP text patching
# def patch_bfp_reads(code, bfp_gpu_names, block_size=32, mantissa_bits=16):
#     """Text-level BFP patching: fix parameter types and replace array reads
#     with bfp_decode calls.  Ported from velocity."""
#     ...


def inject_copy_cast(file_path: Path):
    """Inject copy_cast.h after dace/copy.h to handle mixed-type CopyND calls."""
    with open(file_path, "r") as f:
        content = f.read()

    include = '#include "copy_cast.h"'
    if include in content:
        return

    # Insert after the dace include
    new_content = re.sub(
        r'(#include\s+<dace/dace\.h>)',
        rf'\1\n{include}',
        content,
        count=1,
    )

    if new_content != content:
        with open(file_path, "w") as f:
            f.write(new_content)


def apply_text_patches(codegen_dir: Path, lowprec: str):
    """Apply all text-level patches to generated source files.

    Args:
        codegen_dir: Path to the codegen/ directory with generated .cpp/.h files.
        lowprec: Precision mode string.
    """
    if lowprec in ("fp64", "f64"):
        return  # nothing to patch

    source_files = list(codegen_dir.glob("*.cpp")) + list(codegen_dir.glob("*.h"))

    for f in source_files:
        if lowprec in ("fp16", "f16"):
            fix_mixed_precision_ambiguity(f)
        inject_copy_cast(f)
        floatify_generated_code(f)

    print(f"Applied text patches to {len(source_files)} files in {codegen_dir}")
