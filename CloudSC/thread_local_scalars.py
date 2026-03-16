#!/usr/bin/env python3
"""
rename_threadprivate_zscalars.py

Usage: python rename_threadprivate_zscalars.py cloudsc_py.py > cloudsc_py_renamed.py

Finds all z* scalar variables (not arrays, not array-element assignments) that are
assigned within each top-level loop nest, and renames them with a _<nest_id> suffix
to make their thread-private nature explicit.

Rules:
  - Array variables (allocated via np.ndarray) are NEVER renamed.
  - "Global" scalars assigned before the first loop nest are NEVER renamed.
  - Only z* scalars assigned (bare, not indexed) within a loop nest are renamed.
  - The suffix is _<N> where N is the 0-based index of the top-level loop nest.
"""

import re
import sys


def run(src: str) -> str:
    lines = src.split('\n')

    # --- 1. Collect array names ---
    arrays = set()
    for line in lines:
        m = re.match(r'\s+(\w+)\s*=\s*np\.ndarray\(', line)
        if m:
            arrays.add(m.group(1))

    # --- 2. Find body indent (first np.ndarray allocation or first non-signature line) ---
    body_indent = None
    body_start = None
    past_signature = False
    paren_depth = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('@dace.program'):
            continue
        if stripped.startswith('def '):
            paren_depth += stripped.count('(') - stripped.count(')')
            past_signature = paren_depth <= 0
            continue
        if not past_signature:
            paren_depth += stripped.count('(') - stripped.count(')')
            if paren_depth <= 0 and stripped.endswith(':'):
                past_signature = True
            continue
        if not stripped:
            continue
        body_indent = len(line) - len(line.lstrip())
        body_start = i
        break

    if body_start is None:
        return src

    # --- 3. Find top-level for-loop nests ---
    nests = []  # list of (start, end) line indices
    i = body_start
    while i < len(lines):
        stripped = lines[i].strip()
        if not stripped:
            i += 1
            continue
        indent = len(lines[i]) - len(lines[i].lstrip())
        if indent < body_indent:
            break
        if indent == body_indent and stripped.startswith('for '):
            ns = i
            i += 1
            while i < len(lines):
                s = lines[i].strip()
                if not s:
                    i += 1
                    continue
                if len(lines[i]) - len(lines[i].lstrip()) <= body_indent:
                    break
                i += 1
            nests.append((ns, i))
        else:
            i += 1

    # --- 4. Collect global scalars (assigned before first nest) ---
    global_scalars = set()
    first_nest = nests[0][0] if nests else len(lines)
    for i in range(body_start, first_nest):
        stripped = lines[i].strip()
        # Bare scalar assignment: zfoo = expr  (not zfoo[...] = ..., not np.ndarray)
        m = re.match(r'(z\w+)\s*=\s*(?!\s*np\.ndarray)', stripped)
        if m and not re.match(r'z\w+\s*\[', stripped):
            global_scalars.add(m.group(1))

    print(f"Arrays:          {len(arrays)}", file=sys.stderr)
    print(f"Global scalars:  {len(global_scalars)}: {sorted(global_scalars)}", file=sys.stderr)
    print(f"Loop nests:      {len(nests)}", file=sys.stderr)

    # --- 5. For each nest, identify & rename thread-private scalars ---
    # Process in reverse to preserve line indices
    for idx in range(len(nests) - 1, -1, -1):
        ns, ne = nests[idx]
        # Find all z* scalars assigned (not indexed) in this nest
        nest_scalars = set()
        for li in range(ns, ne):
            stripped = lines[li].strip()
            m = re.match(r'(z\w+)\s*=\s*', stripped)
            if m and not re.match(r'z\w+\s*\[', stripped):
                name = m.group(1)
                if name not in arrays and name not in global_scalars:
                    nest_scalars.add(name)

        if not nest_scalars:
            continue

        suffix = f'_{idx}'
        # Sort longest-first for safe replacement
        by_len = sorted(nest_scalars, key=len, reverse=True)
        # Build single combined pattern for speed
        pattern = re.compile(r'\b(' + '|'.join(re.escape(n) for n in by_len) + r')\b')

        print(f"Nest {idx:2d} (L{ns+1:4d}-{ne:4d}): {len(nest_scalars):2d} scalars {sorted(nest_scalars)}", file=sys.stderr)

        for li in range(ns, ne):
            lines[li] = pattern.sub(lambda m: m.group(0) + suffix, lines[li])

    return '\n'.join(lines)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python rename_threadprivate_zscalars.py INPUT.py [OUTPUT.py]", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        src = f.read()

    result = run(src)

    if len(sys.argv) >= 3:
        with open(sys.argv[2], 'w') as f:
            f.write(result)
        print(f"Written to {sys.argv[2]}", file=sys.stderr)
    else:
        print(result)