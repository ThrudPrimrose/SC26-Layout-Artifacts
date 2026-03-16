# This script detects if it's being sourced in fish and redirects if so.
# We use a syntax that is mostly compatible across bash, zsh, and fish for the first line.
test -n "$fish_version" && source env.fish && return 0 2>/dev/null

# --- Standard bash/zsh logic below ---
VENV_DIR=".venv"

if [ ! -d "$VENV_DIR" ]; then
    echo "Error: Virtual environment directory '$VENV_DIR' not found."
    echo "Please run: python3 -m venv .venv"
    return 1 2>/dev/null || exit 1
fi

source "$VENV_DIR/bin/activate"
export PYTHONPATH=$PYTHONPATH:$(pwd)

echo "Environment activated (bash/zsh)."
echo "Python location: $(which python)"
echo "DaCe version: $(python -c 'import dace; print(dace.__version__)' 2>/dev/null || echo 'Not installed')"
