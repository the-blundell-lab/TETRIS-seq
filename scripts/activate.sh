###############################################################################
# activate.sh
#
# Activates both environments the pipeline needs, in the right order:
#   1. the conda environment, for the command-line tools
#   2. the virtualenv, for the Python the scripts run under
#
# SOURCE it, do not run it - it has to change the current shell:
#
#   source scripts/activate.sh
#
# The virtualenv is looked for next to the repository (../tetris-py), then in
# the repository itself. Override with TETRIS_PY, or pass a path:
#
#   source scripts/activate.sh /path/to/tetris-py
#
# Why two environments: conda-forge's Python has no gdbm extension, so `shelve`
# falls back to dbm.dumb and consensus calling never finishes. See
# docs/INSTALL.md section 1b.
###############################################################################

# --- must be sourced ---------------------------------------------------------
if [ -n "${BASH_SOURCE[0]:-}" ] && [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "activate.sh changes the current shell, so it must be sourced:"
    echo "    source ${0}"
    exit 1
fi

_tetris_repo="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )/.." && pwd )"
_tetris_env_name="${TETRIS_CONDA_ENV:-tetris-seq-pipeline}"

# --- 1. conda ----------------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
    for _c in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/mambaforge" \
              "$HOME/anaconda3" /opt/conda /opt/miniforge3; do
        if [ -f "$_c/etc/profile.d/conda.sh" ]; then
            . "$_c/etc/profile.d/conda.sh"
            break
        fi
    done
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "conda not found. Install Miniforge (see docs/INSTALL.md) or source its"
    echo "profile script first:  . /path/to/miniforge3/etc/profile.d/conda.sh"
    return 1
fi

if ! conda env list | awk '{print $1}' | grep -qx "$_tetris_env_name"; then
    echo "conda environment '$_tetris_env_name' not found. Create it with:"
    echo "    conda env create -f $_tetris_repo/environment_sequencing.yml"
    return 1
fi

conda activate "$_tetris_env_name" || return 1

# --- 2. virtualenv -----------------------------------------------------------
_tetris_py="${1:-${TETRIS_PY:-}}"
if [ -z "$_tetris_py" ]; then
    for _p in "$_tetris_repo/../tetris-py" "$_tetris_repo/tetris-py"; do
        [ -x "$_p/bin/activate" ] || [ -f "$_p/bin/activate" ] && { _tetris_py="$_p"; break; }
    done
fi

if [ -z "$_tetris_py" ] || [ ! -f "$_tetris_py/bin/activate" ]; then
    echo "Virtualenv not found (looked next to the repository for tetris-py/)."
    echo "Create it with:"
    echo "    $_tetris_repo/scripts/setup_python_env.sh"
    echo "or point at an existing one:"
    echo "    source scripts/activate.sh /path/to/tetris-py"
    return 1
fi

# shellcheck disable=SC1091
. "$_tetris_py/bin/activate" || return 1

# --- 3. make the repository easy to refer to ---------------------------------
# TETRIS_SEQ points at the repository, and scripts/ goes on PATH so the wrappers
# can be run by name from any working directory.
export TETRIS_SEQ="$_tetris_repo"

# The pipeline scripts source config/config.sh themselves, but anything run
# outside them - the demo command in demo/README.md, a script of your own -
# needs these too, so export them here as well.
# shellcheck disable=SC1091
[ -f "$_tetris_repo/config/config.sh" ] && . "$_tetris_repo/config/config.sh"
: "${EXTERNAL_TOOLS:=$HOME/Pipeline_tools}"
: "${REF:=$EXTERNAL_TOOLS/Homo_sapiens_assembly19.fasta}"
: "${ANNOVAR_HOME:=$EXTERNAL_TOOLS/annovar}"
: "${PINDEL_DIR:=$EXTERNAL_TOOLS/pindel}"
export EXTERNAL_TOOLS REF ANNOVAR_HOME PINDEL_DIR
case ":$PATH:" in
    *":$_tetris_repo/scripts:"*) ;;
    *) export PATH="$_tetris_repo/scripts:$PATH" ;;
esac

# --- 4. report ---------------------------------------------------------------
_tetris_backend=$(python - <<'EOF' 2>/dev/null
try:
    import dbm.gnu; print("dbm.gnu")
except ImportError:
    try:
        import dbm.ndbm; print("dbm.ndbm")
    except ImportError:
        print("dbm.dumb")
EOF
)

echo "python       $(command -v python)  ($(python --version 2>&1))"
echo "dbm backend  ${_tetris_backend:-unknown}"
echo "tools        $(dirname "$(command -v bwa 2>/dev/null || echo 'bwa-not-found/x')")"
echo "repository   $TETRIS_SEQ  (on PATH: run_sample.sh, check_setup.sh, ...)"
if [ -f "$REF" ]; then
    echo "reference    $REF"
else
    echo "reference    NOT FOUND at $REF"
    echo "             set EXTERNAL_TOOLS (or REF) in config/config.sh - see docs/INSTALL.md"
fi

if [ "$_tetris_backend" = "dbm.dumb" ]; then
    echo
    echo "WARNING: only dbm.dumb is available - consensus calling will not finish."
    echo "         See docs/INSTALL.md section 1b."
fi

unset _tetris_repo _tetris_env_name _tetris_py _tetris_backend _c _p
