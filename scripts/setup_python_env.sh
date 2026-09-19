#!/bin/bash
###############################################################################
# setup_python_env.sh
#
# Creates the Python environment the pipeline's own scripts run in.
#
#   scripts/setup_python_env.sh [destination]
#
# Default destination: ./tetris-py
#
# WHY THIS EXISTS
# ---------------
# The consensus callers store the UMI-grouped reads with `shelve`, which needs a
# real dbm backend (gdbm or ndbm). conda-forge's Python does not ship the gdbm
# extension, so under a conda environment Python silently falls back to
# `dbm.dumb` - a pure-Python store that rewrites its index on every sync. The
# pipeline then appears to run, at low CPU, and never finishes.
#
# So: conda provides the command-line tools (bwa, samtools, Picard, fgbio,
# VarDictJava, Pindel, Java 8) and this virtualenv provides the Python the
# scripts run under. Activate the conda environment first, then this one.
###############################################################################

set -euo pipefail

DEST="${1:-./tetris-py}"
PACKAGES="pysam pyfaidx numpy pandas scipy biopython matplotlib"

# ---- find a Python with a usable dbm backend --------------------------------
have_dbm() {
    "$1" - <<'EOF' >/dev/null 2>&1
try:
    import dbm.gnu
except ImportError:
    import dbm.ndbm
EOF
}

PYTHON=""
for candidate in /usr/bin/python3 /usr/local/bin/python3 python3; do
    cmd=$(command -v "$candidate" 2>/dev/null) || continue
    if have_dbm "$cmd"; then PYTHON="$cmd"; break; fi
done

if [ -z "$PYTHON" ]; then
    cat <<'EOF'
No Python with a working dbm backend (gdbm or ndbm) was found.

The pipeline needs one: `shelve` otherwise falls back to dbm.dumb and consensus
calling never completes. Options, in order of preference:

  * Install your distribution's gdbm bindings, then re-run this script:
        Debian/Ubuntu:  sudo apt install python3-gdbm
        RHEL/Fedora:    sudo dnf install python3-gdbm

  * Use any other Python for which `import dbm.gnu` succeeds:
        scripts/setup_python_env.sh /path/to/venv    # after putting it on PATH

Check a given interpreter with:
    <python> -c "import dbm.gnu"
EOF
    exit 1
fi

echo "python      $PYTHON ($("$PYTHON" --version 2>&1))"
echo "dbm backend $("$PYTHON" -c 'try:
 import dbm.gnu; print("dbm.gnu")
except ImportError:
 import dbm.ndbm; print("dbm.ndbm")')"

# ---- create the virtualenv --------------------------------------------------
if [ -x "$DEST/bin/python" ]; then
    echo "venv        $DEST (already exists)"
else
    echo "create      $DEST"
    if ! "$PYTHON" -m venv "$DEST" 2>/dev/null; then
        # Debian/Ubuntu without python3-venv: ensurepip is missing, so create the
        # environment without pip and bootstrap it.
        echo "            ensurepip unavailable - bootstrapping pip"
        "$PYTHON" -m venv --without-pip "$DEST"
        tmp=$(mktemp -d)
        curl -sSL https://bootstrap.pypa.io/get-pip.py -o "$tmp/get-pip.py"
        "$DEST/bin/python" "$tmp/get-pip.py" >/dev/null
        rm -rf "$tmp"
    fi
fi

# ---- install the packages ---------------------------------------------------
echo "install     $PACKAGES"
# shellcheck disable=SC2086
"$DEST/bin/pip" install --upgrade --quiet pip $PACKAGES

# ---- verify -----------------------------------------------------------------
"$DEST/bin/python" - <<'EOF'
import dbm, pysam, pyfaidx, numpy, pandas, scipy, Bio, matplotlib
try:
    import dbm.gnu as m
except ImportError:
    import dbm.ndbm as m
print("verified    dbm backend " + m.__name__ + ", pysam " + pysam.__version__)
EOF

cat <<EOF

Done. Use it on top of the conda environment, in this order:

    conda activate tetris-seq-pipeline     # tools: bwa, samtools, Picard, fgbio, ...
    source $DEST/bin/activate              # python the pipeline scripts run under

Check with: which python bwa
  python -> $DEST/bin/python
  bwa    -> the conda environment
EOF
