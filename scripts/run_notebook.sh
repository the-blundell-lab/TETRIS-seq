#!/bin/bash
###############################################################################
# run_notebook.sh
#
# Runs one of the repository's notebooks headlessly, against data held
# somewhere else.
#
#   run_notebook.sh <notebook> [data_directory] [output_notebook]
#
#   <notebook>          path to the .ipynb, or just its name - the mCA_caller/
#                       and noise_correction_model/ folders are searched
#   [data_directory]    working directory the notebook runs in (default: .)
#   [output_notebook]   where to write the executed copy
#                       (default: <data_directory>/<name>_run.ipynb)
#
# The notebooks use relative paths from the working directory ("base_path = .").
# nbconvert runs a notebook in the notebook's own directory, which would mean
# copying it next to your data; papermill lets the working directory be chosen,
# so the repository can stay where it is.
#
# Needs the analysis environment:  conda activate tetris-seq-analysis
###############################################################################

set -euo pipefail

P="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$P/.." && pwd)"

[ $# -ge 1 ] || { sed -n '3,22p' "$0"; exit 1; }

NB="$1"
DATA="$(cd "${2:-.}" && pwd)"
if [ ! -f "$NB" ]; then
    for d in "$REPO/mCA_caller" "$REPO/noise_correction_model" "$REPO/chromosomal_rearrangement_caller"; do
        [ -f "$d/$NB" ] && { NB="$d/$NB"; break; }
    done
fi
[ -f "$NB" ] || { echo "Error: notebook not found: $1"; exit 1; }
NB="$(cd "$(dirname "$NB")" && pwd)/$(basename "$NB")"

OUT="${3:-$DATA/$(basename "${NB%.ipynb}")_run.ipynb}"

command -v papermill >/dev/null 2>&1 || {
    echo "Error: papermill not found. Activate the analysis environment:"
    echo "    conda activate tetris-seq-analysis"
    echo "(if it predates papermill:  conda env update -f $REPO/environment_analysis.yml)"
    exit 1
}

echo "notebook   $NB"
echo "data       $DATA"
echo "output     $OUT"
echo

papermill "$NB" "$OUT" --cwd "$DATA" --log-output
