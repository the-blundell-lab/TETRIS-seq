###############################################################################
# require_python_backend.sh
#
# Sourced by the wrapper scripts before they start work. Refuses to run if the
# Python on PATH has no usable dbm backend.
#
# The consensus callers store UMI-grouped reads with shelve. Without gdbm or
# ndbm, Python silently uses dbm.dumb, which rewrites its index on every sync:
# the run appears to proceed, at low CPU, and never finishes. Failing here costs
# a second; not failing here costs hours.
###############################################################################

_backend=$( { command -v python >/dev/null 2>&1 && python - <<'EOF'
try:
    import dbm.gnu; print("dbm.gnu")
except ImportError:
    try:
        import dbm.ndbm; print("dbm.ndbm")
    except ImportError:
        print("dbm.dumb")
EOF
} 2>/dev/null )

case "$_backend" in
    dbm.gnu|dbm.ndbm)
        ;;
    "")
        echo "Error: no 'python' on PATH."
        echo "       Activate the environments first:  source scripts/activate.sh"
        exit 1
        ;;
    *)
        echo "Error: this Python has no usable dbm backend (only dbm.dumb)."
        echo
        echo "       python:  $(command -v python)"
        echo
        echo "       Consensus calling shelves millions of reads and would never"
        echo "       finish. conda-forge's Python does not ship gdbm - use the"
        echo "       virtualenv built by scripts/setup_python_env.sh:"
        echo
        echo "           source scripts/activate.sh"
        echo
        echo "       See docs/INSTALL.md section 1b."
        exit 1
        ;;
esac
unset _backend
