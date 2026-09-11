#!/usr/bin/env bash
# ---------------------------------------------------------------------
# Starts the Prosculpt Dashboard on Linux (e.g. a cluster headnode) -
# meant to run inside a detached `screen` session. Linux counterpart to
# start_dashboard.ps1 (Windows/Task Scheduler); same job, same log format.
#
# Launch (detached):
#   screen -dmS prosculpt_dashboard bash /path/to/start_dashboard.sh
# Reattach later (Ctrl-A D to detach again without killing it):
#   screen -r prosculpt_dashboard
# Stop:
#   screen -X -S prosculpt_dashboard quit
# ---------------------------------------------------------------------

set -u

# --- Adjust these for this machine ---
DASHBOARD_DIR="/home/folivieri/prosculpt_dev/dashboard"   # where app.py lives
PYTHON="python3"                                          # or an absolute path, e.g. a venv/conda env's python
PORT=""                                                   # leave empty to use dashboard_config.yaml's own `defaults: {port: ...}` (or 5000 if that's unset too) - set a value here only to force one regardless of that file, since a PORT env var always wins over it
MAX_WAIT_SECS=60                                          # how long to wait for DASHBOARD_DIR before giving up
LOG_FILE="$HOME/prosculpt_dashboard.log"
# --------------------------------------

log() {
  printf '%s  %s\n' "$(date -Is)" "$1" >>"$LOG_FILE"
}

# Mirrors start_dashboard.ps1's own wait-for-network-drive step - harmless
# here too if DASHBOARD_DIR is on a network mount (NFS/etc.) not yet ready
# at boot.
waited=0
while [ ! -d "$DASHBOARD_DIR" ] && [ "$waited" -lt "$MAX_WAIT_SECS" ]; do
  sleep 2
  waited=$((waited + 2))
done
if [ ! -d "$DASHBOARD_DIR" ]; then
  log "ERROR: $DASHBOARD_DIR never became available after ${MAX_WAIT_SECS}s - dashboard not started."
  exit 1
fi

cd "$DASHBOARD_DIR" || { log "ERROR: could not cd into $DASHBOARD_DIR"; exit 1; }

export PORT
export PYTHONUNBUFFERED=1  # otherwise output only flushes in delayed bursts once stdout isn't a real terminal

if [ -n "$PORT" ]; then
  log "Starting dashboard on port $PORT (cwd: $DASHBOARD_DIR)..."
else
  log "Starting dashboard (port from dashboard_config.yaml, or 5000 if that's unset too; cwd: $DASHBOARD_DIR)..."
fi
# tee, not just >>: also shows live output when you `screen -r` back in.
"$PYTHON" app.py 2>&1 | tee -a "$LOG_FILE"
# $? here would be tee's status (always 0), which hides how the dashboard
# actually ended. PIPESTATUS[0] is python's own, and a status above 128 means
# it was killed by signal (status - 128) rather than exiting on its own.
status=${PIPESTATUS[0]}
if [ "$status" -gt 128 ]; then
  log "Dashboard process was killed by signal $((status - 128))."
else
  log "Dashboard process exited (status $status)."
fi
exit "$status"
