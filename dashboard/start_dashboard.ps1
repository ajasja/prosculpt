# ---------------------------------------------------------------------
# Starts the Prosculpt Dashboard, meant to be launched automatically at
# Windows logon (see the Task Scheduler / Startup-folder setup in
# dashboard/README.md) - not something you'd normally run by hand, though
# it works fine that way too for testing.
#
# Waits for the mounted network drive to actually be reachable first: a
# persistent drive mapping isn't always reconnected the instant a logon
# session starts (the network stack can still be coming up), so launching
# immediately at logon can race that and fail with "path not found" even
# though the mapping is configured correctly.
#
# Runs with no visible console window when triggered from Task Scheduler/
# Startup (see README) - all output goes to $LogFile instead, since
# there's nowhere else for it to go once nothing is watching a terminal.
# ---------------------------------------------------------------------

# Deliberately NOT setting $ErrorActionPreference = "Stop" here: Windows
# PowerShell wraps every line a *native* process (python.exe) writes to
# stderr into an ErrorRecord once it's captured via `*>>` redirection
# below - completely harmless on its own (Flask's dev server always
# prints its "this is a development server" banner to stderr on startup),
# but with ErrorActionPreference set to Stop, that first ordinary banner
# line becomes a terminating error and kills this script before the
# server ever gets a chance to actually serve anything. The explicit
# `if (...) { ...; exit 1 }` checks below don't depend on this preference
# at all, so there's nothing this was actually protecting.

# --- Adjust these for this machine ---
$DashboardDir = "Q:\home\folivieri\prosculpt_dev\dashboard"   # where app.py lives on THIS machine's own drive mapping
$Python       = "python"                                       # or an absolute path, e.g. an isolated venv's python.exe
$Port         = "5000"
$DriveLetter  = "Q:"                                            # the mapped network drive $DashboardDir lives on
$MaxWaitSecs  = 60                                               # how long to wait for it to reconnect before giving up
$LogFile      = Join-Path $env:USERPROFILE "prosculpt_dashboard.log"
# --------------------------------------

function Log($msg) {
    Add-Content -Path $LogFile -Value "$(Get-Date -Format o)  $msg"
}

$waited = 0
while (-not (Test-Path $DriveLetter) -and $waited -lt $MaxWaitSecs) {
    Start-Sleep -Seconds 2
    $waited += 2
}
if (-not (Test-Path $DriveLetter)) {
    Log "ERROR: $DriveLetter never became available after $MaxWaitSecs s - dashboard not started."
    exit 1
}

if (-not (Test-Path $DashboardDir)) {
    Log "ERROR: $DashboardDir does not exist (drive is mounted, but this path under it isn't) - check `$DashboardDir above."
    exit 1
}

Set-Location $DashboardDir
$env:PORT = $Port

Log "Starting dashboard on port $Port (cwd: $DashboardDir)..."
# PYTHONUNBUFFERED=1: Python block-buffers stdout by default once it isn't
# attached to a real console (as it never is here, whether launched from
# Task Scheduler or run by hand) - without this, log lines arrive in
# delayed bursts instead of as they're actually printed, which matters
# for a log you're trying to read/tail live.
$env:PYTHONUNBUFFERED = "1"
# 2>&1 merges stderr into stdout so both streams get captured, piped into
# Out-File with an explicit -Encoding. Deliberately NOT `*>> $LogFile`
# (which looks equivalent and used to be what this line did) - Windows
# PowerShell's raw native-process redirection (`*>>`/`>>`) has a
# long-standing quirk where a native exe's UTF-8 output gets silently
# re-encoded as UTF-16 on the way into the file, inserting a null byte
# after every character. That reads back as a space between every single
# letter in anything that doesn't auto-detect UTF-16 (Notepad does;
# `type`/`cat`-style viewers and many simpler editors don't) - piping
# through Out-File with an explicit -Encoding avoids that redirection
# path entirely instead of fighting it.
& $Python app.py 2>&1 | Out-File -FilePath $LogFile -Append -Encoding utf8
Log "Dashboard process exited."
