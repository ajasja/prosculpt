"""
Prosculpt Dashboard - Flask backend.

Run with:
    python app.py
then open http://localhost:5000 in a browser.

The server needs read access to the log file and to the pipeline's output
directory - run it on the machine (e.g. the cluster login node) that can
see that filesystem, and access it from your laptop via an SSH tunnel if
needed, e.g.:
    ssh -L 5000:localhost:5000 you@cluster
"""

from __future__ import annotations

import csv
import os

from flask import Flask, jsonify, request, send_file, abort

import parser as P

app = Flask(__name__, static_folder="static", template_folder="templates")

# Simple in-memory cache of the last resolved output_dir per log path, so
# helper endpoints (pdb/json fetchers) don't need to re-parse the whole log
# on every click.
_output_dir_cache: dict[str, str] = {}


def _get_output_dir(log_path: str) -> str:
    if log_path in _output_dir_cache:
        return _output_dir_cache[log_path]
    loc = P.resolve_output_dir(log_path)
    if not loc.get("output_dir"):
        abort(400, description=loc.get("error") or "Could not resolve output_dir")
    _output_dir_cache[log_path] = loc["output_dir"]
    return loc["output_dir"]


def _safe_join(base: str, *parts: str) -> str:
    """Join and ensure the result stays inside `base` (avoid path traversal
    from the file-browser / pdb-fetch endpoints)."""
    full = os.path.normpath(os.path.join(base, *parts))
    base_norm = os.path.normpath(base)
    if not (full == base_norm or full.startswith(base_norm + os.sep)):
        abort(400, description="Invalid path")
    return full


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------


@app.route("/")
def index():
    from flask import render_template

    return render_template("index.html")


# ---------------------------------------------------------------------------
# Status / overview
# ---------------------------------------------------------------------------


@app.route("/api/status")
def api_status():
    log_path = request.args.get("log")
    if not log_path:
        return jsonify({"error": "Missing ?log= parameter"}), 400
    status = P.get_status(log_path)
    if status.get("location", {}).get("output_dir"):
        _output_dir_cache[log_path] = status["location"]["output_dir"]
    return jsonify(status)


# ---------------------------------------------------------------------------
# Backbones (rfdiffusion)
# ---------------------------------------------------------------------------


@app.route("/api/backbones")
def api_backbones():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    return jsonify(P.list_backbones(output_dir))


@app.route("/api/backbone_pdb")
def api_backbone_pdb():
    log_path = request.args.get("log")
    name = request.args.get("name")
    status = request.args.get("status", "passed")
    output_dir = _get_output_dir(log_path)
    subdir = "failed_filters" if status == "failed_filter" else ""
    path = _safe_join(output_dir, "1_rfdiff", subdir, f"{name}.pdb")
    if not os.path.isfile(path):
        abort(404, description="pdb not found")
    return send_file(path, mimetype="chemical/x-pdb")


# ---------------------------------------------------------------------------
# Sequences (mpnn)
# ---------------------------------------------------------------------------


@app.route("/api/sequences")
def api_sequences():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    return jsonify(P.list_sequences(output_dir))


# ---------------------------------------------------------------------------
# Models (modeling stage)
# ---------------------------------------------------------------------------


@app.route("/api/models")
def api_models():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    return jsonify(P.list_models(output_dir))


def _resolve_within_output_dir(output_dir: str, path: str) -> str:
    """`path` is normally an absolute path already under output_dir (as
    returned by list_models/list_backbones); re-validate it stays inside
    output_dir regardless of whether it arrived absolute or relative."""
    rel = os.path.relpath(path, output_dir) if os.path.isabs(path) else path
    return _safe_join(output_dir, rel)


@app.route("/api/model_pdb")
def api_model_pdb():
    log_path = request.args.get("log")
    path = request.args.get("path")
    output_dir = _get_output_dir(log_path)
    full = _resolve_within_output_dir(output_dir, path)
    if not os.path.isfile(full):
        abort(404, description="pdb not found")
    return send_file(full, mimetype="chemical/x-pdb")


@app.route("/api/model_confidence")
def api_model_confidence():
    log_path = request.args.get("log")
    path = request.args.get("path")
    output_dir = _get_output_dir(log_path)
    full = _resolve_within_output_dir(output_dir, path)
    if not os.path.isfile(full):
        abort(404, description="confidence file not found")
    return jsonify(P.load_confidence(full))


# ---------------------------------------------------------------------------
# Results (final)
# ---------------------------------------------------------------------------


@app.route("/api/results")
def api_results():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    return jsonify(P.results_summary(output_dir))


@app.route("/api/final_pdb")
def api_final_pdb():
    log_path = request.args.get("log")
    name = request.args.get("name")
    output_dir = _get_output_dir(log_path)
    path = _safe_join(output_dir, "final_pdbs", name)
    if not os.path.isfile(path):
        abort(404, description="pdb not found")
    return send_file(path, mimetype="chemical/x-pdb")


@app.route("/api/final_csv")
def api_final_csv():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    csv_path = os.path.join(output_dir, "final_output.csv")
    if not os.path.isfile(csv_path):
        abort(404, description="final_output.csv not found")
    with open(csv_path, newline="", errors="replace") as f:
        reader = csv.reader(f)
        rows = list(reader)
    if not rows:
        return jsonify({"columns": [], "rows": []})
    columns, data_rows = rows[0], rows[1:]
    return jsonify({"columns": columns, "rows": data_rows})


# ---------------------------------------------------------------------------
# File browser (helper for picking the log file in the UI)
# ---------------------------------------------------------------------------


@app.route("/api/browse")
def api_browse():
    path = request.args.get("path") or os.path.expanduser("~")
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        abort(400, description="Not a directory")
    entries = []
    try:
        for name in sorted(os.listdir(path)):
            full = os.path.join(path, name)
            is_dir = os.path.isdir(full)
            if not is_dir and not (name.endswith(".out") or name.endswith(".log") or name.endswith(".txt")):
                continue
            entries.append({"name": name, "path": full, "is_dir": is_dir})
    except PermissionError:
        abort(403, description="Permission denied")
    parent = os.path.dirname(path) if path != os.path.dirname(path) else None
    return jsonify({"path": path, "parent": parent, "entries": entries})


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    debug = os.environ.get("FLASK_DEBUG", "0") == "1"
    app.run(host="0.0.0.0", port=port, debug=debug)
