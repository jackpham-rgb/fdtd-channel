"""Build and run the telegrapher_fdtd C++ core, and extract the reflection
coefficient from its output.

Keeping the physics in the compiled core and the science here (loading CSVs,
finding pulses, plotting) matches this project's stated split. Nothing in
this file does numerical time-stepping; it only orchestrates and analyzes.
"""
from __future__ import annotations

import csv
import pathlib
import subprocess
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC_FILE = ROOT / "src" / "telegrapher_fdtd.cpp"
BIN_FILE = ROOT / "telegrapher_fdtd.exe" if sys.platform == "win32" else ROOT / "telegrapher_fdtd"


def build(force: bool = False) -> pathlib.Path:
    """Compile the C++ core if it doesn't exist yet or the source changed."""
    if not force and BIN_FILE.exists() and BIN_FILE.stat().st_mtime > SRC_FILE.stat().st_mtime:
        return BIN_FILE
    cmd = ["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-o", str(BIN_FILE), str(SRC_FILE)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"build failed:\n{result.stdout}\n{result.stderr}")
    return BIN_FILE


def run(N, dx, n_steps, L, C, ZL, t0, sigma, obs_frac, out_csv,
        snap_csv=None, n_snapshots=0, binary=None):
    """Run the compiled FDTD core and return (times, v_obs) arrays."""
    binary = binary or build()
    args = [str(binary), str(N), str(dx), str(n_steps), str(L), str(C), str(ZL),
            str(t0), str(sigma), str(obs_frac), str(out_csv)]
    if snap_csv is not None:
        args += [str(snap_csv), str(n_snapshots)]
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"telegrapher_fdtd failed:\n{result.stdout}\n{result.stderr}")
    return load_observation(out_csv)


def load_observation(csv_path):
    times, vals = [], []
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            times.append(float(row[1]))
            vals.append(float(row[2]))
    return np.array(times), np.array(vals)


def load_snapshots(csv_path):
    """Returns (times, X) where X is an (n_snapshots, N+1) array of V(x)."""
    times = []
    rows = []
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            times.append(float(row[0]))
            rows.append([float(v) for v in row[1:]])
    return np.array(times), np.array(rows)


def measure_reflection_coefficient(times, v_obs, sigma, mask_sigmas=3.0):
    """Find the incident pulse's peak, then the reflected pulse's peak
    (searching only outside a +-mask_sigmas*sigma window around the
    incident peak's time, so the search can't just re-find the incident
    pulse's own tail). Returns (gamma, incident_peak, incident_time,
    reflected_peak, reflected_time).

    This only works because the source is a MATCHED (Thevenin) source: it
    absorbs whatever returns to it instead of re-reflecting it, so the
    first two pulses seen at any observation point really are "the pulse
    going out" and "the pulse coming back once," not an uncontrolled mix of
    multiple bounces.
    """
    idx1 = int(np.argmax(np.abs(v_obs)))
    incident_peak = v_obs[idx1]
    t1 = times[idx1]

    mask = np.abs(times - t1) < (mask_sigmas * sigma)
    remaining = v_obs.copy()
    remaining[mask] = 0.0
    idx2 = int(np.argmax(np.abs(remaining)))
    reflected_peak = remaining[idx2]
    t2 = times[idx2]

    gamma = reflected_peak / incident_peak
    return gamma, incident_peak, t1, reflected_peak, t2


def analytic_gamma(ZL, Z0):
    if ZL < 0:
        return 1.0  # open circuit
    return (ZL - Z0) / (ZL + Z0)
