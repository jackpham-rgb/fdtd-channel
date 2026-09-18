"""Build and run the lossy_fdtd C++ core, and extract insertion loss
(|S21|, as a transfer function between two observation points) from the
simulated fields, compared against the analytic complex propagation
constant of a lossy transmission line.
"""
from __future__ import annotations

import csv
import pathlib
import re
import subprocess
import sys

import numpy as np
import scipy.signal.windows

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC_FILE = ROOT / "src" / "lossy_fdtd.cpp"
BIN_FILE = ROOT / "lossy_fdtd.exe" if sys.platform == "win32" else ROOT / "lossy_fdtd"


def build(force: bool = False) -> pathlib.Path:
    if not force and BIN_FILE.exists() and BIN_FILE.stat().st_mtime > SRC_FILE.stat().st_mtime:
        return BIN_FILE
    cmd = ["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-o", str(BIN_FILE), str(SRC_FILE)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"build failed:\n{result.stdout}\n{result.stderr}")
    return BIN_FILE


def run(N, dx, n_steps, L, C, R, G, t0, sigma, obs_near_frac, obs_far_frac, out_csv, binary=None):
    """Run the compiled core. Returns (times, v_near, v_far, dt)."""
    binary = binary or build()
    args = [str(binary), str(N), str(dx), str(n_steps), str(L), str(C), str(R), str(G),
            str(t0), str(sigma), str(obs_near_frac), str(obs_far_frac), str(out_csv)]
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"lossy_fdtd failed:\n{result.stdout}\n{result.stderr}")
    m = re.search(r"dt=([0-9.eE+-]+)", result.stdout)
    dt = float(m.group(1)) if m else None
    times, v_near, v_far = load_two_obs(out_csv)
    return times, v_near, v_far, dt


def load_two_obs(csv_path):
    times, v_near, v_far = [], [], []
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            times.append(float(row[1]))
            v_near.append(float(row[2]))
            v_far.append(float(row[3]))
    return np.array(times), np.array(v_near), np.array(v_far)


def auto_crop_window(*signals, threshold=0.01, margin_frac=0.5):
    """Same idea as twoport.py's version: find a sample range containing
    all the given signals' significant energy, with margin so a window
    function's tapered edges don't land on real content."""
    peak = max(np.max(np.abs(s)) for s in signals)
    cutoff = threshold * peak
    first = len(signals[0])
    last = 0
    for s in signals:
        above = np.flatnonzero(np.abs(s) > cutoff)
        if above.size:
            first = min(first, above[0])
            last = max(last, above[-1])
    span = last - first
    margin = int(margin_frac * span)
    n = len(signals[0])
    start = max(0, first - margin)
    end = min(n, last + margin)
    return start, end


def measured_transfer_function(v_near, v_far, dt, crop=None, nfft=None):
    """H(f) = FFT(v_far) / FFT(v_near): the transfer function between the
    two observation points, i.e. exactly exp(-gamma(f) * distance) for a
    uniform lossy line, IF the crop captures the direct wave at both
    points and nothing else (no reflection from the matched load, which
    is why this core always uses a matched load). Windowed with a Tukey
    window (see docs/02-sparameters.md and twoport.py for why not Hamming:
    it has no flat top and biases signals differently depending on where
    in the record they land).

    `nfft`, if given and larger than the cropped length, zero-pads before
    the FFT. This doesn't add information the short capture didn't already
    have -- it can't separate two frequencies closer than 1/(capture
    duration) -- but the underlying transfer function IS smooth in
    frequency, so zero-padding correctly interpolates that smooth curve
    onto a finer frequency grid instead of only showing it at the native
    (coarse) bin spacing. Useful for seeing a loss transition (a "corner"
    frequency) that falls between two native bins.
    """
    if crop is not None:
        start, end = crop
        v_near = v_near[start:end]
        v_far = v_far[start:end]
    n = len(v_near)
    window = scipy.signal.windows.tukey(n, alpha=0.1)
    nfft = max(n, nfft) if nfft is not None else n
    freqs = np.fft.rfftfreq(nfft, d=dt)

    A_near = np.fft.rfft(v_near * window, n=nfft)
    A_far = np.fft.rfft(v_far * window, n=nfft)

    eps = 1e-12 * np.max(np.abs(A_near))
    denom = np.where(np.abs(A_near) > eps, A_near, np.nan)
    with np.errstate(invalid="ignore"):
        H = A_far / denom
    return freqs, H


def analytic_gamma(freqs, L, C, R, G):
    """Complex propagation constant gamma(f) = alpha(f) + j*beta(f) for a
    lossy line: gamma = sqrt((R + j*2*pi*f*L) * (G + j*2*pi*f*C))."""
    freqs = np.asarray(freqs, dtype=float)
    w = 2.0 * np.pi * freqs
    return np.sqrt((R + 1j * w * L) * (G + 1j * w * C))
