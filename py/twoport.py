"""Build and run the twoport_fdtd C++ core, and turn its output into a
virtual network analyzer: S11(f) and S21(f) for an impedance-step
discontinuity, extracted the same way a real VNA does it (time-domain
gating), cross-checked against the analytic ABCD-matrix result for a single
transmission-line section between two matched lines.
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
SRC_FILE = ROOT / "src" / "twoport_fdtd.cpp"
BIN_FILE = ROOT / "twoport_fdtd.exe" if sys.platform == "win32" else ROOT / "twoport_fdtd"


def build(force: bool = False) -> pathlib.Path:
    if not force and BIN_FILE.exists() and BIN_FILE.stat().st_mtime > SRC_FILE.stat().st_mtime:
        return BIN_FILE
    cmd = ["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-o", str(BIN_FILE), str(SRC_FILE)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"build failed:\n{result.stdout}\n{result.stderr}")
    return BIN_FILE


def run(N, dx, n_steps, L1, C1, Zmid, x1_frac, x2_frac, t0, sigma,
        obs1_frac, obs2_frac, out_csv, binary=None):
    """Run the compiled core. Returns (times, v_obs1, v_obs2, dt)."""
    binary = binary or build()
    args = [str(binary), str(N), str(dx), str(n_steps), str(L1), str(C1), str(Zmid),
            str(x1_frac), str(x2_frac), str(t0), str(sigma),
            str(obs1_frac), str(obs2_frac), str(out_csv)]
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"twoport_fdtd failed:\n{result.stdout}\n{result.stderr}")
    m = re.search(r"dt=([0-9.eE+-]+)", result.stdout)
    dt = float(m.group(1)) if m else None
    times, v1, v2 = load_two_obs(out_csv)
    return times, v1, v2, dt


def load_two_obs(csv_path):
    times, v1, v2 = [], [], []
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            times.append(float(row[1]))
            v1.append(float(row[2]))
            v2.append(float(row[3]))
    return np.array(times), np.array(v1), np.array(v2)


def auto_crop_window(*signals, threshold=0.01, margin_frac=0.5):
    """Find a (start, end) sample range that safely contains all the given
    signals' significant energy, with margin on both sides so a later
    window function's tapered edges fall on near-zero content instead of
    clipping real signal. `threshold` is a fraction of the global peak
    magnitude (across all signals) used to find first/last "significant"
    samples; `margin_frac` extends that span by this fraction of its own
    length on each side.
    """
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


def measured_sparams(v_incident, v_reflected, v_transmitted, dt, crop=None, fmax=None):
    """Hamming-window and FFT the incident/reflected/transmitted time
    series (all the same length, sampled at the same dt), and return
    (freqs, S11, S21) as complex arrays, using S11 = FFT(reflected) /
    FFT(incident) and S21 = FFT(transmitted) / FFT(incident). Valid without
    impedance renormalization because both ports share the same reference
    impedance Z1 here.

    `crop`, if given, is a (start, end) sample-index pair applied to all
    three series before windowing. This matters: a Hamming window tapers
    toward zero at the EDGES of whatever array it's given, not around each
    signal's own center. If the simulation runs much longer than the pulses
    actually need (a lot of empty padding), the incident/reflected/
    transmitted pulses end up sitting at different distances from those
    tapered edges and get suppressed by different amounts, corrupting the
    S11/S21 ratio.

    Cropping alone isn't enough, though: a Hamming window is a raised
    cosine over its ENTIRE length (0.08 at the very edges, but still only
    0.54 a quarter of the way in, reaching 1.0 only exactly at the center).
    It has no flat top. So even with a generous, well-placed crop, the
    incident pulse, the first reflection, and later (delayed) reflections
    each sit at a different distance from the window's center and each get
    a DIFFERENT attenuation factor, which corrupts the ratio just the same
    (found empirically: S11 measured 2-4x too high, consistently, even
    after fixing the crop). A Tukey window with a small taper fraction is
    flat (=1.0) across the middle ~(1-alpha) of its length and only tapers
    over the outer edges, so signals anywhere in that flat region get equal
    treatment; alpha=0.1 here means a 90%-flat, 10%-tapered window.
    """
    if crop is not None:
        start, end = crop
        v_incident = v_incident[start:end]
        v_reflected = v_reflected[start:end]
        v_transmitted = v_transmitted[start:end]
    n = len(v_incident)
    window = scipy.signal.windows.tukey(n, alpha=0.1)
    freqs = np.fft.rfftfreq(n, d=dt)

    A_inc = np.fft.rfft(v_incident * window)
    A_refl = np.fft.rfft(v_reflected * window)
    A_trans = np.fft.rfft(v_transmitted * window)

    # Avoid dividing by near-zero incident spectrum outside the source's
    # bandwidth, where S11/S21 are meaningless anyway.
    eps = 1e-12 * np.max(np.abs(A_inc))
    denom = np.where(np.abs(A_inc) > eps, A_inc, np.nan)
    with np.errstate(invalid="ignore"):
        S11 = A_refl / denom
        S21 = A_trans / denom

    if fmax is not None:
        mask = freqs <= fmax
        return freqs[mask], S11[mask], S21[mask]
    return freqs, S11, S21


def analytic_sparams(freqs, Z1, Zmid, l2, v):
    """S11(f), S21(f) for a single uniform line section of impedance Zmid
    and length l2, embedded between two semi-infinite (matched) lines of
    impedance Z1, via the line section's ABCD matrix (Pozar, Microwave
    Engineering, table of two-port conversions). This accounts for ALL
    orders of internal reflection within the section, not just the first
    bounce, since it's the exact frequency-domain (steady-state) result.
    """
    freqs = np.asarray(freqs, dtype=float)
    beta = 2.0 * np.pi * freqs / v
    theta = beta * l2

    A = np.cos(theta)
    B = 1j * Zmid * np.sin(theta)
    Cm = 1j * np.sin(theta) / Zmid
    D = np.cos(theta)

    denom = A + B / Z1 + Cm * Z1 + D
    S11 = (A + B / Z1 - Cm * Z1 - D) / denom
    S21 = 2.0 / denom
    return S11, S21
