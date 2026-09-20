"""S-parameters driver: turn the FDTD line into a virtual network analyzer.

Runs a reference (uniform) simulation and a discontinuity simulation, uses
their difference to separate incident/reflected/transmitted waves at the
two port reference planes, FFTs them into S11(f)/S21(f), and compares
against the analytic ABCD-matrix result for the same impedance-step slab.

Run from the repo root:
    python py/run_sparams.py
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "py"))
import twoport as tp  # noqa: E402

IMG_DIR = ROOT / "docs" / "imgs"

N = 800
DX = 0.001
L1 = 250e-9
C1 = 100e-12
Z1 = np.sqrt(L1 / C1)
V = 1.0 / np.sqrt(L1 * C1)
T0 = 2.0e-9
SIGMA = 0.15e-9
N_STEPS = 40000
X1_FRAC, X2_FRAC = 0.4, 0.5
OBS1_FRAC, OBS2_FRAC = 0.15, 0.85
ZMID = 100.0
L2 = (X2_FRAC - X1_FRAC) * N * DX


def run_both(binary, Zmid):
    """Reference run (Zmid=Z1, no discontinuity) and actual run (the real
    Zmid), then separate incident/reflected/transmitted. See
    docs/02-sparameters.md for why the reference-run subtraction is needed
    at the port-1 plane but not at port 2."""
    times, v1_ref, v2_ref, dt = tp.run(
        N, DX, N_STEPS, L1, C1, Z1, X1_FRAC, X2_FRAC, T0, SIGMA,
        OBS1_FRAC, OBS2_FRAC, ROOT / "_tmp_ref.csv", binary=binary)
    times2, v1_act, v2_act, dt2 = tp.run(
        N, DX, N_STEPS, L1, C1, Zmid, X1_FRAC, X2_FRAC, T0, SIGMA,
        OBS1_FRAC, OBS2_FRAC, ROOT / "_tmp_act.csv", binary=binary)
    assert dt == dt2, "reference and actual runs must share the same dt"

    v_inc = v1_ref
    v_refl = v1_act - v1_ref
    v_trans = v2_act
    return times, v_inc, v_refl, v_trans, dt


def fig_sparams(binary):
    times, v_inc, v_refl, v_trans, dt = run_both(binary, ZMID)
    crop = tp.auto_crop_window(v_inc, v_refl, v_trans, threshold=0.001, margin_frac=0.4)
    freqs, S11, S21 = tp.measured_sparams(v_inc, v_refl, v_trans, dt, crop=crop)
    S11a, S21a = tp.analytic_sparams(freqs, Z1, ZMID, L2, V)

    bw = 1.0 / (2 * np.pi * SIGMA)
    mask = (freqs > 0.02 * bw) & (freqs < 0.9 * bw)
    f_ghz = freqs[mask] / 1e9

    def to_db(x):
        return 20 * np.log10(np.abs(x))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(f_ghz, to_db(S11[mask]), color="C0", label="|S11| measured (return loss)")
    ax.plot(f_ghz, to_db(S11a[mask]), "--", color="C0", alpha=0.6, label="|S11| analytic")
    ax.plot(f_ghz, to_db(S21[mask]), color="C1", label="|S21| measured (insertion loss)")
    ax.plot(f_ghz, to_db(S21a[mask]), "--", color="C1", alpha=0.6, label="|S21| analytic")
    ax.set_xlabel("frequency (GHz)")
    ax.set_ylabel("magnitude (dB)")
    ax.set_title(f"S11/S21 of a {ZMID:.0f} ohm slab in a {Z1:.0f} ohm line, FDTD vs. analytic")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(IMG_DIR / "sparams.png", dpi=150)
    plt.close(fig)

    max_err_11 = np.max(np.abs(np.abs(S11[mask]) - np.abs(S11a[mask])))
    max_err_21 = np.max(np.abs(np.abs(S21[mask]) - np.abs(S21a[mask])))
    print(f"[check] max |S11| magnitude error: {max_err_11:.4f}")
    print(f"[check] max |S21| magnitude error: {max_err_21:.4f}")
    return times, v_inc, v_refl, v_trans, dt


def fig_time_domain(times, v_inc, v_refl, v_trans):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(times * 1e9, v_inc, label="incident (port 1)", color="C0")
    ax.plot(times * 1e9, v_refl, label="reflected (port 1)", color="C2")
    ax.plot(times * 1e9, v_trans, label="transmitted (port 2)", color="C1")
    ax.set_xlim(0, 12)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("V")
    ax.set_title(f"Time-domain gating: a {ZMID:.0f} ohm slab embedded in a {Z1:.0f} ohm line")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(IMG_DIR / "sparams_time_domain.png", dpi=150)
    plt.close(fig)


def main():
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    binary = tp.build()
    print(f"[build] {binary}")

    times, v_inc, v_refl, v_trans, dt = fig_sparams(binary)
    fig_time_domain(times, v_inc, v_refl, v_trans)
    print(f"\nFigures written to {IMG_DIR}")

    for f in (ROOT / "_tmp_ref.csv", ROOT / "_tmp_act.csv"):
        f.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
