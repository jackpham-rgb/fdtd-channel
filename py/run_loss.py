"""Loss & dispersion driver: add series resistance R and shunt conductance G to the
line, and show insertion loss actually rolling off with frequency instead
of staying flat forever, checked against the analytic lossy propagation
constant.

Run from the repo root:
    python py/run_loss.py
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
import lossy  # noqa: E402

IMG_DIR = ROOT / "docs" / "imgs"

N = 800
DX = 0.001
L = 250e-9
C = 100e-12
Z0 = np.sqrt(L / C)
R = 200.0        # ohm/m: dominant loss mechanism here (conductor-loss-like)
G = 0.0002       # S/m: kept small so its own corner sits far below our band
T0 = 2.0e-9
SIGMA = 0.15e-9
N_STEPS = 300000
OBS_NEAR_FRAC, OBS_FAR_FRAC = 0.1, 0.9
DIST = (OBS_FAR_FRAC - OBS_NEAR_FRAC) * N * DX

FIT_BAND = (1.1e8, 1.0e9)  # see docs/03-loss.md for why the band starts here


def fig_loss_vs_freq(binary):
    times, v_near, v_far, dt = lossy.run(
        N, DX, N_STEPS, L, C, R, G, T0, SIGMA, OBS_NEAR_FRAC, OBS_FAR_FRAC,
        ROOT / "_tmp_lossy.csv", binary=binary)
    crop = lossy.auto_crop_window(v_near, v_far, threshold=1e-4, margin_frac=0.3)
    freqs, H = lossy.measured_transfer_function(v_near, v_far, dt, crop=crop)
    gamma = lossy.analytic_gamma(freqs, L, C, R, G)
    H_analytic = np.exp(-gamma * DIST)

    mask = (freqs > 3e7) & (freqs < 1.2e9)
    f_ghz = freqs[mask] / 1e9
    np_per_m_to_db_per_m = 20 * np.log10(np.e)
    alpha_meas = -np.log(np.abs(H[mask])) / DIST * np_per_m_to_db_per_m
    alpha_ana = np.real(gamma[mask]) * np_per_m_to_db_per_m

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(f_ghz, alpha_meas, color="C0", label="measured (FDTD, two-point transfer function)")
    ax.semilogx(f_ghz, alpha_ana, "--", color="C1", label="analytic: Re(sqrt((R+jwL)(G+jwC)))")
    ax.axvline(R / (2 * np.pi * L) / 1e9, color="gray", linestyle=":", linewidth=1,
               label=f"R/L corner ({R/(2*np.pi*L)/1e6:.0f} MHz)")
    ax.set_xlabel("frequency (GHz)")
    ax.set_ylabel("attenuation (dB/m)")
    ax.set_title(f"Loss rolls off with frequency: R={R:.0f} ohm/m, G={G} S/m line")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(IMG_DIR / "loss_vs_freq.png", dpi=150)
    plt.close(fig)

    band_mask = (freqs[mask] > FIT_BAND[0]) & (freqs[mask] < FIT_BAND[1])
    max_err_db = np.max(np.abs(alpha_meas[band_mask] - alpha_ana[band_mask]))
    print(f"[check] max attenuation error over {FIT_BAND[0]/1e6:.0f}-{FIT_BAND[1]/1e6:.0f} MHz: {max_err_db:.3f} dB/m")

    return times, v_near, v_far


def fig_pulse_dispersion(times, v_near, v_far):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(times * 1e9, v_near, label=f"near (x={OBS_NEAR_FRAC*N*DX:.2f} m)", color="C0")
    ax.plot(times * 1e9, v_far, label=f"far (x={OBS_FAR_FRAC*N*DX:.2f} m)", color="C1")
    ax.set_xlim(0, 15)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("V")
    ax.set_title(f"Loss + dispersion over {DIST:.2f} m: the far pulse is both smaller AND wider")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(IMG_DIR / "loss_dispersion_time_domain.png", dpi=150)
    plt.close(fig)


def main():
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    binary = lossy.build()
    print(f"[build] {binary}")

    times, v_near, v_far = fig_loss_vs_freq(binary)
    fig_pulse_dispersion(times, v_near, v_far)
    print(f"\nFigures written to {IMG_DIR}")

    (ROOT / "_tmp_lossy.csv").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
