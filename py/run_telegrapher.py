"""One-command reflection driver: builds the C++ core, runs the reflection
cases (matched/open/short/arbitrary), measures Gamma against the analytic
formula, and produces the hero figure: snapshots of the pulse traveling
down the line and reflecting off an impedance step.

Run from the repo root:
    python py/run_telegrapher.py
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
import fdtd  # noqa: E402

IMG_DIR = ROOT / "docs" / "imgs"

N = 400
DX = 0.001
L = 250e-9
C = 100e-12
Z0 = np.sqrt(L / C)
T0 = 3e-9
SIGMA = 0.5e-9
N_STEPS = 20000
OBS_FRAC = 0.5


def fig_reflection_table(binary):
    cases = [("matched (ZL=Z0)", Z0), ("open", -1), ("short", 0), ("ZL=25 ohm", 25), ("ZL=100 ohm", 100)]
    rows = []
    for label, ZL in cases:
        out_csv = ROOT / "_tmp_obs.csv"
        times, v_obs = fdtd.run(N, DX, N_STEPS, L, C, ZL, T0, SIGMA, OBS_FRAC, out_csv, binary=binary)
        gamma, incident, t1, reflected, t2 = fdtd.measure_reflection_coefficient(times, v_obs, SIGMA)
        analytic = fdtd.analytic_gamma(ZL, Z0)
        rows.append((label, ZL, gamma, analytic, abs(gamma - analytic)))
        out_csv.unlink(missing_ok=True)

    print(f"{'case':16s} {'ZL':>8s} {'measured':>10s} {'analytic':>10s} {'abs err':>8s}")
    for label, ZL, gamma, analytic, err in rows:
        print(f"{label:16s} {ZL:8.2f} {gamma:10.4f} {analytic:10.4f} {err:8.4f}")
    return rows


def fig_wave_propagation(binary):
    """The hero image: V(x) at several times, showing a pulse launched from
    a matched source, traveling down the line, and reflecting off an open
    end (Gamma=+1, so the reflected pulse comes back right-side up and
    roughly doubles the voltage where it overlaps the tail of the
    incident pulse)."""
    obs_csv = ROOT / "_tmp_obs.csv"
    snap_csv = ROOT / "_tmp_snap.csv"
    fdtd.run(N, DX, N_STEPS, L, C, -1, T0, SIGMA, OBS_FRAC, obs_csv,
             snap_csv=snap_csv, n_snapshots=8, binary=binary)
    times, X = fdtd.load_snapshots(snap_csv)
    x = np.arange(X.shape[1]) * DX * 100  # cm

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, t in enumerate(times):
        ax.plot(x, X[i] + i * 0.3, color="C0", linewidth=1.2)
        ax.text(x[-1] * 1.01, X[i, -1] + i * 0.3, f"{t*1e9:.1f} ns", fontsize=7, va="center")
    ax.axvline(x[-1], color="gray", linestyle="--", linewidth=1, label="open end (load)")
    ax.set_xlabel("position along line (cm)")
    ax.set_ylabel("V(x), offset per snapshot for clarity")
    ax.set_title("Pulse launched from a matched source, reflecting off an open end")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(IMG_DIR / "wave_reflection.png", dpi=150)
    plt.close(fig)

    obs_csv.unlink(missing_ok=True)
    snap_csv.unlink(missing_ok=True)


def fig_observation_trace(binary):
    """A simpler, more classic-looking figure: the observation-point
    voltage over time for the open case, with the incident and reflected
    pulses both visible and labeled."""
    out_csv = ROOT / "_tmp_obs.csv"
    times, v_obs = fdtd.run(N, DX, N_STEPS, L, C, -1, T0, SIGMA, OBS_FRAC, out_csv, binary=binary)
    gamma, incident, t1, reflected, t2 = fdtd.measure_reflection_coefficient(times, v_obs, SIGMA)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(times * 1e9, v_obs, color="C0")
    ax.annotate("incident", xy=(t1 * 1e9, incident), xytext=(t1 * 1e9 - 1, incident + 0.1),
                arrowprops=dict(arrowstyle="->"))
    ax.annotate("reflected", xy=(t2 * 1e9, reflected), xytext=(t2 * 1e9 + 0.3, reflected + 0.1),
                arrowprops=dict(arrowstyle="->"))
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("V at observation point")
    ax.set_title(f"Open-circuit reflection: measured Gamma = {gamma:.3f} (analytic +1.000)")
    fig.tight_layout()
    fig.savefig(IMG_DIR / "observation_trace.png", dpi=150)
    plt.close(fig)
    out_csv.unlink(missing_ok=True)


def main():
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    binary = fdtd.build()
    print(f"[build] {binary}")

    rows = fig_reflection_table(binary)
    max_err = max(r[4] for r in rows)
    print(f"\n[check] max |measured - analytic| Gamma error: {max_err:.4f}")

    fig_wave_propagation(binary)
    fig_observation_trace(binary)
    print(f"\nFigures written to {IMG_DIR}")


if __name__ == "__main__":
    main()
