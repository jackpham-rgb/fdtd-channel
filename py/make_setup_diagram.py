"""Draw the three simulation setups to scale, with the values used in
py/run_telegrapher.py, py/run_sparams.py and py/run_loss.py.

    python py/make_setup_diagram.py
"""
from __future__ import annotations

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

ROOT = pathlib.Path(__file__).resolve().parents[1]

INK = "#1b1f24"
LINE = "#1f5aa6"
SLAB = "#c8421f"
PROBE = "#1a7f72"
GRAY = "#7a828c"
LOSS = "#b07d12"


def source(ax, x, y, label):
    ax.add_patch(Ellipse((x, y - 0.32), 0.018, 0.1, fill=False, color=INK, lw=1.2))
    ax.text(x, y - 0.32, "~", ha="center", va="center", fontsize=9, color=INK)
    ax.plot([x, x], [y - 0.27, y], color=INK, lw=1.2)
    ax.plot([x, x], [y - 0.37, y - 0.5], color=INK, lw=1.2)
    ax.text(x - 0.02, y - 0.32, label, ha="right", va="center", fontsize=8, color=INK)


def resistor(ax, x, y, label, side="right"):
    ax.add_patch(Rectangle((x - 0.012, y - 0.42), 0.024, 0.24, facecolor="white", edgecolor=INK, lw=1.2))
    ax.plot([x, x], [y - 0.18, y], color=INK, lw=1.2)
    ax.plot([x, x], [y - 0.42, y - 0.5], color=INK, lw=1.2)
    ha = "left" if side == "right" else "right"
    dx = 0.02 if side == "right" else -0.02
    ax.text(x + dx, y - 0.3, label, ha=ha, va="center", fontsize=8, color=INK)


def ground(ax, x0, x1, y):
    ax.plot([x0, x1], [y - 0.5, y - 0.5], color=GRAY, lw=1.0)
    ax.text((x0 + x1) / 2, y - 0.6, "return path (ground)", ha="center", va="top", fontsize=7.5, color=GRAY)


def probe(ax, x, y, name):
    ax.plot([x, x], [y, y + 0.28], color=PROBE, lw=1.4)
    ax.plot(x, y, "o", color=PROBE, ms=4)
    ax.text(x, y + 0.32, name, ha="center", va="bottom", fontsize=8, color=PROBE)


def ruler(ax, x0, x1, y, text):
    ax.annotate("", (x0, y), (x1, y), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.9))
    ax.text((x0 + x1) / 2, y - 0.06, text, ha="center", va="top", fontsize=7.5, color=GRAY)


def panel(ax, title, note):
    ax.set_xlim(-0.12, 0.92)
    ax.set_ylim(-0.95, 0.75)
    ax.axis("off")
    ax.text(-0.12, 0.72, title, fontsize=11, fontweight="bold", color=INK, va="top")
    ax.text(-0.12, 0.56, note, fontsize=8, color=INK, va="top")


def main():
    plt.rcParams["font.family"] = ["Bahnschrift", "DejaVu Sans"]
    fig, axes = plt.subplots(3, 1, figsize=(11, 9.6))
    y = 0.0

    # 1. Reflection: 0.4 m, Z0 = 50 ohm, probe at the middle, load ZL
    ax = axes[0]
    panel(ax, "1. Reflection", "A Gaussian pulse travels down a 50 ohm line, hits the load, and part of it comes back. "
          "The probe sees both and gives the reflection coefficient.")
    ax.plot([0, 0.4], [y, y], color=LINE, lw=4, solid_capstyle="butt")
    ax.text(0.2, -0.08, "Z0 = 50 ohm  (L = 250 nH/m, C = 100 pF/m)", ha="center", va="top", fontsize=8, color=LINE)
    source(ax, 0, y, "matched source\nRs = 50 ohm")
    resistor(ax, 0.4, y, "load ZL: open, short,\n25, 50 or 100 ohm")
    probe(ax, 0.2, y, "probe (middle)")
    ground(ax, 0, 0.4, y)
    ruler(ax, 0, 0.4, y - 0.72, "0.4 m (400 cells of 1 mm)")

    # 2. S-parameters: 0.8 m, 100 ohm slab from 0.32 to 0.40 m, probes at 0.12 and 0.68 m
    ax = axes[1]
    panel(ax, "2. S-parameters (a virtual network analyzer)",
          "A 100 ohm section sits inside a 50 ohm line. Port 1 sees what reflects (S11), port 2 what gets through (S21). "
          "Both ends are matched, so nothing else reflects.")
    ax.plot([0, 0.32], [y, y], color=LINE, lw=4, solid_capstyle="butt")
    ax.plot([0.40, 0.8], [y, y], color=LINE, lw=4, solid_capstyle="butt")
    ax.plot([0.32, 0.40], [y, y], color=SLAB, lw=8, solid_capstyle="butt")
    ax.text(0.16, 0.05, "50 ohm", ha="center", va="bottom", fontsize=8, color=LINE)
    ax.text(0.6, 0.05, "50 ohm", ha="center", va="bottom", fontsize=8, color=LINE)
    ax.text(0.36, -0.1, "100 ohm\n80 mm", ha="center", va="top", fontsize=8, color=SLAB)
    source(ax, 0, y, "matched source\n50 ohm")
    resistor(ax, 0.8, y, "matched\n50 ohm")
    probe(ax, 0.12, y, "port 1")
    probe(ax, 0.68, y, "port 2")
    ground(ax, 0, 0.8, y)
    ruler(ax, 0, 0.8, y - 0.72, "0.8 m (800 cells); each case is also run once with no section, to subtract the incident wave")

    # 3. Loss: 0.8 m lossy line, probes 0.64 m apart
    ax = axes[2]
    panel(ax, "3. Loss and dispersion",
          "Resistance along the wire and leakage through the board are added everywhere on the line. "
          "Comparing the two probes gives the loss per meter at each frequency.")
    ax.plot([0, 0.8], [y, y], color=LOSS, lw=4, solid_capstyle="butt")
    ax.text(0.4, 0.05, "50 ohm line with R = 200 ohm/m in series and G = 0.0002 S/m to ground", ha="center",
            va="bottom", fontsize=8, color=LOSS)
    source(ax, 0, y, "matched source\n50 ohm")
    resistor(ax, 0.8, y, "matched\n50 ohm")
    probe(ax, 0.08, y, "near")
    probe(ax, 0.72, y, "far")
    ground(ax, 0, 0.8, y)
    ruler(ax, 0.08, 0.72, y - 0.72, "0.64 m between the probes")

    fig.tight_layout(h_pad=1.5)
    out = ROOT / "docs" / "imgs" / "setups.png"
    fig.savefig(out, dpi=160, facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
