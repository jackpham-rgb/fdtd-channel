"""Builds the C++ lossy_fdtd core and checks the measured attenuation
(alpha(f), from a two-point transfer function) against the analytic lossy
propagation constant Re(sqrt((R+jwL)(G+jwC))), for a few (R, G) values.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "py"))

import lossy  # noqa: E402

N = 800
DX = 0.001
L = 250e-9
C = 100e-12
T0 = 2.0e-9
SIGMA = 0.15e-9
N_STEPS = 300000
OBS_NEAR_FRAC, OBS_FAR_FRAC = 0.1, 0.9
DIST = (OBS_FAR_FRAC - OBS_NEAR_FRAC) * N * DX

# Below this, a ~20-30ns capture can't reliably resolve the transfer
# function: near-DC FFT bins are sensitive to how the Gaussian pulse's tail
# gets truncated by the crop, which shows up as a measured attenuation that
# flattens out too early instead of continuing down toward the true
# sqrt(R*G) plateau. See docs/03-loss.md.
BAND = (1.1e8, 1.0e9)


@pytest.fixture(scope="module")
def binary():
    return lossy.build()


@pytest.mark.parametrize("R,G", [
    (80.0, 0.0001),
    (200.0, 0.0002),
    (400.0, 0.0005),
])
def test_attenuation_matches_analytic(binary, tmp_path, R, G):
    _, v_near, v_far, dt = lossy.run(N, DX, N_STEPS, L, C, R, G, T0, SIGMA,
                                      OBS_NEAR_FRAC, OBS_FAR_FRAC,
                                      tmp_path / "lossy.csv", binary=binary)
    crop = lossy.auto_crop_window(v_near, v_far, threshold=1e-4, margin_frac=0.3)
    freqs, H = lossy.measured_transfer_function(v_near, v_far, dt, crop=crop)
    gamma = lossy.analytic_gamma(freqs, L, C, R, G)

    mask = (freqs > BAND[0]) & (freqs < BAND[1])
    alpha_meas = -np.log(np.abs(H[mask])) / DIST
    alpha_ana = np.real(gamma[mask])

    rel_err = np.max(np.abs(alpha_meas - alpha_ana) / alpha_ana)
    assert rel_err < 0.15, f"R={R}, G={G}: max relative attenuation error {rel_err:.3f}"


def test_lossless_limit_matches_stage_a(binary, tmp_path):
    """R=G=0 should reduce to a plain lossless matched line: no measurable
    attenuation between the two observation points."""
    _, v_near, v_far, dt = lossy.run(N, DX, N_STEPS, L, C, 0.0, 0.0, T0, SIGMA,
                                      OBS_NEAR_FRAC, OBS_FAR_FRAC,
                                      tmp_path / "lossless.csv", binary=binary)
    peak_near = np.max(np.abs(v_near))
    peak_far = np.max(np.abs(v_far))
    assert abs(peak_far - peak_near) / peak_near < 0.02
