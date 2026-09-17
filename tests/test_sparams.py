"""Builds the C++ twoport_fdtd core and checks the measured S11(f)/S21(f)
magnitude (a virtual VNA reading, via time-domain gating + FFT) against the
analytic ABCD-matrix result for an impedance-step slab, across several
(Zmid, slab length) combinations.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "py"))

import twoport as tp  # noqa: E402

N = 800
DX = 0.001
L1 = 250e-9
C1 = 100e-12
Z1 = np.sqrt(L1 / C1)
V = 1.0 / np.sqrt(L1 * C1)
T0 = 2.0e-9
SIGMA = 0.15e-9
N_STEPS = 40000
OBS1_FRAC, OBS2_FRAC = 0.15, 0.85


@pytest.fixture(scope="module")
def binary():
    return tp.build()


@pytest.mark.parametrize("Zmid,x1_frac,x2_frac", [
    (100.0, 0.4, 0.5),
    (20.0, 0.35, 0.55),
    (Z1, 0.4, 0.5),   # a matched "discontinuity" is really no discontinuity at all
])
def test_sparams_match_analytic(binary, tmp_path, Zmid, x1_frac, x2_frac):
    l2 = (x2_frac - x1_frac) * N * DX

    _, v1_ref, v2_ref, dt = tp.run(N, DX, N_STEPS, L1, C1, Z1, x1_frac, x2_frac, T0, SIGMA,
                                    OBS1_FRAC, OBS2_FRAC, tmp_path / "ref.csv", binary=binary)
    _, v1_act, v2_act, dt2 = tp.run(N, DX, N_STEPS, L1, C1, Zmid, x1_frac, x2_frac, T0, SIGMA,
                                     OBS1_FRAC, OBS2_FRAC, tmp_path / "act.csv", binary=binary)
    assert dt == dt2

    v_inc = v1_ref
    v_refl = v1_act - v1_ref
    v_trans = v2_act

    crop = tp.auto_crop_window(v_inc, v_refl, v_trans, threshold=0.001, margin_frac=0.4)
    freqs, S11, S21 = tp.measured_sparams(v_inc, v_refl, v_trans, dt, crop=crop)
    S11a, S21a = tp.analytic_sparams(freqs, Z1, Zmid, l2, V)

    bw = 1.0 / (2 * np.pi * SIGMA)
    mask = (freqs > 0.03 * bw) & (freqs < 0.85 * bw)

    err11 = np.max(np.abs(np.abs(S11[mask]) - np.abs(S11a[mask])))
    err21 = np.max(np.abs(np.abs(S21[mask]) - np.abs(S21a[mask])))
    assert err11 < 0.02, f"Zmid={Zmid}: max |S11| error {err11:.4f}"
    assert err21 < 0.02, f"Zmid={Zmid}: max |S21| error {err21:.4f}"


def test_windowing_choice_matters(binary, tmp_path):
    """A regression guard for the real bug this stage hit: a Hamming window
    is a raised cosine over its WHOLE length (no flat top), so it gives
    the incident/reflected/transmitted pulses different attenuation
    depending on where in the crop they land, corrupting S21 well outside
    tolerance even with a good crop. Tukey (mostly-flat, small taper) does
    not have this problem. See docs/02-sparameters.md."""
    Zmid, x1_frac, x2_frac = 100.0, 0.4, 0.5
    _, v1_ref, v2_ref, dt = tp.run(N, DX, N_STEPS, L1, C1, Z1, x1_frac, x2_frac, T0, SIGMA,
                                    OBS1_FRAC, OBS2_FRAC, tmp_path / "ref.csv", binary=binary)
    _, v1_act, v2_act, dt2 = tp.run(N, DX, N_STEPS, L1, C1, Zmid, x1_frac, x2_frac, T0, SIGMA,
                                     OBS1_FRAC, OBS2_FRAC, tmp_path / "act.csv", binary=binary)
    v_inc = v1_ref
    v_refl = v1_act - v1_ref
    v_trans = v2_act
    crop = tp.auto_crop_window(v_inc, v_refl, v_trans, threshold=0.001, margin_frac=0.4)

    import scipy.signal.windows
    start, end = crop
    n = end - start
    hamming_window = np.hamming(n)
    tukey_window = scipy.signal.windows.tukey(n, alpha=0.1)

    freqs = np.fft.rfftfreq(n, d=dt)
    l2 = (x2_frac - x1_frac) * N * DX
    _, S21a = tp.analytic_sparams(freqs, Z1, Zmid, l2, V)
    bw = 1.0 / (2 * np.pi * SIGMA)
    mask = (freqs > 0.03 * bw) & (freqs < 0.85 * bw)

    def s21_with_window(window):
        A_inc = np.fft.rfft(v_inc[start:end] * window)
        A_trans = np.fft.rfft(v_trans[start:end] * window)
        return A_trans / A_inc

    err_hamming = np.max(np.abs(np.abs(s21_with_window(hamming_window)[mask]) - np.abs(S21a[mask])))
    err_tukey = np.max(np.abs(np.abs(s21_with_window(tukey_window)[mask]) - np.abs(S21a[mask])))
    assert err_tukey < 0.02
    assert err_hamming > 0.1, "expected the Hamming window to visibly fail this comparison"
