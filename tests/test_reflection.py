"""Builds the C++ FDTD core (needs g++ on PATH) and checks the measured
reflection coefficient against the analytic Gamma = (ZL - Z0) / (ZL + Z0)
for several terminations.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "py"))

import fdtd  # noqa: E402

N = 400
DX = 0.001
L = 250e-9
C = 100e-12
Z0 = np.sqrt(L / C)
T0 = 3e-9
SIGMA = 0.5e-9
N_STEPS = 20000
OBS_FRAC = 0.5


@pytest.fixture(scope="module")
def binary(tmp_path_factory):
    return fdtd.build()


@pytest.mark.parametrize("ZL", [Z0, -1, 0, 25, 100])
def test_reflection_matches_analytic(binary, tmp_path, ZL):
    out_csv = tmp_path / "obs.csv"
    times, v_obs = fdtd.run(N, DX, N_STEPS, L, C, ZL, T0, SIGMA, OBS_FRAC,
                             out_csv, binary=binary)
    gamma, _, _, _, _ = fdtd.measure_reflection_coefficient(times, v_obs, SIGMA)
    expected = fdtd.analytic_gamma(ZL, Z0)
    assert abs(gamma - expected) < 0.02, f"ZL={ZL}: measured {gamma:.4f}, expected {expected:.4f}"


def test_discrete_causality_signal_cannot_arrive_early(binary, tmp_path):
    """This explicit local finite-difference stencil can only move
    information one cell per timestep, so a point far from the source MUST
    be exactly zero until at least that many steps have run, no matter how
    early the source's own (long, technically-infinite) Gaussian tail
    turns on. An earlier version of this test placed the observation point
    close to the source and asserted it should read ~0 several sigma
    before the pulse center; that failed, but the sim wasn't wrong, the
    test was: a near-source point genuinely does see the source's own
    small pre-peak tail (a few mV, not a bug), and that's not what
    causality actually forbids. What causality DOES guarantee is checked
    here instead: an exact zero before the light-cone step, not an
    approximately-small value before some arbitrary time."""
    out_csv = tmp_path / "obs.csv"
    obs_frac = 0.9  # far from the source
    v = 1.0 / np.sqrt(L * C)
    dt = 0.9 * DX / v  # matches the core's own Courant dt for ZL<0/open (see src)
    obs_cell = int(obs_frac * N)
    arrival_step = int(obs_cell * DX / v / dt)

    times, v_obs = fdtd.run(N, DX, arrival_step - 5, L, C, -1, T0, SIGMA, obs_frac,
                             out_csv, binary=binary)
    assert np.all(v_obs == 0.0)
