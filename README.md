# fdtd-channel

An FDTD (Finite-Difference Time-Domain) transmission-line/EM solver, built
from the telegrapher's equations, to watch a signal physically propagate
down an interconnect and measure its reflection behavior from the simulated
fields.

- **Stage A**: a 1D lossless line, a matched source, a resistive load, and a
  reflection coefficient measured from the simulated waveform and checked
  against the analytic formula.
- **Stage B**: an impedance-step discontinuity turned into a virtual
  network analyzer: S11(f)/S21(f) extracted from the simulated fields via
  time-domain gating, checked against the analytic two-port (ABCD-matrix)
  result.
- **Stage C** (current): series resistance R and shunt conductance G added
  to the line, so insertion loss actually rolls off with frequency instead
  of staying flat, checked against the analytic lossy propagation
  constant.

## What's here

- `src/telegrapher_fdtd.cpp`: Stage A's time-stepping core (C++17). A
  Yee-staggered leapfrog solver for the 1D telegrapher's equations, with a
  matched (Thevenin) source and a resistive load boundary.
- `py/fdtd.py`: builds Stage A's core, runs it, loads its CSV output, and
  measures the reflection coefficient from the waveform.
- `py/run_telegrapher.py`: one command that regenerates every Stage A
  figure and the reflection-coefficient table.
- `tests/test_reflection.py`: checks measured vs. analytic reflection
  coefficient for 5 terminations, plus a discrete-causality test.
- `src/twoport_fdtd.cpp`: Stage B's core. The same scheme as Stage A,
  generalized to position-dependent L(x)/C(x) so a slab of a different
  impedance can sit in the middle of the line.
- `py/twoport.py`: builds Stage B's core, runs it, separates
  incident/reflected/transmitted waves by time-domain gating, and computes
  S11(f)/S21(f) (with an analytic ABCD-matrix cross-check).
- `py/run_sparams.py`: one command that regenerates Stage B's S-parameter
  and time-domain figures.
- `tests/test_sparams.py`: checks measured vs. analytic S11/S21 for 3
  (impedance, slab length) combinations, plus a regression test for a
  windowing bug this stage hit.
- `docs/00-spec.md`: what this project is and isn't, and the ML/optimization
  scope line (same as [serdes-link](https://github.com/jackpham-rgb/serdes-link)).
- `docs/01-telegrapher.md`: Stage A writeup, including the three real bugs
  hit while building this and how each was found.
- `docs/02-sparameters.md`: Stage B writeup, including the multi-bounce
  amplitude mistake and the windowing bug, and how each was found.
- `src/lossy_fdtd.cpp`: Stage C's core. Adds series resistance R and shunt
  conductance G to Stage A's scheme, with a semi-implicit loss update so it
  doesn't add its own stability limit.
- `py/lossy.py`: builds Stage C's core, runs it, and computes the measured
  transfer function between two points (with an analytic lossy-line
  cross-check).
- `py/run_loss.py`: one command that regenerates Stage C's loss-vs-frequency
  and pulse-dispersion figures.
- `tests/test_loss.py`: checks measured vs. analytic attenuation for 3
  (R, G) combinations, plus a lossless-limit sanity check.
- `docs/03-loss.md`: Stage C writeup, including why constant R/G still
  rolls off with frequency, and the low-frequency resolution limitation
  found while validating it.

## Requirements

- Python 3.10+ (numpy, scipy, matplotlib, pytest)
- g++ with C++17 support, on PATH

## Run it

```bash
python -m venv .venv
.venv\Scripts\activate   # or: source .venv/bin/activate on Linux/Mac
pip install numpy scipy matplotlib pytest

python py/run_telegrapher.py   # Stage A: builds the core, runs all cases, writes docs/imgs/*.png
python py/run_sparams.py       # Stage B: builds the core, runs the VNA sweep, writes docs/imgs/*.png
python py/run_loss.py          # Stage C: builds the core, runs the loss sweep, writes docs/imgs/*.png
pytest tests/ -v
```

## Honesty line

An educational, from-first-principles FDTD field solver. Not a replacement
for commercial EM tools (HFSS/ADS). The value is that channel behavior is
derived from the telegrapher's equations and correlated against analytic
theory, not fit to match a target answer.
