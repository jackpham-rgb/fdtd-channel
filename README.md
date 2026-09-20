# fdtd-channel

An FDTD (Finite-Difference Time-Domain) transmission-line/EM solver, built
from the telegrapher's equations, to watch a signal physically propagate
down an interconnect. Three things it does, each checked against analytic
theory instead of just plotted and trusted:

- **Reflection**: a matched source and a resistive load on a 1D lossless
  line; the reflection coefficient measured from the simulated waveform
  matches the analytic formula.
- **S-parameters**: an impedance-step discontinuity turned into a virtual
  network analyzer. S11(f)/S21(f) are extracted from the simulated fields
  via time-domain gating and match the analytic two-port (ABCD-matrix)
  result.
- **Loss & dispersion**: series resistance and shunt conductance added to
  the line, so insertion loss actually rolls off with frequency instead of
  staying flat, matching the analytic lossy propagation constant.

## What's here

- `src/telegrapher_fdtd.cpp` / `py/fdtd.py` / `py/run_telegrapher.py`:
  the reflection solver. A Yee-staggered leapfrog core (C++17) with a
  matched (Thevenin) source and a resistive load boundary; the Python
  layer builds it, runs it, and measures the reflection coefficient.
  `tests/test_reflection.py` checks 5 terminations plus a
  discrete-causality test.
- `src/twoport_fdtd.cpp` / `py/twoport.py` / `py/run_sparams.py`:
  the S-parameter (virtual VNA) solver. The same core scheme, generalized
  to position-dependent L(x)/C(x) so an impedance slab can sit in the
  line; separates incident/reflected/transmitted waves by time-domain
  gating and computes S11(f)/S21(f). `tests/test_sparams.py` checks 3
  (impedance, slab length) combinations.
- `src/lossy_fdtd.cpp` / `py/lossy.py` / `py/run_loss.py`: the lossy-line
  solver. Adds series resistance R and shunt conductance G with a
  semi-implicit loss update (so it doesn't add its own stability limit),
  and measures the transfer function between two points against the
  analytic lossy propagation constant. `tests/test_loss.py` checks 3
  (R, G) combinations plus a lossless-limit sanity check.
- `docs/00-spec.md`: what this project is and isn't, and the ML/optimization
  scope line (same as [serdes-link](https://github.com/jackpham-rgb/serdes-link)).
- `docs/01-telegrapher.md`, `docs/02-sparameters.md`, `docs/03-loss.md`:
  a writeup per capability above, each including the real bugs hit while
  building it and how they were found.

## Requirements

- Python 3.10+ (numpy, scipy, matplotlib, pytest)
- g++ with C++17 support, on PATH

## Run it

```bash
python -m venv .venv
.venv\Scripts\activate   # or: source .venv/bin/activate on Linux/Mac
pip install numpy scipy matplotlib pytest

python py/run_telegrapher.py   # reflection: builds the core, runs all cases, writes docs/imgs/*.png
python py/run_sparams.py       # S-parameters: builds the core, runs the VNA sweep, writes docs/imgs/*.png
python py/run_loss.py          # loss & dispersion: builds the core, runs the loss sweep, writes docs/imgs/*.png
pytest tests/ -v
```

## Honesty line

An educational, from-first-principles FDTD field solver. Not a replacement
for commercial EM tools (HFSS/ADS). The value is that channel behavior is
derived from the telegrapher's equations and correlated against analytic
theory, not fit to match a target answer.
