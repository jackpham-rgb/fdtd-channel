# fdtd-channel

An FDTD (Finite-Difference Time-Domain) transmission-line/EM solver, built
from the telegrapher's equations, to watch a signal physically propagate
down an interconnect and measure its reflection behavior from the simulated
fields.

Stage A (current): a 1D lossless line, a matched source, a resistive load,
and a reflection coefficient measured from the simulated waveform and
checked against the analytic formula.

## What's here

- `src/telegrapher_fdtd.cpp`: the time-stepping core (C++17). A Yee-staggered
  leapfrog solver for the 1D telegrapher's equations, with a matched
  (Thevenin) source and a resistive load boundary.
- `py/fdtd.py`: builds the core, runs it, loads its CSV output, and measures
  the reflection coefficient from the waveform.
- `py/run_telegrapher.py`: one command that regenerates every figure and the
  reflection-coefficient table.
- `tests/test_reflection.py`: checks measured vs. analytic reflection
  coefficient for 5 terminations, plus a discrete-causality test.
- `docs/00-spec.md`: what this project is and isn't, and the ML/optimization
  scope line (same as [serdes-link](https://github.com/jackpham-rgb/serdes-link)).
- `docs/01-telegrapher.md`: Stage A writeup, including the three real bugs
  hit while building this and how each was found.

## Requirements

- Python 3.10+ (numpy, scipy, matplotlib, pytest)
- g++ with C++17 support, on PATH

## Run it

```bash
python -m venv .venv
.venv\Scripts\activate   # or: source .venv/bin/activate on Linux/Mac
pip install numpy scipy matplotlib pytest

python py/run_telegrapher.py   # builds the core, runs all cases, writes docs/imgs/*.png
pytest tests/ -v
```

## Honesty line

An educational, from-first-principles FDTD field solver. Not a replacement
for commercial EM tools (HFSS/ADS). The value is that channel behavior is
derived from the telegrapher's equations and correlated against analytic
theory, not fit to match a target answer.
