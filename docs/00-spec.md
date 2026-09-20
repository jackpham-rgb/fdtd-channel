# Spec

## What this is

An FDTD (Finite-Difference Time-Domain) electromagnetic/transmission-line
solver, built from Maxwell's/the telegrapher's equations, to watch a signal
physically propagate down an interconnect and extract quantities engineers
actually measure (reflection coefficient, S-parameters, insertion loss)
from the simulated fields. See the README for the three capabilities this
covers.

## Core language

C++17 for the time-stepping cores (`src/*.cpp`), matching this project's
own stated skill mapping (the same numerical-methods muscle as `gw_sim`).
Python (`py/`) builds each core, runs it, and does all analysis and
plotting; no physics or time-stepping happens in Python.

## Honesty line

An educational, from-first-principles FDTD field solver. Not a replacement
for commercial EM tools (HFSS/ADS). The value is that channel behavior is
derived from Maxwell's/the telegrapher's equations and correlated against
analytic theory.

## The ML / optimization decision

Same boundary as `serdes-link`: ML that designs or judges a circuit's
layout is out of scope here. Applied optimization on measurement/simulation
data is not, if a genuine use for it shows up (for example, fitting a
grid-convergence trend in a future extension).
