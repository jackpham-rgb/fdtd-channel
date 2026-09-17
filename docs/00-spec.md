# Stage 0: Spec

## What this is

An FDTD (Finite-Difference Time-Domain) electromagnetic/transmission-line
solver, built from Maxwell's/the telegrapher's equations, to watch a signal
physically propagate down an interconnect and extract its S-parameters from
the simulated fields. Stage A (this stage) is the 1D lossless telegrapher's
equation only: a source, a line, a load, and a reflection coefficient.

## Core language

C++17 for the time-stepping core (`src/telegrapher_fdtd.cpp`), matching
this project's own stated skill mapping (the same numerical-methods muscle
as `gw_sim`). Python (`py/`) builds the core, runs it, and does all
analysis and plotting; no physics or time-stepping happens in Python.

## Honesty line

An educational, from-first-principles FDTD field solver. Not a replacement
for commercial EM tools (HFSS/ADS). The value is that channel behavior is
derived from Maxwell's/the telegrapher's equations and correlated against
analytic theory (and, in later stages, LTspice and a measured channel).

## The ML / optimization decision

Same boundary as `serdes-link`: ML that designs or judges a circuit's
layout is out of scope here. Applied optimization on measurement/simulation
data is not, if a genuine use for it shows up in later stages (for example,
fitting a grid-convergence trend in Stage E).
