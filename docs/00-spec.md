# Spec

## What this is

A field solver for a transmission line, built from the telegrapher's equations with the finite-difference time-domain (FDTD) method. It lets a signal propagate down an interconnect and extracts the quantities engineers measure in the lab: reflection coefficient, S-parameters and insertion loss. See the README for the three solvers.

## Language split

C++17 for the time-stepping cores (`src/*.cpp`). Python (`py/`) builds each core, runs it, and does all analysis and plotting. No time-stepping happens in Python.

## Scope

A from-first-principles solver for learning and validation, not a replacement for commercial EM tools such as HFSS or ADS. Every result is compared against the analytic answer rather than tuned to match it.
