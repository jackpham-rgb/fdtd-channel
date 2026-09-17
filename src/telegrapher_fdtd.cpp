// telegrapher_fdtd: 1D lossless telegrapher's-equation FDTD (Yee leapfrog).
//
// Solves dV/dx = -L dI/dt, dI/dx = -C dV/dt on a staggered grid: V lives on
// integer nodes 0..N, I lives on half-integer nodes 0..N-1 (I[i] is the
// current between V[i] and V[i+1]). A matched (Thevenin, Rs = Z0) source
// drives node 0, so it launches a clean outgoing pulse AND absorbs
// anything that returns to it, instead of re-reflecting it. Node N is
// terminated in a resistive load ZL (ZL < 0 means open circuit, ZL == 0
// means short).
//
// This is the compiled core; py/run_fdtd.py builds it, runs it, and does
// all the analysis and plotting. Keeping the physics here and the science
// in Python matches this project's own stated split (see the master
// planning reference).
//
// Usage:
//   telegrapher_fdtd N dx n_steps L C ZL t0 sigma obs_frac out.csv [snap.csv] [n_snapshots]
//     N            number of grid cells (V has N+1 nodes, I has N)
//     dx           cell size, meters
//     n_steps      number of time steps to run
//     L            inductance per unit length, H/m
//     C            capacitance per unit length, F/m
//     ZL           load resistance, ohms (negative = open circuit, 0 = short)
//     t0           source pulse center time, seconds
//     sigma        source pulse width (std dev), seconds
//     obs_frac     where to record V, as a fraction of the line (0..1)
//     out.csv      output path: step,time,V_obs
//     snap.csv     optional: full V(x) snapshots, one row per snapshot
//     n_snapshots  optional: how many evenly-time-spaced snapshots to write (default 0)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 11) {
        std::fprintf(stderr,
            "usage: %s N dx n_steps L C ZL t0 sigma obs_frac out.csv\n", argv[0]);
        return 1;
    }
    const int N          = std::atoi(argv[1]);
    const double dx      = std::atof(argv[2]);
    const int n_steps    = std::atoi(argv[3]);
    const double L       = std::atof(argv[4]);
    const double C       = std::atof(argv[5]);
    const double ZL      = std::atof(argv[6]);
    const double t0      = std::atof(argv[7]);
    const double sigma   = std::atof(argv[8]);
    const double obsFrac = std::atof(argv[9]);
    const std::string outfile = argv[10];
    const std::string snapfile = (argc > 11) ? argv[11] : "";
    const int nSnapshots = (argc > 12) ? std::atoi(argv[12]) : 0;

    const double v  = 1.0 / std::sqrt(L * C);
    const double Z0 = std::sqrt(L / C);
    double dt = 0.9 * dx / v;  // Courant limit for the interior update: v*dt/dx <= 1

    // Both the source and the load boundaries are lumped half-cell
    // capacitors (C*dx/2) discharging through a resistor (Rs at the
    // source, ZL at the load), each its own explicit-Euler ODE:
    //   dV/dt = (2/(C*dx)) * (I_in - V/R)
    // which is only stable for dt <= C*dx*R. For a low-impedance
    // termination this can be TIGHTER than the interior Courant limit;
    // take whichever limit is smallest. Found by the sim blowing up (NaN)
    // at moderate ZL with only the interior Courant check applied; see
    // docs/01-telegrapher.md.
    const double dtSourceLimit = 0.2 * C * dx * Z0;
    if (dtSourceLimit < dt) dt = dtSourceLimit;
    if (ZL > 0.0) {
        const double dtLoadLimit = 0.2 * C * dx * ZL;
        if (dtLoadLimit < dt) dt = dtLoadLimit;
    }

    std::vector<double> V(N + 1, 0.0);
    std::vector<double> I(N, 0.0);

    const int obsIdx = static_cast<int>(obsFrac * N);

    FILE* f = std::fopen(outfile.c_str(), "w");
    if (!f) {
        std::fprintf(stderr, "could not open %s for writing\n", outfile.c_str());
        return 1;
    }
    std::fprintf(f, "step,time,V_obs\n");

    FILE* snapF = nullptr;
    std::vector<int> snapSteps;
    if (nSnapshots > 0 && !snapfile.empty()) {
        snapF = std::fopen(snapfile.c_str(), "w");
        if (!snapF) {
            std::fprintf(stderr, "could not open %s for writing\n", snapfile.c_str());
            return 1;
        }
        std::fprintf(snapF, "time");
        for (int i = 0; i <= N; ++i) std::fprintf(snapF, ",x%d", i);
        std::fprintf(snapF, "\n");
        for (int s = 0; s < nSnapshots; ++s) {
            snapSteps.push_back(static_cast<int>((double)(s + 1) / (nSnapshots + 1) * n_steps));
        }
    }
    size_t nextSnap = 0;

    for (int n = 0; n < n_steps; ++n) {
        const double t = n * dt;

        // 1. Interior V nodes, using the OLD current (leapfrog: V then I).
        for (int i = 1; i < N; ++i) {
            V[i] -= (dt / (C * dx)) * (I[i] - I[i - 1]);
        }

        // 2. Matched (Thevenin) source at node 0: an ideal source Vs(t) in
        //    series with Rs = Z0, so it both launches a clean pulse and
        //    absorbs anything that comes back (no re-reflection at the
        //    source). An earlier version just added the raw pulse to V[0]
        //    directly; that both re-reflected returning waves AND scaled
        //    with 1/dt (more timesteps under the same pulse width injected
        //    proportionally more total "charge"), and blew up. See
        //    docs/01-telegrapher.md.
        const double Vs = std::exp(-((t - t0) * (t - t0)) / (2.0 * sigma * sigma));
        V[0] += (2.0 * dt / (C * dx)) * ((Vs - V[0]) / Z0 - I[0]);

        // 3. Load boundary at node N: same half-cell treatment through ZL.
        //    ZL<0 is open (no load current); ZL==0 is a short.
        if (ZL < 0.0) {
            V[N] += (2.0 * dt / (C * dx)) * I[N - 1];
        } else if (ZL == 0.0) {
            V[N] = 0.0;
        } else {
            V[N] += (2.0 * dt / (C * dx)) * (I[N - 1] - V[N] / ZL);
        }

        // 4. Current, using the NEW voltage.
        for (int i = 0; i < N; ++i) {
            I[i] -= (dt / (L * dx)) * (V[i + 1] - V[i]);
        }

        std::fprintf(f, "%d,%.9e,%.9e\n", n, t, V[obsIdx]);

        if (snapF && nextSnap < snapSteps.size() && n == snapSteps[nextSnap]) {
            std::fprintf(snapF, "%.9e", t);
            for (int i = 0; i <= N; ++i) std::fprintf(snapF, ",%.9e", V[i]);
            std::fprintf(snapF, "\n");
            ++nextSnap;
        }
    }

    std::fclose(f);
    if (snapF) std::fclose(snapF);
    return 0;
}
