// twoport_fdtd: same 1D telegrapher's-equation FDTD as telegrapher_fdtd, but
// with a middle line section of a DIFFERENT characteristic impedance (an
// impedance-step discontinuity), used as a virtual two-port network.
//
// Layout: region 1 = [0, x1), region 2 = [x1, x2), region 3 = [x2, N*dx].
// Regions 1 and 3 share impedance Z1 = sqrt(L1/C1); region 2 has impedance
// Zmid = sqrt(L2/C2). The propagation velocity v = 1/sqrt(L*C) is kept the
// SAME in all three regions (only the L/C ratio changes), so any reflection
// is purely from the impedance step, not a change in wave speed. Both ends
// (node 0 and node N) are matched (Thevenin, Rs = Z1), so the ONLY
// reflections in the whole line come from the two impedance steps around
// region 2, exactly like a real via or connector transition.
//
// Two observation points: obs1 in region 1 (the port-1 reference plane,
// where py/twoport.py separates incident from reflected), obs2 in region 3
// (the port-2 reference plane, where the transmitted wave is read directly
// since the far end is matched and nothing reflects back into region 3).
//
// Usage:
//   twoport_fdtd N dx n_steps L1 C1 Zmid x1_frac x2_frac t0 sigma
//                obs1_frac obs2_frac out.csv
//
// Prints "dt=<value>" to stdout so Python doesn't have to recompute it.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 14) {
        std::fprintf(stderr,
            "usage: %s N dx n_steps L1 C1 Zmid x1_frac x2_frac t0 sigma "
            "obs1_frac obs2_frac out.csv\n", argv[0]);
        return 1;
    }
    const int N          = std::atoi(argv[1]);
    const double dx      = std::atof(argv[2]);
    const int n_steps    = std::atoi(argv[3]);
    const double L1      = std::atof(argv[4]);
    const double C1      = std::atof(argv[5]);
    const double Zmid    = std::atof(argv[6]);
    const double x1Frac  = std::atof(argv[7]);
    const double x2Frac  = std::atof(argv[8]);
    const double t0      = std::atof(argv[9]);
    const double sigma   = std::atof(argv[10]);
    const double obs1Frac = std::atof(argv[11]);
    const double obs2Frac = std::atof(argv[12]);
    const std::string outfile = argv[13];

    const double v  = 1.0 / std::sqrt(L1 * C1);   // same velocity everywhere
    const double Z1 = std::sqrt(L1 / C1);
    const double L2 = Zmid / v;
    const double C2 = 1.0 / (Zmid * v);

    double dt = 0.9 * dx / v;  // Courant limit (v is the same in every region)
    const double dtSourceLimit = 0.2 * C1 * dx * Z1;
    if (dtSourceLimit < dt) dt = dtSourceLimit;
    const double dtLoadLimit = 0.2 * C1 * dx * Z1;  // node N sits in region 3 (Z1)
    if (dtLoadLimit < dt) dt = dtLoadLimit;

    // Per-position L (attached to I[i], the cell between V[i] and V[i+1])
    // and C (attached to V[i]). A node/cell belongs to region 2 if its
    // position falls in [x1, x2); ties at the exact boundary go to the
    // region on the right, which only matters in the dx -> 0 limit.
    const double x1 = x1Frac * N * dx;
    const double x2 = x2Frac * N * dx;

    std::vector<double> Lc(N), Cc(N + 1);
    for (int i = 0; i < N; ++i) {
        const double xCell = (i + 0.5) * dx;
        Lc[i] = (xCell >= x1 && xCell < x2) ? L2 : L1;
    }
    for (int i = 0; i <= N; ++i) {
        const double xNode = i * dx;
        Cc[i] = (xNode >= x1 && xNode < x2) ? C2 : C1;
    }

    std::vector<double> V(N + 1, 0.0);
    std::vector<double> I(N, 0.0);

    const int obs1Idx = static_cast<int>(obs1Frac * N);
    const int obs2Idx = static_cast<int>(obs2Frac * N);

    FILE* f = std::fopen(outfile.c_str(), "w");
    if (!f) {
        std::fprintf(stderr, "could not open %s for writing\n", outfile.c_str());
        return 1;
    }
    std::fprintf(f, "step,time,V_obs1,V_obs2\n");

    for (int n = 0; n < n_steps; ++n) {
        const double t = n * dt;

        for (int i = 1; i < N; ++i) {
            V[i] -= (dt / (Cc[i] * dx)) * (I[i] - I[i - 1]);
        }

        const double Vs = std::exp(-((t - t0) * (t - t0)) / (2.0 * sigma * sigma));
        V[0] += (2.0 * dt / (Cc[0] * dx)) * ((Vs - V[0]) / Z1 - I[0]);
        V[N] += (2.0 * dt / (Cc[N] * dx)) * (I[N - 1] - V[N] / Z1);

        for (int i = 0; i < N; ++i) {
            I[i] -= (dt / (Lc[i] * dx)) * (V[i + 1] - V[i]);
        }

        std::fprintf(f, "%d,%.9e,%.9e,%.9e\n", n, t, V[obs1Idx], V[obs2Idx]);
    }

    std::fclose(f);
    std::fprintf(stdout, "dt=%.12e\n", dt);
    return 0;
}
