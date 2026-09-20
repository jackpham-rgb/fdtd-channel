// lossy_fdtd: the 1D telegrapher's-equation FDTD, extended with a uniform
// series resistance R (H/m -> ohm/m) and shunt conductance G (F/m -> S/m),
// so a signal actually attenuates and disperses like a real (lossy) trace
// instead of propagating forever unchanged.
//
//   dV/dx = -L dI/dt - R I
//   dI/dx = -C dV/dt - G V
//
// The loss terms are stepped semi-implicitly (Taflove's standard lossy-FDTD
// treatment: average the loss term between the old and new field value,
// same idea as a lossy Yee update's Ca/Cb coefficients). This keeps the
// scheme stable under the SAME Courant condition as the lossless case; a
// naive fully-explicit loss term would add its own (tighter) stability
// limit instead, the same kind of problem the reflection solver hit with the resistive
// load boundary.
//
// Both ends are matched to Z0 = sqrt(L/C), the LOSSLESS characteristic
// impedance. For a low-loss line (R << 2*pi*f*L, G << 2*pi*f*C -- true
// here by choice, see docs/03-loss.md) the real lossy impedance is close
// enough to Z0 that this is a good approximation, not an exact match; any
// small residual reflection this leaves shows up as a small, honestly-
// reported floor in the measurement, same as the other two solvers' boundaries.
//
// This solver ALWAYS uses a matched load (never open/short/other ZL): the
// point here is to measure how a signal decays and spreads while
// propagating, which needs the far end to absorb it cleanly, not reflect
// it back to contaminate the measurement.
//
// Two observation points (near the source, near the load) let Python
// measure the INSERTION LOSS directly, as the ratio of what's seen far
// vs. near: exactly the transfer function exp(-gamma(f) * distance).
//
// Usage:
//   lossy_fdtd N dx n_steps L C R G t0 sigma obsNear_frac obsFar_frac out.csv

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 13) {
        std::fprintf(stderr,
            "usage: %s N dx n_steps L C R G t0 sigma obsNear_frac obsFar_frac out.csv\n",
            argv[0]);
        return 1;
    }
    const int N           = std::atoi(argv[1]);
    const double dx       = std::atof(argv[2]);
    const int n_steps     = std::atoi(argv[3]);
    const double L        = std::atof(argv[4]);
    const double C        = std::atof(argv[5]);
    const double R        = std::atof(argv[6]);
    const double G        = std::atof(argv[7]);
    const double t0       = std::atof(argv[8]);
    const double sigma    = std::atof(argv[9]);
    const double obsNearFrac = std::atof(argv[10]);
    const double obsFarFrac  = std::atof(argv[11]);
    const std::string outfile = argv[12];

    const double v  = 1.0 / std::sqrt(L * C);
    const double Z0 = std::sqrt(L / C);

    double dt = 0.9 * dx / v;
    const double dtSourceLimit = 0.2 * C * dx * Z0;
    if (dtSourceLimit < dt) dt = dtSourceLimit;
    const double dtLoadLimit = 0.2 * C * dx * Z0;
    if (dtLoadLimit < dt) dt = dtLoadLimit;

    // Semi-implicit loss coefficients (interior nodes/cells).
    const double gHalf = G * dt / (2.0 * C);
    const double caV = (1.0 - gHalf) / (1.0 + gHalf);
    const double cbV = (dt / (C * dx)) / (1.0 + gHalf);
    const double rHalf = R * dt / (2.0 * L);
    const double caI = (1.0 - rHalf) / (1.0 + rHalf);
    const double cbI = (dt / (L * dx)) / (1.0 + rHalf);

    std::vector<double> V(N + 1, 0.0);
    std::vector<double> I(N, 0.0);

    const int obsNearIdx = static_cast<int>(obsNearFrac * N);
    const int obsFarIdx  = static_cast<int>(obsFarFrac * N);

    FILE* f = std::fopen(outfile.c_str(), "w");
    if (!f) {
        std::fprintf(stderr, "could not open %s for writing\n", outfile.c_str());
        return 1;
    }
    std::fprintf(f, "step,time,V_near,V_far\n");

    for (int n = 0; n < n_steps; ++n) {
        const double t = n * dt;

        // 1. Interior V, semi-implicit in G.
        for (int i = 1; i < N; ++i) {
            V[i] = caV * V[i] - cbV * (I[i] - I[i - 1]);
        }

        // 2. Matched source boundary (node 0), semi-implicit in G.
        const double Vs = std::exp(-((t - t0) * (t - t0)) / (2.0 * sigma * sigma));
        V[0] = caV * V[0] + 2.0 * cbV * ((Vs - V[0]) / Z0 - I[0]);

        // 3. Matched load boundary (node N), semi-implicit in G.
        V[N] = caV * V[N] + 2.0 * cbV * (I[N - 1] - V[N] / Z0);

        // 4. Current, semi-implicit in R.
        for (int i = 0; i < N; ++i) {
            I[i] = caI * I[i] - cbI * (V[i + 1] - V[i]);
        }

        std::fprintf(f, "%d,%.9e,%.9e,%.9e\n", n, t, V[obsNearIdx], V[obsFarIdx]);
    }

    std::fclose(f);
    std::fprintf(stdout, "dt=%.12e\n", dt);
    return 0;
}
