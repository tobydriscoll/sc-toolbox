#include "sctoolbox/ode45.hpp"

#include <algorithm>
#include <cmath>

namespace sctoolbox {

namespace {
// Dormand-Prince RK5(4) coefficients.
constexpr double kC[7] = {0.0, 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1.0, 1.0};
constexpr double kA[7][6] = {
    {0, 0, 0, 0, 0, 0},
    {1.0 / 5, 0, 0, 0, 0, 0},
    {3.0 / 40, 9.0 / 40, 0, 0, 0, 0},
    {44.0 / 45, -56.0 / 15, 32.0 / 9, 0, 0, 0},
    {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729, 0, 0},
    {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656, 0},
    {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84},
};
constexpr double kB5[7] = {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84, 0};
constexpr double kB4[7] = {5179.0 / 57600, 0, 7571.0 / 16695, 393.0 / 640, -92097.0 / 339200, 187.0 / 2100, 1.0 / 40};
}  // namespace

Eigen::VectorXd ode45(const OdeFun& f, double t0, double t1, const Eigen::VectorXd& y0, double abstol,
                       double reltol) {
    const double dir = (t1 >= t0) ? 1.0 : -1.0;
    double t = t0;
    Eigen::VectorXd y = y0;
    double h = dir * std::min(std::abs(t1 - t0) / 10.0, 0.1);
    if (h == 0.0) return y;

    Eigen::VectorXd k[7];

    int safety = 0;
    constexpr int kMaxSteps = 100000;
    while (dir * (t1 - t) > 0.0) {
        if (++safety > kMaxSteps) break;  // avoid hanging if NaN persists at every step size
        if (dir * (t + h - t1) > 0.0) h = t1 - t;

        k[0] = f(t, y);
        for (int s = 1; s < 7; ++s) {
            Eigen::VectorXd ys = y;
            for (int j = 0; j < s; ++j) ys += h * kA[s][j] * k[j];
            k[s] = f(t + kC[s] * h, ys);
        }

        Eigen::VectorXd y5 = y, y4 = y;
        for (int s = 0; s < 7; ++s) {
            y5 += h * kB5[s] * k[s];
            y4 += h * kB4[s] * k[s];
        }

        double errNorm = 0.0;
        for (int i = 0; i < y.size(); ++i) {
            const double sc = abstol + reltol * std::max(std::abs(y5(i)), std::abs(y(i)));
            errNorm = std::max(errNorm, std::abs(y5(i) - y4(i)) / sc);
        }

        const bool stepIsNan = std::isnan(errNorm);
        if (!stepIsNan && (errNorm <= 1.0 || std::abs(h) < 1e-14)) {
            t += h;
            y = y5;
        }

        // A NaN errNorm means the RHS produced NaN/Inf somewhere in this
        // trial step (e.g. a transient RK stage evaluation landed outside
        // the function's valid domain) -- treat it as a severe error and
        // shrink aggressively. Without this check, `errNorm > 0.0` is false
        // for NaN exactly like the "perfect, errNorm==0" case below, so the
        // step size would incorrectly *grow* (factor=5.0) right when it
        // needs to shrink, and (per the line above) a NaN step could even
        // get silently accepted once h underflows below 1e-14.
        double factor;
        if (stepIsNan) {
            factor = 0.2;
        } else if (errNorm > 0.0) {
            factor = 0.9 * std::pow(1.0 / errNorm, 0.2);
        } else {
            factor = 5.0;
        }
        factor = std::min(5.0, std::max(0.2, factor));
        h *= factor;
        if (dir * h <= 0.0) h = dir * 1e-10;
    }

    return y;
}

}  // namespace sctoolbox
