#include "sctoolbox/scaddvtx.hpp"

#include <cmath>
#include <complex>
#include <vector>

namespace sctoolbox {

namespace {
bool isInf(const std::complex<double>& z) { return std::isinf(z.real()) || std::isinf(z.imag()); }
}  // namespace

ScaddvtxResult scaddvtx(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, int pos,
                         std::array<double, 4> window) {
    const int n = static_cast<int>(w.size());
    const int pos1 = (pos + 1) % n;
    std::complex<double> newv;

    if (!isInf(w(pos)) && !isInf(w(pos1))) {
        newv = 0.5 * (w(pos) + w(pos1));
    } else {
        // Find a pair of adjacent finite vertices as a basis for absolute angle.
        std::vector<double> ang(n, 0.0);
        int base = -1;
        for (int i = 0; i < n; ++i) {
            if (!isInf(w(i)) && !isInf(w((i + 1) % n))) {
                base = i;
                break;
            }
        }
        ang[base] = std::arg(w((base + 1) % n) - w(base));

        for (int step = 1; step < n; ++step) {
            const int j = (base + step) % n;
            const int prev = (j - 1 + n) % n;
            ang[j] = ang[prev] - M_PI * beta(j);
            if (j == pos) break;
        }

        double lensum = 0.0;
        int lencount = 0;
        for (int i = 0; i < n; ++i) {
            const double l = std::abs(w((i + 1) % n) - w(i));
            if (!std::isinf(l)) {
                lensum += l;
                ++lencount;
            }
        }
        const double avglen0 = lensum / lencount;

        std::complex<double> basept;
        double dirang;
        if (isInf(w(pos))) {
            basept = w(pos1);
            dirang = ang[pos] + M_PI;
        } else {
            basept = w(pos);
            dirang = ang[pos];
        }
        const std::complex<double> dir = std::polar(1.0, dirang);

        double avglen = avglen0;
        newv = basept + avglen * dir;
        while (newv.real() < window[0] || newv.real() > window[1] || newv.imag() < window[2] ||
               newv.imag() > window[3]) {
            avglen /= 2.0;
            newv = basept + avglen * dir;
        }
    }

    Eigen::VectorXcd wn(n + 1);
    Eigen::VectorXd betan(n + 1);
    for (int i = 0; i <= pos; ++i) {
        wn(i) = w(i);
        betan(i) = beta(i);
    }
    wn(pos + 1) = newv;
    betan(pos + 1) = 0.0;
    for (int i = pos + 1; i < n; ++i) {
        wn(i + 1) = w(i);
        betan(i + 1) = beta(i);
    }

    return {wn, betan};
}

}  // namespace sctoolbox
