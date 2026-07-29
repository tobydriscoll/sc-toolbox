#include "sctoolbox/rcorners.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sctoolbox {

RCorners rcorners(const Eigen::VectorXcd& wIn, const Eigen::VectorXd& betaIn, const Eigen::VectorXcd& zIn) {
    const int n = static_cast<int>(wIn.size());
    const double eps = std::numeric_limits<double>::epsilon();

    double minRe = zIn(0).real(), maxRe = zIn(0).real();
    double minIm = zIn(0).imag(), maxIm = zIn(0).imag();
    for (int i = 1; i < n; ++i) {
        minRe = std::min(minRe, zIn(i).real());
        maxRe = std::max(maxRe, zIn(i).real());
        minIm = std::min(minIm, zIn(i).imag());
        maxIm = std::max(maxIm, zIn(i).imag());
    }

    // 0-indexed corner candidates (more than one boundary condition holds).
    std::vector<int> corners;  // 0-indexed
    for (int i = 0; i < n; ++i) {
        const bool left = std::abs(zIn(i).real() - minRe) < eps;
        const bool right = std::abs(zIn(i).real() - maxRe) < eps;
        const bool top = std::abs(zIn(i).imag() - maxIm) < eps;
        const bool bot = std::abs(zIn(i).imag() - minIm) < eps;
        const int sum = (left ? 1 : 0) + (right ? 1 : 0) + (top ? 1 : 0) + (bot ? 1 : 0);
        if (sum - 1 > 0) corners.push_back(i);
    }

    int c1 = 0;
    {
        double best = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            const double d = std::abs(zIn(i) - std::complex<double>(maxRe, 0.0));
            if (d < best) {
                best = d;
                c1 = i;
            }
        }
    }

    int offset = 0;
    for (size_t k = 0; k < corners.size(); ++k) {
        if (corners[k] == c1) {
            offset = static_cast<int>(k);
            break;
        }
    }

    std::vector<int> cornersRot;
    for (size_t k = offset; k < corners.size(); ++k) cornersRot.push_back(corners[k]);
    for (int k = 0; k < offset; ++k) cornersRot.push_back(corners[k]);
    corners = cornersRot;

    // renum = [corners(1):n, 1:corners(1)-1] in MATLAB's 1-indexed terms;
    // here corners[0] is 0-indexed, so renum starts at corners[0].
    const int start = corners[0];
    std::vector<int> renum;
    for (int i = start; i < n; ++i) renum.push_back(i);
    for (int i = 0; i < start; ++i) renum.push_back(i);

    RCorners r;
    r.w.resize(n);
    r.beta.resize(n);
    r.z.resize(n);
    for (int i = 0; i < n; ++i) {
        r.w(i) = wIn(renum[i]);
        r.beta(i) = betaIn(renum[i]);
        r.z(i) = zIn(renum[i]);
    }

    // corners = rem(corners - corners(1) + 1 + n - 1, n) + 1  (MATLAB, 1-indexed)
    // In 0-indexed terms this is simply (corners[k] - start) mod n.
    r.corners.resize(corners.size());
    for (size_t k = 0; k < corners.size(); ++k) {
        int v = (corners[k] - start) % n;
        if (v < 0) v += n;
        r.corners[k] = v;
    }

    return r;
}

}  // namespace sctoolbox
