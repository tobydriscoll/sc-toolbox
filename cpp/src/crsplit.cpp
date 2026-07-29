#include "sctoolbox/crsplit.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

#include "sctoolbox/crcdt.hpp"
#include "sctoolbox/crpsdist.hpp"
#include "sctoolbox/scangle.hpp"

namespace sctoolbox {

namespace {

// Port of crsplit.m's nested crpsgd: geodesic distance from two point sets
// to a directed segment (points "inside" i.e. left of the segment use plain
// point-segment distance; points "outside" use distance to the nearer
// segment endpoint).
void crpsgd(std::complex<double> segA, std::complex<double> segB, const Eigen::VectorXcd& pts1,
            const Eigen::VectorXcd& pts2, Eigen::VectorXd& d1, Eigen::VectorXd& d2) {
    const int n1 = static_cast<int>(pts1.size());
    const int n2 = static_cast<int>(pts2.size());
    d1.resize(n1);
    d2.resize(n2);
    const std::complex<double> diffSeg = segB - segA;
    const std::complex<double> rot = std::conj(diffSeg / std::abs(diffSeg));
    const double epsD = std::numeric_limits<double>::epsilon();

    for (int i = 0; i < n1; ++i) {
        const std::complex<double> p = (pts1(i) - segA) * rot;
        if (p.imag() > epsD) {
            Eigen::VectorXcd one(1);
            one(0) = pts1(i);
            d1(i) = crpsdist({segA, segB}, one)(0);
        } else {
            d1(i) = std::abs(pts1(i) - segA);
        }
    }
    for (int i = 0; i < n2; ++i) {
        const std::complex<double> p = (pts2(i) - segA) * rot;
        if (p.imag() > epsD) {
            Eigen::VectorXcd one(1);
            one(0) = pts2(i);
            d2(i) = crpsdist({segA, segB}, one)(0);
        } else {
            d2(i) = std::abs(pts2(i) - segB);
        }
    }
}

}  // namespace

CrSplitResult crsplit(const Eigen::VectorXcd& wIn) {
    const int n0 = static_cast<int>(wIn.size());
    const Eigen::VectorXd beta0 = scangle(wIn);
    std::vector<bool> sharp0(n0);
    for (int i = 0; i < n0; ++i) sharp0[i] = beta0(i) < -0.75;

    const std::complex<double> nanc(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    Eigen::MatrixXcd neww = Eigen::MatrixXcd::Constant(3, n0, nanc);
    for (int j = 0; j < n0; ++j) neww(1, j) = wIn(j);

    for (int j = 0; j < n0; ++j) {
        if (!sharp0[j]) continue;
        const int jm1 = (j - 1 + n0) % n0;
        const int jp1 = (j + 1) % n0;

        Eigen::VectorXd ang(n0);
        for (int i = 0; i < n0; ++i) {
            const std::complex<double> v = (wIn(i) - wIn(j)) / (wIn(jp1) - wIn(j));
            ang(i) = std::fmod(std::arg(v) + 2.0 * M_PI, 2.0 * M_PI);
        }
        std::vector<bool> vis(n0);
        for (int i = 0; i < n0; ++i) vis[i] = (ang(i) > 0.0) && (ang(i) < ang(jm1));
        vis[jm1] = true;
        vis[j] = false;
        vis[jp1] = true;

        double mindist = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n0; ++i)
            if (vis[i]) mindist = std::min(mindist, std::abs(wIn(i) - wIn(j)));

        const std::complex<double> signJm1 = (wIn(jm1) - wIn(j)) / std::abs(wIn(jm1) - wIn(j));
        const std::complex<double> signJp1 = (wIn(jp1) - wIn(j)) / std::abs(wIn(jp1) - wIn(j));
        neww(0, j) = wIn(j) + 0.5 * mindist * signJm1;
        neww(2, j) = wIn(j) + 0.5 * mindist * signJp1;
    }

    std::vector<std::complex<double>> wVec;
    std::vector<bool> sharpVec, origVec;
    for (int j = 0; j < n0; ++j) {
        for (int r = 0; r < 3; ++r) {
            const std::complex<double> val = neww(r, j);
            if (std::isnan(val.real())) continue;
            wVec.push_back(val);
            sharpVec.push_back(r == 1 ? sharp0[j] : false);
            origVec.push_back(r == 1);
        }
    }

    Eigen::VectorXcd w(wVec.size());
    for (size_t i = 0; i < wVec.size(); ++i) w(i) = wVec[i];
    std::vector<bool> sharp = sharpVec;
    std::vector<bool> orig = origVec;

    CrTriangulation tri = crtriang(w);
    tri = crcdt(w, tri);

    // Phase 2: detect (but do not yet act on) narrow channels needing a
    // split. Build vertex adjacency from the current CDT.
    const int n = static_cast<int>(w.size());
    Eigen::VectorXd beta = scangle(w);
    Eigen::MatrixXi V = Eigen::MatrixXi::Zero(n, n);
    for (int c = 0; c < tri.edge.cols(); ++c) {
        const int a = tri.edge(0, c) - 1, b = tri.edge(1, c) - 1;
        V(a, b) = 1;
        V(b, a) = 1;
    }

    for (int j = 0; j < n; ++j) {
        const int jp1Sharp = (j + 1) % n;
        if (sharp[j] || sharp[jp1Sharp]) continue;

        const int jp1 = (j + 1) % n;
        const int jm1 = (j - 1 + n) % n;
        const int jp2 = (jp1 + 1) % n;

        std::vector<int> idx1, idx2;
        for (int i = 0; i < n; ++i) {
            if (i == j || i == jp1) continue;
            if (V(i, j) != 0) idx1.push_back(i);
        }
        if (std::abs(beta(j) - 1.0) < std::numeric_limits<double>::epsilon()) {
            idx1.erase(std::remove(idx1.begin(), idx1.end(), jm1), idx1.end());
        }
        for (int i = 0; i < n; ++i) {
            if (i == j || i == jp1) continue;
            if (V(i, jp1) != 0) idx2.push_back(i);
        }
        if (std::abs(beta(jp1) - 1.0) < std::numeric_limits<double>::epsilon()) {
            idx2.erase(std::remove(idx2.begin(), idx2.end(), jp2), idx2.end());
        }

        Eigen::VectorXcd w1(idx1.size()), w2(idx2.size());
        for (size_t i = 0; i < idx1.size(); ++i) w1(i) = w(idx1[i]);
        for (size_t i = 0; i < idx2.size(); ++i) w2(i) = w(idx2[i]);

        Eigen::VectorXd dist1, dist2;
        crpsgd(w(j), w(jp1), w1, w2, dist1, dist2);

        double mind = std::numeric_limits<double>::infinity();
        for (int i = 0; i < dist1.size(); ++i) mind = std::min(mind, dist1(i));
        for (int i = 0; i < dist2.size(); ++i) mind = std::min(mind, dist2(i));

        if (mind / std::abs(w(jp1) - w(j)) < 1.0 / (std::sqrt(2.0) * 3.0)) {
            throw std::runtime_error(
                "crsplit: polygon requires narrow-channel edge splitting, which is not "
                "implemented in this C++ port (no existing golden test exercises this path)");
        }
    }

    return CrSplitResult{w, orig, tri};
}

}  // namespace sctoolbox
