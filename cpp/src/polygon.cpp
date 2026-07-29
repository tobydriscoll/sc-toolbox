#include "sctoolbox/polygon.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "sctoolbox/isinpoly.hpp"

namespace sctoolbox {

namespace {

const double kEps = std::numeric_limits<double>::epsilon();

bool isInfC(const std::complex<double>& z) { return std::isinf(z.real()) || std::isinf(z.imag()); }

struct AngleResult {
    Eigen::VectorXd alpha;
    bool isccw;
    int index;
};

// Port of @polygon/angle.m.
AngleResult computeAngleInfo(const Eigen::VectorXcd& w, const Eigen::VectorXd& providedAlpha) {
    const int n = static_cast<int>(w.size());
    Eigen::VectorXd alpha;

    if (providedAlpha.size() > 0) {
        alpha = providedAlpha;
    } else {
        for (int i = 0; i < n; ++i)
            if (isInfC(w(i))) throw std::runtime_error("Cannot compute angles for unbounded polygons.");

        alpha.resize(n);
        for (int i = 0; i < n; ++i) {
            const std::complex<double> incoming = w(i) - w((i - 1 + n) % n);
            const std::complex<double> outgoing = w((i + 1) % n) - w(i);
            double a = std::arg(-incoming * std::conj(outgoing)) / M_PI;
            a = std::fmod(a, 2.0);
            if (a < 0.0) a += 2.0;
            alpha(i) = a;
        }

        std::vector<bool> mask(n);
        bool allMask = true;
        for (int i = 0; i < n; ++i) {
            mask[i] = (alpha(i) < 100.0 * kEps) || (2.0 - alpha(i) < 100.0 * kEps);
            if (!mask[i]) allMask = false;
        }
        if (allMask) {
            // All vertices collinear -- degenerate polygon. MATLAB returns
            // early here, bypassing the index/isccw computation below.
            alpha.setZero();
            return {alpha, true, 1};
        }

        std::vector<std::complex<double>> suspicious, rest;
        std::vector<int> suspiciousIdx;
        for (int i = 0; i < n; ++i) {
            if (mask[i]) {
                suspicious.push_back(w(i));
                suspiciousIdx.push_back(i);
            } else {
                rest.push_back(w(i));
            }
        }
        Eigen::VectorXcd zmask(suspicious.size()), wsub(rest.size());
        for (size_t i = 0; i < suspicious.size(); ++i) zmask(i) = suspicious[i];
        for (size_t i = 0; i < rest.size(); ++i) wsub(i) = rest[i];
        const Eigen::VectorXd slit = isinpoly(zmask, wsub);
        for (size_t i = 0; i < suspiciousIdx.size(); ++i) alpha(suspiciousIdx[i]) = (slit(i) != 0.0) ? 2.0 : 0.0;
    }

    auto computeIndex = [&]() {
        double idx = 0.0;
        for (int i = 0; i < n; ++i) idx += (alpha(i) - 1.0);
        return idx / 2.0;
    };

    double index = computeIndex();
    if (std::abs(index - std::round(index)) > 100.0 * std::sqrt(static_cast<double>(n)) * kEps) {
        for (int i = 0; i < n; ++i) {
            const bool nearDegenerate = (alpha(i) < 2.0 * kEps) || (2.0 - alpha(i) < 2.0 * kEps);
            if (!nearDegenerate) alpha(i) = 2.0 - alpha(i);
        }
        index = computeIndex();
        if (std::abs(index - std::round(index)) > 100.0 * std::sqrt(static_cast<double>(n)) * kEps)
            throw std::runtime_error("Invalid polygon.");
    }

    const int idxRounded = static_cast<int>(std::lround(index));
    return {alpha, idxRounded < 0, idxRounded};
}

}  // namespace

Polygon::Polygon(Eigen::VectorXcd vertices) { init(std::move(vertices), Eigen::VectorXd()); }

Polygon::Polygon(Eigen::VectorXcd vertices, Eigen::VectorXd angles) {
    init(std::move(vertices), std::move(angles));
}

void Polygon::init(Eigen::VectorXcd w, Eigen::VectorXd alpha) {
    int n0 = static_cast<int>(w.size());
    if (n0 > 0 && std::abs(w(n0 - 1) - w(0)) < 3.0 * kEps) {
        w.conservativeResize(n0 - 1);
        n0 -= 1;
    }
    vertex_ = w;
    angle_ = alpha;

    if (n0 == 0) return;

    const AngleResult r = computeAngleInfo(vertex_, angle_);
    if (!r.isccw) {
        vertex_ = vertex_.reverse().eval();
        angle_ = (Eigen::VectorXd::Constant(n0, 2.0) - r.alpha.reverse().eval()).eval();
    } else {
        angle_ = r.alpha;
    }
    // abs(index) > 1 would warn "Polygon is multiple-sheeted" in MATLAB;
    // not fatal, so no action needed here.
}

bool Polygon::isInf() const {
    for (int i = 0; i < vertex_.size(); ++i)
        if (isInfC(vertex_(i))) return true;
    return false;
}

}  // namespace sctoolbox
