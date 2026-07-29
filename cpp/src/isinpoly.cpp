#include "sctoolbox/isinpoly.hpp"

#include <cmath>
#include <complex>

#include "sctoolbox/scangle.hpp"

namespace sctoolbox {

namespace {
std::complex<double> csign(std::complex<double> v) {
    const double a = std::abs(v);
    return a == 0.0 ? std::complex<double>(0.0, 0.0) : v / a;
}
}  // namespace

Eigen::VectorXd isinpoly(const Eigen::VectorXcd& z, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                          double tol) {
    const int n = static_cast<int>(w.size());
    const int np = static_cast<int>(z.size());
    Eigen::VectorXd index = Eigen::VectorXd::Zero(np);
    const double eps_ = std::numeric_limits<double>::epsilon();

    double scale = 0.0;
    for (int i = 0; i < n; ++i) scale += std::abs(w((i + 1) % n) - w(i));
    scale /= n;
    if (!(scale > eps_)) return index;

    const Eigen::VectorXcd ws = w / scale;
    const Eigen::VectorXcd zs = z / scale;

    Eigen::MatrixXcd d(n, np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) {
            std::complex<double> dij = ws(i) - zs(j);
            if (std::abs(dij) < eps_) dij = std::complex<double>(eps_, 0.0);
            d(i, j) = dij;
        }
    }

    Eigen::MatrixXd ang(n, np);
    for (int i = 0; i < n; ++i) {
        const int ni = (i + 1) % n;
        for (int j = 0; j < np; ++j) ang(i, j) = std::arg(d(ni, j) / d(i, j)) / M_PI;
    }

    Eigen::VectorXcd tangents(n);
    for (int i = 0; i < n; ++i) tangents(i) = csign(ws((i + 1) % n) - ws(i));
    for (int p = 0; p < n; ++p) {
        if (tangents(p) == std::complex<double>(0.0, 0.0)) {
            for (int step = 1; step <= n; ++step) {
                const int idx = (p + step) % n;
                if (ws(idx) != ws(p)) {
                    tangents(p) = csign(ws(idx) - ws(p));
                    break;
                }
            }
        }
    }

    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> onbdy(n, np), onvtx(n, np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) {
            const std::complex<double> rot = d(i, j) / tangents(i);
            onvtx(i, j) = std::abs(d(i, j)) < tol;
            onbdy(i, j) = std::abs(rot.imag()) < 10.0 * tol;
        }
    }
    for (int i = 0; i < n; ++i) {
        const int ni = (i + 1) % n;
        for (int j = 0; j < np; ++j) {
            const bool onSegment = (std::abs(ang(i, j)) > 0.9) || onvtx(i, j) || onvtx(ni, j);
            onbdy(i, j) = onbdy(i, j) && onSegment;
        }
    }

    for (int j = 0; j < np; ++j) {
        bool anyOnbdy = false;
        for (int i = 0; i < n; ++i) {
            if (onbdy(i, j)) {
                anyOnbdy = true;
                break;
            }
        }
        if (!anyOnbdy) {
            double s = 0.0;
            for (int i = 0; i < n; ++i) s += ang(i, j);
            index(j) = std::round(s / 2.0);
        } else {
            double S = 0.0;
            for (int i = 0; i < n; ++i)
                if (!onbdy(i, j)) S += ang(i, j);
            double bsum = 0.0;
            int onbdyCount = 0, onvtxCount = 0;
            for (int i = 0; i < n; ++i) {
                if (onbdy(i, j)) ++onbdyCount;
                if (onvtx(i, j)) {
                    ++onvtxCount;
                    bsum += beta(i);
                }
            }
            const double augment = onbdyCount - onvtxCount - bsum;
            const double signS = (S > 0.0) - (S < 0.0);
            index(j) = std::round(augment * signS + S) / 2.0;
        }
    }

    for (int j = 0; j < np; ++j) index(j) = (index(j) != 0.0) ? 1.0 : 0.0;
    return index;
}

Eigen::VectorXd isinpoly(const Eigen::VectorXcd& z, const Eigen::VectorXcd& w) {
    return isinpoly(z, w, scangle(w), std::numeric_limits<double>::epsilon());
}

Eigen::VectorXd isinpoly(const Eigen::VectorXcd& z, const Eigen::VectorXcd& w, double tol) {
    return isinpoly(z, w, scangle(w), tol);
}

}  // namespace sctoolbox
