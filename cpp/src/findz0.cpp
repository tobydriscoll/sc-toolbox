#include "sctoolbox/findz0.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace sctoolbox {

namespace {
bool isInf(std::complex<double> v) { return !std::isfinite(v.real()) || !std::isfinite(v.imag()); }

// 1/(sigma_max/sigma_min) for a 2x2 real matrix, matching MATLAB's rcond
// well enough to detect near-singular (collinear) systems.
double rcond2x2(const Eigen::Matrix2d& A) {
    const Eigen::JacobiSVD<Eigen::Matrix2d> svd(A);
    const double smax = svd.singularValues()(0);
    const double smin = svd.singularValues()(1);
    if (smax == 0.0) return 0.0;
    return smin / smax;
}
}  // namespace

FindZ0Result findz0(const std::string& prefix, const Eigen::VectorXcd& wpIn,
                     const std::function<Eigen::VectorXcd(const Eigen::VectorXcd&)>& mapfun,
                     const Eigen::VectorXcd& wIn, const Eigen::VectorXd& betaIn, const Eigen::VectorXcd& zIn,
                     std::complex<double> /*c*/, const Eigen::MatrixXd& qdatIn) {
    const int n = static_cast<int>(wIn.size());
    const bool fromDisk = !prefix.empty() && prefix[0] == 'd';
    const bool fromHp = prefix == "hp";
    const bool fromStrip = prefix == "st";
    const bool fromRect = prefix == "r";

    Eigen::VectorXcd w = wIn;
    Eigen::VectorXcd z = zIn;
    Eigen::VectorXd beta = betaIn;
    Eigen::MatrixXd qdat = qdatIn;

    int kinf0 = -1;  // 0-indexed, only set for from_strip
    if (fromStrip) {
        int atinf0 = -1;
        for (int i = 0; i < n; ++i)
            if (isInf(z(i))) {
                atinf0 = i;
                break;
            }
        std::vector<int> renum;
        for (int i = atinf0; i < n; ++i) renum.push_back(i);
        for (int i = 0; i < atinf0; ++i) renum.push_back(i);

        Eigen::VectorXcd w2(n), z2(n);
        Eigen::VectorXd beta2(n);
        for (int i = 0; i < n; ++i) {
            w2(i) = w(renum[i]);
            z2(i) = z(renum[i]);
            beta2(i) = beta(renum[i]);
        }
        w = w2;
        z = z2;
        beta = beta2;

        Eigen::MatrixXd qdat2 = qdat;
        for (int i = 0; i < n; ++i) {
            qdat2.col(i) = qdat.col(renum[i]);
            qdat2.col(n + 1 + i) = qdat.col(n + 1 + renum[i]);
        }
        qdat = qdat2;

        for (int i = n - 1; i >= 0; --i)
            if (isInf(z(i))) {
                kinf0 = i;
                break;
            }
    }

    // argw(j) (0-indexed), j=0..n-1
    Eigen::VectorXd argw(n);
    if (fromStrip) {
        argw(0) = std::arg(w(2) - w(1));
        for (int j = 1; j < n; ++j) {
            const int srcBeta = (j + 2) % n;  // beta([3:n,1]) 1-indexed -> 0-indexed (j+2 mod n)
            argw(j) = argw(j - 1) - M_PI * beta(srcBeta);
        }
        // rotate: argw = argw([n,1:n-1]) (1-indexed) -> 0-indexed: [n-1,0,1,...,n-2]
        Eigen::VectorXd rotated(n);
        rotated(0) = argw(n - 1);
        for (int j = 1; j < n; ++j) rotated(j) = argw(j - 1);
        argw = rotated;
    } else {
        argw(0) = std::arg(w(1) - w(0));
        for (int j = 1; j < n; ++j) argw(j) = argw(j - 1) - M_PI * beta(j);
    }

    std::vector<bool> infty(n);
    for (int i = 0; i < n; ++i) infty[i] = isInf(w(i));
    std::vector<int> fwd(n);
    for (int i = 0; i < n; ++i) fwd[i] = (i + 1) % n;

    Eigen::VectorXcd anchor(n);
    for (int i = 0; i < n; ++i) anchor(i) = infty[i] ? w(fwd[i]) : w(i);

    Eigen::VectorXcd direcn(n);
    for (int i = 0; i < n; ++i) {
        direcn(i) = std::exp(std::complex<double>(0.0, argw(i)));
        if (infty[i]) direcn(i) = -direcn(i);
    }

    Eigen::VectorXd len(n);
    for (int i = 0; i < n; ++i) len(i) = std::abs(w(fwd[i]) - w(i));

    Eigen::VectorXd argz;
    if (fromDisk) {
        argz.resize(n);
        for (int i = 0; i < n; ++i) {
            argz(i) = std::arg(z(i));
            if (argz(i) <= 0.0) argz(i) += 2.0 * M_PI;
        }
    }

    double factor = 0.5;
    const int m0 = static_cast<int>(wpIn.size());
    std::vector<bool> done(m0, false);
    const double tol = 1000.0 * std::pow(10.0, -static_cast<double>(qdat.rows()));

    Eigen::VectorXcd zbase(n), wbase(n);
    Eigen::VectorXcd z0 = wpIn, w0 = wpIn;
    std::vector<int> idx(m0, -1);  // 0-indexed; -1 sentinel for "not yet assigned"
    bool idxAssigned = false;

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    int iter = 0;
    int mLeft = m0;
    while (mLeft > 0) {
        for (int j = 0; j < n; ++j) {
            if (fromDisk) {
                if (j < n - 1) {
                    zbase(j) = std::exp(std::complex<double>(0.0, factor * argz(j) + (1 - factor) * argz(j + 1)));
                } else {
                    zbase(j) = std::exp(std::complex<double>(0.0, factor * argz(n - 1) + (1 - factor) * (2 * M_PI + argz(0))));
                }
            } else if (fromHp) {
                if (j < n - 2) {
                    zbase(j) = z(j) + factor * (z(j + 1) - z(j));
                } else if (j == n - 2) {
                    zbase(j) = std::max(10.0, z(n - 2).real()) / factor;
                } else {
                    zbase(j) = std::min(-10.0, z(0).real()) / factor;
                }
            } else if (fromStrip) {
                if (j == 0) {
                    zbase(j) = std::min(-1.0, z(1).real()) / factor;
                } else if (j == kinf0 - 1) {
                    zbase(j) = std::max(1.0, z(kinf0 - 1).real()) / factor;
                } else if (j == kinf0) {
                    zbase(j) = std::complex<double>(0.0, 1.0) + std::max(1.0, z(kinf0 + 1).real()) / factor;
                } else if (j == n - 1) {
                    zbase(j) = std::complex<double>(0.0, 1.0) + std::min(-1.0, z(n - 1).real()) / factor;
                } else {
                    zbase(j) = z(j) + factor * (z(j + 1) - z(j));
                }
            } else if (fromRect) {
                zbase(j) = z(j) + factor * (z((j + 1) % n) - z(j));
                if (std::abs(zbase(j)) < 1e-4) {
                    zbase(j) += std::complex<double>(0.0, 0.2);
                } else if (std::abs(zbase(j) - std::complex<double>(0.0, z.array().imag().maxCoeff())) < 1e-4) {
                    zbase(j) -= std::complex<double>(0.0, 0.2);
                }
            }

            Eigen::VectorXcd zbv(1);
            zbv(0) = zbase(j);
            wbase(j) = mapfun(zbv)(0);

            const double proj = ((wbase(j) - anchor(j)) * std::conj(direcn(j))).real();
            wbase(j) = anchor(j) + proj * direcn(j);
        }

        if (!idxAssigned) {
            for (int p = 0; p < m0; ++p) {
                if (done[p]) continue;
                double best = std::numeric_limits<double>::infinity();
                int bestJ = 0;
                for (int j = 0; j < n; ++j) {
                    const double d = std::abs(wpIn(p) - wbase(j));
                    if (d < best) {
                        best = d;
                        bestJ = j;
                    }
                }
                idx[p] = bestJ;
            }
            idxAssigned = true;
        } else {
            for (int p = 0; p < m0; ++p) {
                if (!done[p]) idx[p] = (idx[p] + 1) % n;
            }
        }

        for (int p = 0; p < m0; ++p) {
            if (!done[p]) {
                z0(p) = zbase(idx[p]);
                w0(p) = wbase(idx[p]);
            }
        }

        for (int j = 0; j < n; ++j) {
            std::vector<int> active;
            for (int p = 0; p < m0; ++p)
                if (!done[p] && idx[p] == j) active.push_back(p);
            if (active.empty()) continue;

            for (int p : active) done[p] = true;

            for (int k = 0; k < n; ++k) {
                if (k == j) continue;
                Eigen::Matrix2d A;
                A(0, 0) = direcn(k).real();
                A(1, 0) = direcn(k).imag();
                for (int p : active) {
                    const std::complex<double> dif = w0(p) - wpIn(p);
                    A(0, 1) = dif.real();
                    A(1, 1) = dif.imag();

                    if (rcond2x2(A) < std::numeric_limits<double>::epsilon()) {
                        const double wpx = ((wpIn(p) - anchor(k)) / direcn(k)).real();
                        const double w0x = ((w0(p) - anchor(k)) / direcn(k)).real();
                        if (wpx * w0x < 0.0 || (wpx - len(k)) * (w0x - len(k)) < 0.0) done[p] = false;
                    } else {
                        const std::complex<double> dif2 = w0(p) - anchor(k);
                        Eigen::Vector2d rhs(dif2.real(), dif2.imag());
                        const Eigen::Vector2d s = A.fullPivLu().solve(rhs);
                        if (s(0) >= 0.0 && s(0) <= len(k)) {
                            if (std::abs(s(1) - 1.0) < tol) {
                                z0(p) = zbase(k);
                                w0(p) = wbase(k);
                            } else if (std::abs(s(1)) < tol) {
                                const std::complex<double> n1 = std::conj(wpIn(p) - w0(p)) * std::complex<double>(0.0, 1.0) * direcn(k);
                                if (n1.real() > 0.0) done[p] = false;
                            } else if (s(1) > 0.0 && s(1) < 1.0) {
                                done[p] = false;
                            }
                        }
                    }
                }
            }

            mLeft = 0;
            for (int p = 0; p < m0; ++p)
                if (!done[p]) ++mLeft;
            if (mLeft == 0) break;
        }

        if (iter > 2 * n) {
            throw std::runtime_error("findz0: can't seem to choose starting points");
        }
        ++iter;
        factor = uni(rng);
    }

    return FindZ0Result{z0, w0};
}

}  // namespace sctoolbox
