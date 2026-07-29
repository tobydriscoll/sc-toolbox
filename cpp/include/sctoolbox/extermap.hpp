#pragma once
#include <Eigen/Dense>

#include "sctoolbox/dinvmap.hpp"
#include "sctoolbox/polygon.hpp"

namespace sctoolbox {

// Port of @extermap (disk -> polygon-exterior Schwarz-Christoffel map
// class), restricted as DiskMap is. Internally, every method re-derives the
// "negated/flipped" (w,beta) convention deparam.m's residual needs
// (w=flipud(vertex(poly)), beta=flipud(1-angle(poly))) -- matching the
// @extermap/*.m methods exactly, and consistent with how the deparam
// golden-test generator had to be fixed (see CPP_PLAN.md 9.4) to call
// deparam with this same convention rather than the normal alpha-1 one.
class ExterMap {
public:
    explicit ExterMap(Polygon poly, double tol = 1e-8);

    // Port of @extermap/eval.m. Points with abs(zp) > 1+eps map to NaN.
    Eigen::VectorXcd eval(const Eigen::VectorXcd& zp) const;
    // Port of @extermap/evalinv.m.
    InvMapResult evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0 = {}) const;
    // Port of @extermap/evaldiff.m.
    Eigen::VectorXcd evaldiff(const Eigen::VectorXcd& zp) const;
    // Port of @extermap/accuracy.m. Memoized after first call.
    double accuracy() const;

    const Polygon& polygon() const { return poly_; }
    const Eigen::VectorXcd& prevertex() const { return z_; }
    std::complex<double> constant() const { return c_; }
    const Eigen::MatrixXd& qdata() const { return qdat_; }

private:
    Polygon poly_;
    Eigen::VectorXcd z_;
    std::complex<double> c_;
    Eigen::MatrixXd qdat_;
    double tol_;
    mutable double acc_ = -1.0;
};

}  // namespace sctoolbox
