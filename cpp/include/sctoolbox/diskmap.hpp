#pragma once
#include <Eigen/Dense>

#include "sctoolbox/dinvmap.hpp"
#include "sctoolbox/polygon.hpp"

namespace sctoolbox {

// Port of @diskmap (the disk -> polygon-interior Schwarz-Christoffel map
// class), restricted to the "construct from polygon, solve the parameter
// problem with an automatic initial guess" form -- the continuation-map,
// given-prevertex, and tol-override call forms of the MATLAB class are not
// ported, matching every XXparam restriction made earlier in this port.
class DiskMap {
public:
    explicit DiskMap(Polygon poly, double tol = 1e-8);

    // Port of @diskmap/eval.m. Points with abs(zp) > 1+eps map to NaN.
    Eigen::VectorXcd eval(const Eigen::VectorXcd& zp) const;
    // Port of @diskmap/evalinv.m.
    InvMapResult evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0 = {}) const;
    // Port of @diskmap/evaldiff.m.
    Eigen::VectorXcd evaldiff(const Eigen::VectorXcd& zp) const;
    // Port of @diskmap/accuracy.m. Memoized after first call.
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
