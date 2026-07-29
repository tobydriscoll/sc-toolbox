#pragma once
#include <Eigen/Dense>

#include "sctoolbox/dinvmap.hpp"
#include "sctoolbox/polygon.hpp"

namespace sctoolbox {

// Port of @hplmap (half-plane -> polygon-interior Schwarz-Christoffel map
// class), restricted as DiskMap is.
class HplMap {
public:
    explicit HplMap(Polygon poly, double tol = 1e-8);

    // Port of @hplmap/eval.m. Points with imag(zp) <= -eps map to NaN.
    Eigen::VectorXcd eval(const Eigen::VectorXcd& zp) const;
    // Port of @hplmap/evalinv.m.
    InvMapResult evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0 = {}) const;
    // Port of @hplmap/evaldiff.m.
    Eigen::VectorXcd evaldiff(const Eigen::VectorXcd& zp) const;
    // Port of @hplmap/accuracy.m. Memoized after first call.
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
