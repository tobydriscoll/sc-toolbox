#pragma once
#include <Eigen/Dense>
#include <array>

#include "sctoolbox/dinvmap.hpp"
#include "sctoolbox/polygon.hpp"

namespace sctoolbox {

// Port of @rectmap (rectangle -> polygon-interior Schwarz-Christoffel map
// class), restricted as DiskMap is. `corners` is the 1-indexed 4-vector of
// corner positions, as rparam.m requires (first two entries describe a
// long side, CCW order).
class RectMap {
public:
    RectMap(Polygon poly, std::array<int, 4> corners, double tol = 1e-8);

    // Port of @rectmap/eval.m. Points outside the source rectangle map to NaN.
    Eigen::VectorXcd eval(const Eigen::VectorXcd& zp) const;
    // Port of @rectmap/evalinv.m.
    InvMapResult evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0 = {}) const;
    // Port of @rectmap/evaldiff.m.
    Eigen::VectorXcd evaldiff(const Eigen::VectorXcd& zp) const;
    // Port of @rectmap/accuracy.m. Memoized after first call.
    double accuracy() const;
    // Port of @rectmap/corners.m: indices of the rectangle/generalized
    // quadrilateral corners within prevertex().
    std::array<int, 4> corners() const;

    const Polygon& polygon() const { return poly_; }
    const Eigen::VectorXcd& prevertex() const { return z_; }
    std::complex<double> constant() const { return c_; }
    double stripL() const { return L_; }
    const Eigen::MatrixXd& qdata() const { return qdat_; }

private:
    Polygon poly_;
    Eigen::VectorXcd z_;
    std::complex<double> c_;
    double L_;
    Eigen::MatrixXd qdat_;
    double tol_;
    mutable double acc_ = -1.0;
};

}  // namespace sctoolbox
