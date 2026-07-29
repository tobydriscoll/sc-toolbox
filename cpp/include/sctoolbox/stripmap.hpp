#pragma once
#include <Eigen/Dense>
#include <array>

#include "sctoolbox/dinvmap.hpp"
#include "sctoolbox/polygon.hpp"

namespace sctoolbox {

// Port of @stripmap (infinite strip -> polygon-interior Schwarz-Christoffel
// map class), restricted as DiskMap is. `ends` is the 1-indexed 2-vector of
// vertex positions mapping to the strip's ends, as stparam.m requires (no
// interactive mouse-click fallback here).
class StripMap {
public:
    StripMap(Polygon poly, std::array<int, 2> ends, double tol = 1e-8);

    // Port of @stripmap/eval.m. Points outside 0<=Im(zp)<=1 map to NaN.
    Eigen::VectorXcd eval(const Eigen::VectorXcd& zp) const;
    // Port of @stripmap/evalinv.m.
    InvMapResult evalinv(const Eigen::VectorXcd& wp, const Eigen::VectorXcd& z0 = {}) const;
    // Port of @stripmap/evaldiff.m.
    Eigen::VectorXcd evaldiff(const Eigen::VectorXcd& zp) const;
    // Port of @stripmap/accuracy.m. Memoized after first call.
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
