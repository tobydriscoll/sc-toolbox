#pragma once
#include <Eigen/Dense>
#include <complex>

#include "sctoolbox/annulusmap_internal.hpp"
#include "sctoolbox/polygon.hpp"

namespace sctoolbox {

// Port of @annulusmap (Schwarz-Christoffel map to a doubly connected
// polygonal region), restricted to bounded outer polygons -- the MATLAB
// constructor's 'truncate' option for unbounded outer polygons is not
// ported, matching this codebase's scope.
class AnnulusMap {
public:
    AnnulusMap(const Polygon& outerPolygon, const Polygon& innerPolygon);

    // Forward map: canonical annulus point -> doubly connected polygon point.
    std::complex<double> eval(std::complex<double> w) const;
    // Inverse map: doubly connected polygon point -> canonical annulus point.
    std::complex<double> evalinv(std::complex<double> z) const;

    int M() const { return data_.M; }
    int N() const { return data_.N; }
    double u() const { return params_.u; }
    std::complex<double> c() const { return params_.c; }
    const Eigen::VectorXcd& w0() const { return params_.w0; }
    const Eigen::VectorXcd& w1() const { return params_.w1; }
    const Eigen::VectorXd& phi0() const { return params_.phi0; }
    const Eigen::VectorXd& phi1() const { return params_.phi1; }
    const Eigen::VectorXd& qwork() const { return qwork_; }

private:
    AnnulusData data_;
    Eigen::VectorXd qwork_;
    DscParams params_;
    static constexpr int kNptq = 8;
};

}  // namespace sctoolbox
