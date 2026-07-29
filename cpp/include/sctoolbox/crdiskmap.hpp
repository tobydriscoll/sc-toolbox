#pragma once
#include <Eigen/Dense>

#include "sctoolbox/polygon.hpp"
#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

// Port of @crdiskmap (disk -> polygon-interior Schwarz-Christoffel map
// class, cross-ratio formulation), restricted as DiskMap is.
class CrDiskMap {
public:
    explicit CrDiskMap(Polygon poly, double tol = 1e-8);

    // Port of @crdiskmap/eval.m. Points with abs(zp) > 1+eps map to NaN.
    Eigen::VectorXcd eval(const Eigen::VectorXcd& zp) const;
    // Port of @crdiskmap/evalinv.m. Points outside the polygon map to NaN.
    Eigen::VectorXcd evalinv(const Eigen::VectorXcd& wp) const;
    // Port of @crdiskmap/evaldiff.m.
    Eigen::VectorXcd evaldiff(const Eigen::VectorXcd& zp) const;
    // Port of @crdiskmap/accuracy.m. Memoized after first call.
    double accuracy() const;

    const Polygon& polygon() const { return poly_; }
    const Eigen::VectorXd& crossratio() const { return cr_; }
    const Eigen::MatrixXcd& affine() const { return aff_; }
    const QGraph& qlgraph() const { return Q_; }
    const Eigen::VectorXcd& wcfix() const { return wcfix_; }
    const Eigen::MatrixXd& qdata() const { return qdat_; }

private:
    Polygon poly_;
    Eigen::VectorXd cr_;
    Eigen::MatrixXcd aff_;
    QGraph Q_;
    Eigen::VectorXcd wcfix_;
    Eigen::MatrixXd qdat_;
    double tol_;
    mutable double acc_ = -1.0;
};

}  // namespace sctoolbox
