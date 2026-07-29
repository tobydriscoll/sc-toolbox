#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Port of @polygon/polygon.m + @polygon/angle.m. Vertices are stored
// counterclockwise (reversing the input order if necessary, matching the
// MATLAB constructor); `angle` holds interior angles normalized by pi.
class Polygon {
public:
    // Bounded polygon: angles computed automatically from vertex geometry.
    explicit Polygon(Eigen::VectorXcd vertices);
    // Angles supplied explicitly (required for unbounded polygons; also
    // accepted for bounded ones, mirroring MATLAB's POLYGON(W,ALPHA) form).
    Polygon(Eigen::VectorXcd vertices, Eigen::VectorXd angles);

    const Eigen::VectorXcd& vertex() const { return vertex_; }
    const Eigen::VectorXd& angle() const { return angle_; }
    int length() const { return static_cast<int>(vertex_.size()); }
    bool isInf() const;

private:
    void init(Eigen::VectorXcd w, Eigen::VectorXd alpha);

    Eigen::VectorXcd vertex_;
    Eigen::VectorXd angle_;
};

}  // namespace sctoolbox
