#pragma once
#include <Eigen/Dense>
#include <limits>

namespace sctoolbox {

// Winding number (collapsed to 0/1, matching MATLAB's `index = logical(index)`
// final line) of each point in z with respect to the polygon w, beta.
// Port of +sctool/isinpoly.m's (z,w,beta,tol) form.
Eigen::VectorXd isinpoly(const Eigen::VectorXcd& z, const Eigen::VectorXcd& w, const Eigen::VectorXd& beta,
                          double tol = std::numeric_limits<double>::epsilon());

// Convenience overload matching the 2-argument MATLAB call site used by
// polygon angle computation: beta defaults to scangle(w).
Eigen::VectorXd isinpoly(const Eigen::VectorXcd& z, const Eigen::VectorXcd& w);

// Port of isinpoly.m's (z,w,tol) 3-argument form (MATLAB reinterprets a
// scalar 3rd argument as tol, not beta): beta defaults to scangle(w).
Eigen::VectorXd isinpoly(const Eigen::VectorXcd& z, const Eigen::VectorXcd& w, double tol);

}  // namespace sctoolbox
