#pragma once
#include <Eigen/Dense>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

struct CrParamResult {
    Eigen::VectorXcd w;     // possibly-subdivided polygon vertices
    Eigen::VectorXd beta;   // matching turning angles
    Eigen::VectorXd cr;     // n-3 prevertex crossratios
    Eigen::MatrixXcd aff;   // (n-3) x 2 affine table
    QGraph Q;
    std::vector<bool> orig;  // true at entries of w corresponding to original vertices
    Eigen::MatrixXd qdat;
};

// Port of @crdiskmap/private/crparam.m, restricted to the cr0=[] (automatic
// initial guess) calling form and the full (non-abbreviated) output form,
// matching the other XXparam ports' restrictions.
CrParamResult crparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, double tol = 1e-8, int method = 2);

}  // namespace sctoolbox
