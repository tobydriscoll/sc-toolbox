#pragma once
#include <Eigen/Dense>
#include <complex>

namespace sctoolbox {

// Port of @stripmap/private/stderiv.m. z/beta must include the strip's
// two infinite "end" prevertices (stderiv strips them out internally).
// j follows MATLAB's 1-indexed convention directly: j == 0 means no
// Gauss-Jacobi normalization is applied, otherwise j is the 1-indexed
// position (within the full z, including the infinite ends) of the
// singularity to normalize.
Eigen::VectorXcd stderiv(const Eigen::VectorXcd& zp, const Eigen::VectorXcd& z, const Eigen::VectorXd& beta,
                         std::complex<double> c = 1.0, int j = 0);

}  // namespace sctoolbox
