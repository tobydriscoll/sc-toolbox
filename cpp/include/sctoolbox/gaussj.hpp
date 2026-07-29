#pragma once
#include <Eigen/Dense>

namespace sctoolbox {

// Nodes and weights for Gauss-Jacobi integration on [-1, 1] with weight
// (1-x)^alf * (1+x)^bet. Direct port of +sctool/gaussj.m (Lanczos
// recurrence -> symmetric tridiagonal eigenproblem).
void gaussj(int n, double alf, double bet, Eigen::VectorXd& z, Eigen::VectorXd& w);

}  // namespace sctoolbox
