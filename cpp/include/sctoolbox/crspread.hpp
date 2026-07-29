#pragma once
#include <Eigen/Dense>

#include "sctoolbox/qgraph.hpp"

namespace sctoolbox {

struct CrSpreadResult {
    Eigen::MatrixXcd ul;  // n3 x m
    Eigen::MatrixXcd dl;  // n3 x m
};

// Port of @crdiskmap/private/crspread.m. quadnum is 0-indexed.
CrSpreadResult crspread(const Eigen::VectorXcd& u, int quadnum, const Eigen::VectorXd& cr, const QGraph& Q);

}  // namespace sctoolbox
