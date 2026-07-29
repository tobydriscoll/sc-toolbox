#pragma once
#include <Eigen/Dense>
#include <functional>

namespace sctoolbox {

using Fvec = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

// Finite-difference Jacobian approximation (Dennis & Schnabel Algorithm
// A5.4.1). details(12) (0-based index 12, MATLAB's details(13)) supplies
// sqrt(eta) for the step size; a caller passing an unprocessed details
// vector with that entry left at 0 will get a zero step size and hence a
// NaN-filled Jacobian, exactly as the MATLAB version does.
Eigen::MatrixXd nefdjac(const Fvec& fvec, const Eigen::VectorXd& fc, Eigen::VectorXd xc,
                         const Eigen::VectorXd& sx, const Eigen::VectorXd& details, int& nofun);

}  // namespace sctoolbox
