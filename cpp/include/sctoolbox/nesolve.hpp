#pragma once
#include <Eigen/Dense>

#include "sctoolbox/nefdjac.hpp"

namespace sctoolbox {

struct NesolveResult {
    Eigen::VectorXd xf;
    int termcode;
};

// Port of +sctool/nesolve.m restricted to the (fvec, x0, details) calling
// form -- no fparam / analytic jacobian / scale arguments, since no real
// caller in this codebase uses them and CPP_PLAN.md's golden cases only
// exercise the plain 3-argument form. `details` is the 16-element Dennis &
// Schnabel options vector; entries keep the same meaning as the MATLAB
// version's details(1..16), just reindexed to details[0..15] here.
//
// `identityInitialJacobian` (false by default) selects +sctool/nesolvei.m's
// sole behavioral difference from nesolve.m: the very first Jacobian is the
// identity matrix instead of a finite-difference approximation (subsequent
// restart Jacobians, if any, still use nefdjac, matching nesolvei.m).
// crparam.m calls nesolvei specifically because its residual's Jacobian is
// already well-approximated by the identity near a good starting guess.
NesolveResult nesolve(const Fvec& fvec, const Eigen::VectorXd& x0, Eigen::VectorXd details,
                      bool identityInitialJacobian = false);

}  // namespace sctoolbox
