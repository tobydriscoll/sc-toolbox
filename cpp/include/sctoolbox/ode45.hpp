#pragma once
#include <Eigen/Dense>
#include <functional>

namespace sctoolbox {

using OdeFun = std::function<Eigen::VectorXd(double t, const Eigen::VectorXd& y)>;

// Adaptive Dormand-Prince RK5(4) integrator from t0 to t1. Substitutes for
// MATLAB's ode23/ode113 in the XXinvmap ODE-continuation step: only the
// final state is needed (the result is polished by Newton iteration
// afterward), so bit-exact replication of MATLAB's stepper is unnecessary —
// any integrator accurate to abstol/reltol lands in the same Newton basin.
Eigen::VectorXd ode45(const OdeFun& f, double t0, double t1, const Eigen::VectorXd& y0, double abstol,
                       double reltol);

}  // namespace sctoolbox
