#pragma once
#include <Eigen/Dense>
#include <functional>
#include <string>

namespace sctoolbox {

struct FindZ0Result {
    Eigen::VectorXcd z0;
    Eigen::VectorXcd w0;
};

// Port of +sctool/findz0.m. prefix is one of "d" (disk), "de" (exterior,
// also disk-based), "hp" (half-plane), "st" (strip), "r" (rectangle).
// MATLAB's `from_rect` branch repurposes the `qdat` argument slot to carry
// the strip-conformal-modulus `L` (passing the real qdat as an extra `aux`
// argument) purely so its positional `feval(mapfun, ...)` calling convention
// works; `L` itself is never read again inside findz0.m. In C++ the mapfun
// closure already captures `L` directly, so no equivalent extra parameter
// is needed here -- just pass the real qdat as `qdat` as usual.
FindZ0Result findz0(const std::string& prefix, const Eigen::VectorXcd& wp,
                     const std::function<Eigen::VectorXcd(const Eigen::VectorXcd&)>& mapfun,
                     const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, const Eigen::VectorXcd& z,
                     std::complex<double> c, const Eigen::MatrixXd& qdat);

}  // namespace sctoolbox
