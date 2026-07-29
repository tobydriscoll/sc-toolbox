#include "sctoolbox/scinvopt.hpp"

namespace sctoolbox {

InvOpt scinvopt(const std::vector<double>& options) {
    double opt[3] = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < options.size() && i < 3; ++i) opt[i] = options[i];
    if (opt[0] == 0.0) opt[0] = 0.0;
    if (opt[1] == 0.0) opt[1] = 1e-8;
    if (opt[2] == 0.0) opt[2] = 10.0;

    InvOpt out;
    out.ode = (opt[0] == 0.0) || (opt[0] == 1.0);
    out.newton = (opt[0] == 0.0) || (opt[0] == 2.0);
    out.tol = opt[1];
    out.maxiter = static_cast<int>(opt[2]);
    return out;
}

}  // namespace sctoolbox
