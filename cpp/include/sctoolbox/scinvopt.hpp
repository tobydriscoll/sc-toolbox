#pragma once
#include <vector>

namespace sctoolbox {

struct InvOpt {
    bool ode = true;
    bool newton = true;
    double tol = 1e-8;
    int maxiter = 10;
};

// Port of +sctool/scinvopt.m.
InvOpt scinvopt(const std::vector<double>& options = {});

}  // namespace sctoolbox
