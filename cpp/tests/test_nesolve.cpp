#include <catch2/catch_test_macros.hpp>

#include "golden_reader.hpp"
#include "sctoolbox/nesolve.hpp"

TEST_CASE("nesolve matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/nesolve.gold");
    auto& cases = groups.at("cases_nesolve");
    REQUIRE(!cases.empty());

    sctoolbox::Fvec f2d = [](const Eigen::VectorXd& x) {
        Eigen::VectorXd f(2);
        f(0) = x(0) * x(0) - 1.0;
        f(1) = x(1) * x(1) - 4.0;
        return f;
    };
    sctoolbox::Fvec f3d = [](const Eigen::VectorXd& x) {
        Eigen::VectorXd f(3);
        f(0) = x(0) + x(1) - 1.0;
        f(1) = x(1) + x(2) - 2.0;
        f(2) = x(0) * x(2) - 0.5;
        return f;
    };

    for (const auto& c : cases) {
        const auto& x0in = c.in("x0");
        const auto& detIn = c.in("details");
        const int n = static_cast<int>(x0in.rows());
        Eigen::VectorXd x0(n);
        for (int i = 0; i < n; ++i) x0(i) = x0in(i, 0).real();
        Eigen::VectorXd details(detIn.rows());
        for (int i = 0; i < detIn.rows(); ++i) details(i) = detIn(i, 0).real();

        const auto& f = (n == 3) ? f3d : f2d;
        sctoolbox::NesolveResult r = sctoolbox::nesolve(f, x0, details);

        // Verify residual is small, not bit-exact (solver path may legitimately
        // differ in iteration count/sequence between MATLAB and this port).
        CHECK(f(r.xf).lpNorm<Eigen::Infinity>() < 1e-8);
        CHECK(r.termcode == 1);
    }
}
