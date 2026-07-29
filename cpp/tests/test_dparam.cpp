#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/dparam.hpp"

TEST_CASE("dparam matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/diskmap_private.gold");
    auto& cases = groups.at("cases_dparam");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& win = c.in("w");
        const auto& betain = c.in("beta");

        const int n = static_cast<int>(win.rows());
        Eigen::VectorXcd w(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            w(i) = win(i, 0);
            beta(i) = betain(i, 0).real();
        }

        const auto result = sctoolbox::dparam(w, beta, 1e-12);

        const auto& expectedZ = c.out("z");
        for (int i = 0; i < n; ++i) {
            CHECK(std::abs(result.z(i) - expectedZ(i, 0)) < c.tol);
        }
        CHECK(std::abs(result.c - c.out("c")(0, 0)) < c.tol);
    }
}
