#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/rparam.hpp"

TEST_CASE("rparam matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/rectmap_private.gold");
    auto& cases = groups.at("cases_rparam");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& win = c.in("w");
        const auto& betain = c.in("beta");
        const auto& cornersin = c.in("corners");

        const int n = static_cast<int>(win.rows());
        Eigen::VectorXcd w(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            w(i) = win(i, 0);
            beta(i) = betain(i, 0).real();
        }
        std::array<int, 4> cnr;
        for (int i = 0; i < 4; ++i) cnr[i] = static_cast<int>(cornersin(0, i).real());

        const auto result = sctoolbox::rparam(w, beta, cnr, 1e-9);

        // rparam runs two sequential iterative solves (nesolve, then a
        // Newton refinement onto the rectangle boundary); cross-
        // implementation floating-point path differences between this and
        // MATLAB's nesolve compound slightly across both stages, so allow
        // a small multiple of the case's own tolerance (matching the
        // pattern already used for the most iteration-heavy invmap-style
        // comparisons elsewhere in this suite).
        const double tol = 5.0 * c.tol;
        const auto& expectedZ = c.out("z");
        for (int i = 0; i < n; ++i) {
            CHECK(std::abs(result.z(i) - expectedZ(i, 0)) < tol);
        }
        CHECK(std::abs(result.c - c.out("c")(0, 0)) < tol);
        CHECK(std::abs(result.L - c.out("L")(0, 0).real()) < tol);
    }
}
