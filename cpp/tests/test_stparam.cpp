#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/stparam.hpp"

TEST_CASE("stparam matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/stripmap_private.gold");
    auto& cases = groups.at("cases_stparam");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& win = c.in("w");
        const auto& betain = c.in("beta");
        const auto& endsin = c.in("ends");

        const int n = static_cast<int>(win.rows());
        Eigen::VectorXcd w(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            w(i) = win(i, 0);
            beta(i) = betain(i, 0).real();
        }
        std::array<int, 2> ends{static_cast<int>(endsin(0, 0).real()), static_cast<int>(endsin(0, 1).real())};

        const auto result = sctoolbox::stparam(w, beta, ends, 1e-12);

        const auto& expectedZ = c.out("z");
        for (int i = 0; i < n; ++i) {
            const auto e = expectedZ(i, 0);
            if (std::isinf(e.real())) {
                CHECK(std::signbit(e.real()) == std::signbit(result.z(i).real()));
                CHECK(std::isinf(result.z(i).real()));
            } else {
                CHECK(std::abs(result.z(i) - e) < c.tol);
            }
        }
        CHECK(std::abs(result.c - c.out("c")(0, 0)) < c.tol);
    }
}
