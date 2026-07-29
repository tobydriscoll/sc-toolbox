#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/ellipkkp.hpp"

TEST_CASE("ellipkkp matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/rectmap_private.gold");
    auto& cases = groups.at("cases_ellipkkp");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const double L = c.in("L")(0, 0).real();
        const auto result = sctoolbox::ellipkkp(L);

        CHECK(std::abs(result.K - c.out("K")(0, 0).real()) < c.tol);
        CHECK(std::abs(result.Kp - c.out("Kp")(0, 0).real()) < c.tol);
    }
}
