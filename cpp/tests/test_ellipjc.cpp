#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/ellipjc.hpp"

TEST_CASE("ellipjc matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/rectmap_private.gold");
    auto& cases = groups.at("cases_ellipjc");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const std::complex<double> u = c.in("u")(0, 0);
        const double L = c.in("L")(0, 0).real();

        Eigen::VectorXcd uv(1);
        uv(0) = u;
        const auto result = sctoolbox::ellipjc(uv, L);

        CHECK(std::abs(result.sn(0) - c.out("sn")(0, 0)) < c.tol);
        CHECK(std::abs(result.cn(0) - c.out("cn")(0, 0)) < c.tol);
        CHECK(std::abs(result.dn(0) - c.out("dn")(0, 0)) < c.tol);
    }
}
