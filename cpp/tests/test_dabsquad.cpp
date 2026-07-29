#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/dabsquad.hpp"

TEST_CASE("dabsquad matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/diskmap_private.gold");
    auto& cases = groups.at("cases_dabsquad");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& z1in = c.in("z1");
        const auto& z2in = c.in("z2");
        const auto& sing1in = c.in("sing1");
        const auto& zin = c.in("z");
        const auto& betain = c.in("beta");
        const auto& qdat = c.in("qdat");

        const int m = static_cast<int>(z1in.rows());
        const int n = static_cast<int>(zin.rows());

        Eigen::VectorXcd z1(m), z2(m);
        std::vector<int> sing1(m);
        for (int i = 0; i < m; ++i) {
            z1(i) = z1in(i, 0);
            z2(i) = z2in(i, 0);
            sing1[i] = static_cast<int>(sing1in(i, 0).real());
        }

        Eigen::VectorXcd z(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            z(i) = zin(i, 0);
            beta(i) = betain(i, 0).real();
        }

        const auto result = sctoolbox::dabsquad(z1, z2, sing1, z, beta, qdat.real());

        const auto& expected = c.out("I");
        for (int i = 0; i < m; ++i) {
            CHECK(std::abs(result(i) - expected(i, 0).real()) < c.tol);
        }
    }
}
