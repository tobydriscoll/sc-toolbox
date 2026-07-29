#include <catch2/catch_test_macros.hpp>

#include "golden_reader.hpp"
#include "sctoolbox/neqrdcmp.hpp"

TEST_CASE("neqrdcmp matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/nesolve.gold");
    auto& cases = groups.at("cases_neqrdcmp");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& Ain = c.in("A");
        const int n = static_cast<int>(Ain.rows());
        Eigen::MatrixXd M(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) M(i, j) = Ain(i, j).real();

        Eigen::VectorXd M1, M2;
        int sing;
        sctoolbox::neqrdcmp(M, M1, M2, sing);

        const auto& Mg = c.out("M");
        const auto& M1g = c.out("M1");
        const auto& M2g = c.out("M2");
        CHECK(sing == static_cast<int>(c.out("sing")(0, 0).real()));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                CHECK(std::abs(M(i, j) - Mg(i, j).real()) < c.tol);
            }
            CHECK(std::abs(M1(i) - M1g(i, 0).real()) < c.tol);
            CHECK(std::abs(M2(i) - M2g(i, 0).real()) < c.tol);
        }
    }
}
