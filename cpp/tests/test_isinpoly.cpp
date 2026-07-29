#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/isinpoly.hpp"

TEST_CASE("isinpoly matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/isinpoly.gold");
    auto& cases = groups.at("cases_isinpoly");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& zin = c.in("z");
        const auto& win = c.in("w");
        const int nz = static_cast<int>(zin.rows());
        const int nw = static_cast<int>(win.rows());

        Eigen::VectorXcd z(nz), w(nw);
        for (int i = 0; i < nz; ++i) z(i) = zin(i, 0);
        for (int i = 0; i < nw; ++i) w(i) = win(i, 0);

        const Eigen::VectorXd idx = sctoolbox::isinpoly(z, w);

        const auto& idxg = c.out("idx");
        REQUIRE(idx.size() == idxg.rows());
        for (int i = 0; i < idxg.rows(); ++i) CHECK(std::abs(idx(i) - idxg(i, 0).real()) < c.tol);
    }
}
