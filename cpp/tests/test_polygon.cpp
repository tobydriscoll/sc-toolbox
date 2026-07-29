#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/polygon.hpp"

TEST_CASE("polygon constructor matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/polygon.gold");
    auto& cases = groups.at("cases_polygon");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& win = c.in("w");
        const auto& alphain = c.in("alpha");
        const int n = static_cast<int>(win.rows());

        Eigen::VectorXcd w(n);
        for (int i = 0; i < n; ++i) w(i) = win(i, 0);

        sctoolbox::Polygon p = [&]() {
            if (alphain.rows() == 0) return sctoolbox::Polygon(w);
            Eigen::VectorXd alpha(alphain.rows());
            for (int i = 0; i < alphain.rows(); ++i) alpha(i) = alphain(i, 0).real();
            return sctoolbox::Polygon(w, alpha);
        }();

        const auto& vg = c.out("vertex");
        const auto& ag = c.out("angle");
        const auto& ig = c.out("isinf");

        REQUIRE(p.vertex().size() == vg.rows());
        for (int i = 0; i < vg.rows(); ++i) {
            const auto expected = vg(i, 0);
            if (std::isinf(expected.real()) || std::isinf(expected.imag())) {
                CHECK(p.vertex()(i) == expected);
            } else {
                CHECK(std::abs(p.vertex()(i) - expected) < c.tol);
            }
        }
        REQUIRE(p.angle().size() == ag.rows());
        for (int i = 0; i < ag.rows(); ++i) CHECK(std::abs(p.angle()(i) - ag(i, 0).real()) < c.tol);
        CHECK(p.isInf() == (ig(0, 0).real() != 0.0));
    }
}
