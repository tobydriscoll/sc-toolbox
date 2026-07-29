#include <catch2/catch_test_macros.hpp>
#include <array>
#include <cmath>
#include <stdexcept>

#include "golden_reader.hpp"
#include "sctoolbox/composite.hpp"
#include "sctoolbox/diskmap.hpp"
#include "sctoolbox/moebius.hpp"

namespace {

Eigen::VectorXcd column(const golden::Matrix& m) {
    Eigen::VectorXcd v(m.rows());
    for (int i = 0; i < m.rows(); ++i) v(i) = m(i, 0);
    return v;
}

std::array<std::complex<double>, 4> coeff(const golden::Matrix& m) {
    REQUIRE(m.rows() == 4);
    return {m(0, 0), m(1, 0), m(2, 0), m(3, 0)};
}

sctoolbox::DiskMap buildMap(const golden::Case& c) {
    const auto& win = c.in("w");
    const auto& betain = c.in("beta");
    const int n = static_cast<int>(win.rows());
    Eigen::VectorXcd w(n);
    Eigen::VectorXd beta(n);
    for (int i = 0; i < n; ++i) {
        w(i) = win(i, 0);
        beta(i) = betain(i, 0).real();
    }
    return sctoolbox::DiskMap(sctoolbox::Polygon(w, beta.array() + 1.0), 1e-12);
}

// Rebuild the composition the generator recorded. `order` names which
// members are present and in what sequence; see gen_composite.m.
sctoolbox::Composite buildComposite(const golden::Case& c, const std::string& order) {
    const sctoolbox::DiskMap map = buildMap(c);
    const sctoolbox::Moebius pre(coeff(c.in("coeff_pre")));
    const sctoolbox::Moebius post(coeff(c.in("coeff_post")));

    sctoolbox::Composite f;
    if (order == "map_then_mob") {
        f.append(sctoolbox::Composite::member(map)).append(sctoolbox::Composite::member(post));
    } else if (order == "mob_then_map") {
        f.append(sctoolbox::Composite::member(pre)).append(sctoolbox::Composite::member(map));
    } else if (order == "mob_map_mob") {
        f.append(sctoolbox::Composite::member(pre))
            .append(sctoolbox::Composite::member(map))
            .append(sctoolbox::Composite::member(post));
    } else {
        FAIL("unknown composite order: " << order);
    }
    return f;
}

golden::GroupMap load() {
    return golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/composite.gold");
}

}  // namespace

TEST_CASE("Composite::eval matches MATLAB goldens") {
    auto groups = load();
    auto& cases = groups.at("cases_eval");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        INFO(c.desc);
        const auto f = buildComposite(c, c.inStr("order"));
        const Eigen::VectorXcd wp = f.eval(column(c.in("zp")));
        const auto& expected = c.out("wp");
        REQUIRE(wp.size() == expected.rows());
        for (int i = 0; i < expected.rows(); ++i) CHECK(std::abs(wp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("Composite::inverse matches MATLAB goldens") {
    auto groups = load();
    auto& cases = groups.at("cases_inv");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        INFO(c.desc);
        const auto fi = buildComposite(c, c.inStr("order")).inverse();
        const Eigen::VectorXcd zp = fi.eval(column(c.in("wp")));
        const auto& expected = c.out("zp");
        REQUIRE(zp.size() == expected.rows());
        for (int i = 0; i < expected.rows(); ++i) CHECK(std::abs(zp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("Composite flattens a nested composite into its members") {
    auto groups = load();
    auto& cases = groups.at("cases_flatten");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        INFO(c.desc);
        const sctoolbox::DiskMap map = buildMap(c);
        const sctoolbox::Moebius pre(coeff(c.in("coeff_pre")));
        const sctoolbox::Moebius post(coeff(c.in("coeff_post")));

        sctoolbox::Composite inner;
        inner.append(sctoolbox::Composite::member(pre)).append(sctoolbox::Composite::member(map));
        sctoolbox::Composite outer;
        outer.append(inner).append(sctoolbox::Composite::member(post));

        CHECK(outer.length() == static_cast<int>(std::lround(c.out("nmaps")(0, 0).real())));

        const Eigen::VectorXcd wp = outer.eval(column(c.in("zp")));
        const auto& expected = c.out("wp");
        REQUIRE(wp.size() == expected.rows());
        for (int i = 0; i < expected.rows(); ++i) CHECK(std::abs(wp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("Composite::inverse rejects a member with no inverse") {
    // Stands in for MATLAB's "Can't invert INLINE maps." error.
    sctoolbox::Composite f;
    f.append(sctoolbox::Composite::function(
        [](const Eigen::VectorXcd& z) { return Eigen::VectorXcd(z.array().square()); }, "square"));
    CHECK_THROWS_AS(f.inverse(), std::runtime_error);
}
