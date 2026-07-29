#include <catch2/catch_test_macros.hpp>
#include <array>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/moebius.hpp"

namespace {

Eigen::VectorXcd column(const golden::Matrix& m) {
    Eigen::VectorXcd v(m.rows());
    for (int i = 0; i < m.rows(); ++i) v(i) = m(i, 0);
    return v;
}

Eigen::Vector3cd triple(const golden::Matrix& m) {
    REQUIRE(m.rows() == 3);
    return Eigen::Vector3cd(m(0, 0), m(1, 0), m(2, 0));
}

std::array<std::complex<double>, 4> coeff(const golden::Matrix& m) {
    REQUIRE(m.rows() == 4);
    return {m(0, 0), m(1, 0), m(2, 0), m(3, 0)};
}

bool isInf(std::complex<double> v) { return std::isinf(v.real()) || std::isinf(v.imag()); }
bool isNan(std::complex<double> v) { return std::isnan(v.real()) || std::isnan(v.imag()); }

// How to compare a non-finite entry.
enum class NonFinite {
    // Inf must match Inf. eval.m funnels every degenerate result to Inf and
    // never returns NaN, so this is the strict mode used for eval outputs;
    // only the fact of being infinite is compared, since the sign and
    // imaginary part of an infinite result carry no meaning.
    kInfExact,
    // Any non-finite matches any other. Needed for diff.m, which -- unlike
    // eval.m -- has no infinity handling at all: at z = Inf it evaluates
    // 0*Inf and returns NaN. C++ agrees on this platform, but the standard
    // explicitly permits a complex multiplication to recover an infinity
    // where MATLAB's IEEE arithmetic produces NaN, so pinning down which
    // flavor of non-finite comes out would be testing the standard library
    // rather than this port.
    kAnyNonFinite,
};

void checkClose(std::complex<double> got, std::complex<double> want, double tol,
                NonFinite mode = NonFinite::kInfExact) {
    if (mode == NonFinite::kAnyNonFinite) {
        const bool wantBad = isInf(want) || isNan(want);
        const bool gotBad = isInf(got) || isNan(got);
        if (wantBad || gotBad) {
            CHECK(wantBad == gotBad);
            return;
        }
    } else if (isInf(want) || isInf(got)) {
        CHECK(isInf(want) == isInf(got));
        return;
    }
    CHECK(std::abs(got - want) < tol);
}

void checkVector(const Eigen::VectorXcd& got, const golden::Matrix& want, double tol,
                 NonFinite mode = NonFinite::kInfExact) {
    REQUIRE(got.size() == want.rows());
    for (int i = 0; i < want.rows(); ++i) checkClose(got(i), want(i, 0), tol, mode);
}

void checkCoeff(const std::array<std::complex<double>, 4>& got, const golden::Matrix& want,
                double tol) {
    REQUIRE(want.rows() == 4);
    for (int i = 0; i < 4; ++i) checkClose(got[i], want(i, 0), tol);
}

golden::GroupMap load() {
    return golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/moebius.gold");
}

}  // namespace

TEST_CASE("Moebius three-point constructor matches MATLAB goldens") {
    auto groups = load();
    auto& cases = groups.at("cases_construct");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        INFO(c.desc);
        const sctoolbox::Moebius m(triple(c.in("z")), triple(c.in("w")));
        checkCoeff(m.coeff(), c.out("coeff"), c.tol);
        // source/image are stored renumbered by the constructor.
        checkVector(m.source(), c.out("source"), c.tol);
        checkVector(m.image(), c.out("image"), c.tol);
    }
}

TEST_CASE("Moebius eval and diff match MATLAB goldens") {
    auto groups = load();
    auto& cases = groups.at("cases_eval");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        INFO(c.desc);
        const sctoolbox::Moebius m(coeff(c.in("coeff")));
        const Eigen::VectorXcd zp = column(c.in("zp"));
        checkVector(m.eval(zp), c.out("fp"), c.tol);
        checkVector(m.diff(zp), c.out("dp"), c.tol, NonFinite::kAnyNonFinite);
    }
}

TEST_CASE("Moebius operators match MATLAB goldens") {
    auto groups = load();
    auto& cases = groups.at("cases_ops");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        INFO(c.desc);
        const sctoolbox::Moebius m(coeff(c.in("coeff")));
        const std::complex<double> s = c.in("s")(0, 0);
        checkCoeff(m.inverse().coeff(), c.out("coeff_inv"), c.tol);
        checkCoeff(m.normal().coeff(), c.out("coeff_normal"), c.tol);
        checkCoeff((-m).coeff(), c.out("coeff_neg"), c.tol);
        checkCoeff((m + s).coeff(), c.out("coeff_plus"), c.tol);
        checkCoeff((m - s).coeff(), c.out("coeff_minus"), c.tol);
        checkCoeff((m * s).coeff(), c.out("coeff_times"), c.tol);
        checkCoeff((m / s).coeff(), c.out("coeff_rdiv"), c.tol);
        checkCoeff((s / m).coeff(), c.out("coeff_ldiv"), c.tol);
    }
}

TEST_CASE("Moebius inverse swaps source and image") {
    auto groups = load();
    auto& cases = groups.at("cases_invpoints");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        INFO(c.desc);
        const sctoolbox::Moebius mi =
            sctoolbox::Moebius(triple(c.in("z")), triple(c.in("w"))).inverse();
        checkCoeff(mi.coeff(), c.out("coeff_inv"), c.tol);
        checkVector(mi.source(), c.out("source_inv"), c.tol);
        checkVector(mi.image(), c.out("image_inv"), c.tol);
    }
}

TEST_CASE("Moebius composition matches MATLAB goldens") {
    auto groups = load();
    auto& cases = groups.at("cases_compose");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        INFO(c.desc);
        const sctoolbox::Moebius m1(coeff(c.in("coeff1")));
        const sctoolbox::Moebius m2(coeff(c.in("coeff2")));
        const sctoolbox::Moebius m = m1.compose(m2);
        checkCoeff(m.coeff(), c.out("coeff"), c.tol);
        checkVector(m.eval(column(c.in("zp"))), c.out("fp"), c.tol);
    }
}

TEST_CASE("Moebius three-point map carries its source points to its image points") {
    // A property check independent of the goldens: the defining property of
    // the (z,w) constructor, verified for the all-finite cases (the ones
    // where every source point has a well-defined finite comparison).
    const Eigen::Vector3cd z(std::complex<double>(0, 0), std::complex<double>(1, 0),
                             std::complex<double>(0, 1));
    const Eigen::Vector3cd w(std::complex<double>(2, 1), std::complex<double>(3, 0),
                             std::complex<double>(0, -1));
    const sctoolbox::Moebius m(z, w);
    for (int i = 0; i < 3; ++i) CHECK(std::abs(m.eval(z(i)) - w(i)) < 1e-12);

    // ...and the inverse carries them back.
    const sctoolbox::Moebius mi = m.inverse();
    for (int i = 0; i < 3; ++i) CHECK(std::abs(mi.eval(w(i)) - z(i)) < 1e-12);
}
