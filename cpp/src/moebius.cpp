#include "sctoolbox/moebius.hpp"

#include <cmath>
#include <limits>

namespace sctoolbox {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

bool isInf(std::complex<double> v) { return std::isinf(v.real()) || std::isinf(v.imag()); }
bool isNan(std::complex<double> v) { return std::isnan(v.real()) || std::isnan(v.imag()); }

// Index of the first infinite entry, or -1. Mirrors MATLAB's find(isinf(v))
// followed by an implicit "use the first" (the colon expression built from it
// only ever consumes one element).
int firstInf(const Eigen::Vector3cd& v) {
    for (int i = 0; i < 3; ++i)
        if (isInf(v(i))) return i;
    return -1;
}

// MATLAB: renum = rem((j-1:j+1)+2,3)+1, with j 1-indexed. Returned 0-indexed:
// the permutation that moves entry j into the middle slot.
std::array<int, 3> renumber(int j0) {
    const int j = j0 + 1;
    std::array<int, 3> p{};
    for (int k = 0; k < 3; ++k) p[k] = (j - 1 + k + 2) % 3;
    return p;
}

Eigen::Vector3cd permute(const Eigen::Vector3cd& v, const std::array<int, 3>& p) {
    return Eigen::Vector3cd(v(p[0]), v(p[1]), v(p[2]));
}

}  // namespace

Moebius::Moebius() : coeff_{0.0, 0.0, 0.0, 0.0}, source_(Eigen::Vector3cd::Zero()), image_(Eigen::Vector3cd::Zero()) {}

Moebius::Moebius(const std::array<std::complex<double>, 4>& coeff)
    : coeff_(coeff), source_(Eigen::Vector3cd::Zero()), image_(Eigen::Vector3cd::Zero()) {}

Moebius::Moebius(std::complex<double> c1, std::complex<double> c2, std::complex<double> c3,
                 std::complex<double> c4)
    : Moebius(std::array<std::complex<double>, 4>{c1, c2, c3, c4}) {}

Moebius::Moebius(const Eigen::Vector3cd& zIn, const Eigen::Vector3cd& wIn) {
    Eigen::Vector3cd z = zIn, w = wIn;
    std::complex<double> t1 = 0.0, t2 = 0.0;
    bool haveA = false;  // MATLAB tests isnan(A(1)) to detect this branch
    std::array<std::complex<double>, 4> A{0.0, 0.0, 0.0, 0.0};

    const int jw = firstInf(w);
    const int jz = firstInf(z);

    if (jw >= 0) {
        // Renumber so that w(2) == Inf (1-indexed), then branch on z.
        const auto p = renumber(jw);
        z = permute(z, p);
        w = permute(w, p);

        const int jz2 = firstInf(z);
        if (jz2 < 0) {
            t1 = z(1) - z(0);
            t2 = z(1) - z(2);
        } else if (jz2 == 1) {
            t1 = 1.0;
            t2 = 1.0;
        } else {
            if (jz2 != 0) {
                // Move Inf to the beginning of z.
                z = Eigen::Vector3cd(z(2), z(1), z(0));
                w = Eigen::Vector3cd(w(2), w(1), w(0));
            }
            A = {w(2) * (z(2) - z(1)) - w(0) * z(2), w(0), -z(1), 1.0};
            haveA = true;
        }
    } else if (jz >= 0) {
        // All w finite; renumber so that z(2) == Inf (1-indexed).
        const auto p = renumber(jz);
        z = permute(z, p);
        w = permute(w, p);
        t1 = w(1) - w(2);
        t2 = w(1) - w(0);
    } else {
        t1 = -(z(1) - z(0)) * (w(2) - w(1));
        t2 = -(z(2) - z(1)) * (w(1) - w(0));
    }

    if (!haveA) {
        A = {w(2) * z(0) * t2 - w(0) * z(2) * t1, w(0) * t1 - w(2) * t2, z(0) * t2 - z(2) * t1,
             t1 - t2};
    }

    coeff_ = A;
    source_ = z;
    image_ = w;
    hasPoints_ = true;
}

std::complex<double> Moebius::eval(std::complex<double> z) const {
    std::complex<double> f;
    if (isInf(z)) {
        // MATLAB computes coeff(2)/coeff(4) and lets the trailing
        // isnan -> Inf sweep absorb a zero denominator.
        f = (coeff_[3] == std::complex<double>(0.0, 0.0)) ? std::complex<double>(kInf, 0.0)
                                                          : coeff_[1] / coeff_[3];
    } else {
        std::complex<double> num = coeff_[1] * z + coeff_[0];
        std::complex<double> den = coeff_[3] * z + coeff_[2];
        if (std::abs(den) < 3.0 * std::numeric_limits<double>::epsilon()) {
            return std::complex<double>(kInf, 0.0);
        }
        f = num / den;
    }
    if (isNan(f)) f = std::complex<double>(kInf, 0.0);
    return f;
}

Eigen::VectorXcd Moebius::eval(const Eigen::VectorXcd& z) const {
    Eigen::VectorXcd f(z.size());
    for (int i = 0; i < z.size(); ++i) f(i) = eval(z(i));
    return f;
}

std::complex<double> Moebius::diff(std::complex<double> z) const {
    const std::complex<double> a = coeff_[0], b = coeff_[1], c = coeff_[2], d = coeff_[3];
    const std::complex<double> den = d * z + c;
    return (b * c - a * d) / (den * den);
}

Eigen::VectorXcd Moebius::diff(const Eigen::VectorXcd& z) const {
    Eigen::VectorXcd f(z.size());
    for (int i = 0; i < z.size(); ++i) f(i) = diff(z(i));
    return f;
}

Moebius Moebius::inverse() const {
    Moebius m;
    m.coeff_ = {coeff_[0], -coeff_[2], -coeff_[1], coeff_[3]};
    m.source_ = image_;
    m.image_ = source_;
    m.hasPoints_ = hasPoints_;
    return m;
}

Moebius Moebius::normal() const {
    std::complex<double> c = coeff_[2];
    if (c == std::complex<double>(0.0, 0.0)) c = coeff_[3];
    Moebius m = *this;
    for (auto& v : m.coeff_) v /= c;
    return m;
}

Moebius Moebius::compose(const Moebius& inner) const {
    const auto& C = coeff_;
    const auto& D = inner.coeff_;
    return Moebius(std::array<std::complex<double>, 4>{
        C[0] * D[2] + C[1] * D[0], C[0] * D[3] + C[1] * D[1], C[2] * D[2] + C[3] * D[0],
        C[2] * D[3] + C[3] * D[1]});
}

Moebius Moebius::plus(std::complex<double> s) const {
    auto C = coeff_;
    C[0] += s * C[2];
    C[1] += s * C[3];
    return Moebius(C);
}

Moebius Moebius::minus(std::complex<double> s) const { return plus(-s); }

Moebius Moebius::times(std::complex<double> s) const {
    auto C = coeff_;
    C[0] *= s;
    C[1] *= s;
    return Moebius(C);
}

Moebius Moebius::dividedBy(std::complex<double> s) const {
    auto C = coeff_;
    C[2] *= s;
    C[3] *= s;
    return Moebius(C);
}

Moebius Moebius::divide(std::complex<double> s, const Moebius& m) {
    // MATLAB exchanges numerator and denominator, then scales.
    const auto& C = m.coeff_;
    return Moebius(std::array<std::complex<double>, 4>{C[2], C[3], C[0], C[1]}).times(s);
}

Moebius Moebius::negated() const {
    Moebius m = *this;
    m.coeff_[0] = -m.coeff_[0];
    m.coeff_[1] = -m.coeff_[1];
    return m;
}

}  // namespace sctoolbox
