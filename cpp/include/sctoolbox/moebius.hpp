#pragma once
#include <Eigen/Dense>
#include <array>
#include <complex>

namespace sctoolbox {

// Port of @moebius: the Moebius transformation
//
//     f(z) = (c1 + c2*z) / (c3 + c4*z)
//
// stored as the coefficient 4-vector [c1 c2 c3 c4]. Infinity is a valid
// input and output, both for the three-point constructor and for eval.
//
// `moebius3` (moebius3.hpp) is the all-finite fast path of the three-point
// constructor, kept as a separate free function because crspread/crgather
// call it on hot paths and never need the infinite branches.
//
// The pretty-printing methods of the MATLAB class (char/disp/display) are
// not ported.
class Moebius {
public:
    // MATLAB MOEBIUS with no arguments: empty coefficients. Evaluating such
    // a map is not meaningful; this exists so `inv`/arithmetic can build up
    // a result the way the .m sources do.
    Moebius();

    // MATLAB MOEBIUS([C1 C2 C3 C4]) and MOEBIUS(C1,C2,C3,C4).
    explicit Moebius(const std::array<std::complex<double>, 4>& coeff);
    Moebius(std::complex<double> c1, std::complex<double> c2, std::complex<double> c3,
            std::complex<double> c4);

    // MATLAB MOEBIUS(Z,W): the transformation carrying the 3-vector z to w.
    // Infinite entries are allowed in either (but the MATLAB source assumes
    // at most one per vector, and so does this port).
    Moebius(const Eigen::Vector3cd& z, const Eigen::Vector3cd& w);

    // Port of @moebius/eval.m. Infinite inputs are mapped to c2/c4, and any
    // result that comes out NaN (including a vanishing denominator, which
    // MATLAB detects with the |den| < 3*eps test) is reported as Inf.
    Eigen::VectorXcd eval(const Eigen::VectorXcd& z) const;
    std::complex<double> eval(std::complex<double> z) const;

    // Port of @moebius/diff.m: (c2*c3 - c1*c4) / (c4*z + c3)^2.
    Eigen::VectorXcd diff(const Eigen::VectorXcd& z) const;
    std::complex<double> diff(std::complex<double> z) const;

    // Port of @moebius/inv.m.
    Moebius inverse() const;

    // Port of @moebius/normal.m: divide through by the denominator's
    // constant term, or by its linear coefficient when the constant is zero.
    Moebius normal() const;

    // Port of @moebius/subsref.m's M1(M2) form: the map z -> this(inner(z)).
    Moebius compose(const Moebius& inner) const;

    // Ports of @moebius/{plus,minus,mtimes,mrdivide,uminus,uplus}.m. All are
    // scalar operations; MATLAB errors on any other operand type.
    Moebius plus(std::complex<double> s) const;
    Moebius minus(std::complex<double> s) const;
    Moebius times(std::complex<double> s) const;
    Moebius dividedBy(std::complex<double> s) const;
    // MATLAB's `s / M`: reciprocate the map, then scale by s.
    static Moebius divide(std::complex<double> s, const Moebius& m);
    Moebius negated() const;

    // Port of @moebius/double.m.
    const std::array<std::complex<double>, 4>& coeff() const { return coeff_; }

    // The three-point data, when the map was built from a (z,w) pair. Stored
    // in the renumbered order the MATLAB constructor leaves them in, since
    // inv() swaps them.
    bool hasPoints() const { return hasPoints_; }
    const Eigen::Vector3cd& source() const { return source_; }
    const Eigen::Vector3cd& image() const { return image_; }

private:
    std::array<std::complex<double>, 4> coeff_;
    Eigen::Vector3cd source_;
    Eigen::Vector3cd image_;
    bool hasPoints_ = false;
};

inline Moebius operator+(const Moebius& m, std::complex<double> s) { return m.plus(s); }
inline Moebius operator+(std::complex<double> s, const Moebius& m) { return m.plus(s); }
inline Moebius operator-(const Moebius& m, std::complex<double> s) { return m.minus(s); }
inline Moebius operator*(const Moebius& m, std::complex<double> s) { return m.times(s); }
inline Moebius operator*(std::complex<double> s, const Moebius& m) { return m.times(s); }
inline Moebius operator/(const Moebius& m, std::complex<double> s) { return m.dividedBy(s); }
inline Moebius operator/(std::complex<double> s, const Moebius& m) { return Moebius::divide(s, m); }
inline Moebius operator-(const Moebius& m) { return m.negated(); }

}  // namespace sctoolbox
