#include "sctoolbox/deparam.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "sctoolbox/dabsquad.hpp"
#include "sctoolbox/dequad.hpp"
#include "sctoolbox/nesolve.hpp"
#include "sctoolbox/scqdata.hpp"

namespace sctoolbox {

namespace {

// cs = cumsum(cumprod([1;exp(-y)])); theta = 2*pi*cs(1:n-1)/cs(n).
Eigen::VectorXd yToTheta(const Eigen::VectorXd& y, int n) {
    Eigen::VectorXd v(n);
    v(0) = 1.0;
    for (int i = 0; i < n - 1; ++i) v(i + 1) = std::exp(-y(i));
    double prod = 1.0;
    for (int i = 0; i < n; ++i) {
        prod *= v(i);
        v(i) = prod;
    }
    Eigen::VectorXd cs(n);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += v(i);
        cs(i) = sum;
    }
    Eigen::VectorXd theta(n - 1);
    for (int i = 0; i < n - 1; ++i) theta(i) = 2.0 * M_PI * cs(i) / cs(n - 1);
    return theta;
}

Eigen::VectorXcd thetaToZ(const Eigen::VectorXd& theta, int n) {
    Eigen::VectorXcd z = Eigen::VectorXcd::Ones(n);
    for (int i = 0; i < n - 1; ++i) z(i) = std::exp(std::complex<double>(0.0, theta(i)));
    return z;
}

// Port of deparam.m's nested depfun residual function.
Eigen::VectorXd depfun(const Eigen::VectorXd& y, int n, const Eigen::VectorXd& beta, const Eigen::VectorXd& nmlen,
                       const Eigen::MatrixXd& qdat) {
    const Eigen::VectorXd theta = yToTheta(y, n);
    const Eigen::VectorXcd z = thetaToZ(theta, n);

    const int k = n - 2;
    Eigen::VectorXcd mid(k), z1(k), z2(k);
    std::vector<int> sing1(k), sing2(k);
    for (int i = 0; i < k; ++i) {
        mid(i) = std::exp(std::complex<double>(0.0, (theta(i) + theta(i + 1)) / 2.0));
        z1(i) = z(i);
        z2(i) = z(i + 1);
        sing1[i] = i + 1;
        sing2[i] = i + 2;
    }

    const Eigen::VectorXd I1 = dabsquad(z1, mid, sing1, z, beta, qdat);
    const Eigen::VectorXd I2 = dabsquad(z2, mid, sing2, z, beta, qdat);
    const Eigen::VectorXd ints = I1 + I2;

    const int nmain = (n > 3) ? (n - 3) : 0;
    Eigen::VectorXd F(nmain + 2);
    for (int i = 0; i < nmain; ++i) F(i) = std::abs(ints(i + 1)) / std::abs(ints(0)) - nmlen(i);

    std::complex<double> sumBetaOverZ(0.0, 0.0);
    for (int i = 0; i < n; ++i) sumBetaOverZ += beta(i) / z(i);
    const std::complex<double> res = -sumBetaOverZ / ints(0);

    F(nmain) = res.real();
    F(nmain + 1) = res.imag();
    return F;
}

}  // namespace

DeParamResult deparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, double tol, int method) {
    const int n = static_cast<int>(w.size());
    const int nqpts = std::max(static_cast<int>(std::ceil(-std::log10(tol))), 2);
    const Eigen::MatrixXd qdat = scqdata(beta, nqpts);

    Eigen::VectorXcd z(n);

    if (n == 2) {
        z(0) = -1.0;
        z(1) = 1.0;
    } else {
        Eigen::VectorXd len(n);
        len(0) = std::abs(w(0) - w(n - 1));
        for (int i = 1; i < n; ++i) len(i) = std::abs(w(i) - w(i - 1));

        Eigen::VectorXd nmlen(n - 3);
        for (int i = 0; i < n - 3; ++i) nmlen(i) = std::abs(len(i + 2) / len(1));

        const Eigen::VectorXd y0 = Eigen::VectorXd::Zero(n - 1);

        const Fvec fvec = [&](const Eigen::VectorXd& y) { return depfun(y, n, beta, nmlen, qdat); };

        Eigen::VectorXd details = Eigen::VectorXd::Zero(16);
        details(0) = 0.0;
        details(1) = static_cast<double>(method);
        details(5) = 100.0 * (n - 3);
        details(7) = tol;
        details(8) = std::min(std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0), tol / 10.0);
        details(11) = nqpts;

        const NesolveResult r = nesolve(fvec, y0, details);
        const Eigen::VectorXd theta = yToTheta(r.xf, n);
        z = thetaToZ(theta, n);
    }

    const std::complex<double> mid = z(0) * std::exp(std::complex<double>(0.0, std::arg(z(1) / z(0)) / 2.0));
    const std::vector<int> singAt1{1}, singAt2{2};
    const Eigen::VectorXcd I1 = dequad(z.segment(0, 1), Eigen::VectorXcd::Constant(1, mid), singAt1, z, beta, qdat);
    const Eigen::VectorXcd I2 = dequad(z.segment(1, 1), Eigen::VectorXcd::Constant(1, mid), singAt2, z, beta, qdat);
    const std::complex<double> c = (w(1) - w(0)) / (I1(0) - I2(0));

    return DeParamResult{z, c, qdat};
}

}  // namespace sctoolbox
