#pragma once
#include <Eigen/Dense>
#include <complex>

namespace sctoolbox {

// Plain-data bundle mirroring the `dataz` struct threaded through
// @annulusmap/private/*.m (M, N, Z0, Z1, ALFA0, ALFA1). Z0/Z1/ALFA0/ALFA1
// are 0-indexed Eigen vectors of length M (resp. N), matching the MATLAB
// 1-indexed arrays of the same length.
struct AnnulusData {
    int M = 0;
    int N = 0;
    Eigen::VectorXcd Z0;
    Eigen::VectorXcd Z1;
    Eigen::VectorXd ALFA0;
    Eigen::VectorXd ALFA1;
};

// Port of @annulusmap/private/qinit.m. Returns the flat 0-indexed work array
// (length 2*nptq*(M+N+1)) of Gauss-Jacobi nodes/weights.
Eigen::VectorXd qinit(const AnnulusData& dataz, int nptq);

// Port of @annulusmap/private/thdata.m. Precomputed theta-function data that
// depends only on the inner radius u.
struct ThetaParam4 {
    bool small = true;  // true iff u < 0.63 (series-sum branch in wtheta)
    bool closedForm = false;  // true iff u >= 0.94 (closed-form sub-branch)
    double VARY = 0.0;
    double DLAM = 0.0;
    int IU = 0;
    Eigen::VectorXd UARY;  // length IU, UARY(k) = u^(k^2), 1-indexed k stored 0-indexed
};
ThetaParam4 thdata(double u);

// Port of @annulusmap/private/wtheta.m, evaluated at a single point (the
// MATLAB version is array-valued; callers here loop themselves).
std::complex<double> wtheta(double u, const ThetaParam4& param4, std::complex<double> w);

// Port of @annulusmap/private/wprod.m (the D-SC integrand).
Eigen::VectorXcd wprod(const Eigen::VectorXcd& w, double u, const Eigen::VectorXcd& uw0,
                       const Eigen::VectorXcd& u_w1, const AnnulusData& dataz, const ThetaParam4& param4);

// Port of @annulusmap/private/wqsum.m.
std::complex<double> wqsum(std::complex<double> wa, double phia, int kwa, int ic, std::complex<double> wb,
                           double phib, double radius, double u, const Eigen::VectorXcd& w0,
                           const Eigen::VectorXcd& w1, int nptq, const Eigen::VectorXd& qwork, int linearc,
                           const AnnulusData& dataz, const ThetaParam4& param4);

// Port of @annulusmap/private/wquad1.m.
std::complex<double> wquad1(std::complex<double> wa, double phia, int kwa, int ic, std::complex<double> wb,
                            double phib, double radius, double u, const Eigen::VectorXcd& w0,
                            const Eigen::VectorXcd& w1, int nptq, const Eigen::VectorXd& qwork, int linearc,
                            const AnnulusData& dataz, const ThetaParam4& param4);

// Port of @annulusmap/private/wquad.m.
std::complex<double> wquad(std::complex<double> wa, double phia, int kwa, int ica, std::complex<double> wb,
                           double phib, int kwb, int icb, double radius, double u, const Eigen::VectorXcd& w0,
                           const Eigen::VectorXcd& w1, int nptq, const Eigen::VectorXd& qwork, int linearc,
                           int ievl, const AnnulusData& dataz, const ThetaParam4& param4);

// D-SC parameters, mirroring annulusmap.m's u/c/w0/w1/phi0/phi1 properties.
struct DscParams {
    double u = 0.0;
    std::complex<double> c;
    Eigen::VectorXcd w0;
    Eigen::VectorXcd w1;
    Eigen::VectorXd phi0;
    Eigen::VectorXd phi1;
};

// Port of @annulusmap/private/xwtran.m. x is the unconstrained-parameter
// vector (1-indexed in MATLAB; here a 0-indexed Eigen::VectorXd of the same
// length, x(k) in MATLAB == x(k-1) here). w0/w1/phi0/phi1 are updated in
// place (matching MATLAB's in/out parameter reuse for the fixed M-th / N-th
// entries that xwtran never recomputes).
void xwtran(const Eigen::VectorXd& x, Eigen::VectorXcd& w0, Eigen::VectorXcd& w1, Eigen::VectorXd& phi0,
            Eigen::VectorXd& phi1, const AnnulusData& dataz, double& u, std::complex<double>& c);

// Port of @annulusmap/private/dscsolv.m. ishape must be false (bounded outer
// polygon); the unbounded/truncated-polygon branch (ishape==1) is not
// ported, matching this codebase's scope (no `truncate` constructor path).
DscParams dscsolv(int nptq, const Eigen::VectorXd& qwork, bool ishape, int linearc, const AnnulusData& dataz);

// Port of @annulusmap/private/nearw.m. Returns (knear, inear), 0-indexed
// knear into w0 (inear==0) or w1 (inear==1).
std::pair<int, int> nearw(std::complex<double> w, const Eigen::VectorXcd& w0, const Eigen::VectorXcd& w1,
                          const AnnulusData& dataz);

// Port of @annulusmap/private/nearz.m. Returns (knz, inz), 0-indexed knz
// into Z0 (inz==0) or Z1 (inz==1); inz==2 signals "no vertex found".
std::pair<int, int> nearz(std::complex<double> z, const AnnulusData& dataz);

// Port of @annulusmap/private/zdsc.m (forward map).
std::complex<double> zdsc(std::complex<double> ww, int kww, int ic, double u, std::complex<double> c,
                          const Eigen::VectorXcd& w0, const Eigen::VectorXcd& w1, const Eigen::VectorXd& phi0,
                          const Eigen::VectorXd& phi1, int nptq, const Eigen::VectorXd& qwork, int iopt,
                          const AnnulusData& dataz);

// Port of @annulusmap/private/wdsc.m (inverse map).
std::complex<double> wdsc(std::complex<double> zz, double u, std::complex<double> c, const Eigen::VectorXcd& w0,
                          const Eigen::VectorXcd& w1, const Eigen::VectorXd& phi0, const Eigen::VectorXd& phi1,
                          int nptq, const Eigen::VectorXd& qwork, double eps, int iopt, const AnnulusData& dataz);

}  // namespace sctoolbox
