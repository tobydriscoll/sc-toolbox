#include "sctoolbox/annulusmap_internal.hpp"

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>

#include "sctoolbox/gaussj.hpp"
#include "sctoolbox/nesolve.hpp"

namespace sctoolbox {

namespace {

constexpr double kPi = 3.14159265358979323846;

// 0-indexed [nodeStart, weightStart) extraction matching wqsum.m's
// iwt1/iwt2/ioffst index arithmetic (see annulusmap_internal.hpp).
std::pair<int, int> qworkBlockStart(int nptq, int M, int N, int kwa, int ic) {
    int nodeStart;
    if (kwa == 0) {
        nodeStart = nptq * (M + N);
    } else {
        nodeStart = nptq * (ic * M + kwa - 1);
    }
    const int ioffst = nptq * (M + N + 1);
    return {nodeStart, nodeStart + ioffst};
}

}  // namespace

Eigen::VectorXd qinit(const AnnulusData& dataz, int nptq) {
    const int M = dataz.M;
    const int N = dataz.N;
    Eigen::VectorXd qwork = Eigen::VectorXd::Zero(2 * nptq * (M + N + 1));

    for (int K = 1; K <= M + N; ++K) {
        const int inodes = nptq * (K - 1);  // 0-indexed
        const int iwts = nptq * (M + N + K);
        double alpha;
        if (K <= M) {
            alpha = dataz.ALFA0(K - 1) - 1.0;
            if (dataz.ALFA0(K - 1) > 0.0) {
                Eigen::VectorXd nodes, wts;
                gaussj(nptq, 0.0, alpha, nodes, wts);
                qwork.segment(inodes, nptq) = nodes;
                qwork.segment(iwts, nptq) = wts;
            } else {
                qwork.segment(inodes, nptq).setZero();
                qwork.segment(iwts, nptq).setZero();
            }
        } else {
            alpha = dataz.ALFA1(K - M - 1) - 1.0;
            Eigen::VectorXd nodes, wts;
            gaussj(nptq, 0.0, alpha, nodes, wts);
            qwork.segment(inodes, nptq) = nodes;
            qwork.segment(iwts, nptq) = wts;
        }
        for (int j = 0; j < nptq; ++j) {
            qwork(iwts + j) = qwork(iwts + j) * std::pow(1.0 + qwork(inodes + j), -alpha);
        }
    }

    const int inodes = nptq * (M + N);
    const int iwts = nptq * (2 * (M + N) + 1);
    Eigen::VectorXd nodes, wts;
    gaussj(nptq, 0.0, 0.0, nodes, wts);
    qwork.segment(inodes, nptq) = nodes;
    qwork.segment(iwts, nptq) = wts;

    return qwork;
}

ThetaParam4 thdata(double u) {
    ThetaParam4 p;
    if (u >= 0.63) {
        p.small = false;
        p.closedForm = (u >= 0.94);
        p.VARY = std::exp(kPi * kPi / std::log(u));
        p.DLAM = -std::log(u) / kPi;
        return p;
    }
    p.small = true;
    p.closedForm = false;
    if (u < 0.06) p.IU = 3;
    else if (u < 0.19) p.IU = 4;
    else if (u < 0.33) p.IU = 5;
    else if (u < 0.45) p.IU = 6;
    else if (u < 0.55) p.IU = 7;
    else p.IU = 8;

    p.UARY.resize(p.IU);
    for (int k = 1; k <= p.IU; ++k) p.UARY(k - 1) = std::pow(u, static_cast<double>(k) * k);
    return p;
}

std::complex<double> wtheta(double u, const ThetaParam4& param4, std::complex<double> w) {
    if (param4.small) {
        std::complex<double> sum(1.0, 0.0);
        std::complex<double> negw = -w;
        std::complex<double> pk = negw;
        for (int k = 1; k <= param4.IU; ++k) {
            sum += param4.UARY(k - 1) * (pk + 1.0 / pk);
            pk *= negw;
        }
        return sum;
    }
    const std::complex<double> wt = std::complex<double>(0.0, -1.0) * std::log(-w);
    if (param4.closedForm) {
        return std::exp(-0.25 * (wt * wt) / (kPi * param4.DLAM)) / std::sqrt(param4.DLAM);
    }
    const std::complex<double> bracket = 1.0 + 2.0 * param4.VARY * std::cosh(wt / param4.DLAM);
    return std::exp(-0.25 * (wt * wt) / (kPi * param4.DLAM)) * (bracket / std::sqrt(param4.DLAM));
}

Eigen::VectorXcd wprod(const Eigen::VectorXcd& w, double u, const Eigen::VectorXcd& uw0,
                       const Eigen::VectorXcd& u_w1, const AnnulusData& dataz, const ThetaParam4& param4) {
    const int n = static_cast<int>(w.size());
    Eigen::VectorXcd result(n);
    for (int i = 0; i < n; ++i) {
        std::complex<double> acc(0.0, 0.0);
        for (int j = 0; j < uw0.size(); ++j) {
            acc += std::log(wtheta(u, param4, w(i) / uw0(j))) * (dataz.ALFA0(j) - 1.0);
        }
        for (int j = 0; j < u_w1.size(); ++j) {
            acc += std::log(wtheta(u, param4, w(i) * u_w1(j))) * (dataz.ALFA1(j) - 1.0);
        }
        result(i) = std::exp(acc);
    }
    return result;
}

std::complex<double> wqsum(std::complex<double> wa, double phia, int kwa, int ic, std::complex<double> wb,
                           double phib, double radius, double u, const Eigen::VectorXcd& w0,
                           const Eigen::VectorXcd& w1, int nptq, const Eigen::VectorXd& qwork, int linearc,
                           const AnnulusData& dataz, const ThetaParam4& param4) {
    const auto [nodeStart, weightStart] = qworkBlockStart(nptq, dataz.M, dataz.N, kwa, ic);
    const Eigen::VectorXcd uw0 = u * w0;
    Eigen::VectorXcd u_w1(w1.size());
    for (int j = 0; j < w1.size(); ++j) u_w1(j) = u / w1(j);

    if (linearc == 1) {
        const double pwh = (phib - phia) / 2.0;
        const double pwc = (phib + phia) / 2.0;
        Eigen::VectorXcd w(nptq);
        for (int j = 0; j < nptq; ++j) {
            w(j) = radius * std::exp(std::complex<double>(0.0, pwc + pwh * qwork(nodeStart + j)));
        }
        const Eigen::VectorXcd prod = wprod(w, u, uw0, u_w1, dataz, param4);
        std::complex<double> sum(0.0, 0.0);
        for (int j = 0; j < nptq; ++j) sum += qwork(weightStart + j) * w(j) * prod(j);
        return std::complex<double>(0.0, pwh) * sum;
    } else {
        const std::complex<double> wh = (wb - wa) / 2.0;
        const std::complex<double> wc = (wa + wb) / 2.0;
        Eigen::VectorXcd w(nptq);
        for (int j = 0; j < nptq; ++j) w(j) = wc + wh * qwork(nodeStart + j);
        const Eigen::VectorXcd prod = wprod(w, u, uw0, u_w1, dataz, param4);
        std::complex<double> sum(0.0, 0.0);
        for (int j = 0; j < nptq; ++j) sum += qwork(weightStart + j) * prod(j);
        return wh * sum;
    }
}

namespace {

double nearestVertexDistance(std::complex<double> w, const Eigen::VectorXcd& w0, const Eigen::VectorXcd& w1) {
    double dmin = 2.0;
    for (int i = 0; i < w0.size(); ++i) {
        double d = std::abs(w - w0(i));
        if (d < 1e-15) d = 0.0;
        if (d != 0.0 && d < dmin) dmin = d;
    }
    for (int i = 0; i < w1.size(); ++i) {
        double d = std::abs(w - w1(i));
        if (d < 1e-15) d = 0.0;
        if (d != 0.0 && d < dmin) dmin = d;
    }
    return dmin;
}

}  // namespace

std::complex<double> wquad1(std::complex<double> wa, double phia, int kwa, int ic, std::complex<double> wb,
                            double phib, double radius, double u, const Eigen::VectorXcd& w0,
                            const Eigen::VectorXcd& w1, int nptq, const Eigen::VectorXd& qwork, int linearc,
                            const AnnulusData& dataz, const ThetaParam4& param4) {
    if (std::abs(wa - wb) == 0.0) return std::complex<double>(0.0, 0.0);

    if (linearc == 1) {
        const double wdist = std::min(2.0, nearestVertexDistance(wa, w0, w1));
        double r = std::min(1.0, wdist / std::abs(wb - wa));
        double phaa = phia + r * (phib - phia);
        std::complex<double> waa = radius * std::exp(std::complex<double>(0.0, phaa));
        std::complex<double> result = wqsum(wa, phia, kwa, ic, waa, phaa, radius, u, w0, w1, nptq, qwork, linearc,
                                            dataz, param4);
        int iters = 0;
        while (r != 1.0) {
            if (++iters > 10000) throw std::runtime_error("wquad1: circular-arc subdivision did not terminate");
            const double wdist2 = std::min(2.0, nearestVertexDistance(waa, w0, w1));
            r = std::min(1.0, wdist2 / std::abs(waa - wb));
            const double phbb = phaa + r * (phib - phaa);
            const std::complex<double> wbb = radius * std::exp(std::complex<double>(0.0, phbb));
            result += wqsum(waa, phaa, 0, 2, wbb, phbb, radius, u, w0, w1, nptq, qwork, linearc, dataz, param4);
            phaa = phbb;
            waa = wbb;
        }
        return result;
    } else {
        const double wdist = std::min(2.0, nearestVertexDistance(wa, w0, w1));
        double r = std::min(1.0, wdist / std::abs(wb - wa));
        std::complex<double> waa = wa + r * (wb - wa);
        std::complex<double> result = wqsum(wa, 0.0, kwa, ic, waa, 0.0, 0.0, u, w0, w1, nptq, qwork, linearc, dataz,
                                            param4);
        if (r != 1.0) {
            // MATLAB's continuation loop here (wquad1.m line 69) contains a
            // self-indexing bug (`d(d(d ~= 0))`) that would error if MATLAB
            // ever executed it; this path is therefore unreachable in
            // practice and is not ported, matching this codebase's
            // established precedent (e.g. crsplit's deferred mesh surgery).
            throw std::runtime_error(
                "wquad1: line-segment continuation not supported (matches an unreachable/buggy MATLAB branch)");
        }
        return result;
    }
}

std::complex<double> wquad(std::complex<double> wa, double phia, int kwa, int ica, std::complex<double> wb,
                           double phib, int kwb, int icb, double radius, double u, const Eigen::VectorXcd& w0,
                           const Eigen::VectorXcd& w1, int nptq, const Eigen::VectorXd& qwork, int linearc,
                           int ievl, const AnnulusData& dataz, const ThetaParam4& param4) {
    std::complex<double> wmid, wmida, wmidb;
    double phmid = 0.0, phmida = 0.0, phmidb = 0.0;

    if (linearc == 0) {
        wmid = (wa + wb) / 2.0;
        wmida = (wa + wmid) / 2.0;
        wmidb = (wb + wmid) / 2.0;
    } else {
        if (ievl != 1) {
            if (phib < phia) phia -= 2.0 * kPi;
        }
        phmid = (phia + phib) / 2.0;
        wmid = radius * std::exp(std::complex<double>(0.0, phmid));
        phmida = (phia + phmid) / 2.0;
        wmida = radius * std::exp(std::complex<double>(0.0, phmida));
        phmidb = (phib + phmid) / 2.0;
        wmidb = radius * std::exp(std::complex<double>(0.0, phmidb));
    }

    const std::complex<double> wa_to_wmida =
        wquad1(wa, phia, kwa, ica, wmida, phmida, radius, u, w0, w1, nptq, qwork, linearc, dataz, param4);
    const std::complex<double> wmid_to_wmida =
        wquad1(wmid, phmid, 0, 2, wmida, phmida, radius, u, w0, w1, nptq, qwork, linearc, dataz, param4);
    const std::complex<double> wqa = wa_to_wmida - wmid_to_wmida;

    const std::complex<double> wb_to_wmidb =
        wquad1(wb, phib, kwb, icb, wmidb, phmidb, radius, u, w0, w1, nptq, qwork, linearc, dataz, param4);
    const std::complex<double> wmid_to_wmidb =
        wquad1(wmid, phmid, 0, 2, wmidb, phmidb, radius, u, w0, w1, nptq, qwork, linearc, dataz, param4);
    const std::complex<double> wqb = wb_to_wmidb - wmid_to_wmidb;

    return wqa - wqb;
}

void xwtran(const Eigen::VectorXd& x, Eigen::VectorXcd& w0, Eigen::VectorXcd& w1, Eigen::VectorXd& phi0,
            Eigen::VectorXd& phi1, const AnnulusData& dataz, double& u, std::complex<double>& c) {
    const int M = dataz.M;
    const int N = dataz.N;
    // x is 0-indexed here; x(k) in MATLAB (1-indexed) == x(k-1) here.
    const double x1 = x(0);
    if (std::abs(x1) <= 1e-14) {
        u = 0.50;
    } else {
        double uu = (x1 - 2.0 - std::sqrt(0.9216 * x1 * x1 + 4.0)) / (2.0 * x1);
        u = (0.0196 * x1 - 1.0) / (uu * x1);
    }

    c = std::complex<double>(x(1), x(2));

    const double xn3 = x(N + 2);  // MATLAB x(N+3)
    if (std::abs(xn3) <= 1e-14) {
        phi1(N - 1) = 0.0;
    } else {
        const double ph = (1.0 + std::sqrt(1.0 + kPi * kPi * xn3 * xn3)) / xn3;
        phi1(N - 1) = (kPi * kPi) / ph;
    }

    double dph = 1.0;
    double phsum = dph;
    for (int k = 1; k <= N - 1; ++k) {
        dph = dph / std::exp(x(2 + k));  // MATLAB x(3+k)
        phsum += dph;
    }
    dph = 2.0 * kPi / phsum;
    phi1(0) = phi1(N - 1) + dph;
    w1(0) = u * std::complex<double>(std::cos(phi1(0)), std::sin(phi1(0)));
    w1(N - 1) = u * std::complex<double>(std::cos(phi1(N - 1)), std::sin(phi1(N - 1)));
    phsum = phi1(0);
    for (int k = 1; k <= N - 2; ++k) {
        dph = dph / std::exp(x(2 + k));
        phsum += dph;
        phi1(k) = phsum;
        w1(k) = u * std::complex<double>(std::cos(phsum), std::sin(phsum));
    }

    dph = 1.0;
    phsum = dph;
    for (int k = 1; k <= M - 1; ++k) {
        dph = dph / std::exp(x(N + 2 + k));  // MATLAB x(N+3+k)
        phsum += dph;
    }
    dph = 2.0 * kPi / phsum;
    phsum = dph;
    phi0(0) = dph;
    w0(0) = std::complex<double>(std::cos(dph), std::sin(dph));
    for (int k = 1; k <= M - 2; ++k) {
        dph = dph / std::exp(x(N + 2 + k));
        phsum += dph;
        phi0(k) = phsum;
        w0(k) = std::complex<double>(std::cos(phsum), std::sin(phsum));
    }
}

std::pair<int, int> nearw(std::complex<double> w, const Eigen::VectorXcd& w0, const Eigen::VectorXcd& w1,
                          const AnnulusData& dataz) {
    double dist = 2.0;
    int knear = -1;
    for (int i = 0; i < w0.size(); ++i) {
        if (dataz.ALFA0(i) > 0.0) {
            const double d = std::abs(w - w0(i));
            if (d < dist) {
                dist = d;
                knear = i;
            }
        }
    }
    int inear = 0;
    for (int i = 0; i < w1.size(); ++i) {
        const double d = std::abs(w - w1(i));
        if (d < dist) {
            dist = d;
            knear = i;
            inear = 1;
        }
    }
    return {knear, inear};
}

std::pair<int, int> nearz(std::complex<double> z, const AnnulusData& dataz) {
    int inz = 2;
    double dist = 99.0;
    int knz = -1;
    for (int i = 0; i < dataz.Z0.size(); ++i) {
        if (dataz.ALFA0(i) > 0.0) {
            const double d = std::abs(z - dataz.Z0(i));
            if (d < dist) {
                dist = d;
                knz = i;
                inz = 0;
            }
        }
    }
    for (int i = 0; i < dataz.Z1.size(); ++i) {
        const double d = std::abs(z - dataz.Z1(i));
        if (d < dist) {
            dist = d;
            knz = i;
            inz = 1;
        }
    }
    return {knz, inz};
}

std::complex<double> zdsc(std::complex<double> ww, int kww, int ic, double u, std::complex<double> c,
                          const Eigen::VectorXcd& w0, const Eigen::VectorXcd& w1, const Eigen::VectorXd& phi0,
                          const Eigen::VectorXd& phi1, int nptq, const Eigen::VectorXd& qwork, int iopt,
                          const AnnulusData& dataz) {
    const ThetaParam4 param4 = thdata(u);
    int ibd = 1;
    std::complex<double> ww0 = ww;

    if (iopt == 1) {
        if (std::abs(std::abs(ww0) - 1.0) == 0.0) {
            ww0 = (1.0 + u) * ww0 / 2.0;
            ibd = 2;
        }
    }

    const auto [knear, inear] = nearw(ww0, w0, w1, dataz);
    std::complex<double> za, wa;
    if (inear == 0) {
        za = dataz.Z0(knear);
        wa = w0(knear);
    } else {
        za = dataz.Z1(knear);
        wa = w1(knear);
    }

    if (iopt != 1) {
        return za + c * wquad(wa, 0.0, knear + 1, inear, ww, 0.0, kww, ic, 0.0, u, w0, w1, nptq, qwork, 0, 1, dataz,
                              param4);
    }

    double phiww0 = (ww0.imag() >= 0.0) ? std::arg(ww0) : std::arg(ww0) + 2.0 * kPi;
    const double dww0 = std::abs(ww0);
    double phiwb = (inear == 0) ? phi0(knear) : phi1(knear);
    std::complex<double> wb = dww0 * std::exp(std::complex<double>(0.0, phiwb));

    std::complex<double> wint1(0.0, 0.0);
    if (std::abs(wb - wa) != 0.0) {
        wint1 = wquad(wa, 0.0, knear + 1, inear, wb, 0.0, 0, 2, 0.0, u, w0, w1, nptq, qwork, 0, 1, dataz, param4);
    }

    if (std::abs(wb - ww0) == 0.0) {
        std::complex<double> z_dsc = za + c * wint1;
        if (ibd == 2) {
            z_dsc += c * wquad(ww0, 0.0, 0, 2, ww, 0.0, kww, ic, 0.0, u, w0, w1, nptq, qwork, 0, 1, dataz, param4);
        }
        return z_dsc;
    }

    if (std::abs(phiwb - 2.0 * kPi - phiww0) < std::abs(phiwb - phiww0)) phiwb -= 2.0 * kPi;
    if (std::abs(phiww0 - 2.0 * kPi - phiwb) < std::abs(phiwb - phiww0)) phiww0 -= 2.0 * kPi;

    std::complex<double> wint2;
    if (std::abs(wb - wa) == 0.0) {
        wint2 = wquad(wb, phiwb, knear + 1, inear, ww0, phiww0, kww, ic, dww0, u, w0, w1, nptq, qwork, 1, 1, dataz,
                      param4);
        std::complex<double> z_dsc = za + c * (wint1 + wint2);
        if (ibd == 2) {
            z_dsc += c * wquad(ww0, 0.0, 0, 2, ww, 0.0, kww, ic, 0.0, u, w0, w1, nptq, qwork, 0, 1, dataz, param4);
        }
        return z_dsc;
    }

    wint2 = wquad(wb, phiwb, 0, 2, ww0, phiww0, 0, 2, dww0, u, w0, w1, nptq, qwork, 1, 1, dataz, param4);
    std::complex<double> z_dsc = za + c * (wint1 + wint2);
    if (ibd == 2) {
        z_dsc += c * wquad(ww0, 0.0, 0, 2, ww, 0.0, kww, ic, 0.0, u, w0, w1, nptq, qwork, 0, 1, dataz, param4);
    }
    return z_dsc;
}

std::complex<double> wdsc(std::complex<double> zz, double u, std::complex<double> c, const Eigen::VectorXcd& w0,
                          const Eigen::VectorXcd& w1, const Eigen::VectorXd& phi0, const Eigen::VectorXd& phi1,
                          int nptq, const Eigen::VectorXd& qwork, double eps, int iopt, const AnnulusData& dataz) {
    const ThetaParam4 param4 = thdata(u);
    const Eigen::VectorXcd uw0 = u * w0;
    Eigen::VectorXcd u_w1(w1.size());
    for (int j = 0; j < w1.size(); ++j) u_w1(j) = u / w1(j);

    const auto [knz, inz] = nearz(zz, dataz);
    if (inz == 2) {
        throw std::runtime_error("wdsc: inverse evaluation failed (no nearby vertex found)");
    }

    std::complex<double> wi;
    if (inz == 0) {
        if (knz >= 1) {
            wi = (w0(knz) + w0(knz - 1)) / 2.0;
        } else {
            wi = (w0(knz) + w0(dataz.M - 1)) / 2.0;
        }
    } else {
        if (knz >= 1) {
            wi = (w1(knz) + w1(knz - 1)) / 2.0;
        } else {
            wi = (w0(knz) + w1(dataz.N - 1)) / 2.0;
        }
        wi = u * wi / std::abs(wi);
    }

    const std::complex<double> zs = zdsc(wi, 0, 2, u, c, w0, w1, phi0, phi1, nptq, qwork, 1, dataz);

    for (int k = 0; k < 20; ++k) {
        wi = wi + (zz - zs) / (20.0 * c * wprod(Eigen::VectorXcd::Constant(1, wi), u, uw0, u_w1, dataz, param4)(0));
    }
    if (std::abs(wi) > 1.0) wi = wi / (std::abs(wi) + std::abs(1.0 - std::abs(wi)));
    if (std::abs(wi) < u) wi = u * wi * (std::abs(wi) + std::abs(u - std::abs(wi))) / std::abs(wi);

    int it = 1;
    std::complex<double> wfn = zz - zdsc(wi, 0, 2, u, c, w0, w1, phi0, phi1, nptq, qwork, iopt, dataz);
    while (std::abs(wfn) >= eps && it <= 15) {
        wi = wi + wfn / (c * wprod(Eigen::VectorXcd::Constant(1, wi), u, uw0, u_w1, dataz, param4)(0));
        if (std::abs(wi) > 1.0) wi = wi / (std::abs(wi) + std::abs(1.0 - std::abs(wi)));
        if (std::abs(wi) < u) wi = u * wi * (std::abs(wi) + std::abs(u - std::abs(wi))) / std::abs(wi);
        ++it;
        wfn = zz - zdsc(wi, 0, 2, u, c, w0, w1, phi0, phi1, nptq, qwork, iopt, dataz);
    }

    std::complex<double> zz1;
    std::complex<double> w_dsc;
    if (std::abs(wfn) < eps) {
        w_dsc = wi;
        zz1 = zdsc(wi, 0, 2, u, c, w0, w1, phi0, phi1, nptq, qwork, 1, dataz);
    } else {
        zz1 = std::complex<double>(1.0, 1.0) + zz;
    }

    if (std::abs(wi) >= u && std::abs(wi) <= 1.0 && std::abs(zz - zz1) <= 1e-3) {
        return w_dsc;
    }

    AnnulusData temp = dataz;
    if (inz == 0) {
        temp.Z0(knz) = zz;
    } else {
        temp.Z1(knz) = zz;
    }

    return wdsc(zz, u, c, w0, w1, phi0, phi1, nptq, qwork, eps, iopt, temp);
}

namespace {

// Port of @annulusmap/private/dscfun.m, restricted to the bounded outer
// polygon case (ishape==0); the unbounded/truncated-polygon branches
// (ishape==1, the `nshape`/`ind` bookkeeping) are not ported since the
// C++ AnnulusMap constructor only supports bounded outer polygons.
Eigen::VectorXd dscfun(const Eigen::VectorXd& x, int nptq, const Eigen::VectorXd& qwork, int linearc,
                       const AnnulusData& dataz) {
    const int M = dataz.M;
    const int N = dataz.N;

    // xwtran never touches w0(M-1)/phi0(M-1) (the one prevertex fixed at the
    // start of dscsolv); they must be pre-set here too, matching MATLAB's
    // dscfun.m which captures the same fixed w0/phi0 array via its closure.
    Eigen::VectorXcd w0 = Eigen::VectorXcd::Zero(M);
    Eigen::VectorXcd w1 = Eigen::VectorXcd::Zero(N);
    Eigen::VectorXd phi0 = Eigen::VectorXd::Zero(M);
    Eigen::VectorXd phi1 = Eigen::VectorXd::Zero(N);
    w0(M - 1) = std::complex<double>(1.0, 0.0);
    phi0(M - 1) = 0.0;

    double u;
    std::complex<double> c;
    xwtran(x, w0, w1, phi0, phi1, dataz, u, c);
    const ThetaParam4 param4 = thdata(u);

    Eigen::VectorXd fval(M + N + 2);

    const std::complex<double> win1 =
        dataz.Z1(0) - dataz.Z1(N - 1) -
        c * wquad(w1(N - 1), phi1(N - 1), N, 1, w1(0), phi1(0), 1, 1, u, u, w0, w1, nptq, qwork, linearc, 2, dataz,
                  param4);
    fval(0) = win1.real();
    fval(1) = win1.imag();

    for (int I = 1; I <= N - 1; ++I) {
        const std::complex<double> wint1 = wquad(w1(I - 1), phi1(I - 1), I, 1, w1(I), phi1(I), I + 1, 1, u, u, w0,
                                                 w1, nptq, qwork, linearc, 2, dataz, param4);
        fval(I + 1) = std::abs(dataz.Z1(I) - dataz.Z1(I - 1)) - std::abs(c * wint1);
    }

    const double test1 = std::cos(phi1(N - 1));
    std::complex<double> win2;
    if (test1 >= u) {
        win2 = wquad(w0(M - 1), 0.0, M, 0, w1(N - 1), 0.0, N, 1, 0.0, u, w0, w1, nptq, qwork, 0, 2, dataz, param4);
    } else {
        const double wx = u;
        const std::complex<double> wline =
            wquad(w0(M - 1), 0.0, M, 0, wx, 0.0, 0, 2, 0.0, u, w0, w1, nptq, qwork, 0, 2, dataz, param4);
        std::complex<double> warc;
        if (phi1(N - 1) <= 0.0) {
            warc = wquad(w1(N - 1), phi1(N - 1), N, 1, wx, 0.0, 0, 2, u, u, w0, w1, nptq, qwork, 1, 2, dataz, param4);
            win2 = wline - warc;
        } else {
            warc = wquad(wx, 0.0, 0, 2, w1(N - 1), phi1(N - 1), N, 1, u, u, w0, w1, nptq, qwork, 1, 2, dataz, param4);
            win2 = wline + warc;
        }
    }
    const std::complex<double> win2full = dataz.Z1(N - 1) - dataz.Z0(M - 1) - c * win2;
    fval(N + 1) = win2full.real();
    fval(N + 2) = win2full.imag();

    const std::complex<double> win3 =
        dataz.Z0(0) - dataz.Z0(M - 1) -
        c * wquad(w0(M - 1), phi0(M - 1), M, 0, w0(0), phi0(0), 1, 0, 1.0, u, w0, w1, nptq, qwork, linearc, 2, dataz,
                  param4);
    fval(N + 3) = win3.real();
    fval(N + 4) = win3.imag();

    if (M == 3) return fval;

    for (int J = 1; J <= M - 3; ++J) {
        const std::complex<double> wint2 = wquad(w0(J - 1), phi0(J - 1), J, 0, w0(J), phi0(J), J + 1, 0, 1.0, u, w0,
                                                 w1, nptq, qwork, linearc, 2, dataz, param4);
        fval(N + 4 + J) = std::abs(dataz.Z0(J) - dataz.Z0(J - 1)) - std::abs(c * wint2);
    }
    return fval;
}

}  // namespace

DscParams dscsolv(int nptq, const Eigen::VectorXd& qwork, bool ishape, int linearc, const AnnulusData& dataz) {
    if (ishape) {
        throw std::runtime_error("dscsolv: unbounded/truncated outer polygon (ishape==1) not supported");
    }
    const int M = dataz.M;
    const int N = dataz.N;

    Eigen::VectorXcd w0 = Eigen::VectorXcd::Zero(M);
    Eigen::VectorXcd w1 = Eigen::VectorXcd::Zero(N);
    Eigen::VectorXd phi0 = Eigen::VectorXd::Zero(M);
    Eigen::VectorXd phi1 = Eigen::VectorXd::Zero(N);
    w0(M - 1) = std::complex<double>(1.0, 0.0);
    phi0(M - 1) = 0.0;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(M + N + 2);
    x(1) = 0.0;
    x(2) = 0.0;

    // INITIAL GUESS (iguess==0) -- the only branch the constructor uses.
    x(0) = 1.0 / 0.5 - 1.0 / 0.46;
    const double ave_in = 2.0 * kPi / N;
    for (int k = 1; k <= N - 2; ++k) {
        x(2 + k) = std::log((ave_in + 0.0001 * k) / (ave_in + 0.0001 * (k + 1)));
    }
    x(N + 1) = std::log((ave_in + 0.0001 * (N - 1)) / (2.0 * kPi - (N - 1) * (ave_in + N * 0.00005)));
    x(N + 2) = 1.0 / (4.0 - 0.1) - 1.0 / (4.0 + 0.1);
    const double ave_out = 2.0 * kPi / M;
    for (int t = 0; t <= M - 3; ++t) {
        x(N + 3 + t) = std::log((ave_out + 0.0001 * t) / (ave_out + 0.0001 * (t + 1)));
    }
    x(M + N + 1) = std::log((ave_out + 0.0001 * (M - 2)) / (2.0 * kPi - (M - 1) * (ave_out + (M - 2) * 0.00005)));

    // Calculate the initial guess x(2) & x(3) (0-indexed x(1),x(2)) to match
    // the choice of the other entries:
    double u;
    std::complex<double> c;
    xwtran(x, w0, w1, phi0, phi1, dataz, u, c);
    ThetaParam4 param4 = thdata(u);
    const std::complex<double> wint =
        wquad(w0(M - 1), 0.0, M, 0, w1(N - 1), 0.0, N, 1, 0.0, u, w0, w1, nptq, qwork, 0, 2, dataz, param4);
    const std::complex<double> c1 = (dataz.Z1(N - 1) - dataz.Z0(M - 1)) / wint;
    x(1) = c1.real();
    x(2) = c1.imag();

    const Fvec fvec = [&](const Eigen::VectorXd& xx) { return dscfun(xx, nptq, qwork, linearc, dataz); };

    Eigen::VectorXd details = Eigen::VectorXd::Zero(16);
    details(0) = 2.0;
    details(1) = 2.0;
    details(5) = 100.0 * (16 - 3);
    details(7) = 1e-8;
    details(8) = std::min(std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0), 1e-8 / 10.0);
    details(11) = nptq;

    const NesolveResult r = nesolve(fvec, x, details);
    x = r.xf;

    DscParams result;
    xwtran(x, w0, w1, phi0, phi1, dataz, result.u, result.c);
    result.w0 = w0;
    result.w1 = w1;
    result.phi0 = phi0;
    result.phi1 = phi1;
    return result;
}

}  // namespace sctoolbox
