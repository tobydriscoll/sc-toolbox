#include "sctoolbox/nesolve.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sctoolbox/nechdcmp.hpp"
#include "sctoolbox/neqrdcmp.hpp"

namespace sctoolbox {
namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();

double sign1(double x) { return (x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0)) + (x == 0.0 ? 1.0 : 0.0); }

// +sctool/nefn.m
double nefn(const Eigen::VectorXd& xplus, const Eigen::VectorXd& SF, const Fvec& fvec,
            int& nofun, Eigen::VectorXd& FVplus) {
    FVplus = fvec(xplus);
    ++nofun;
    return 0.5 * (SF.array() * FVplus.array()).square().sum();
}

// +sctool/neinck.m (scale argument always absent in this port)
void neinck(const Eigen::VectorXd& x0, Eigen::VectorXd F0, Eigen::VectorXd& dout,
            Eigen::VectorXd& Sx, Eigen::VectorXd& SF, int& termcode) {
    const int n = static_cast<int>(x0.size());
    termcode = 0;

    if (dout(15) == 0.0) {
        Sx = Eigen::VectorXd::Ones(n);
        SF = Eigen::VectorXd::Ones(n);
    } else {
        Eigen::VectorXd x0safe = x0;
        for (int i = 0; i < n; ++i) if (x0safe(i) == 0.0) x0safe(i) = 1.0;
        for (int i = 0; i < n; ++i) if (F0(i) == 0.0) F0(i) = 1.0;
        Sx = x0safe.cwiseAbs().cwiseInverse();
        SF = F0.cwiseAbs().cwiseInverse();
    }

    dout(12) = (dout(11) <= 0.0) ? kEps : std::max(kEps, std::pow(10.0, -dout(11)));
    if (dout(12) > 0.01) {
        termcode = -2;
        return;
    }

    if (dout(1) <= 0.0) dout(1) = 1.0;
    if ((dout(1) == 2.0 || dout(1) == 3.0) && dout(6) <= 0.0) dout(6) = -1.0;

    if (dout(5) <= 1.0) dout(5) = 100.0;
    if (dout(7) <= 0.0) dout(7) = std::pow(kEps, 1.0 / 3.0);
    if (dout(8) <= 0.0) dout(8) = std::pow(kEps, 2.0 / 3.0);
    if (dout(9) <= 0.0) dout(9) = std::pow(kEps, 2.0 / 3.0);
    if (dout(10) <= 0.0) {
        // norm(diag(Sx)) (the matrix 2-norm) equals max(abs(Sx)).
        dout(10) = 1000.0 * std::max((Sx.array() * x0.array()).matrix().norm(), Sx.cwiseAbs().maxCoeff());
    }
}

// +sctool/nersolv.m
void nersolv(const Eigen::MatrixXd& M, const Eigen::VectorXd& M2, Eigen::VectorXd& b) {
    const int n = static_cast<int>(M.rows());
    b(n - 1) = b(n - 1) / M2(n - 1);
    for (int i = n - 2; i >= 0; --i) {
        b(i) = (b(i) - M.row(i).segment(i + 1, n - 1 - i) * b.segment(i + 1, n - 1 - i)) / M2(i);
    }
}

// +sctool/neqrsolv.m
void neqrsolv(const Eigen::MatrixXd& M, const Eigen::VectorXd& M1, const Eigen::VectorXd& M2,
              Eigen::VectorXd& b) {
    const int n = static_cast<int>(M.rows());
    for (int j = 0; j < n - 1; ++j) {
        const double tau = (M.col(j).segment(j, n - j).dot(b.segment(j, n - j))) / M1(j);
        b.segment(j, n - j) -= tau * M.col(j).segment(j, n - j);
    }
    nersolv(M, M2, b);
}

// +sctool/neconest.m
double neconest(const Eigen::MatrixXd& M, const Eigen::VectorXd& M2) {
    const int n = static_cast<int>(M.rows());
    Eigen::VectorXd p = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd pm = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

    Eigen::MatrixXd Rtri = M.triangularView<Eigen::Upper>();
    Rtri.diagonal() = M2;
    double est = Rtri.lpNorm<1>();

    x(0) = 1.0 / M2(0);
    for (int i = 1; i < n; ++i) p(i) = M(0, i) * x(0);

    for (int j = 1; j < n; ++j) {
        const double xp = (1.0 - p(j)) / M2(j);
        const double xm = (-1.0 - p(j)) / M2(j);
        double temp = std::abs(xp);
        double tempm = std::abs(xm);
        for (int i = j + 1; i < n; ++i) {
            pm(i) = p(i) + M(j, i) * xm;
            tempm += std::abs(p(i)) / std::abs(M2(i));
            p(i) = p(i) + M(j, i) * xp;
            temp += std::abs(p(i)) / std::abs(M2(i));
        }
        if (temp > tempm) {
            x(j) = xp;
        } else {
            x(j) = xm;
            for (int i = j + 1; i < n; ++i) p(i) = pm(i);
        }
    }

    est = est / x.lpNorm<1>();
    nersolv(M, M2, x);
    est = est * x.lpNorm<1>();
    return est;
}

struct ModelResult {
    Eigen::MatrixXd m;  // encoded QR (line-search) or Cholesky factor L (hook step)
    Eigen::MatrixXd h;  // approximate Hessian J'J, only populated for hook step
    Eigen::VectorXd sn;
};

// +sctool/nemodel.m
ModelResult nemodel(const Eigen::VectorXd& fc, const Eigen::MatrixXd& J, const Eigen::VectorXd& g,
                     const Eigen::VectorXd& sf, const Eigen::VectorXd& sx, int globmeth) {
    const int n = static_cast<int>(J.rows());
    ModelResult r;
    Eigen::MatrixXd m = sf.asDiagonal() * J;

    Eigen::VectorXd m1, m2;
    int sing = 0;
    neqrdcmp(m, m1, m2, sing);

    double est = 0.0;
    if (sing == 0) {
        for (int j = 1; j < n; ++j) m.col(j).head(j) /= sx(j);
        m2 = m2.cwiseQuotient(sx);
        est = neconest(m, m2);
    }

    if (sing == 1 || est > 1.0 / kEps || std::isnan(est)) {
        Eigen::MatrixXd h = J.transpose() * sf.asDiagonal();
        h = h * h.transpose();
        Eigen::VectorXd invSx = sx.cwiseInverse();
        double hnorm = (1.0 / sx(0)) * (h.row(0).cwiseAbs() * invSx)(0);
        for (int i = 1; i < n; ++i) {
            const double tem1 = (h.col(i).cwiseAbs().array() / sx.array()).sum();
            const double tem2 = (h.row(i).transpose().cwiseAbs().array() / sx.array()).sum();
            const double temp = (1.0 / sx(i)) / (tem1 + tem2);
            hnorm = std::max(temp, hnorm);
        }
        h += std::sqrt(n * kEps) * hnorm * sx.array().square().matrix().asDiagonal().toDenseMatrix();
        Eigen::MatrixXd L;
        double maxadd;
        nechdcmp(h, 0.0, L, maxadd);
        r.m = L;
        r.h = h;
        r.sn = -(L.transpose().triangularView<Eigen::Upper>().solve(
            L.triangularView<Eigen::Lower>().solve(g)));
    } else {
        for (int j = 1; j < n; ++j) m.col(j).head(j) *= sx(j);
        m2 = m2.cwiseProduct(sx);
        Eigen::VectorXd sn = -sf.cwiseProduct(fc);
        neqrsolv(m, m1, m2, sn);
        if (globmeth == 2 || globmeth == 3) {
            Eigen::MatrixXd mm = m.triangularView<Eigen::Upper>();
            mm = mm + Eigen::MatrixXd(mm.transpose());
            mm.diagonal() = m2;
            m = mm;
        }
        r.sn = sn;
        r.m = m;
        if (globmeth == 2) {
            // The Cholesky factor (for later use) is the same as tril(m).
            Eigen::MatrixXd L = m.triangularView<Eigen::Lower>();
            r.h = L * L.transpose();  // J'J approximation of H
        } else {
            r.h.resize(0, 0);
        }
    }
    return r;
}

// +sctool/nebroyuf.m
void nebroyuf(Eigen::MatrixXd& A, const Eigen::VectorXd& xc, const Eigen::VectorXd& xp,
              const Eigen::VectorXd& fc, const Eigen::VectorXd& fp, const Eigen::VectorXd& sx,
              double eta) {
    const Eigen::VectorXd s = xp - xc;
    const double denom = (sx.array() * s.array()).matrix().squaredNorm();
    Eigen::VectorXd tempi = fp - fc - A * s;
    for (int i = 0; i < tempi.size(); ++i) {
        if (std::abs(tempi(i)) < eta * (std::abs(fp(i)) + std::abs(fc(i)))) tempi(i) = 0.0;
    }
    A += (tempi / denom) * (s.array() * sx.array().square()).matrix().transpose();
}

// +sctool/nestop.m
void nestop(const Eigen::VectorXd& xc, const Eigen::VectorXd& xp, const Eigen::VectorXd& F,
            double Fnorm, const Eigen::VectorXd& g, const Eigen::VectorXd& sx,
            const Eigen::VectorXd& sf, int retcode, const Eigen::VectorXd& details, int itncount,
            bool maxtaken, int& consecmax, int& termcode) {
    const int n = static_cast<int>(xc.size());
    termcode = 0;
    Eigen::VectorXd invSx = sx.cwiseInverse();

    if (retcode == 1) {
        termcode = 3;
    } else if ((sf.array() * F.array().abs()).maxCoeff() <= details(7)) {
        termcode = 1;
    } else if (((xp - xc).array().abs() / xp.cwiseAbs().cwiseMax(invSx).array()).maxCoeff() <= details(8)) {
        termcode = 2;
    } else if (itncount >= details(5)) {
        termcode = 4;
    } else if (maxtaken) {
        ++consecmax;
        if (consecmax == 5) termcode = 5;
    } else {
        consecmax = 0;
        if (details(3) != 0.0 || details(2) != 0.0) {
            if ((g.cwiseAbs().array() * xp.cwiseAbs().cwiseMax(invSx).array()).maxCoeff() /
                    std::max(Fnorm, n / 2.0) <=
                details(9)) {
                termcode = 6;
            }
        }
    }
}

struct LnsrchResult {
    int retcode;
    Eigen::VectorXd xp;
    double fp;
    Eigen::VectorXd Fp;
    bool maxtaken;
};

// +sctool/nelnsrch.m (NE form only -- umflag branch not needed here)
LnsrchResult nelnsrch(const Eigen::VectorXd& xc, double fc, const Fvec& fvec,
                       const Eigen::VectorXd& g, Eigen::VectorXd p, const Eigen::VectorXd& sx,
                       const Eigen::VectorXd& sf, const Eigen::VectorXd& details, int& nofun,
                       std::vector<int>& btrack) {
    const int n = static_cast<int>(xc.size());
    LnsrchResult r;
    r.xp = Eigen::VectorXd::Zero(n);
    r.fp = 0.0;
    bool maxtaken = false;
    int retcode = 2;
    const double alpha = 1e-4;

    double newtlen = (sx.array() * p.array()).matrix().norm();
    if (newtlen > details(10)) {
        p = p * (details(10) / newtlen);
        newtlen = details(10);
    }

    const double initslope = g.dot(p);
    const double rellength = (p.cwiseAbs().array() / xc.cwiseAbs().cwiseMax(sx.cwiseInverse()).array()).maxCoeff();
    const double minlambda = details(8) / rellength;

    double lambda = 1.0;
    double lambdaprev = 0.0, fpprev = 0.0;
    int bt = 0;
    Eigen::VectorXd xp;
    double fp = 0.0;
    Eigen::VectorXd Fp;

    while (retcode >= 2) {
        xp = xc + lambda * p;
        fp = nefn(xp, sf, fvec, nofun, Fp);
        if (fp <= fc + alpha * lambda * initslope) {
            retcode = 0;
            maxtaken = (lambda == 1.0) && (newtlen > 0.99 * details(10));
        } else if (lambda < minlambda) {
            retcode = 1;
            xp = xc;
        } else {
            double lambdatemp;
            if (lambda == 1.0) {
                ++bt;
                lambdatemp = -initslope / (2.0 * (fp - fc - initslope));
            } else {
                ++bt;
                Eigen::Matrix2d Acoef;
                Acoef << 1.0 / (lambda * lambda), -1.0 / (lambdaprev * lambdaprev),
                    -lambdaprev / (lambda * lambda), lambda / (lambdaprev * lambdaprev);
                Eigen::Vector2d rhs(fp - fc - lambda * initslope, fpprev - fc - lambdaprev * initslope);
                Eigen::Vector2d a = (1.0 / (lambda - lambdaprev)) * (Acoef * rhs);
                const double disc = a(1) * a(1) - 3.0 * a(0) * initslope;
                if (a(0) == 0.0) {
                    lambdatemp = -initslope / (2.0 * a(1));
                } else {
                    lambdatemp = (-a(1) + std::sqrt(disc)) / (3.0 * a(0));
                }
                if (lambdatemp > 0.5 * lambda) lambdatemp = 0.5 * lambda;
            }
            lambdaprev = lambda;
            fpprev = fp;
            lambda = (lambdatemp <= 0.1 * lambda) ? 0.1 * lambda : lambdatemp;
        }
    }

    if (bt < static_cast<int>(btrack.size())) {
        btrack[bt] += 1;
    } else {
        btrack.resize(bt + 1, 0);
        btrack[bt] = 1;
    }

    r.retcode = retcode;
    r.xp = xp;
    r.fp = fp;
    r.Fp = Fp;
    r.maxtaken = maxtaken;
    return r;
}

struct TrustResult {
    Eigen::VectorXd xp;
    double fp;
    Eigen::VectorXd Fp;
    bool maxtaken;
    int retcode;
    Eigen::VectorXd xpprev;
    double fpprev;
    Eigen::VectorXd Fpprev;
};

// +sctool/netrust.m
TrustResult netrust(int retcode, Eigen::VectorXd xpprev, double fpprev, Eigen::VectorXd Fpprev,
                     const Eigen::VectorXd& xc, double fc, const Fvec& fvec,
                     const Eigen::VectorXd& g, const Eigen::MatrixXd& L, const Eigen::VectorXd& s,
                     const Eigen::VectorXd& sx, const Eigen::VectorXd& sf, bool newttaken,
                     Eigen::VectorXd& details, int steptype, const Eigen::MatrixXd& H, int& nofun) {
    TrustResult r;
    bool maxtaken = false;
    const double alpha = 1e-4;
    const double steplen = (sx.array() * s.array()).matrix().norm();
    Eigen::VectorXd xp = xc + s;

    Eigen::VectorXd Fp;
    const double fp = nefn(xp, sf, fvec, nofun, Fp);

    const double deltaf = fp - fc;
    const double initslope = g.dot(s);
    if (retcode != 3) fpprev = 0.0;

    if (retcode == 3 && (fp >= fpprev || deltaf > alpha * initslope)) {
        retcode = 0;
        xp = xpprev;
        // fp/Fp reassigned from prev below
        details(6) = details(6) / 2.0;
        r.xp = xp;
        r.fp = fpprev;
        r.Fp = Fpprev;
        r.maxtaken = maxtaken;
        r.retcode = retcode;
        r.xpprev = xpprev;
        r.fpprev = fpprev;
        r.Fpprev = Fpprev;
        return r;
    }

    double fpOut = fp;
    Eigen::VectorXd FpOut = Fp;
    Eigen::VectorXd xpOut = xp;

    if (deltaf >= alpha * initslope) {
        const double rellength = (xp.cwiseAbs().cwiseMax(sx.cwiseInverse()).array()).matrix().cwiseInverse().cwiseProduct(s.cwiseAbs()).maxCoeff();
        if (rellength < details(8)) {
            retcode = 1;
            xpOut = xc;
        } else {
            retcode = 2;
            const double deltatemp = -initslope * steplen / (2.0 * (deltaf - initslope));
            if (deltatemp < 0.1 * details(6)) {
                details(6) = 0.1 * details(6);
            } else if (deltatemp > 0.5 * details(6)) {
                details(6) = 0.5 * details(6);
            } else {
                details(6) = deltatemp;
            }
        }
    } else {
        double deltafpred = initslope;
        if (steptype == 1) {
            deltafpred += 0.5 * (s.transpose() * H * s)(0, 0);
        } else {
            Eigen::VectorXd ttemp = L.transpose() * s;
            deltafpred += 0.5 * ttemp.squaredNorm();
        }
        if (retcode != 2 &&
            ((std::abs(deltafpred - deltaf) <= 0.1 * std::abs(deltaf)) || (deltaf <= initslope)) &&
            !newttaken && details(6) <= 0.99 * details(10)) {
            retcode = 3;
            xpprev = xp;
            fpprev = fp;
            Fpprev = Fp;
            details(6) = std::min(2.0 * details(6), details(10));
        } else {
            retcode = 0;
            if (steplen > 0.99 * details(10)) maxtaken = true;
            if (deltaf >= 0.1 * deltafpred) {
                details(6) = details(6) / 2.0;
            } else if (deltaf <= 0.75 * deltafpred) {
                details(6) = std::min(2.0 * details(6), details(10));
            }
        }
    }

    r.xp = xpOut;
    r.fp = fpOut;
    r.Fp = FpOut;
    r.maxtaken = maxtaken;
    r.retcode = retcode;
    r.xpprev = xpprev;
    r.fpprev = fpprev;
    r.Fpprev = Fpprev;
    return r;
}

struct HookResult {
    int retcode;
    Eigen::VectorXd xp;
    double fp;
    Eigen::VectorXd Fp;
    bool maxtaken;
};

// +sctool/nehook.m (incorporates A6.4.1 + A6.4.2)
HookResult nehook(const Eigen::VectorXd& xc, double fc, const Fvec& fvec, const Eigen::VectorXd& g,
                   const Eigen::MatrixXd& L, const Eigen::MatrixXd& H, const Eigen::VectorXd& sN,
                   const Eigen::VectorXd& sx, const Eigen::VectorXd& sf, Eigen::VectorXd& details,
                   int itn, Eigen::VectorXd& trustvars, int& nofun) {
    int retcode = 4;
    bool firsthook = true;
    const double newtlen = (sx.array() * sN.array()).matrix().norm();

    if (itn == 1 || details(6) == -1.0) {
        trustvars(0) = 0.0;
        if (details(6) == -1.0) {
            Eigen::VectorXd alphav = g.cwiseQuotient(sx);
            const double alpha = alphav.dot(alphav);
            Eigen::VectorXd betav = L.transpose() * (g.array() / (sx.array() * sx.array())).matrix();
            const double beta = betav.dot(betav);
            details(6) = std::pow(alpha, 1.5) / beta;
            if (details(6) > details(10)) details(6) = details(10);
        }
    }

    Eigen::VectorXd xpprev = Eigen::VectorXd::Zero(xc.size());
    double fpprev = 0.0;
    Eigen::VectorXd Fpprev = Eigen::VectorXd::Zero(xc.size());

    Eigen::VectorXd xp, Fp, s;
    double fp = 0.0;
    bool maxtaken = false;
    double phiprimeinit = 0.0;

    while (retcode >= 2) {
        const double hi = 1.5, lo = 0.75;
        bool newttaken;
        if (newtlen <= hi * details(6)) {
            newttaken = true;
            s = sN;
            trustvars(0) = 0.0;
            details(6) = std::min(details(6), newtlen);
        } else {
            newttaken = false;
            if (trustvars(0) > 0.0) {
                trustvars(0) -= ((trustvars(2) + trustvars(1)) / details(6)) *
                                (((trustvars(1) - details(6)) + trustvars(2)) / trustvars(3));
            }
            trustvars(2) = newtlen - details(6);
            if (firsthook) {
                firsthook = false;
                Eigen::VectorXd rhs = (sx.array().square() * sN.array()).matrix();
                Eigen::VectorXd tempvec = L.triangularView<Eigen::Lower>().solve(rhs);
                phiprimeinit = -(tempvec.dot(tempvec)) / newtlen;
            }
            double mulow = -trustvars(2) / phiprimeinit;
            double muup = (g.array() / sx.array()).matrix().norm() / details(6);
            bool done = false;
            while (!done) {
                if (trustvars(0) < mulow || trustvars(0) > muup) {
                    trustvars(0) = std::max(std::sqrt(mulow * muup), muup * 1e-3);
                }
                Eigen::MatrixXd Hp =
                    H + trustvars(0) * sx.array().square().matrix().asDiagonal().toDenseMatrix();
                Eigen::MatrixXd L642;
                double maxadd;
                nechdcmp(Hp, 0.0, L642, maxadd);
                Eigen::VectorXd s642 =
                    -(L642.transpose().triangularView<Eigen::Upper>().solve(
                        L642.triangularView<Eigen::Lower>().solve(g)));
                s = s642;
                const double steplen = (sx.array() * s.array()).matrix().norm();
                trustvars(2) = steplen - details(6);
                Eigen::VectorXd rhs2 = (sx.array().square() * s.array()).matrix();
                Eigen::VectorXd tempvec2 = L642.triangularView<Eigen::Lower>().solve(rhs2);
                trustvars(3) = -(tempvec2.dot(tempvec2)) / steplen;
                if ((steplen >= lo * details(6) && steplen <= hi * details(6)) || (muup - mulow <= 0.0)) {
                    done = true;
                } else {
                    mulow = std::max(mulow, trustvars(0) - (trustvars(2) / trustvars(3)));
                    if (trustvars(2) < 0.0) muup = trustvars(0);
                    trustvars(0) -= (steplen / details(6)) * (trustvars(2) / trustvars(3));
                }
            }
        }
        trustvars(1) = details(6);

        TrustResult tr = netrust(retcode, xpprev, fpprev, Fpprev, xc, fc, fvec, g, L, s, sx, sf,
                                  newttaken, details, 1, H, nofun);
        xp = tr.xp;
        fp = tr.fp;
        Fp = tr.Fp;
        maxtaken = tr.maxtaken;
        retcode = tr.retcode;
        xpprev = tr.xpprev;
        fpprev = tr.fpprev;
        Fpprev = tr.Fpprev;
    }

    HookResult r;
    r.retcode = retcode;
    r.xp = xp;
    r.fp = fp;
    r.Fp = Fp;
    r.maxtaken = maxtaken;
    return r;
}

}  // namespace

NesolveResult nesolve(const Fvec& fvec, const Eigen::VectorXd& x0in, Eigen::VectorXd details,
                      bool identityInitialJacobian) {
    if (details.size() < 16) {
        Eigen::VectorXd padded = Eigen::VectorXd::Zero(16);
        padded.head(details.size()) = details;
        details = padded;
    }
    details(14) = 0.0;  // no fparam support (this port never accepts fparam)
    details(3) = 0.0;   // no analytic jacobian support (this port never accepts jac)
    if (details(15) == 2.0) details(15) = 1.0;  // no scale matrix support

    const Eigen::VectorXd& x0 = x0in;
    int nofun = 0;
    std::vector<int> btrack;
    Eigen::VectorXd trustvars = Eigen::VectorXd::Zero(4);

    Eigen::VectorXd FVplus = Eigen::VectorXd::Zero(x0.size());

    Eigen::VectorXd Sx, SF;
    int termcode = 0;
    neinck(x0, FVplus, details, Sx, SF, termcode);

    NesolveResult result;
    if (termcode < 0) {
        result.xf = x0;
        result.termcode = termcode;
        return result;
    }

    int itncount = 0;

    double fc = nefn(x0, SF, fvec, nofun, FVplus);

    int consecmax = 0;
    if ((SF.array() * FVplus.array().abs()).maxCoeff() <= 1e-2 * details(7)) {
        termcode = 1;
    } else {
        termcode = 0;
    }

    Eigen::VectorXd xc = x0;
    Eigen::MatrixXd Jc;
    Eigen::VectorXd gc;
    Eigen::VectorXd FVc;

    if (termcode > 0) {
        result.xf = x0;
        result.termcode = termcode;
        return result;
    } else {
        Jc = identityInitialJacobian ? Eigen::MatrixXd::Identity(x0.size(), x0.size())
                                      : nefdjac(fvec, FVplus, x0, Sx, details, nofun);
        gc = Jc.transpose() * (FVplus.array() * SF.array().square()).matrix();
        FVc = FVplus;
    }

    bool restart = true;
    double xplus_fp = 0.0;
    Eigen::VectorXd xplus, Fplus;

    while (termcode == 0) {
        ++itncount;

        if (details(3) != 0.0 || details(2) != 0.0 || (1.0 - details(4)) != 0.0) {
            ModelResult mr = nemodel(FVc, Jc, gc, SF, Sx, static_cast<int>(details(1)));
            Eigen::MatrixXd M = mr.m;
            Eigen::MatrixXd Hc = mr.h;
            Eigen::VectorXd sN = mr.sn;

            int retcode;
            double fplus;
            bool maxtaken;
            if (details(1) == 1.0) {
                LnsrchResult lr = nelnsrch(xc, fc, fvec, gc, sN, Sx, SF, details, nofun, btrack);
                retcode = lr.retcode;
                xplus = lr.xp;
                fplus = lr.fp;
                Fplus = lr.Fp;
                maxtaken = lr.maxtaken;
            } else if (details(1) == 2.0) {
                Eigen::MatrixXd L = M.triangularView<Eigen::Lower>();
                HookResult hr = nehook(xc, fc, fvec, gc, L, Hc, sN, Sx, SF, details, itncount,
                                        trustvars, nofun);
                retcode = hr.retcode;
                xplus = hr.xp;
                fplus = hr.fp;
                Fplus = hr.Fp;
                maxtaken = hr.maxtaken;
            } else {
                throw std::runtime_error("nesolve: Dogleg not implemented.");
            }

            if (retcode != 1 || restart || details(3) != 0.0 || details(2) != 0.0) {
                if (details(3) != 0.0) {
                    throw std::runtime_error("nesolve: analytic jacobian not supported in this port.");
                } else if (details(2) != 0.0) {
                    Jc = nefdjac(fvec, Fplus, xplus, Sx, details, nofun);
                } else if (details(4) != 0.0) {
                    throw std::runtime_error("nesolve: factored secant method not implemented.");
                } else {
                    nebroyuf(Jc, xc, xplus, FVc, Fplus, Sx, details(12));
                }
                if (details(4) != 0.0) {
                    throw std::runtime_error("nesolve: gradient calc for factored method not implemented.");
                } else {
                    gc = Jc.transpose() * (Fplus.array() * SF.array().square()).matrix();
                }
                nestop(xc, xplus, Fplus, fplus, gc, Sx, SF, retcode, details, itncount, maxtaken,
                       consecmax, termcode);
            }

            if ((retcode == 1 || termcode == 2) && !restart && details(3) == 0.0 && details(2) == 0.0) {
                Jc = nefdjac(fvec, FVc, xc, Sx, details, nofun);
                gc = Jc.transpose() * (FVc.array() * SF.array().square()).matrix();
                if (details(1) == 2.0 || details(1) == 3.0) details(6) = -1.0;
                restart = true;
                if (termcode == 2) termcode = 0;
            } else {
                if (termcode > 0) {
                    result.xf = xplus;
                } else {
                    restart = false;
                }
                xc = xplus;
                fc = fplus;
                FVc = Fplus;
            }
        } else {
            throw std::runtime_error("nesolve: Factored model not implemented.");
        }
    }

    result.xf = xc;
    result.termcode = termcode;
    return result;
}

}  // namespace sctoolbox
