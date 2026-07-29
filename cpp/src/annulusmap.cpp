#include "sctoolbox/annulusmap.hpp"

#include <stdexcept>

namespace sctoolbox {

AnnulusMap::AnnulusMap(const Polygon& outerPolygon, const Polygon& innerPolygon) {
    if (outerPolygon.isInf()) {
        throw std::runtime_error("AnnulusMap: unbounded outer polygon ('truncate' option) not supported");
    }
    data_.M = outerPolygon.length();
    data_.N = innerPolygon.length();
    data_.Z0 = outerPolygon.vertex();
    data_.Z1 = innerPolygon.vertex();
    data_.ALFA0 = outerPolygon.angle();
    data_.ALFA1 = 2.0 - innerPolygon.angle().array();

    qwork_ = qinit(data_, kNptq);
    params_ = dscsolv(kNptq, qwork_, /*ishape=*/false, /*linearc=*/1, data_);
}

std::complex<double> AnnulusMap::eval(std::complex<double> w) const {
    int kww = 0;
    int ic = 2;
    for (int i = 0; i < params_.w0.size(); ++i) {
        if (params_.w0(i) == w) {
            kww = i + 1;
            ic = 0;
            break;
        }
    }
    if (ic == 2) {
        for (int i = 0; i < params_.w1.size(); ++i) {
            if (params_.w1(i) == w) {
                kww = i + 1;
                ic = 1;
                break;
            }
        }
    }
    return zdsc(w, kww, ic, params_.u, params_.c, params_.w0, params_.w1, params_.phi0, params_.phi1, kNptq, qwork_,
                1, data_);
}

std::complex<double> AnnulusMap::evalinv(std::complex<double> z) const {
    for (int i = 0; i < data_.Z0.size(); ++i) {
        if (data_.Z0(i) == z) throw std::runtime_error("AnnulusMap::evalinv: the point calculated is a vertex.");
    }
    for (int i = 0; i < data_.Z1.size(); ++i) {
        if (data_.Z1(i) == z) throw std::runtime_error("AnnulusMap::evalinv: the point calculated is a vertex.");
    }
    constexpr double eps = 1e-9;
    return wdsc(z, params_.u, params_.c, params_.w0, params_.w1, params_.phi0, params_.phi1, kNptq, qwork_, eps, 1,
                data_);
}

}  // namespace sctoolbox
