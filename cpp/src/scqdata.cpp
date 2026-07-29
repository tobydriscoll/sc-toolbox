#include "sctoolbox/scqdata.hpp"

#include "sctoolbox/gaussj.hpp"

namespace sctoolbox {

Eigen::MatrixXd scqdata(const Eigen::VectorXd& beta, int nqpts) {
    const int n = static_cast<int>(beta.size());
    Eigen::MatrixXd qnode = Eigen::MatrixXd::Zero(nqpts, n + 1);
    Eigen::MatrixXd qwght = Eigen::MatrixXd::Zero(nqpts, n + 1);

    for (int j = 0; j < n; ++j) {
        if (beta(j) > -1.0) {  // false for NaN too, matching MATLAB's beta>-1
            Eigen::VectorXd z, w;
            gaussj(nqpts, 0.0, beta(j), z, w);
            qnode.col(j) = z;
            qwght.col(j) = w;
        }
    }
    Eigen::VectorXd z0, w0;
    gaussj(nqpts, 0.0, 0.0, z0, w0);
    qnode.col(n) = z0;
    qwght.col(n) = w0;

    Eigen::MatrixXd qdat(nqpts, 2 * (n + 1));
    qdat << qnode, qwght;
    return qdat;
}

}  // namespace sctoolbox
