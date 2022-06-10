#ifndef LANCZOS_HH_
#define LANCZOS_HH_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iostream>

template <class DeMat>
struct result_lanczos {
  DeMat a;
  DeMat b;
  DeMat v;
  Eigen::VectorXcd ev;
};

//  give matrix A and number of eigenvalues you want m, and start vector aV1
template <class aMat, typename aData, int aN,
          template <typename, int, int> class aVec>
result_lanczos<Eigen::Matrix<aData, -1, -1>> lanczos(aMat A, const int m,
                                                     aVec<aData, aN, 1> aV1) {
  // We define a Dense Matrix of Dynamik Size used in the calculations based on
  // the Template Datatype
  typedef Eigen::Matrix<aData, -1, -1> tmpMat;
  assert(A.cols() == A.rows());
  assert(A.cols() == aV1.rows());
  assert(m <= aV1.rows());
  assert(sqrt(aV1.dot(aV1)) == 1);

  tmpMat TMat = tmpMat::Zero(m, m);
  tmpMat v = tmpMat::Zero(aV1.rows(), m + 2);
  v.col(1) = aV1;
  aData b_i = 0;

  // the algorithm
  for (int i = 1; i <= m; i++) {
    Eigen::Matrix<aData, aN, 1> w_i =
        A * v.col(i) - b_i * v.col(i - 1);  // line2
    aData a_i = w_i.dot(v.col(i));          // line3
    w_i = w_i - a_i * v.col(i);             // line4
    aData b_ip1 = sqrt(w_i.dot(w_i));       // line5

    TMat(i - 1, i - 1) = a_i;
    if (i < m) {
      TMat(i, i - 1) = b_ip1;
      TMat(i - 1, i) = b_ip1;
    }
    if (b_ip1 == 0) break;  // line 6
    v.col(i + 1) = w_i / b_ip1;
  }
  // std::cout << "TMat = " << TMat << std::endl;
  result_lanczos<tmpMat> res;
  res.v = v;
  Eigen::EigenSolver<tmpMat> es(TMat);
  // std::cout << "the evs are " << es.eigenvalues() << std::endl;
  res.ev = es.eigenvalues();
  return res;
}
#endif
