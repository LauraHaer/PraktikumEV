#include <iostream>

#include <Eigen/Dense>

#include "Lanczos.hh"
#include "standard_include.hh"

int main() {
  int n = 5;
  // A in 5 dimensions, include function here
  // using dense for now
  Eigen::MatrixXd A{{4, -1, 0, 0, -1},
                    {-1, 4, -1, 0, 0},
                    {0, -1, 4, -1, 0},
                    {0, 0, -1, 4, 0},
                    {-1, 0, 0, 0, 4}};

  //std::cout << "A = " << A << std::endl;

  DeVec v1 = DeVec::Zero(n);
  v1(0) = 1;
  v1 << -1, 0, 1, -1, 1;
  v1 = v1 / sqrt(v1.dot(v1));

  //mos << PRINT_REFLECTION(A*v1) << std::endl;

  result res = lanczos(A, n, v1);

//  Eigen::MatrixXd B = Eigen::MatrixXd::Random(101,101);
//  Eigen::MatrixXd C = B + B.transpose();
//  Eigen::EigenSolver<DeMat> es(C);
//  //mos << PRINT_REFLECTION(C) << std::endl;
//  res = lanczos(C, 10);
//
//  B = Eigen::MatrixXd::Random(101,101);
//  C = B + B.transpose();
//  res = lanczos(C, 10);

  //mos << PRINT_REFLECTION(es.eigenvalues()) << std::endl;;
  //std::cout << "v = " << res.v << std::endl;
  //std::cout << "ev = " << res.ev << std::endl;
}
