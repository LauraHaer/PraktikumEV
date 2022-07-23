#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ctime>
#include <cstdlib>

#include "Lanczos.hh"
#include "helpfunctions.hh"
#include "standard_include.hh"
#include "tridiag_ev_solver.hh"

TEST(QR_EV_SOLVER, CalculateEVFromSmallDim) {
    std::srand(std::time(nullptr));
    int n = std::rand() % 5 + 5;
    Eigen::MatrixXd A(n,n); 
    Eigen::VectorXd diag = Eigen::VectorXd::Random(n);
    Eigen::VectorXd sdiag = Eigen::VectorXd::Random(n-1);
    A.diagonal(1) = sdiag;
    A.diagonal(-1) = sdiag;
    A.diagonal() = diag;

    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    Eigen::VectorXd evs = tridiag_ev_solver(diag, sdiag);
    std::cout << "Eigen evs = " << es.eigenvalues() << std::endl;
  
    int e = 1;
    EXPECT_LE(e, 1e-8 );

}
