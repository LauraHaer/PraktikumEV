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

<<<<<<< Updated upstream
=======

TEST(QR_EV_SOLVER, CalculateEVFromFixedMat) {
Eigen::MatrixXd A{{1.6 , 8.7 , 0 , 0 , 0},
                    {8.7 , 8.8 , 5.4 , 0 , 0},
                    {0.0 , 5.4 , 1.5 , 9.4 , 0},
                    {0 , 0 , 9.4 , 9 , 10.2},
                    {0 , 0 , 0 , 10.2, 2.1}};

    Eigen::VectorXd diag{{1.6,8.8,1.5,9,2.1}};
    Eigen::VectorXd sdiag{{8.7,5.4,9.4,10.2}};

    
    /*Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    Eigen::VectorXd evs = tridiag_ev_solver(diag,sdiag);
    std::cout << "Eigen evs = " << es.eigenvalues() << std::endl;*/

    std::vector<double> result = runTridEVSolver(diag, sdiag, true);
}


>>>>>>> Stashed changes
TEST(QR_EV_SOLVER, CalculateEVFromSmallDim) {
    std::srand(std::time(nullptr));
    int n = std::rand() % 5 + 5;
    Eigen::MatrixXd A(n,n); 
    Eigen::VectorXd diag = Eigen::VectorXd::Random(n);
    Eigen::VectorXd sdiag = Eigen::VectorXd::Random(n-1);
    A.diagonal(1) = sdiag;
    A.diagonal(-1) = sdiag;
    A.diagonal() = diag;

<<<<<<< Updated upstream
    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    Eigen::VectorXd evs = tridiag_ev_solver(diag, sdiag);
    std::cout << "Eigen evs = " << es.eigenvalues() << std::endl;
  
    int e = 1;
    EXPECT_LE(e, 1e-8 );

}
=======
    std::vector<double> result = runTridEVSolver(diag, sdiag, true);
}

//only necessary to calculate a small number of evs, 
//so not testing too large matrices
TEST(QR_EV_SOLVER, CalculateEVFromLargeDim) {
    std::srand(std::time(nullptr));
    int n = 10;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n,n);
    
    Eigen::VectorXd diag = Eigen::VectorXd::Random(n) * 10;
    Eigen::VectorXd sdiag = Eigen::VectorXd::Random(n-1) * 10;
    A.diagonal(1) = sdiag;
    A.diagonal(-1) = sdiag;
    A.diagonal() = diag;

    std::vector<double> result = runTridEVSolver(diag, sdiag, true);
}
>>>>>>> Stashed changes
