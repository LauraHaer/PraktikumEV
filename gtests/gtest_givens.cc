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

TEST(GIVENS, CalculateQFromSmallDense) {
    int n = 5;
    Eigen::MatrixXd A{{1.6,8.7,3.6,1.5,2.4},
                        {5.5,8.8,5.4,3.5,10},
                        {6,2.8,1.5,9.4,8.7},
                        {5.9,7.9,2,9,10.2},
                        {6,4.5,5.5,3.7,2.1}};

    q = givens_q(A);
    Q_exact = {{-0.1,0.8,0.1,0},
                {-0.5,0.3,0.3,-05,0.6},
                {-0.5,-0.4,-0.3,0.5,0.5},
                {-0.5,0.2,-0.6,-0.3,-0.5},
                {-0.5,-0.2,0.7,0.1,-0.5}}
    DiffQ = q - Q_exact;
    Eigen::MatrixXd I = Eigen::MatrixX::Identity(n,n);
    DiffId = (q.transpose() * q) - I;
  
    EXPECT_LE((DiffQ).norm(), 1e-8 );
    EXPECT_LE((DiffId).norm(), 1e-8 );
}

TEST(GIVENS, CalculateQFromLargeRandom) {
    std::srand(std::time(nullptr));
    int n = std::rand() % 400 + 100;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    A = A + A.transpose().eval();
    Eigen::MatrixXd I = Eigen::MatrixX::Identity(n,n);
    DiffId = (q.transpose() * q) - I;
  
    EXPECT_LE((DiffId).norm(), 1e-8 );
}