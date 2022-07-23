#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ctime>
#include <cstdlib>

#include "gtest_helpfunctions.hh"
#include "Lanczos.hh"
#include "helpfunctions.hh"
#include "standard_include.hh"
#include "tridiag_ev_solver.hh"

TEST(GIVENS, CalculateQFromSmallDense) {
    int n = 5;
    Eigen::MatrixXd A{{1.6,8.7,0,0,0},
                        {5.5,8.8,5.4,0,0},
                        {0.0,2.8,1.5,9.4,0},
                        {0,0,2,9,10.2},
                        {0,0,0,3.7,2.1}};

    Eigen::MatrixXd q = givens_q(A);
    //test against exact q from internet
    Eigen::MatrixXd Q_exact{{-0.3,0.9,0.3,0,0.3},
                {-0.1,-0.3,-0.1,0,-0.1},
                {0,0.4,-0.6,0.1,-0.6},
                {0,0,-0.7,-0.1,0.7},
                {0,0,0,-1,-0.1}};
    Eigen::MatrixXd DiffQ = q - Q_exact;

    //test if q is orthogonal
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);
    Eigen::MatrixXd DiffId = (q.transpose() * q) - I;
    std::cout << "givens q = " << q << std::endl;

    //test if R is upper triangular
    Eigen::MatrixXd R = q.transpose() * A;
    std::cout << "givens r = " << R << std::endl;

    double e = 1;
    double e_largest = 1;
    for (int i = 1; i< n; i++) {
        e = std::max(R.diagonal(-i).maxCoeff(), - R.diagonal(-i).minCoeff());
        if (e_largest < e) {
            e_largest = e;
        }
    }
  
    EXPECT_LE(e, 1e-8 );
    //EXPECT_LE((DiffQ).norm(), 1e-8 );
    EXPECT_LE((DiffId).norm(), 1e-8 );
}

TEST(GIVENS, CalculateQFromLargeRandom) {
    std::srand(std::time(nullptr));
    int n = std::rand() % 400 + 100;
    Eigen::MatrixXd A = CreateTridiagMatrix(n);

    Eigen::MatrixXd q = givens_q(A);

    //test if q is orthogonal
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);
    Eigen::MatrixXd DiffId = (q.transpose() * q) - I;

    //std::cout << "givens q = " << q << std::endl;

    //test if R is upper triangular
    Eigen::MatrixXd R = q.transpose() * A;
    //std::cout << "givens r = " << R << std::endl;

    double e = 1;
    double e_largest = 1;
    for (int i = 1; i< n; i++) {
        e = std::max(R.diagonal(-i).maxCoeff(), - R.diagonal(-i).minCoeff());
        if (e_largest < e) {
            e_largest = e;
        }
    }
  
    EXPECT_LE(e, 1e-8 );
  
    EXPECT_LE((DiffId).norm(), 1e-8 );
}