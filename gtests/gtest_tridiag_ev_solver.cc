#include <gtest/gtest.h>

#include <Eigen/Sparse>
#include "tridiag_ev_solver.hh"
#include <iostream>


TEST(EvSolver, calculateCorrectEigenvalues) {
    //just some vectors that build tridiag matrix
    Eigen::VectorXd d {{0.0,4.0,3.0,2.0,1.0}};
    Eigen::VectorXd s {{0.0,0.0,5.0,6.0,7.0}};

    //wolfram alpha: evs are 12.2102, -7.7614, 6.32971, -0.77852
    Eigen::VectorXd res {{12.2102, -7.7614, 6.32971, -0.77852}};
    Eigen::VectorXd ev = tridiag_ev_solver(d,s);

    EXPECT_LE((res - ev).norm(), 1e-4);
    //std::cout << "Difference in norm= " << (res - ev).norm() << std::endl;

}
