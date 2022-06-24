#ifndef TRIDIAG_EV_SOLVER_HH
#define TRIDIAG_EV_SOLVER_HH

//#include "Lanczos.hh"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>

typedef Eigen::VectorXd aVec;   //change into template
typedef Eigen::MatrixXd aMat;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

aVec tridiag_ev_solver(aVec diag, aVec sdiag) {
    int n = diag.rows() - 1;

    //build T with given diagonal and subdiagonal
    aMat TMat = aMat::Identity(n,n);
    SpMat tempMat(n,n);
    std::vector<T> tripletList;
    for(int i = 0; i < n; i++) {
        tripletList.push_back(T(i,i,diag(i+1)));
        if (i < n-1) {
            tripletList.push_back(T(i,i+1,sdiag(i+2)));
            tripletList.push_back(T(i+1,i,sdiag(i+2)));
        }
    }
    tempMat.setFromTriplets(tripletList.begin(), tripletList.end());
    TMat = TMat * tempMat;

    // build q, is identity matrix in beginning (faster way?)
    std::vector<T> tripletList3;    //triplets for identiy
    for(int i = 0; i < n; i++) {
        tripletList3.push_back(T(i,i,1));
    }

    //document error as maximum of subdiagonal,
    double e = 1;

    //qr iteration until T converged to diagonal matrix
    while (e > 1e-8) {
        //compute q of qr decomposition of T
        for ( int i = 0; i < n-1; ++i) {
            //get givens matrix entries c and s
            double c = 1;
            double s = 0;
            double x = TMat(i, i);
            double y = TMat(i+1,i);
            double r = std::sqrt((std::pow(x,2) + std::pow(y,2)));
            if (r > 1e-8) {
                c = x/r;
                s = y/r;
            }

            //q_i has givens matrix inside identity matrix
            SpMat q_i(n,n);
            std::vector<T> tripletList2;
            tripletList2 = tripletList3;
            tripletList2.push_back(T(i,i,c-1));
            tripletList2.push_back(T(i+1,i+1,c-1));
            tripletList2.push_back(T(i,i+1,s));
            tripletList2.push_back(T(i+1,i,-s));
            q_i.setFromTriplets(tripletList2.begin(), tripletList2.end());

            // new T = r*q
            TMat = q_i * TMat * q_i.transpose();

            //update error
            e = TMat.diagonal<-1>().maxCoeff();
        }
    }
    return TMat.diagonal();
}

#endif
