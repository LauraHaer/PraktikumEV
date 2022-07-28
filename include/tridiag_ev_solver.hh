#ifndef TRIDIAG_EV_SOLVER_HH
#define TRIDIAG_EV_SOLVER_HH

//#include "Lanczos.hh"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>

#include "standard_include.hh"
#include "helpfunctions.hh"

// typedef Eigen::VectorXd aVec;  // change into template
// typedef Eigen::MatrixXd aMat;

typedef Eigen::SparseMatrix<double>
    SpMat;  // declares a column-major sparse matrix type of double

template <class aVec>
Eigen::MatrixXd createTMatrix(aVec diag, aVec sdiag) {
  int n = diag.rows() - 1;

  // build T with given diagonal and subdiagonal
  Eigen::MatrixXd TMat = Eigen::MatrixXd::Zero(n, n);
  TMat.diagonal(1) = sdiag(Eigen::seqN(1, n - 1));
  TMat.diagonal(-1) = sdiag(Eigen::seqN(1, n - 1));
  TMat.diagonal() = diag(Eigen::lastN(n));

  return TMat;
}

template <class aVec>
Eigen::MatrixXd new_createTMatrix(aVec diag, aVec sdiag) {
  int n = diag.rows();

  // build T with given diagonal and subdiagonal
  Eigen::MatrixXd TMat = Eigen::MatrixXd::Zero(n, n);
  TMat.diagonal(1) = sdiag;
  TMat.diagonal(-1) = sdiag;
  TMat.diagonal() = diag;

  return TMat;
}

template<class aMat>
Eigen::MatrixXd givens_q(aMat A) {
  int n = A.cols();
  Eigen::MatrixXd q = Eigen::MatrixXd::Identity(n, n);
  for (int i = 0; i < n - 1; ++i) {
    // get givens matrix entries c and s
    double c = 1;
    double s = 0;
    double x = A(i, i);
    double y = A(i + 1, i);
    double r = std::sqrt((std::pow(x, 2) + std::pow(y, 2)));
    if (r > 1e-8) {
      c = x / r;
      s = y / r;
    }

    // q_i has givens matrix inside identity matrix
<<<<<<< Updated upstream
    SpMat q_i(n, n);
    for(int ii=0; ii < n; ++ii) q_i.coeffRef(ii,ii) = 1;
    q_i.coeffRef(i, i) = c;
    q_i.coeffRef(i + 1, i + 1) = c;
    q_i.coeffRef(i, i + 1) = s;
    q_i.coeffRef(i + 1, i) = -s;
=======
    aMat q_i = aMat::Identity(n,n);
    q_i(i, i) = c;
    q_i(i + 1, i + 1) = c;
    q_i(i, i + 1) = s;
    q_i(i + 1, i) = -s;
>>>>>>> Stashed changes

    A = q_i * A * q_i.transpose();

    q = q_i * q;

  }
  return q.transpose();
}

template <class aMat>
void elimErrorVal (aMat& TMat, aMat& q, int h, int j) {

  double r = std::sqrt((std::pow(TMat(h,h), 2) + std::pow(TMat(j,h), 2)));
  double c = TMat(h,h)/r;
  double s = TMat(j,h)/r;

  aMat q_plus = aMat::Identity(TMat.cols(),TMat.cols());
  q_plus(h,h) = c;
  q_plus(j,j) = c;
  q_plus(h,j) = s;
  q_plus(j,h) = -s;

  TMat = q_plus * TMat * q_plus.transpose();

  q = q_plus * q;
}

template <class aVec>
aVec tridiag_ev_solver(aVec diag, aVec sdiag) {
<<<<<<< Updated upstream
  // build T with given diagonal and subdiagonal
  Eigen::MatrixXd TMat = new_createTMatrix(diag, sdiag);
  Eigen::EigenSolver<Eigen::MatrixXd> es(TMat);
  // document error as maximum of subdiagonal,
<<<<<<< HEAD
  double e = 1;

  int i = 0;
  // qr iteration until T converged to diagonal matrix
  while (e > 1e-6) {
    Eigen::MatrixXd q = givens_q(TMat);
    TMat = q.transpose() * TMat * q;
  //  // update error
    Eigen::MatrixXd test = Eigen::MatrixXd::Zero(TMat.rows(), TMat.rows());
=======
  typedef Eigen::MatrixXd aMat;
  // build T with given diagonal and subdiagonal
  Eigen::MatrixXd TMat = new_createTMatrix(diag, sdiag);
  std::cout << "Created Matrix: " << std::endl;
  mos << PRINT_REFLECTION(TMat) << std::endl;

  // document error as maximum of subdiagonal,
  double e = 1;
  int i = 0;
  // qr iteration until T converged to diagonal matrix
  while (e > 1e-3) {
    if (i>5000) {
      std::cout << " i > 5000 " << std::endl;
      break;
    }
    i++;
    aMat q = givens_q(TMat);

    if ((q-aMat::Identity(TMat.rows(),TMat.rows())).norm() < 1e-8) {
      std::cout << "q is Identity, stopped converging at i = " << i <<
                    "with error e = " << e << std::endl;
      break;
    }
    TMat = q.transpose() * TMat * q;

  // update error
    aMat test = Eigen::MatrixXd::Zero(TMat.rows(), TMat.rows());
>>>>>>> Stashed changes
    test.diagonal() = TMat.diagonal();
    test = TMat - test;
      //TMat - Eigen::MatrixXd::Identity(TMat.rows(), TMat.rows()) * TMat.diagonal().eval();
    e = std::max(test.maxCoeff(), - test.minCoeff());
<<<<<<< Updated upstream
    if (i % 100 == 0) mos << PRINT_REFLECTION(e) << PRINT_REFLECTION(TMat) << std::endl;
    ++i;
  }
  //mos << PRINT_REFLECTION(es.eigenvalues()) << std::endl;
  mos << PRINT_REFLECTION(i) << std::endl;
  mos << PRINT_REFLECTION(e) << std::endl;
  mos << PRINT_REFLECTION(TMat.diagonal()) << std::endl;
  //mos << PRINT_REFLECTION(TMat) << std::endl;
  //return es.eigenvalues().real();
=======

    //catch large values outside of the diag and subdiags
    //use givens rotation
    /*for (int k=0; k<test.rows(); k++) {
      for (int l=0; l<test.rows(); l++) {
        if(std::abs(k-l) > 2) {
          if (abs(test(k,l)) > 1e-5) {
            mos << TMat << " at " << "k = " << k << "," 
                                  << "l = " << l << "," 
                                  << "has value " << TMat(k,l) << std::endl;
            elimErrorVal(TMat,q,k,l);
            mos << TMat << " after elminination of "<< "k = " << k << "," 
                        << "l = " << l 
                        << "is now " <<  TMat(k,l) << std::endl;
            test.diagonal() = TMat.diagonal();
            test = TMat - test;
          }
        }
      }
    }*/
    //e = std::max(test.maxCoeff(), - test.minCoeff());
    //e = test.norm();
    e = ( q- aMat::Identity(TMat.rows(),TMat.rows())).norm();
     mos << PRINT_REFLECTION(e) << std::endl;
    if (i % 1000 == 0) {
      if (i < 100 && i%10 == 0) {
        mos << PRINT_REFLECTION(i) << PRINT_REFLECTION(e) << std::endl;
        mos << PRINT_REFLECTION(TMat) << std::endl;
        mos << PRINT_REFLECTION(q) << std::endl;
      }
      mos << PRINT_REFLECTION(i) << PRINT_REFLECTION(e) << std::endl;
      mos << PRINT_REFLECTION(TMat) << std::endl;
    }
  }
   if (i<5000) {
      std::cout <<"Converged after " << i << " loops." << std::endl;
    }
>>>>>>> Stashed changes
  return TMat.diagonal();
=======
//  double e = 1;
//
//  int i = 0;
//  // qr iteration until T converged to diagonal matrix
//  while (e > 1e-6) {
//    Eigen::MatrixXd q = givens_q(TMat);
//    TMat = q.transpose() * TMat * q;
//    TMat = tridiagonalize(TMat);
//  //  // update error
//    Eigen::MatrixXd test = Eigen::MatrixXd::Zero(TMat.rows(), TMat.rows());
//   test.diagonal() = TMat.diagonal();
//    test = TMat - test;
//      //TMat - Eigen::MatrixXd::Identity(TMat.rows(), TMat.rows()) * TMat.diagonal().eval();
//    e = std::max(test.maxCoeff(), - test.minCoeff());
//    //if (i % 100 == 0) mos << PRINT_REFLECTION(e) << PRINT_REFLECTION(TMat) << std::endl;
//    ++i;
//  }
//  //mos << PRINT_REFLECTION(es.eigenvalues()) << std::endl;
//  mos << PRINT_REFLECTION(i) << std::endl;
//  mos << PRINT_REFLECTION(e) << std::endl;
//  mos << PRINT_REFLECTION(TMat.diagonal()) << std::endl;
  //mos << PRINT_REFLECTION(TMat) << std::endl;
  aVec res = es.eigenvalues().real();
  return res;
  //return TMat.diagonal();
>>>>>>> 2866b7b502d2583d0f300a04fc513fec85e2c2a4
}


template <class aMat>
Eigen::VectorXd inverseIteration(const aMat T, const double mu, const double eps = 1e-8) {
  Eigen::VectorXd x_prev = Eigen::VectorXd::Zero(T.cols());
  double e1 = 1;
  double e2 = 1;

  Eigen::VectorXd x = Eigen::VectorXd::Random(T.cols());
  while( e1 > eps && e2 > eps ) {
    x = (T - (mu * Eigen::MatrixXd::Identity(T.cols(), T.cols()))).inverse() * x;
    x = x.normalized();

    e1 = (x - x_prev).norm();
    e2 = (x - (x_prev * (-1))).norm();

    x_prev = x;
  }
  return x;
}

//template <class aMat, class aVec>
//Eigen::MatrixXd eigenvectorsA(aMat T, aVec ev, aMat v) {
//  v = v(Eigen::all, Eigen::lastN(ev.cols()));
//  Eigen::MatrixXd vec(ev.cols(), ev.cols());
//
//  for (int i = 0; i < res.ev.cols(); i++) {
////    res.vec.col(i) = inverseIteration(TMat, res.ev(i).real());
//    //compute eigenvectors of T
//    vec.col(i) = inverseIteration(t, ev(i));
//    //compute eigenvectors of A from eigenvectors of T
//    vec.col(i) = v*vec.col(i);
//  }
//
//  for
//}

#endif
