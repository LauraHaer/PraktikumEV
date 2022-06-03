#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <iostream>
#include <math.h>  
#include <algorithm>    // std::min

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::SparseVector<double> SpVec; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

typedef Eigen::MatrixXd DeMat; // declares a column-major sparse matrix type of double
typedef Eigen::VectorXd DeVec; // declares a column-major sparse matrix type of double

struct result {
  DeMat a;
  DeMat b;
  DeMat v;
  DeVec ev;
};

// give matrix A and number of eigenvalues you want m
result lanczos(Eigen::MatrixXd A, int m)
{
  int n = A.cols(); //dimension of A
  DeMat T(n,n);
  // create matrix to store vectors v with v0=0, norm(v1)=1
  DeMat v(n, n); 
  for (int i = 0; i < n; i++)
    {
      // v0 = 0
      v(i,0) = 0;
      v(i,1) = 1;
    }

  // b1 = 0
  double b_i = 0;
  double b_ip1;
  double a_i;

  // the algorithm
  for (int i = 1; i < std::min(m,n-1); i++)
    {
      // access relevant columns of v (easier way?)
      DeVec v_i(n);
      DeVec v_im1(n);
      for (int j = 0; j < n; j++)
        {
          v_i(j) = v(j,i);
          v_im1(j) = v(j,i-1);
        }

      // line2
      DeVec w_i(n);
      w_i = A * v_i - b_i * v_im1;    //rappend(wi)
      // line3
      a_i = w_i.dot(v_i);
      // line4
      DeVec prod;
      prod = a_i * v_i;
      for (int j = 0; j < n; j++)
        { // easier way?
          w_i(j) = w_i(j) - prod(j);
        }
      // line5
      double root = sqrt(w_i.dot(w_i));
      b_ip1 = root;
      //enter values into T
      T(i-1,i-1) = a_i;
      T(i,i-1) = b_i;
      T(i-1,i) = b_i;
      // line6
      if (b_ip1 == 0) {
        result res;
        res.v = v; 
        Eigen::EigenSolver<DeMat> es(T);
        std::cout << "the evs are " << es.eigenvalues() << std::endl;
        //res.ev = es.eigenvalues();
        //res.ev = T.eigenvalues(); 
        return res;
      }    
      DeVec div = w_i/b_ip1;
      for (int j = 0; j < n; j++) {
        v(j,i+1) = div(j);    //q append(vi)
      }
    }
  std::cout << "T = " << T << std::endl;
  result res;
  res.v = v; 
  Eigen::EigenSolver<DeMat> es(T);
  std::cout << "the evs are " << es.eigenvalues() << std::endl;
  //res.ev = es.eigenvalues();
  //res.ev = T.eigenvalues(); 
  return res;
}

int 
main()
{
  int n = 5;
  // A in 5 dimensions, include function here
  //using dense for now
  Eigen::MatrixXd A{{4, -1, 0, 0, -1},
                    {-1, 4, -1, 0, 0},
                    {0, -1, 4, -1, 0},
                    {0, 0, -1, 4, 0},
                    {-1, 0, 0, 0, 4}};
                  

  std::cout << "A = " << A << std::endl;

  result res = lanczos(A,n);
  std::cout << "v = " << res.v << std::endl;
 // std::cout << "ev = " << res.ev << std::endl;
}