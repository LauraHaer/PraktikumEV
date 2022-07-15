#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <ctime>

//#include "Lanczos.hh"
#include "helpfunctions.hh"
#include "standard_include.hh"
#include "tridiag_ev_solver.hh"


int main () {
    std::cout << "Random Full Dense:" << ", ";

    std::srand(std::time(nullptr));
    int m = std::rand() % 5 + 5;
    int n = m * 2;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    A = A + A.transpose().eval();
    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    auto eigenvalues = es.eigenvalues();
    auto eigenvectors = es.eigenvectors();

    Eigen::MatrixXd vecs(n,n);
    for (int i = 0; i < n; i++) {
        vecs.col(i) = inverseIteration(A, eigenvalues(i).real());
    }
    mos << "Eigen calculated eigenvectors" << std::endl;
    std::cout << eigenvectors << ", ";
    std::cout << std::endl;
    mos << "inverse iteration eigenvectors" << std::endl;
    std::cout << vecs << ", ";
    std::cout << std::endl;

    /*std::cout << "Dense Diagonal:" << ", ";

    std::srand(std::time(nullptr));
    m = std::rand() % 5 + 5;
    n = m * 2;
    A = Eigen::MatrixXd::Identity(n,n);
    A.diagonal().setRandom();
    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    eigenvalues = es.eigenvalues();;
    eigenvectors = es.eigenvectors();

    Eigen::MatrixXd vecs(eigenvalues.cols(),eigenvalues.cols());
    for (int i = 0; i < eigenvalues.cols(); i++) {
        vecs.col(i) = inverseIteration(A, eigenvalues(i).real());
    }
    mos << "Eigen calculated eigenvectors" << std::endl;
    std::cout << eigenvectors << ", ";
    std::cout << std::endl;
    mos << "inverse iteration eigenvectors" << std::endl;
    std::cout << vecs << ", ";
    std::cout << std::endl;*/

    return 0;

}
