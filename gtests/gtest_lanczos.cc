/*
TEST(LANCZOS, DISABLED_TestRandomMatrixWithEigenvector) {
  for(int i = 0; i < 10; ++i) {
//    int n = rand() % 15 + 5;
//    int m = rand() % 16 + 10;
    int n = rand() % 5 + 5;
    int m = rand() % 6 + 10;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n*n, n*n);
    A = A + A.transpose().eval();
    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    Eigen::VectorXd v1 = es.eigenvectors().col(0).real();
    result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-17);
    EXPECT_LE(std::abs(res.ev[0].real() - es.eigenvalues()[0].real()), 1e-8 );
  }
}
*/
