//
// ------------------------------------------------------------------------------------------------
// Test Cholesky decomposition of SuperK CCQE ring-counting efficiency covariance matrix
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_cholesky_decomp_sk_ccqe_rc_covariance(void)
{
  using namespace std;
  using namespace t2kana;

  TMatrixD V(3,3);

  const double sig1 = 0.017; // Ereco: < 0.4 GeV
  const double sig2 = 0.035; // Ereco: 0.4 - 1.1 GeV
  const double sig3 = 0.093; // Ereco: > 1.1 GeV

  V(0,0) = sig1 * sig1;
  V(1,1) = sig2 * sig2;
  V(2,2) = sig3 * sig3;
  V(0,1) = 0.59 * sig1 * sig2;
  V(0,2) = 0.33 * sig1 * sig3;
  V(1,2) = 0.61 * sig2 * sig3;
  V(1,0) = V(0,1);
  V(2,0) = V(0,2);
  V(2,1) = V(1,2);

  cout << "Covariance matrix: " << endl;

  V.Print();
//V.Invert().Print();

  cout << "Running Cholesky decomposition: " << endl;

  TMatrixD L = util::CholeskyDecomposition(V);
  L.Print();
}

