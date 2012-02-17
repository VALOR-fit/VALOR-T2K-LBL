//
// ------------------------------------------------------------------------------------------------
// Test Cholesky decomposition of BANFF covariance matrix
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_cholesky_decomp_banff_covariance
(
  const char * config_file = "analysis_configuration.xml", /* analysis configuration file */
)

{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  TMatrixD & C = skexp->GetBANFF()->CovMatrix();

  cout << "BANFF covariance matrix: " << endl;
  C.Invert().Print();

  cout << "Performing Cholesky decomposition of BANFF covariance matrix... " << endl;
  TMatrixD L = util::CholeskyDecomposition(C);

  cout << "Lower-triangular Cholesky matrix L: " << endl;
  L.Print();
}

