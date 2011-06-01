//
// ------------------------------------------------------------------------------------------------
// Test Cholesky decomposition of SuperK flux covariance matrix
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_cholesky_decomp_sk_flux_covariance
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

  TMatrixD & sk_flux_cov = skexp->GetNuFlux()->SuperkNuMuFluxCov();

  sk_flux_cov.Invert().Print();

  TMatrixD L = util::CholeskyDecomposition(sk_flux_cov);
  L.Print();
}

