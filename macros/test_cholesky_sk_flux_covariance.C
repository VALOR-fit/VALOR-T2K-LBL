//
// ------------------------------------------------------------------------------------------------
// Test Cholesky decomposition of SuperK flux covariance matrix
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

{
  using namespace std;
  using namespace t2kana;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  TMatrixD & sk_flux_cov = skexp->GetNuFlux()->SuperkNuMuFluxCov();

  sk_flux_cov.Invert().Print();

  TMatrixD L = util::CholeskyDecomposition(sk_flux_cov);
  L.Print();
}

