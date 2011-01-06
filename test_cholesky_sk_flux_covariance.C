//
// ------------------------------------------------------------------------------------------------
// Test Cholesky decomposition of SuperK flux covariance matrix
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("libT2KDisapp2010a");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::numu10a;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  TMatrixD & sk_flux_cov = skexp->GetBeamInputs()->SuperkNuMuFluxCov();

  sk_flux_cov.Invert().Print();

  TMatrixD L = util::CholeskyDecomposition(sk_flux_cov);
  L.Print();
}

