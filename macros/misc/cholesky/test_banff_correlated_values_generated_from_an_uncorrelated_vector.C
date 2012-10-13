#include <iomanip.h>

void test_banff_correlated_values_generated_from_an_uncorrelated_vector(const char * config_file = "analysis_configuration.xml")
{

  using namespace std;
  using namespace t2kana;
  
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
  
  TMatrixD C = skexp->GetBANFF()->CovMatrix(); 
  
  int ncols = C.GetNcols();
  int nrows = C.GetNrows();
  assert(ncols==nrows);
  
  cout << "BANFF covariance matrix: " << endl;
  C.Print();

  cout << "Performing Cholesky decomposition of BANFF covariance matrix... " << endl;
  TMatrixD L = util::CholeskyDecomposition(C);

  cout << "Lower-triangular Cholesky matrix L: " << endl;
  L.Print();

  TVectorD bestfit = skexp->GetBANFF()->BestFitValues();
  bestfit.Print();

  // try with a null uncorrelated vector - should return something equivalent to the bestfit
  TVectorD uncorrelateddeviations(ncols);
  for(int ic=0; ic < ncols; ic++) {
    uncorrelateddeviations[ic] = 0;
  }
  TVectorD variated = util::CholeskyGenerateCorrelatedParams(L, bestfit, uncorrelateddeviations);
  variated.Print();

  // now try with random gaussian vector - check for banff current values going negative
  TRandom3 gR;
  int negbanff = 0, ntries = 1000000;
  for(int i = 0; i < ntries; i++) {
    for(int ic=0; ic < ncols; ic++) {
      uncorrelateddeviations[ic] = gR.Gaus(0,1);
    }
    TVectorD uncorrelateddeviationscopy(uncorrelateddeviations);
    variated = util::CholeskyGenerateCorrelatedParams(L, bestfit, uncorrelateddeviations);
    for(int ic=0; ic < ncols; ic++) {
      if(variated[ic] < 0) {
	//variated.Print();
	//uncorrelateddeviationscopy.Print();
	negbanff++;
      }
    }
  }
  cout << negbanff << " out of " << ntries << " times there are negative results in the banff current values" << endl;

}
