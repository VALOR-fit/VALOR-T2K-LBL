//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits for mock-data 
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_mock_data_fit()
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

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  const double input_sinsq2theta23 = 0.97;
  const double input_dm23sq        = 2.5E-3;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->CreateMockData(*skexp, *fparm);

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs = 
      new OscillationFitter::FitInputs(skobs,skexp,fparm);

  // get fitter
  OscillationFitter * fitter = new OscillationFitter(fit_inputs,"minuit");
  fitter->ActualOscParams(input_sinsq2theta23, input_dm23sq);

  // do the fit
  fitter.DoTheFit();

  // draw results
  fitter->Draw(new TCanvas(), "Results");
}
