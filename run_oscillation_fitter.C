//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance oscillation fit
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

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitShapeOnlyNoSyst();
//PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read("inputs/analysis_configuration.xml");

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs = new OscillationFitter::FitInputs(skobs,skexp,fparm);

  // get fitter
  OscillationFitter * fitter = new OscillationFitter(fit_inputs,"minuit");

  fitter->ActualOscParams(0.97,2.4E-3);

  // do the fit
  fitter.DoTheFit();

  // see inputs
  fitter->Draw(new TCanvas(), "Results");
}
