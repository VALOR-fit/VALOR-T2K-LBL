//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance oscillation analysis
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
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation skobs;
  skobs.Read("inputs/analysis_configuration.xml);

  // create fit inputs
  OscillationFitter::FitInputs fit_inputs(skobs,skexp,*fparm);

  // get fitter
  OscillationFitter fitter(fit_inputs,"minuit");

  // do the fit
  fitter.DoTheFit();
}
