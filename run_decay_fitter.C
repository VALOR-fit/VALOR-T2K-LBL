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
  PhysParamList * fparm = PhysParamList::NumuDecayFitNoSyst();
  cout << *fparm << endl;

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read("inputs/data.xml",0);

  // create fit inputs
  DecayFitter::FitInputs * fit_inputs = new DecayFitter::FitInputs(skobs,skexp,fparm);

  // get fitter
  //DecayFitter * fitter = new DecayFitter(fit_inputs,"minuit");
  DecayFitter * fitter = new DecayFitter(fit_inputs,"custom");

  // do the fit
  fitter.DoTheFit();

  // see inputs
  fitter->Draw(new TCanvas(), "StandardDecFitPlots");
  fitter->Draw(new TCanvas(), "ChisqSurface");
}
