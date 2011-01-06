//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance oscillation fit
//
// Inputs: 
//	- dataset_id (see inputs/data.xml)
//	- fitter name ("minuit", "custom" (2-D grid scan))
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_oscillation_fitter(int dataset_id=0, const char * fitter_name="custom")
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

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read("inputs/data.xml",dataset_id);

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs = new OscillationFitter::FitInputs(skobs,skexp,fparm);

  // get fitter
  OscillationFitter * fitter = new OscillationFitter(fit_inputs,fitter_name);

  // do the fit
  fitter->DoTheFit();

  // see inputs
  fitter->Draw(new TCanvas(), "StandardOscFitPlots+Minos90CL");
  fitter->Draw(new TCanvas(), "ChisqSurface");
  fitter->Draw(new TCanvas(), "MarginalizedChisq");
}
