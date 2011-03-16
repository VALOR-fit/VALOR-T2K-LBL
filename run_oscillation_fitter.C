//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance oscillation fit
//
// Inputs: 
//	- dataset_id (number in <data_set_%d> XML tag, see inputs/data.xml)
//      - data_file  
//	- fitter_name ("minuit", "custom" (2-D grid scan))
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_oscillation_fitter(
  int dataset_id=0, const char * data_file="inputs/data.xml", const char * fitter_name="custom")
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

  gSystem->Load("libT2KNuOsc");

  using namespace std;
  using namespace t2kana;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs = new OscillationFitter::FitInputs(skobs,skexp,fparm);

  // get fitter
  OscillationFitter * fitter = new OscillationFitter(fit_inputs,fitter_name);

  // do the fit
  fitter->DoTheFit();

  // see inputs
  const char * outfile = Form("fit%d",dataset_id);
  fitter->Draw(new TCanvas(), "StandardOscFitPlots+Minos90CL",outfile);
  fitter->Draw(new TCanvas(), "ChisqSurface");
  fitter->Draw(new TCanvas(), "MarginalizedChisq");
}
