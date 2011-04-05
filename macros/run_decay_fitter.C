//
// ------------------------------------------------------------------------------------------------
// Run neutrino decay fit
//
// Inputs:
//      - dataset_id (see inputs/data.xml)
//      - fitter name ("minuit", "custom" (2-D grid scan))
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_decay_fitter(int dataset_id=0, const char * fitter_name="custom")
{
  using namespace std;
  using namespace t2kana;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDecayFitNoSyst();
  cout << *fparm << endl;

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read("inputs/data.xml",dataset_id);

  // create fit inputs
  FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

  // get fitter
  NuMuDecayFitter * fitter = new NuMuDecayFitter(fit_inputs,fitter_name);

  // do the fit
  fitter->DoTheFit();

  // see inputs
  fitter->Draw(new TCanvas(), "StandardDecFitPlots");
  fitter->Draw(new TCanvas(), "ChisqSurface");
}
