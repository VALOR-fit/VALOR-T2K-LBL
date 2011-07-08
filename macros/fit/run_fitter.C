//
// ------------------------------------------------------------------------------------------------
// Run fitter on data
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_fitter
(
  int          dataset_id   = 0,                            /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file    = "inputs/data.xml",            /* input SK data */
  const char * config_file  = "analysis_configuration.xml", /* analysis configuration file */
  const char * hypothesis   = "numudisapp2f",               /* physics hypothesis ("numudisapp2f", "decay", "decoh")  */
  const char * fitparamlist = "numudisapp2f",               /* fit parameter list ("numudisapp2f", "numudisapp2f+nuisance", "decay", "decoh")  */
  const char * fitter_name  = "minuit",                     /* fitter name ("minuit", "custom" (2-D grid scan) */
  bool         plot         = false
)
{
  using namespace std;
  using namespace t2kana;

  cout << "Fit configration:" << endl;
  cout << "* Hypothesis...........: " << hypothesis   << endl;
  cout << "* Fit parameter list...: " << fitparamlist << endl;
  cout << "* Fitter...............: " << fitter_name  << endl;

  // get superk spectrum predictor & configure

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // get superk spectrum observation

  const char * data_file_fullpath = Form("%s/%s", gSystem->Getenv("T2KOSCANA"), data_file);
  cout << "Reading dataset: " << dataset_id << " from: " << data_file_fullpath << endl;

  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);

  // get specified fit parameter list & build fit inputs
  FitParamList * fparm = FitParamList::NamedParameterList(fitparamlist);
  FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

  // get specified fitter
  FitterABC * fitter = 0;
  //... 2-flavour oscillation / numu disappearance
  if(hypothesis == "numudisapp2f") { fitter = new NuMuDisappFitter(fit_inputs,fitter_name); }
  // ... neutrino decay
  if(hypothesis == "decay")        { fitter = new NuMuDecayFitter (fit_inputs,fitter_name); }
  // ... neutrino decoherence
  if(hypothesis == "decoh")        { fitter = new NuMuDecohFitter (fit_inputs,fitter_name); }
  assert(fitter);

  // do the fit

  TStopwatch w;
  w.Start();

  fitter->DoTheFit();

  w.Stop();
  cout << "Fit performed in " << w.CpuTime() << " sec." << endl;

  // draw/save outputs
  if(plot) {
    const char * outfile = Form("fit_%d",dataset_id);
    fitter->Draw(new TCanvas(), "StandardOscFitPlots",outfile);
  //fitter->Draw(new TCanvas(), "StandardOscFitPlots+Minos90CL",outfile);
    fitter->Draw(new TCanvas(), "ChisqSurface");
    fitter->Draw(new TCanvas(), "MarginalizedChisq");
  }


}
