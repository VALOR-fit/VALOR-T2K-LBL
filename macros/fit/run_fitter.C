//
// ------------------------------------------------------------------------------------------------
// Run fitter on data
//
// The following fit parameter list / fitter options may be used:
// ..................................................................
// Fit parameter list              Fitter option
// (see FitParamList.h)
// ..................................................................
// nuosc3f                         minuit
// nuosc3f_12constrain             minuit
// nuosc3f_numudisapp_fix4         minuit
// nuosc3f_nueapp_fix4             minuit
// numudisapp2f                    minuit
// numudisapp2f+nuisance           minuit
// decay                           minuit
// decoh                           minuit
// ..................................................................
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_fitter
(
  int          dataset_id     = 0,                            /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file      = "inputs/data.xml",            /* input SK data */
  const char * config_file    = "analysis_configuration.xml", /* analysis configuration file */
  const char * fit_param_list = "numudisapp2f",               /* fit parameter list  */
  const char * fit_option     = "minuit",                     /* fitter option */
  const char * output_option  = "Observation+BestFitSpectra+ChisqSurface+Canvases"  /* what to write-out once the fit is done */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Data file
  const char * data_file_fullpath = Form("%s/%s", gSystem->Getenv("T2KOSCANA"), data_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Print settings
  cout << "[Job inputs]" << endl;
  cout << "  Reading analysis configuration from: " << config << endl;
  cout << "  Reading dataset: " << dataset_id << " from: " << data_file_fullpath << endl;
  cout << "  Fit parameter list: " << fit_param_list << endl;
  cout << "  Fitter option: " << fit_option << endl;
  cout << "  Output option: " << output_option  << endl;
  cout << *Opt::Instance() << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);

  // Get specified fit parameter list & build fit inputs
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  assert(fparm);
  FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

  // Get specified fitter
  FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);
  assert(fitter);

  // Do the fit
  fitter->DoTheFit();

  // Post-fit analysis
  fitter->Analyze(output_option);

  // Save results
  const char * outfile = Form("fit_%d.root",dataset_id);
  fitter->Save(outfile);
}
