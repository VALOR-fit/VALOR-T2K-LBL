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
// nuosc3f_12fix                   minuit
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
// ChangeLog:                                                                                                                              
// @ 01/02/2012, Lorena
// Added the possibility of selecting one of the two different fitters
// depending on the initial sample we define
//
void run_fitter
(
  int          dataset_id     = 0,                                                  /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file      = "inputs/data.xml",                                  /* input SK data */
  const char * config_file    = "analysis_configuration.xml",                       /* analysis configuration file */
  const char * fit_param_list = "nuosc3f_numudisapp_fix4",                          /* fit parameter list - options above  */
  const char * fit_option     = "minuit",                                           /* fitter option */
  const string sample         = "1Rmu",                                             /* options: "1Rmu", "1Re" */
  //set initial values to different from default ones
  const char * hypothesis         = "3f_osc",                                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12     = 0.8704,                                         /* value from nue appearance TN52 */
  double       input_mixing13     = 0.10,                                           /* */
  double       input_mixing23     = 1.00,                                           /* value from nue appearance TN52 */
  double       input_deltacp      = 0.00,                                           /* */
  double       input_dm12sq       = 7.6E-5,                                         /* value from nue appearance TN52 */
  double       input_dm23sq       = 2.4E-3,                                         /* value from nue appearance TN52 */
  const char * output_option  = "Observation+BestFitSpectra+ChisqSurface+Canvases"  /* what to write-out once the fit is done */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = util::SKSampleFromString(sample);

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


  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  } else
    if(ph.Hypothesis() == conv::kHypOsc3F) {
      ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
      ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
      ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
      ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
      ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
      ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
    }
  
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
  switch(s) {
  case (conv::kSKsFCFV1RingMuLike): {
    FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);
    break;
  }
  case (conv::kSKsFCFV1RingELike): {
    FitterABC * fitter = new Superk1ReNFitter(fit_inputs,fit_option);
    break;
     }
  }

  assert(fitter);
  
  // Do the fit
  fitter->DoTheFit();

  // Post-fit analysis
  fitter->Analyze(output_option);

  // Save results
  const char * outfile = Form("fit_%s_%d.root",sample.c_str(),dataset_id);
  fitter->Save(outfile);
}
