//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are read from the input XML file.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_n_mock_data_fits_xml
(
  const char * tag,                                                /* */
  int          run,                                                /* */
  const char * xmlfile,                                            /* input XML file with toy MC experiments */
  int          first_experiment,                                   /* */
  int          last_experiment,                                    /* */
  const char * config_file         = "analysis_configuration.xml", /* analysis configuration file */
  const char * fit_param_list      = "numudisapp2f",               /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option          = "minuit"                      /* fitter option (see run_fitter.C) */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiments from: " << xmlfile << endl; 
  cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl; 
  cout << "- Fit" << endl;
  cout << "   - Fit param list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit param list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Store oscillation fit results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","", 
     "iexp:mixing12:mixing13:mixing23:deltacp:dm12sq:dm23sq:chisq:fitstatus");

  // Loop over toy-MC experiments
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // Re-seed fitter (reset current values to specified init values)
     fparm->Reset();

     // Get superk spectrum observation (toy experiment)
     SuperkObservation * skobs = new SuperkObservation;
     skobs->Read(xmlfile,iexperiment);

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

     // Do the fit
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();

     // Get best-fit chisquare
     bestfit_chisq = fit_inputs->BestFitChisquare();

     // Get best-fit params
     double bestfit_mixing12 = -1;
     double bestfit_mixing13 = -1;
     double bestfit_mixing23 = -1;
     double bestfit_deltacp  = -1;
     double bestfit_dm12sq   = -1;
     double bestfit_dm23sq   = -1;
     double bestfit_chisq    = -1;
     if(hyp == conv::kHypOsc3F) {
        bestfit_mixing12 = fparm->BestFitValue(t2kana::conv::kPhypMixing12);
        bestfit_mixing13 = fparm->BestFitValue(t2kana::conv::kPhypMixing13);
        bestfit_mixing23 = fparm->BestFitValue(t2kana::conv::kPhypMixing23);
        bestfit_deltacp  = fparm->BestFitValue(t2kana::conv::kPhypDeltaCP);
        bestfit_dm12sq   = fparm->BestFitValue(t2kana::conv::kPhypDm12sq);
        bestfit_dm23sq   = fparm->BestFitValue(t2kana::conv::kPhypDm23sq);
     } 
     else
     if(hyp == conv::kHypNuMuDisapp2F) {
        bestfit_mixing23 = fparm->BestFitValue(t2kana::conv::kPhypMixing23);
        bestfit_dm23sq   = fparm->BestFitValue(t2kana::conv::kPhypDm23sq);
     }

     // Print-out best-fit params
     fparm->ShowOnlyBestFitValuesInNextPrint();
     cout << "- Best-fit oscillation params (below) yield chisq = " << bestfit_chisq 
          << *fparm;

     // Store best-fit values
     nfitnt->Fill(
           iexperiment, 
           bestfit_mixing12, bestfit_mixing13, bestfit_mixing23,
           bestfit_deltacp, bestfit_dm12sq, bestfit_dm23sq, 
           bestfit_chisq, fitstatus);

     delete skobs;
     delete fit_inputs;
     delete fitter;
  }

  //
  // Store results
  //
  TFile f(Form("%s-run_%d-expts_%d_%d.root",
               tag, run, first_experiment, last_experiment), "recreate");
  nfitnt->Write("nfitnt");
  f.Close();
}
