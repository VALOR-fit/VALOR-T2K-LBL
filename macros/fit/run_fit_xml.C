//
// ------------------------------------------------------------------------------------------------
// Run fits to experiments read from an input XML file and store best-fit parameters and spectra.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_fit_xml
(
  const char * tag,                                                    /* */
  int          run,                                                    /* */
  const char * xmlfile,                                                /* input XML file with toy MC experiments */
  int          first_experiment = 0,                                   /* */
  int          last_experiment  = 0,                                   /* */
  const char * config_file      = "analysis_configuration.xml",        /* analysis configuration file */
  const char * fit_param_list   = "nuosc3f_numudisapp_fix4+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option       = "minuit",                            /* fitter option */
  bool         store_toy_expts         = true,                         /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = true                          /* disable this if you don't want to store the best-fit spectra */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gOscDmsqScheme        = conv::kNuDmsq12andAbs23;
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

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each expt.
  TObjArray * dat_array = new TObjArray;
  TObjArray * pdf_array = new TObjArray;
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);

  // Create file to store results
  TFile f(Form("%s-run_%d-expts_%d_%d.root", tag, run, first_experiment, last_experiment), "recreate");

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params", "t2kana::FitParamList",  &fparm,     1280000,0);
//fitnt->Branch("input_params",    "t2kana::PhysParamList", &iparm,     1280000,0);
  fitnt->Branch("data",            "TObjArray",             &dat_array, 1280000,0);
  fitnt->Branch("best_fit_pdf",    "TObjArray",             &pdf_array, 1280000,0);

  // Loop over toy-MC experiments
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // Reset output spectra
     dat_array->Clear();
     pdf_array->Clear();

     // Re-seed fitter (reset current values to specified init values)
     fparm->ResetCurrToInitValues();

     // Get superk spectrum observation (toy experiment)
     SuperkObservation * skobs = new SuperkObservation;
     skobs->Read(xmlfile,iexperiment);

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

     // Do the fit
     skexp->PrintSystematicVariations();
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();
     double bestfit_chisq = fit_inputs->BestFitChisquare();

     // Print-out best-fit params
     fparm->ShowOnlyBestFitValuesInNextPrint();
     cout << "- Best-fit oscillation params (below) yield chisq = " << bestfit_chisq 
          << *fparm;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array, "pdf");
     }

     // Add all results in the tree
     fitnt->Fill();

     delete skobs;
     delete fit_inputs;
     delete fitter;

  }//iexperiment

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();
}
