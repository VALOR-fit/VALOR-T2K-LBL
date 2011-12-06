//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are generated on the fly.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_n_mock_data_fits
(
  const char * tag                     = "toyexptfits",                  /* tag for output file */
  int          run                     = 0,                              /* run number */
  int          nexperiments            = 1000,                           /* number of toy MC experiment fits */
  double       input_mixing12          = 0.86,                           /* 12-mixing, sinq2theta or sinsqtheta depending on config choice, if -999 then randomize */
  double       input_mixing13          = 0.10,                           /* 13-mixing, sinq2theta or sinsqtheta depending on config choice, if -999 then randomize */
  double       input_mixing23          = 1.00,                           /* 23-mixing, sinq2theta or sinsqtheta depending on config choice */
  double       input_deltacp           = 0.00,                           /* deltacp, if -999 then randomize */
  double       input_dm12sq            = 7.58E-5,                        /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq            = 2.32E-3,                        /* |Dm^2_{23}| */
  const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
  bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics       = true,                           /* generate toy MC with randomized systematics? */
  const char * fit_param_list          = "numudisapp2f",                 /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option              = "minuit",                       /* fitter option (see run_fitter.C) */
  bool         store_chisq_surf        = false,                          /* by default storing only best-fit values, enable this to store the chisq surfaces   */
  bool         store_toy_expts         = false,                          /* by default storing only best-fit values, enable this to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false                           /* by default storing only best-fit values, enable this to store the best-fit spectra */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Check whether some osc. params need to be randomized
  bool randomize_mixing12 = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13 = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_deltacp  = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq   = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - Randomize systematics : "      << ((tweak_systematics) ?      "YES " : "NO") << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ?      "YES " : "NO") << endl;
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << ((randomize_mixing12) ? "* randomized *" : input_mixing12) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << ((randomize_mixing13) ? "* randomized *" : input_mixing13) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << input_mixing23 << endl;
  cout << "       > " << util::AsString(kPhypDeltaCP)  << " = " << ((randomize_deltacp) ? "* randomized *" : input_deltacp) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing chisq surfaces: "      << ((store_chisq_surf) ?       "YES " : "NO") << endl;
  cout << "     - Storing toy MC experiments:  " << ((store_toy_expts) ?        "YES " : "NO") << endl;
  cout << "     - Storing best-fit spectra: "    << ((store_bestfit_spectrum) ? "YES " : "NO") << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Create file to store results
  TFile f(Form("%s-mixing23_%.2f-dm23sq_%.5feV2-run_%d.root",
          tag, input_mixing23,input_dm23sq,run), "recreate");

  // Create ntuple to store oscillation fit results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","",
     "iexp:bestfit_mixing12:bestfit_mixing13:bestfit_mixing23:bestfit_deltacp:bestfit_dm12sq:bestfit_dm23sq:chisq:fitstatus:input_mixing12:input_mixing13:input_mixing23:input_deltacp:input_dm12sq:input_dm23sq");

  // Fit parameter list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Loop over toy-MC experiments
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // randomize specified osc params for current toy MC expt
     if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
     if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
     if (randomize_deltacp ) { input_deltacp  = util::RandomizeDeltaCP();  }
     if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }

     // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
     if(hyp == kHypOsc3F) {
       fparm->Update(kPhypMixing12, input_mixing12);
       fparm->Update(kPhypMixing13, input_mixing13);
       fparm->Update(kPhypMixing23, input_mixing23);
       fparm->Update(kPhypDeltaCP,  input_deltacp);
       fparm->Update(kPhypDm12sq,   input_dm12sq);
       fparm->Update(kPhypDm23sq,   input_dm23sq);
     }
     else
     if(hyp == kHypNuMuDisapp2F) {
       fparm->Update(kPhypMixing23, input_mixing23);
       fparm->Update(kPhypDm23sq,   input_dm23sq);
     }

     //
     // get superk spectrum observation (toy experiment)
     //
     SuperkObservation * skobs = new SuperkObservation;
     if(tweak_systematics) {
        // tweak all sources of systematic err (gaussian err)
        cout << "- Randomizing all systematic params in pdf..." << endl;
        skexp->TweakAllSystematics();     
        skexp->PrintSystematicVariations();
        // generate mock data using the current PDF (incl. systematic variations) 
        // + statistical fluctuations if requested (Poisson fluctuation in each bin)
        cout << "- Generating toy-MC..." << endl;
        skobs->CreateMockData(*skexp, *fparm, stat_fluctuations);
        // reset all systematics so that the PDF used in the fitter is the nominal one
        cout << "- Reseting all systematic params in pdf..." << endl;
        skexp->ResetSystematicVariations();
        //skexp->PrintSystematicVariations();
     } else {       
        // generate mock data using the current PDF (without systematic variations) 
        // with statistical fluctuations if requested (Poisson fluctuation in each bin)
        cout << "- Generating toy-MC..." << endl;
        skobs->CreateMockData(*skexp, *fparm, stat_fluctuations); 
     }

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

     // Chisq for true (input) params
     double true_chisq = fit_inputs->Chisquare();

     // Print-out true params
     fparm->ShowOnlyCurrValuesInNextPrint();
     cout << "- True oscillation params (below) yield chisq = " << true_chisq
          << *fparm;

     // Do the fit
     fparm->ResetCurrToInitValues(); // reset osc params
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();
     bestfit_chisq = fit_inputs->BestFitChisquare();

     // Print-out best-fit params
     fparm->ShowOnlyBestFitValuesInNextPrint();
     cout << "- Best-fit oscillation params (below) yield chisq = " << bestfit_chisq 
          << *fparm;

     // Get & store best-fit params
     // Init best-fit values
     double bestfit_mixing12 = -1;
     double bestfit_mixing13 = -1;
     double bestfit_mixing23 = -1;
     double bestfit_deltacp  = -1;
     double bestfit_dm12sq   = -1;
     double bestfit_dm23sq   = -1;
     double bestfit_chisq    = -1;
     if(hyp == kHypOsc3F) {
        bestfit_mixing12 = fparm->BestFitValue(kPhypMixing12);
        bestfit_mixing13 = fparm->BestFitValue(kPhypMixing13);
        bestfit_mixing23 = fparm->BestFitValue(kPhypMixing23);
        bestfit_deltacp  = fparm->BestFitValue(kPhypDeltaCP);
        bestfit_dm12sq   = fparm->BestFitValue(kPhypDm12sq);
        bestfit_dm23sq   = fparm->BestFitValue(kPhypDm23sq);
     } else
     if(hyp == kHypNuMuDisapp2F) {
        bestfit_mixing23 = fparm->BestFitValue(kPhypMixing23);
        bestfit_dm23sq   = fparm->BestFitValue(kPhypDm23sq);
     }
     nfitnt->Fill(iexperiment, 
                  bestfit_mixing12, bestfit_mixing13, bestfit_mixing23, 
                  bestfit_deltacp, bestfit_dm12sq, bestfit_dm23sq, 
                  bestfit_chisq, fitstatus,
		  input_mixing12, input_mixing13, input_mixing23, 
                  input_deltacp, input_dm12sq, input_dm23sq);

     string output_option = "";
     if ( store_toy_expts        ) { output_option += "Observation";     }
     if ( store_bestfit_spectrum ) { output_option += "+BestFitSpectra"; }
     if ( store_chisq_surf       ) { output_option += "+ChisqSurface";   }
     fitter->Analyze(output_option.c_str());

     // Store toy-MC expt (if asked)
     if(store_toy_expts) {
        TH1D * h = (TH1D*) fitter->GetResults().FindObject("data");
        TH1D * toyexpspectrum = new TH1D(*h);
        f.cd();
        toyexpspectrum->SetName(Form("toyexp_ereco_spectrum_1Rmu_%d",iexperiment));
        toyexpspectrum->Write();
     }

     // Store best-fit spectrum (if asked)
     if(store_bestfit_spectrum) {
        TH1D * h = (TH1D*) fitter->GetResults().FindObject("bestfit_tot");
        TH1D * bestfitspectrum = new TH1D(*h);
        bestfitspectrum->SetName(Form("bestfit_ereco_spectrum_1Rmu_%d",iexperiment));
        f.cd();
        bestfitspectrum->Write();
     }
     // Store chisq surface (if asked)
     if(store_chisq_surf) {
        TH2D * h2 = (TH2D*) fitter->GetResults().FindObject("chisq_surface");
        TH2D * chisqsurf = new TH2D(*h2); 
	chisqsurf->SetName(Form("chisqsurf_%d",iexperiment)); 
	f.cd();
	chisqsurf->Write();
     }

     delete skobs;
     delete fit_inputs;
     delete fitter;
  }

  // Store ntuple with best-fit results
  f.cd();
  nfitnt->Write("nfitnt");

  f.Close();
}
