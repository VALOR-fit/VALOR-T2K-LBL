//
// ------------------------------------------------------------------------------------------------
// Run fits to toy MC experiments generated on-the-fly and store the toy MC, the best-fit parameter
// values and the best-fit spectra, and the true input parameter values.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_fit_dcp
(
  const char * tag                     = "toyexptfits",                       /* tag for output file */
  int          run                     = 0,                                   /* run number */
  int          nexperiments            = 1,                                   /* number of toy MC experiment fits */
  double       input_deltacp           = 0,                                   /* deltacp */
  double       input_mixing12          = -99,                                 /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize. If -99, then bestfit from Conventions.h */
  double       input_mixing13          = -99,                                 /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize. If -99, then bestfit from Conventions.h  */
  double       input_mixing23          = -99,                                 /* 23-mixing, sinq2theta or sinsqtheta depending on user choice. If -99, then bestfit from Conventions.h  */
  double       input_dm12sq            = -99,                                 /* Dm^2_{12}, if -999 then randomize. If -99, then bestfit from Conventions.h  */
  double       input_dm23sq            = -99,                                 /* |Dm^2_{23}|, if -999 then randomize. If -99, then bestfit from Conventions.h  */
  const char * config_file             = "analysis_configuration.xml",        /* analysis configuration file */
  bool         stat_fluctuations       = true,                                /* include statistical fluctuation in toy MC expts? */
  int          tweak_systematics       = 1,                                   /* generate toy MC with randomized systematics? */
  const char * fit_param_list          = "nuosc3f_dcpsensitivity_fix5",       /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option              = "minuit",                            /* fitter option */
  bool         store_toy_expts         = true,                                /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = true,                                /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one"                                /* one==kNuOctantOne two==kNuOctantTwo */
 )
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

  // Check whether some osc. params need to be randomized
  bool randomize_mixing23 = (TMath::Abs(input_mixing23 + 999) < kASmallNum);
  bool randomize_mixing12 = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13 = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_dm12sq   = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);
  bool randomize_dm23sq   = (TMath::Abs(input_dm23sq   + 999) < kASmallNum);

  // Check whether some osc. params need to be randomized
  bool randomize_future_mixing23 = (TMath::Abs(input_mixing23 + 9999) < kASmallNum);
  bool randomize_future_mixing12 = (TMath::Abs(input_mixing12 + 9999) < kASmallNum);
  bool randomize_future_mixing13 = (TMath::Abs(input_mixing13 + 9999) < kASmallNum);
  bool randomize_future_dm12sq   = (TMath::Abs(input_dm12sq   + 9999) < kASmallNum);
  bool randomize_future_dm23sq   = (TMath::Abs(input_dm23sq   + 9999) < kASmallNum);

  // Check whether some osc. params need to take best fit from Conventions.h
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(TMath::Abs(input_mixing13 + 99) < kASmallNum) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(TMath::Abs(input_mixing23 + 99) < kASmallNum) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(TMath::Abs(input_mixing12 + 99) < kASmallNum) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(TMath::Abs(input_mixing13 + 99) < kASmallNum) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(TMath::Abs(input_mixing23 + 99) < kASmallNum) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(TMath::Abs(input_mixing12 + 99) < kASmallNum) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(TMath::Abs(input_dm23sq + 99) < kASmallNum)  input_dm23sq = conv::kGlobFitDm23sq;
  if(TMath::Abs(input_dm12sq + 99) < kASmallNum)  input_dm12sq = conv::kGlobFitDm12sq;


  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - Randomize systematics : "      << ((tweak_systematics) ?      "YES " : "NO");
  if(tweak_systematics) cout << tweak_systematics << endl;
  else cout << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ?      "YES " : "NO") << endl;
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypDeltaCP)  << " = " << input_deltacp << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << ((randomize_mixing12) ? "* randomized *" : input_mixing12) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << ((randomize_mixing13) ? "* randomized *" : input_mixing13) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << ((randomize_mixing23) ? "* randomized *" : input_mixing23) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << ((randomize_dm23sq)  ? "* randomized *" : input_dm23sq ) << " eV^2/c^4" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing toy MC experiments:  " << ((store_toy_expts) ?        "YES " : "NO") << endl;
  cout << "     - Storing best-fit spectra: "    << ((store_bestfit_spectrum) ? "YES " : "NO") << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;

  // Fit parameter list
  FitParamList * fparm  = FitParamList::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = FitParamList::NamedParameterList(fit_param_list);
  FitParamList * fparm2 = FitParamList::NamedParameterList(fit_param_list);
  FitParamList * fparm_dcp0 = FitParamList::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1_dcp0 = FitParamList::NamedParameterList(fit_param_list);
  FitParamList * fparm2_dcp0 = FitParamList::NamedParameterList(fit_param_list);
  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each toy-MC.
  TObjArray * dat_array = new TObjArray; //only used to initialise branch
  TObjArray * pdf_array = new TObjArray; //only used to initialise branch
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);
  TObjArray * pdf_array1 = new TObjArray;
  pdf_array1->SetOwner(false);
  TObjArray * pdf_array2 = new TObjArray;
  pdf_array2->SetOwner(false);
  TObjArray * pdf_array_dcp0 = new TObjArray; //only used to initialise branch
  pdf_array_dcp0->SetOwner(false);
  TObjArray * pdf_array1_dcp0 = new TObjArray;
  pdf_array1_dcp0->SetOwner(false);
  TObjArray * pdf_array2_dcp0 = new TObjArray;
  pdf_array2_dcp0->SetOwner(false);

  // Create file to store results
  string oscstr = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("%s-run_%d-%s.root", tag, run, oscstr.c_str()), "recreate");

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",           "t2kana::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("best_fit_params_fit1",      "t2kana::FitParamList",  &fparm1,     1280000,0);
  fitnt->Branch("best_fit_params_fit2",      "t2kana::FitParamList",  &fparm2,     1280000,0);
  fitnt->Branch("best_fit_params_dcp0",      "t2kana::FitParamList",  &fparm_dcp0, 1280000,0);
  fitnt->Branch("best_fit_params_dcp0_fit1", "t2kana::FitParamList",  &fparm1_dcp0,1280000,0);
  fitnt->Branch("best_fit_params_dcp0_fit2", "t2kana::FitParamList",  &fparm2_dcp0,1280000,0);
  fitnt->Branch("input_params",              "t2kana::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("data",                      "TObjArray",             &dat_array,  1280000,0);
  fitnt->Branch("best_fit_pdf",              "TObjArray",             &pdf_array,  1280000,0);
  fitnt->Branch("best_fit_pdf_fit1",         "TObjArray",             &pdf_array1, 1280000,0);
  fitnt->Branch("best_fit_pdf_fit2",         "TObjArray",             &pdf_array2, 1280000,0);
  fitnt->Branch("best_fit_pdf_dcp0",         "TObjArray",             &pdf_array_dcp0,  1280000,0);
  fitnt->Branch("best_fit_pdf_dcp0_fit1",    "TObjArray",             &pdf_array1_dcp0, 1280000,0);
  fitnt->Branch("best_fit_pdf_dcp0_fit2",    "TObjArray",             &pdf_array2_dcp0, 1280000,0);
  TBranch * br_fparm          = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array       = fitnt->GetBranch("best_fit_pdf");
  TBranch * br_fparm_dcp0     = fitnt->GetBranch("best_fit_params_dcp0");
  TBranch * br_bf_array_dcp0  = fitnt->GetBranch("best_fit_pdf_dcp0");

  // Loop over toy-MC experiments
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // randomize specified osc params for current toy MC expt
     if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
     if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
     if (randomize_mixing23) { input_mixing23 = util::RandomizeMixing23(); }
     if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }
     if (randomize_dm23sq  ) { input_dm23sq   = util::RandomizeDm23sq();   }

     // randomize specified osc params for current toy MC expt (future error case)
     if (randomize_future_mixing12) {
       cerr << "No future error on mixing12 defined, using current best fit value" << endl;
       if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
	 input_mixing12 = conv::kGlobFitSinsq2Theta12;
       else
	 input_mixing12 = conv::kGlobFitSinsqTheta12;
     }
     if (randomize_future_mixing13) {
       cerr << "No future error on mixing13 defined, using current best fit value" << endl;
       if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
	 input_mixing13 = conv::kGlobFitSinsq2Theta13;
       else
	 input_mixing13 = conv::kGlobFitSinsqTheta13;
     }
     if (randomize_future_dm12sq) {
       cerr << "No future error on dm21sq defined, using current best fit value" << endl;
       input_mixing12 = conv::kGlobFitSinsqTheta12;
     }
     if (randomize_future_mixing23) {
       input_mixing23 = -1;
       while(input_mixing23 < 0 || input_mixing23 > 1) {
	 if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
	   double mean  = conv::kGlobFitSinsq2Theta23;
	   double sigma = 0.172;
	 }
	 else {
	   double mean  = conv::kGlobFitSinsqTheta23;
	   double sigma = 0.045;
	 }
	 input_mixing23 = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
       }//while
     }
     if (randomize_future_dm23sq) {
       input_dm23sq = -1;
       while(input_dm23sq < 0 || input_dm23sq > 1) {
	 double mean  = conv::kGlobFitDm23sq;
	 double sigma = 0.045;
	 input_dm23sq = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
       }//while
     }

     // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
     if(hyp == kHypOsc3F) {
       fparm1->SetInitValue(kPhypMixing12, input_mixing12);
       fparm1->SetInitValue(kPhypMixing13, input_mixing13);
       fparm1->SetInitValue(kPhypMixing23, input_mixing23);
       fparm1->SetInitValue(kPhypDeltaCP,  input_deltacp);
       fparm1->SetInitValue(kPhypDm12sq,   input_dm12sq);
       fparm1->SetInitValue(kPhypDm23sq,   input_dm23sq);
       fparm1_dcp0->SetInitValue(kPhypMixing12, input_mixing12);
       fparm1_dcp0->SetInitValue(kPhypMixing13, input_mixing13);
       fparm1_dcp0->SetInitValue(kPhypMixing23, input_mixing23);
       fparm1_dcp0->SetInitValue(kPhypDeltaCP,  input_deltacp);
       fparm1_dcp0->SetInitValue(kPhypDm12sq,   input_dm12sq);
       fparm1_dcp0->SetInitValue(kPhypDm23sq,   input_dm23sq);
     }
     else
     if(hyp == kHypNuMuDisapp2F) {
       fparm1->SetInitValue(kPhypMixing23, input_mixing23);
       fparm1->SetInitValue(kPhypDm23sq,   input_dm23sq);
       fparm1_dcp0->SetInitValue(kPhypMixing23, input_mixing23);
       fparm1_dcp0->SetInitValue(kPhypDm23sq,   input_dm23sq);
     }

     // Reset all output spectra
     dat_array ->Clear();
     pdf_array ->Clear();
     pdf_array1->Clear();
     pdf_array2->Clear();
     pdf_array_dcp0 ->Clear();
     pdf_array1_dcp0->Clear();
     pdf_array2_dcp0->Clear();

     // Reset input true params and add true oscillation params
     iparm->Clear();
     PhysParamList true_oscillation_params = fparm1->AsPhysParamList("current");
     iparm->Append(true_oscillation_params);

     //
     // get toy MC with/without systematics
     //
     SuperkObservation * skobs = new SuperkObservation;
     bool applied_negevents_condition = true;
     while(applied_negevents_condition) {
       if(tweak_systematics) {
	 // tweak all sources of systematic error
	 if(tweak_systematics == 1) {
	   cout << "- Randomizing all systematic params in pdf..." << endl;
	   skexp->TweakAllSystematics(); 
	 }
	 else if(tweak_systematics == 2) {
	   cout << "- Randomizing BANFF systematic params in pdf..." << endl;
	   skexp->TweakAllBANFFSystematics();
	 }
	 else if(tweak_systematics == 3) {
	   cout << "- Randomizing uncorrelated-xsec systematic params in pdf..." << endl;
	   skexp->TweakAllUncorrelatedXSecSystematics();
	 }
	 else if(tweak_systematics == 4) {
	   cout << "- Randomizing SKDetFSI systematic params in pdf..." << endl;
	   skexp->TweakAllSKDetectorFSISystematics();
	 }
	 skexp->PrintSystematicVariations();
	 // generate mock data using the current PDF (incl. systematic variations) 
	 // + statistical fluctuations if requested (Poisson fluctuation in each bin)
	 cout << "- Generating toy-MC..." << endl;
	 //note that fparm does not contain the tweaks, so use only the PhysHyp here
	 applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), stat_fluctuations);
	 if(!applied_negevents_condition) {
	   // keep track of all nuisance params set so as to add in output tree
	   PhysParamList true_nuisance_params = skexp->Tweaks();
	   iparm->Append(true_nuisance_params);
	 }
	 // reset all systematics so that the PDF used in the fitter is the nominal one
	 cout << "- Reseting all systematic params in pdf..." << endl;
	 skexp->ResetSystematicVariations();
	 //skexp->PrintSystematicVariations();
       }
       else {
	 //make sure systematics are reset
	 skexp->ResetSystematicVariations();
	 // generate mock data using the current PDF (without systematic variations) 
	 // with statistical fluctuations if requested (Poisson fluctuation in each bin)
	 cout << "- Generating toy-MC..." << endl;
	 applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), stat_fluctuations); 
       }
     }//applied_negevents_condition

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

     // Chisq for true (input) params
     if(fparm1->Exists(conv::kNuisXSecSpecFunc))
	fparm1->Update(conv::kNuisXSecSpecFunc, 0); //init value of SpecFunc is 0.5
     double true_chisq = fit_inputs->Chisquare();

     // Print-out true params
     fparm1->ShowOnlyCurrValuesInNextPrint();
     cout << "- True oscillation params (below) yield chisq = " << true_chisq
          << *fparm1;

     // Do the fit
     fparm1->ResetCurrToInitValues(); // reset osc params
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();
     double bestfit_chisq  = fit_inputs->BestFitChisquare();
     double bestfit_chisq2 = 999999;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }

     //Do the fit again at the mirror point
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta && fparm1->GetNumberFreeParameters() > 1) {
       //get the correct intial values (best fit points from first fit)
       fparm2->CopyBestFitToInitValues(*fparm1);
       //set the mixing23 seed to be the mirror point of the initial best fit
       double sinsqtheta23_bestfit = fparm1->BestFitValue(conv::kPhypMixing23);
       double sinsqtheta13_bestfit = fparm1->BestFitValue(conv::kPhypMixing13);
       double centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       double reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
	 reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
	 reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm2->SetInitValue(conv::kPhypMixing23, reseed);
       // Get fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
       FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
       fparm2->ResetCurrToInitValues(); // reset osc params
       fitter2->DoTheFit();
       bestfit_chisq2 = fit_inputs2->BestFitChisquare();
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
	 skexp->Write(pdf_array2, "pdf");
       }
     }

     //Fit again, fixing deltacp to 0. Start from nominal (not previous best fits) for other parameters
     fparm1_dcp0->SetInitValue(kPhypDeltaCP,  0);
     fparm1_dcp0->FixParameter(kPhypDeltaCP);
     // Get fitter
     FitInputs * fit_inputs_dcp0 = new FitInputs(skobs,skexp,fparm1_dcp0);
     FitterABC * fitter_dcp0 = util::GetFitter(fit_inputs_dcp0,fit_option);
     // Chisq for true (input) params
     double true_chisq_dcp0 = fit_inputs_dcp0->Chisquare();
     // Print-out true params
     fparm1_dcp0->ShowOnlyCurrValuesInNextPrint();
     cout << "- True oscillation params (with dcp=0) (below) yield chisq = " << true_chisq_dcp0
          << *fparm1_dcp0;
     // Do the fit
	 cout << *fparm1_dcp0;
     fparm1_dcp0->ResetCurrToInitValues(); // reset osc params
     double bestfit_chisq_dcp0 = 0;
     fitter_dcp0->DoTheFit();
     int fitstatus_dcp0 = fitter_dcp0->FitStatusCode();
     bestfit_chisq_dcp0  = fit_inputs_dcp0->BestFitChisquare();
     double bestfit_chisq2_dcp0 = 999999;
     // Save spectra (data and best-fit prediction)
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1_dcp0, "pdf");
     }

     //Do the fit again at the mirror point
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta && fparm1_dcp0->GetNumberFreeParameters() > 1) {
       //get the correct intial values (best fit points from first fit)
       fparm2_dcp0->CopyBestFitToInitValues(*fparm1_dcp0);
       //set the mixing23 seed to be the mirror point of the initial best fit
       double sinsqtheta23_bestfit = fparm1_dcp0->BestFitValue(conv::kPhypMixing23);
       double sinsqtheta13_bestfit = fparm1_dcp0->BestFitValue(conv::kPhypMixing13);
       double centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       double reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
	 reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
	 reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm2_dcp0->SetInitValue(conv::kPhypMixing23, reseed);
       // Get fitter
       FitInputs * fit_inputs2_dcp0 = new FitInputs(skobs,skexp,fparm2_dcp0);
       FitterABC * fitter2_dcp0 = util::GetFitter(fit_inputs2_dcp0,fit_option);
       fparm2_dcp0->ResetCurrToInitValues(); // reset osc params
       fitter2_dcp0->DoTheFit();
       bestfit_chisq2_dcp0 = fit_inputs2_dcp0->BestFitChisquare();
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
	 skexp->Write(pdf_array2_dcp0, "pdf");
       }
     }

     // Print-out best-fit params, and set the "best of the two fits" branches correctly
     cout << "- Best-fit oscillation params (below) yield chisq = ";
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if(bestfit_chisq2 < bestfit_chisq) {
	 cout << bestfit_chisq2 << " from the second seed fit" << endl;
	 fparm2->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2;
	 br_fparm->SetAddress(&fparm2);
	 br_bf_array->SetAddress(&pdf_array2);
       }
       else {
	 cout << bestfit_chisq << " from the first seed fit" << endl;
	 fparm1->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1;
	 br_fparm->SetAddress(&fparm1);
	 br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << bestfit_chisq << endl;
       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1;
       br_fparm->SetAddress(&fparm1);
       br_bf_array->SetAddress(&pdf_array1);
     }

     // Print-out best-fit params, assuming dcp=0, and set the "best of the two fits" branches correctly
     cout << "- Best-fit oscillation params with dcp=0 (below) yield chisq = ";
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if(bestfit_chisq2_dcp0 < bestfit_chisq_dcp0) {
	 cout << bestfit_chisq2_dcp0 << " from the second seed fit" << endl;
	 fparm2_dcp0->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2_dcp0;
	 br_fparm_dcp0->SetAddress(&fparm2_dcp0);
	 br_bf_array_dcp0->SetAddress(&pdf_array2_dcp0);
       }
       else {
	 cout << bestfit_chisq_dcp0 << " from the first seed fit" << endl;
	 fparm1_dcp0->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1_dcp0;
	 br_fparm_dcp0->SetAddress(&fparm1_dcp0);
	 br_bf_array_dcp0->SetAddress(&pdf_array1_dcp0);
       }
     }
     else {
       cout << bestfit_chisq_dcp0 << endl;
       fparm1_dcp0->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1_dcp0;
       br_fparm_dcp0->SetAddress(&fparm1_dcp0);
       br_bf_array_dcp0->SetAddress(&pdf_array1_dcp0);
     }

     // Add all results in the tree
     fitnt->Fill();

     // Clean-up
     delete skobs;
     delete fit_inputs;
     delete fitter;
     delete fit_inputs_dcp0;
     delete fitter_dcp0;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta && fparm1->GetNumberFreeParameters() > 1) {
       delete fit_inputs2;
       delete fitter2;
     }
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta && fparm1_dcp0->GetNumberFreeParameters() > 1) {
       delete fit_inputs2_dcp0;
       delete fitter2_dcp0;
     }

  }//iexperiment

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();
}
