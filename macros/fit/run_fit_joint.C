//
// ------------------------------------------------------------------------------------------------
// Run fits to toy MC experiments generated on-the-fly and store the toy MC, the best-fit parameter
// values and the best-fit spectra, and the true input parameter values.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_fit_joint
(
  const char * tag                     = "toyexptfits",                       /* tag for output file */
  int          run                     = 0,                                   /* run number */
  int          nexperiments            = 1,                                   /* number of toy MC experiment fits */
  double       input_mixing12          = -99, //0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing13          = -99, //0.02146, //0.0251,                                /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing23          = -99, //0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
  double       input_deltacp           = -99, //0.,                                /* deltacp, if -999 then randomize x TMath::Pi()!!!!*/
  double       input_dm12sq            = -99, //7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq            = -99, //2.32E-3,                             /* |Dm^2_{23}| */

  const char * config_file             = "analysis_configuration_run1to5_test.xml",        /* analysis configuration file */
  //const char * config_file             = "analysis_configuration_run1to6.xml",        /* analysis configuration file */
  
  //const char * config_file             = "analysis_configuration_test.xml",        /* analysis configuration file */
  //const char * config_file             = "analysis_configuration_test_rhc.xml",        /* analysis configuration file */

  //const char * config_file             = "analysis_configuration_test_multiring.xml",        /* analysis configuration file */
  //const char * config_file             = "analysis_configuration_test_multiring_mrmc.xml",        /* analysis configuration file */

  bool         stat_fluctuations       = true,                                /* include s
tatistical fluctuation in toy MC expts? */
  int          tweak_systematics       = true,                                /* generate toy MC with randomized systematics? */

  //const char * fit_param_list          = "nuosc3f_12fix",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */ 
  const char * fit_param_list          = "nuosc3f_12fix+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  //const char * fit_param_list          = "nuosc3fbar_fixnonbar_12fix",
  //const char * fit_param_list          = "nuosc3fbar_fixnonbar_12fix+nuisance",

  const char * fit_option              = "minuit",                            /* fitter option */
  //const char * fit_option              = "minuit+hesse",                            /* fitter option */

  bool         store_toy_expts         = true,                                /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false,                                /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one",                               /* one==kNuOctantOne two==kNuOctantTwo */
  //bool         apply_2nu_boundary      = true                                 /* if using the 2nu flavour osc prob, apply the physical boundary? */
  double       dcp_value               = 0.,                                    /* value of deltacp (in units of PI) to assume (when deltacp is not a free parameter) */
  bool         reactor                 = true,                               /* enable this if you want to apply the reactor constraint */
  double       tolerance               = 0.01,                                 /* value for the fitter tolerance */
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  double bestfit_chisq,bestfit_chisq2;

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
  
  Opt::Instance()->SetDefaultParameterValues();

  // Check whether some osc. params need to be randomized
  bool randomize_mixing12 = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13 = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_deltacp  = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq   = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);

  //Set to defaults
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1 && input_dm23sq > -900)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1 && input_dm12sq > -900)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90 && input_deltacp > -900) input_deltacp = conv::kGlobFitDeltaCP;

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
  cout << "       > " << util::AsString(kPhypCPV    )  << " = " << ((randomize_deltacp) ? "* randomized *" : input_deltacp) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "   - Fitter tolerance: " << tolerance << endl;
  cout << "  Applying reactor constraint: " << ((reactor)  ? "YES" : "NO" ) << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing toy MC experiments:  " << ((store_toy_expts) ?        "YES " : "NO") << endl;
  cout << "     - Storing best-fit spectra: "    << ((store_bestfit_spectrum) ? "YES " : "NO") << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Fit parameter list
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = util::NamedParameterList(fit_param_list);
  FitParamList * fparm2 = util::NamedParameterList(fit_param_list);
  
  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  //make an assumption for deltacp
  if(fparm->IsFixed(conv::kPhypCPV    )) {
    fparm ->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm1->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm2->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
   }





  //Set RPA value (discrete, so fixed in the fit)
  double RPAvalue = 0.;
  //#ifdef __FHC_RUNNING__
  fparm->Update(conv::kNuisXSecRPA, RPAvalue); //init value of RPA to 0
  fparm1->Update(conv::kNuisXSecRPA, RPAvalue); //init value of RPA to 0
  fparm2->Update(conv::kNuisXSecRPA, RPAvalue);
  //#endif
  //#ifdef __RHC_RUNNING
  //fparm->Update(conv::kRHCNuisXSecRPA, RPAvalue); //init value of RPA to 0
  //fparm1->Update(conv::kRHCNuisXSecRPA, RPAvalue); //init value of RPA to 0
  //fparm2->Update(conv::kRHCNuisXSecRPA, RPAvalue);
  //#endif





  // // Get superk spectrum predictor & configure
  // SuperkExpectation * skexp = new SuperkExpectation;
  // skexp->Read(config.Data());

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

  // Create file to store results
  string oscstr = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("%s-run_%d-%s.root", tag, run, oscstr.c_str()), "recreate");

  double fittime1, fittime2;
  int fitstatus1, fitstatus2, fitstatus;
  int fitstatus_err1, fitstatus_err2, fitstatus_err;
  double chisq1, chisq2, chisq;
  
  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("best_fit_params_fit1", "valor::FitParamList",  &fparm1,     1280000,0);
  fitnt->Branch("best_fit_params_fit2", "valor::FitParamList",  &fparm2,     1280000,0);
  fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("data",                 "TObjArray",             &dat_array,  1280000,0);
  fitnt->Branch("best_fit_pdf",         "TObjArray",             &pdf_array,  1280000,0);
  fitnt->Branch("best_fit_pdf_fit1",    "TObjArray",             &pdf_array1, 1280000,0);
  fitnt->Branch("best_fit_pdf_fit2",    "TObjArray",             &pdf_array2, 1280000,0);
  fitnt->Branch("fittime1",&fittime1,"fittime1/D");
  fitnt->Branch("fittime2",&fittime2,"fittime2/D");
  fitnt->Branch("fitstatus1",&fitstatus1,"fitstatus1/I");
  fitnt->Branch("fitstatus2",&fitstatus2,"fitstatus2/I");
  fitnt->Branch("fitstatus",&fitstatus,"fitstatus/I");
  fitnt->Branch("fitstatus_err1",&fitstatus_err1,"fitstatus_err1/I");
  fitnt->Branch("fitstatus_err2",&fitstatus_err2,"fitstatus_err2/I");
  fitnt->Branch("fitstatus_err",&fitstatus_err,"fitstatus_err/I");
  fitnt->Branch("chisq1",&chisq1,"chisq1/D");
  fitnt->Branch("chisq2",&chisq2,"chisq2/D");
  fitnt->Branch("chisq",&chisq,"chisq/D");

  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  // Loop over toy-MC experiments
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // randomize specified osc params for current toy MC expt
     if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
     if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
     if (randomize_deltacp ) { input_deltacp  = util::RandomizeCPV();  }
     if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }

     // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
     if(hyp == kHypOsc3F) {
       fparm1->Update(kPhypMixing12, input_mixing12);
       fparm1->Update(kPhypMixing13, input_mixing13);
       fparm1->Update(kPhypMixing23, input_mixing23);
       fparm1->Update(kPhypCPV    ,  input_deltacp*conv::kpi);
       fparm1->Update(kPhypDm12sq,   input_dm12sq);
       fparm1->Update(kPhypDm23sq,   input_dm23sq);
       /////////////////////////////////////////
       fparm2->Update(kPhypMixing12, input_mixing12);
       fparm2->Update(kPhypMixing13, input_mixing13);
       fparm2->Update(kPhypMixing23, input_mixing23);
       fparm2->Update(kPhypCPV    ,  input_deltacp*conv::kpi);
       fparm2->Update(kPhypDm12sq,   input_dm12sq);
       fparm2->Update(kPhypDm23sq,   input_dm23sq);
       /////////////////////////////////////////
     }
     else
     if(hyp == kHypNuMuDisapp2F) {
       fparm1->Update(kPhypMixing23, input_mixing23);
       fparm1->Update(kPhypDm23sq,   input_dm23sq);
       /////////////////////////////////////////       
       fparm2->Update(kPhypMixing23, input_mixing23);
       fparm2->Update(kPhypDm23sq,   input_dm23sq);
       /////////////////////////////////////////
     }
   
     cout << *fparm1 << endl;

     // Reset all output spectra
     dat_array ->Clear();
     pdf_array ->Clear();
     pdf_array1->Clear();
     pdf_array2->Clear();

     // Reset input true params and add true oscillation params
     iparm->Clear();
     PhysParamList true_oscillation_params = fparm1->AsPhysParamList("current");
     iparm->Append(true_oscillation_params);

     // get toy MC with/without systematics
     
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
     
    //if theta23 is not fixed then initialize to 0.4 for the first fit (first octant) and 0.6 for the second one
     if(!(fparm1->IsFixed(conv::kPhypMixing23))) {     
       fparm1->SetInitValue(conv::kPhypMixing23,0.4);
       fparm1->Update(conv::kPhypMixing23,0.4);
       fparm2->SetInitValue(conv::kPhypMixing23,0.6);
       fparm2->Update(conv::kPhypMixing23,0.6);
     }
     //if reactor constraint enabled, set nuisance
     if(reactor) {
       fparm1->SetNuisance(conv::kPhypMixing13,true);
       fparm2->SetNuisance(conv::kPhypMixing13,true);
     }
     //reset systematic variations, otherwise fit is initialized at
     //best fit values of systematics from previous fit
     skexp->ResetSystematicVariations();

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

     // Set SF Initial value to 0
     if(fparm1->Exists(conv::kNuisXSecSpecFunc))
	fparm1->Update(conv::kNuisXSecSpecFunc, 0); //init value of SpecFunc is 0.5
     
     // Chisq for true (input) params
     double true_chisq = fit_inputs->Statistic();

     // Print-out true params
     fparm1->ShowOnlyCurrValuesInNextPrint();
     cout << "- True oscillation params (below) yield chisq = " << true_chisq
          << *fparm1;

     // Do the fit

     cout << endl;
     cout << endl;
     cout << "First Octant" << endl;
     cout << endl;
     cout << endl;

     fparm1->ResetCurrToInitValues(); // reset osc params
     //skexp->PrintSystematicVariations();
     fitter->SetTol(tolerance);
     fitter->DoTheFit();
     fitstatus1 = fitter->FitStatusCode("migrad_min");
     fitstatus_err1 = fitter->FitStatusCode("hesse_err");     
     fittime1 = fitter->FitTime();
    
     double bestfit_chisq  = fparm1->GetBestFitStatistic();
     chisq1 = bestfit_chisq;
     double bestfit_chisq2 = 0;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }


     //Do the fit again at the mirror point
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta &&
	(!(fparm2->IsFixed(conv::kPhypMixing23))) ) {

       cout << endl;
       cout << endl;
       cout << "Second Octant" << endl;
       cout << endl;
       cout << endl;

       skexp->ResetSystematicVariations();

       //get the correct intial values (best fit points from first fit)
       //fparm2->CopyBestFitToInitValues(*fparm1); // SF=0 not needed!!!
       fparm2->ResetCurrToInitValues(); // reset osc params
       if(fparm2->Exists(conv::kNuisXSecSpecFunc))
	 fparm2->Update(conv::kNuisXSecSpecFunc, 0); //init value of SpecFunc is 0.5

       cout << *fparm2 << endl;

       // //set the mixing23 seed to be the mirror point of the initial best fit
       // double sinsqtheta23_bestfit = fparm1->BestFitValue(conv::kPhypMixing23);
       // double sinsqtheta13_bestfit = fparm1->BestFitValue(conv::kPhypMixing13);
       // double centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       // double reseed = sinsqtheta23_bestfit;
       // if     (reseed < centre)
       // 	 reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       // else if(reseed > centre)
       // 	 reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       // fparm2->SetInitValue(conv::kPhypMixing23, reseed);

       // Get fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
       FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
       fitter2->SetTol(tolerance);
       fitter2->DoTheFit();
       fitstatus2 = fitter->FitStatusCode("migrad_min");
       fitstatus_err2 = fitter->FitStatusCode("hesse_err");     
       bestfit_chisq2 =  fparm2->GetBestFitStatistic();
       chisq2 = bestfit_chisq2;

       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
	 skexp->Write(pdf_array2, "pdf");
       }
     }

     // Print-out best-fit params, and set the "best of the two fits" branches correctly
     cout << "- Best-fit oscillation params (below) yield chisq = ";
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if(bestfit_chisq2 < bestfit_chisq) {

	 fitstatus = fitstatus2;
	 fitstatus_err = fitstatus_err2;
	 chisq = chisq2;
	 
	 cout << bestfit_chisq2 << " from the second seed fit" << endl;
	 fparm2->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2;
	 br_fparm->SetAddress(&fparm2);
	 br_bf_array->SetAddress(&pdf_array2);
       }
       else {

	 fitstatus = fitstatus1;
	 fitstatus_err = fitstatus_err1;
	 chisq = chisq1;

	 cout << bestfit_chisq << " from the first seed fit" << endl;
	 fparm1->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1;
	 br_fparm->SetAddress(&fparm1);
	 br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << bestfit_chisq << endl;

       fitstatus = fitstatus1;
       fitstatus_err = fitstatus_err1;
       chisq = chisq1;

       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1;
       br_fparm->SetAddress(&fparm1);
       br_bf_array->SetAddress(&pdf_array1);
     }

     // Add all results in the tree
     fitnt->Fill();

     // Clean-up
     delete skobs;
     delete fit_inputs;
     delete fitter;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       delete fit_inputs2;
       delete fitter2;
     }

  }

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();
}
