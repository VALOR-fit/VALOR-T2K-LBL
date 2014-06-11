//
// ------------------------------------------------------------------------------------------------
// Run fits to toy MC experiments generated on-the-fly and store the toy MC, the best-fit parameter
// values and the best-fit spectra, and the true input parameter values.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_sensitivity_joint_toys
(
  const char * tag                     = "toyexptsensi",                       /* tag for output file */
  int          run                     = 0,                                   /* run number */
  int          nexperiments            = 3,                                   /* number of toy MC experiment fits */
  double       input_mixing12          = 0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing13          = 0.0251,                                /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing23          = 0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
  double       input_deltacp           = 0.00,                                /* deltacp, if -999 then randomize */
  double       input_dm12sq            = 7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq            = 2.32E-3,                             /* |Dm^2_{23}| */
  const char * config_file             = "analysis_configuration_test.xml",        /* analysis configuration file */
  bool         stat_fluctuations       = true,                                /* include statistical fluctuation in toy MC expts? */
  int          tweak_systematics       = 1,                                /* generate toy MC with randomized systematics? */
  //const char * fit_param_list        = "nuosc3f_numudisapp_nuis3fixdcp+nuisance",  /* fit parameter list for the given physics hypothesis (see UtlFitParamListFactory.h) */
  //const char * fit_param_list          = "nuosc3f_12constrain+nuisance",  /* fit parameter list for the given physics hypothesis (see UtlFitParamListFactory.h) */
  const char * fit_param_list          = "nuosc3f_12constrain",  /* fit parameter list for the given physics hypothesis (see UtlFitParamListFactory.h) */
  const char * fit_option              = "minuit+hesse",                            /* fitter option */
  bool         store_toy_expts         = true,                                /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false,                                /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one",                               /* one==kNuOctantOne two==kNuOctantTwo */
  //bool         apply_2nu_boundary      = true                                 /* if using the 2nu flavour osc prob, apply the physical boundary? */
  double       dcp_value               = 0,                                    /* value of deltacp (in units of PI) to assume (when deltacp is not a free parameter) */

  //

  bool         reactor                 = true,                               /* enable this if you want to apply the reactor constraint */
  double       tolerance               = 0.01,                                 /* value for the fitter tolerance */
  const char * param1                  = "deltacp",
  const char * param2                  = "none",                  /* second parameter to be fitting vs (for 2D) - default: none for a single fit or 1D */ 
  int          min_bin_param1          = 1 ,                                /* minimum bin for param1 */
  int          max_bin_param1          = 51 ,                                /* maximum bin for param1 */
  int          min_bin_param2          = 1 ,                                /* minimum bin for param2 */
  int          max_bin_param2          = 1                                  /* maximum bin for param2 */
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

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
  cout << "- Fit vs " << param1 << " and " << param2 << endl; 
  cout << "- Applying reactor constraint: " << ((reactor)  ? "YES" : "NO" ) << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;

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
  FitParamList * fparm0 = util::NamedParameterList("nuosc3f_12fix"); //only used to initialise branch
  FitParamList * fparm00= util::NamedParameterList("nuosc3f_12fix"); //only used to initialise branch
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = util::NamedParameterList(fit_param_list); //first fit (first octant)
  FitParamList * fparm2 = util::NamedParameterList(fit_param_list); //second fit (second octant)
  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  //make an assumption for deltacp
  if(fparm->IsFixed(conv::kPhypCPV    )) {
    fparm ->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm1->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm2->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
  }

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

  double fittime1, fittime2;
  int fitstatus1, fitstatus2, fitstatus;
  int fitstatus_err1, fitstatus_err2, fitstatus_err;
  double chisq1, chisq2, chisq;
  double inparam1, inparam2;
  int nbins_param1, nbins_param2;

  // Create file to store results
  string oscstr = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("%s-run_%d-%s.root", tag, run, oscstr.c_str()), "recreate");

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
  fitnt->Branch("inparam1",&inparam1,"inparam1/D");
  fitnt->Branch("inparam2",&inparam2,"inparam2/D");
  fitnt->Branch("nbins_param1",&nbins_param1,"nbins_param1/I");
  fitnt->Branch("nbins_param2",&nbins_param2,"nbins_param2/I");
  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");
  
  /////////////////////////////////////////////////////
  //                                                 //
  //               MAKE TOYS FOR SINGLE HP           //
  //                                                 //
  /////////////////////////////////////////////////////
  
  // Loop over toy-MC experiments
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     nbins_param1 = max_bin_param1;
     nbins_param2 = max_bin_param2;

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
     }
     else
     if(hyp == kHypNuMuDisapp2F) {
       fparm1->Update(kPhypMixing23, input_mixing23);
       fparm1->Update(kPhypDm23sq,   input_dm23sq);
     }

     // Reset all output spectra
     dat_array ->Clear();
     pdf_array ->Clear();
     pdf_array1->Clear();
     pdf_array2->Clear();

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
     
     /////////////////////////////////////////////////////
     //                                                 //
     //               FIT IN EACH GRID POINT            //
     //                                                 //
     /////////////////////////////////////////////////////
     
     double mixing13_min = 0.0;
     double mixing13_max = 0.1;
     double mixing13_step = 0.005;
     
     if(reactor){
       mixing13_min = 0.015;
       mixing13_max = 0.035;
       mixing13_step = 0.001;
     }
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
       mixing13_min = 0.0;
       mixing13_max = 0.4;
       mixing13_step = 0.02;
     }
     if((Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) && reactor){
       mixing13_min = 0.06;
       mixing13_max = 0.22;//0.14;
       mixing13_step = 0.004;
     }
     
     double mixing23_min = 0.3;
     double mixing23_max = 0.7;
     //double mixing23_step = 0.01;
     double mixing23_step = 0.01;

     double dmsq_min = 2.0E-3;
     double dmsq_max = 3.0E-3;
     //double dmsq_step = 0.05E-3;
     double dmsq_step = (dmsq_max-dmsq_min)/(max_bin_param1 - min_bin_param1);
     
     double deltacp_min = -conv::kpi;
     double deltacp_max = conv::kpi;
     //double deltacp_step = (2*conv::kpi)/50.;
     double deltacp_step = (deltacp_max-deltacp_min)/(max_bin_param1 - min_bin_param1);

     int parvalue1 = -999;
     int parvalue2 = -999;
     
     //loop for the first parameter (if any)
     for (int i = min_bin_param1; i <= max_bin_param1; i++){
       if(param1 == "mixing13") {
	 parvalue1 = conv::kPhypMixing13;
       	 fparm0->SetInitValue(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
	 fparm0->Update(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
	 fparm0->FixParameter(conv::kPhypMixing13);
	 fparm00->SetInitValue(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
	 fparm00->Update(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
	 fparm00->FixParameter(conv::kPhypMixing13);
	 fparm1->SetInitValue(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
	 fparm1->Update(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
	 fparm1->FixParameter(conv::kPhypMixing13);
	 fparm2->SetInitValue(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
	 fparm2->Update(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
	 fparm2->FixParameter(conv::kPhypMixing13);
       }
       if(param1 == "mixing23") {
	 parvalue1 = conv::kPhypMixing23;
	 fparm0->SetInitValue(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
	 fparm0->Update(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
	 fparm0->FixParameter(conv::kPhypMixing23);
	 fparm00->SetInitValue(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
	 fparm00->Update(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
	 fparm00->FixParameter(conv::kPhypMixing23);
	 fparm1->SetInitValue(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
	 fparm1->Update(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
	 fparm1->FixParameter(conv::kPhypMixing23);
	 fparm2->SetInitValue(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
	 fparm2->Update(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
	 fparm2->FixParameter(conv::kPhypMixing23);
       }
       if(param1 == "dmsq") {
	 parvalue1 = conv::kPhypDm23sq;
	 fparm0->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
	 fparm0->Update(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
	 fparm0->FixParameter(conv::kPhypDm23sq);
	 fparm00->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
	 fparm00->Update(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
	 fparm00->FixParameter(conv::kPhypDm23sq);
	 fparm1->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
	 fparm1->Update(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
	 fparm1->FixParameter(conv::kPhypDm23sq);
	 fparm2->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
	 fparm2->Update(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
	 fparm2->FixParameter(conv::kPhypDm23sq);
       }
       if(param1 == "deltacp") {
	 parvalue1 = conv::kPhypCPV;
	 fparm0->SetInitValue(conv::kPhypCPV,(deltacp_min + (i-1)*deltacp_step));
	 fparm0->Update(conv::kPhypCPV,(deltacp_min + (i-1)*deltacp_step));
	 fparm0->FixParameter(conv::kPhypCPV);
	 fparm00->SetInitValue(conv::kPhypCPV,(deltacp_min + (i-1)*deltacp_step));
	 fparm00->Update(conv::kPhypCPV,(deltacp_min + (i-1)*deltacp_step));
	 fparm00->FixParameter(conv::kPhypCPV);
	 fparm1->SetInitValue(conv::kPhypCPV,(deltacp_min + (i-1)*deltacp_step));
	 fparm1->Update(conv::kPhypCPV,(deltacp_min + (i-1)*deltacp_step));
	 fparm1->FixParameter(conv::kPhypCPV);
	 fparm2->SetInitValue(conv::kPhypCPV,(deltacp_min + (i-1)*deltacp_step));
	 fparm2->Update(conv::kPhypCPV,(deltacp_min + (i-1)*deltacp_step));
	 fparm2->FixParameter(conv::kPhypCPV);
       }
       //loop for the second parameter (if any)
       for (int j = min_bin_param2; j <= max_bin_param2; j++){
	 if(param2 == "mixing13") {
	   parvalue2 = conv::kPhypMixing13;
	   fparm0->SetInitValue(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	   fparm0->Update(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	   fparm0->FixParameter(conv::kPhypMixing13);
	   fparm00->SetInitValue(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	   fparm00->Update(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	   fparm00->FixParameter(conv::kPhypMixing13);
	   fparm1->SetInitValue(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	   fparm1->Update(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	   fparm1->FixParameter(conv::kPhypMixing13);
	   fparm2->SetInitValue(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	   fparm2->Update(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	   fparm2->FixParameter(conv::kPhypMixing13);
	 }
	 if(param2 == "mixing23") {
	   parvalue2 = conv::kPhypMixing23;
	   fparm0->SetInitValue(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	   fparm0->Update(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	   fparm0->FixParameter(conv::kPhypMixing23);
	   fparm00->SetInitValue(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	   fparm00->Update(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	   fparm00->FixParameter(conv::kPhypMixing23);
	   fparm1->SetInitValue(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	   fparm1->Update(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	   fparm1->FixParameter(conv::kPhypMixing23);
	   fparm2->SetInitValue(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	   fparm2->Update(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	   fparm2->FixParameter(conv::kPhypMixing23);
	 }
	 if(param2 == "dmsq") {
	   parvalue2 = conv::kPhypDm23sq;
	   fparm0->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	   fparm0->Update(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	   fparm0->FixParameter(conv::kPhypDm23sq);
	   fparm00->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	   fparm00->Update(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	   fparm00->FixParameter(conv::kPhypDm23sq);
	   fparm1->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	   fparm1->Update(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	   fparm1->FixParameter(conv::kPhypDm23sq);
	   fparm2->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	   fparm2->Update(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	   fparm2->FixParameter(conv::kPhypDm23sq);
	 }
	 if(param2 == "deltacp") {
	   parvalue2 = conv::kPhypCPV;
	   fparm0->SetInitValue(conv::kPhypCPV,(deltacp_min + (j-1)*deltacp_step));
	   fparm0->Update(conv::kPhypCPV,(deltacp_min + (j-1)*deltacp_step));
	   fparm0->FixParameter(conv::kPhypCPV);
	   fparm00->SetInitValue(conv::kPhypCPV,(deltacp_min + (j-1)*deltacp_step));
	   fparm00->Update(conv::kPhypCPV,(deltacp_min + (j-1)*deltacp_step));
	   fparm00->FixParameter(conv::kPhypCPV);
	   fparm1->SetInitValue(conv::kPhypCPV,(deltacp_min + (j-1)*deltacp_step));
	   fparm1->Update(conv::kPhypCPV,(deltacp_min + (j-1)*deltacp_step));
	   fparm1->FixParameter(conv::kPhypCPV);
	   fparm2->SetInitValue(conv::kPhypCPV,(deltacp_min + (j-1)*deltacp_step));
	   fparm2->Update(conv::kPhypCPV,(deltacp_min + (j-1)*deltacp_step));
	   fparm2->FixParameter(conv::kPhypCPV);
	 }
	 
	 //if theta23 is not fixed then initialize to 0.4 for the first fit (first octant) and 0.6 for the second one
	 if(!(fparm1->IsFixed(conv::kPhypMixing23))) {     
	   fparm0->SetInitValue(conv::kPhypMixing23,0.4);
	   fparm0->Update(conv::kPhypMixing23,0.4);
	   fparm00->SetInitValue(conv::kPhypMixing23,0.6);
	   fparm00->Update(conv::kPhypMixing23,0.6);
	   fparm1->SetInitValue(conv::kPhypMixing23,0.4);
	   fparm1->Update(conv::kPhypMixing23,0.4);
	   fparm2->SetInitValue(conv::kPhypMixing23,0.6);
	   fparm2->Update(conv::kPhypMixing23,0.6);
	 }
	 //if reactor constraint enabled, set nuisance
	 if(reactor) {
	   fparm0->SetNuisance(conv::kPhypMixing13,true);
	   fparm00->SetNuisance(conv::kPhypMixing13,true);
	   fparm1->SetNuisance(conv::kPhypMixing13,true);
	   fparm2->SetNuisance(conv::kPhypMixing13,true);
	 }

	 // Print out the current grid point

	 cout << endl;
	 cout << endl;
	 cout << endl;
	 cout << "/////////////////////////////////////////" << endl;
	 cout << endl;
	 cout << "Grid point:" << endl;
	 cout << endl;

	 conv::PhysParamId_t p1 = fparm1->WhichParam(conv::kPhypCPV);
	 conv::PhysParamId_t p2 = fparm1->WhichParam(conv::kPhypCPV);
	 if(param1!="none"){
	   TString namepar1 = fparm1->Name(p1);
	   cout << namepar1.Data() << " = " << fparm1->CurrValue(p1) << " ";
	 }
	 if(param2!="none"){
	   TString namepar2 = fparm1->Name(p2);
	   cout << namepar2.Data() << " = " << fparm1->CurrValue(p2) << " ";
	 }
	 cout << endl << endl;
	 
	 // Save the grid point in the tree
	 inparam1 = fparm1->CurrValue(p1);  
	 inparam2 = fparm1->CurrValue(p2);

	 //reset systematic variations, otherwise fit is initialized at
	 //best fit values of systematics from previous fit
	 skexp->ResetSystematicVariations();
	 
     	 //first fit with only osc
	 fparm0->ResetCurrToInitValues(); // reset osc params
	 FitInputs * fit_inputs0 = new FitInputs(skobs,skexp,fparm0);
	 FitterABC * fitter0 = util::GetFitter(fit_inputs0,fit_option);
	 fitter0->SetTol(tolerance);
	 fitter0->DoTheFit();
	 
	 skexp->ResetSystematicVariations();
	 
	 fparm1->ResetCurrToInitValues(); // reset osc params
	 // Get fitter
	 FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
	 FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
      
	 // Chisq for true (input) params
	 if(fparm1->Exists(conv::kNuisXSecSpecFunc))
	   fparm1->Update(conv::kNuisXSecSpecFunc, 0); //init value of SpecFunc is 0.5
	 double true_chisq = fit_inputs->Statistic();
	 
	 // Print-out true params
	 cout << "- True oscillation params (below) yield chisq = " << true_chisq
	      << *fparm1;
	 
	 // Do the fit
	 fitter->SetTol(tolerance);
	 fitter->DoTheFit();
	 fitstatus1 = fitter->FitStatusCode("migrad_min");
	 fittime1 = fitter->FitTime();
	 chisq1  = fparm1->GetBestFitStatistic();
	 cout << "chisq = " << chisq1 << endl;
	 // Save spectra for best-fit prediction
	 if(store_bestfit_spectrum) {
	   skexp->Write(pdf_array1, "pdf");
	 }
            
	 //now, only if theta23 is not fixed, repeat the fit at mirror point...
	 if(((Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta)||(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta2Theta)) && (!(fparm1->IsFixed(conv::kPhypMixing23))) ) {
	   skexp->ResetSystematicVariations();
	   //first fit with only osc
	   fparm00->ResetCurrToInitValues(); // reset osc params
	   FitInputs * fit_inputs00 = new FitInputs(skobs,skexp,fparm00);
	   FitterABC * fitter00 = util::GetFitter(fit_inputs00,fit_option);
	   fitter00->SetTol(tolerance);
	   fitter00->DoTheFit();
	   
	   skexp->ResetSystematicVariations();

	   fparm2->ResetCurrToInitValues(); // reset osc params

	   // Get fitter
	   FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
	   FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
	   
	   // Print-out true params
	   cout << "- Second fit param list: " 
		<< *fparm2;
	   
	   // Do the fit
	   fparm2->ResetCurrToInitValues(); // reset osc params
	   fitter2->SetTol(tolerance);
	   fitter2->DoTheFit();
	   fitstatus2 = fitter2->FitStatusCode("migrad_min");
	   fittime2 = fitter2->FitTime();
	   chisq2  = fparm2->GetBestFitStatistic();
	   
	   // Save spectra (best-fit prediction)
	   if(store_bestfit_spectrum) {
	     skexp->Write(pdf_array2, "pdf");
	   }
	 }
	 
	 // Print-out best-fit params, and set the "best of the two fits" branches correctly
	 cout << "- Best-fit oscillation params (below) yield chisq = ";
	 if((Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta || Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)  &&
	    (!(fparm1->IsFixed(conv::kPhypMixing23)))) {
	   
	   if(chisq2 < chisq1) {
	     fitstatus = fitstatus2;
	     fitstatus_err = fitstatus_err2;
	     chisq = chisq2;	     
	     cout << chisq2 << " from the second seed fit" << endl;
	     fparm2->ShowOnlyBestFitValuesInNextPrint();
	     cout << *fparm2;
	     br_fparm->SetAddress(&fparm2);
	     br_bf_array->SetAddress(&pdf_array2);
	   }
	   else {
	     fitstatus = fitstatus1;
	     fitstatus_err = fitstatus_err1;
	     chisq = chisq1;	     
	     cout << chisq1 << " from the first seed fit" << endl;
	     fparm1->ShowOnlyBestFitValuesInNextPrint();
	     cout << *fparm1;
	     br_fparm->SetAddress(&fparm1);
	     br_bf_array->SetAddress(&pdf_array1);
	   }
	 }
	 else {
	   fitstatus = fitstatus1;
	   fitstatus_err = fitstatus_err1;
	   chisq = chisq1;	   
	   cout << chisq1 << endl;
	   fparm1->ShowOnlyBestFitValuesInNextPrint();
	   cout << *fparm1;
	   br_fparm->SetAddress(&fparm1);
	   br_bf_array->SetAddress(&pdf_array1);
	 }
	 	 
	 // Add all results in the tree
	 fitnt->Fill();
       
	 // Clean-up
	 delete fit_inputs;
	 delete fitter;
	 delete fit_inputs0;
	 delete fitter0;
	 if((Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta ||Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)  &&
	    (!(fparm1->IsFixed(conv::kPhypMixing23)))) {
	   delete fit_inputs2;
	   delete fitter2;
	   delete fit_inputs00;
	   delete fitter00;
	 }      
       
       } // loop grid param2
     } // loop grid param1
     
  }
  f.cd();
  fitnt->Write();
  f.Close();
  
  return;
}
