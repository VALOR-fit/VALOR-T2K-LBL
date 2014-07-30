//
// ------------------------------------------------------------------------------------------------
// Run fits to experiments read from an input XML file and store best-fit parameters and spectra.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//
 
void make_toygrid_toys_joint
(
 const char * xmlfile = "Asimov/toymc_Asimov_dcp0_NH-0.xml",             /* input XML file with toy MC experiments */
 //const char * xmlfile = "data_newstruct.xml",                                                /* input XML file with toy MC experiments */
 //const char * xmlfile = "toymc_FHCpRHC_DataRun1to4_MCRun5and6_NH-0.xml",                 /* input XML file with toy MC experiments */
 
 const char * tag = "prova",                                                    /* */
 int          run = 0,                                                    /* */
 int          first_experiment = 0,                                   /* */
 int          last_experiment  = 0,                                   /* */
 
 double       input_mixing12          = -99, //0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13          = -99, //0.0251,                                /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23          = -99, //0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp           = -99, //-3.1,                                /* deltacp, if -999 then randomize */
 double       input_dm12sq            = -99, //7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
 double       input_dm23sq            = -99, //2.32E-3,                             /* |Dm^2_{23}| */

 //const char * config_file      = "analysis_configuration_run1to6.xml",        /* analysis configuration file */
 //const char * config_file      = "analysis_configuration_run5and6.xml",        /* analysis configuration file */
 //const char * config_file      = "analysis_configuration_test.xml",        /* analysis configuration file */
 const char * config_file        = "analysis_configuration_test_multiring_mrmc.xml",   /* analysis configuration file */

 bool         stat_fluctuations       = true,                                /* include statistical fluctuation in toy MC expts? */
 int          tweak_systematics       = 0,                                /* generate toy MC with randomized systematics? */

 const char * fit_param_list   = "nuosc3f_12fix",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
 //const char * fit_param_list   = "nuosc3f_12fix+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
 
 //const char * fit_option       = "minuit+hesse",                            /* fitter option */
 const char * fit_option       = "minuit",                            /* fitter option */
 
 bool         store_toy_expts         = false,                         /* disable this if you don't want to store the toy MC expts     */
 bool         store_bestfit_spectrum  = false,                         /* disable this if you don't want to store the best-fit spectra */
 
 const char * dmsqscheme          = "abs23",                               /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy           = "normal",                              /* normal or inverted */
 const char * mixparamchoice      = "sinsqtheta",                          /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant      = "one",                                 /* one==kNuOctantOne two==kNuOctantTwo */
 double       dcp_value           = 0,                                      /* value of deltacp (in units of PI) to assume (when deltacp is not a free parameter) */
 bool         reactor             = true,                               /* enable this if you want to apply the reactor constraint */
 double       tolerance           = 0.1,                                 /* value for the fitter tolerance */
 bool         using_xml_file      = true,                                   /* if false, skobs reads xmlfile as a root file */
 
 //const char * param1                  = "mixing23",                        /* first parameter to be fitting vs - default: none for a single fit */ 
 //const char * param2                  = "dmsq",                        /* second parameter to be fitting vs (for 2D) - default: none for a single fit or 1D */ 
 //const char * param3                  = "mixing13",                        /* second parameter to be fitting vs (for 2D) - default: none for a single fit or 1D */ 
 const char * param1                  = "deltacp",                        /* first parameter to be fitting vs - default: none for a single fit */ 
 const char * param2                  = "none",                        /* second parameter to be fitting vs (for 2D) - default: none for a single fit or 1D */ 
 const char * param3                  = "none",                        /* second parameter to be fitting vs (for 2D) - default: none for a single fit or 1D */ 
 
 int          nbins_tot_param1            = 10
,                              /* tot # grid values for param1 (start from bin=1) --> part of the grid --> step = nbins-1 */
 int          nbins_tot_param2            = 1,                              /* tot # grid values for param2 (start from bin=1) --> part of the grid --> step = nbins-1 */
 int          nbins_tot_param3            = 1,                              /* tot # grid values for param3 (start from bin=1) --> part of the grid --> step = nbins-1 */

 int          min_bin_param1          = 1,                                /* minimum bin for param1 */
 int          max_bin_param1          = 10,                                /* maximum bin for param1 */
 int          min_bin_param2          = 1,                                /* minimum bin for param2 */
 int          max_bin_param2          = 1,                                  /* maximum bin for param2 */
 int          min_bin_param3          = 1,                                /* minimum bin for param2 */
 int          max_bin_param3          = 1                                  /* maximum bin for param2 */
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  double bestfit_chisq,bestfit_chisq2;

  // Analysis configuration
  TString config   = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);
  
  FileStat_t dummy;
  
  // Check for config file, and XML data file.
  // If not found in inputs or config, see if it works as a normal path
  if (gSystem->GetPathInfo(config.Data(), dummy)) 
    {
      config = config_file;
      if (gSystem->GetPathInfo(config.Data(), dummy)) {
	cerr<<"Config file not found: "<<config_file<<endl;
	exit(-1);
      }
    }
  
  if(using_xml_file){
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
      {
	xmlfile2 = xmlfile;
	if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
	  cerr<<"Data input file not found: "<<xmlfile<<endl;
	  exit(-1);
	}
      }
  }

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

  cout << "- Toy-MC generation" << endl;
  //cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - Number of toy experiments: " << last_experiment-first_experiment+1 << endl; 
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
  cout << "   - Fitter tolerance: " << tolerance << endl;
  cout << "  Applying reactor constraint: " << ((reactor)  ? "YES" : "NO" ) << endl;
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

  // Fit param list
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = util::NamedParameterList(fit_param_list);
  FitParamList * fparm2 = util::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // ///// TO REMOVE ????
  // //make an assumption for deltacp
  // if(fparm->IsFixed(conv::kPhypCPV    )) {
  //   fparm ->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
  //   fparm1->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
  //   fparm2->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
  // }
  // /////////////////////

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());
    
  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each expt.
  TObjArray * dat_array = new TObjArray; //only used to initialise branch
  TObjArray * pdf_array = new TObjArray; //only used to initialise branch
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);
  TObjArray * pdf_array1 = new TObjArray;
  pdf_array1->SetOwner(false);
  TObjArray * pdf_array2 = new TObjArray;
  pdf_array2->SetOwner(false);

  // Create file to store results
  //TFile f(Form("fitgrid-%s_expts_%d_%d.root", tag, run, first_experiment, last_experiment), "recreate");
  char outtreename[200];
  sprintf(outtreename,"fitgrid_%s_%s_%s_%s_BinsTot_%d_BinsTot_%d_BinsTot_%d_BinsPart_%d_%d_BinsPart_%d_%d_BinsPart_%d_%d.root",
	  tag,param1,param2,param3,
	  nbins_tot_param1,nbins_tot_param2,nbins_tot_param3,
	  min_bin_param1,max_bin_param1,min_bin_param2,max_bin_param2,min_bin_param3,max_bin_param3
	  );
  TFile f(outtreename, "recreate");
  
  double fittime1, fittime2;
  int fitstatus1, fitstatus2, fitstatus;
  int fitstatus_err1, fitstatus_err2, fitstatus_err;
  double chisq1, chisq2, chisq;
  double inparam1, inparam2, inparam3;
  int nbins_param1, nbins_param2, nbins_param3;
  int max_bin_param1, min_bin_param1, max_bin_param2, min_bin_param2, max_bin_param3, min_bin_param3;
  double max_val_param1, min_val_param1, max_val_param2, min_val_param2, max_val_param3, min_val_param3;
  double step1, step2, step3;

  TH1D *hLike1D = new TH1D(); // dcp
  TH2D *hLike2D_12 = new TH2D(); // 23, |dm2|
  TH2D *hLike2D_13 = new TH2D(); // 23, 13
  TH2D *hLike2D_23 = new TH2D(); // |dm2|, 13
  TH3D *hLike3D = new TH3D(); // 13, 23, |dm2|

  //string parname1, parname2;
  std::string *parname1 = new std::string();
  std::string *parname2 = new std::string();
  std::string *parname3 = new std::string();
  parname1->append(param1);
  parname2->append(param2);
  parname3->append(param3);

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("best_fit_params_fit1", "valor::FitParamList",  &fparm1,     1280000,0);
  fitnt->Branch("best_fit_params_fit2", "valor::FitParamList",  &fparm2,     1280000,0);
  //fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);
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
  if(param1!="none" || param2!="none" || param3!="none"){
    fitnt->Branch("inparam1",&inparam1,"inparam1/D");
    fitnt->Branch("inparam2",&inparam2,"inparam2/D");
    fitnt->Branch("inparam3",&inparam3,"inparam3/D");
    fitnt->Branch("nbins_param1",&nbins_param1,"nbins_param1/I");
    fitnt->Branch("nbins_param2",&nbins_param2,"nbins_param2/I");
    fitnt->Branch("nbins_param3",&nbins_param3,"nbins_param3/I");
    fitnt->Branch("nbins_tot_param1",&nbins_tot_param1,"nbins_tot_param1/I");
    fitnt->Branch("nbins_tot_param2",&nbins_tot_param2,"nbins_tot_param2/I");
    fitnt->Branch("nbins_tot_param3",&nbins_tot_param3,"nbins_tot_param3/I");
    fitnt->Branch("max_val_param1",&max_val_param1,"max_val_param1/D");
    fitnt->Branch("min_val_param1",&min_val_param1,"min_val_param1/D");
    fitnt->Branch("max_val_param2",&max_val_param2,"max_val_param2/D");
    fitnt->Branch("min_val_param2",&min_val_param2,"min_val_param2/D");
    fitnt->Branch("max_val_param3",&max_val_param3,"max_val_param3/D");
    fitnt->Branch("min_val_param3",&min_val_param3,"min_val_param3/D");
    fitnt->Branch("step1",&step1,"step1/D");
    fitnt->Branch("step2",&step2,"step2/D");
    fitnt->Branch("step3",&step3,"step3/D");
    fitnt->Branch("parname1", &parname1); // name of the grid parameter 1
    fitnt->Branch("parname2", &parname2); // name of the grid parameter 2    
    fitnt->Branch("parname3", &parname3); // name of the grid parameter 3    
    fitnt->Branch("hLike1D","TH1D",&hLike1D);
    fitnt->Branch("hLike2D_12","TH2D",&hLike2D_12);
    fitnt->Branch("hLike2D_13","TH2D",&hLike2D_13);
    fitnt->Branch("hLike2D_23","TH2D",&hLike2D_23);
    fitnt->Branch("hLike3D","TH3D",&hLike3D);
  }
  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");
  

  //Open the ROOT file that has the input PhysParamList (if applicable)
  
  // TFile * f_sysparm = 0;
  // if(!using_xml_file) {
  //   f_sysparm = new TFile(xmlfile2.Data(),"READ");
  //   if(f_sysparm->IsZombie()) {
  //     cerr << "ROOT file not found. Exiting..." << endl;
  //     exit(1);
  //   }
  // }

  // Loop over toy-MC experiments
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {
         
    nbins_param1 = (max_bin_param1-min_bin_param1)+1;
    nbins_param2 = (max_bin_param2-min_bin_param2)+1;
    nbins_param3 = (max_bin_param3-min_bin_param3)+1;
     
    if( nbins_param1 > nbins_tot_param1 ||
	nbins_param2 > nbins_tot_param2 ||
	nbins_param3 > nbins_tot_param3 ) {

      cout << "nbins_param1 = " << nbins_param1 << endl;
      cout << "nbins_param2 = " << nbins_param2 << endl;
      cout << "nbins_param3 = " << nbins_param3 << endl;
      cout << endl;
      cout << "nbins_tot_param1 = " << nbins_tot_param1 << endl;
      cout << "nbins_tot_param2 = " << nbins_tot_param2 << endl;
      cout << "nbins_tot_param3 = " << nbins_tot_param3 << endl;
      
      cout << "Partial # of grid bins > Total # of grid bins!!!" << endl;
      exit(1);
    }

    cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;
      
    if(!using_xml_file){
      // randomize specified osc params for current toy MC expt
      if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
      if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
      if (randomize_deltacp ) { input_deltacp  = util::RandomizeCPV();  }
      if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }
    }

     // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
     if(hyp == kHypOsc3F) {
       fparm1->Update(kPhypMixing12, input_mixing12);
       fparm1->Update(kPhypMixing13, input_mixing13);
       fparm1->Update(kPhypMixing23, input_mixing23);
       fparm1->Update(kPhypCPV    ,  input_deltacp);
       fparm1->Update(kPhypDm12sq,   input_dm12sq);
       fparm1->Update(kPhypDm23sq,   input_dm23sq);
       //////////////////////////////////////////////
       fparm2->Update(kPhypMixing12, input_mixing12);
       fparm2->Update(kPhypMixing13, input_mixing13);
       fparm2->Update(kPhypMixing23, input_mixing23);
       fparm2->Update(kPhypCPV    ,  input_deltacp);
       fparm2->Update(kPhypDm12sq,   input_dm12sq);
       fparm2->Update(kPhypDm23sq,   input_dm23sq);
       //////////////////////////////////////////////
     }
     else
     if(hyp == kHypNuMuDisapp2F) {
       fparm1->Update(kPhypMixing23, input_mixing23);
       fparm1->Update(kPhypDm23sq,   input_dm23sq);
       //////////////////////////////////////////////
       fparm2->Update(kPhypMixing23, input_mixing23);
       fparm2->Update(kPhypDm23sq,   input_dm23sq);
       //////////////////////////////////////////////
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
     SuperkObservation *skobs = new SuperkObservation;

     if(!using_xml_file){
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
     }// if(!using_xml_file)
     // Get superk spectrum observation (toy experiment)
     //SuperkObservation * skobs = new SuperkObservation;
     else 
       skobs->Read(xmlfile2.Data(),iexperiment,using_xml_file);

    // Set exclusion grid points
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
    if(param1 == "mixing13")      mixing13_step = (mixing13_max-mixing13_min)/(nbins_tot_param1-1);
    else if(param2 == "mixing13") mixing13_step = (mixing13_max-mixing13_min)/(nbins_tot_param2-1);     
    else if(param3 == "mixing13") mixing13_step = (mixing13_max-mixing13_min)/(nbins_tot_param3-1);     
    
    double mixing23_min = 0.3;
    double mixing23_max = 0.7;
    double mixing23_step = 0.01;
    if(param1 == "mixing23")      mixing23_step = (mixing23_max-mixing23_min)/(nbins_tot_param1-1);
    else if(param2 == "mixing23") mixing23_step = (mixing23_max-mixing23_min)/(nbins_tot_param2-1);
    else if(param3 == "mixing23") mixing23_step = (mixing23_max-mixing23_min)/(nbins_tot_param3-1);
    
    double dmsq_min = 2.0E-3;
    double dmsq_max = 3.0E-3;
    double dmsq_step = 0.05E-3;
    if(param1 == "dmsq")      dmsq_step = (dmsq_max-dmsq_min)/(nbins_tot_param1-1);
    else if(param2 == "dmsq") dmsq_step = (dmsq_max-dmsq_min)/(nbins_tot_param2-1);          
    else if(param3 == "dmsq") dmsq_step = (dmsq_max-dmsq_min)/(nbins_tot_param3-1);          

    double deltacp_min = -conv::kpi;
    double deltacp_max = conv::kpi;
    double deltacp_step = (2*conv::kpi)/50.;
    if(param1 == "deltacp")      deltacp_step = (deltacp_max-deltacp_min)/(nbins_tot_param1-1);
    else if(param2 == "deltacp") deltacp_step = (deltacp_max-deltacp_min)/(nbins_tot_param2-1);     
    else if(param3 == "deltacp") deltacp_step = (deltacp_max-deltacp_min)/(nbins_tot_param3-1);     

    if(param1=="mixing13"){
      step1 = mixing13_step;
      max_val_param1 = mixing13_max;
      min_val_param1 = mixing13_min;
    }
    else if(param1=="mixing23"){
      step1 = mixing23_step;
      max_val_param1 = mixing23_max;
      min_val_param1 = mixing23_min;
    }
    else if(param1=="dmsq"){
      step1 = dmsq_step;
      max_val_param1 = dmsq_max;
      min_val_param1 = dmsq_min;
    }
    else if(param1=="deltacp"){
      step1 = deltacp_step;
      max_val_param1 = deltacp_max;
      min_val_param1 = deltacp_min;
    }
    else if(param1=="none"){
      step1 = 0;
      max_val_param1 = 0;
      min_val_param1 = 0;
    }
    else{
      cout << "You did not choose the right oscillation parameter 1 !!!" << endl;
      exit(0);
    }
   
    if(param2=="mixing13"){
      step2 = mixing13_step;
      max_val_param2 = mixing13_max;
      min_val_param2 = mixing13_min;
    }
    else if(param2=="mixing23"){
      step2 = mixing23_step;
      max_val_param2 = mixing23_max;
      min_val_param2 = mixing23_min;
    }
    else if(param2=="dmsq"){
      step2 = dmsq_step;
      max_val_param2 = dmsq_max;
      min_val_param2 = dmsq_min;
    }
    else if(param2=="deltacp"){
      step2 = deltacp_step;
      max_val_param2 = deltacp_max;
      min_val_param2 = deltacp_min;
    }
    else if(param2=="none"){
      step2 = 0;
      max_val_param2 = 0;
      min_val_param2 = 0;
    }
    else{
      cout << "You did not choose the right oscillation parameter 2 !!!" << endl;
      exit(0);
    }
    
    if(param3=="mixing13"){
      step3 = mixing13_step;
      max_val_param3 = mixing13_max;
      min_val_param3 = mixing13_min;
    }
    else if(param3=="mixing23"){
      step3 = mixing23_step;
      max_val_param3 = mixing23_max;
      min_val_param3 = mixing23_min;
    }
    else if(param3=="dmsq"){
      step3 = dmsq_step;
      max_val_param3 = dmsq_max;
      min_val_param3 = dmsq_min;
    }
    else if(param3=="deltacp"){
      step3 = deltacp_step;
      max_val_param3 = deltacp_max;
      min_val_param3 = deltacp_min;
    }
    else if(param3=="none"){
      step3 = 0;
      max_val_param3 = 0;
      min_val_param3 = 0;
    }
    else{
      cout << "You did not choose the right oscillation parameter 3 !!!" << endl;
      exit(0);
    }

    
    // Set histograms of likelihood grid
    if(param1!="none"){
      int n_param1 = (max_val_param1 - min_val_param1)/step1 + 1; // +1 for last point
      double min1 = min_val_param1;
      double max1 = max_val_param1+step1;
      hLike1D = new TH1D("hLike1D","hLike1D",n_param1,min1,max1);
    }
    if(param1!="none" && param2!="none"){
      hLike2D_12->Clear();
      int n_param1 = (max_val_param1 - min_val_param1)/step1 + 1; // +1 for last point
      double min1 = min_val_param1;
      double max1 = max_val_param1+step1;
      int n_param2 = (max_val_param2 - min_val_param2)/step2 + 1; // +1 for last point
      double min2 = min_val_param2;
      double max2 = max_val_param2+step2;
      hLike2D_12 = new TH2D("hLike2D_12","hLike2D_12",n_param1,min1,max1,n_param2,min2,max2);
    }
    if(param1!="none" && param3!="none"){
      hLike2D_13->Clear();
      int n_param1 = (max_val_param1 - min_val_param1)/step1 + 1; // +1 for last point
      double min1 = min_val_param1;
      double max1 = max_val_param1+step1;
      int n_param3 = (max_val_param3 - min_val_param3)/step3 + 1; // +1 for last point
      double min3 = min_val_param3;
      double max3 = max_val_param3+step3;
      hLike2D_13 = new TH2D("hLike2D_13","hLike2D_13",n_param1,min1,max1,n_param3,min3,max3);
    }
    if(param2!="none" && param3!="none"){
      hLike2D_23->Clear();
      int n_param2 = (max_val_param2 - min_val_param2)/step2 + 1; // +1 for last point
      double min2 = min_val_param2;
      double max2 = max_val_param2+step2;
      int n_param3 = (max_val_param3 - min_val_param3)/step3 + 1; // +1 for last point
      double min3 = min_val_param3;
      double max3 = max_val_param3+step3;
      hLike2D_23 = new TH2D("hLike2D_23","hLike2D_23",n_param2,min2,max2,n_param3,min3,max3);
    }
    if(param1!="none" && param2!="none" && param3!="none"){
      hLike3D->Clear();
      int n_param1 = (max_val_param1 - min_val_param1)/step1 + 1; // +1 for last point
      double min1 = min_val_param1;
      double max1 = max_val_param1+step1;
      int n_param2 = (max_val_param2 - min_val_param2)/step2 + 1; // +1 for last point
      double min2 = min_val_param2;
      double max2 = max_val_param2+step2;
      int n_param3 = (max_val_param3 - min_val_param3)/step3 + 1; // +1 for last point
      double min3 = min_val_param3;
      double max3 = max_val_param3+step3;
      hLike3D = new TH3D("hLike3D","hLike3D",n_param1,min1,max1,n_param2,min2,max2,n_param3,min3,max3);
    }
    
    //
    
    int parvalue1 = -999;
    int parvalue2 = -999;
    int parvalue3 = -999;
    //loop for the first parameter (if any)
    for (int i = min_bin_param1; i <= max_bin_param1; i++){
      if(param1 == "mixing13") {
	parvalue1 = conv::kPhypMixing13;
	double currval = mixing13_min + (i-1)*mixing13_step;
	fparm1->SetInitValue(conv::kPhypMixing13,currval);
	fparm1->Update(conv::kPhypMixing13,currval);
	fparm1->FixParameter(conv::kPhypMixing13);
	fparm2->SetInitValue(conv::kPhypMixing13,currval);
	fparm2->Update(conv::kPhypMixing13,currval);
	fparm2->FixParameter(conv::kPhypMixing13);
      }
      if(param1 == "mixing23") {
	parvalue1 = conv::kPhypMixing23;
	double currval = mixing23_min + (i-1)*mixing23_step;
	fparm1->SetInitValue(conv::kPhypMixing23,currval);
	fparm1->Update(conv::kPhypMixing23,currval);
	fparm1->FixParameter(conv::kPhypMixing23);
	fparm2->SetInitValue(conv::kPhypMixing23,currval);
	fparm2->Update(conv::kPhypMixing23,currval);
	fparm2->FixParameter(conv::kPhypMixing23);
      }
      if(param1 == "dmsq") {
	parvalue1 = conv::kPhypDm23sq;
	double currval = dmsq_min + (i-1)*dmsq_step;
      	fparm1->SetInitValue(conv::kPhypDm23sq,currval);
      	fparm1->Update(conv::kPhypDm23sq,currval);
      	fparm1->FixParameter(conv::kPhypDm23sq);
      	fparm2->SetInitValue(conv::kPhypDm23sq,currval);
      	fparm2->Update(conv::kPhypDm23sq,currval);
      	fparm2->FixParameter(conv::kPhypDm23sq);
      }
      if(param1 == "deltacp") {
	parvalue1 = conv::kPhypCPV;
	double currval = deltacp_min + (i-1)*deltacp_step;
	fparm1->SetInitValue(conv::kPhypCPV,currval);
	fparm1->Update(conv::kPhypCPV,currval);
	fparm1->FixParameter(conv::kPhypCPV);
	fparm2->SetInitValue(conv::kPhypCPV,currval);
	fparm2->Update(conv::kPhypCPV,currval);
	fparm2->FixParameter(conv::kPhypCPV);
      }
      //loop for the second parameter (if any)
      for (int j = min_bin_param2; j <= max_bin_param2; j++){
	if(param2 == "mixing13") {
	  parvalue2 = conv::kPhypMixing13;
	  double currval = mixing13_min + (j-1)*mixing13_step;
	  fparm1->SetInitValue(conv::kPhypMixing13,currval);
	  fparm1->Update(conv::kPhypMixing13,currval);
	  fparm1->FixParameter(conv::kPhypMixing13);
	  fparm2->SetInitValue(conv::kPhypMixing13,currval);
	  fparm2->Update(conv::kPhypMixing13,currval);
	  fparm2->FixParameter(conv::kPhypMixing13);
	}
	if(param2 == "mixing23") {
	  parvalue2 = conv::kPhypMixing23;
	  double currval = mixing23_min + (j-1)*mixing23_step;
	  fparm1->SetInitValue(conv::kPhypMixing23,currval);
	  fparm1->Update(conv::kPhypMixing23,currval);
	  fparm1->FixParameter(conv::kPhypMixing23);
	  fparm2->SetInitValue(conv::kPhypMixing23,currval);
	  fparm2->Update(conv::kPhypMixing23,currval);
	  fparm2->FixParameter(conv::kPhypMixing23);
	}
	if(param2 == "dmsq") {
	  parvalue2 = conv::kPhypDm23sq;
	  double currval = dmsq_min + (j-1)*dmsq_step;
	  fparm1->SetInitValue(conv::kPhypDm23sq,currval);
	  fparm1->Update(conv::kPhypDm23sq,currval);
	  fparm1->FixParameter(conv::kPhypDm23sq);
	  fparm2->SetInitValue(conv::kPhypDm23sq,currval);
	  fparm2->Update(conv::kPhypDm23sq,currval);
	  fparm2->FixParameter(conv::kPhypDm23sq);
	}
	if(param2 == "deltacp")	 {
	  parvalue2 = conv::kPhypCPV; 
	  double currval = deltacp_min + (j-1)*deltacp_step;
	  fparm1->SetInitValue(conv::kPhypCPV,currval);
	  fparm1->Update(conv::kPhypCPV,currval);
	  fparm1->FixParameter(conv::kPhypCPV);
	  fparm2->SetInitValue(conv::kPhypCPV,currval);
	  fparm2->Update(conv::kPhypCPV,currval);
	  fparm2->FixParameter(conv::kPhypCPV);
	}
	//loop for the third parameter (if any)
	for (int k = min_bin_param3; k <= max_bin_param3; k++){
	  if(param3 == "mixing13") {
	    parvalue3 = conv::kPhypMixing13;
	    double currval = mixing13_min + (k-1)*mixing13_step;
	    fparm1->SetInitValue(conv::kPhypMixing13,currval);
	    fparm1->Update(conv::kPhypMixing13,currval);
	    fparm1->FixParameter(conv::kPhypMixing13);
	    fparm2->SetInitValue(conv::kPhypMixing13,currval);
	    fparm2->Update(conv::kPhypMixing13,currval);
	    fparm2->FixParameter(conv::kPhypMixing13);
	  }
	  if(param3 == "mixing23") {
	    parvalue3 = conv::kPhypMixing23;
	    double currval = mixing23_min + (k-1)*mixing23_step;
	    fparm1->SetInitValue(conv::kPhypMixing23,currval);
	    fparm1->Update(conv::kPhypMixing23,currval);
	    fparm1->FixParameter(conv::kPhypMixing23);
	    fparm2->SetInitValue(conv::kPhypMixing23,currval);
	    fparm2->Update(conv::kPhypMixing23,currval);
	    fparm2->FixParameter(conv::kPhypMixing23);
	  }
	  if(param3 == "dmsq") {
	    parvalue3 = conv::kPhypDm23sq;
	    double currval = dmsq_min + (k-1)*dmsq_step;
	    fparm1->SetInitValue(conv::kPhypDm23sq,currval);
	    fparm1->Update(conv::kPhypDm23sq,currval);
	    fparm1->FixParameter(conv::kPhypDm23sq);
	    fparm2->SetInitValue(conv::kPhypDm23sq,currval);
	    fparm2->Update(conv::kPhypDm23sq,currval);
	    fparm2->FixParameter(conv::kPhypDm23sq);
	  }
	  if(param3 == "deltacp") {
	    parvalue3 = conv::kPhypCPV;
	    double currval = deltacp_min + (k-1)*deltacp_step;
	    fparm1->SetInitValue(conv::kPhypCPV,currval);
	    fparm1->Update(conv::kPhypCPV,currval);
	    fparm1->FixParameter(conv::kPhypCPV);
	    fparm2->SetInitValue(conv::kPhypCPV,currval);
	    fparm2->Update(conv::kPhypCPV,currval);
	    fparm2->FixParameter(conv::kPhypCPV);
	  }

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

	  conv::PhysParamId_t p1;
	  conv::PhysParamId_t p2;
	  conv::PhysParamId_t p3;
	  if(param1 != "none" || param2 != "none" || param3 != "none"){
	    cout << endl;
	    cout << endl;
	    cout << endl;
	    cout << "/////////////////////////////////////////" << endl;
	    cout << endl;
	    cout << "Grid point:" << endl;
	    cout << endl;
	    if(param1!="none"){
	      p1 = fparm1->WhichParam(parvalue1);
	      TString namepar1 = fparm1->Name(p1);
	      cout << namepar1.Data() << " = " << fparm1->CurrValue(p1) << " ";
	      cout << "step1 = " << step1 
		   << " - min_val_param1 = " << min_val_param1 
		   << " - max_val_param1 = " << max_val_param1 
		   << " - parvalue1 = " << parvalue1 
		   << endl;
	    }
	    if(param2!="none"){
	      p2 = fparm1->WhichParam(parvalue2);
	      TString namepar2 = fparm1->Name(p2);
	      cout << namepar2.Data() << " = " << fparm1->CurrValue(p2) << " ";
	      cout << "step2 = " << step2 
		   << " - min_val_param2 = " << min_val_param2 
		   << " - max_val_param2 = " << max_val_param2 
		   << " - parvalue2 = " << parvalue2 
		   << endl;
		}
	    if(param3!="none"){
	      p3 = fparm1->WhichParam(parvalue3);
	      TString namepar3 = fparm1->Name(p3);
	      cout << namepar3.Data() << " = " << fparm1->CurrValue(p3) << " ";
	      cout << "step3 = " << step3 
		   << " - min_val_param3 = " << min_val_param3 
		   << " - max_val_param3 = " << max_val_param3 
		   << " - parvalue3 = " << parvalue3 
		   << endl;
	    }
	    cout << endl;
	    cout << endl;

	    // Save the grid point in the tree
	    if(param1!="none") inparam1 = fparm1->CurrValue(p1);  
	    if(param2!="none") inparam2 = fparm1->CurrValue(p2);
	    if(param3!="none") inparam3 = fparm1->CurrValue(p3);
	  }
	
	  //reset systematic variations, otherwise fit is initialized at
	  //best fit values of systematics from previous fit
	  skexp->ResetSystematicVariations();
	  fparm1->ResetCurrToInitValues(); // reset osc params
	
	  // Get fitter
	  FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
	  FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

	  // Set SF Initial value to 0
	  if(fparm1->Exists(conv::kNuisXSecSpecFunc))
	    fparm1->Update(conv::kNuisXSecSpecFunc, 0); //init value of SpecFunc is 0.5

	  // Chisq for true (input) params
	  double true_chisq = fit_inputs->Statistic();

	  // Print-out true params
	  //fparm1->ShowOnlyCurrValuesInNextPrint();
	  cout << "- True oscillation params (below) yield chisq = " << true_chisq
	       << *fparm1;

	  // Do the fit

	  cout << endl;
	  cout << endl;
	  cout << "First Octant" << endl;
	  cout << endl;
	  cout << endl;

	  //skexp->PrintSystematicVariations();
	  fitter->SetTol(tolerance);
	  fitter->DoTheFit();
	  fitstatus1 = fitter->FitStatusCode("migrad_min");
	  fitstatus_err1 = fitter->FitStatusCode("hesse_err");     

	  double bestfit_chisq  = fparm1->GetBestFitStatistic();
	  chisq1 = bestfit_chisq;
	  double bestfit_chisq2 = 9999.;

	  // Save spectra (data and best-fit prediction)
	  if(store_toy_expts) {
	    skobs->Write(dat_array, "obs");
	  }
	  if(store_bestfit_spectrum) {
	    skexp->Write(pdf_array1, "pdf");
	  }

	  //Do the fit again
	  if(((Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta)||(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta2Theta)) 
	     && (!(fparm1->IsFixed(conv::kPhypMixing23))) ) {

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
       
	    // Get fitter
	    FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
	    FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
	    fitter2->SetTol(tolerance);
	    fitter2->DoTheFit();
	    fitstatus2 = fitter->FitStatusCode("migrad_min");
	    fitstatus_err2 = fitter->FitStatusCode("hesse_err");     
	    bestfit_chisq2 = fparm2->GetBestFitStatistic();
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

	  cout << endl;
	  cout << "Best-fit chisq1 = " << chisq1 << endl;
	  cout << "Best-fit chisq2 = " << chisq2 << endl;
	  cout << "Best-fit chisq = " << chisq << endl;
	  cout << endl;
	
	
	  // Fill hLike grid with chisq

	  if(param1!="none"){
	    double x = fparm1->CurrValue(p1);
	    double y = fparm1->CurrValue(p2);
	    int bin = hLike1D->FindBin(x+conv::kASmallNum,y+conv::kASmallNum);
	    hLike1D->SetBinContent(bin,chisq);
	  }
	  if(param1!="none" && param2!="none"){
	    double x = fparm1->CurrValue(p1);
	    double y = fparm2->CurrValue(p2);
	    int bin = hLike2D_12->FindBin(x+conv::kASmallNum,y+conv::kASmallNum);
	    hLike2D_12->SetBinContent(bin,chisq);
	  }
	  if(param1!="none" && param3!="none"){
	    double x = fparm1->CurrValue(p1);
	    double y = fparm1->CurrValue(p3);
	    int bin = hLike2D_13->FindBin(x+conv::kASmallNum,y+conv::kASmallNum);
	    hLike2D_13->SetBinContent(bin,chisq);
	  }
	  if(param2!="none" && param3!="none"){
	    double x = fparm1->CurrValue(p2);
	    double y = fparm1->CurrValue(p3);
	    int bin = hLike2D_23->FindBin(x+conv::kASmallNum,y+conv::kASmallNum);
	    hLike2D_23->SetBinContent(bin,chisq);
	  }
	  if(param1!="none" && param2!="none" && param3!="none"){
	    double x = fparm1->CurrValue(p1);
	    double y = fparm1->CurrValue(p2);
	    double z = fparm1->CurrValue(p3);
	    int bin = hLike3D->FindBin(x+conv::kASmallNum,y+conv::kASmallNum,z+conv::kASmallNum);
	    hLike3D->SetBinContent(bin,chisq);
	  }
	
	  // Add all results in the tree
	  fitnt->Fill(); 
     
	  // Clean-up
	  //if(using_xml_file) delete skobs;
	  delete fit_inputs;
	  delete fitter;
	  if((Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta ||Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)  &&
	     (!(fparm1->IsFixed(conv::kPhypMixing23)))) {
	    delete fit_inputs2;
	    delete fitter2;
	  }    
	
	} // k 
      } // j
    } // i
    
    // Store histogram in root file for data only --> only 1 toy!!!
    if(last_experiment-first_experiment==0){
      if(param1!="none" || param2!="none" || param3!="none"){
	char foutname[200];
	sprintf(foutname,"LikeGrid_%s_%s_%s_%s_BinsTot_%d_BinsTot_%d_BinsTot_%d_BinsPart_%d_%d_BinsPart_%d_%d_BinsPart_%d_%d.root",
		tag,param1,param2,param3,
		nbins_tot_param1,nbins_tot_param2,nbins_tot_param3,
		min_bin_param1,max_bin_param1,min_bin_param2,max_bin_param2,min_bin_param3,max_bin_param3
		);
	TFile *fout = new TFile(foutname,"RECREATE");
	cout << endl;
	cout << "Store likelihood grid in " << foutname << endl;
	cout << endl;
	hLike1D->Write();
	hLike2D_12->Write();
	hLike2D_13->Write();
	hLike2D_23->Write();
	hLike3D->Write();
	fout->Close();
      }
    }
    
  }//iexperiment
  
  // if(f_sysparm) {
  //   if(f_sysparm->IsOpen())
  //     f_sysparm->Close();
  //   delete f_sysparm;
  // }

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();  
  
  return;

}
