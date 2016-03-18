//
// ------------------------------------------------------------------------------------------------
// Run fits to toy MC experiments generated on-the-fly and store the toy MC, the best-fit parameter
// values and the best-fit spectra, and the true input parameter values.
//
// Costas Andreopoulos
//
// Fit for nubar mode analysis: different oscillation parameters and probabilities are used to
// oscillate neutrinos and anti-neutrinos. Used in general for the RHC-only mode to search for 
// anti-nue appearance and anti-numu disappearance 
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------
//

#include <vector>

void run_fit_joint_bar_Dm2LocalMin
(
  const char * tag                     = "toyexptfits_Dm2Loc_10",                       /* tag for output file */
  int          run                     = 0,                                   /* run number */
  int          nexperiments            = 2
,                               /* number of toy MC experiment fits */
  
  double       input_mixing12          = -99, //0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing13          = -99, //0.02146, //0.0251,                                /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing23          = -99, //0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
  double       input_deltacp           = -99, //0.,                                /* deltacp, if -999 then randomize x TMath::Pi()!!!!*/
  double       input_dm12sq            = -99, //7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq            = -99, //2.32E-3,                             /* |Dm^2_{23}| */

  double       input_mixing12_bar      = -99, //0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing13_bar      = -99, //0.02146, //0.0251,                    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing23_bar      = -99, //0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
  double       input_deltacp_bar       = -99, //0.,                                /* deltacp, if -999 then randomize x TMath::Pi()!!!!*/
  double       input_dm12sq_bar        = -99, //7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq_bar        = -99, //2.32E-3,                             /* |Dm^2_{23}| */

  const char * config_file             = "analysis_configuration_run1to6_test.xml",        /* analysis configuration file */

  bool         stat_fluctuations       = true,                                /* include statistical fluctuation in toy MC expts? */
  int          tweak_systematics       = false,                                /* generate toy MC with randomized systematics? */

  //const char * fit_param_list          = "nuosc3f_12fix",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */ 
  //const char * fit_param_list          = "nuosc3f_12fix+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_param_list          = "nuosc3fbar_fixnonbar_12fix",
  //const char * fit_param_list          = "nuosc3fbar_fixnonbar_12fix+nuisance",

  const char * fit_option              = "minuit",                            /* fitter option */
  //const char * fit_option              = "minuit+hesse",                            /* fitter option */

  bool         store_toy_expts         = false,                                /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false,                                /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one",                               /* one==kNuOctantOne two==kNuOctantTwo */

  double       dcp_value               = 0.,                                    /* value of deltacp (in units of PI) to assume (when deltacp is not a free parameter) */
  bool         reactor                 = false,                               /* enable this if you want to apply the reactor constraint */
  double       tolerance               = 0.001,                                 /* value for the fitter tolerance */

  char *optfixpar        = "mixing13, dmsq12, deltacp"  // list of the params to be fixed in the fit   
 //TString optfixparams         = "mixing23, dmsq32, dmsq12, mixing13, deltacp"  // list of the params to be fixed in the fit   
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;


  const int Nmin = 3;
  double Dm2LocalMin[Nmin] = {
    0.0025, 0.020, 0.035
    //0.0025
    //0.04
  };

  // const int Nmin = 9;
  // const double Dm2LocalMin[Nmin] = {
  //   0.0025, 0.007, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040
  // };
  

  double Dm2LocalMinChisq;
  double Dm2LocalMinChisqVec[Nmin];
  double Dm2LocalMinVal;
  double Dm2LocalMinValVec[Nmin];
  
  int Dm2LocalMinIndex;
  TH1F *hDm2LocalMinIndex = new TH1F("hDm2LocalMinIndex","hDm2LocalMinIndex",Nmin,0,Nmin);
  


  double bestfit_chisq,bestfit_chisq2;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  TString optfixparams(optfixpar);

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

  bool randomize_mixing12_bar = (TMath::Abs(input_mixing12_bar + 999) < kASmallNum);
  bool randomize_mixing13_bar = (TMath::Abs(input_mixing13_bar + 999) < kASmallNum);
  bool randomize_deltacp_bar  = (TMath::Abs(input_deltacp_bar  + 999) < kASmallNum);
  bool randomize_dm12sq_bar   = (TMath::Abs(input_dm12sq_bar   + 999) < kASmallNum);

  //Set to defaults
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsq2Theta12;

    if(input_mixing13_bar < -1 && input_mixing13_bar > -900) input_mixing13_bar = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_bar < -1 && input_mixing23_bar > -900) input_mixing23_bar = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_bar < -1 && input_mixing12_bar > -900) input_mixing12_bar = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsqTheta12;

    if(input_mixing13_bar < -1 && input_mixing13_bar > -900) input_mixing13_bar = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_bar < -1 && input_mixing23_bar > -900) input_mixing23_bar = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_bar < -1 && input_mixing12_bar > -900) input_mixing12_bar = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1 && input_dm23sq > -900)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1 && input_dm12sq > -900)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90 && input_deltacp > -900) input_deltacp = conv::kGlobFitDeltaCP;

  if(input_dm23sq_bar < -1 && input_dm23sq_bar > -900)  input_dm23sq_bar = conv::kGlobFitDm23sq;
  if(input_dm12sq_bar < -1 && input_dm12sq_bar > -900)  input_dm12sq_bar = conv::kGlobFitDm12sq;
  if(input_deltacp_bar < -90 && input_deltacp_bar > -900) input_deltacp_bar = conv::kGlobFitDeltaCP;


  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Fit parameter list
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = util::NamedParameterList(fit_param_list);
  FitParamList * fparm2 = util::NamedParameterList(fit_param_list);

  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();

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
  if(hyp == kHypOsc3FBar){
    cout << "   - True oscillation parameters bar: " << endl;
    cout << "       > " << util::RawParamAsString(kPhypMixing12Bar) << " = " << ((randomize_mixing12_bar) ? "* randomized *" : input_mixing12_bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing13Bar) << " = " << ((randomize_mixing13_bar) ? "* randomized *" : input_mixing13_bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing23Bar) << " = " << input_mixing23_bar << endl;
    cout << "       > " << util::AsString(kPhypCPVBar    )  << " = " << ((randomize_deltacp_bar) ? "* randomized *" : input_deltacp_bar) << endl;
    cout << "       > " << util::AsString(kPhypDm12sqBar)   << " = " << ((randomize_dm12sq_bar)  ? "* randomized *" : input_dm12sq_bar ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)   << " = " << input_dm23sq_bar << " eV^2/c^4" << endl;
  }
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

  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  //make an assumption for deltacp
  if(hyp==kHypOsc3F && fparm->IsFixed(conv::kPhypCPV)) {
    fparm ->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm1->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm2->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
   }
  else if(hyp==kHypOsc3FBar && fparm->IsFixed(conv::kPhypCPVBar)) {
    fparm ->SetInitValue(conv::kPhypCPVBar    , dcp_value_bar * conv::kpi);
    fparm1->SetInitValue(conv::kPhypCPVBar    , dcp_value_bar * conv::kpi);
    fparm2->SetInitValue(conv::kPhypCPVBar    , dcp_value_bar * conv::kpi);
   }

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

  fitnt->Branch("Nmin",&Nmin,"Nmin/I");
  fitnt->Branch("Dm2LocalMinValVec",&Dm2LocalMinValVec,"Dm2LocalMinValVec[Nmin]/D");
  fitnt->Branch("Dm2LocalMinChisqVec",&Dm2LocalMinChisqVec,"Dm2LocalMinChisqVec[Nmin]/D");
  fitnt->Branch("Dm2LocalMin",&Dm2LocalMin,"Dm2LocalMin[Nmin]/D");

  //fitnt->Branch("Dm2LocalMinIndex",hDm2LocalMinIndex,"Dm2LocalMinIndex/I");
  //fitnt->Branch("Dm2LocalMinVal",hDm2LocalMinVal,"Dm2LocalMinVal/D");

  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");




  ////////////////////////
  //                    //
  //   Fix parameters   //
  //                    //
  ////////////////////////
  
  if(hyp == kHypOsc3F){	
    if(optfixparams.Contains("mixing23")){ 
      fparm->FixParameter(kPhypMixing23);
      fparm1->FixParameter(kPhypMixing23);
      fparm2->FixParameter(kPhypMixing23);
    }
    if(optfixparams.Contains("mixing13")){ 
      fparm->FixParameter(kPhypMixing13);
      fparm1->FixParameter(kPhypMixing13);
      fparm2->FixParameter(kPhypMixing13);
    }
    if(optfixparams.Contains("dmsq32")){ 
      fparm->FixParameter(kPhypDm23sq);
      fparm1->FixParameter(kPhypDm23sq);
      fparm2->FixParameter(kPhypDm23sq);
    }
    if(optfixparams.Contains("dmsq12")){ 
      fparm->FixParameter(kPhypDm12sq);
      fparm1->FixParameter(kPhypDm12sq);
      fparm2->FixParameter(kPhypDm12sq);
    }
    if(optfixparams.Contains("deltacp")){ 
      fparm->FixParameter(kPhypCPV);
      fparm1->FixParameter(kPhypCPV);
      fparm2->FixParameter(kPhypCPV);
    }
  }
  else if(hyp == kHypOsc3FBar){	
    if(optfixparams.Contains("mixing23")){ 
      fparm->FixParameter(kPhypMixing23Bar);
      fparm1->FixParameter(kPhypMixing23Bar);
      fparm2->FixParameter(kPhypMixing23Bar);
    }
    if(optfixparams.Contains("mixing13")){ 
      fparm->FixParameter(kPhypMixing13Bar);
      fparm1->FixParameter(kPhypMixing13Bar);
      fparm2->FixParameter(kPhypMixing13Bar);
    }
    if(optfixparams.Contains("dmsq32")){ 
      fparm->FixParameter(kPhypDm23sqBar);
      fparm1->FixParameter(kPhypDm23sqBar);
      fparm2->FixParameter(kPhypDm23sqBar);
    }
    if(optfixparams.Contains("dmsq12")){ 
      fparm->FixParameter(kPhypDm12sqBar);
      fparm1->FixParameter(kPhypDm12sqBar);
      fparm2->FixParameter(kPhypDm12sqBar);
    }
    if(optfixparams.Contains("deltacp")){ 
      fparm->FixParameter(kPhypCPVBar);
      fparm1->FixParameter(kPhypCPVBar);
      fparm2->FixParameter(kPhypCPVBar);
    }
  }
  




  // Loop over toy-MC experiments
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
    
    cout << endl;
    cout << endl;
    cout << "********************************************" << endl;
    cout << endl;
    cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;
    cout << endl;
    cout << "********************************************" << endl;
    cout << endl;

     // randomize specified osc params for current toy MC expt
     if(hyp == kHypOsc3F){
       if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
       if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
       if (randomize_deltacp ) { input_deltacp  = util::RandomizeCPV();  }
       if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }
     }
     else if(hyp == kHypOsc3FBar){ // randomize nu osc params only if not fixed
       if (randomize_mixing12_bar) { input_mixing12_bar = util::RandomizeMixing12Bar(); }
       if (randomize_mixing13_bar) { input_mixing13_bar = util::RandomizeMixing13Bar(); }
       if (randomize_deltacp_bar ) { input_deltacp_bar  = util::RandomizeCPVBar();  }
       if (randomize_dm12sq_bar  ) { input_dm12sq_bar   = util::RandomizeDm12sqBar();   }
     }

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
     else if(hyp == kHypOsc3FBar) {
       fparm1->Update(kPhypMixing12, input_mixing12);
       fparm1->Update(kPhypMixing13, input_mixing13);
       fparm1->Update(kPhypMixing23, input_mixing23);
       fparm1->Update(kPhypCPV    ,  input_deltacp*conv::kpi);
       fparm1->Update(kPhypDm12sq,   input_dm12sq);
       fparm1->Update(kPhypDm23sq,   input_dm23sq);

       fparm1->Update(kPhypMixing12Bar, input_mixing12_bar);
       fparm1->Update(kPhypMixing13Bar, input_mixing13_bar);
       fparm1->Update(kPhypMixing23Bar, input_mixing23_bar);
       fparm1->Update(kPhypCPVBar    ,  input_deltacp_bar*conv::kpi);
       fparm1->Update(kPhypDm12sqBar,   input_dm12sq_bar);
       fparm1->Update(kPhypDm23sqBar,   input_dm23sq_bar);
       /////////////////////////////////////////
       fparm2->Update(kPhypMixing12, input_mixing12);
       fparm2->Update(kPhypMixing13, input_mixing13);
       fparm2->Update(kPhypMixing23, input_mixing23);
       fparm2->Update(kPhypCPV    ,  input_deltacp*conv::kpi);
       fparm2->Update(kPhypDm12sq,   input_dm12sq);
       fparm2->Update(kPhypDm23sq,   input_dm23sq);

       fparm2->Update(kPhypMixing12Bar, input_mixing12_bar);
       fparm2->Update(kPhypMixing13Bar, input_mixing13_bar);
       fparm2->Update(kPhypMixing23Bar, input_mixing23_bar);
       fparm2->Update(kPhypCPVBar    ,  input_deltacp_bar*conv::kpi);
       fparm2->Update(kPhypDm12sqBar,   input_dm12sq_bar);
       fparm2->Update(kPhypDm23sqBar,   input_dm23sq_bar);
     }
     else if(hyp == kHypNuMuDisapp2F) {
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

     if(hyp == kHypOsc3F) {
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
     }
     else if(hyp == kHypOsc3FBar){
       if(!(fparm1->IsFixed(conv::kPhypDm23sqBar))) {   
	 fparm1->SetInitValue(conv::kPhypMixing23Bar,0.4);
	 fparm1->Update(conv::kPhypMixing23Bar,0.4);
	 fparm2->SetInitValue(conv::kPhypMixing23Bar,0.6);
	 fparm2->Update(conv::kPhypMixing23Bar,0.6);
       }
       //if reactor constraint enabled, set nuisance
       if(reactor) {
	 fparm1->SetNuisance(conv::kPhypMixing13Bar,true);
	 fparm2->SetNuisance(conv::kPhypMixing13Bar,true);
       }
     }

     //reset systematic variations, otherwise fit is initialized at
     //best fit values of systematics from previous fit
     skexp->ResetSystematicVariations();

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

     // Get fitter
     FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
     FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);

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


     // Loop over local minima
     for(int imin=0;imin<Nmin;imin++){

       cout << endl;
       cout << "/////////////////////////////////" << endl;
       cout << endl;
       cout << "Start Fit from " << Dm2LocalMin[imin] << endl;
       cout << endl;
       
       fparm1->SetInitValue(conv::kPhypDm23sqBar,Dm2LocalMin[imin] );
       fparm1->Update(conv::kPhypDm23sqBar,Dm2LocalMin[imin] );
       fparm2->SetInitValue(conv::kPhypDm23sqBar,Dm2LocalMin[imin] );
       fparm2->Update(conv::kPhypDm23sqBar,Dm2LocalMin[imin] );
       
       skexp->ResetSystematicVariations();
       fparm1->ResetCurrToInitValues(); // reset osc params
 	    
       cout << endl;
       cout << endl;
       cout << "First Octant" << endl;
       cout << endl;
       cout << endl;

       // Print-out true params
       fparm1->ShowOnlyCurrValuesInNextPrint();
       cout << "- True oscillation params (below) yield chisq = " << true_chisq
	    << *fparm1;      

       //skexp->PrintSystematicVariations();
       fitter->SetTol(tolerance);
       fitter->DoTheFit();
       fitstatus1 = fitter->FitStatusCode("migrad_min");
       fitstatus_err1 = fitter->FitStatusCode("hesse_err");     
       fittime1 = fitter->FitTime();
       
       double bestfit_chisq  = fparm1->GetBestFitStatistic();
       chisq1 = bestfit_chisq;
       double bestfit_chisq2 = 0;
      
       //Do the fit again at the mirror point
       if( Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta &&
	   (!(fparm2->IsFixed(conv::kPhypMixing23)) && hyp == kHypOsc3F) || 
	   (!(fparm2->IsFixed(conv::kPhypMixing23Bar)) && hyp == kHypOsc3FBar)
	   ){
	 
	 cout << endl;
	 cout << endl;
	 cout << "Second Octant" << endl;
	 cout << endl;
	 cout << endl;
	 
	 skexp->ResetSystematicVariations();
	 
	 fparm2->ResetCurrToInitValues(); // reset osc params
	 
	 // Chisq for true (input) params
	 double true_chisq = fit_inputs2->Statistic();
	 // Print-out true params
	 fparm2->ShowOnlyCurrValuesInNextPrint();
	 cout << "- True oscillation params (below) yield chisq = " << true_chisq
	      << *fparm2;
	 
	 fitter2->SetTol(tolerance);
	 fitter2->DoTheFit();
	 fitstatus2 = fitter->FitStatusCode("migrad_min");
	 fitstatus_err2 = fitter->FitStatusCode("hesse_err");     
	 bestfit_chisq2 =  fparm2->GetBestFitStatistic();
	 chisq2 = bestfit_chisq2;

       } // octant
     
       cout << "chisq1 = " << chisq1 << " - chisq2 = " << chisq2 << endl;
       if(chisq1<chisq2){ 
	 Dm2LocalMinChisqVec[imin] = chisq1;
	 Dm2LocalMinValVec[imin] = fparm1->BestFitValue(conv::kPhypDm23sqBar);
       }
       else{
	 Dm2LocalMinChisqVec[imin] = chisq2;
      	 Dm2LocalMinValVec[imin] = fparm2->BestFitValue(conv::kPhypDm23sqBar);
       }   
  
       //cout << "Dm2LocalMinChisqVec[" << imin << "] = " << Dm2LocalMinChisqVec[imin] << endl;
       //cout << "Dm2LocalMinValVec[" << imin << "] = " << Dm2LocalMinValVec[imin] << endl;
       //getchar();
     
     } // end loop over Dm2

     cout << endl;
     cout << endl;
     cout << endl;
     cout << endl;
     cout << endl;
     cout << endl;

     double chisqval = 99999.;

     for(int imin=0;imin<Nmin;imin++){
       
       // if(Dm2LocalMinChisqVec[imin] < chisqval ){
       // 	 chisqval = Dm2LocalMinChisqVec[imin];
       // 	 Dm2LocalMinChisq = chisqval;
       // 	 Dm2LocalMinVal = Dm2LocalMinValVec[imin];
       // 	 Dm2LocalMinIndex = imin;
       // }

       cout << "Dm2LocalMinChisq[" << imin << "] = " << Dm2LocalMinChisqVec[imin] << endl;
       cout << "Dm2LocalMinValVec[" << imin << "] = " << Dm2LocalMinValVec[imin] << endl;
     }

     // cout << endl;
     // cout << "Dm2LocalMinIndex = " << Dm2LocalMinIndex << endl;
     // cout << "Dm2LocalMinVal = " << Dm2LocalMinVal << endl;
     // hDm2LocalMinIndex->Fill(Dm2LocalMinIndex);
 



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
  } // end toy experiments
  


  // cout << hDm2LocalMinIndex->Integral() << endl;
  // hDm2LocalMinIndex->Draw("hist");
  

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  //hDm2LocalMinIndex->Write();
  f.Close();
}
