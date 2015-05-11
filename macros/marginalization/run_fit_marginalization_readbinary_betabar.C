
// ------------------------------------------------------------------------------------------------           
//
// Macro to perform fit marginalizing the systematic and oscillation parameters
// Run toys for marginalization for all or sub-range of the free parameter and fill the Chisq.
// In the end the minimum of the marginalized Chisq distribution is found.
// Toys for different sub-ranges of the free paramters can be run in parallel jobs.
// The other possibility is to produce toys with build_marg_toys.C and store in binary files.
// Then this macro reads the toys from binary files, fill the marginalized Chisq and minimize. 
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------




//
// TODO:
// - set correct oscillation parameters for reactor constraint
// 

using namespace std;
using namespace valor;
using namespace valor::conv;

void run_fit_marginalization_readbinary_betabar
(

 ////////////////////////// GENERAL PARAMETERS /////////////////////

 const char * tag                      = "fitmarg_betabar",    /* tag for output file */

 int          first_experiment         = 0, //4404, //4739,                       /* */
 int          last_experiment          = 0, //4404, //4739,                       /* */

 //const char * xmlfile = "toymc/toymc_nuebarapp_fixosc_betabar_1_rethrowunphysical_4E20POT-0.xml",
 //const char * xmlfile = "toymc/toymc_nuebarapp_fixosc_betabar_0_rethrowunphysical_4E20POT-0.xml",
 const char * xmlfile = "toymc/toymc_nuebarapp_fixosc_betabar_1_rethrowunphysical_7p8E21POT-0.xml",
 
 bool         dotoys                   = false,

 // Default
 double       input_mixing12           = -99,           
 double       input_mixing13           = -99, 
 double       input_mixing23           = -99, 
 double       input_deltacp            = -99, 
 double       input_dm12sq             = -99,
 double       input_dm23sq             = -99,
 double       input_mixing12_bar       = -99, 
 double       input_mixing13_bar       = -99,
 double       input_mixing23_bar       = -99,
 double       input_deltacp_bar        = -99,
 double       input_dm12sq_bar         = -99,
 double       input_dm23sq_bar         = -99,
 double       input_beta_bar           = -99,

 //const char * config_file              = "analysis_configuration_run1to6_test_app.xml",
 const char * config_file              = "analysis_configuration_run1to6_test_app_fullstat.xml",
 
 bool         stat_fluctuations        = false,    /* include statistical fluctuation in toy MC expts? (not for marginalization) */
 int          tweak_systematics        = 0,       /* generate toy (both fakedata sets and marginalization) MC with randomized systematics? */
 
 //const char * fit_param_list       = "nuosc3fbar_normbar_fixall+nuisance",
 const char * fit_param_list       = "nuosc3fbar_normbar_fixosc+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_normbar_fitdcpbar+nuisance",

 const char * dmsqscheme               = "abs23",   /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamchoice           = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * hierarchy                = "normal",     /* normal or inverted */
 const char * mixparamoctant           = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

 ////////////////////// MARGINALIZATION PARAMETERS ///////////////////

 //const char * infilename           = "../../inputs/marginalization/toys/nuebarapp_1D_T2Kreact_betabar_fixosc_10kToys/BetaBar_4p0E20POT",  // T2K+reactor (betabar=0,1)
 //const char * infilename           = "../../inputs/marginalization/toys/nuebarapp_1D_T2Kreact_betabar_fixosc_10kToys_grid/build_marg",  // T2K+reactor (betabar grid)
 
 const char * infilename           = "../../inputs/marginalization/toys/nuebarapp_1D_T2KFullStat_T2Kreact_betabar_fixosc_10kToys/BetaBar_7p8E21POT",  // T2K+reactor T2K Full Stat (betabar=0,1) 


 const int    nfiles               = 1, //401, 1        // # of files to read

 const int    NtoysMarg            = 10000,   // # toys marginalization 
 
 const char * freeparname0        = "betaBar",  // free parameter
 const double freeparmin0         = 0.,         // minimum value of free parameter
 const double freeparmax0         = 1, // 10., 1         // maxinum value of free parameter
 const int    NptScan0            = 2, // 801, 2          // # pts of free parameter to scan
 
 const bool  doMargMH             = true,    // marginalize Mass Hierarchy
 
 const TString TypeError          = "",  // "HESSE-MINOS" (contained by the string)  
 const bool DrawFullChisq         = false,

 const string teststat            = "rate-only"   // "rate-only" or "rate+shape"
)
{

  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  // Set the second parameter (1D fit)
  const char * freeparname1        = "";  // free parameter
  const double freeparmin1         = 0.; // minimum value of free parameter
  const double freeparmax1         = 0.; // maxinum value of free parameter
  const int    NptScan1            = 1; //201, //21,      // # pts of free parameter to scan
  //

  TString opt(fit_param_list); 

  // Analysis configuration
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file));
  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile));

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
  
  if (!dotoys && gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
  {
    xmlfile2 = xmlfile;
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
      cerr<<"Data input file not found: "<<xmlfile<<endl;
      exit(-1);
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
  bool randomize_mixing23 = (TMath::Abs(input_mixing23 + 999) < kASmallNum);
  bool randomize_deltacp  = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq   = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);
  bool randomize_dm23sq   = (TMath::Abs(input_dm23sq   + 999) < kASmallNum);

  bool randomize_mixing12_bar = (TMath::Abs(input_mixing12_bar + 999) < kASmallNum);
  bool randomize_mixing13_bar = (TMath::Abs(input_mixing13_bar + 999) < kASmallNum);
  bool randomize_mixing23_bar = (TMath::Abs(input_mixing23_bar + 999) < kASmallNum);
  bool randomize_deltacp_bar  = (TMath::Abs(input_deltacp_bar  + 999) < kASmallNum);
  bool randomize_dm12sq_bar   = (TMath::Abs(input_dm12sq_bar   + 999) < kASmallNum);
  bool randomize_dm23sq_bar   = (TMath::Abs(input_dm23sq_bar   + 999) < kASmallNum);

  bool randomize_beta_bar   = (TMath::Abs(input_beta_bar   + 999) < kASmallNum);
    
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

  if(input_beta_bar < -90 && input_beta_bar > -900) input_beta_bar = conv::kGlobFitBetaBar;


  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());
 
  // Fit param list for fakedata sets
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); 
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used for fakedata sets: " << util::AsString(hyp) << endl;

  // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
  if(hyp == kHypOsc3F) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
  }
  else if(hyp == kHypOsc3FBar) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
    
    fparm->SetInitValue(kPhypMixing12Bar, input_mixing12_bar);
    fparm->SetInitValue(kPhypMixing13Bar, input_mixing13_bar);
    fparm->SetInitValue(kPhypMixing23Bar, input_mixing23_bar);
    fparm->SetInitValue(kPhypCPVBar    ,  input_deltacp_bar);
    fparm->SetInitValue(kPhypDm12sqBar,   input_dm12sq_bar);
    fparm->SetInitValue(kPhypDm23sqBar,   input_dm23sq_bar);
  }
  else if(hyp == kHypOsc3FBarNormBar) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
    
    fparm->SetInitValue(kPhypMixing12Bar, input_mixing12_bar);
    fparm->SetInitValue(kPhypMixing13Bar, input_mixing13_bar);
    fparm->SetInitValue(kPhypMixing23Bar, input_mixing23_bar);
    fparm->SetInitValue(kPhypCPVBar    ,  input_deltacp_bar);
    fparm->SetInitValue(kPhypDm12sqBar,   input_dm12sq_bar);
    fparm->SetInitValue(kPhypDm23sqBar,   input_dm23sq_bar);

    fparm->SetInitValue(kPhypBetaBar,     input_beta_bar);
  }

  else if(hyp == kHypNuMuDisapp2F) {
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
  }

   
  int nexperiments = last_experiment - first_experiment + 1; 
  
  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << endl;
  cout << "- Reading analysis configuration from: " << config << endl;

  if(dotoys){
    cout << "- Toy-MC generation" << endl;
    cout << "   - Number of toy experiments: " << nexperiments << endl; 
    cout << "   - Randomize systematics : "      << ((tweak_systematics) ?      "YES " : "NO") << endl;
    cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ?      "YES " : "NO") << endl;
    cout << "   - True oscillation parameters: " << endl;
    cout << "       > " << util::AsString(kPhypMixing12) << " = " << ((randomize_mixing12) ? "* randomized *" : input_mixing12) << endl;
    cout << "       > " << util::AsString(kPhypMixing13) << " = " << ((randomize_mixing13) ? "* randomized *" : input_mixing13) << endl;
    cout << "       > " << util::AsString(kPhypMixing23) << " = " << ((randomize_mixing23) ? "* randomized *" : input_mixing23) << endl;
    cout << "       > " << util::AsString(kPhypCPV    )  << " = " << ((randomize_deltacp) ? "* randomized *" : input_deltacp) << endl;
    cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << ((randomize_dm23sq)  ? "* randomized *" : input_dm23sq ) << " eV^2/c^4" << endl;
    if(hyp == kHypOsc3FBar){
      cout << "   - True oscillation parameters bar: " << endl;
      cout << "       > " << util::RawParamAsString(kPhypMixing12Bar) << " = " << ((randomize_mixing12_bar) ? "* randomized *" : input_mixing12_bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing13Bar) << " = " << ((randomize_mixing13_bar) ? "* randomized *" : input_mixing13_bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing23Bar) << " = " << ((randomize_mixing23_bar) ? "* randomized *" : input_mixing23_bar) << endl;
      cout << "       > " << util::AsString(kPhypCPVBar    )  << " = " << ((randomize_deltacp_bar) ? "* randomized *" : input_deltacp_bar) << endl;
      cout << "       > " << util::AsString(kPhypDm12sqBar)   << " = " << ((randomize_dm12sq_bar)  ? "* randomized *" : input_dm12sq_bar ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypDm23sqBar)   << " = " << ((randomize_dm23sq_bar)  ? "* randomized *" : input_dm23sq_bar ) << " eV^2/c^4" << endl;
    }
    else if(hyp == kHypOsc3FBarNormBar){
      cout << "   - True oscillation parameters bar: " << endl;
      cout << "       > " << util::RawParamAsString(kPhypMixing12Bar) << " = " << ((randomize_mixing12_bar) ? "* randomized *" : input_mixing12_bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing13Bar) << " = " << ((randomize_mixing13_bar) ? "* randomized *" : input_mixing13_bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing23Bar) << " = " << ((randomize_mixing23_bar) ? "* randomized *" : input_mixing23_bar) << endl;
      cout << "       > " << util::AsString(kPhypCPVBar    )  << " = " << ((randomize_deltacp_bar) ? "* randomized *" : input_deltacp_bar) << endl;
      cout << "       > " << util::AsString(kPhypDm12sqBar)   << " = " << ((randomize_dm12sq_bar)  ? "* randomized *" : input_dm12sq_bar ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypDm23sqBar)   << " = " << ((randomize_dm23sq_bar)  ? "* randomized *" : input_dm23sq_bar ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypBetaBar)     << " = " << ((randomize_beta_bar)  ? "* randomized *" : input_beta_bar ) << endl;
    }
  }
  else{
      cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiments from: " << xmlfile2 << endl; 
  cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl; 
  }
  
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "- Outputs" << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;

  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each toy-MC.
  TObjArray * dat_array = new TObjArray; //only used to initialise branch
  TObjArray * pdf_array = new TObjArray; //only used to initialise branch
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);
  
  // Create file to store results
  string oscstr = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);  
  TFile f(Form("%s-%s_%d_%d.root", tag, oscstr.c_str(), first_experiment, last_experiment), "recreate");

  double truepar0=-999;
  double truepar1=-999;    
  int nfits=-999;
  int nchisqcalls=-999;
  int niterations=-999;
  int fitstatus=-999;
  double chisqmin=conv::kALargeNum;
  double chisqtrue0=conv::kALargeNum;
  double chisqtrue1=conv::kALargeNum;
  double deltachisqtrue=conv::kALargeNum;
  double deltachisq01=conv::kALargeNum;
  double deltachisq0BF=conv::kALargeNum;
  TH2F *hchisq2D_out;
  TH1F *hchisq1D_out;
  TH2F *hstep2D_out;
  TH1F *hstep1D_out;
  double bestfit0=-999;
  double bestfit1=-999;
  double bestfiterrNeg0=-999;
  double bestfiterrNeg1=-999;
  double bestfiterrPos0=-999;
  double bestfiterrPos1=-999;
  double bestfiterrSym0=-999;
  double bestfiterrSym1=-999;
  double bestfitFullChisq0=-999;
  double bestfitFullChisq1=-999;
  double chisqminFullChisq=conv::kALargeNum;
  int ToyID=-999;
  TH1D *hErecoToy;
  
  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("truepar0", &truepar0 );
  fitnt->Branch("truepar1", &truepar1 );
  fitnt->Branch("bestfit0", &bestfit0 );
  fitnt->Branch("bestfit1", &bestfit1 );
  fitnt->Branch("hchisq2D", &hchisq2D_out );
  fitnt->Branch("hchisq1D", &hchisq1D_out );
  fitnt->Branch("hstep2D", &hstep2D_out );
  fitnt->Branch("hstep1D", &hstep1D_out );
  fitnt->Branch("fitstatus",&fitstatus,"fitstatus/I");
  fitnt->Branch("chisqmin",&chisqmin,"chisqmin/D");
  fitnt->Branch("chisqtrue0",&chisqtrue0,"chisqtrue0/D");
  fitnt->Branch("chisqtrue1",&chisqtrue1,"chisqtrue1/D");
  fitnt->Branch("deltachisqtrue",&deltachisqtrue,"deltachisqtrue/D");
  fitnt->Branch("nfits",&nfits,"nfits/I");
  fitnt->Branch("nchisqcalls",&nchisqcalls,"nchisqcalls/I");
  fitnt->Branch("niterations",&niterations,"niterations/I");
  fitnt->Branch("bestfiterrNeg0",&bestfiterrNeg0);
  fitnt->Branch("bestfiterrNeg1",&bestfiterrNeg1);
  fitnt->Branch("bestfiterrPos0",&bestfiterrPos0);
  fitnt->Branch("bestfiterrPos1",&bestfiterrPos1);
  fitnt->Branch("bestfiterrSym0",&bestfiterrSym0);
  fitnt->Branch("bestfiterrSym1",&bestfiterrSym1);
  fitnt->Branch("bestfitFullChisq0", &bestfitFullChisq0 );
  fitnt->Branch("bestfitFullChisq1", &bestfitFullChisq1 );
  fitnt->Branch("chisqminFullChisq",&chisqminFullChisq);
  fitnt->Branch("ToyID",&ToyID);
  fitnt->Branch("hErecoToy",&hErecoToy);
  
  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  // Loop over toy-MC experiments
  //for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {
    cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;
    
    ToyID = iexperiment; // Store the ToyID in the output tree

    // randomize specified osc params for current toy MC expt
    if(hyp == kHypOsc3F){
      if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
      if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
      if (randomize_mixing23) { input_mixing23 = util::RandomizeMixing23(); }
      if (randomize_deltacp ) { input_deltacp  = util::RandomizeCPV();  }
      if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }
      if (randomize_dm23sq  ) { input_dm23sq   = util::RandomizeDm23sq();   }
    }
    else if(hyp == kHypOsc3FBar){ // randomize nu osc params only if not fixed
      if (randomize_mixing12_bar) { input_mixing12_bar = util::RandomizeMixing12Bar(); }
      if (randomize_mixing13_bar) { input_mixing13_bar = util::RandomizeMixing13Bar(); }
      if (randomize_mixing23_bar) { input_mixing23_bar = util::RandomizeMixing23Bar(); }
      if (randomize_deltacp_bar ) { input_deltacp_bar  = util::RandomizeCPVBar();  }
      if (randomize_dm12sq_bar  ) { input_dm12sq_bar   = util::RandomizeDm12sqBar();   }
      if (randomize_dm23sq_bar  ) { input_dm23sq_bar   = util::RandomizeDm23sqBar();   }
    }
    else if(hyp == kHypOsc3FBarNormBar){ // randomize nu osc params only if not fixed
      if (randomize_mixing12_bar) { input_mixing12_bar = util::RandomizeMixing12Bar(); }
      if (randomize_mixing13_bar) { input_mixing13_bar = util::RandomizeMixing13Bar(); }
      if (randomize_mixing23_bar) { input_mixing23_bar = util::RandomizeMixing23Bar(); }
      if (randomize_deltacp_bar ) { input_deltacp_bar  = util::RandomizeCPVBar();  }
      if (randomize_dm12sq_bar  ) { input_dm12sq_bar   = util::RandomizeDm12sqBar();   }
      if (randomize_dm23sq_bar  ) { input_dm23sq_bar   = util::RandomizeDm23sqBar();   }
      if (randomize_beta_bar  )   { input_beta_bar     = util::RandomizeBetaBar();   }
    }

    // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
    if(hyp == kHypOsc3F) {
      fparm->SetInitValue(kPhypMixing12, input_mixing12);
      fparm->SetInitValue(kPhypMixing13, input_mixing13);
      fparm->SetInitValue(kPhypMixing23, input_mixing23);
      fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
      fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
      fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
    }
    else if(hyp == kHypOsc3FBar) {
      fparm->SetInitValue(kPhypMixing12, input_mixing12);
      fparm->SetInitValue(kPhypMixing13, input_mixing13);
      fparm->SetInitValue(kPhypMixing23, input_mixing23);
      fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
      fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
      fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);

      fparm->SetInitValue(kPhypMixing12Bar, input_mixing12_bar);
      fparm->SetInitValue(kPhypMixing13Bar, input_mixing13_bar);
      fparm->SetInitValue(kPhypMixing23Bar, input_mixing23_bar);
      fparm->SetInitValue(kPhypCPVBar    ,  input_deltacp_bar);
      fparm->SetInitValue(kPhypDm12sqBar,   input_dm12sq_bar);
      fparm->SetInitValue(kPhypDm23sqBar,   input_dm23sq_bar);
    }
    else if(hyp == kHypOsc3FBarNormBar) {
      fparm->SetInitValue(kPhypMixing12, input_mixing12);
      fparm->SetInitValue(kPhypMixing13, input_mixing13);
      fparm->SetInitValue(kPhypMixing23, input_mixing23);
      fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
      fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
      fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);

      fparm->SetInitValue(kPhypMixing12Bar, input_mixing12_bar);
      fparm->SetInitValue(kPhypMixing13Bar, input_mixing13_bar);
      fparm->SetInitValue(kPhypMixing23Bar, input_mixing23_bar);
      fparm->SetInitValue(kPhypCPVBar    ,  input_deltacp_bar);
      fparm->SetInitValue(kPhypDm12sqBar,   input_dm12sq_bar);
      fparm->SetInitValue(kPhypDm23sqBar,   input_dm23sq_bar);

      fparm->SetInitValue(kPhypBetaBar,     input_beta_bar);
    }
    else if(hyp == kHypNuMuDisapp2F) {
      fparm->SetInitValue(kPhypMixing23, input_mixing23);
      fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
    }
       
    // Reset all output spectra
    dat_array->Clear();
    pdf_array->Clear();
 
    // Reset input true params and add true oscillation params
    iparm->Clear();
    PhysParamList true_oscillation_params = fparm->AsPhysParamList("current");
    iparm->Append(true_oscillation_params);

    // get toy MC with/without systematics
    
    SuperkObservation * skobs_fakedata = new SuperkObservation;
    bool using_xml_file = true;
    if(!dotoys && using_xml_file){
      skobs_fakedata->Read(xmlfile2.Data(),iexperiment,using_xml_file);
    }
    else{
      bool applied_negevents_condition = true;
      bool unphysical_tweak = true;
      bool rethrow_if_parameters_unphysical = true;// both marginalization and fakedata sets
      while(applied_negevents_condition || unphysical_tweak) {
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
          
          PhysParamList true_nuisance_params = skexp->Tweaks();
          if (rethrow_if_parameters_unphysical) {
            unphysical_tweak = !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
            if (unphysical_tweak) {cout<<"Parameters are unphysical - rethrowing"<<std::endl;}
          }
          else {
            unphysical_tweak = false;
          }
          
          // generate mock data using the current PDF (incl. systematic variations) 
          // + statistical fluctuations if requested (Poisson fluctuation in each bin)
          cout << "- Generating toy-MC..." << endl;
          //note that fparm does not contain the tweaks, so use only the PhysHyp here
          applied_negevents_condition = skobs_fakedata->CreateMockData(*skexp, PhysHypothesis(*fparm), stat_fluctuations);
	  cout << "applied_negevents_condition = " << applied_negevents_condition << endl;

          
          if(!(applied_negevents_condition || unphysical_tweak)) {
            // keep track of all nuisance params set so as to add in output tree
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
          applied_negevents_condition = skobs_fakedata->CreateMockData(*skexp, PhysHypothesis(*fparm), stat_fluctuations); 
	  cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
	  
	  unphysical_tweak = false; // Otherwise it doesn't exit the loop
        }
      }//applied_negevents_condition
    }
        
    cout << *fparm << endl;

    /////////////////// SET MARGINALIZATION FIT ////////////////
    
    // Reset all the variables
    truepar0=-999;
    truepar1=-999;    
    nfits=-999;
    nchisqcalls=-999;
    niterations=-999;
    fitstatus=-999;
    chisqmin=conv::kALargeNum;
    chisqtrue0=conv::kALargeNum;
    chisqtrue1=conv::kALargeNum;
    deltachisqtrue=conv::kALargeNum;
    bestfit0=-999;
    bestfit1=-999;
    bestfiterrNeg0=-999;
    bestfiterrNeg1=-999;
    bestfiterrPos0=-999;
    bestfiterrPos1=-999;
    bestfiterrSym0=-999;
    bestfiterrSym1=-999;  
    bestfitFullChisq0=-999;
    bestfitFullChisq1=-999;
    chisqminFullChisq=conv::kALargeNum;

    //

    PhysParamId_t freepar0;
    string sfreeparname0(freeparname0);  
    // NuBar mode
    if(sfreeparname0=="mixing23Bar")      freepar0 = conv::kPhypMixing23Bar;
    else if(sfreeparname0=="mixing13Bar") freepar0 = conv::kPhypMixing13Bar;
    else if(sfreeparname0=="dmsq23Bar")   freepar0 = conv::kPhypDm23sqBar;
    else if(sfreeparname0=="dmsq12Bar")   freepar0 = conv::kPhypDm12sqBar;
    else if(sfreeparname0=="mixing12Bar") freepar0 = conv::kPhypMixing12Bar;
    else if(sfreeparname0=="deltacpBar")  freepar0 = conv::kPhypCPVBar;
    else if(sfreeparname0=="betaBar")     freepar0 = conv::kPhypBetaBar;
    // Nu mode
    else if(sfreeparname0=="mixing23")    freepar0 = conv::kPhypMixing23;
    else if(sfreeparname0=="mixing13")    freepar0 = conv::kPhypMixing13;
    else if(sfreeparname0=="dmsq23")      freepar0 = conv::kPhypDm23sq;
    else if(sfreeparname0=="dmsq12")      freepar0 = conv::kPhypDm12sq;
    else if(sfreeparname0=="mixing12")    freepar0 = conv::kPhypMixing12;
    else if(sfreeparname0=="deltacp")     freepar0 = conv::kPhypCPV;
    //
    else{
      freepar0 = conv::kPhypUndefined;
      cout << "Param 0 is not defined" << endl;
      exit(1);
    }
    
    PhysParamId_t freepar1 = conv::kPhypUndefined;
    string sfreeparname1(freeparname1);  
    // NuBar mode
    if(sfreeparname1=="mixing23Bar")      freepar1 = conv::kPhypMixing23Bar;
    else if(sfreeparname1=="mixing13Bar") freepar1 = conv::kPhypMixing13Bar;
    else if(sfreeparname1=="dmsq23Bar")   freepar1 = conv::kPhypDm23sqBar;
    else if(sfreeparname1=="dmsq12Bar")   freepar1 = conv::kPhypDm12sqBar;
    else if(sfreeparname1=="mixing12Bar") freepar1 = conv::kPhypMixing12Bar;
    else if(sfreeparname1=="deltacpBar")  freepar1 = conv::kPhypCPVBar;
    else if(sfreeparname1=="betaBar")     freepar1 = conv::kPhypBetaBar;
    // Nu mode
    else if(sfreeparname1=="mixing23")    freepar1 = conv::kPhypMixing23;
    else if(sfreeparname1=="mixing13")    freepar1 = conv::kPhypMixing13;
    else if(sfreeparname1=="dmsq23")      freepar1 = conv::kPhypDm23sq;
    else if(sfreeparname1=="dmsq12")      freepar1 = conv::kPhypDm12sq;
    else if(sfreeparname1=="mixing12")    freepar1 = conv::kPhypMixing12;
    else if(sfreeparname1=="deltacp")     freepar1 = conv::kPhypCPV;
    //
    else{
      freepar1 = conv::kPhypUndefined;
      cout << "1D fit" << endl;
    }
    
    
    // Set free param, ntoys, # pt to scan, min par, max par (if min and max not passed, use default ones for that par)
    
    // free parameter 
    PhysParamId_t parfree[2];
    parfree[0] = freepar0;
    parfree[1] = freepar1;
    
    // # of toys for marginalization
    int ntoys = NtoysMarg;
    
    // # of free parameter pts to scan 
    int nbinspar[2];
    nbinspar[0] = NptScan0;
    nbinspar[1] = NptScan1;
    
    // free parameter range
    double minpar[2];
    double maxpar[2];
    minpar[0] = freeparmin0;
    maxpar[0] = freeparmax0;
    minpar[1] = freeparmin1;
    maxpar[1] = freeparmax1;
    
    // free parameter sub-range
    
    const int    first_bin_freepar0   = -999; // first freepar of the sub-range
    const int    last_bin_freepar0    = -999; // last freepar of the sub-range
    const int    first_bin_freepar1   = -999; // first freepar of the sub-range
    const int    last_bin_freepar1    = -999; // last freepar of the sub-range
 
    int subminbin[2];
    int submaxbin[2];
    subminbin[0] = first_bin_freepar0; submaxbin[0] = last_bin_freepar0;
    subminbin[1] = first_bin_freepar1; submaxbin[1] = last_bin_freepar1;
    
    // Implement the constructor giving the following input parameters:
    FitMarginalization *fitmarg = new FitMarginalization();   
    fitmarg->SetPhysParFree(parfree,ntoys,nbinspar,minpar,maxpar);
    fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
    fitmarg->SetTestStatistic(teststat);     // Set Test Statistic ("rate-only" or "rate+shape")
    fitmarg->Init(skexp); // Before need SetPhysParFree() to initialize the arrays
    
    
    // Set Mass Hierarchy marginalization
    fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
    
    // Set fit parameter list
    fitmarg->SetFitParamList(fit_param_list);
    fitmarg->PrintPhysHypothesis();
   

    // Build the map fileID - sub-range
    
    fitmarg->SetMaxIterationSingleFit(2);  // Max # of iterations per fit
    fitmarg->SetMaxFitNumber(3); // Max # of fits max # of iterations

    fitmarg->SetTagFileName(infilename);
    fitmarg->SetFileID(nfiles);
    
    
    // Fit 1

    cout << endl;
    cout << "/////////////// Fit 1 (all range) ////////////////" << endl;

    fitmarg->SetMinimIntervalRead(freeparmin0,freeparmax0,freeparmin1,freeparmax1,true);
    fitmarg->SetStartBin(0.,1.); 
    fitmarg->SetNumStartBin(1,1);
    fitmarg->SetGridStartRange(freeparmin0,freeparmax0,freeparmin1,freeparmax1);
    fitmarg->SetExtractStartPt("GridFix");
    fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

    double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
    double bestfit0_1 = fitmarg->GetParamBestFit(0);
    double bestfit1_1 = -999.;
    if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
    int nfits_1 = fitmarg->GetNumFit();
    int fitstatus_1 = fitmarg->GetFitStatus();
    int niterations_1 = fitmarg->GetNumIterations();
    if(chisqmin_1<chisqmin){
      chisqmin = chisqmin_1;
      bestfit0 = bestfit0_1; 
      bestfit1 = bestfit1_1; 
      nfits = nfits_1;    
      fitstatus = fitstatus_1;
      niterations = niterations_1;
    }


    nchisqcalls = fitmarg->GetNCallsCalcChisqMargPt();
    fitmarg->BuildDeltaChisqMarg();

    // Set best-fit parameters in FitMarginalization object
    // Otherwise error computation is done with respect to a wrong minimum
    fitmarg->SetChisqBestFit(chisqmin);
    fitmarg->SetParamBestFit(0,bestfit0);
    if(fitmarg->IsFit2D()) fitmarg->SetParamBestFit(1,bestfit1);

    //
    // Get fit errors
    //
    // Since the chisq was minimized reading binary files at each step you 
    // must adopt the same approach for error computations as well.
    //
    
    if(TypeError.Contains("MINOS")) fitmarg->CalcBestFitErrors(skobs_fakedata,"MINOS",true);
    if(TypeError.Contains("HESSE")) fitmarg->CalcBestFitErrors(skobs_fakedata,"HESSE",true);

    bestfiterrNeg0 = fitmarg->GetParamBestFitErrNeg(0);
    bestfiterrNeg1 = fitmarg->GetParamBestFitErrNeg(1);
    bestfiterrPos0 = fitmarg->GetParamBestFitErrPos(0);
    bestfiterrPos1 = fitmarg->GetParamBestFitErrPos(1);
    
    bestfiterrSym0 = fitmarg->GetParamBestFitErrSym(0);
    bestfiterrSym1 = fitmarg->GetParamBestFitErrSym(1);


    // Get true chisq 

    truepar0 = fparm->InitValue(freepar0);
    if(fitmarg->IsFit2D()) truepar1 = fparm->InitValue(freepar1);	
    chisqtrue0 = (double)fitmarg->GetChisqMargTrueRead(0.,truepar1,skobs_fakedata,true);
    chisqtrue1 = (double)fitmarg->GetChisqMargTrueRead(1.,truepar1,skobs_fakedata,true);
    
    // Draw Full Chisq
    if(DrawFullChisq){
      fitmarg->FillChisqMarg(skobs_fakedata,true);
      fitmarg->MinimizeHisto();
      bestfitFullChisq0 = fitmarg->GetParamBestFit(0);
      if(fitmarg->IsFit2D()) bestfitFullChisq1 = fitmarg->GetParamBestFit(1); 
      chisqminFullChisq = (double)fitmarg->GetChisqBestFit();
      fitmarg->BuildDeltaChisqMarg();
    }

    // Get delta chisq
    deltachisq01 = chisqtrue0 - chisqtrue1;
    deltachisq0BF = chisqtrue0 - chisqmin;
 
    // Get the fakedata histogram

    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
      else if(s == kSKsRHCFCFV1RingELike){
	//cout << "    " << util::AsStringCompact(s) << endl;
	hErecoToy = new TH1D(*skobs_fakedata->ErecoSpectrum(s));
	//cout <<  skobs_fakedata->ErecoSpectrum(s)->Integral() << endl;
      }
    }
    //cout << fitmarg->GetNSKsamples() << endl;


    // Fill the output tree
   
    if(fitmarg->IsFit2D()){
      hchisq2D_out = (TH2F*)fitmarg->GetChisqMarg2D()->Clone("hchisq2D");
      //if(minimizeHisto) hstep2D_out = (TH2F*)fitmarg->GetChisqStep_2D()->Clone("hstep2D");
    }
    else{
      hchisq1D_out = (TH1F*)fitmarg->GetChisqMarg1D()->Clone("hchisq1D");
      //if(minimizeHisto) hstep1D_out = (TH1F*)fitmarg->GetChisqStep_1D()->Clone("hstep1D");
    }
            
    fitnt->Fill();

         
    // Print output
    
    cout << endl;
    cout << "Minimum Marginalized Chisq = " << chisqmin << endl;        
    cout << "bestfit 0 = " << bestfit0;
    cout << " +" << bestfiterrPos0 
	 << " -" << bestfiterrNeg0 
	 << " (MINOS)"
	 << endl;
    cout << "bestfit 0 = " << bestfit0
	 << " +" << bestfiterrSym0 
	 << " -" << bestfiterrSym0 
	 << " (HESSE)"
	 << endl;
    if(fitmarg->IsFit2D()){
      cout << "bestfit 1 = " << bestfit1; 
      cout << " +" << bestfiterrPos1
	   << " -" << bestfiterrNeg1 
	   << " (MINOS)"
	   << endl;
      cout << "bestfit 1 = " << bestfit1
	   << " +" << bestfiterrSym1
	   << " -" << bestfiterrSym1 
	   << " (HESSE)"
	   << endl;
    }
    cout << endl;
    cout << "# of calls of CalcChisqMargPt(): " 
	 << nchisqcalls << endl;
    cout << endl;
    cout << "Chisq True (beta=0) = " << chisqtrue0 << endl;
    cout << "Chisq True (beta=1) = " << chisqtrue1 << endl;
    cout << "Delta Chisq = Chisq(beta=0) - Chisq(beta=1) = " << deltachisq01 << endl;
    cout << "Delta Chisq = Chisq(beta=0) - Chisq(betaBF) = " << deltachisq0BF << endl;
    cout << endl;
    cout << "Minimum of full Chisq histogram:" << endl;
    cout << "  - chisqmin = " << chisqminFullChisq << endl; 
    cout << "  - bestfit0 = " << bestfitFullChisq0 << endl; 
    if(fitmarg->IsFit2D())
      cout << "  - bestfit1 = " << bestfitFullChisq1 << endl; 
    cout << endl;

     



    TCanvas *cErecoToy = new TCanvas("cErecoToy","cErecoToy");
    hErecoToy->DrawCopy();

    TCanvas *cChisqMarg = new TCanvas("cChisqMarg","cChisqMarg");
    fitmarg->GetChisqMarg1D()->DrawCopy();    
    
    TCanvas *cChisqStep = new TCanvas("cChisqStep","cChisqStep");
    fitmarg->GetChisqStep_1D()->DrawCopy();
    
    TCanvas *cDeltaChisqMarg =  new TCanvas("cDeltaChisqMarg","cDeltaChisqMarg");
    fitmarg->GetDeltaChisqMarg1D()->DrawCopy();


    // Clean-up
    delete skobs_fakedata;  

    fitmarg->CleanUp();
    delete fitmarg;

  } // for iexperiment

   



  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();


  return;
}
