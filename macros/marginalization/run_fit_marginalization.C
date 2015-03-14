
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

void run_fit_marginalization
(

 ////////////////////////// GENERAL PARAMETERS /////////////////////

 const char * tag                      = "toyexptmarg",    /* tag for output file */
 int          run                      = 0,                /* run number */
 int          first_experiment         = 1,                       /* */
 int          last_experiment          = 1,                       /* */

 //const char * xmlfile        = "fakedata/201411.fakedata/fakedata_201411_0p170e21_A.xml",  
 const char * xmlfile        = "toymc/toymc_nueapp-0.xml",  
 
 bool using_xml_file                   = true,
 bool         dotoys                   = false,
 
 double       input_mixing12           = -99, //0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13           = -99, //0.02146, //0.0251,                                /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23           = -99, //0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp            = -99, //0.,                                /* deltacp, if -999 then randomize x TMath::Pi()!!!!*/
 double       input_dm12sq             = -99, //7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
 double       input_dm23sq             = -99, //2.32E-3,                             /* |Dm^2_{23}| */
 
 double       input_mixing12_bar       = -99, //0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13_bar       = -99, //0.02146, //0.0251,                    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23_bar       = -99, //-99, //0.7, //0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp_bar        = -99, //0.,                                /* deltacp, if -999 then randomize x TMath::Pi()!!!!*/
 double       input_dm12sq_bar         = -99, //7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
 double       input_dm23sq_bar         = -99, //-99, //2.7E-3, //2.32E-3,                             /* |Dm^2_{23}| */

 //const char * config_file              = "analysis_configuration_1.7E20_prefit.xml",  
 //const char * config_file              = "analysis_configuration_1.7E20.xml",  
 const char * config_file              = "analysis_configuration_run1to6_test_app.xml",
 
 bool         stat_fluctuations        = true,                                /* include statistical fluctuation in toy MC expts? (not for marginalization) */
 int          tweak_systematics        = 1,       /* generate toy (both fakedata sets and marginalization) MC with randomized systematics? */
 
 //const char * fit_param_list           = "nuosc3fbar_fixnonbar_fix5+nuisance",
 const char * fit_param_list           = "nuosc3fbar_fixnonbar_12fix+nuisance",
 //const char * fit_param_list           = "nuosc3fbar_numudisapp_constrainall+nuisance",  
 //const char * fit_param_list           = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

 bool rethrow_if_parameters_unphysical = true,   // for both marginalization and fakedata sets

 bool         store_toy_expts          = true,  // store the toy MC expts     
 bool         store_bestfit_spectrum   = false,  // store the best-fit spectra 
 const char * dmsqscheme               = "abs23",   /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy                = "normal",     /* normal or inverted */
 const char * mixparamchoice           = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant           = "one",       /* one==kNuOctantOne two==kNuOctantTwo */
 bool         reactor                  = false,   /* reactor constraint */

 const char *optfixpar                 = "mixing12,mixing13,deltacp,dmsq12",




 ////////////////////// MARGINALIZATION PARAMETERS ///////////////////

 const bool   readtoys             = true,

 // Parameter for chisq minimization
 const bool   minimizeHisto        = true, // if not all the chisq histo is filled
 const int    MaxIterFit           = 2, 
 const int    MaxNfit              = 3,
 
 const int    strategy             = 0,    // 0 --> use steps, 1 --> use intervals

 // Step strategy
 const double start_X              = 0.102,
 const double start_Y              = -1.66504,
 const string StrategyStartPt      = "GridRnd", //"Single", "GridRnd", "GridFix"
 const int    NStartPt_X           = 12,  // # of throws for starting pt extraction
 const int    NStartPt_Y           = 4,  // # of throws for starting pt extraction

 // Interval strategy
 const bool   doUncorrXY           = false, // DON'T USE FOR sin2th13 Vs dcp
 const bool   RedoBorderFit        = true,
 const int    MaxNfitBorder        = 1,

 //const double MinStep0             = 0.,     // min x value for fit steps
 //const double MaxStep0             = 1.,     // max x value for fit steps
 const double MinStep0             = 0.,     // min x value for fit steps
 const double MaxStep0             = 0.3,     // max x value for fit steps
 
 //const double MinStep1             = 0.0002, // min y value for fit steps
 //const double MaxStep1             = 0.02,   // max y value for fit steps
 //const double MinStep1             = -TMath::Pi(), // min y value for fit steps
 const double MinStep1             = -TMath::Pi(), // min y value for fit steps
 const double MaxStep1             = +TMath::Pi(),   // max y value for fit steps
 //

 //const char * outfilename          = "ChisqMarg_mixing23_Vs_dmsq23_Nbins201x201",
 const char * outfilename          = "ChisqMarg_mixing23_Vs_dmsq23_Nbins21x21",
 //const char * outfilename          = "ChisqMarg_prova",
 
 const int    IDoutfile            = 0,
 const int    first_bin_freepar0   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar0    = -999, // last freepar of the sub-range
 const int    first_bin_freepar1   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar1    = -999, // last freepar of the sub-range

 //const char * infilename           = "../../inputs/marginalization/toys/mixing23_1p7E20_postfit_NoRethrow_Merged",
 //const char * infilename           = "submit_mixing23_1p7E20_postfit_NoRethrow",
 //const char * infilename           = "submit_mixing23_1p7E20_postfit_NoRethrow_Merged", // test
 //const char * infilename           = "../../inputs/marginalization/toys/mixing23_dmsq32_1p7E20_postfit_NoRethrow_Merged",
 //const char * infilename           = "MH_FitMarg_mixing13Vsdeltacp_1p7E20_prefit",
 const char * infilename           = "MH_FitMarg_mixing13Vsdeltacp_1p7E20_201x201_prefit",

 const int    nfiles               = 1,         // # of files to read

 const int    NtoysMarg            = 1,   // # toys marginalization 
 //const int    NtoysMarg            = 2,   // # toys marginalization 
 
 //const char * freeparname0        = "mixing23Bar",  // free parameter
 //const double freeparmin0         = 0.,              // minimum value of free parameter
 //const double freeparmax0         = 1.,            // maxinum value of free parameter
 //const int    NptScan0            = 21, //2, //21,        // # pts of free parameter to scan
 const char * freeparname0        = "mixing13Bar",  // free parameter
 const double freeparmin0         = 0.,              // minimum value of free parameter
 const double freeparmax0         = 0.3, //0.05, //0.3,   // maxinum value of free parameter
 const int    NptScan0            = 201, //2, //21,        // # pts of free parameter to scan
 
 //const char * freeparname1        = "",  // free parameter
 const char * freeparname1        = "deltacpBar",  // free parameter
 const double freeparmin1         = -TMath::Pi(), // minimum value of free parameter
 const double freeparmax1         = +TMath::Pi(), // maxinum value of free parameter
 const int    NptScan1            = 201, //21,      // # pts of free parameter to scan
 //const char * freeparname1        = "dmsq23Bar",  // free parameter
 //const double freeparmin1         = 0.0002, // minimum value of free parameter
 //const double freeparmax1         = 0.02, // maxinum value of free parameter
 //const int    NptScan1            = 21, //21,      // # pts of free parameter to scan

 const bool  doMargMH             = true    // marginalize Mass Hierarchy
)
{

  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  TString opt(fit_param_list); 
  TString optfixparams(optfixpar);

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
  else if(hyp == kHypNuMuDisapp2F) {
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
  }

   
  int nexperiments = last_experiment - first_experiment + 1; 
  
  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;

  if(!dotoys){
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
  }
  else{
      cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiments from: " << xmlfile2 << endl; 
  cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl; 
  }
  
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "  Applying reactor constraint: " << ((reactor)  ? "YES" : "NO" ) << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing toy MC experiments:  " << ((store_toy_expts) ?        "YES " : "NO") << endl;
  cout << "     - Storing best-fit spectra: "    << ((store_bestfit_spectrum) ? "YES " : "NO") << endl;
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
  TFile f(Form("%s-run_%d-%s_%d_%d.root", tag, run, oscstr.c_str(), first_experiment, last_experiment), "recreate");

  int nfits=-999;
  int nchisqcalls=-999;
  int niterations=-999;
  int fitstatus=-999;
  double chisqmin=-999;
  double chisqtrue=-999;
  TH2F *hchisq2D_out;
  TH1F *hchisq1D_out;
  TH2F *hstep2D_out;
  TH1F *hstep1D_out;
  float bestfit0=-999;
  float bestfit1=-999;
  
  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("bestfit0", &bestfit0 );
  fitnt->Branch("bestfit1", &bestfit1 );
  fitnt->Branch("hchisq2D", &hchisq2D_out );
  fitnt->Branch("hchisq1D", &hchisq1D_out );
  fitnt->Branch("hstep2D", &hstep2D_out );
  fitnt->Branch("hstep1D", &hstep1D_out );
  fitnt->Branch("fitstatus",&fitstatus,"fitstatus/I");
  fitnt->Branch("chisqmin",&chisqmin,"chisqmin/D");
  fitnt->Branch("chisqtrue",&chisqtrue,"chisqtrue/D");
  fitnt->Branch("nfits",&nfits,"nfits/I");
  fitnt->Branch("nchisqcalls",&nchisqcalls,"nchisqcalls/I");
  fitnt->Branch("niterations",&niterations,"niterations/I");

  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  // Loop over toy-MC experiments
  //for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {
    cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;
    
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
    else if(hyp == kHypNuMuDisapp2F) {
      fparm->SetInitValue(kPhypMixing23, input_mixing23);
      fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
    }
   
    cout << *fparm << endl;
    

    // Reset all output spectra
    dat_array->Clear();
    pdf_array->Clear();
 
    // Reset input true params and add true oscillation params
    iparm->Clear();
    PhysParamList true_oscillation_params = fparm->AsPhysParamList("current");
    iparm->Append(true_oscillation_params);

    // get toy MC with/without systematics
    
    SuperkObservation * skobs_fakedata = new SuperkObservation;
    if(!dotoys && using_xml_file){
      skobs_fakedata->Read(xmlfile2.Data(),iexperiment,using_xml_file);
    }
    else{
      bool applied_negevents_condition = true;
      bool unphysical_tweak = true;
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
        }
      }//applied_negevents_condition
    }

    


  
    /////////////////// SET MARGINALIZATION FIT ////////////////
    
    PhysParamId_t freepar0;
    string sfreeparname0(freeparname0);  
    // NuBar mode
    if(sfreeparname0=="mixing23Bar")      freepar0 = conv::kPhypMixing23Bar;
    else if(sfreeparname0=="mixing13Bar") freepar0 = conv::kPhypMixing13Bar;
    else if(sfreeparname0=="dmsq23Bar")   freepar0 = conv::kPhypDm23sqBar;
    else if(sfreeparname0=="dmsq12Bar")   freepar0 = conv::kPhypDm12sqBar;
    else if(sfreeparname0=="mixing12Bar") freepar0 = conv::kPhypMixing12Bar;
    else if(sfreeparname0=="deltacpBar")  freepar0 = conv::kPhypCPVBar;
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
    int subminbin[2];
    int submaxbin[2];
    subminbin[0] = first_bin_freepar0; submaxbin[0] = last_bin_freepar0;
    subminbin[1] = first_bin_freepar1; submaxbin[1] = last_bin_freepar1;
    
    // Implement the constructor giving the following input parameters:
    FitMarginalization *fitmarg = new FitMarginalization();   
    fitmarg->SetPhysParFree(parfree,ntoys,nbinspar,minpar,maxpar);
    fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
    fitmarg->Init(); // Before need SetPhysParFree() to initialize the arrays
    
    // Set Mass Hierarchy marginalization
    fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
    
    // Set fit parameter list
    fitmarg->SetFitParamList(fit_param_list);
    fitmarg->PrintPhysHypothesis();
    
    if(!readtoys){
      
      // Decide whether to rethrow or not the unphysical systematic toys 
      fitmarg->SetRethrowSystUnphysical(doRethrowUnphysical);
      
      // Produce toys and fill marginalized chisq
      
      // Set priors for not fixed oscillation parameters
      fitmarg->SetPrior(kPhypMixing23Bar,"flat",0,1);
      if(reactor) fitmarg->SetPrior(kPhypMixing13Bar,"gaus",0.02146;0.059); // Daya Bay 2014
      else        fitmarg->SetPrior(kPhypMixing13Bar,"flat",0.,1.);
      fitmarg->SetPrior(kPhypMixing12Bar,"flat",0,1);
      fitmarg->SetPrior(kPhypDm23sqBar,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
      fitmarg->SetPrior(kPhypDm12sqBar,"flat",conv::kMinFitDm12sq,conv::kMaxFitDm12sq);
      fitmarg->SetPrior(kPhypCPVBar,"flat",-TMath::Pi(),+TMath::Pi());
      
      fitmarg->MakeToysMarginalization(tweak_systematics,skexp);     
      
      if(!minimizeHisto){
	fitmarg->FillChisqMarg(skobs_fakedata);
	fitmarg->MinimizeHisto();
      }
      else{
	// Max # of iterations per fit
	fitmarg->SetMaxIterationSingleFit(MaxIterFit); 
	// Max # of fits max # of iterations
	fitmarg->SetMaxFitNumber(MaxNfit);   
	
	// Use chisq minimization with steps
	if(strategy==0){
	  fitmarg->SetMinimInterval(MinStep0,MaxStep0,MinStep1,MaxStep1);
	  fitmarg->SetStartBin(start_X,start_Y); 
	  fitmarg->SetNumStartBin(NStartPt_X,NStartPt_Y);
	  fitmarg->SetExtractStartPt(StrategyStartPt);
	  fitmarg->MinimizeChisqStep(skobs_fakedata);      
	}
	
	// Use chisq minimization with intervals
	else if(strategy==1){
	  // Redo the fit if minimum on the border  
	  fitmarg->SetBorderDoFit(RedoBorderFit); // true by default
	  // Max # of fits max # of iterations if the first minimum is on the border
	  fitmarg->SetBorderMaxFitNumber(MaxNfitBorder); // 2*Ndim by default
	  fitmarg->SetMinimInterval(MinStep0,MaxStep0,MinStep1,MaxStep1);
	  fitmarg->SetMinimUncorrXY(doUncorrXY);
	  fitmarg->MinimizeChisqInterval(skobs_fakedata);      
	}
	
      }
    }
    
    else{
      
      // Read toys from binary files and fill marginalized chisq
      
      for(int ifile=0;ifile<nfiles;ifile++){
	
	cout << endl;
	cout << "////////////////////////////////" << endl;
	cout << endl;
	cout << "Read file ID " << ifile << endl; 
	cout << endl;
	
	fitmarg->ReadToysFromBinary(ifile,infilename);
	
	if(!minimizeHisto){
	  fitmarg->FillChisqMarg(skobs_fakedata);
	  fitmarg->MinimizeHisto();
	}
	else{
	  // Max # of iterations per fit
	  fitmarg->SetMaxIterationSingleFit(MaxIterFit); 
	  // Max # of fits max # of iterations
	  fitmarg->SetMaxFitNumber(MaxNfit);   
	  
	  // Use chisq minimization with steps
	  if(strategy==0){
	    fitmarg->SetMinimInterval(MinStep0,MaxStep0,MinStep1,MaxStep1);
	    fitmarg->SetStartBin(start_X,start_Y); 
	    fitmarg->SetNumStartBin(NStartPt_X,NStartPt_Y);
	    fitmarg->SetExtractStartPt(StrategyStartPt);
	    fitmarg->MinimizeChisqStep(skobs_fakedata);      
	  }
	  
	  // Use chisq minimization with intervals
	  else if(strategy==1){
	    // Redo the fit if minimum on the border  
	    fitmarg->SetBorderDoFit(RedoBorderFit); // true by default
	    // Max # of fits max # of iterations if the first minimum is on the border
	    fitmarg->SetBorderMaxFitNumber(MaxNfitBorder); // 2*Ndim by default
	    fitmarg->SetMinimInterval(MinStep0,MaxStep0,MinStep1,MaxStep1);
	    fitmarg->SetMinimUncorrXY(doUncorrXY);
	    fitmarg->MinimizeChisqInterval(skobs_fakedata);      
	  }
	  
	}
	
      }
    }
    

    
    // Find minimum of marginalized chisq and find best-fit values
    cout << endl;
    cout << "# of calls of CalcChisqMargPt(): " 
	 << fitmarg->GetNCallsCalcChisqMargPt() << endl;
    cout << "Fit status " << fitmarg->GetFitStatus() << endl;
    fitmarg->PrintBestFit();
    

    // Fill the output tree

    if(fitmarg->IsFit2D()){
      hchisq2D_out = (TH2F*)fitmarg->GetChisqMarg2D()->Clone("hchisq2D");
      if(minimizeHisto) hstep2D_out = (TH2F*)fitmarg->GetChisqStep_2D()->Clone("hstep2D");
    }
    else{
      hchisq1D_out = (TH1F*)fitmarg->GetChisqMarg1D()->Clone("hchisq1D");
      if(minimizeHisto) hstep1D_out = (TH1F*)fitmarg->GetChisqStep_1D()->Clone("hstep1D");
    }

    chisqmin = (double)fitmarg->GetChisqBestFit();
    nfits = fitmarg->GetNumFit();
    nchisqcalls = fitmarg->GetNCallsCalcChisqMargPt();
    niterations = fitmarg->GetNumIterations();
    fitstatus = fitmarg->GetFitStatus();

    bestfit0 = fitmarg->GetParamBestFit(0); 
    if(fitmarg->IsFit2D()) bestfit1 = fitmarg->GetParamBestFit(0); 
    

    // Get true marginalized chisq

    double truepar0=-999;
    double truepar1=-999;
    
    truepar0 = fparm->InitValue(freepar0);
    if(fitmarg->IsFit2D()) truepar1 = fparm->InitValue(freepar1);
    chisqtrue = (double)fitmarg->GetChisqMargTrue(truepar0,truepar1,skobs_fakedata);
    
    cout << endl;
    cout << "True Marginalized Chisq for: " << endl;
    cout << "   - " << freeparname0 << " = " << truepar0 << endl;
    cout << "   - " << freeparname1 << " = " << truepar1 << endl;
    cout << "   - chisq = " << chisqtrue << endl;  
    cout << endl;

    fitnt->Fill();




    TCanvas *cChisqMarg = new TCanvas("cChisqMarg","cChisqMarg");
    if(fitmarg->IsFit2D()){
      fitmarg->GetChisqMarg2D()->DrawCopy("colz");
      //fitmarg->GetChisqMarg2D()->DrawCopy("colz text");
      //fitmarg->GetChisqMarg2D()->DrawCopy("LEGO");
    }
    else{
      double min = fitmarg->GetChisqBestFit();
      fitmarg->GetChisqMarg()->SetMinimum(min);
      fitmarg->GetChisqMarg()->DrawCopy();
    }


    if(minimizeHisto){
      TCanvas *cChisqStep = new TCanvas("cChisqStep","cChisqStep");
      if(fitmarg->IsFit2D())
	fitmarg->GetChisqStep_2D()->DrawCopy("colz");
      else
	fitmarg->GetChisqStep_1D()->DrawCopy("");
    }


    // Clean-up
    delete skobs_fakedata;  

  } // for iexperiment

    

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();

  return;
}
