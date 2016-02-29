
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

#include "Build.h"
#include "AnalysisOptions.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

void run_fit_marginalization_readbinary_joint_2d
(

 ////////////////////////// GENERAL PARAMETERS /////////////////////

 const char * tag                      = "fitmarg_joint_PROVA",    /* tag for output file */

 int          first_experiment         = 0, //4404, //4739,                       /* */
 int          last_experiment          = 0, //4404, //4739,                       /* */

 //const char * xmlfile = "data_run1to6/skdata_run1to6.xml",
 const char * xmlfile = "../outputs/test/toymc_prova_joint-0.xml",

 const char * rootfile = "../outputs/test/toymc_prova_joint-0.root",

 bool using_xml_file = false,

 bool         dotoys               = 1,

 // Mach3, T2K Run1-4 + reactors best fit (A)
 double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp        = -1.601, //-99,    /* deltacp, */
 double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */

 double       input_mixing23Bar    = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_dm23sqBar      = 2.509E-3, //-99,    /* |Dm^2_{23}| */

 // // Mach3, T2K Run1-4 + reactors best fit - dcp=0 (B)
 // double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 0., //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */

 // // Mach3, T2K only Run1-4 best fit (C)
 // double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0422,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 1.91, //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */

 // // Mach3, T2K Run1-4 + reactors best fit - th23 non maximal (D)
 // double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.4, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = -1.601, //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */

 const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_NIWG2015c_ExpSummer2016_75_25.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_FullStat.xml",
 //const char * config_file      = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml",
 //const char * config_file      = "analysis_configuration_run1to4_joint_Nu6.57E20POT.xml",
 
 bool         stat_fluctuations        = false,    /* include statistical fluctuation in toy MC expts? (not for marginalization) */
 int          tweak_systematics        = 0,       /* generate toy (both fakedata sets and marginalization) MC with randomized systematics? */
 
 const char * fit_param_list        = "nuosc3fbar_numudissap_allsyst", // dcp
 //const char * fit_param_list        = "nuosc3f_12fix+nuisance", // dcp
 
 const char * dmsqscheme               = "auto",   /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamchoice           = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * hierarchy                = "normal",     /* normal or inverted */
 const char * mixparamoctant           = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

 ////////////////////// MARGINALIZATION PARAMETERS ///////////////////
 
 //const char * infilename = "submit_prova",
 //const char * infilename = "submit_prova_sinsqth13Vsdcp",
 //const char * infilename = "submit_prova_sinsqth13Vsdcp_react_manyfiles",
 const char * infilename = "Run14",
 //const char * infilename = "submit_prova_sinsqth23Vsdmsq23",
 //const char * infilename = "submit_prova_sinsqth23Vsdmsq23_manyfiles",
 
 const int    nfiles               = 1, //4, //1, //201        // # of files to read

 const int    NtoysMarg            = 10, //500, //10000, //1,  // # toys marginalization
 
 // const char * freeparname0        = "mixing13",  // free parameter
 // const double freeparmin0         = 0.015, // minimum value of free parameter
 // const double freeparmax0         = 0.036, // maxinum value of free parameter
 // const int    NptScan0            = 81, //201, //21,      // # pts of free parameter to scan
 // const char * freeparname1        = "deltacp",  // free parameter
 // const double freeparmin1         = -3.14159265358979312e+00,         // minimum value of free parameter
 // const double freeparmax1         = +3.14159265358979312e+00, // 10., 1         // maxinum value of free parameter
 // const int    NptScan1            = 51, //201, //21, // 801,          // # pts of free parameter to scan
 
 // const char * freeparname0        = "mixing13",  // free parameter
 // const double freeparmin0         = 0., // minimum value of free parameter
 // const double freeparmax0         = 0.113, // maxinum value of free parameter
 // const int    NptScan0            = 3,       // # pts of free parameter to scan
 // const char * freeparname1        = "deltacp",  // free parameter
 // const double freeparmin1         = -3.14159265358979312e+00,         // minimum value of free parameter
 // const double freeparmax1         = +3.14159265358979312e+00, // 10., 1         // maxinum value of free parameter
 // const int    NptScan1            = 3, //51, //201, //3, //21, // 801, 2          // # pts of free parameter to scan

// const char * freeparname0        = "mixing13",  // free parameter
// const double freeparmin0         = 0., // minimum value of free parameter
// const double freeparmax0         = 0.113, // maxinum value of free parameter
// const int    NptScan0            = 81,       // # pts of free parameter to scan
// const char * freeparname1        = "deltacp",  // free parameter
// const double freeparmin1         = -3.14159265358979312e+00,         // minimum value of free parameter
// const double freeparmax1         = +3.14159265358979312e+00, // 10., 1         // maxinum value of free parameter
// const int    NptScan1            = 51, //51, //201, //3, //21, // 801, 2          // # pts of free parameter to scan

  const char * freeparname0        = "mixing23Bar",  // free parameter
  const double freeparmin0         = 0.25, // minimum value of free parameter
  const double freeparmax0         = 0.75, // maxinum value of free parameter
  const int    NptScan0            = 81, //81, //201, //21,      // # pts of free parameter to scan
  const char * freeparname1        = "dmsq23Bar",  // free parameter
  const double freeparmin1         = 1.8E-3,         // minimum value of free parameter
  const double freeparmax1         = 3.2E-3, // 10., 1         // maxinum value of free parameter
  const int    NptScan1            = 21, //51, //201, //3, //21, // 801, 2          // # pts of free parameter to scan

 const bool  doMargMH             = false,    // marginalize Mass Hierarchy
 const char * hierarchy_marg      = "normal",     /* normal or inverted */
 
 const TString TypeError          = "",  // "HESSE-MINOS" (contained by the string)  
 const bool DrawFullChisq         = true,

 const bool StoreLikeWeights      = false,

 const string teststat            = "rate+shape",   // "rate-only" or "rate+shape",

 const bool constrain_reactor     = false
)
{

  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  TString opt(fit_param_list); 

  // Analysis configuration
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file));
  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile));
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));

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
  bool randomize_mixing23Bar = (TMath::Abs(input_mixing23Bar + 999) < kASmallNum);
  bool randomize_dm23sqBar   = (TMath::Abs(input_dm23sqBar   + 999) < kASmallNum);

  //Set to defaults
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsq2Theta12;
    if(input_mixing23Bar < -1 && input_mixing23Bar > -900) input_mixing23Bar = conv::kGlobFitSinsq2Theta23;

  }
  else {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsqTheta12;
    if(input_mixing23Bar < -1 && input_mixing23Bar > -900) input_mixing23Bar = conv::kGlobFitSinsqTheta23;

  }
  if(input_dm23sq < -1 && input_dm23sq > -900)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm23sqBar < -1 && input_dm23sqBar > -900)  input_dm23sqBar = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1 && input_dm12sq > -900)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90 && input_deltacp > -900) input_deltacp = conv::kGlobFitDeltaCP;

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
  if(hyp == kHypOsc3FBarDissap) {
    fparm->SetInitValue(kPhypMixing23Bar, input_mixing23Bar);
    fparm->SetInitValue(kPhypDm23sqBar,   input_dm23sqBar);
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
  }
  else{
      cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  if(using_xml_file) cout << "   - Reading experiments from: " << xmlfile2 << endl; 
  else cout << "   - Reading experiments from: " << rootfile2 << endl; 
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
  string oscstr = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq,0,input_mixing23Bar,input_dm23sqBar);
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
  TH2D *hchisq2D_out;
  TH1D *hchisq1D_out;
  TH2D *hpenaltyterm2D_out;
  TH1D *hpenaltyterm1D_out;
  TH2D *hdeltachisq2D_out;
  TH1D *hdeltachisq1D_out;
  TH2D *hstep2D_out;
  TH1D *hstep1D_out;
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
  TH1F *hErecoToy_1Rmu;
  TH1F *hErecoToy_1Re;
  TH1F *hErecoToy_RHC1Rmu;
  TH1F *hErecoToy_RHC1Re;
  TH2F *hErecoToy2D_1Rmu;
  TH2F *hErecoToy2D_1Re;
  TH2F *hErecoToy2D_RHC1Rmu;
  TH2F *hErecoToy2D_RHC1Re;
  TH1F *hErecoExp_1Rmu;
  TH1F *hErecoExp_1Re;
  TH1F *hErecoExp_RHC1Rmu;
  TH1F *hErecoExp_RHC1Re;
  TH2F *hErecoExp2D_1Rmu;
  TH2F *hErecoExp2D_1Re;
  TH2F *hErecoExp2D_RHC1Rmu;
  TH2F *hErecoExp2D_RHC1Re;
  
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
  fitnt->Branch("hpenaltyterm2D", &hpenaltyterm2D_out );
  fitnt->Branch("hpenaltyterm1D", &hpenaltyterm1D_out );
  fitnt->Branch("hdeltachisq2D", &hdeltachisq2D_out );
  fitnt->Branch("hdeltachisq1D", &hdeltachisq1D_out );
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
  fitnt->Branch("hErecoToy_1Rmu",&hErecoToy_1Rmu);
  fitnt->Branch("hErecoToy_1Re",&hErecoToy_1Re);
  fitnt->Branch("hErecoToy_RHC1Rmu",&hErecoToy_RHC1Rmu);
  fitnt->Branch("hErecoToy_RHC1Re",&hErecoToy_RHC1Re);
  fitnt->Branch("hErecoToy2D_1Rmu",&hErecoToy2D_1Rmu);
  fitnt->Branch("hErecoToy2D_1Re",&hErecoToy2D_1Re);
  fitnt->Branch("hErecoToy2D_RHC1Rmu",&hErecoToy2D_RHC1Rmu);
  fitnt->Branch("hErecoToy2D_RHC1Re",&hErecoToy2D_RHC1Re);
  fitnt->Branch("hErecoExp_1Rmu",&hErecoExp_1Rmu);
  fitnt->Branch("hErecoExp_1Re",&hErecoExp_1Re);
  fitnt->Branch("hErecoExp_RHC1Rmu",&hErecoExp_RHC1Rmu);
  fitnt->Branch("hErecoExp_RHC1Re",&hErecoExp_RHC1Re);
  fitnt->Branch("hErecoExp2D_1Rmu",&hErecoExp2D_1Rmu);
  fitnt->Branch("hErecoExp2D_1Re",&hErecoExp2D_1Re);
  fitnt->Branch("hErecoExp2D_RHC1Rmu",&hErecoExp2D_RHC1Rmu);
  fitnt->Branch("hErecoExp2D_RHC1Re",&hErecoExp2D_RHC1Re);  

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

    if(hyp = kHypOsc3FBarDissap){
      if (randomize_mixing23Bar) { input_mixing23Bar = util::RandomizeMixing23(); }
      if (randomize_dm23sqBar  ) { input_dm23sqBar   = util::RandomizeDm23sq();   }
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
    
    if(hyp = kHypOsc3FBarDissap){
      fparm->SetInitValue(kPhypMixing23Bar, input_mixing23Bar);
      fparm->SetInitValue(kPhypDm23sqBar,   input_dm23sqBar);
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
    //bool using_xml_file = true;
    if(!dotoys){
      if(using_xml_file){
	cout << "Use xml file" << endl;
	skobs_fakedata->Read(xmlfile2.Data(),iexperiment,using_xml_file);
      }
      else{
	cout << "Use root file" << endl;
	skobs_fakedata->Read(rootfile2.Data(),iexperiment,using_xml_file);
      }
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
    
    // Set if to store full likelihood for each toy (after fitmarg->Init())
    fitmarg->SetStoreLikeWeights(StoreLikeWeights);
    
    // Set Mass Hierarchy marginalization
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy_marg); // needed to build the right binary file name 
    fitmarg->SetMargMassHierarchy(doMargMH); // false by default     

    // Set fit parameter list
    fitmarg->SetFitParamList(fit_param_list);
    fitmarg->PrintPhysHypothesis();
   

    // Build the map fileID - sub-range
    
    fitmarg->SetMaxIterationSingleFit(2);  // Max # of iterations per fit
    fitmarg->SetMaxFitNumber(3); // Max # of fits max # of iterations

    fitmarg->SetTagFileName(infilename);
    fitmarg->SetFileID(nfiles);
    

    // Penalty term
    
    if(constrain_reactor){
      
      //double sinsq2th13_mean = 0.093; // PDG 2014 (sinsq2th13)
      //double sinsq2th13_sigma = 0.008; // PDG 2014 (sinsq2th13)
      double sinsq2th13_mean = conv::kPrior_Sinsq2th13_Used_Mean; // Used (sinsq2th13)
      double sinsq2th13_sigma = conv::kPrior_Sinsq2th13_Used_RMS; // Used (sinsq2th13)
      
      fitmarg->UsePenaltyTerm(true,false);
      fitmarg->SetPenaltyTerm(sinsq2th13_mean, sinsq2th13_sigma); // to be used after UsePenaltyTerm()
    }

    //
    //////////////// Implement the minimization ///////////////
    //

    //
    ///////////////////////////////////////////////////////////
    //

    
    // Draw Full Chisq
    if(DrawFullChisq){
      fitmarg->FillChisqMarg(skobs_fakedata,true);
      fitmarg->MinimizeHisto();
      bestfitFullChisq0 = fitmarg->GetParamBestFit(0);
      if(fitmarg->IsFit2D()) bestfitFullChisq1 = fitmarg->GetParamBestFit(1); 
      chisqminFullChisq = (double)fitmarg->GetChisqBestFit();
      fitmarg->BuildDeltaChisqMarg();
    }

    nchisqcalls = fitmarg->GetNCallsCalcChisqMargPt();

    // Set best-fit parameters in FitMarginalization object
    // Otherwise error computation is done with respect to a wrong minimum
    fitmarg->SetChisqBestFit(chisqminFullChisq);
    fitmarg->SetParamBestFit(0,bestfitFullChisq0);
    if(fitmarg->IsFit2D()) fitmarg->SetParamBestFit(1,bestfitFullChisq1);
    
    // Build the DeltaChisqMarg after the best-fit Chisq has been set
    // Otherwise DeltaChisqMarg becomes negative if shifted wrt a local minimum
    fitmarg->BuildDeltaChisqMarg();
	
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
    chisqtrue0 = (double)fitmarg->GetChisqMargTrueRead(truepar0,truepar1,skobs_fakedata,true);

    // Store toy likelihoods (needed to draw best-fit spectra w/ marginalization)
    if(StoreLikeWeights){
      fitmarg->WriteLikeToys2Binary(tag);
    }

    // Get delta chisq
    if(chisqmin==conv::kBadResult) deltachisq0BF = chisqtrue0 - chisqmin;
    else deltachisq0BF = chisqtrue0 - chisqminFullChisq;

    // Get the fakedata histogram

    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;

#ifndef __USE_2D_TEMPLATES__
      else if(s == kSKsRHCFCFV1RingMuLike){
	hErecoToy_RHC1Rmu = (TH1F*)skobs_fakedata->ErecoSpectrum(s);
      	hErecoExp_RHC1Rmu = (TH1F*)skexp->CurrPredictedErecoSpectrum(s);
      }
      else if(s == kSKsRHCFCFV1RingELike){
	hErecoToy_RHC1Re = (TH1F*)skobs_fakedata->ErecoSpectrum(s);
	hErecoExp_RHC1Re = (TH1F*)skexp->CurrPredictedErecoSpectrum(s);  
      }
      else if(s == kSKsFCFV1RingMuLike){
	hErecoToy_1Rmu = (TH1F*)skobs_fakedata->ErecoSpectrum(s);
	hErecoExp_1Rmu = (TH1F*)skexp->CurrPredictedErecoSpectrum(s);  
      }
      else if(s == kSKsFCFV1RingELike){
	hErecoToy_1Re = (TH1F*)skobs_fakedata->ErecoSpectrum(s);
      	hErecoExp_1Re = (TH1F*)skexp->CurrPredictedErecoSpectrum(s);  
      }
      //cout << "    " << util::AsStringCompact(s) << endl;
#else
      else if(s == kSKsRHCFCFV1RingMuLike){
	hErecoToy2D_RHC1Rmu = (TH2F*)skobs_fakedata->TwoDSpectrum(s);
	hErecoExp2D_RHC1Rmu = (TH2F*)skexp->CurrPredictedErecoVsQ2Spectrum(s);  
      }
      else if(s == kSKsRHCFCFV1RingELike){
	hErecoToy2D_RHC1Re = (TH2F*)skobs_fakedata->TwoDSpectrum(s);
	hErecoExp2D_RHC1Re = (TH2F*)skexp->CurrPredictedErecoVsQ2Spectrum(s);  
      }
      else if(s == kSKsFCFV1RingMuLike){
	hErecoToy2D_1Rmu = (TH2F*)skobs_fakedata->TwoDSpectrum(s);
    	hErecoExp2D_1Rmu = (TH2F*)skexp->CurrPredictedErecoVsQ2Spectrum(s);  
      }
      else if(s == kSKsFCFV1RingELike){
	hErecoToy2D_1Re = (TH2F*)skobs_fakedata->TwoDSpectrum(s);
    	hErecoTExp2D_1Re = (TH2F*)skexp->CurrPredictedErecoVsQ2Spectrum(s);  
      }
      //cout << "    " << util::AsStringCompact(s) << endl;      
#endif
    }


    // Fill the output tree
   
    if(fitmarg->IsFit2D()){
      hchisq2D_out = (TH2D*)fitmarg->GetChisqMarg2D()->Clone("hchisq2D");
      hdeltachisq2D_out = (TH2D*)fitmarg->GetDeltaChisqMarg2D()->Clone("hdeltachisq2D");
      hpenaltyterm2D_out = (TH2D*)fitmarg->GetPenaltyTerm2D()->Clone("hpenaltyterm2D");
      //if(minimizeHisto) hstep2D_out = (TH2D*)fitmarg->GetChisqStep_2D()->Clone("hstep2D");
    }
    else{
      hchisq1D_out = (TH1D*)fitmarg->GetChisqMarg1D()->Clone("hchisq1D");
      hdeltachisq1D_out = (TH1D*)fitmarg->GetDeltaChisqMarg1D()->Clone("hcdeltahisq1D");
      hpenaltyterm1D_out = (TH1D*)fitmarg->GetPenaltyTerm1D()->Clone("hpenaltyterm1D");
      //if(minimizeHisto) hstep1D_out = (TH1D*)fitmarg->GetChisqStep_1D()->Clone("hstep1D");
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
    cout << "Chisq True = " << chisqtrue0 << endl;
    cout << "Delta Chisq = Chisq_true - Chisq_bf = " << deltachisq0BF << endl;
    cout << endl;
    cout << "Minimum of full Chisq histogram:" << endl;
    cout << "  - chisqmin = " << chisqminFullChisq << endl; 
    cout << "  - bestfit0 = " << bestfitFullChisq0 << endl; 
    if(fitmarg->IsFit2D())
      cout << "  - bestfit1 = " << bestfitFullChisq1 << endl; 
    cout << endl;

    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
	
      cout << "Sample " << util::AsStringCompact(s) << endl;
      cout << " - # of observed events = ";       
#ifndef __USE_2D_TEMPLATES__
      cout << skobs_fakedata->ErecoSpectrum(s)->Integral() << endl;
#else
      cout << skobs_fakedata->TwoDSpectrum(s)->Integral() << endl;
#endif
      cout << " - # of expected events = ";       
#ifndef __USE_2D_TEMPLATES__
      cout << skexp->CurrPredictedErecoSpectrum(s)->Integral() << endl;
#else
      cout << skexp->CurrPredictedErecoVsQ2Spectrum(s)->Integral() << endl;
#endif
      cout << endl;
    }
    cout << endl;
 
    // Draw observed spectra

#ifndef __USE_2D_TEMPLATES__
    TCanvas *cErecoToy_1Rmu = new TCanvas("cErecoToy_1Rmu","cErecoToy_1Rmu");
    hErecoToy_1Rmu->DrawCopy();

    TCanvas *cErecoToy_1Re = new TCanvas("cErecoToy_1Re","cErecoToy_1Re");
    hErecoToy_1Re->DrawCopy();

    TCanvas *cErecoToy_RHC1Rmu = new TCanvas("cErecoToy_RHC1Rmu","cErecoToy_RHC1Rmu");
    hErecoToy_RHC1Rmu->DrawCopy();

    TCanvas *cErecoToy_RHC1Re = new TCanvas("cErecoToy_RHC1Re","cErecoToy_RHC1Re");
    hErecoToy_RHC1Re->DrawCopy();
#else

    TCanvas *cErecoToy2D_1Rmu = new TCanvas("cErecoToy2D_1Rmu","cErecoToy2D_1Rmu");
    hErecoToy2D_1Rmu->DrawCopy("colz");

    TCanvas *cErecoToy2D_1Re = new TCanvas("cErecoToy2D_1Re","cErecoToy2D_1Re");
    hErecoToy2D_1Re->DrawCopy("colz");

    TCanvas *cErecoToy2D_RHC1Rmu = new TCanvas("cErecoToy2D_RHC1Rmu","cErecoToy2D_RHC1Rmu");
    hErecoToy2D_RHC1Rmu->DrawCopy("colz");

    TCanvas *cErecoToy2D_RHC1Re = new TCanvas("cErecoToy2D_RHC1Re","cErecoToy2D_RHC1Re");
    hErecoToy2D_RHC1Re->DrawCopy("colz");

#endif

    // Draw expected spectra
    
#ifndef __USE_2D_TEMPLATES__
    TCanvas *cErecoExp_1Rmu = new TCanvas("cErecoExp_1Rmu","cErecoExp_1Rmu");
    hErecoExp_1Rmu->DrawCopy();

    TCanvas *cErecoExp_1Re = new TCanvas("cErecoExp_1Re","cErecoExp_1Re");
    hErecoExp_1Re->DrawCopy();

    TCanvas *cErecoExp_RHC1Rmu = new TCanvas("cErecoExp_RHC1Rmu","cErecoExp_RHC1Rmu");
    hErecoExp_RHC1Rmu->DrawCopy();

    TCanvas *cErecoExp_RHC1Re = new TCanvas("cErecoExp_RHC1Re","cErecoExp_RHC1Re");
    hErecoExp_RHC1Re->DrawCopy();
#else

    TCanvas *cErecoExp2D_1Rmu = new TCanvas("cErecoExp2D_1Rmu","cErecoExp2D_1Rmu");
    hErecoExp2D_1Rmu->DrawCopy("colz");

    TCanvas *cErecoExp2D_1Re = new TCanvas("cErecoExp2D_1Re","cErecoExp2D_1Re");
    hErecoExp2D_1Re->DrawCopy("colz");

    TCanvas *cErecoExp2D_RHC1Rmu = new TCanvas("cErecoExp2D_RHC1Rmu","cErecoExp2D_RHC1Rmu");
    hErecoExp2D_RHC1Rmu->DrawCopy("colz");

    TCanvas *cErecoExp2D_RHC1Re = new TCanvas("cErecoExp2D_RHC1Re","cErecoExp2D_RHC1Re");
    hErecoExp2D_RHC1Re->DrawCopy("colz");

#endif



   TCanvas *cChisqMarg = new TCanvas("cChisqMarg","cChisqMarg");
   fitmarg->GetChisqMarg2D()->DrawCopy("colz");    
   
   TCanvas *cDeltaChisqMarg =  new TCanvas("cDeltaChisqMarg","cDeltaChisqMarg");
   fitmarg->GetDeltaChisqMarg2D()->DrawCopy("colz");


    // Clean-up
    delete skobs_fakedata;  

    fitmarg->CleanUp();
    delete fitmarg;
    delete hErecoToy_1Rmu;
    delete hErecoToy_1Re;
    delete hErecoToy_RHC1Rmu;
    delete hErecoToy_RHC1Re;
    delete hErecoToy2D_1Rmu;
    delete hErecoToy2D_1Re;
    delete hErecoToy2D_RHC1Rmu;
    delete hErecoToy2D_RHC1Re;
    delete hErecoExp_1Rmu;
    delete hErecoExp_1Re;
    delete hErecoExp_RHC1Rmu;
    delete hErecoExp_RHC1Re;
    delete hErecoExp2D_1Rmu;
    delete hErecoExp2D_1Re;
    delete hErecoExp2D_RHC1Rmu;
    delete hErecoExp2D_RHC1Re;

  } // for iexperiment


  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();


  return;
}
