
//
// ------------------------------------------------------------------------------------------------           
//
// Run toys for marginalization and store into binary files. 
// Toys for different sub-ranges of the free paramters can be run in parallel jobs.
// A file .bin contains the array, while a twin file .bin contains the sub-range of the free parameter
// Each sub-range contains the same number of toys, coming from the same random seed.
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;


void build_marg_toys_joint
(
 const int    fileid               = 0,   // file id

 const int    first_bin_freepar0   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar0    = -999, // last freepar of the sub-range
 const int    first_bin_freepar1   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar1    = -999, // last freepar of the sub-range

 const char * tag                  = "PROVA_SystAll_TemplVArray",
 //const char * tag                  = "PROVA_SystAll",
 //const char * tag                  = "PROVA_SystOsc",
 //const char * tag                  = "PROVA_NoSyst",
 
 // const char * config_file          = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml",
  const char * config_file          = "analysis_configuration.xml",
 //const char * config_file          = "analysis_configuration_run1to6_test_app_4.0108E20POT_fluxJune.xml",
 
 const char * fit_param_list        = "nuosc3f+nuisance", // all osc params are free
 //const char * fit_param_list        = "nuosc3f_12fix+nuisance", // dcp
 //const char * fit_param_list        = "nuosc3f_fixall_fitdcp+nuisance", // fix all osc params (stat-only)
 //const char * fit_param_list        = "nuosc3f_nueapp_fix4+nuisance", // fix all osc params except sinsqth13 and dcp (stat-only)
 
 // Mach3, T2K Run1-4 + reactors best fit (A)
 double       input_mixing12_nh       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13_nh       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23_nh       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp_nh        = -1.601, //-99,    /* deltacp, */
 double       input_dm12sq_nh         = 7.53E-5,    /* Dm^2_{12},  */
 double       input_dm23sq_nh         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 double       input_mixing12_ih       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13_ih       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23_ih       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp_ih        = -1.601, //-99,    /* deltacp, */
 double       input_dm12sq_ih         = 7.53E-5,    /* Dm^2_{12},  */
 double       input_dm23sq_ih         = 2.509E-3, //-99,    /* |Dm^2_{23}| */ 
 double       input_beta_bar          = -99,

 // // Mach3, T2K only Run1-4 best fit (C)
 // double       input_mixing12_nh       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13_nh       = 0.0422,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23_nh       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp_nh        = 1.91, //-99,    /* deltacp, */
 // double       input_dm12sq_nh         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq_nh         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_mixing12_ih       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13_ih       = 0.0422,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23_ih       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp_ih        = 1.91, //-99,    /* deltacp, */
 // double       input_dm12sq_ih         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq_ih         = 2.509E-3, //-99,    /* |Dm^2_{23}| */ 
 // double       input_beta_bar          = -99,

 // const char * freeparname0          = "deltacp",  // free parameter
 // const double freeparmin0           = -3.14159265358979312e+00,     // minimum value of free parameter
 // const double freeparmax0           = +3.14159265358979312e+00,     // maxinum value of free parameter
 // const int    NptScan0              = 3, //3, //21,  // # pts of free parameter to scan
 // const char * freeparname1          = "",  // free parameter
 // const double freeparmin1           = 0.,     // minimum value of free parameter
 // const double freeparmax1           = 0.,     // maxinum value of free parameter
 // const int    NptScan1              = 1,  // # pts of free parameter to scan
 
 const char * freeparname0          = "mixing13",  // free parameter
 const double freeparmin0           = 0.,     // minimum value of free parameter
 const double freeparmax0           = 0.4,     // maxinum value of free parameter
 const int    NptScan0              = 2, //81,  // # pts of free parameter to scan
 const char * freeparname1          = "deltacp",  // free parameter
 const double freeparmin1           = -3.14159265358979312e+00,     // minimum value of free parameter
 const double freeparmax1           = +3.14159265358979312e+00,     // maxinum value of free parameter
 const int    NptScan1              = 2, //51, //3, //21,  // # pts of free parameter to scan

 const int    NtoysMarg            = 100, //3, //500  // # toys marginalization 
 const int tweak_systematics       = 1,      /* generate toy MC with randomized systematics? */

 const bool   doMargMH             = 0, //true,      // marginalize Mass Hierarchy
 const double PriorProbNH          = 0.5, // Prior for MH --> give the prob. of NH as input

 const bool   doRethrowUnphysical  = true,    // rethrow unphysical systematic toys
 
 const char * dmsqscheme           = "auto",   /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy            = "normal",    /* normal or inverted */
 const char * mixparamchoice       = "sinsqtheta",  /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant       = "two",         /* one==kNuOctantOne two==kNuOctantTwo */
 
 //TString file_priors = "../../inputs/priors/joint_T2Kreact_run14/Priors_T2Kreact_run14_Mach3_22Maj15.root", // Set to "" if no user priors are used (VALOR)
 TString file_priors = "", // Set to "" if no user priors are used (VALOR)
 
 bool reactors       = false // reactor constraint
 )
{  

  FileStat_t dummy;
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

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

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

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
    // Parameter 0 must be always defined
    cout << "Parameter 0 is not defined" << endl;
    exit(1);
  }

  PhysParamId_t freepar1;
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
    cout << "Parameter 1 is not defined" << endl;
    cout << "1D fit" << endl;
    freepar1 = conv::kPhypUndefined;
  }



  // Implement the constructor giving the following input parameters:
  // free param, ntoys, # pt to scan, min par, max par.
  // If min and max not passed, use default ones for that parameter (see FitMarginalization.h)
  
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
  fitmarg->Init(skexp); // Before need SetPhysParFree() to initialize the arrays
  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  fitmarg->PrintPhysHypothesis();
  //Saving param lists test
  fitmarg->SetSaveParamLists(1);
  //Set marginalization oscillation parameters to defaults  

  FitParamList *fparm    = fitmarg->GetFitParamList();
  FitParamList *fparm_nh = fitmarg->GetFitParamList_NH();
  FitParamList *fparm_ih = fitmarg->GetFitParamList_IH();
  PhysHypothesis_t hyp   = fparm->PhysicsHypothesis();

  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13_nh == -99) input_mixing13_nh = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_nh == -99) input_mixing23_nh = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_nh == -99) input_mixing12_nh = conv::kGlobFitSinsq2Theta12;

    if(input_mixing13_ih == -99) input_mixing13_ih = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_ih == -99) input_mixing23_ih = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_ih == -99) input_mixing12_ih = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13_nh == -99) input_mixing13_nh = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_nh == -99) input_mixing23_nh = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_nh == -99) input_mixing12_nh = conv::kGlobFitSinsqTheta12;

    if(input_mixing13_ih == -99) input_mixing13_ih = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_ih == -99) input_mixing23_ih = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_ih == -99) input_mixing12_ih = conv::kGlobFitSinsqTheta12;
  }

  if(input_dm23sq_nh == -99)  input_dm23sq_nh = conv::kGlobFitDm23sq;
  if(input_dm12sq_nh == -99)  input_dm12sq_nh = conv::kGlobFitDm12sq;
  if(input_deltacp_nh == -99) input_deltacp_nh = conv::kGlobFitDeltaCP;

  if(input_dm23sq_ih == -99)  input_dm23sq_ih = conv::kGlobFitDm23sq;
  if(input_dm12sq_ih == -99)  input_dm12sq_ih = conv::kGlobFitDm12sq;
  if(input_deltacp_ih == -99) input_deltacp_ih = conv::kGlobFitDeltaCP;

  if(input_beta_bar == -99) input_beta_bar = conv::kGlobFitBetaBar;

  // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
  if(hyp == kHypOsc3F) {
    fparm_nh->SetInitValue(kPhypMixing12, input_mixing12_nh);
    fparm_nh->SetInitValue(kPhypMixing13, input_mixing13_nh);
    fparm_nh->SetInitValue(kPhypMixing23, input_mixing23_nh);
    fparm_nh->SetInitValue(kPhypCPV    ,  input_deltacp_nh);
    fparm_nh->SetInitValue(kPhypDm12sq,   input_dm12sq_nh);
    fparm_nh->SetInitValue(kPhypDm23sq,   input_dm23sq_nh);

    fparm_ih->SetInitValue(kPhypMixing12, input_mixing12_ih);
    fparm_ih->SetInitValue(kPhypMixing13, input_mixing13_ih);
    fparm_ih->SetInitValue(kPhypMixing23, input_mixing23_ih);
    fparm_ih->SetInitValue(kPhypCPV    ,  input_deltacp_ih);
    fparm_ih->SetInitValue(kPhypDm12sq,   input_dm12sq_ih);
    fparm_ih->SetInitValue(kPhypDm23sq,   input_dm23sq_ih);
  }
  else if(hyp == kHypOsc3FBar) {
    fparm_nh->SetInitValue(kPhypMixing12, input_mixing12_nh);
    fparm_nh->SetInitValue(kPhypMixing13, input_mixing13_nh);
    fparm_nh->SetInitValue(kPhypMixing23, input_mixing23_nh);
    fparm_nh->SetInitValue(kPhypCPV    ,  input_deltacp_nh);
    fparm_nh->SetInitValue(kPhypDm12sq,   input_dm12sq_nh);
    fparm_nh->SetInitValue(kPhypDm23sq,   input_dm23sq_nh);    
    fparm_nh->SetInitValue(kPhypMixing12Bar, input_mixing12_nh);
    fparm_nh->SetInitValue(kPhypMixing13Bar, input_mixing13_nh);
    fparm_nh->SetInitValue(kPhypMixing23Bar, input_mixing23_nh);
    fparm_nh->SetInitValue(kPhypCPVBar    ,  input_deltacp_nh);
    fparm_nh->SetInitValue(kPhypDm12sqBar,   input_dm12sq_nh);
    fparm_nh->SetInitValue(kPhypDm23sqBar,   input_dm23sq_nh);

    fparm_ih->SetInitValue(kPhypMixing12, input_mixing12_ih);
    fparm_ih->SetInitValue(kPhypMixing13, input_mixing13_ih);
    fparm_ih->SetInitValue(kPhypMixing23, input_mixing23_ih);
    fparm_ih->SetInitValue(kPhypCPV    ,  input_deltacp_ih);
    fparm_ih->SetInitValue(kPhypDm12sq,   input_dm12sq_ih);
    fparm_ih->SetInitValue(kPhypDm23sq,   input_dm23sq_ih);    
    fparm_ih->SetInitValue(kPhypMixing12Bar, input_mixing12_ih);
    fparm_ih->SetInitValue(kPhypMixing13Bar, input_mixing13_ih);
    fparm_ih->SetInitValue(kPhypMixing23Bar, input_mixing23_ih);
    fparm_ih->SetInitValue(kPhypCPVBar    ,  input_deltacp_ih);
    fparm_ih->SetInitValue(kPhypDm12sqBar,   input_dm12sq_ih);
    fparm_ih->SetInitValue(kPhypDm23sqBar,   input_dm23sq_ih);
  }  
  else if(hyp == kHypOsc3FNormBar) {
    fparm_nh->SetInitValue(kPhypMixing12, input_mixing12_nh);
    fparm_nh->SetInitValue(kPhypMixing13, input_mixing13_nh);
    fparm_nh->SetInitValue(kPhypMixing23, input_mixing23_nh);
    fparm_nh->SetInitValue(kPhypCPV    ,  input_deltacp_nh);
    fparm_nh->SetInitValue(kPhypDm12sq,   input_dm12sq_nh);
    fparm_nh->SetInitValue(kPhypDm23sq,   input_dm23sq_nh);    
    fparm_nh->SetInitValue(kPhypBetaBar,  input_beta_bar);

    fparm_ih->SetInitValue(kPhypMixing12, input_mixing12_ih);
    fparm_ih->SetInitValue(kPhypMixing13, input_mixing13_ih);
    fparm_ih->SetInitValue(kPhypMixing23, input_mixing23_ih);
    fparm_ih->SetInitValue(kPhypCPV    ,  input_deltacp_ih);
    fparm_ih->SetInitValue(kPhypDm12sq,   input_dm12sq_ih);
    fparm_ih->SetInitValue(kPhypDm23sq,   input_dm23sq_ih);    
    fparm_ih->SetInitValue(kPhypBetaBar,  input_beta_bar);
  }
  else if(hyp == kHypOsc3FBarNormBar) {
    fparm_nh->SetInitValue(kPhypMixing12, input_mixing12_nh);
    fparm_nh->SetInitValue(kPhypMixing13, input_mixing13_nh);
    fparm_nh->SetInitValue(kPhypMixing23, input_mixing23_nh);
    fparm_nh->SetInitValue(kPhypCPV    ,  input_deltacp_nh);
    fparm_nh->SetInitValue(kPhypDm12sq,   input_dm12sq_nh);
    fparm_nh->SetInitValue(kPhypDm23sq,   input_dm23sq_nh);    
    fparm_nh->SetInitValue(kPhypMixing12Bar, input_mixing12_nh);
    fparm_nh->SetInitValue(kPhypMixing13Bar, input_mixing13_nh);
    fparm_nh->SetInitValue(kPhypMixing23Bar, input_mixing23_nh);
    fparm_nh->SetInitValue(kPhypCPVBar    ,  input_deltacp_nh);
    fparm_nh->SetInitValue(kPhypDm12sqBar,   input_dm12sq_nh);
    fparm_nh->SetInitValue(kPhypDm23sqBar,   input_dm23sq_nh);
    fparm_nh->SetInitValue(kPhypBetaBar,     input_beta_bar);

    fparm_ih->SetInitValue(kPhypMixing12, input_mixing12_ih);
    fparm_ih->SetInitValue(kPhypMixing13, input_mixing13_ih);
    fparm_ih->SetInitValue(kPhypMixing23, input_mixing23_ih);
    fparm_ih->SetInitValue(kPhypCPV    ,  input_deltacp_ih);
    fparm_ih->SetInitValue(kPhypDm12sq,   input_dm12sq_ih);
    fparm_ih->SetInitValue(kPhypDm23sq,   input_dm23sq_ih);    
    fparm_ih->SetInitValue(kPhypMixing12Bar, input_mixing12_ih);
    fparm_ih->SetInitValue(kPhypMixing13Bar, input_mixing13_ih);
    fparm_ih->SetInitValue(kPhypMixing23Bar, input_mixing23_ih);
    fparm_ih->SetInitValue(kPhypCPVBar    ,  input_deltacp_ih);
    fparm_ih->SetInitValue(kPhypDm12sqBar,   input_dm12sq_ih);
    fparm_ih->SetInitValue(kPhypDm23sqBar,   input_dm23sq_ih);
    fparm_ih->SetInitValue(kPhypBetaBar,     input_beta_bar);
  }
  else if(hyp == kHypOsc3FBarDissap) {
    fparm_nh->SetInitValue(kPhypMixing12, input_mixing12_nh);
    fparm_nh->SetInitValue(kPhypMixing13, input_mixing13_nh);
    fparm_nh->SetInitValue(kPhypMixing23, input_mixing23_nh);
    fparm_nh->SetInitValue(kPhypCPV    ,  input_deltacp_nh);
    fparm_nh->SetInitValue(kPhypDm12sq,   input_dm12sq_nh);
    fparm_nh->SetInitValue(kPhypDm23sq,   input_dm23sq_nh);
    fparm_nh->SetInitValue(kPhypMixing23Bar, input_mixing23_nh);
    fparm_nh->SetInitValue(kPhypDm23sqBar,   input_dm23sq_nh);
    
    fparm_ih->SetInitValue(kPhypMixing12, input_mixing12_ih);
    fparm_ih->SetInitValue(kPhypMixing13, input_mixing13_ih);
    fparm_ih->SetInitValue(kPhypMixing23, input_mixing23_ih);
    fparm_ih->SetInitValue(kPhypCPV    ,  input_deltacp_ih);
    fparm_ih->SetInitValue(kPhypDm12sq,   input_dm12sq_ih);
    fparm_ih->SetInitValue(kPhypDm23sq,   input_dm23sq_ih);
    fparm_ih->SetInitValue(kPhypMixing23Bar, input_mixing23_ih);
    fparm_ih->SetInitValue(kPhypDm23sqBar,   input_dm23sq_ih);
  }
  else if(hyp == kHypNuMuDisapp2F) {
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
  }
  
  cout << "   - True oscillation parameters (Normal Hierarchy): " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing12) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing13) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing23) << endl;
  cout << "       > " << util::AsString(kPhypCPV     ) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypCPV     ) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq  ) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm12sq  ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq  ) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm23sq  ) << " eV^2/c^4" << endl;
  if(hyp == kHypOsc3FBar){
    cout << "       > " << util::AsString(kPhypMixing12Bar) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing12Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing13Bar) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing13Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing23Bar) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing23Bar) << endl;
    cout << "       > " << util::AsString(kPhypCPVBar    )  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypCPVBar     ) << endl;
    cout << "       > " << util::AsString(kPhypDm12sqBar)   << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)   << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
  }
  else if(hyp == kHypOsc3FNormBar) {
    cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypBetaBar    ) << endl;
  }  
  else if(hyp == kHypOsc3FBarNormBar) {
    cout << "       > " << util::AsString(kPhypMixing12Bar)  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing12Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing13Bar)  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing13Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing23Bar) << endl;
    cout << "       > " << util::AsString(kPhypCPVBar    )   << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypCPVBar     ) << endl;
    cout << "       > " << util::AsString(kPhypDm12sqBar)    << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypBetaBar    ) << endl;
  }  
  else if(hyp == kHypOsc3FBarDissap) {
    cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing23Bar) << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
  }

  cout << "   - True oscillation parameters (Inverted Hierarchy): " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing12) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing13) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing23) << endl;
  cout << "       > " << util::AsString(kPhypCPV     ) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypCPV     ) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq  ) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm12sq  ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq  ) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm23sq  ) << " eV^2/c^4" << endl;
  if(hyp == kHypOsc3FBar){
    cout << "       > " << util::AsString(kPhypMixing12Bar) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing12Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing13Bar) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing13Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing23Bar) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing23Bar) << endl;
    cout << "       > " << util::AsString(kPhypCPVBar    )  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypCPVBar     ) << endl;
    cout << "       > " << util::AsString(kPhypDm12sqBar)   << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)   << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
  }
  else if(hyp == kHypOsc3FNormBar) {
    cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypBetaBar    ) << endl;
  }  
  else if(hyp == kHypOsc3FBarNormBar) {
    cout << "       > " << util::AsString(kPhypMixing12Bar)  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing12Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing13Bar)  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing13Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing23Bar) << endl;
    cout << "       > " << util::AsString(kPhypCPVBar    )   << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypCPVBar     ) << endl;
    cout << "       > " << util::AsString(kPhypDm12sqBar)    << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypBetaBar    ) << endl;
  }
  else if(hyp == kHypOsc3FBarDissap) {
    cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing23Bar) << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
  }
  
  cout << *fitmarg->GetFitParamList_NH() << endl;
  //cout << *fitmarg->GetFitParamList_IH() << endl;

  // Set Prior for Mass Hierarchy (set it for NH and it's automatically done for IH as well)
  fitmarg->SetProbNormalHierarchy(PriorProbNH);
  
  // Set Mass Hierarchy marginalization
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
  
  // Decide whether to rethrow or not the unphysical systematic toys 
  fitmarg->SetRethrowSystUnphysical(doRethrowUnphysical);

  // Initialize the tree to store the marg toy values of oscillation parameters
  fitmarg->InitTreeOscToys();

  // Take prior distributions
 
  TString name_prior_th13_dcp_nh = "hPrior_th13_dcp_nh";
  TString name_prior_th13_dcp_ih = "hPrior_th13_dcp_ih";
  TString name_prior_23_nh = "hPrior_L23_nh";
  TString name_prior_23_ih = "hPrior_L23_ih";
 
  TH2D *hPrior_th13_dcp_nh = 0;
  TH2D *hPrior_23_nh = 0;
  TH2D *hPrior_th13_dcp_ih = 0;
  TH2D *hPrior_23_ih = 0;
 
  bool UsePrior_dcp = false;
  bool UsePrior_23 = false;

  if(file_priors!=""){

    UsePrior_dcp = true;
    UsePrior_23 = true;
    
    cout << endl;
    cout << "Take priors from: " << file_priors << endl;
    TFile *fPrior = new TFile(file_priors.Data(),"READ");
    if(!fPrior->IsOpen() || !fPrior){
      cerr << endl;
      cerr << "Error: the file " << file_priors << " is not open!" << endl;
      cerr << endl;
      exit(1);
    }

    hPrior_th13_dcp_nh = (TH2D*) fPrior->Get(name_prior_th13_dcp_nh.Data());        
    if(hPrior_th13_dcp_nh->IsZombie() || !hPrior_th13_dcp_nh){
      cerr << endl;
      cerr << "Error: the histogram theta13 Vs deltaCP normal hierarchy is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
    hPrior_th13_dcp_ih = (TH2D*) fPrior->Get(name_prior_th13_dcp_ih.Data());       
    if(hPrior_th13_dcp_ih->IsZombie() || !hPrior_th13_dcp_ih){
      cerr << endl;
      cerr << "Error: the histogram theta13 Vs deltaCP inverted hierarchy is not taken!" << endl;
      cerr << endl;
      exit(1);
    }

    hPrior_23_nh = (TH2D*) fPrior->Get(name_prior_23_nh.Data()); 
    if(hPrior_23_nh->IsZombie() || !hPrior_23_nh){
      cerr << endl;
      cerr << "Error: the histogram par23_hist_normal is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
    hPrior_23_ih = (TH2D*) fPrior->Get(name_prior_23_ih.Data()); 
    if(hPrior_23_ih->IsZombie() || !hPrior_23_ih){
      cerr << endl;
      cerr << "Error: the histogram par23_hist_inverted is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
    
    // TCanvas *c0 = new TCanvas("c0","c0");
    // hPrior_th13_dcp_nh->Draw("colz");
    // TCanvas *c1 = new TCanvas("c1","c1");
    // hPrior_th13_dcp_ih->Draw("colz");
    // TCanvas *c2 = new TCanvas("c2","c2");
    // hPrior_23_nh->Draw("colz");
    // TCanvas *c3 = new TCanvas("c3","c3");
    // hPrior_23_ih->Draw("colz");
    // return;
    
  } // if(file_priors)
  

  
  // Set theta13 and deltaCP prior
  
  if(UsePrior_dcp){
    fitmarg->SetPrior_NH(kPhypCPV,"user",hPrior_th13_dcp_nh);  // mix13 and dcp must have the same prior
    fitmarg->SetPrior_IH(kPhypCPV,"user",hPrior_th13_dcp_ih);

    fitmarg->SetPrior_NH(kPhypMixing13,"user",hPrior_th13_dcp_nh); 
    fitmarg->SetPrior_IH(kPhypMixing13,"user",hPrior_th13_dcp_ih);

    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar){ 
      fitmarg->SetPrior_NH(kPhypCPVBar,"user",hPrior_th13_dcp_nh);
      fitmarg->SetPrior_IH(kPhypCPVBar,"user",hPrior_th13_dcp_ih);

      fitmarg->SetPrior_NH(kPhypMixing13Bar,"user",hPrior_th13_dcp_nh);
      fitmarg->SetPrior_IH(kPhypMixing13Bar,"user",hPrior_th13_dcp_ih);
    }
  }
  else{
    // deltaCP
    fitmarg->SetPrior_NH(kPhypCPV,"flat",
			 conv::kPrior_deltaCP_Flat_Low,
			 conv::kPrior_deltaCP_Flat_Up
			 );
    fitmarg->SetPrior_IH(kPhypCPV,"flat",
			 conv::kPrior_deltaCP_Flat_Low,
			 conv::kPrior_deltaCP_Flat_Up
			 );
    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar){ 
      fitmarg->SetPrior_NH(kPhypCPVBar,"flat",
			   conv::kPrior_deltaCP_Flat_Low,
			   conv::kPrior_deltaCP_Flat_Up
			   );
      fitmarg->SetPrior_IH(kPhypCPVBar,"flat",
			   conv::kPrior_deltaCP_Flat_Low,
			   conv::kPrior_deltaCP_Flat_Up
			   );
    }
    
    // mixing13
    if(reactors){
      fitmarg->SetPrior_NH(kPhypMixing13,"gaus",
      			   conv::kPrior_Sinsq2th13_Used_Mean,
      			   conv::kPrior_Sinsq2th13_Used_RMS
      			   ); // prior on sinsq2th
      fitmarg->SetPrior_IH(kPhypMixing13,"gaus",
      			   conv::kPrior_Sinsq2th13_Used_Mean,
      			   conv::kPrior_Sinsq2th13_Used_RMS
      			   );   
      if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar){
      	fitmarg->SetPrior_NH(kPhypMixing13Bar,"gaus",
      			     conv::kPrior_Sinsq2th13_Used_Mean,
      			     conv::kPrior_Sinsq2th13_Used_RMS
      			     );
      	fitmarg->SetPrior_IH(kPhypMixing13Bar,"gaus",
      			     conv::kPrior_Sinsq2th13_Used_Mean,
      			     conv::kPrior_Sinsq2th13_Used_RMS
      			     );
      }      
    }
    else{
      fitmarg->SetPrior_NH(kPhypMixing13,"flat",
			   conv::kPrior_Sinsq2th13_Flat_Low,
			   conv::kPrior_Sinsq2th13_Flat_Up
			   ); // prior on sinsq2th
      fitmarg->SetPrior_IH(kPhypMixing13,"flat",
			   conv::kPrior_Sinsq2th13_Flat_Low,
			   conv::kPrior_Sinsq2th13_Flat_Up
			   );   
      if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar){ 
	fitmarg->SetPrior_NH(kPhypMixing13Bar,"flat",
			     conv::kPrior_Sinsq2th13_Flat_Low,
			     conv::kPrior_Sinsq2th13_Flat_Up
			     );
	fitmarg->SetPrior_IH(kPhypMixing13Bar,"flat",
			     conv::kPrior_Sinsq2th13_Flat_Low,
			     conv::kPrior_Sinsq2th13_Flat_Up
			     );
      }
    }
    
  } // if(UsePrior_dcp)

  // Set 23 prior
  if(UsePrior_23){
    fitmarg->SetPrior_NH(kPhypMixing23,"user",hPrior_23_nh);  // mix23 and dmsq must have the same prior
    fitmarg->SetPrior_IH(kPhypMixing23,"user",hPrior_23_ih);
    fitmarg->SetPrior_NH(kPhypDm23sq,"user",hPrior_23_nh);
    fitmarg->SetPrior_IH(kPhypDm23sq,"user",hPrior_23_ih);
    cout << "setting 23bar priors " << endl;
    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar || hyp == kHypOsc3FBarDissap){
      cout << "in here" << endl;
      fitmarg->SetPrior_NH(kPhypMixing23Bar,"user",hPrior_23_nh);
      fitmarg->SetPrior_IH(kPhypMixing23Bar,"user",hPrior_23_ih);
      fitmarg->SetPrior_NH(kPhypDm23sqBar,"user",hPrior_23_nh);
      fitmarg->SetPrior_IH(kPhypDm23sqBar,"user",hPrior_23_ih);
    }
  }
  else{
    // mixing23 
    fitmarg->SetPrior_NH(kPhypMixing23,"flat",
			 conv::kPrior_Sinsq2th23_Flat_Low,
			 conv::kPrior_Sinsq2th23_Flat_Up
			 ); // prior on sinsqth
    fitmarg->SetPrior_IH(kPhypMixing23,"flat",
			 conv::kPrior_Sinsq2th23_Flat_Low,
			 conv::kPrior_Sinsq2th23_Flat_Up
			 );
    cout << "setting 23bar priors " << endl;
    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar || hyp == kHypOsc3FBarDissap){
      cout << "in here" << endl;
      fitmarg->SetPrior_NH(kPhypMixing23Bar,"flat",
			   conv::kPrior_Sinsq2th23_Flat_Low,
			   conv::kPrior_Sinsq2th23_Flat_Up
			   );
      fitmarg->SetPrior_IH(kPhypMixing23Bar,"flat",
			   conv::kPrior_Sinsq2th23_Flat_Low,
			   conv::kPrior_Sinsq2th23_Flat_Up
			   );
    }
    
    // dmsq23
    fitmarg->SetPrior_NH(kPhypDm23sq,"flat",
			 conv::kPrior_Dm23sq_Flat_Low,
			 conv::kPrior_Dm23sq_Flat_Up
			 );
    fitmarg->SetPrior_IH(kPhypDm23sq,"flat",
			 conv::kPrior_Dm23sq_Flat_Low,
			 conv::kPrior_Dm23sq_Flat_Up
			 );
    cout << "setting 23bar priors " << endl;
    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar || hyp == kHypOsc3FBarDissap){
      cout << "in here" << endl;
      fitmarg->SetPrior_NH(kPhypDm23sqBar,"flat",
			   conv::kPrior_Dm23sq_Flat_Low,
			   conv::kPrior_Dm23sq_Flat_Up
			   );
      fitmarg->SetPrior_IH(kPhypDm23sqBar,"flat",
			   conv::kPrior_Dm23sq_Flat_Low,
			   conv::kPrior_Dm23sq_Flat_Up
			   );
    }
  }
  
  // dmsq12
  fitmarg->SetPrior_NH(kPhypDm12sq,"gaus",
		       conv::kPrior_Dm12sq_Used_Mean,
		       conv::kPrior_Dm12sq_Used_RMS
		       );
  fitmarg->SetPrior_IH(kPhypDm12sq,"gaus",
		       conv::kPrior_Dm12sq_Used_Mean,
		       conv::kPrior_Dm12sq_Used_RMS
		       );

  // mixing12
  fitmarg->SetPrior_NH(kPhypMixing12,"gaus",
		       conv::kPrior_Sinsq2th12_Used_Mean,
		       conv::kPrior_Sinsq2th12_Used_RMS
		       ); // prior on sinsq2th
  fitmarg->SetPrior_IH(kPhypMixing12,"gaus",
		       conv::kPrior_Sinsq2th12_Used_Mean,
		       conv::kPrior_Sinsq2th12_Used_RMS
		       );
 

  // Make toys and store in binary files
    
  // Produce toys for different sub-ranges of the free parameter and store in different files.
  // Useful for running on cluster when many toys are needed (e.g. 10k x 200)

  fitmarg->MakeToysMarginalization(tweak_systematics,skexp); 
  fitmarg->WriteToys2Binary(fileid,tag);  

  // Write tree with stored oscillation parameter values 
  fitmarg->WriteOscToys2Tree(fileid,tag);
  

  // Draw prior and toy histograms for marginalized oscillation parameters
 
  TH1D *hToyHisto_dcp_nh;
  TH1D *hToyHisto_Mix13_nh;

  TH1D *hToyHisto_Mix23_nh;
  TH1D *hToyHisto_Dm23_nh;
  
  TH1D *hToyHisto_Mix12_nh;
  TH1D *hToyHisto_Dm12_nh;
  
  TH1D *hToyHisto_Mix23Bar_nh;
  TH1D *hToyHisto_Dm23Bar_nh;
  
  TH1D *hToyHisto_Mix13_ih;
  TH1D *hToyHisto_dcp_ih;

  TH1D *hToyHisto_Mix23_ih;
  TH1D *hToyHisto_Dm23_ih;

  TH1D *hToyHisto_Mix12_ih;
  TH1D *hToyHisto_Dm12_ih;
  
  TH1D *hToyHisto_Mix23Bar_ih;
  TH1D *hToyHisto_Dm23Bar_ih;
  
  TH1D *hToyHisto_MH;

  PhysHypothesis ph(hyp);

  if(ph.Contains(conv::kPhypMixing12) && fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypMixing12)){    
    hToyHisto_Mix12_nh = new TH1D(*fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypMixing12));
    TCanvas *cToyHisto_Mix12_nh = new TCanvas("cToyHisto_Mix12_nh","cToyHisto_Mix12_nh");
    hToyHisto_Mix12_nh->Draw("");
  }


  if(ph.Contains(conv::kPhypDm12sq) && fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypDm12sq)){
    hToyHisto_Dm12_nh = new TH1D(*fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypDm12sq));
    TCanvas *cToyHisto_Dm12_nh = new TCanvas("cToyHisto_Dm12_nh","cToyHisto_Dm12_nh");
    hToyHisto_Dm12_nh->Draw("");
  }
  
  if(ph.Contains(conv::kPhypMixing13) && fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypMixing13)){    
    hToyHisto_Mix13_nh = new TH1D(*fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypMixing13));
    TCanvas *cToyHisto_Mix13_nh = new TCanvas("cToyHisto_Mix13_nh","cToyHisto_Mix13_nh");
    hToyHisto_Mix13_nh->Draw("");
  }
  
  if(ph.Contains(conv::kPhypCPV) && fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypCPV)){
    hToyHisto_dcp_nh = new TH1D(*fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypCPV));
    TCanvas *cToyHisto_dcp_nh = new TCanvas("cToyHisto_dcp_nh","cToyHisto_dcp_nh");
     hToyHisto_dcp_nh->Draw("");
  }

  if(ph.Contains(conv::kPhypMixing23) && fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypMixing23)){
    hToyHisto_Mix23_nh = new TH1D(*fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypMixing23));
    TCanvas *cToyHisto_Mix23_nh = new TCanvas("cToyHisto_Mix23_nh","cToyHisto_Mix23_nh");
    hToyHisto_Mix23_nh->Draw("");
  }
  
  if(ph.Contains(conv::kPhypDm23sq) && fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypDm23sq)){
    hToyHisto_Dm23_nh = new TH1D(*fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypDm23sq));
    TCanvas *cToyHisto_Dm23_nh = new TCanvas("cToyHisto_Dm23_nh","cToyHisto_Dm23_nh");
    hToyHisto_Dm23_nh->Draw("");
  }
  
  if(ph.Contains(conv::kPhypMixing23Bar) && fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypMixing23Bar)){
    hToyHisto_Mix23Bar_nh = new TH1D(*fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypMixing23Bar));
    TCanvas *cToyHisto_Mix23Bar_nh = new TCanvas("cToyHisto_Mix23Bar_nh","cToyHisto_Mix23Bar_nh");
    hToyHisto_Mix23Bar_nh->Draw("");
  }
  
  if(ph.Contains(conv::kPhypDm23sqBar) && fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypDm23sqBar)){
    hToyHisto_Dm23Bar_nh = new TH1D(*fitmarg->GetFitParamList_NH()->GetToyHisto(conv::kPhypDm23sqBar));
    TCanvas *cToyHisto_Dm23Bar_nh = new TCanvas("cToyHisto_Dm23Bar_nh","cToyHisto_Dm23Bar_nh");
    hToyHisto_Dm23Bar_nh->Draw("");
  }

  if(ph.Contains(conv::kPhypMixing12) && fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypMixing12)){
    hToyHisto_Mix12_ih = new TH1D(*fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypMixing12));
    TCanvas *cToyHisto_Mix12_ih = new TCanvas("cToyHisto_Mix12_ih","cToyHisto_Mix12_ih");
    hToyHisto_Mix12_ih->Draw("");
  }

  if(ph.Contains(conv::kPhypDm12sq) && fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypDm12sq)){
    hToyHisto_Dm12_ih = new TH1D(*fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypDm12sq));
    TCanvas *cToyHisto_Dm12_ih = new TCanvas("cToyHisto_Dm12_ih","cToyHisto_Dm12_ih");
    hToyHisto_Dm12_ih->Draw("");
  }

  if(ph.Contains(conv::kPhypMixing13) && fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypMixing13)){
    hToyHisto_Mix13_ih = new TH1D(*fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypMixing13));
    TCanvas *cToyHisto_Mix13_ih = new TCanvas("cToyHisto_Mix13_ih","cToyHisto_Mix13_ih");
    hToyHisto_Mix13_ih->Draw("");
  }

  if(ph.Contains(conv::kPhypCPV) && fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypCPV)){
    hToyHisto_dcp_ih = new TH1D(*fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypCPV));
    TCanvas *cToyHisto_dcp_ih = new TCanvas("cToyHisto_dcp_ih","cToyHisto_dcp_ih");
    hToyHisto_dcp_ih->Draw("");
  }

  if(ph.Contains(conv::kPhypMixing23) && fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypMixing23)){
    hToyHisto_Mix23_ih = new TH1D(*fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypMixing23));
    TCanvas *cToyHisto_Mix23_ih = new TCanvas("cToyHisto_Mix23_ih","cToyHisto_Mix23_ih");
    hToyHisto_Mix23_ih->Draw("");
  }

  if(ph.Contains(conv::kPhypDm23sq) && fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypDm23sq)){
    hToyHisto_Dm23_ih = new TH1D(*fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypDm23sq));
    TCanvas *cToyHisto_Dm23_ih = new TCanvas("cToyHisto_Dm23_ih","cToyHisto_Dm23_ih");
    hToyHisto_Dm23_ih->Draw("");
  }
  
  if(ph.Contains(conv::kPhypMixing23Bar) && fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypMixing23Bar)){
    hToyHisto_Mix23Bar_ih = new TH1D(*fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypMixing23Bar));
    TCanvas *cToyHisto_Mix23Bar_ih = new TCanvas("cToyHisto_Mix23Bar_ih","cToyHisto_Mix23Bar_ih");
    hToyHisto_Mix23Bar_ih->Draw("");
  }
  
  if(ph.Contains(conv::kPhypDm23sqBar) && fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypDm23sqBar)){
    hToyHisto_Dm23Bar_ih = new TH1D(*fitmarg->GetFitParamList_IH()->GetToyHisto(conv::kPhypDm23sqBar));
    TCanvas *cToyHisto_Dm23Bar_ih = new TCanvas("cToyHisto_Dm23Bar_ih","cToyHisto_Dm23Bar_ih");
    hToyHisto_Dm23Bar_ih->Draw("");
  }

  if(fitmarg->GetHistoPriorMH()){
    hPriorHisto_MH = new TH1D(*fitmarg->GetHistoPriorMH());
    TCanvas *cPriorHisto_MH = new TCanvas("cPriorHisto_MH","cPriorHisto_MH");
    hPriorHisto_MH->Draw();
  }

  if(fitmarg->GetHistoToyMH()){
    hToyHisto_MH = new TH1D(*fitmarg->GetHistoToyMH());
    TCanvas *cToyHisto_MH = new TCanvas("cToyHisto_MH","cToyHisto_MH");
    hToyHisto_MH->Draw();
  }

  
  // Write histograms in output file
 
  cout << "Write histograms in output file" << endl;

  char name[200];
  sprintf(name,"%s_histotoys_file%d.root",tag,fileid);
  TFile *fout = new TFile(name,"RECREATE");

  if(hToyHisto_Mix12_nh){
    hToyHisto_Mix12_nh->SetName("hToyHisto_Mix12_nh");
    hToyHisto_Mix12_nh->Write();
  }
  if(hToyHisto_Dm12_nh){
    hToyHisto_Dm12_nh->SetName("hToyHisto_Dm12_nh");
    hToyHisto_Dm12_nh->Write();
  }
  if(hToyHisto_Mix13_nh){
    hToyHisto_Mix13_nh->SetName("hToyHisto_Mix13_nh");
    hToyHisto_Mix13_nh->Write();
  }
  if(hToyHisto_dcp_nh){
  hToyHisto_dcp_nh->SetName("hToyHisto_dcp_nh");
  hToyHisto_dcp_nh->Write();
  }
  if(hToyHisto_Mix23_nh){
    hToyHisto_Mix23_nh->SetName("hToyHisto_Mix23_nh");
    hToyHisto_Mix23_nh->Write();
  }
  if(hToyHisto_Dm23_nh){
    hToyHisto_Dm23_nh->SetName("hToyHisto_Dm23_nh");
    hToyHisto_Dm23_nh->Write();
  }
  if(hToyHisto_Mix23Bar_nh){
    hToyHisto_Mix23Bar_nh->SetName("hToyHisto_Mix23Bar_nh");
    hToyHisto_Mix23Bar_nh->Write();
  }
  if(hToyHisto_Dm23Bar_nh){
    hToyHisto_Dm23Bar_nh->SetName("hToyHisto_Dm23Bar_nh");
    hToyHisto_Dm23Bar_nh->Write();
  }

  if(hToyHisto_Mix12_ih){
    hToyHisto_Mix12_ih->SetName("hToyHisto_Mix12_ih");
    hToyHisto_Mix12_ih->Write();
  }
  if(hToyHisto_Dm12_ih){
    hToyHisto_Dm12_ih->SetName("hToyHisto_Dm12_ih");
    hToyHisto_Dm12_ih->Write();
  }
  if(hToyHisto_Mix13_ih){
    hToyHisto_Mix13_ih->SetName("hToyHisto_Mix13_ih");
    hToyHisto_Mix13_ih->Write();
  }
  if(hToyHisto_dcp_ih){
    hToyHisto_dcp_ih->SetName("hToyHisto_dcp_ih");
    hToyHisto_dcp_ih->Write();
  }
  if(hToyHisto_Mix23_ih){
    hToyHisto_Mix23_ih->SetName("hToyHisto_Mix23_ih");
    hToyHisto_Mix23_ih->Write();
  }
  if(hToyHisto_Dm23_ih){
    hToyHisto_Dm23_ih->SetName("hToyHisto_Dm23_ih");
    hToyHisto_Dm23_ih->Write();
  }
  if(hToyHisto_Mix23Bar_ih){
    hToyHisto_Mix23Bar_ih->SetName("hToyHisto_Mix23Bar_ih");
    hToyHisto_Mix23Bar_ih->Write();
  }
  if(hToyHisto_Dm23Bar_ih){
    hToyHisto_Dm23Bar_ih->SetName("hToyHisto_Dm23Bar_ih");
    hToyHisto_Dm23Bar_ih->Write();
  }

  
  if(hPriorHisto_MH){
    hPriorHisto_MH->SetName("hPriorHisto_MH");
    hPriorHisto_MH->Write();
  }
  if(hToyHisto_MH){
    hToyHisto_MH->SetName("hToyHisto_MH");
    hToyHisto_MH->Write();
  }
  
  
  fitmarg->SaveTreeOfToyParamLists(fout);

  
  fout->Close();

  return;
}
