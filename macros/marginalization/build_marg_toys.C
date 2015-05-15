
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

void build_marg_toys
(
 const int    fileid               = 0,   // file id

 const int    first_bin_freepar0   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar0    = -999, // last freepar of the sub-range
 const int    first_bin_freepar1   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar1    = -999, // last freepar of the sub-range

 const char * tag                  = "BetaBar_4p0E20POT_Toy23dcp",

 //const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",   
 //const char * config_file          = "analysis_configuration_1.7E20.xml",   
 const char * config_file          = "analysis_configuration_run1to6_test_app.xml",

 //const char * fit_param_list     = "nuosc3fbar_numudisapp_constrainall+nuisance",  
 //const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_fixnonbar_12fix+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_fixnonbar_fix5+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_nueapp_fixnonbar_fix4+nuisance",
 
 //const char * fit_param_list       = "nuosc3fbar_normbar_fixall+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_normbar_fixosc+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_normbar_fitdcpbar+nuisance",
 //const char * fit_param_list       = "nuosc3f_normbar_fitdcp+nuisance",
 const char * fit_param_list       = "nuosc3f_normbar_fit23_fitdcp+nuisance",

 double       input_mixing12           = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13           = -99,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23           = -99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp            = -99,    /* deltacp, */
 double       input_dm12sq             = -99,    /* Dm^2_{12},  */
 double       input_dm23sq             = -99,    /* |Dm^2_{23}| */
 double       input_mixing12_bar       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13_bar       = -99,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23_bar       = -99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp_bar        = -99,    /* deltacp */
 double       input_dm12sq_bar         = -99,    /* Dm^2_{12} */
 double       input_dm23sq_bar         = -99,    /* |Dm^2_{23}| */
 double       input_beta_bar           = -99,    // numubar->nuebar normalization

 // const char * freeparname0          = "mixing23Bar",  // free parameter
 // const double freeparmin0           = 0,     // minimum value of free parameter
 // const double freeparmax0           = 1,     // maxinum value of free parameter
 // const char * freeparname0          = "mixing13Bar",  // free parameter
 // const double freeparmin0           = 0.,     // minimum value of free parameter
 // const double freeparmax0           = 0.3,     // maxinum value of free parameter
 const char * freeparname0          = "betaBar",  // free parameter
 const double freeparmin0           = 0.,     // minimum value of free parameter
 const double freeparmax0           = 1.,     // maxinum value of free parameter
 const int   NptScan0               = 2.,  // # pts of free parameter to scan

 const char * freeparname1          = "",  // free parameter
 //const char * freeparname1          = "dmsq23Bar",  // free parameter
 //const double freeparmin1           = 0.0002,     // minimum value of free parameter
 //const double freeparmax1           = 0.02,     // maxinum value of free parameter
 //const char * freeparname1          = "deltacpBar",  // free parameter
 const double freeparmin1           = -3.14159265358979312e+00,     // minimum value of free parameter
 const double freeparmax1           = +3.14159265358979312e+00,     // maxinum value of free parameter
 const int    NptScan1              = 1,  // # pts of free parameter to scan

 const int    NtoysMarg            = 1,   // # toys marginalization 
 const int tweak_systematics       = 0,      /* generate toy MC with randomized systematics? */

 const bool   doMargMH             = true,    // marginalize Mass Hierarchy
 const bool   doRethrowUnphysical  = true,    // rethrow unphysical systematic toys

 const char * dmsqscheme           = "auto",   /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy            = "normal",    /* normal or inverted */
 const char * mixparamchoice       = "sinsqtheta",  /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant       = "two",         /* one==kNuOctantOne two==kNuOctantTwo */
 
 //TString file_prior_dcp    = "../../inputs/priors/joint_T2Kreact_run14/Priors_deltaCP_run14_T2Kreact_NH.root", // Set to "" if no user priors are used (VALOR)
 TString file_prior_dcp    = "../../inputs/priors/joint_T2Kreact_run14/deltacp_T2Kreact_run14_Mach3.root", // Set to "" if no user priors are used (Mach3)

 TString file_prior_23     = "../../inputs/priors/joint_T2Kreact_run14/parameters23_T2Kreact_run14_Mach3.root" // Set to "" if no user priors are used (Mach3)

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

  //Set marginalization oscillation parameters to defaults  

  FitParamList *fparm = fitmarg->GetFitParamList();
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();

  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 == -99) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 == -99) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 == -99) input_mixing12 = conv::kGlobFitSinsq2Theta12;
    if(input_mixing13_bar == -99) input_mixing13_bar = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_bar == -99) input_mixing23_bar = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_bar == -99) input_mixing12_bar = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 == -99) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 == -99) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 == -99) input_mixing12 = conv::kGlobFitSinsqTheta12;
    if(input_mixing13_bar == -99) input_mixing13_bar = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_bar == -99) input_mixing23_bar = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_bar == -99) input_mixing12_bar = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq == -99)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq == -99)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp == -99) input_deltacp = conv::kGlobFitDeltaCP;
  if(input_dm23sq_bar == -99)  input_dm23sq_bar = conv::kGlobFitDm23sq;
  if(input_dm12sq_bar == -99)  input_dm12sq_bar = conv::kGlobFitDm12sq;
  if(input_deltacp_bar == -99) input_deltacp_bar = conv::kGlobFitDeltaCP;
  if(input_beta_bar == -99) input_beta_bar = conv::kGlobFitBetaBar;
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
  else if(hyp == kHypOsc3FNormBar) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);    
    fparm->SetInitValue(kPhypBetaBar,  input_beta_bar);
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
  
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing12) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing13) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing23) << endl;
  cout << "       > " << util::AsString(kPhypCPV     ) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypCPV     ) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq  ) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypDm12sq  ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq  ) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypDm23sq  ) << " eV^2/c^4" << endl;
  if(hyp == kHypOsc3FBar){
    cout << "       > " << util::AsString(kPhypMixing12Bar) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing12Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing13Bar) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing13Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing23Bar) << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing23Bar) << endl;
    cout << "       > " << util::AsString(kPhypCPVBar    )  << " = " << fitmarg->GetFitParamList()->InitValue(kPhypCPVBar     ) << endl;
    cout << "       > " << util::AsString(kPhypDm12sqBar)   << " = " << fitmarg->GetFitParamList()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)   << " = " << fitmarg->GetFitParamList()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
  }
  else if(hyp == kHypOsc3FNormBar) {
    cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList()->InitValue(kPhypBetaBar    ) << endl;
  }  
  else if(hyp == kHypOsc3FBarNormBar) {
    cout << "       > " << util::AsString(kPhypMixing12Bar)  << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing12Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing13Bar)  << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing13Bar) << endl;
    cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList()->InitValue(kPhypMixing23Bar) << endl;
    cout << "       > " << util::AsString(kPhypCPVBar    )   << " = " << fitmarg->GetFitParamList()->InitValue(kPhypCPVBar     ) << endl;
    cout << "       > " << util::AsString(kPhypDm12sqBar)    << " = " << fitmarg->GetFitParamList()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
    cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList()->InitValue(kPhypBetaBar    ) << endl;
  }  
  

  // Set Mass Hierarchy marginalization
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
  
  // Decide whether to rethrow or not the unphysical systematic toys 
  fitmarg->SetRethrowSystUnphysical(doRethrowUnphysical);

  
  // Take prior distributions
  
  TH1D *hPrior_dcp;
  TH2D *hPrior_23;
 
  if(file_prior_dcp!=""){
    cout << endl;
    cout << "Take prior for deltaCP from: " << file_prior_dcp << endl;
    TFile *fPrior_dcp = new TFile(file_prior_dcp.Data(),"READ");
    if(!fPrior_dcp->IsOpen() || !fPrior_dcp){
      cerr << endl;
      cerr << "Error: the file is not open!" << endl;
      cerr << endl;
      exit(1);
    }
    hPrior_dcp = (TH1D*) fPrior_dcp->Get("dcp_hist");       // Mach3 
    //hPrior_dcp = (TH1D*) fPrior_dcp->Get("hPrior_Ldcp");  // VALOR
    if(hPrior_dcp->IsZombie() || !hPrior_dcp){
      cerr << endl;
      cerr << "Error: the histogram hPrior_Ldcp is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
  }

  if(file_prior_23!=""){
    cout << endl;
    cout << "Take prior for 23 parameters from: " << file_prior_23 << endl;
    TFile *fPrior_23 = new TFile(file_prior_23.Data(),"READ");
    if(!fPrior_23->IsOpen() || !fPrior_23){
      cerr << endl;
      cerr << "Error: the file is not open!" << endl;
      cerr << endl;
      exit(1);
    }
    hPrior_23 = (TH2D*) fPrior_23->Get("par23_hist_normal"); 
    if(hPrior_23->IsZombie() || !hPrior_23){
      cerr << endl;
      cerr << "Error: the histogram par23_hist_normal is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
  }
  

  // Set priors for not fixed and not free oscillation parameters
    
  // fitmarg->SetPrior(kPhypMixing23,"flat",0,1);
  // fitmarg->SetPrior(kPhypMixing13,"flat",0,0.3);
  // fitmarg->SetPrior(kPhypMixing12,"flat",0,1);
  // fitmarg->SetPrior(kPhypDm23sq,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
  // fitmarg->SetPrior(kPhypDm12sq,"flat",conv::kMinFitDm12sq,conv::kMaxFitDm12sq);
  
  // if(hyp == kHypOsc3FBar){
  //   fitmarg->SetPrior(kPhypMixing23Bar,"flat",0,1);
  //   fitmarg->SetPrior(kPhypMixing13Bar,"flat",0,0.3);
  //   fitmarg->SetPrior(kPhypMixing12Bar,"flat",0,1);
  //   fitmarg->SetPrior(kPhypDm23sqBar,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
  //   fitmarg->SetPrior(kPhypDm12sqBar,"flat",conv::kMinFitDm12sq,conv::kMaxFitDm12sq);
  // }
  // else if(hyp == kHypOsc3FBarNormBar){
  //   fitmarg->SetPrior(kPhypMixing23Bar,"flat",0,1);
  //   fitmarg->SetPrior(kPhypMixing13Bar,"flat",0,0.3);
  //   fitmarg->SetPrior(kPhypMixing12Bar,"flat",0,1);
  //   fitmarg->SetPrior(kPhypDm23sqBar,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
  //   fitmarg->SetPrior(kPhypDm12sqBar,"flat",conv::kMinFitDm12sq,conv::kMaxFitDm12sq);
  //   fitmarg->SetPrior(kPhypBetaBar,"flat",conv::kMinFitBetaBar,conv::kMaxFitBetaBar);
  // }
  // else if(hyp == kHypOsc3FNormBar){
  //   fitmarg->SetPrior(kPhypBetaBar,"flat",conv::kMinFitBetaBar,conv::kMaxFitBetaBar);
  // }


  // Set deltaCP prior
  if(file_prior_dcp!=""){
    fitmarg->SetPrior(kPhypCPV,"user",hPrior_dcp);
    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar) 
      fitmarg->SetPrior(kPhypCPVBar,"user",hPrior_dcp);
  }
  else{
    cout << endl;
    cout << "Use flat prior for deltaCP" << endl;
    fitmarg->SetPrior(kPhypCPV,"flat",-TMath::Pi(),+TMath::Pi());
    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar) 
      fitmarg->SetPrior(kPhypCPVBar,"flat",-TMath::Pi(),+TMath::Pi());
  }


  // Set 23 prior
  if(file_prior_23!=""){

    fitmarg->SetPrior(kPhypMixing23,"user",hPrior_23);
    fitmarg->SetPrior(kPhypDm23sq,"user",hPrior_23);
   
    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar){
      fitmarg->SetPrior(kPhypMixing23Bar,"user",hPrior_23);
      fitmarg->SetPrior(kPhypDm23sqBar,"user",hPrior_23);
    }
  }
  else{
    cout << endl;
    cout << "Use flat prior for deltaCP" << endl;
    fitmarg->SetPrior(kPhypMixing23,"flat",0,1);
    fitmarg->SetPrior(kPhypDm23sq,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
    if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar){
      fitmarg->SetPrior(kPhypMixing23Bar,"flat",0,1);
      fitmarg->SetPrior(kPhypDm23sqBar,"flat",conv::kMinFitDm23sqBar,conv::kMaxFitDm23sqBar);
    }
  }
  
  
  // Make toys and store in binary files
    
  // Produce toys for different sub-ranges of the free parameter and store in different files.
  // Useful for running on cluster when many toys are needed (e.g. 10k x 200)

  fitmarg->MakeToysMarginalization(tweak_systematics,skexp); 
  fitmarg->WriteToys2Binary(fileid,tag);  

  
  // Draw prior and toy histograms for marginalized oscillation parameters
 
  TH1D *hPriorHisto_dcp = new TH1D(*fitmarg->GetFitParamList()->GetPriorHisto(conv::kPhypCPV));
  TH1D *hToyHisto_dcp = new TH1D(*fitmarg->GetFitParamList()->GetToyHisto(conv::kPhypCPV));

  TH2D *hPriorHisto_Mix23 = new TH2D(*fitmarg->GetFitParamList()->GetPriorHisto2D(conv::kPhypMixing23));
  TH2D *hToyHisto_Mix23 = new TH2D(*fitmarg->GetFitParamList()->GetToyHisto2D(conv::kPhypMixing23));

  TH2D *hPriorHisto_Dm23 = new TH2D(*fitmarg->GetFitParamList()->GetPriorHisto2D(conv::kPhypDm23sq));
  TH2D *hToyHisto_Dm23 = new TH2D(*fitmarg->GetFitParamList()->GetToyHisto2D(conv::kPhypDm23sq));

  
  TCanvas *cPriorHisto_dcp = new TCanvas("cPriorHisto_dcp","cPriorHisto_dcp");
  hPriorHisto_dcp->Draw();
  TCanvas *cToyHisto_dcp = new TCanvas("cToyHisto_dcp","cToyHisto_dcp");
  hToyHisto_dcp->Draw();

  TCanvas *cPriorHisto_Mix23 = new TCanvas("cPriorHisto_Mix23","cPriorHisto_Mix23");
  hPriorHisto_Mix23->Draw("colz");
  TCanvas *cToyHisto_Mix23 = new TCanvas("cToyHisto_Mix23","cToyHisto_Mix23");
  hToyHisto_Mix23->Draw("colz");

  TCanvas *cPriorHisto_Dm23 = new TCanvas("cPriorHisto_Dm23","cPriorHisto_Dm23");
  hPriorHisto_Dm23->Draw("colz");
  TCanvas *cToyHisto_Dm23 = new TCanvas("cToyHisto_Dm23","cToyHisto_Dm23");
  hToyHisto_Dm23->Draw("colz");

  char name[200];
  sprintf(name,"%s.root",tag);
  TFile *fout = new TFile(name,"RECREATE");
  hPriorHisto_dcp->Write();
  hToyHisto_dcp->Write();
  hPriorHisto_Mix23->Write();
  hToyHisto_Mix23->Write();
  hPriorHisto_Dm23->Write();
  hToyHisto_Dm23->Write();

  return;
}
