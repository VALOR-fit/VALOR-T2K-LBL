
//
// ------------------------------------------------------------------------------------------------
//
// Check the extraction of the systematic and oscillation parameters done for marginalization
// as well as the Ereco distributions.
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------

#include "Build.h"
#include "AnalysisOptions.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

void Check_toys_marginalization_readbinary
(
 //const char * config_file      = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml",
 //const char * config_file      = "analysis_configuration_run1to4_joint_Nu6.57E20POT.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_Nu6.57E20POT.xml",
 const char * config_file      = "analysis_configuration_Run1-7cPOT_theta.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_FullStat.xml",

 const char * fit_param_list        = "nuosc3f+nuisance", // dcp

 const char * sk_sample = "RHC1Rmu",   /* "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */

 const char * hypothesis           = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_oscbar_normbar, 3f_osc_normbar" */

 double       input_mixing12       = 0.304,                         
 double       input_mixing13       = 0.0243,                         
 double       input_mixing23       = 0.528,                         
 double       input_deltacp        = -1.601,                         
 double       input_dm12sq         = 7.53E-5,                         
 double       input_dm23sq         = 2.509E-3,                       
 double       input_betabar        = 1.,
 double       input_mixing23Bar    = 0.528,
 double       input_dm23sqBar      = 2.509E-3,

 const char * dmsqscheme          = "auto",   /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy           = "normal"   ,
 const char * mixparamchoice      = "sinsqtheta",
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

 const int InputBin0  = 0,
 const int InputBin1  = 0,

 ////////////////////// MARGINALIZATION PARAMETERS ///////////////////

 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_NH_10kToys_2Dtempl_ErecoVsTheta/joint_dcp_t2kreact_constrall_IH_10kToys_2Dtempl",             
 //const char * infilename = "submit_prova",
 //const char * infilename = "PROVA_SystAll",
 //const char * infilename = "submit_prova_TemplVArray_2d",
 const char * infilename = "/data/t2k/raj/studies/final_2015_joint/toys/ccpi/dcp/nh/Final_ccpi",

 const int    nfiles               = 21, //101, //201, 21, 1        // # of files to read

 const int    NtoysMarg            = 10000, //100, //10000, //1,  // # toys marginalization
 
 const char * freeparname0        = "deltacp",  // free parameter
 const double freeparmin0         = -3.14159265358979312e+00,         // minimum value of free parameter
 const double freeparmax0         = +3.14159265358979312e+00, // 10., 1         // maxinum value of free parameter
 const int    NptScan0            = 201, //9, //201, //3, //21, // 801, 2          // # pts of free parameter to scan
// 
 const char * freeparname1        = "",  // free parameter
 const double freeparmin1         = 0., // minimum value of free parameter
 const double freeparmax1         = 0., // maxinum value of free parameter
 const int    NptScan1            = 1, //201, //21,      // # pts of free parameter to scan

 //const char * freeparname0        = "mixing23Bar",  // free parameter
 //const double freeparmin0         = 0.25,         // minimum value of free parameter
 //const double freeparmax0         = 0.75, // 10., 1         // maxinum value of free parameter
 //const int    NptScan0            = 81, //9, //201, //3, //21, // 801, 2          // # pts of free parameter to scan
 
 //const char * freeparname1        = "dmsq23Bar",  // free parameter
 //const double freeparmin1         = 1.8E-3, // minimum value of free parameter
 //const double freeparmax1         = 3.2E-3, // maxinum value of free parameter
 //const int    NptScan1            = 21, //201, //21,      // # pts of free parameter to scan

 
 const bool  doMargMH             = false,     // marginalize Mass Hierarchy

 //const char *FileToyOscParam       = "../../inputs/marginalization/toys//nuebarapp_1D_NuOsc_T2Kreact_betabar_Toy_23_mix13dcp_Joint_MH_10kToys/BetaBar_4p0E20POT_Toy_23_th13dcp_Corr_Mach3.root"
 const char *FileToyOscParam       = ""
 )
{
  
  // TFile *fToyOscParam = new TFile(FileToyOscParam,"READ");

  // TH2D *hPriorHisto_th13_nh = (TH2D*) fToyOscParam->Get("hPriorHisto_Mix13_nh");
  // TH2D *hPriorHisto_th13_ih = (TH2D*) fToyOscParam->Get("hPriorHisto_Mix13_ih");
  // TH2D *hToyHisto_th13_nh = (TH2D*) fToyOscParam->Get("hToyHisto_Mix13_nh");
  // TH2D *hToyHisto_th13_ih = (TH2D*) fToyOscParam->Get("hToyHisto_Mix13_ih");
  
  // TH2D *hPriorHisto_dcp_nh = (TH2D*) fToyOscParam->Get("hPriorHisto_dcp_nh");
  // TH2D *hPriorHisto_dcp_ih = (TH2D*) fToyOscParam->Get("hPriorHisto_dcp_ih");
  // TH2D *hToyHisto_dcp_nh = (TH2D*) fToyOscParam->Get("hToyHisto_dcp_nh");
  // TH2D *hToyHisto_dcp_ih = (TH2D*) fToyOscParam->Get("hToyHisto_dcp_ih");
  
  // TH2D *hPriorHisto_Mix23_nh = (TH2D*) fToyOscParam->Get("hPriorHisto_Mix23_nh");
  // TH2D *hPriorHisto_Mix23_ih = (TH2D*) fToyOscParam->Get("hPriorHisto_Mix23_ih");
  // TH2D *hToyHisto_Mix23_nh = (TH2D*) fToyOscParam->Get("hToyHisto_Mix23_nh");
  // TH2D *hToyHisto_Mix23_ih = (TH2D*) fToyOscParam->Get("hToyHisto_Mix23_ih");

  // TH2D *hPriorHisto_Dm23_nh = (TH2D*) fToyOscParam->Get("hPriorHisto_Dm23_nh");
  // TH2D *hPriorHisto_Dm23_ih = (TH2D*) fToyOscParam->Get("hPriorHisto_Dm23_ih");
  // TH2D *hToyHisto_Dm23_nh = (TH2D*) fToyOscParam->Get("hToyHisto_Dm23_nh");
  // TH2D *hToyHisto_Dm23_ih = (TH2D*) fToyOscParam->Get("hToyHisto_Dm23_ih");
  
  // TH1D *hPriorHisto_MH = (TH1D*) fToyOscParam->Get("hPriorHisto_MH");
  // TH1D *hToyHisto_MH = (TH1D*) fToyOscParam->Get("hToyHisto_MH");

  
  // TCanvas *cPriorHisto_th13_dcp_nh = new TCanvas("cPriorHisto_th13_dcp_nh","cPriorHisto_th13_dcp_nh");
  // hPriorHisto_th13_nh->GetXaxis()->SetTitle("sin^2 #theta_{13}");
  // hPriorHisto_th13_nh->GetYaxis()->SetTitle("#delta_{CP}");
  // hPriorHisto_th13_nh->Draw("colz");

  // TCanvas *cToyHisto_th13_dcp_nh = new TCanvas("cToyHisto_th13_dcp_nh","cToyHisto_th13_dcp_nh");
  // hToyHisto_th13_nh->GetXaxis()->SetTitle("sin^2 #theta_{13}");
  // hToyHisto_th13_nh->GetYaxis()->SetTitle("#delta_{CP}");
  // hToyHisto_th13_nh->Draw("colz");

  // TCanvas *cPriorHisto_th13_dcp_ih = new TCanvas("cPriorHisto_th13_dcp_ih","cPriorHisto_th13_dcp_ih");
  // hPriorHisto_th13_ih->GetXaxis()->SetTitle("sin^2 #theta_{13}");
  // hPriorHisto_th13_ih->GetYaxis()->SetTitle("#delta_{CP}");
  // hPriorHisto_th13_ih->Draw("colz");

  // TCanvas *cToyHisto_th13_dcp_ih = new TCanvas("cToyHisto_th13_dcp_ih","cToyHisto_th13_dcp_ih");
  // hToyHisto_th13_ih->GetXaxis()->SetTitle("sin^2 #theta_{13}");
  // hToyHisto_th13_ih->GetYaxis()->SetTitle("#delta_{CP}");
  // hToyHisto_th13_ih->Draw("colz");

  // TCanvas *cPriorHisto_Mix23_nh = new TCanvas("cPriorHisto_Mix23_nh","cPriorHisto_Mix23_nh");
  // hPriorHisto_Mix23_nh->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  // hPriorHisto_Mix23_nh->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
  // hPriorHisto_Mix23_nh->Draw("colz");

  // TCanvas *cToyHisto_Mix23_nh = new TCanvas("cToyHisto_Mix23_nh","cToyHisto_Mix23_nh");
  // hToyHisto_Mix23_nh->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  // hToyHisto_Mix23_nh->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
  // hToyHisto_Mix23_nh->Draw("colz");

  // TCanvas *cPriorHisto_Mix23_ih = new TCanvas("cPriorHisto_Mix23_ih","cPriorHisto_Mix23_ih");
  // hPriorHisto_Mix23_ih->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  // hPriorHisto_Mix23_ih->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
  // hPriorHisto_Mix23_ih->Draw("colz");

  // TCanvas *cToyHisto_Mix23_ih = new TCanvas("cToyHisto_Mix23_ih","cToyHisto_Mix23_ih");
  // hToyHisto_Mix23_ih->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  // hToyHisto_Mix23_ih->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
  // hToyHisto_Mix23_ih->Draw("colz");

  //return;



  //

  FileStat_t dummy;
  
  TString opt(fit_param_list); 
  
  // Analysis configuration
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
 
  // Fit param list for fakedata sets
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); 
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used for fakedata sets: " << util::AsString(hyp) << endl;

  // Set Marginalization

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
  else if(sfreeparname0=="betaBar")     freepar0 = conv::kPhypBetaBar;
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
  else if(sfreeparname1=="betaBar")     freepar1 = conv::kPhypBetaBar;
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
    
  // Get Ereco toys
  
  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  
  int n_ereco_bins = util::NumErecoBins(s);  
  double * sumsq_per_bin = (double*)calloc(n_ereco_bins,sizeof(double));
  VALORArray * nsk_per_toy = new VALORArray(NtoysMarg);
  double * nsk_per_toy_arr = nsk_per_toy->GetArray();
  TH1D * hNsk = new TH1D("hNsk","hNsk",10000,0,1000);
  
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  

  int ipar0 = InputBin0;
  int ipar1 = InputBin1;
  
  double parvalue = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(InputBin0+1);

  cout << "par0 = " << ipar0 << " --> " << freeparname0 << " = " << parvalue << endl;
  

  TH1F* h_marg_toy;
  TH1F* h_up = new TH1F("h_up","h_up",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1F* h_down = new TH1F("h_down","h_down",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1D* h_nominal;
  TH2F* h_2d = new TH2F("h_2d","Envelope",util::NumErecoBins(s),util::ErecoBinEdges(s),2000,0,10);
  
  // Get tuned MC - with BANFF tuning

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

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
    if(input_mixing23Bar < -1) input_mixing23Bar = conv::kGlobFitSinsq2Theta23;

  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
    if(input_mixing23 < -1) input_mixing23Bar = conv::kGlobFitSinsqTheta23;

  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm23sq < -1)  input_dm23sqBar = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
  if(input_betabar < -1) input_betabar = conv::kGlobFitBetaBar;
 
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FBarNormBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypBetaBar,    input_betabar);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FNormBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypBetaBar,    input_betabar);
  } else
    if(ph.Hypothesis() == conv::kHypOsc3FBarDissap){
      ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
      ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
      ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
      ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
      ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
      ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
      ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23Bar);
      ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sqBar);
    }
 
  skexp->Generate(ph);

#ifndef __USE_2D_TEMPLATES__
  h_nominal = new TH1D(*skexp->CurrPredictedErecoSpectrum(s));
#else
  h_nominal = new TH1D(*skexp->CurrPredictedErecoVsQ2Spectrum(s)->ProjectionY());
#endif

  double nominal_total_nsk = h_nominal->Integral();

  // Loop over marginalization toys

  for(int itoy=0;itoy<NtoysMarg;itoy++){
  
#ifndef __USE_2D_TEMPLATES__
    h_marg_toy = (TH1F*)fitmarg->GetToyMargHisto(s,itoy,ipar0,ipar1);
#else
    h_marg_toy = (TH1F*)fitmarg->GetToyMargHisto2D(s,itoy,ipar0,ipar1)->ProjectionY();
#endif

    cout << "# of events = " << h_marg_toy->Integral() << endl;
    
    double NNormPOT = fitmarg->GetNormPOT(s);
    double NPOT = fitmarg->GetPOT(s);
    //cout << "NNormPOT = " << NNormPOT << endl;
    //cout << "NPOT = " << NPOT << endl;

    for(int i=1;i<=util::NumErecoBins(s);i++){
      h_2d->Fill(h_marg_toy->GetBinCenter(i),h_marg_toy->GetBinContent(i));
      //cout << h_marg_toy->GetBinContent(i) / NNormPOT << " "; // normalized to # 1x10^20 POT
    }
    //cout << endl;
	
    //store Nsk
    double nskint = h_marg_toy->Integral();
    hNsk->Fill(nskint);
    nsk_per_toy_arr[itoy] = nskint;
    
    // Store tweaks.
    for (int ibin = 0 ; ibin < n_ereco_bins ; ibin++) {
      double nom   = h_nominal->GetBinContent(ibin+1);
      double N = h_marg_toy->GetBinContent(ibin+1);
      sumsq_per_bin[ibin] += (N-nom)*(N-nom);
    }
  
  } // loop itoy
 

  // Assume symmetric errors

  // Find min/max and calculate fractional deviation from nominal
  TH1D * hErecoMax = new TH1D(*h_nominal);
  TH1D * hErecoMin = new TH1D(*h_nominal);
  hErecoMax->SetDirectory(0);
  hErecoMin->SetDirectory(0);
  hErecoMax->GetYaxis()->SetTitle("Ratio of events at 1 #sigma to nominal");
  hErecoMin->GetYaxis()->SetTitle("Ratio of events at 1 #sigma to nominal");
  TH1D * hErecoMaxEnv = new TH1D(*h_nominal);
  TH1D * hErecoMinEnv = new TH1D(*h_nominal);
  hErecoMaxEnv->SetDirectory(0);
  hErecoMinEnv->SetDirectory(0);
  for(int ibin=1; ibin <= util::NumErecoBins(s); ibin++) {
    double nom   = h_nominal->GetBinContent(ibin);
    double sumsq = 0;
    double sig = TMath::Sqrt(sumsq_per_bin[ibin-1]/NtoysMarg);
    double min = nom-sig;
    double max = nom+sig;
    if(nom>0){
      hErecoMax -> SetBinContent (ibin,  (max-nom)/nom);
      hErecoMin -> SetBinContent (ibin,  (min-nom)/nom);
      hErecoMaxEnv -> SetBinContent (ibin,  max);
      hErecoMinEnv -> SetBinContent (ibin,  min);
    }
  }//ibin




  // Assume symmetric errors wrt the average (not nominal)

  for(int i=1;i<=util::NumErecoBins(s);i++){
    TH1D* h_proj = h_2d->ProjectionY("",i,i);
    double average = average(h_proj);
    double rms = rms(h_proj,average);
    h_up->SetBinContent(i,average+ rms);
    h_down->SetBinContent(i,average-rms);  
  }
  
  //

  double dn_over_n = hNsk->GetRMS() / hNsk->GetMean();
  cout << endl;
  cout << "Nominal = " << h_nominal->Integral() << endl;
  cout << "Mean = " << hNsk->GetMean() << endl;
  cout << "RMS = " << hNsk->GetRMS() << endl;
  cout << "RMS / Mean  = " << dn_over_n << endl;
  cout << endl;

  double ave = average(hNsk);
  double rms = rms(hNsk,ave);
  cout << "ave = " << ave << " - rms = " << rms << endl;

  //if (use_new_definition_of_dn_over_n) {
  if(1){
    double sum_nsk = 0;
    for (int iexp = 0 ; iexp < NtoysMarg; iexp++) {
      sum_nsk += nsk_per_toy_arr[iexp];
    }
    double mean_nsk = sum_nsk / NtoysMarg;
    
    double sumsq = 0;
    for (int iexp = 0 ; iexp < NtoysMarg; iexp++) {
      //cerr<<nsk_per_toy_arr[iexp]<<"\t"<<mean_nsk<<endl;
      sumsq += pow(nsk_per_toy_arr[iexp]-mean_nsk,2);
    }
    
    double expected_deviation = sqrt(sumsq/NtoysMarg);
    cerr<<"Expected deviation: "<<expected_deviation<<endl;
    cerr<<"Nominal total SK: "<<nominal_total_nsk<<endl;
    
    dn_over_n = expected_deviation/nominal_total_nsk; 
  }
  cout << "RMS / Nominal = " << dn_over_n << endl;
  cout << endl;
  




  //Draw & save
  TCanvas * c = new TCanvas();
  hErecoMaxEnv->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
  hErecoMaxEnv->SetLineStyle(3);
  hErecoMinEnv->SetLineStyle(3);
  hErecoMaxEnv->SetLineColor(kRed);
  hErecoMinEnv->SetLineColor(kRed);
  hErecoMaxEnv->SetFillColor(kRed);
  hErecoMinEnv->SetFillColor(10); //white, to unfill the lower area
  hErecoMaxEnv->SetFillStyle(3244); // hatched
  hErecoMinEnv->SetFillStyle(1001); // solid (white, to unfill the lower area)
  hErecoMaxEnv->DrawCopy();
  hErecoMinEnv->DrawCopy("same");
  h_nominal->DrawCopy("same");

  
  // Draw histograms

  TCanvas * c = new TCanvas("c","",20,20,1000,1000);
  //h_2d->Draw("colz");
  h_up->Draw("");
  h_down->Draw("same");
  h_nominal->SetLineColor(kYellow);
  h_nominal->Draw("same");
      
  return;
}

double average(TH1D* h1){
  double avg(0);
  
  for(int i=1; i<=h1->GetNbinsX();i++){
    avg += h1->GetBinCenter(i)*h1->GetBinContent(i);
  }
  avg = avg/h1->GetEntries();
  
  return avg;
}

double rms(TH1D* h1,double avg){
  double rms(0);
  
  for(int i=1; i<=h1->GetNbinsX();i++){
    rms += ((h1->GetBinCenter(i) - avg) **2) * (h1->GetBinContent(i));
  }
  rms = rms/h1->GetEntries();
  rms = rms**0.5;
  return rms;
}

TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}

