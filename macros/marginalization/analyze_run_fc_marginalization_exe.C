
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
//
// Davide Sgalaberna                                                                                             
//                                                                                           
// ------------------------------------------------------------------------------------------------

#include "Build.h"
#include "AnalysisOptions.h"
#include <vector>
#include "TString.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

const double CLfrac_1sigma = 0.6827;  // 1sigma% CL
const double CLfrac_90CL   = 0.9;  // 90% CL
const double CLfrac_95CL   = 0.95;  // 95% CL
const double CLfrac_2sigma = 0.9545;  // 2sigma% CL
const double CLfrac_3sigma = 0.9973;  // 3sigma% CL

const double CRITVAL_1SIGMA_1D = 1.;
const double CRITVAL_90CL_1D =  2.71;
const double CRITVAL_95CL_1D = 3.84;
const double CRITVAL_2SIGMA_1D = 4.;
const double CRITVAL_3SIGMA_1D = 9.;

const double CRITVAL_1SIGMA_2D = 2.3;
const double CRITVAL_90CL_2D = 4.61;
const double CRITVAL_95CL_2D = 5.99;
const double CRITVAL_2SIGMA_2D = 6.18;
const double CRITVAL_3SIGMA_2D = 11.83;

void analyze_run_fc_marginalization_exe
( 
 string namepar0 = "#delta_{CP}",
 string namepar1 = "sin^{2} #theta_{13}",

 const bool UseFullChisqFit = false,
 const bool CheckSingleToy = false,
 const int  ToySel = 0,
 const bool PrintToy = false,
 const char * rootfile = "data_run1to8/SKdata_Run18.root", 

 const bool use_posterior_predictive  = false,
 const bool use_post_pred_fhc_mulike  = false,
 const bool use_post_pred_fhc_elike   = false,
 const bool use_post_pred_rhc_mulike  = false,
 const bool use_post_pred_rhc_elike   = false,
 const bool use_post_pred_fhc_cc1pilike = false, 
// const char infilename[300] = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17c_FLUXMERGE/fit_fitmarg_fc_nh_10k_run17c_fluxmerge_normal_3.142.root",
 const char infilename[300] = "/scratch/cbarry/tmp/tn_toys/FC/dcp/fc_toys/dcp_0_nh_update/dcp_0_nh_merge.root", //"/scratch/cbarry/tmp/tn_toys/FC/dcp/fc_toys/dcp_pi_ih/fitmarg_fc_dcp_pi_ih_merge.root",
 string mytag = "Run1to8_final",
 const double mydeltachisq = 7.692815042, // dcp=pi // Don't draw data deltachisq
  )
{

  cout << "Input data deltachisq: " << mydeltachisq << endl;
  const double MaxHisto = 20.;
  const double MinHisto = -5.;
  const int    NBinHisto = 250.;
  const double NBinHistoMany = 2500000.;
  const double newbinwidth = 10000.;

  char name[300];
  FitParamList * fparm  = util::NamedParameterList("nuosc3f"); 
  PhysParamList * iparm = new PhysParamList;
  vector<double> myvector;
  vector<double> myvector_nh;
  vector<double> myvector_ih;

  vector<double> vec_like_penalty;

  //TH1D *hdeltachisq_posterior = new TH1D("hdeltachisq_posterior","hdeltachisq_posterior",NBinHisto,MinHisto,MaxHisto);
  TH1D *hdeltachisq_posterior = new TH1D("hdeltachisq_posterior","hdeltachisq_posterior",NBinHistoMany,MinHisto,MaxHisto);

  //TH2D *hdeltachisq_Vs_chisqweight_posterior = new TH2D("hdeltachisq_Vs_chisqweight_posterior","hdeltachisq_Vs_chisqweight_posterior",NBinHisto,MinHisto,MaxHisto,NBinHisto,250,350);
  TH2D *hdeltachisq_Vs_chisqweight_posterior = new TH2D("hdeltachisq_Vs_chisqweight_posterior","hdeltachisq_Vs_chisqweight_posterior",NBinHisto,MinHisto,MaxHisto,NBinHisto,0,40);
  
  TH1D *hWeight_posterior = new TH1D("hWeight_posterior","hWeight_posterior",300,200,500);
  
  //TH2D *hbestfit_Vs_chisqweight_posterior = new TH2D("hbestfit_Vs_chisqweight_posterior","hbestfit_Vs_chisqweight_posterior",97,-3.15,3.15,NBinHisto,250,350);
  TH2D *hbestfit_Vs_chisqweight_posterior = new TH2D("hbestfit_Vs_chisqweight_posterior","hbestfit_Vs_chisqweight_posterior",97,-3.15,3.15,NBinHisto,0,40);
  TH1D *hbestfit_posterior_1D = new TH1D("hbestfit_posterior_1D","hbestfit_posterior_1D",97,-3.15,3.15);
  
  //TH2D *hMHbf_Vs_chisqweight_posterior = new TH2D("hMHbf_Vs_chisqweight_posterior","hMHbf_Vs_chisqweight_posterior",2,0,2,NBinHisto,250,350);
  TH2D *hMHbf_Vs_chisqweight_posterior = new TH2D("hMHbf_Vs_chisqweight_posterior","hMHbf_Vs_chisqweight_posterior",2,0,2,NBinHisto,0,40);
  TH1D *hMHbf_posterior_1D = new TH1D("hMHbf_posterior_1D","hMHbf_posterior_1D",2,0,2);

  double deltachisq0BF_nh = conv::kBadResult; 
  double deltachisq0BF_ih = conv::kBadResult; 

  TFile *filein = new TFile(infilename,"READ");
  if(!filein->IsOpen()){
    cerr << endl;
    cerr << "File " << infilename << " is not open" << endl;
    cerr << endl;
  }
  
  cout << endl;
  cout << "Reading file " << infilename << endl;
  cout << endl;

  TTree *treein = (TTree*)filein->Get("fitnt");
  if(!treein){
    cerr << endl;
    cerr << "Tree object is empty" << endl;
    cerr << endl;
  }

  // Set Tree toy/fit

  double CPUtime_Minim_nh=conv::kBadResult;
  double CPUtime_Minim_ih=conv::kBadResult;
  double CPUtime_FullChisq_nh=conv::kBadResult;
  double CPUtime_FullChisq_ih=conv::kBadResult;
  double REALtime_Minim_nh=conv::kBadResult;
  double REALtime_Minim_ih=conv::kBadResult;
  double REALtime_FullChisq_nh=conv::kBadResult;
  double REALtime_FullChisq_ih=conv::kBadResult;
  double CPUtime_end=conv::kBadResult;
  double REALtime_end=conv::kBadResult;

  bool FlagFit2D=false;
  bool IsFit2D_nh=false;
  bool IsFit2D_ih=false;
  int run=-999;
  string *tag;
  int nexpts=-999;
  int ToyID=-999;
  
  // fit status variables
  int nfits_nh=-999;
  int nfits_ih=-999;
  int nchisqcalls_nh=-999;
  int nchisqcalls_ih=-999;
  int niterations_nh=-999;
  int niterations_ih=-999;
  int fitstatus_nh=-999;
  int fitstatus_ih=-999;
  
  // true variables
  double truepar0=-999;
  double truepar1=-999;    
  double chisqtrue0=conv::kALargeNum;
  string *trueMH;
  string *trueparname0;
  string *trueparname1;
    
  // global fit output
  string *MH_glob;
  double deltachisq0BF_glob=conv::kALargeNum;  
  double chisqmin_glob=conv::kALargeNum;
  double bestfit0_glob=-999;
  double bestfit1_glob=-999;
  double bestfiterrNeg0_glob=-999;
  double bestfiterrNeg1_glob=-999;
  double bestfiterrPos0_glob=-999;
  double bestfiterrPos1_glob=-999;
  double bestfiterrSym0_glob=-999;
  double bestfiterrSym1_glob=-999;
  double bestfitFullChisq0_glob=-999;
  double bestfitFullChisq1_glob=-999;
  double chisqminFullChisq_glob=conv::kALargeNum;

  // normal hierarchy fit output
  double chisqmin_nh=conv::kALargeNum;
  double bestfit0_nh=-999;
  double bestfit1_nh=-999;
  double bestfiterrNeg0_nh=-999;
  double bestfiterrNeg1_nh=-999;
  double bestfiterrPos0_nh=-999;
  double bestfiterrPos1_nh=-999;
  double bestfiterrSym0_nh=-999;
  double bestfiterrSym1_nh=-999;
  double bestfitFullChisq0_nh=-999;
  double bestfitFullChisq1_nh=-999;
  double chisqminFullChisq_nh=conv::kALargeNum;

  // inverted hierarchy fit output
  double chisqmin_ih=conv::kALargeNum;
  double bestfit0_ih=-999;
  double bestfit1_ih=-999;
  double bestfiterrNeg0_ih=-999;
  double bestfiterrNeg1_ih=-999;
  double bestfiterrPos0_ih=-999;
  double bestfiterrPos1_ih=-999;
  double bestfiterrSym0_ih=-999;
  double bestfiterrSym1_ih=-999;
  double bestfitFullChisq0_ih=-999;
  double bestfitFullChisq1_ih=-999;
  double chisqminFullChisq_ih=conv::kALargeNum;

  // histograms
  TH2D *hchisq2D_out_nh = new TH2D("hchisq2D_out_nh","hchisq2D_out_nh",1,0,1,1,0,1);
  TH2D *hchisq2D_out_ih = new TH2D("hchisq2D_out_ih","hchisq2D_out_ih",1,0,1,1,0,1);
  TH1D *hchisq1D_out_nh = new TH1D("hchisq1D_out_nh","hchisq1D_out_nh",1,0,1);
  TH1D *hchisq1D_out_ih = new TH1D("hchisq1D_out_ih","hchisq1D_out_ih",1,0,1);
  TH2D *hdeltachisq2D_out_nh = new TH2D("hdeltachisq2D_out_nh","hdeltachisq2D_out_nh",1,0,1,1,0,1);
  TH2D *hdeltachisq2D_out_ih = new TH2D("hdeltachisq2D_out_ih","hdeltachisq2D_out_ih",1,0,1,1,0,1);
  TH1D *hdeltachisq1D_out_nh = new TH1D("hdeltachisq1D_out_nh","hdeltachisq1D_out_nh",1,0,1);
  TH1D *hdeltachisq1D_out_ih = new TH1D("hdeltachisq1D_out_ih","hdeltachisq1D_out_ih",1,0,1);
  TH2D *hstep2D_out_nh = new TH2D("hstep2D_out_nh","hstep2D_out_nh",1,0,1,1,0,1);
  TH2D *hstep2D_out_ih = new TH2D("hstep2D_out_ih","hstep2D_out_ih",1,0,1,1,0,1);
  TH1D *hstep1D_out_nh = new TH1D("hstep1D_out_nh","hstep1D_out_nh",1,0,1);
  TH1D *hstep1D_out_ih = new TH1D("hstep1D_out_ih","hstep1D_out_ih",1,0,1);

#ifndef __USE_2D_TEMPLATES__
  TH1D *hErecoToy_1Rmu = new TH1D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1);
  TH1D *hErecoToy_1Re = new TH1D("hErecoToy_1Re","hErecoToy_1Re",1,0,1);
  TH1D *hErecoToy_RHC1Rmu = new TH1D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1);
  TH1D *hErecoToy_RHC1Re = new TH1D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1);
  TH1D *hErecoToy_MultiRe = new TH1D("hErecoToy_MultiRe","hErecoToy_MultiRe",1,0,1);
  TH1D *hErecoToy_sample[conv::kns] = {};
  TH1D *hErecoExp_1Rmu = new TH1D("hErecoExp_1Rmu","hErecoExp_1Rmu",1,0,1);
  TH1D *hErecoExp_1Re = new TH1D("hErecoExp_1Re","hErecoExp_1Re",1,0,1);
  TH1D *hErecoExp_RHC1Rmu = new TH1D("hErecoExp_RHC1Rmu","hErecoExp_RHC1Rmu",1,0,1);
  TH1D *hErecoExp_RHC1Re = new TH1D("hErecoExp_RHC1Re","hErecoExp_RHC1Re",1,0,1);
  TH1D *hErecoExp_MultiRe = new TH1D("hErecoExp_MultiRe","hErecoExp_MultiRe",1,0,1);
  TH1D *hErecoExp_sample[conv::kns] = {};
#else
  TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_MultiRe = new TH2D("hErecoToy_MultiRe","hErecoToy_MultiRe",1,0,1,1,0,1);
  TH2D *hErecoToy_sample[conv::kns] = {};
  TH2D *hErecoExp_1Rmu = new TH2D("hErecoExp2D_1Rmu","hErecoExp2D_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_1Re = new TH2D("hErecoExp2D_1Re","hErecoExp2D_1Re",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Rmu = new TH2D("hErecoExp2D_RHC1Rmu","hErecoExp2D_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_MultiRe = new TH2D("hErecoExp2D_MultiRe","hErecoExp2D_MultiRe",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Re = new TH2D("hErecoExp2D_RHC1Re","hErecoExp2D_RHC1Re",1,0,1,1,0,1); 
  cout << "2D templates" << endl;
#endif


  // Create tree to store outputs
  
  treein->SetBranchAddress("best_fit_params",&fparm);
  treein->SetBranchAddress("input_params",&iparm);
  
  treein->SetBranchAddress("CPUtime_Minim_nh",&CPUtime_Minim_nh);
  treein->SetBranchAddress("CPUtime_Minim_ih",&CPUtime_Minim_ih);
  treein->SetBranchAddress("CPUtime_FullChisq_nh",&CPUtime_FullChisq_nh);
  treein->SetBranchAddress("CPUtime_FullChisq_ih",&CPUtime_FullChisq_ih);
  treein->SetBranchAddress("REALtime_Minim_nh",&REALtime_Minim_nh);
  treein->SetBranchAddress("REALtime_Minim_ih",&REALtime_Minim_ih);
  treein->SetBranchAddress("REALtime_FullChisq_nh",&REALtime_FullChisq_nh);
  treein->SetBranchAddress("REALtime_FullChisq_ih",&REALtime_FullChisq_ih);
  treein->SetBranchAddress("CPUtime_end",&CPUtime_end);
  treein->SetBranchAddress("REALtime_end",&REALtime_end);

  treein->SetBranchAddress("ToyID",&ToyID);
  treein->SetBranchAddress("Job",&run);
  treein->SetBranchAddress("Tag",&tag);
  treein->SetBranchAddress("nexpts",&nexpts);

  //fit status
  treein->SetBranchAddress("IsFit2D_nh",&IsFit2D_nh);
  treein->SetBranchAddress("IsFit2D_ih",&IsFit2D_ih);
  treein->SetBranchAddress("fitstatus_nh",&fitstatus_nh);
  treein->SetBranchAddress("fitstatus_ih",&fitstatus_ih);
  treein->SetBranchAddress("nfits_nh",&nfits_nh);
  treein->SetBranchAddress("nfits_ih",&nfits_ih);
  treein->SetBranchAddress("nchisqcalls_nh",&nchisqcalls_nh);
  treein->SetBranchAddress("nchisqcalls_ih",&nchisqcalls_ih);
  treein->SetBranchAddress("niterations_nh",&niterations_nh);
  treein->SetBranchAddress("niterations_ih",&niterations_ih);

  //true variables
  treein->SetBranchAddress("truepar0", &truepar0 );
  treein->SetBranchAddress("truepar1", &truepar1 );
  treein->SetBranchAddress("trueMH",   &trueMH   );
  treein->SetBranchAddress("trueparname0",   &trueparname0   );
  treein->SetBranchAddress("trueparname1",   &trueparname1   );
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  
  //global fit output
  treein->SetBranchAddress("deltachisq0BF_glob",&deltachisq0BF_glob);
  treein->SetBranchAddress("MH_glob",   &MH_glob   );
  treein->SetBranchAddress("chisqmin_glob",&chisqmin_glob);
  treein->SetBranchAddress("bestfit0_glob", &bestfit0_glob );
  treein->SetBranchAddress("bestfit1_glob", &bestfit1_glob );
  treein->SetBranchAddress("bestfiterrNeg0_glob",&bestfiterrNeg0_glob);
  treein->SetBranchAddress("bestfiterrNeg1_glob",&bestfiterrNeg1_glob);
  treein->SetBranchAddress("bestfiterrPos0_glob",&bestfiterrPos0_glob);
  treein->SetBranchAddress("bestfiterrPos1_glob",&bestfiterrPos1_glob);
  treein->SetBranchAddress("bestfiterrSym0_glob",&bestfiterrSym0_glob);
  treein->SetBranchAddress("bestfiterrSym1_glob",&bestfiterrSym1_glob);
  treein->SetBranchAddress("bestfitFullChisq0_glob", &bestfitFullChisq0_glob);
  treein->SetBranchAddress("chisqminFullChisq_glob",&chisqminFullChisq_glob);

  //normal hierarchy fit output
  treein->SetBranchAddress("chisqmin_nh",&chisqmin_nh);
  treein->SetBranchAddress("bestfit0_nh", &bestfit0_nh );
  treein->SetBranchAddress("bestfit1_nh", &bestfit1_nh );
  treein->SetBranchAddress("bestfiterrNeg0_nh",&bestfiterrNeg0_nh);
  treein->SetBranchAddress("bestfiterrNeg1_nh",&bestfiterrNeg1_nh);
  treein->SetBranchAddress("bestfiterrPos0_nh",&bestfiterrPos0_nh);
  treein->SetBranchAddress("bestfiterrPos1_nh",&bestfiterrPos1_nh);
  treein->SetBranchAddress("bestfiterrSym0_nh",&bestfiterrSym0_nh);
  treein->SetBranchAddress("bestfiterrSym1_nh",&bestfiterrSym1_nh);
  treein->SetBranchAddress("bestfitFullChisq0_nh", &bestfitFullChisq0_nh );
  treein->SetBranchAddress("chisqminFullChisq_nh",&chisqminFullChisq_nh);

  //inverted hierarchy fit output
  treein->SetBranchAddress("chisqmin_ih",&chisqmin_ih);
  treein->SetBranchAddress("bestfit0_ih", &bestfit0_ih );
  treein->SetBranchAddress("bestfit1_ih", &bestfit1_ih );
  treein->SetBranchAddress("bestfiterrNeg0_ih",&bestfiterrNeg0_ih);
  treein->SetBranchAddress("bestfiterrNeg1_ih",&bestfiterrNeg1_ih);
  treein->SetBranchAddress("bestfiterrPos0_ih",&bestfiterrPos0_ih);
  treein->SetBranchAddress("bestfiterrPos1_ih",&bestfiterrPos1_ih);
  treein->SetBranchAddress("bestfiterrSym0_ih",&bestfiterrSym0_ih);
  treein->SetBranchAddress("bestfiterrSym1_ih",&bestfiterrSym1_ih);
  treein->SetBranchAddress("bestfitFullChisq0_ih", &bestfitFullChisq0_ih );
  treein->SetBranchAddress("chisqminFullChisq_ih",&chisqminFullChisq_ih);
    
  // histograms
  treein->SetBranchAddress("hchisq1D_nh",     &hchisq1D_out_nh      );
  treein->SetBranchAddress("hchisq1D_ih",     &hchisq1D_out_ih      );
  treein->SetBranchAddress("hstep1D_nh",      &hstep1D_out_nh       );
  treein->SetBranchAddress("hstep1D_ih",      &hstep1D_out_ih       );
  treein->SetBranchAddress("hdeltachisq1D_nh",&hdeltachisq1D_out_nh );
  treein->SetBranchAddress("hdeltachisq1D_ih",&hdeltachisq1D_out_ih );
  treein->SetBranchAddress("hchisq2D_nh",     &hchisq2D_out_nh      );
  treein->SetBranchAddress("hchisq2D_ih",     &hchisq2D_out_ih      );
  treein->SetBranchAddress("hdeltachisq2D_nh",&hdeltachisq2D_out_nh );
  treein->SetBranchAddress("hdeltachisq2D_ih",&hdeltachisq2D_out_ih );
  treein->SetBranchAddress("hstep2D_nh",      &hstep2D_out_nh       );
  treein->SetBranchAddress("hstep2D_ih",      &hstep2D_out_ih       );

#ifndef __USE_2D_TEMPLATES__
  treein->SetBranchAddress("hErecoToy_1Rmu",     &hErecoToy_1Rmu     );
  treein->SetBranchAddress("hErecoToy_1Re",      &hErecoToy_1Re      );
  treein->SetBranchAddress("hErecoToy_RHC1Rmu",  &hErecoToy_RHC1Rmu  );
  treein->SetBranchAddress("hErecoToy_RHC1Re",   &hErecoToy_RHC1Re   );
  treein->SetBranchAddress("hErecoExp_1Rmu",     &hErecoExp_1Rmu     );
  treein->SetBranchAddress("hErecoExp_1Re",      &hErecoExp_1Re      );
  treein->SetBranchAddress("hErecoExp_RHC1Rmu",  &hErecoExp_RHC1Rmu  );
  treein->SetBranchAddress("hErecoExp_RHC1Re",   &hErecoExp_RHC1Re   );
  treein->SetBranchAddress("hErecoExp_MultiRe",  &hErecoExp_MultiRe  );  
#else  
  cout << "HERE" << endl;
  treein->SetBranchAddress("hErecoToy_1Rmu",     &hErecoToy_1Rmu     );
  treein->SetBranchAddress("hErecoToy_1Re",      &hErecoToy_1Re      );
  treein->SetBranchAddress("hErecoToy_RHC1Rmu",  &hErecoToy_RHC1Rmu  );
  treein->SetBranchAddress("hErecoToy_RHC1Re",   &hErecoToy_RHC1Re   );
  treein->SetBranchAddress("hErecoToy_MultiRe",  &hErecoToy_MultiRe  );
  treein->SetBranchAddress("hErecoExp_1Rmu",     &hErecoExp_1Rmu     );
  treein->SetBranchAddress("hErecoExp_1Re",      &hErecoExp_1Re      );
  treein->SetBranchAddress("hErecoExp_RHC1Rmu",  &hErecoExp_RHC1Rmu  );
  treein->SetBranchAddress("hErecoExp_RHC1Re",   &hErecoExp_RHC1Re   );
  treein->SetBranchAddress("hErecoExp_MultiRe",  &hErecoExp_MultiRe  );

#endif



  // Data files

  TH2D *h_data_1rmu;
  TH2D *h_data_1re;
  TH2D *h_data_multire;
  TH2D *h_data_rhc1rmu;
  TH2D *h_data_rhc1re;
  TH2D *h_data_multire;  
  
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));
  SuperkObservation * skobs_fakedata = new SuperkObservation;
  TFile *infile;
  if(use_posterior_predictive){
    cout << endl;
    cout << "PERFORM A POSTERIOR PREDICTIVE ANALYSIS:" << endl;
    cout << " - use root file " << rootfile2 << endl;
    cout << endl;
    infile = new TFile(rootfile2.Data(),"READ");
    if(!infile->IsOpen()){
      cout << endl;
      cout << "The data root file doesn't exist: don't draw data spectra" << endl;
      cout << endl;
    }
    infile->Close();
    
    conv::SuperkSample_t s_1rmu    = util::SKSampleFromString("1Rmu"); 
    conv::SuperkSample_t s_1re     = util::SKSampleFromString("1Re");
    conv::SuperkSample_t s_multire = util::SKSampleFromString("MultiRe");
    conv::SuperkSample_t s_rhc1rmu = util::SKSampleFromString("RHC1Rmu"); 
    conv::SuperkSample_t s_rhc1re  = util::SKSampleFromString("RHC1Re");
    conv::SuperkSample_t s_multire = util::SKSampleFromString("MultiRe"); 
    
    //int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
    bool using_xml_file = false;
    int iexperiment = 0;
    skobs_fakedata->Read(rootfile2.Data(),iexperiment,using_xml_file);

    h_data_1rmu    = (TH2D*) skobs_fakedata->TwoDSpectrum(s_1rmu)   ->Clone("h_data_1rmu");
    h_data_1re     = (TH2D*) skobs_fakedata->TwoDSpectrum(s_1re)    ->Clone("h_data_1re");
    h_data_multire = (TH2D*) skobs_fakedata->TwoDSpectrum(s_1re)    ->Clone("h_data_multire");
    h_data_rhc1rmu = (TH2D*) skobs_fakedata->TwoDSpectrum(s_rhc1rmu)->Clone("h_data_rhc1rmu");
    h_data_rhc1re  = (TH2D*) skobs_fakedata->TwoDSpectrum(s_rhc1re) ->Clone("h_data_rhc1re");
    h_data_multire = (TH2D*) skobs_fakedata->TwoDSpectrum(s_multire)->Clone("h_data_multire");     
    
    cout << "# of events 1Rmu = "    << h_data_1rmu->Integral()    << endl;
    cout << "# of events 1Re = "     << h_data_1re->Integral()     << endl;
    cout << "# of events MultiRe = " << h_data_multire->Integral() << endl;
    cout << "# of events RHC1Rmu = " << h_data_rhc1rmu->Integral() << endl;
    cout << "# of events RHC1Re = "  << h_data_rhc1re->Integral()  << endl;
  }
  

  ///////////// Declare histograms //////////////

  TH1D *htrue_MH = new TH1D("htrue_MH","htrue_MH",2,0,2);
  TH1D *htrue_1D = new TH1D("htrue_1D","htrue_1D",97,-3.15,3.15);
  TH2D *htrue_2D = new TH2D("htrue_2D","htrue_2D",97,0,0.3,97,-3.15,3.15);

  TH1D *hbestfitFullChisq_1D = new TH1D("hbestfitFullChisq_1D","hbestfitFullChisq_1D",97,-3.15,3.15);
  TH2D *hbestfitFullChisq_2D = new TH2D("hbestfitFullChisq_2D","hbestfitFullChisq_2D",97,0,0.3,97,-3.15,3.15);

  TH1D *hbestfit_1D = new TH1D("hbestfit_1D","hbestfit_1D",97,-3.15,3.15);
  TH1D *hbestfit_1D_nh = new TH1D("hbestfit_1D_nh","hbestfit_1D_nh",97,-3.15,3.15);
  TH1D *hbestfit_1D_ih = new TH1D("hbestfit_1D_ih","hbestfit_1D_ih",97,-3.15,3.15);

  TH1D *hbestfit_MH = new TH1D("hbestfit_MH","hbestfit_MH",2,0,2);
  
  TH2D *hbestfit_2D = new TH2D("hbestfit_2D","hbestfit_2D",97,0,0.3,97,-3.15,3.15);
  TH2D *hbestfit_2D_nh = new TH2D("hbestfit_2D_nh","hbestfit_2D_nh",97,0,0.3,97,-3.15,3.15);
  TH2D *hbestfit_2D_ih = new TH2D("hbestfit_2D_ih","hbestfit_2D_ih",97,0,0.3,97,-3.15,3.15);

  TH1D *hdeltachisq = new TH1D("hdeltachisq","hdeltachisq",NBinHisto,MinHisto,MaxHisto);
  //TH1D *hdeltachisq = new TH1D("hdeltachisq","hdeltachisq",100,-5,20);
  TH1D *hdeltachisq_nh = new TH1D("hdeltachisq_nh","hdeltachisq_nh",100,-5,20);
  TH1D *hdeltachisq_ih = new TH1D("hdeltachisq_ih","hdeltachisq_ih",100,-5,20);


  TH1D *hchisqmin = new TH1D("hchisqmin","hchisqmin",200,-10,400); 
  TH1D *hchisqminDiff = new TH1D("hchisqminDiff","hchisqminDiff",100,-0.3,3.3);
  TH1D *hchisqminFullChisq = new TH1D("hchisqminFullChisq","hchisqminFullChisq",100,-10,100);
  TH1D *hnchisqcalls = new TH1D("hnchisqcalls","hnchisqcalls",100,0,700);
  TH1D *hfitstatus = new TH1D("hfitstatus","hfitstatus",5,0,5);

  TH1D *hNevents_fhc_mulike = new TH1D("hNevents_fhc_mulike","hNevents_fhc_mulike",300,0,500);
  TH1D *hNevents_fhc_elike  = new TH1D("hNevents_fhc_elike", "hNevents_fhc_elike" ,100,0,250);
  TH1D *hNevents_fhc_cc1pilike  = new TH1D("hNevents_fhc_cc1pilike", "hNevents_fhc_cc1pilike" ,100,0,250);
  TH1D *hNevents_rhc_mulike = new TH1D("hNevents_rhc_mulike","hNevents_rhc_mulike",100,0,500);
  TH1D *hNevents_rhc_elike  = new TH1D("hNevents_rhc_elike", "hNevents_rhc_elike" ,100x,0,250);

  TH1D *hNevents_exp_fhc_mulike = new TH1D("hNevents_exp_fhc_mulike","hNevents_exp_fhc_mulike",100,0,500);
  TH1D *hNevents_exp_fhc_elike  = new TH1D("hNevents_exp_fhc_elike", "hNevents_exp_fhc_elike" ,100,0,250);
  TH1D *hNevents_exp_fhc_cc1pilike  = new TH1D("hNevents_exp_fhc_cc1pilike", "hNevents_exp_fhc_cc1pilike" ,100,0,250);
  TH1D *hNevents_exp_rhc_mulike = new TH1D("hNevents_exp_rhc_mulike","hNevents_exp_rhc_mulike",100,0,500);
  TH1D *hNevents_exp_rhc_elike  = new TH1D("hNevents_exp_rhc_elike", "hNevents_exp_rhc_elike" ,100,0,250);
 
  ///////////////////////////////////////////////
  
  // Check if 1D or 2D fit
  bool IsFit2D;
  treein->GetEntry(0);
  IsFit2D = FlagFit2D;  
  if(IsFit2D){
    cout << endl;
    cout << "The fit is 2D" << endl;
    cout << endl;
  }
  else{
    cout << endl;
    cout << "The fit is 1D" << endl;
    cout << endl;
  }

  // loop over the tree  
  
  int Nin_1sigma = 0;
  int Nin_90CL   = 0;
  int Nin_95CL   = 0;
  int Nin_2sigma = 0;
  int Nin_3sigma = 0;

  int Nin_1sigma_nh = 0;
  int Nin_90CL_nh   = 0;
  int Nin_95CL_nh   = 0;
  int Nin_2sigma_nh = 0;
  int Nin_3sigma_nh = 0;
  
  int Nin_1sigma_ih = 0;
  int Nin_90CL_ih   = 0;
  int Nin_95CL_ih   = 0;
  int Nin_2sigma_ih = 0;
  int Nin_3sigma_ih = 0;

  int NdiffChisq = 0;
  
  int Nentries = treein->GetEntries();
  //Nentries = 5000;

  cout << endl;
  cout << "Total # of toys = " << Nentries << endl;
  cout << endl;

  int empty_nev_exp_rhc1re   = 0;
  int empty_nev_exp_1re      = 0;
  int empty_nev_exp_rhc1rmu  = 0;
  int empty_nev_exp_1rmu     = 0;
  int empty_nev_exp_multire  = 0;


  
  for(int ientry=0;ientry<Nentries;ientry++){
    
    if(CheckSingleToy){ 
      Nentries = 1;
      ientry = ToySel;
    }

    treein->GetEntry(ientry);
 
    if(ientry%1000==0){
      cout << "Entry " << ientry << endl;
    }

    //
    
    //if(deltachisq0BF_glob<2) continue;
    //if(*MH_glob=="normal") continue; // plots a sinistra
    //if(*MH_glob=="inverted") continue; // plots a sinistra
    //if(bestfit0_glob > 0) continue;
    //if(bestfit0_glob < 0) continue;

    //

    double nevents_rhc1re = hErecoToy_RHC1Re->Integral();    
    double nevents_1re = hErecoToy_1Re->Integral();    
    double nevents_rhc1rmu = hErecoToy_RHC1Rmu->Integral();    
    double nevents_1rmu = hErecoToy_1Rmu->Integral();    
    double nevents_multire = hErecoToy_MultiRe->Integral();

    double nevents_exp_rhc1re = hErecoExp_RHC1Re->Integral();    
    double nevents_exp_1re = hErecoExp_1Re->Integral();    
    double nevents_exp_rhc1rmu = hErecoExp_RHC1Rmu->Integral();    
    double nevents_exp_1rmu = hErecoExp_1Rmu->Integral();
    double nevents_exp_multire = hErecoExp_MultiRe->Integral();

    if(ientry%1000==0){
      cout << "Toy 1K 1Rmu Exp 1Rmu " << nevents_1rmu << " " << nevents_exp_1rmu << endl;
    }
   
    if(
      (nevents_rhc1re<0) ||
      (nevents_1re<0) ||
      (nevents_rhc1rmu<0) ||
      (nevents_1rmu<0) ||
      (nevents_multire<0)
      ){
     cerr << endl;
     cerr << "ToyID: " << ToyID << endl;
     cerr << "# of events RHC1Re  = " << nevents_rhc1re  << " < 0" << endl;
     cerr << "# of events 1Re     = " << nevents_1re     << " < 0" << endl;
     cerr << "# of events RHC1Rmu = " << nevents_rhc1rmu << " < 0" << endl;
     cerr << "# of events 1Rmu    = " << nevents_1rmu    << " < 0" << endl;
     cerr << endl;
     exit(1);
    }
    
    hNevents_rhc_elike->Fill(nevents_rhc1re);    
    hNevents_fhc_elike->Fill(nevents_1re);    
    hNevents_fhc_cc1pilike->Fill(nevents_multire);
    hNevents_rhc_mulike->Fill(nevents_rhc1rmu);    
    hNevents_fhc_mulike->Fill(nevents_1rmu);

    hNevents_exp_rhc_elike->Fill(nevents_exp_rhc1re);    
    hNevents_exp_fhc_elike->Fill(nevents_exp_1re);    
    hNevents_exp_fhc_cc1pilike->Fill(nevents_exp_multire);
    hNevents_exp_rhc_mulike->Fill(nevents_exp_rhc1rmu);    
    hNevents_exp_fhc_mulike->Fill(nevents_exp_1rmu);

     // cout << "ToyID: " << ToyID << endl;
     // cout << "# of events RHC1Re  = " << nevents_rhc1re  << endl;
     // cout << "# of events 1Re     = " << nevents_1re     << endl;
     // cout << "# of events RHC1Rmu = " << nevents_rhc1rmu << endl;
     // cout << "# of events 1Rmu    = " << nevents_1rmu    << endl;
     // cout << "# of events MultiRe = " << nevents_multire << endl;

     // cout << "-------------------------" << endl; 
     // cout << "Nevents Exp: " << endl; 
     // cout << "# of events RHC1Re  = " << nevents_exp_rhc1re   << endl;
     // cout << "# of events 1Re     = " << nevents_exp_1re      << endl;
     // cout << "# of events RHC1Rmu = " << nevents_exp_rhc1rmu  << endl;
     // cout << "# of events 1Rmu    = " << nevents_exp_1rmu     << endl;
     // cout << "# of events MultiRe = " << nevents_exp_multire  << endl;

     // cout << endl;

         
    if (nevents_exp_rhc1re  <= 0)  empty_nev_exp_rhc1re++ ;
    if (nevents_exp_1re     <= 0)  empty_nev_exp_1re++    ;
    if (nevents_exp_rhc1rmu <= 0)  empty_nev_exp_rhc1rmu++;
    if (nevents_exp_1rmu    <= 0)  empty_nev_exp_1rmu++   ;
    if (nevents_exp_multire <= 0)  empty_nev_exp_multire++;
                                  
    
    //




    if(ToyID==0) cout << "ToyID = " << ToyID << endl;

    
     // Print output
    
    if(PrintToy){
      //cout << *fparm << endl;
      cout << std::setprecision(9);
      cout << "/////////////////" << endl;
      cout << endl;
      cout << "Toy " << ientry << endl;
      cout << "ToyID " << ToyID << endl;
      cout << endl;      
      cout << "/// True parameters ///" << endl;
      cout << " - truepar0: " << *trueparname0 << " = " << truepar0 << endl;
      cout << " - truepar1: " << *trueparname1 << " = " << truepar1 << endl;
      cout << " - Mass Hierarchy: " << *trueMH << endl;
      cout << " - true chisq: " << chisqtrue0 << endl;
      cout << endl;      
      cout << "/// Normal Hierarchy ///" << endl;
      cout << " - Minimum Marginalized Chisq = " << chisqmin_nh << endl;        
      cout << " - bestfit 0 = " << bestfit0_nh << endl;
      if(IsFit2D) cout << " - bestfit 1 = " << bestfit1_nh << endl;
      cout << " - # of calls of CalcChisqMargPt(): " << nchisqcalls_nh << endl;
      cout << "Minimum of full Chisq histogram:" << endl;
      cout << " - chisqmin = " << chisqminFullChisq_nh << endl; 
      cout << " - bestfit0 = " << bestfitFullChisq0_nh << endl; 
      if(IsFit2D) cout << " - bestfit1 = " << bestfitFullChisq1_nh << endl; 
      cout << endl;      
      cout << "/// Inverted Hierarchy ///" << endl;
      cout << " - Minimum Marginalized Chisq = " << chisqmin_ih << endl;        
      cout << " - bestfit 0 = " << bestfit0_ih << endl;
      if(IsFit2D) cout << " - bestfit 1 = " << bestfit1_ih << endl;
      cout << " - # of calls of CalcChisqMargPt(): " << nchisqcalls_ih << endl;
      cout << "Minimum of full Chisq histogram:" << endl;
      cout << " - chisqmin = " << chisqminFullChisq_ih << endl; 
      cout << " - bestfit0 = " << bestfitFullChisq0_ih << endl; 
      if(IsFit2D) cout << " - bestfit1 = " << bestfitFullChisq1_ih << endl; 
      cout << endl;      
      cout << "/// Global minimum ///" << endl;
      cout << " - DeltaChisq = " << deltachisq0BF_glob << endl;
      cout << " - Minimum Marginalized Chisq = " << chisqmin_glob << endl;  
      cout << " - Mass Hierarchy: " << *MH_glob << endl;
      cout << " - bestfit 0 = " << bestfit0_glob << endl;
      if(IsFit2D) cout << " - bestfit 1 = " << bestfit1_glob << endl;
      cout << "Minimum of full Chisq histogram:" << endl;
      cout << " - chisqmin = " << chisqminFullChisq_glob << endl; 
      cout << " - bestfit0 = " << bestfitFullChisq0_glob << endl; 
      if(IsFit2D) cout << " - bestfit1 = " << bestfitFullChisq1_glob << endl; 
      cout << endl;
    }
   
    // True variables

    if(*trueMH=="normal")        htrue_MH->Fill(0.5);
    else if(*trueMH=="inverted") htrue_MH->Fill(1.5);
    else{
      cerr << "No true Mass Hierarchy" << endl;
      exit(1);
    }
    if(IsFit2D) htrue_2D->Fill(truepar0,truepar1);
    else        htrue_1D->Fill(truepar0);

    // Find bestfit chisq and true chisq
        
    double deltachisq0BF_used = deltachisq0BF_glob;

    // cout << "deltachisq0BF_used = " << deltachisq0BF_used;
    deltachisq0BF_used = chisqtrue0-chisqmin_glob;
    // cout << "  " << deltachisq0BF_used;
    deltachisq0BF_used = chisqtrue0-chisqminFullChisq_glob;
    // cout << "  " << deltachisq0BF_used << endl;
    if(!UseFullChisqFit) deltachisq0BF_used = chisqtrue0-chisqmin_glob;
    else                 deltachisq0BF_used = chisqtrue0-chisqminFullChisq_glob;

    // Print chisq 
    if(*trueMH=="inverted"){   
      int npt_chisq_x = hchisq1D_out_ih->GetNbinsX();
      for(int ix=0;ix<npt_chisq_x;ix++){
        // cout << hchisq1D_out_ih->GetBinContent(ix+1) << " ";
      }
    }
    if(*trueMH=="normal"){
      int npt_chisq_x = hchisq1D_out_nh->GetNbinsX();
      for(int ix=0;ix<npt_chisq_x;ix++){
        // cout << hchisq1D_out_nh->GetBinContent(ix+1) << " ";
      }
    }

    // Fill deltaChisq histogram
    
    if(*trueMH=="normal"){
      deltachisq0BF_nh = chisqtrue0 - chisqmin_nh;
      hdeltachisq_nh->Fill(deltachisq0BF_nh);
      myvector_nh.push_back(deltachisq0BF_nh);
      hbestfit_2D_nh->Fill(bestfit0_nh,bestfit1_nh);
      hbestfit_1D_nh->Fill(bestfit0_nh);
      // Fill variables for coverage study
      if(deltachisq0BF_nh <= CRITVAL_1SIGMA_1D ) Nin_1sigma_nh++;
      if(deltachisq0BF_nh <= CRITVAL_90CL_1D   ) Nin_90CL_nh++;
      if(deltachisq0BF_nh <= CRITVAL_95CL_1D   ) Nin_95CL_nh++;
      if(deltachisq0BF_nh <= CRITVAL_2SIGMA_1D ) Nin_2sigma_nh++;
      if(deltachisq0BF_nh <= CRITVAL_3SIGMA_1D ) Nin_3sigma_nh++;
    }
    else if(*trueMH=="inverted"){
      deltachisq0BF_ih = chisqtrue0 - chisqmin_ih;
      hdeltachisq_ih->Fill(deltachisq0BF_ih);
      myvector_ih.push_back(deltachisq0BF_ih);
      hbestfit_2D_ih->Fill(bestfit0_ih,bestfit1_ih);
      hbestfit_1D_ih->Fill(bestfit0_ih);
      // Fill variables for coverage study
      if(deltachisq0BF_ih <= CRITVAL_1SIGMA_1D ) Nin_1sigma_ih++;
      if(deltachisq0BF_ih <= CRITVAL_90CL_1D   ) Nin_90CL_ih++;
      if(deltachisq0BF_ih <= CRITVAL_95CL_1D   ) Nin_95CL_ih++;
      if(deltachisq0BF_ih <= CRITVAL_2SIGMA_1D ) Nin_2sigma_ih++;
      if(deltachisq0BF_ih <= CRITVAL_3SIGMA_1D ) Nin_3sigma_ih++;
    }
    else{
      cout << "WARNING:" << endl;
      cout << "ToyID = " << ToyID << endl;
      cout << "True Mass Hierarchy is not stored" << endl;
    }
    
    // Global minimum
    
    hdeltachisq->Fill(deltachisq0BF_used);
    myvector.push_back(deltachisq0BF_used);
    hchisqmin->Fill(chisqmin_glob);
    hbestfit_2D->Fill(bestfit0_glob,bestfit1_glob);
    hbestfit_1D->Fill(bestfit0_glob);
    if(*MH_glob=="normal") hbestfit_MH->Fill(0.5);
    else if(*MH_glob=="inverted") hbestfit_MH->Fill(1.5);
    // Fill variables for coverage study
    if(deltachisq0BF_used <= CRITVAL_1SIGMA_1D ) Nin_1sigma++;
    if(deltachisq0BF_used <= CRITVAL_90CL_1D   ) Nin_90CL++;
    if(deltachisq0BF_used <= CRITVAL_95CL_1D   ) Nin_95CL++;
    if(deltachisq0BF_used <= CRITVAL_2SIGMA_1D ) Nin_2sigma++;
    if(deltachisq0BF_used <= CRITVAL_3SIGMA_1D ) Nin_3sigma++;
        
    // Get Likelihood weights --> penalty term
    long double toy_weight = 1.;
    long double toy_chisq = 0;
    if(use_posterior_predictive){
      if(use_post_pred_fhc_mulike) 
	toy_chisq += calc_chisq(h_data_1rmu,   hErecoExp_1Rmu);    // FHC 1Rmu
      if(use_post_pred_fhc_elike)  
	toy_chisq += calc_chisq(h_data_1re,    hErecoExp_1Re);     // FHC 1Re
      if(use_post_pred_fhc_elike)  
	toy_chisq += calc_chisq(h_data_multire,    hErecoExp_MultiRe);     // FHC 1Re
      if(use_post_pred_rhc_mulike) 
	toy_chisq += calc_chisq(h_data_rhc1rmu,hErecoExp_RHC1Rmu); // RHC 1Rmu
      if(use_post_pred_rhc_elike) 
	toy_chisq += calc_chisq(h_data_rhc1re, hErecoExp_RHC1Re);  // RHC 1Re
      if(use_post_pred_fhc_cc1pilike)
        toy_chisq += calc_chisq(h_data_multire, hErecoExp_MultiRe);  // RHC 1Re
      toy_weight = TMath::Exp(-1*toy_chisq/2.);

    }
    vec_like_penalty.push_back(toy_weight);

    hdeltachisq_posterior->Fill(deltachisq0BF_used,toy_weight);
    hdeltachisq_Vs_chisqweight_posterior->Fill(deltachisq0BF_used,toy_chisq);
    hWeight_posterior->Fill(toy_weight);
    hbestfit_posterior_1D->Fill(bestfit0_glob,toy_weight);
    hbestfit_Vs_chisqweight_posterior->Fill(bestfit0_glob,toy_chisq);

    if(*MH_glob=="normal"){
      hMHbf_posterior_1D->Fill(0.1,toy_weight);
      //hMHbf_Vs_chisqweight_posterior->Fill(0.1,toy_chisq);
    }
    else if(*MH_glob=="inverted"){
      hMHbf_posterior_1D->Fill(1.1,toy_weight);
      //hMHbf_Vs_chisqweight_posterior->Fill(1.1,toy_chisq);
    }
    
    //cout << "toy_chisq = " << toy_chisq << " - toy_weight = " << toy_weight << endl;
 

    // Full chisq distribution 
    hchisqminFullChisq->Fill(chisqminFullChisq_glob);
    hbestfitFullChisq_2D->Fill(bestfitFullChisq0_glob,bestfitFullChisq1_glob);
    hbestfitFullChisq_1D->Fill(bestfitFullChisq0_glob);

    // Check that full chisq min is equal to minimized best-fit
    if(chisqminFullChisq_glob!=conv::kALargeNum){
      if( fabs(chisqmin_glob-chisqminFullChisq_glob)>1.E-6 ){
	double diff = chisqmin_glob-chisqminFullChisq_glob;
	hchisqminDiff->Fill(diff);
	cout << std::setprecision(8);
	cout << endl;
	cout << "Toy " << ientry << " --> ToyID = " << ToyID << endl;
	cout << "chisqmin_glob = " << chisqmin_glob
	     << " - bestfit0_glob = " << bestfit0_glob
	     << " - bestfit1_glob = " << bestfit1_glob << endl;
	cout << "chisqminFullChisq_glob = " << chisqminFullChisq_glob 
	     << " - bestfitFullChisq0_glob = " << bestfitFullChisq0_glob 
	     << " - bestfitFullChisq1_glob = " << bestfitFullChisq1_glob << endl;
	cout << "Difference = " << diff << endl;
	cout << endl; 
	NdiffChisq++;
      }
    }
    //   // Check that full chisq min is equal to minimized best-fit
    //   if(bestfitFullChisq0_glob!=conv::kBadResult){
    //     if( fabs(bestfit0_glob-bestfitFullChisq0_glob)>1.E-6 ){
    // 	double diff = bestfit0_glob-bestfitFullChisq0_glob;
    // 	//hchisqminDiff->Fill(diff);
    // 	cout << std::setprecision(8);
    // 	cout << endl;
    // 	cout << "Toy " << ientry << " --> ToyID = " << ToyID << endl;
    // 	cout << "chisqmin_glob = " << chisqmin_glob
    // 	     << " - bestfit0_glob = " << bestfit0_glob
    // 	     << " - bestfit1_glob = " << bestfit1_glob << endl;
    // 	cout << "chisqminFullChisq_glob = " << chisqminFullChisq_glob 
    // 	     << " - bestfitFullChisq0_glob = " << bestfitFullChisq0_glob 
    // 	     << " - bestfitFullChisq1_glob = " << bestfitFullChisq1_glob << endl;
    // 	cout << "Difference = " << diff << endl;
    // 	cout << endl; 
    //     }
    //   }
    
    
  } // end loop over tree entries


  // Calculate coverage
  double FracIn_1sigma = (double) Nin_1sigma / (double) Nentries;
  double FracIn_90CL   = (double) Nin_90CL   / (double) Nentries;
  double FracIn_95CL   = (double) Nin_95CL   / (double) Nentries;
  double FracIn_2sigma = (double) Nin_2sigma / (double) Nentries;  
  double FracIn_3sigma = (double) Nin_3sigma / (double) Nentries;  

  double FracIn_1sigma_nh = (double) Nin_1sigma_nh / (double) Nentries;
  double FracIn_90CL_nh   = (double) Nin_90CL_nh   / (double) Nentries;
  double FracIn_95CL_nh   = (double) Nin_95CL_nh   / (double) Nentries;
  double FracIn_2sigma_nh = (double) Nin_2sigma_nh / (double) Nentries;  
  double FracIn_3sigma_nh = (double) Nin_3sigma_nh / (double) Nentries;  

  double FracIn_1sigma_ih = (double) Nin_1sigma_ih / (double) Nentries;
  double FracIn_90CL_ih   = (double) Nin_90CL_ih   / (double) Nentries;
  double FracIn_95CL_ih   = (double) Nin_95CL_ih   / (double) Nentries;
  double FracIn_2sigma_ih = (double) Nin_2sigma_ih / (double) Nentries;  
  double FracIn_3sigma_ih = (double) Nin_3sigma_ih / (double) Nentries;  
  
  // Calculate p-value and quantile (p-value = gauss area of both tails --> divide by 2)
  double pvalue = GetPvalue_rightside(myvector,mydeltachisq);  
  TF1 f2("f2","ROOT::Math::gaussian_quantile(1.-(x/2.),1.) ",0.,1.);
  double zscore = f2.Eval(pvalue); 

  // Calculate critical value
  
  double FCCritVal_1sigma = CalcDLcrit(myvector,CLfrac_1sigma*100.,0.); 
  double FCCritVal_90CL   = CalcDLcrit(myvector,CLfrac_90CL*100.,0.);
  double FCCritVal_95CL   = CalcDLcrit(myvector,CLfrac_95CL*100.,0.); 
  double FCCritVal_2sigma = CalcDLcrit(myvector,CLfrac_2sigma*100.,0.); 
  double FCCritVal_3sigma = CalcDLcrit(myvector,CLfrac_3sigma*100.,0.); 

  double FCCritVal_1sigma_nh = conv::kBadResult;
  double FCCritVal_90CL_nh = conv::kBadResult;
  double FCCritVal_95CL_nh = conv::kBadResult;
  double FCCritVal_2sigma_nh = conv::kBadResult;
  double FCCritVal_3sigma_nh = conv::kBadResult;
  double FCCritVal_1sigma_ih = conv::kBadResult;
  double FCCritVal_90CL_ih = conv::kBadResult;
  double FCCritVal_95CL_ih = conv::kBadResult;
  double FCCritVal_2sigma_ih = conv::kBadResult;
  double FCCritVal_3sigma_ih = conv::kBadResult;

  if(*trueMH=="normal"){
    FCCritVal_1sigma_nh = CalcDLcrit(myvector_nh,CLfrac_1sigma*100.,0.); 
    FCCritVal_90CL_nh   = CalcDLcrit(myvector_nh,CLfrac_90CL*100.,0.); 
    FCCritVal_95CL_nh   = CalcDLcrit(myvector_nh,CLfrac_95CL*100.,0.); 
    FCCritVal_2sigma_nh = CalcDLcrit(myvector_nh,CLfrac_2sigma*100.,0.); 
    FCCritVal_3sigma_nh = CalcDLcrit(myvector_nh,CLfrac_3sigma*100.,0.); 
  }
  else if(*trueMH=="inverted"){
    FCCritVal_1sigma_ih = CalcDLcrit(myvector_ih,CLfrac_1sigma*100.,0.); 
    FCCritVal_90CL_ih   = CalcDLcrit(myvector_ih,CLfrac_90CL*100.,0.); 
    FCCritVal_95CL_ih   = CalcDLcrit(myvector_ih,CLfrac_95CL*100.,0.); 
    FCCritVal_2sigma_ih = CalcDLcrit(myvector_ih,CLfrac_2sigma*100.,0.); 
    FCCritVal_3sigma_ih = CalcDLcrit(myvector_ih,CLfrac_3sigma*100.,0.); 
  }

  
  // Compute p-value and quantile (use likelihood weights from data sets)
  double pvalue_posterior = GetPvalue_posterior_rightside(myvector,vec_like_penalty,mydeltachisq);  
  double pvalue_posterior_histo = GetPvalue_Histo_rightside(hdeltachisq_posterior,mydeltachisq);  
  double zscore_posterior = f2.Eval(pvalue_posterior); 
 
  double FCCritVal_posterior_1sigma = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_1sigma*100.,0.);
  double FCCritVal_posterior_90CL   = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_90CL*100.,0.);
  double FCCritVal_posterior_95CL   = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_95CL*100.,0.);
  double FCCritVal_posterior_2sigma = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_2sigma*100.,0.);
  double FCCritVal_posterior_3sigma = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_3sigma*100.,0.);
  
  double FCCritVal_posterior_1sigma_histo = CalcDLcrit(hdeltachisq_posterior,CLfrac_1sigma);
  double FCCritVal_posterior_90CL_histo   = CalcDLcrit(hdeltachisq_posterior,CLfrac_90CL);
  double FCCritVal_posterior_95CL_histo   = CalcDLcrit(hdeltachisq_posterior,CLfrac_95CL);
  double FCCritVal_posterior_2sigma_histo = CalcDLcrit(hdeltachisq_posterior,CLfrac_2sigma);
  double FCCritVal_posterior_3sigma_histo = CalcDLcrit(hdeltachisq_posterior,CLfrac_3sigma);

  // Compute errors on the posterior predictive p-value and critical values

  double err_1sigma = get_pvalue_error_vector(myvector,vec_like_penalty,FCCritVal_posterior_1sigma,1.-CLfrac_1sigma); // pass the "p-value"
  double FCCritVal_posterior_1sigma_perr = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_1sigma*100.,+err_1sigma);
  FCCritVal_posterior_1sigma_perr = fabs(FCCritVal_posterior_1sigma - FCCritVal_posterior_1sigma_perr);
  double FCCritVal_posterior_1sigma_merr = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_1sigma*100.,-err_1sigma);
  FCCritVal_posterior_1sigma_merr = fabs(FCCritVal_posterior_1sigma - FCCritVal_posterior_1sigma_merr);
    
  double err_90CL = get_pvalue_error_vector(myvector,vec_like_penalty,FCCritVal_posterior_90CL,1.-CLfrac_90CL); // pass the "p-value"
  double FCCritVal_posterior_90CL_perr = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_90CL*100.,+err_90CL);
  FCCritVal_posterior_90CL_perr = fabs(FCCritVal_posterior_90CL - FCCritVal_posterior_90CL_perr);
  double FCCritVal_posterior_90CL_merr = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_90CL*100.,-err_90CL);
  FCCritVal_posterior_90CL_merr = fabs(FCCritVal_posterior_90CL - FCCritVal_posterior_90CL_merr);

  double err_2sigma = get_pvalue_error_vector(myvector,vec_like_penalty,FCCritVal_posterior_2sigma,1.-CLfrac_2sigma); // pass the "p-value"
  double FCCritVal_posterior_2sigma_perr = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_2sigma*100.,+err_2sigma);
  FCCritVal_posterior_2sigma_perr = fabs(FCCritVal_posterior_2sigma - FCCritVal_posterior_2sigma_perr);
  double FCCritVal_posterior_2sigma_merr = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_2sigma*100.,-err_2sigma);
  FCCritVal_posterior_2sigma_merr = fabs(FCCritVal_posterior_2sigma - FCCritVal_posterior_2sigma_merr);

  double err_3sigma = get_pvalue_error_vector(myvector,vec_like_penalty,FCCritVal_posterior_3sigma,1.-CLfrac_3sigma); // pass the "p-value"
  double FCCritVal_posterior_3sigma_perr = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_3sigma*100.,+err_3sigma);
  FCCritVal_posterior_3sigma_perr = fabs(FCCritVal_posterior_3sigma - FCCritVal_posterior_3sigma_perr);
  double FCCritVal_posterior_3sigma_merr = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_3sigma*100.,-err_3sigma);
  FCCritVal_posterior_3sigma_merr = fabs(FCCritVal_posterior_3sigma - FCCritVal_posterior_3sigma_merr);

  double err_pvalue = get_pvalue_error_vector(myvector,vec_like_penalty,mydeltachisq,pvalue_posterior); 
  double zscore_posterior_perr = f2.Eval(pvalue_posterior-err_pvalue); // smaller p-value gives larger significance 
  zscore_posterior_perr = fabs(zscore_posterior-zscore_posterior_perr); 
  double zscore_posterior_merr = f2.Eval(pvalue_posterior+err_pvalue); // larger p-value gives lower significance 
  zscore_posterior_merr = fabs(zscore_posterior-zscore_posterior_merr); 


  // Draw histograms

  TCanvas *ctrue = new TCanvas("ctrue","ctrue");
  if(IsFit2D){
    htrue_2D->GetXaxis()->SetTitle(namepar0.c_str());
    htrue_2D->GetYaxis()->SetTitle(namepar1.c_str());
    htrue_2D->DrawClone("colz");
  }
  else{
    htrue_1D->GetYaxis()->SetTitle("# entries");
    htrue_1D->GetXaxis()->SetTitle(namepar0.c_str());
    htrue_1D->DrawClone();
  }

  TCanvas *ctrue_MH = new TCanvas("ctrue_MH","ctrue_MH");
  htrue_MH->GetXaxis()->SetBinLabel(1,"normal");
  htrue_MH->GetXaxis()->SetLabelSize(0.09);
  htrue_MH->GetXaxis()->SetBinLabel(2,"inverted");
  htrue_MH->GetXaxis()->SetTitle("True Mass Hierarchy");
  htrue_MH->GetYaxis()->SetTitle("# entries");
  htrue_MH->DrawClone();


  TCanvas *cNevents_fhc_mulike = new TCanvas("cNevents_fhc_mulike","cNevents_fhc_mulike");
  hNevents_fhc_mulike->GetXaxis()->SetTitle("# of events");
  hNevents_fhc_mulike->GetYaxis()->SetTitle("# of entries");
  hNevents_fhc_mulike->DrawClone();

  TCanvas *cNevents_fhc_elike = new TCanvas("cNevents_fhc_elike","cNevents_fhc_elike");
  hNevents_fhc_elike->GetXaxis()->SetTitle("# of events");
  hNevents_fhc_elike->GetYaxis()->SetTitle("# of entries");
  hNevents_fhc_elike->DrawClone();

  TCanvas *cNevents_rhc_mulike = new TCanvas("cNevents_rhc_mulike","cNevents_rhc_mulike");
  hNevents_rhc_mulike->GetXaxis()->SetTitle("# of events");
  hNevents_rhc_mulike->GetYaxis()->SetTitle("# of entries");
  hNevents_rhc_mulike->DrawClone();

  TCanvas *cNevents_rhc_elike = new TCanvas("cNevents_rhc_elike","cNevents_rhc_elike");
  hNevents_rhc_elike->GetXaxis()->SetTitle("# of events");
  hNevents_rhc_elike->GetYaxis()->SetTitle("# of entries");
  hNevents_rhc_elike->DrawClone();

  TCanvas *cNevents_fhc_cc1pilike = new TCanvas("cNevents_fhc_cc1pilike","cNevents_fhc_cc1pilike");
  hNevents_fhc_cc1pilike->GetXaxis()->SetTitle("# of events");
  hNevents_fhc_cc1pilike->GetYaxis()->SetTitle("# of entries");
  hNevents_fhc_cc1pilike->DrawClone();


  TCanvas *cNevents_exp_fhc_mulike = new TCanvas("cNevents_exp_fhc_mulike","cNevents_exp_fhc_mulike");
  hNevents_exp_fhc_mulike->GetXaxis()->SetTitle("# of events");
  hNevents_exp_fhc_mulike->GetYaxis()->SetTitle("# of entries");
  hNevents_exp_fhc_mulike->DrawClone();

  TCanvas *cNevents_exp_fhc_elike = new TCanvas("cNevents_exp_fhc_elike","cNevents_exp_fhc_elike");
  hNevents_exp_fhc_elike->GetXaxis()->SetTitle("# of events");
  hNevents_exp_fhc_elike->GetYaxis()->SetTitle("# of entries");
  hNevents_exp_fhc_elike->DrawClone();

  TCanvas *cNevents_exp_rhc_mulike = new TCanvas("cNevents_exp_rhc_mulike","cNevents_exp_rhc_mulike");
  hNevents_exp_rhc_mulike->GetXaxis()->SetTitle("# of events");
  hNevents_exp_rhc_mulike->GetYaxis()->SetTitle("# of entries");
  hNevents_exp_rhc_mulike->DrawClone();

  TCanvas *cNevents_exp_rhc_elike = new TCanvas("cNevents_exp_rhc_elike","cNevents_exp_rhc_elike");
  hNevents_exp_rhc_elike->GetXaxis()->SetTitle("# of events");
  hNevents_exp_rhc_elike->GetYaxis()->SetTitle("# of entries");
  hNevents_exp_rhc_elike->DrawClone();

  TCanvas *cNevents_exp_fhc_cc1pilike = new TCanvas("cNevents_exp_fhc_cc1pilike","cNevents_exp_fhc_cc1pilike");
  hNevents_exp_fhc_cc1pilike->GetXaxis()->SetTitle("# of events");
  hNevents_exp_fhc_cc1pilike->GetYaxis()->SetTitle("# of entries");
  hNevents_exp_fhc_cc1pilike->DrawClone();


  TCanvas *cbestfit_nh = new TCanvas("cbestfit_nh","cbestfit_nh");
  if(IsFit2D){
    hbestfit_2D_nh->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_2D_nh->GetYaxis()->SetTitle(namepar1.c_str());
    hbestfit_2D_nh->DrawClone("colz");
  }
  else{
    hbestfit_1D_nh->GetYaxis()->SetTitle("# entries");
    hbestfit_1D_nh->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_1D_nh->DrawClone();
  }

  TCanvas *cbestfit_ih = new TCanvas("hbestfit_ih","hbestfit_ih");
  if(IsFit2D){
    hbestfit_2D_ih->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_2D_ih->GetYaxis()->SetTitle(namepar1.c_str());
    hbestfit_2D_ih->DrawClone("colz");
  }
  else{
    hbestfit_1D_ih->GetYaxis()->SetTitle("# entries");
    hbestfit_1D_ih->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_1D_ih->DrawClone();
  }

  TCanvas *cbestfit = new TCanvas("hbestfit","hbestfit");
  if(IsFit2D){
    hbestfit_2D->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_2D->GetYaxis()->SetTitle(namepar1.c_str());
    hbestfit_2D->DrawClone("colz");
  }
  else{
    hbestfit_1D->GetYaxis()->SetTitle("# entries");
    hbestfit_1D->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_1D->DrawClone();
  }

  TCanvas *cbestfitFullChisq = new TCanvas("hbestfitFullChisq","hbestfitFullChisq");
  if(IsFit2D){
    hbestfitFullChisq_2D->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfitFullChisq_2D->GetYaxis()->SetTitle(namepar1.c_str());
    hbestfitFullChisq_2D->DrawClone("colz");
  }
  else{
    hbestfitFullChisq_1D->GetYaxis()->SetTitle("# entries");
    hbestfitFullChisq_1D->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfitFullChisq_1D->DrawClone();
  }

  TCanvas *cbestfit_MH = new TCanvas("hbestfit_MH","hbestfit_MH");
  hbestfit_MH->GetXaxis()->SetBinLabel(1,"normal");
  hbestfit_MH->GetXaxis()->SetLabelSize(0.09);
  hbestfit_MH->GetXaxis()->SetBinLabel(2,"inverted");
  hbestfit_MH->GetXaxis()->SetTitle("Mass Hierarchy");
  hbestfit_MH->GetYaxis()->SetTitle("# entries");
  hbestfit_MH->DrawClone();

  TCanvas *cchisq = new TCanvas("cchisq","cchisq");
  hchisqmin->GetYaxis()->SetTitle("# entries");
  hchisqmin->GetXaxis()->SetTitle("#chi^{2}");
  hchisqmin->DrawClone();

  TCanvas *cchisqminDiff = new TCanvas("hchisqminDiff","hchisqminDiff");
  hchisqminDiff->GetYaxis()->SetTitle("# entries");
  hchisqminDiff->GetXaxis()->SetTitle("#chi^{2}_{Descent} - #chi^{2}_{Global}");
  hchisqminDiff->DrawClone();


  TCanvas *cdeltachisq = new TCanvas("cdeltachisq","cdeltachisq");
  hdeltachisq->GetYaxis()->SetTitle("# entries");
  hdeltachisq->GetXaxis()->SetTitle("#Delta #chi^{2}");
  hdeltachisq->DrawClone();
  TLine *lcritval_90CL = new TLine(FCCritVal_90CL,0,FCCritVal_90CL,hdeltachisq->GetMaximum());
  TLine *lcritval_2sigma = new TLine(FCCritVal_2sigma,0,FCCritVal_2sigma,hdeltachisq->GetMaximum());
  TLine *lcritval_3sigma = new TLine(FCCritVal_3sigma,0,FCCritVal_3sigma,hdeltachisq->GetMaximum());
  lcritval_90CL->SetLineColor(kRed);
  lcritval_2sigma->SetLineColor(kBlue);
  lcritval_3sigma->SetLineColor(kGreen);
  lcritval_90CL->SetLineWidth(2);
  lcritval_2sigma->SetLineWidth(2);
  lcritval_3sigma->SetLineWidth(2);
  lcritval_90CL->DrawClone();
  lcritval_2sigma->DrawClone();
  //lcritval_3sigma->DrawClone();
  TLine *ldeltachisq_data = new TLine(mydeltachisq,0,mydeltachisq,hdeltachisq->GetMaximum());
  ldeltachisq_data->SetLineColor(kGreen);
  ldeltachisq_data->SetLineWidth(2);
  if(mydeltachisq>0.) ldeltachisq_data->DrawClone();
  TLegend *ldeltachisq = new TLegend(0.4,0.6,0.6,0.9);
  ldeltachisq->SetFillColor(0);
  ldeltachisq->SetLineColor(0);
  ldeltachisq->SetTextFont(132);
  ldeltachisq->SetTextSize(0.07);
  ldeltachisq->AddEntry(lcritval_90CL,TString::Format("#Delta #chi^{2}_{crit}(90%%CL) = %1.2f",FCCritVal_90CL).Data(),"l");
  ldeltachisq->AddEntry(lcritval_2sigma,TString::Format("#Delta #chi^{2}_{crit} (2#sigma) = %1.2f",FCCritVal_2sigma).Data(),"l");
  //ldeltachisq->AddEntry(lcritval_3sigma,TString::Format("#Delta #chi^{2}_{crit} (3#sigma) = %1.2f",FCCritVal_3sigma).Data(),"l");
  if(mydeltachisq>0.) ldeltachisq->AddEntry(ldeltachisq_data,TString::Format("#Delta #chi^{2}_{data} = %1.3f",mydeltachisq).Data(),"l");
  ldeltachisq->DrawClone();
  TLatex *latPvalue = new TLatex(hdeltachisq->GetXaxis()->GetBinLowEdge(50),hdeltachisq->GetMaximum()*0.3, TString::Format("P-value = %1.3f",pvalue).Data() );
  if(mydeltachisq>0.) latPvalue->Draw("same");
  TLatex *latZscore = new TLatex(hdeltachisq->GetXaxis()->GetBinLowEdge(50),hdeltachisq->GetMaximum()*0.1, TString::Format("Z = %1.2f#sigma",zscore).Data() );
  if(mydeltachisq>0.) latZscore->Draw("same");


  TCanvas *cdeltachisq_nh = new TCanvas("cdeltachisq_nh","cdeltachisq_nh");
  hdeltachisq_nh->GetYaxis()->SetTitle("# entries");
  hdeltachisq_nh->GetXaxis()->SetTitle("#Delta #chi^{2}");
  hdeltachisq_nh->DrawClone();
  TLine *lcritval_90CL_nh;
  TLine *lcritval_2sigma_nh;
  TLine *lcritval_3sigma_nh;
  if(*trueMH=="normal"){
    lcritval_90CL_nh = new TLine(FCCritVal_90CL_nh,0,FCCritVal_90CL_nh,hdeltachisq_nh->GetMaximum());
    lcritval_2sigma_nh = new TLine(FCCritVal_2sigma_nh,0,FCCritVal_2sigma_nh,hdeltachisq_nh->GetMaximum());
    lcritval_3sigma_nh = new TLine(FCCritVal_3sigma_nh,0,FCCritVal_3sigma_nh,hdeltachisq_nh->GetMaximum());
    lcritval_90CL_nh->SetLineColor(kRed);
    lcritval_2sigma_nh->SetLineColor(kBlue);
    lcritval_3sigma_nh->SetLineColor(kGreen);
    lcritval_90CL_nh->DrawClone();
    lcritval_2sigma_nh->DrawClone();
    //lcritval_3sigma_nh->DrawClone();
  }
  if(mydeltachisq>0.) ldeltachisq_data->DrawClone();
  TLegend *ldeltachisq_nh = new TLegend(0.6,0.45,0.8,0.75);
  ldeltachisq_nh->SetFillColor(0);
  ldeltachisq_nh->SetLineColor(0);
  ldeltachisq_nh->SetTextFont(132);
  ldeltachisq_nh->SetTextSize(0.07);
  ldeltachisq_nh->AddEntry(lcritval_90CL_nh,TString::Format("#Delta #chi^{2}_{crit}(90%%CL) = %1.2f",FCCritVal_90CL_nh).Data(),"l");
  ldeltachisq_nh->AddEntry(lcritval_2sigma_nh,TString::Format("#Delta #chi^{2}_{crit} (2#sigma) = %1.2f",FCCritVal_2sigma_nh).Data(),"l");
  //ldeltachisq_nh->AddEntry(lcritval_3sigma_nh,TString::Format("#Delta #chi^{2}_{crit} (3#sigma) = %1.2f",FCCritVal_3sigma_nh).Data(),"l");
  if(mydeltachisq>0.) ldeltachisq_nh->AddEntry(ldeltachisq_data,TString::Format("#Delta #chi^{2}_{data} = %1.3f",mydeltachisq).Data(),"l");
  ldeltachisq_nh->DrawClone();


  TCanvas *cdeltachisq_ih = new TCanvas("cdeltachisq_ih","cdeltachisq_ih");
  hdeltachisq_ih->GetYaxis()->SetTitle("# entries");
  hdeltachisq_ih->GetXaxis()->SetTitle("#Delta #chi^{2}");
  hdeltachisq_ih->DrawClone();
  TLine *lcritval_90CL_ih;
  TLine *lcritval_2sigma_ih;
  TLine *lcritval_3sigma_ih;
  if(*trueMH=="inverted"){
    lcritval_90CL_ih = new TLine(FCCritVal_90CL_ih,0,FCCritVal_90CL_ih,hdeltachisq_ih->GetMaximum());
    lcritval_2sigma_ih = new TLine(FCCritVal_2sigma_ih,0,FCCritVal_2sigma_ih,hdeltachisq_ih->GetMaximum());
    lcritval_3sigma_ih = new TLine(FCCritVal_3sigma_ih,0,FCCritVal_3sigma_ih,hdeltachisq_ih->GetMaximum());
    lcritval_90CL_ih->SetLineColor(kRed);
    lcritval_2sigma_ih->SetLineColor(kBlue);
    lcritval_3sigma_ih->SetLineColor(kGreen);
    lcritval_90CL_ih->DrawClone();
    //lcritval_2sigma_ih->DrawClone();
    lcritval_3sigma_ih->DrawClone();
  }
  if(mydeltachisq>0.) ldeltachisq_data->DrawClone();
  TLegend *ldeltachisq_ih = new TLegend(0.6,0.45,0.8,0.75);
  ldeltachisq_ih->SetFillColor(0);
  ldeltachisq_ih->SetLineColor(0);
  ldeltachisq_ih->SetTextFont(132);
  ldeltachisq_ih->SetTextSize(0.07);
  ldeltachisq_ih->AddEntry(lcritval_90CL_ih,TString::Format("#Delta #chi^{2}_{crit}(90%%CL) = %1.2f",FCCritVal_90CL_ih).Data(),"l");
  ldeltachisq_ih->AddEntry(lcritval_2sigma_ih,TString::Format("#Delta #chi^{2}_{crit} (2#sigma) = %1.2f",FCCritVal_2sigma_ih).Data(),"l");
  //ldeltachisq_ih->AddEntry(lcritval_3sigma_ih,TString::Format("#Delta #chi^{2}_{crit} (3#sigma) = %1.2f",FCCritVal_3sigma_ih).Data(),"l");
  if(mydeltachisq>0.) ldeltachisq_ih->AddEntry(ldeltachisq_data,TString::Format("#Delta #chi^{2}_{data} = %1.3f",mydeltachisq).Data(),"l");
  ldeltachisq_ih->DrawClone();



  // Posterior

  TCanvas *cdeltachisq_posterior = new TCanvas("cdeltachisq_posterior","cdeltachisq_posterior");
  hdeltachisq_posterior->Scale(1. / hdeltachisq_posterior->Integral());
  double mergebin = newbinwidth; 
  hdeltachisq_posterior->Rebin(mergebin);
  hdeltachisq_posterior->GetXaxis()->SetTitle("#Delta #chi^{2} posterior");
  hdeltachisq_posterior->GetYaxis()->SetTitle("# of entries");
  hdeltachisq_posterior->DrawClone();
  TLine *lcritval_1sigma_post = new TLine(FCCritVal_posterior_1sigma,0,FCCritVal_posterior_1sigma,hdeltachisq_posterior->GetMaximum());
  TLine *lcritval_90CL_post = new TLine(FCCritVal_posterior_90CL,0,FCCritVal_posterior_90CL,hdeltachisq_posterior->GetMaximum());
  TLine *lcritval_2sigma_post = new TLine(FCCritVal_posterior_2sigma,0,FCCritVal_posterior_2sigma,hdeltachisq_posterior->GetMaximum());
  TLine *lcritval_3sigma_post = new TLine(FCCritVal_posterior_3sigma,0,FCCritVal_posterior_3sigma,hdeltachisq_posterior->GetMaximum());
  lcritval_90CL_post->SetLineColor(kRed);
  lcritval_2sigma_post->SetLineColor(kBlue);
  lcritval_3sigma_post->SetLineColor(kGreen);
  lcritval_90CL_post->SetLineWidth(2);
  lcritval_2sigma_post->SetLineWidth(2);
  lcritval_3sigma_post->SetLineWidth(2);
  lcritval_90CL_post->DrawClone();
  lcritval_2sigma_post->DrawClone();
  //lcritval_3sigma_post->DrawClone();

  ldeltachisq_data->SetLineColor(kGreen);
  ldeltachisq_data->SetLineWidth(2);
  ldeltachisq_data->DrawClone();

  TLegend *ldeltachisq_post = new TLegend(0.4,0.6,0.6,0.9);
  ldeltachisq_post->SetFillColor(0);
  ldeltachisq_post->SetLineColor(0);
  ldeltachisq_post->SetTextFont(132);
  ldeltachisq_post->SetTextSize(0.07);
  ldeltachisq_post->AddEntry(lcritval_90CL_post,TString::Format("#Delta #chi^{2}_{crit}(90%%CL) = %1.2f",FCCritVal_posterior_90CL).Data(),"l");
  ldeltachisq_post->AddEntry(lcritval_2sigma_post,TString::Format("#Delta #chi^{2}_{crit} (2#sigma) = %1.2f",FCCritVal_posterior_2sigma).Data(),"l");
  //ldeltachisq_post->AddEntry(lcritval_3sigma_post,TString::Format("#Delta #chi^{2}_{crit} (3#sigma) = %1.2f",FCCritVal_posterior_3sigma).Data(),"l");
  if(mydeltachisq>0.) ldeltachisq_post->AddEntry(ldeltachisq_data,TString::Format("#Delta #chi^{2}_{data} = %1.3f",mydeltachisq).Data(),"l");
  ldeltachisq_post->DrawClone();
  TLatex *latPvalue_post = new TLatex(hdeltachisq_posterior->GetXaxis()->GetBinLowEdge(120),hdeltachisq_posterior->GetMaximum()*0.3, TString::Format("P-value = %1.3f",pvalue_posterior).Data() );
  if(mydeltachisq>0.) latPvalue_post->Draw("same");
  TLatex *latZscore_post = new TLatex(hdeltachisq_posterior->GetXaxis()->GetBinLowEdge(120),hdeltachisq_posterior->GetMaximum()*0.1, TString::Format("Z = %1.2f#sigma",zscore_posterior).Data() );
  if(mydeltachisq>0.) latZscore_post->Draw("same");


  TCanvas *cWeight_posterior = new TCanvas("cWeight_posterior","cWeight_posterior");
  hWeight_posterior->GetXaxis()->SetTitle("L (data,fake exp)");
  hWeight_posterior->GetYaxis()->SetTitle("# of entries");
  hWeight_posterior->DrawClone();

  TCanvas *cdeltachisq_Vs_chisqweight_posterior = new TCanvas("cdeltachisq_Vs_chisqweight_posterior","cdeltachisq_Vs_chisqweight_posterior");
  hdeltachisq_Vs_chisqweight_posterior->GetXaxis()->SetTitle("test statistic");
  hdeltachisq_Vs_chisqweight_posterior->GetYaxis()->SetTitle("#chi^{2} (data,fake exp)");
  hdeltachisq_Vs_chisqweight_posterior->Draw("colz");

  TCanvas *cbestfit_posterior_1D = new TCanvas("cbestfit_posterior_1D","cbestfit_posterior_1D");
  hbestfit_posterior_1D->GetXaxis()->SetTitle("#delta_{CP} best-fit posterior");
  hbestfit_posterior_1D->GetYaxis()->SetTitle("# of entries");
  hbestfit_posterior_1D->DrawClone();

  TCanvas *cbestfit_Vs_chisqweight_posterior = new TCanvas("cbestfit_Vs_chisqweight_posterior","cbestfit_Vs_chisqweight_posterior");
  hbestfit_Vs_chisqweight_posterior->GetXaxis()->SetTitle("#delta_{CP} best-fit posterior");
  hbestfit_Vs_chisqweight_posterior->GetYaxis()->SetTitle("#chi^{2} (data,fake exp)");
  hbestfit_Vs_chisqweight_posterior->Draw("colz");

  TCanvas *cMHbf_posterior_1D = new TCanvas("cMHbf_posterior_1D","cMHbf_posterior_1D");
  hMHbf_posterior_1D->GetXaxis()->SetTitle("Mass Hierarchy best-fit posterior");
  hMHbf_posterior_1D->GetYaxis()->SetTitle("# of entries");
  hMHbf_posterior_1D->DrawClone();

  TCanvas *cMHbf_Vs_chisqweight_posterior = new TCanvas("cMHbf_Vs_chisqweight_posterior","cMHbf_Vs_chisqweight_posterior");
  hMHbf_Vs_chisqweight_posterior->GetXaxis()->SetTitle("Mass Hierarchy best-fit posterior");
  hMHbf_Vs_chisqweight_posterior->GetYaxis()->SetTitle("#chi^{2} (data,fake exp)");
  hMHbf_Vs_chisqweight_posterior->Draw("colz");


  TCanvas *cComp_deltachisq = new TCanvas("cComp_deltachisq","cComp_deltachisq");
  hdeltachisq->Scale(1./hdeltachisq->Integral());
  hdeltachisq->SetLineWidth(2);
  hdeltachisq_posterior->Scale(1./hdeltachisq_posterior->Integral());
  hdeltachisq_posterior->SetLineColor(kRed);
  hdeltachisq_posterior->SetLineWidth(2);
  hdeltachisq_posterior->DrawClone();
  hdeltachisq->DrawClone("same");
  TLegend *ldeltachisq_comp = new TLegend(0.4,0.6,0.6,0.9);
  ldeltachisq_comp->SetFillColor(0);
  ldeltachisq_comp->SetLineColor(0);
  ldeltachisq_comp->SetTextFont(132);
  ldeltachisq_comp->SetTextSize(0.07);
  ldeltachisq_comp->AddEntry(hdeltachisq,"#Delta #chi^{2} prior","l");
  ldeltachisq_comp->AddEntry(hdeltachisq_posterior,"#Delta #chi^{2} posterior","l");
  ldeltachisq_comp->DrawClone();


  // Check a single toy

  if(CheckSingleToy){
    
    treein->GetEntry(ToySel);

    int npt_chisq_x = conv::kBadResult;
    int npt_chisq_y = conv::kBadResult;

    TCanvas *cchisq_out_GetEntry = new TCanvas("cchisq_out_GetEntry","cchisq_out_GetEntry");
    if(IsFit2D){
      hchisq2D_out->GetXaxis()->SetTitle(namepar0.c_str());
      hchisq2D_out->GetYaxis()->SetTitle(namepar1.c_str());
      hchisq2D_out->Draw("colz");
      npt_chisq_x = hchisq2D_out->GetNbinsX();
      npt_chisq_y = hchisq2D_out->GetNbinsY();
    }
    else{
      if(*trueMH=="inverted"){
        hchisq1D_out_ih->GetXaxis()->SetTitle(namepar0.c_str());
        hchisq1D_out_ih->Draw("histo");
        npt_chisq_x = hchisq1D_out_ih->GetNbinsX();
        npt_chisq_y = hchisq1D_out_ih->GetNbinsY();
       }
      if(*trueMH=="normal"){
        hchisq1D_out_nh->GetXaxis()->SetTitle(namepar0.c_str());
        hchisq1D_out_nh->Draw("histo");
        npt_chisq_x = hchisq1D_out_nh->GetNbinsX();
        npt_chisq_y = hchisq1D_out_nh->GetNbinsY();
      }
    }
    TCanvas *cdeltachisq_out_GetEntry = new TCanvas("cdeltachisq_out_GetEntry","cdeltachisq_out_GetEntry");
    if(IsFit2D){
      hdeltachisq2D_out->GetXaxis()->SetTitle(namepar0.c_str());
      hdeltachisq2D_out->GetYaxis()->SetTitle(namepar1.c_str());
      hdeltachisq2D_out->Draw("colz");
    }    
    else{
      if(*trueMH=="inverted"){
        hdeltachisq1D_out_ih->GetXaxis()->SetTitle(namepar0.c_str());
        hdeltachisq1D_out_ih->Draw("histo");
      }
      if(*trueMH=="normal"){
        hdeltachisq1D_out_nh->GetXaxis()->SetTitle(namepar0.c_str());
        hdeltachisq1D_out_nh->Draw("histo");
      }
    }
    cout << endl;
    
    cout << "npt_chisq_x = " << npt_chisq_x << endl;
    cout << "npt_chisq_y = " << npt_chisq_y << endl;
    
    for(int iy=0;iy<npt_chisq_y;iy++){
      for(int ix=0;ix<npt_chisq_x;ix++){
       	if(IsFit2D) cout << hchisq2D_out->GetBinContent(ix+1,iy+1) << " ";

        else
          if(*trueMH=="normal"){
            cout << hchisq1D_out_nh->GetBinContent(ix+1) << " ";
          }
        if(*trueMH=="inverted"){
          cout << hchisq1D_out_ih->GetBinContent(ix+1) << " ";
        }
      }
      cout << endl;
    }
    cout << endl;
    TCanvas *cErecoToy_1Rmu = new TCanvas("cErecoToy_1Rmu","cErecoToy_1Rmu");
    hErecoToy_1Rmu->Draw("colz");
    TCanvas *cErecoToy_1Re = new TCanvas("cErecoToy_1Re","cErecoToy_1Re");
    hErecoToy_1Re->Draw("colz");
    TCanvas *cErecoToy_MultiRe = new TCanvas("cErecoToy_MultiRe","cErecoToy_MultiRe");
    hErecoToy_MultiRe->Draw("colz");
    TCanvas *cErecoToy_RHC1Rmu = new TCanvas("cErecoToy_RHC1Rmu","cErecoToy_RHC1Rmu");
    hErecoToy_RHC1Rmu->Draw("colz");
    TCanvas *cErecoToy_RHC1Re = new TCanvas("cErecoToy_RHC1Re","cErecoToy_RHC1Re");
    hErecoToy_RHC1Re->Draw("colz");
    TCanvas *cErecoToy_MultiRe = new TCanvas("cErecoToy_MultiRe","cErecoToy_MultiRe");
    hErecoToy_MultiRe->Draw("colz");

    cout << endl;
    cout << "Events in fake data spectra:"<< endl;
    cout << " - # of 1Rmu events    = "   << hErecoToy_1Rmu->Integral() << endl;
    cout << " - # of RHC1Rmu events = "   << hErecoToy_RHC1Rmu->Integral() << endl;
    cout << " - # of 1Re events     = "   << hErecoToy_1Re->Integral(0,15,0,25,"") << endl;
    cout << " - # of RHC1Re events  = "   << hErecoToy_RHC1Re->Integral(0,15,0,25,"") << endl;
    cout << " - # of MultiRe events = "   << hErecoToy_MultiRe->Integral(0,15,0,25,"") << endl;
    cout << " - # of MultiRe events = "     << hErecoToy_MultiRe->Integral() << endl;
    cout << endl;
  }

  

  // Print plots
  
  string storedparname0 = *trueparname0;
  string storedparname1 = *trueparname1;
  if(storedparname1=="") storedparname1 = "1D";
  string storedmh       = *trueMH;
  double storetrue0      = truepar0;
  double storetrue1      = truepar1;
  
  //cout << storetrue0 << " " << storetrue1 << endl;

  cNevents_fhc_mulike->Print( TString::Format("%s_Nevents_fhc_mulike_%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  
  cNevents_fhc_elike->Print( TString::Format("%s_Nevents_fhc_elike_%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  
  cNevents_rhc_mulike->Print( TString::Format("%s_Nevents_rhc_mulike_%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  
  cNevents_rhc_elike->Print( TString::Format("%s_Nevents_rhc_elike_%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cNevents_fhc_cc1pilike->Print( TString::Format("%s_Nevents_fhc_cc1pilike_%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  
  cbestfit_nh->Print( TString::Format("%s_cbestfit_nh_%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  
  cbestfit_ih->Print( TString::Format("%s_cbestfit_ih_%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str())  );

  cbestfit->Print( TString::Format("%s_cbestfit_%s_%1.3f_%s_%1.3f_%s.pdf",
				   mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cbestfit_MH->Print( TString::Format("%s_cbestfit_MH_%s_%1.3f_%s_%1.3f_%s.pdf",
				      mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cchisq->Print( TString::Format("%s_chisq_glob_%s_%1.3f_%s_%1.3f_%s.pdf",
				 mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cdeltachisq->Print( TString::Format("%s_deltachisq_glob_%s_%1.3f_%s_%1.3f_%s.pdf",
				      mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  
  cdeltachisq_nh->Print( TString::Format("%s_deltachisq_nh_%s_%1.3f_%s_%1.3f_%s.pdf",
					 mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cdeltachisq_ih->Print( TString::Format("%s_deltachisq_ih_%s_%1.3f_%s_%1.3f_%s.pdf",
					 mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  
  cdeltachisq_posterior->Print( TString::Format("%s_deltachisq_posterior_%s_%1.3f_%s_%1.3f_%s.pdf",
						mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cWeight_posterior->Print( TString::Format("%s_Weight_posterior_%s_%1.3f_%s_%1.3f_%s.pdf",
  					    mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  
  cdeltachisq_Vs_chisqweight_posterior->Print( TString::Format("%s_deltachisq_Vs_chisqweight_posterior_%s_%1.3f_%s_%1.3f_%s.pdf",
  							       mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cbestfit_posterior_1D->Print( TString::Format("%s_bestfit_posterior_1D_%s_%1.3f_%s_%1.3f_%s.pdf",
  						mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cbestfit_Vs_chisqweight_posterior->Print( TString::Format("%s_bestfit_Vs_chisqweight_posterior_%s_%1.3f_%s_%1.3f_%s.pdf",
  							    mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cMHbf_Vs_chisqweight_posterior->Print( TString::Format("%s_MHbf_Vs_chisqweight_posterior_%s_%1.3f_%s_%1.3f_%s.pdf",
  							 mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  cMHbf_posterior_1D->Print( TString::Format("%s_MHbf_posterior_1D_%s_%1.3f_%s_%1.3f_%s.pdf",
  							 mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  TFile f_out("fc_analyze.root","recreate");
  cdeltachisq->Write("cdeltachisq");
  hdeltachisq->Write("hdeltachisq");

  cNevents_fhc_mulike->Write("cNevents_fhc_mulike");  
  hNevents_fhc_mulike->Write("hNevents_fhc_mulike");
  
  cNevents_fhc_elike->Write("cNevents_fhc_elike");  
  hNevents_fhc_elike->Write("hNevents_fhc_elike");
  cNevents_rhc_mulike->Write("cNevents_rhc_mulike");  
  hNevents_rhc_mulike->Write("hNevents_rhc_mulike");
  cNevents_rhc_elike->Write("cNevents_rhc_elike");  
  hNevents_rhc_elike->Write("hNevents_rhc_elike");
  cNevents_fhc_cc1pilike->Write("cNevents_fhc_cc1pilike");  
  hNevents_fhc_cc1pilike->Write("hNevents_fhc_cc1pilike");

  cbestfit_nh->Write("cbestfit_nh");  
  cbestfit_ih->Write("cbestfit_ih");  
  cbestfit->Write("cbestfit");  
  
  cbestfit_MH->Write("cbestfit_MH");  
  hbestfit_MH->Write("hbestfit_MH");
  
  cchisq->Write("cchisq");  
  
  cdeltachisq->Write("cdeltachisq");  
  hdeltachisq->Write("hdeltachisq");
  cdeltachisq_nh->Write("cdeltachisq_nh");  
  hdeltachisq_nh->Write("hdeltachisq_nh");
  cdeltachisq_ih->Write("cdeltachisq_ih");  
  hdeltachisq_ih->Write("hdeltachisq_ih"); 
  
  f_out.Close();


  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  // ->Print( TString::Format("%s__%s_%1.3f_%s_%1.3f_%s.pdf",mytag.c_str(),storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );

  
  


  // Print results 
  
  cout << endl;
  cout << endl;
  cout << "/// True parameters ///" << endl;
  cout << " - truepar0: " << *trueparname0 << " = " << truepar0 << endl;
  cout << " - truepar1: " << *trueparname1 << " = " << truepar1 << endl;
  cout << " - Mass Hierarchy: " << *trueMH << endl;
  cout << endl;
  cout << " - Difference between descent gradient and global minimum" << endl;
  cout << " - # of toys with different minim and full chisq minimum: "
       << NdiffChisq << " / " << Nentries << " toys" << endl;
  cout << " - Mean difference = " << hchisqminDiff->GetMean() << endl;
  cout << " - RMS difference = " << hchisqminDiff->GetRMS() << endl;
  cout << endl;
  cout << "/// Stat+Syst event rate ///" << endl;
  cout << " - FHC mu-like:" << endl;
  cout << "      Mean = " << hNevents_fhc_mulike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_fhc_mulike->GetRMS() << endl; 
  if(hNevents_fhc_mulike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_fhc_mulike->GetRMS() / hNevents_fhc_mulike->GetMean() << endl; 
  cout << " - FHC e-like:" << endl;
  cout << "      Mean = " << hNevents_fhc_elike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_fhc_elike->GetRMS() << endl; 
  if(hNevents_fhc_elike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_fhc_elike->GetRMS() / hNevents_fhc_elike->GetMean() << endl; 
  cout << " - FHC cc1pi-like:" << endl;
  cout << "      Mean = " << hNevents_fhc_cc1pilike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_fhc_cc1pilike->GetRMS() << endl; 
  if(hNevents_fhc_cc1pilike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_fhc_cc1pilike->GetRMS() / hNevents_fhc_cc1pilike->GetMean() << endl; 
  cout << " - RHC mu-like:" << endl;
  cout << "      Mean = " << hNevents_rhc_mulike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_rhc_mulike->GetRMS() << endl; 
  if(hNevents_rhc_mulike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_rhc_mulike->GetRMS() / hNevents_rhc_mulike->GetMean() << endl; 
  cout << " - RHC e-like:" << endl;
  cout << "      Mean = " << hNevents_rhc_elike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_rhc_elike->GetRMS() << endl; 
  if(hNevents_rhc_elike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_rhc_elike->GetRMS() / hNevents_rhc_elike->GetMean() << endl;
  cout << " - CC1 pi e-like:" << endl;
  cout << "      Mean = " << hNevents_fhc_cc1pilike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_fhc_cc1pilike->GetRMS() << endl; 
  if(hNevents_fhc_cc1pilike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_fhc_cc1pilike->GetRMS() / hNevents_fhc_cc1pilike->GetMean() << endl; 
  cout << endl;
  
  cout << "/// Syst only event rate ///" << endl;
  cout << " - FHC mu-like:" << endl;
  cout << "      Mean = " << hNevents_exp_fhc_mulike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_exp_fhc_mulike->GetRMS() << endl; 
  if(hNevents_exp_fhc_mulike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_exp_fhc_mulike->GetRMS() / hNevents_exp_fhc_mulike->GetMean() << endl; 
  cout << " - FHC e-like:" << endl;
  cout << "      Mean = " << hNevents_exp_fhc_elike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_exp_fhc_elike->GetRMS() << endl; 
  if(hNevents_exp_fhc_elike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_exp_fhc_elike->GetRMS() / hNevents_exp_fhc_elike->GetMean() << endl; 
  cout << " - FHC cc1pi-like:" << endl;
  cout << "      Mean = " << hNevents_exp_fhc_cc1pilike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_exp_fhc_cc1pilike->GetRMS() << endl; 
  if(hNevents_exp_fhc_cc1pilike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_exp_fhc_cc1pilike->GetRMS() / hNevents_exp_fhc_cc1pilike->GetMean() << endl; 
  cout << " - RHC mu-like:" << endl;
  cout << "      Mean = " << hNevents_exp_rhc_mulike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_exp_rhc_mulike->GetRMS() << endl; 
  if(hNevents_exp_rhc_mulike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_exp_rhc_mulike->GetRMS() / hNevents_exp_rhc_mulike->GetMean() << endl; 
  cout << " - RHC e-like:" << endl;
  cout << "      Mean = " << hNevents_exp_rhc_elike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_exp_rhc_elike->GetRMS() << endl; 
  if(hNevents_exp_rhc_elike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_exp_rhc_elike->GetRMS() / hNevents_exp_rhc_elike->GetMean() << endl;
  cout << " - CC1pi e-like:" << endl;
  cout << "      Mean = " << hNevents_exp_fhc_cc1pilike->GetMean() << endl;
  cout << "      RMS  = " << hNevents_exp_fhc_cc1pilike->GetRMS() << endl; 
  if(hNevents_exp_fhc_cc1pilike->GetMean()>0) cout << "      RMS / Mean = " << hNevents_exp_fhc_cc1pilike->GetRMS() / hNevents_exp_fhc_cc1pilike->GetMean() << endl; 


  cout << "Empty exp events in: " << Nentries << " Toys" << endl;
  cout << " - FHC mu-like"  << empty_nev_exp_1rmu    << endl;
  cout << " - FHC e-like"   << empty_nev_exp_1re     << endl;
  cout << " - RHC mu-like"  << empty_nev_exp_rhc1rmu << endl;
  cout << " - RHC e-like"   << empty_nev_exp_rhc1re  << endl;
  cout << " - CC1pi e-like" << empty_nev_exp_multire << endl;


  cout << endl;

  cout << "/// Global Result ///" << endl;
  cout << " - P-value = " << pvalue << endl;
  cout << " - Z = " << zscore << " sigma" << endl;
  cout << endl;
  cout << " - Coverage:" << endl;
  cout << "     Fraction of toys inside 1sigma = " << FracIn_1sigma << endl;
  cout << "     Fraction of toys inside 90%CL   = " << FracIn_90CL   << endl;
  cout << "     Fraction of toys inside 95%CL   = " << FracIn_95CL   << endl;
  cout << "     Fraction of toys inside 2sigma = " << FracIn_2sigma << endl;
  cout << "     Fraction of toys inside 3sigma = " << FracIn_3sigma << endl;
  cout << endl;
  cout << " - Feldman-Cousins critical values:" << endl;
  cout << "     Crit Value 1sigma = " << FCCritVal_1sigma << endl; 
  cout << "     Crit Value 90%CL   = " << FCCritVal_90CL << endl; 
  cout << "     Crit Value 95%CL   = " << FCCritVal_95CL << endl; 
  cout << "     Crit Value 2sigma = " << FCCritVal_2sigma << endl; 
  cout << "     Crit Value 3sigma = " << FCCritVal_3sigma << endl; 
  cout << endl;

  cout << "/// Posterior predictive ///" << endl;
  cout << "   Use the data samples:" << endl;
  if(use_post_pred_fhc_mulike) cout << "     - FHC mu-like" << endl; 
  if(use_post_pred_fhc_elike) cout  << "     - FHC e-like"  << endl; 
  if(use_post_pred_fhc_mulike) cout << "     - RHC mu-like" << endl; 
  if(use_post_pred_rhc_elike) cout  << "     - RHC e-like"  << endl;
  if(use_post_pred_fhc_cc1pilike) cout  << "   - CC1pi e-like"  << endl; 
  cout << "   p-value = " << pvalue_posterior 
       << " +- " << err_pvalue
       << endl;
  cout << "   p-value (histo) = " << pvalue_posterior_histo << endl;
  cout << "   z-score = " << zscore_posterior
       << " + " << zscore_posterior_perr
       << " - " << zscore_posterior_merr
       << " sigma" << endl;
  //cout << " - Coverage:" << endl;
  //cout << "     Fraction of toys inside 1sigma = " << FracIn_1sigma << endl;
  //cout << "     Fraction of toys inside 90%CL   = " << FracIn_90CL   << endl;
  //cout << "     Fraction of toys inside 95%CL   = " << FracIn_95CL   << endl;
  //cout << "     Fraction of toys inside 2sigma = " << FracIn_2sigma << endl;
  //cout << "     Fraction of toys inside 3sigma = " << FracIn_3sigma << endl;
  cout << endl;
  cout << " - Feldman-Cousins critical values:" << endl;
  cout << "     Crit Value 1sigma = " << FCCritVal_posterior_1sigma 
       << " + " << FCCritVal_posterior_1sigma_perr
       << " - " << FCCritVal_posterior_1sigma_merr
       << endl; 
  cout << "     Crit Value 90%CL   = " << FCCritVal_posterior_90CL 
       << " + " << FCCritVal_posterior_90CL_perr
       << " - " << FCCritVal_posterior_90CL_merr
       << endl; 
  //cout << "     Crit Value 95%CL   = " << FCCritVal_posterior_95CL << endl; 
  cout << "     Crit Value 2sigma = " << FCCritVal_posterior_2sigma
       << " + " << FCCritVal_posterior_2sigma_perr
       << " - " << FCCritVal_posterior_2sigma_merr
       << endl;
  cout << "     Crit Value 3sigma = " << FCCritVal_posterior_3sigma
       << " + " << FCCritVal_posterior_3sigma_perr
       << " - " << FCCritVal_posterior_3sigma_merr
       << endl;
  cout << " - Feldman-Cousins critical values (histo):" << endl;
  cout << "     Crit Value 1sigma = " << FCCritVal_posterior_1sigma_histo << endl; 
  cout << "     Crit Value 90%CL   = " << FCCritVal_posterior_90CL_histo << endl; 
  cout << "     Crit Value 95%CL   = " << FCCritVal_posterior_95CL_histo << endl; 
  cout << "     Crit Value 2sigma = " << FCCritVal_posterior_2sigma_histo << endl; 
  cout << "     Crit Value 3sigma = " << FCCritVal_posterior_3sigma_histo << endl; 
  cout << endl;

  cout << "/// Normal Hierarchy Result ///" << endl;
  cout << " - Coverage:" << endl;
  cout << "     Fraction of toys inside 1sigma = " << FracIn_1sigma_nh << endl;
  cout << "     Fraction of toys inside 90%CL   = " << FracIn_90CL_nh   << endl;
  cout << "     Fraction of toys inside 95%CL   = " << FracIn_95CL_nh   << endl;
  cout << "     Fraction of toys inside 2sigma = " << FracIn_2sigma_nh << endl;
  cout << "     Fraction of toys inside 3sigma = " << FracIn_3sigma_nh << endl;
  cout << endl;
  cout << " - Feldman-Cousins critical values:" << endl;
  cout << "     Crit Value 1sigma = " << FCCritVal_1sigma_nh << endl; 
  cout << "     Crit Value 90%CL = "   << FCCritVal_90CL_nh << endl; 
  cout << "     Crit Value 95%CL = "   << FCCritVal_95CL_nh << endl; 
  cout << "     Crit Value 2sigma = " << FCCritVal_2sigma_nh << endl; 
  cout << "     Crit Value 3sigma = " << FCCritVal_3sigma_nh << endl; 
  cout << endl;
  cout << "/// Inverted Hierarchy Result ///" << endl;
  cout << " - Coverage:" << endl;
  cout << "     Fraction of toys inside 1sigma = " << FracIn_1sigma_ih << endl;
  cout << "     Fraction of toys inside 90%CL   = " << FracIn_90CL_ih   << endl;
  cout << "     Fraction of toys inside 95%CL   = " << FracIn_95CL_ih   << endl;
  cout << "     Fraction of toys inside 2sigma = " << FracIn_2sigma_ih << endl;
  cout << "     Fraction of toys inside 3sigma = " << FracIn_3sigma_ih << endl;
  cout << endl;
  cout << " - Feldman-Cousins critical values:" << endl;
  cout << "     Crit Value 1sigma = " << FCCritVal_1sigma_ih << endl; 
  cout << "     Crit Value 90%CL   = " << FCCritVal_90CL_ih << endl; 
  cout << "     Crit Value 95%CL   = " << FCCritVal_95CL_ih << endl; 
  cout << "     Crit Value 2sigma = " << FCCritVal_2sigma_ih << endl; 
  cout << "     Crit Value 3sigma = " << FCCritVal_3sigma_ih << endl; 
  cout << endl;
  
  
  ofstream * outfilep = 0;
  outfilep = new ofstream(TString::Format("critval_%s_%1.3f_%s_%1.3f_%s.txt",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()),ios_base::trunc);  
  ofstream & outfile = *outfilep;
  outfile.close();
    
  return;
}


inline double GetConstChisq(bool Fit2D,string CL){
  
  const double CRITVAL_1SIGMA_1D = 1.;
  const double CRITVAL_90CL_1D =  2.71;
  const double CRITVAL_95CL_1D = 3.84;
  const double CRITVAL_2SIGMA_1D = 4.;
  const double CRITVAL_3SIGMA_1D = 9.;

  const double CRITVAL_1SIGMA_2D = 2.3;
  const double CRITVAL_90CL_2D = 4.61;
  const double CRITVAL_95CL_2D = 5.99;
  const double CRITVAL_2SIGMA_2D = 6.18;
  const double CRITVAL_3SIGMA_2D = 11.83;

  const string CL1sigma = "1sigma";
  const string CL90 = "90CL";
  const string CL95 = "95CL";
  const string CL2sigma = "2sigma";
  const string CL3sigma = "3sigma";
  
  if(!Fit2D){ // 1D fit
    if(CL==CL1sigma) return CRITVAL_1SIGMA_1D;
    if(CL==CL90)     return CRITVAL_90CL_1D;
    if(CL==CL95)     return CRITVAL_95CL_1D;
    if(CL==CL2sigma) return CRITVAL_2SIGMA_1D;
    if(CL==CL3sigma) return CRITVAL_3SIGMA_1D;
    else{
      cerr << endl;
      cerr << "The critical value for " << string << " is not available" << endl;
      cerr << "Possible choices are: " << endl;
      cerr << CL1sigma << endl;
      cerr << CL90 << endl;
      cerr << CL95 << endl;
      cerr << CL2sigma << endl;
      cerr << CL3sigma << endl;
      cerr << endl;
      exit(1);
    }
  }
  else{ // 2D fit
    if(CL==CL1sigma) return CRITVAL_1SIGMA_2D;
    if(CL==CL90)     return CRITVAL_90CL_2D;
    if(CL==CL95)     return CRITVAL_95CL_2D;
    if(CL==CL2sigma) return CRITVAL_2SIGMA_2D;
    if(CL==CL3sigma) return CRITVAL_3SIGMA_2D;
    else{
      cerr << endl;
      cerr << "The critical value for " << string << " is not available" << endl;
      cerr << "Possible choices are: " << endl;
      cerr << CL1sigma << endl;
      cerr << CL90 << endl;
      cerr << CL95 << endl;
      cerr << CL2sigma << endl;
      cerr << CL3sigma << endl;
      cerr << endl;
      exit(1);
    }
  }
  return -999;
}


double GetPvalue_Histo_rightside(TH1D *histo,double data){
  int bin_data = histo->FindBin(data);
  double pvalue = histo->Integral(bin_data,histo->GetNbinsX());
  pvalue /= histo->Integral();
  
  return pvalue;
}

double GetPvalue_Histo_leftside(TH1D *histo,double data){
  int bin_data = histo->FindBin(data);
  double pvalue = histo->Integral(1,bin_data);
  pvalue /= histo->Integral();
  
  return pvalue;
}


double GetPvalue_rightside(vector<double> vec_toys,double data){

  double ntot_toys = vec_toys.size();

  //cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
  //cout << "Input total number of toys = " << ntot_toys << endl; 
  //cout << "Set total number of toys to vector size: " << ntot_toys << endl;
  //cout << endl;

  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it >= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}



double GetPvalue_posterior_rightside(vector<double> vec_toys,vector<double> weight,double data){

  double ntot_toys = vec_toys.size();

  //cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
  //cout << "Input total number of toys = " << ntot_toys << endl; 
  //cout << "Set total number of toys to vector size: " << ntot_toys << endl;
  //cout << endl;

  // double nfrac=0;
  // for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
  //   if(*it >= data) nfrac += 1.; // If toy=data include the entry in pvalue integral range 
  // }
  // double pvalue =  nfrac / (double)ntot_toys;
  
  double tot = 0.; 
  double nfrac=0;
  for(int i=0;i<vec_toys.size();i++){
    tot += weight[i];
    double cont = vec_toys[i];
    if(cont >= data) nfrac += weight[i];
  }
  double pvalue =  nfrac / tot;

  return pvalue;
}

double CalcDLcrit(vector<double> dchisq, double conflev, double sigmacl){

  // Get error on critical value: sigmacl = sqrt(N*p*(1-p))

  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector   
  vector<double>::iterator dchisqIter;
  dchisqIter = dchisq.begin() + (int)(0.01*conflev*dchisq.size() + sigmacl);
  double critval = (*dchisqIter);
  return critval;
}

double CalcDLcrit_posterior(vector<double> dchisq, vector<double> weight, double conflev, double sigmacl){ 

  // Get error on coverage: sigmacl = sqrt(N*p*(1-p))/N

  int nelem = dchisq.size();
  
  double *new_dchisq = new double[nelem];
  double *new_weight = new double[nelem];
  int *index = new int[nelem];
  
  // copy vector into array
  for(int i=0;i<nelem;i++){
    new_dchisq[i] = dchisq[i];
    new_weight[i] = weight[i];
    index[i] = i;
  }
  
  sort(dchisq.begin(), dchisq.end());
  
  // store the sorted indices --> use for both dchisq and weights
  //index = TMath::Sort(nelem,new_dchisq,index,0);
  TMath::Sort(nelem,new_dchisq,index,0);

  double tot_weight = 0.;
  for(int i=0;i<nelem;i++) tot_weight += new_weight[index[i]]; 

  //cout << "conflev = " << conflev << endl;
  
  double sum = 0.;
  double critval = -999.;
  for(int i=0;i<nelem;i++){
    sum += new_weight[index[i]]; 
    
    //cout << std::setprecision(16);
    //cout << sum << " ";
    
    double epsilon = 0.000001; // for the double precision in the break condition

    if( sum >= ( (0.01 * conflev + sigmacl)*tot_weight - epsilon ) ){
          
      // cout << endl;
      // cout << "sum = " << sum << endl;
      // cout << "new_weight[index[i]] = " << new_weight[index[i]] << endl;
      // cout << "before: sum = " << (sum - new_weight[index[i]]) << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      // cout << "current: sum = " << sum << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      // cout << "after: sum = " << (sum + new_weight[index[i]]) << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      
      critval = new_dchisq[index[i]]; 
      break;
    }
  }
    
  return critval;
}




double CalcDLcrit(TH1D *h_delta_true,double FCfrac_CL){
  
  //Feldman-Cousins: Calculate chi2 value at 90%CL************************
  double ntoys_total = 0;
  double ntoys = 0;
  double ntoys_frac = 0;
  int chi2_bin = 0;
  double chi2_val = 0;
  double dlval = 0;
  
  ntoys_total = h_delta_true->Integral();

  for(int ibin=0;ibin<h_delta_true->GetNbinsX();ibin++){
    ntoys += h_delta_true->GetBinContent(ibin+1);
    ntoys_frac = ntoys/ntoys_total;
    
    if(ntoys_frac>=FCfrac_CL) {
      
      // chi2_bin = ibin+1;  // histgram rule for the # of bins!!!
      // chi2_val = chi2_bin*h_delta_true->GetBinWidth(1) + MinDL;
      // cout << "chi2_val = " << chi2_val << endl;

      // chi2_bin = ibin;  // histogram rule for the # of bins!!!
      // chi2_val = chi2_bin*h_delta_true->GetBinWidth(1)  
      // 	+ h_delta_true->GetXaxis()->GetBinLowEdge(1);

      chi2_val = h_delta_true->GetBinCenter(ibin+1);
 
      break;
    }
  }
  return chi2_val;
}

double calc_chisq(TH2D* hobs, TH2D* hexp){
  int nbins_x = hobs->GetNbinsX();
  int nbins_y = hobs->GetNbinsY();
  double chisq = 0;
  
  for(int i=1;i<= nbins_x;i++)
    {
      for(int j=1;j<= nbins_y;j++)
	{
	  double obs = hobs->GetBinContent(i,j);
	  double exp = hexp->GetBinContent(i,j);
	  chisq += valor::util::ChisqPoisson(obs,exp);
	}
    } 
  return chisq;
}



double get_error(double obs, TH1D* _hbeta, TH1D* _hbetaUnscaled)
{
  int obs_bin = _hbeta->FindBin(obs);
  double pvalue = _hbeta->Integral(obs_bin,_hbeta->GetXaxis()->GetNbins());
  double err_sqrd = 0;
  //  for(int i=1;i<obs_bin;i++)
  //  {
  //    double pd = _hbeta->Integral();
  //    double wl = _hbeta->GetBinContent(i)/_hbetaUnscaled->GetBinContent(i);
  //    double bl = _hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wl* (pvalue) / pd, 2) * bl;
  //    err_sqrd = err_sqrd + err;
  //  }
  
  //  for(int i=obs_bin;i<=_hbeta->GetXaxis()->GetNbins();i++)
  //  {
  //    double pd = _hbeta->Integral();
  //    double wr = _hbeta->GetBinContent(i)/_hbetaUnscaled->GetBinContent(i);
  //    double br = _hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wr* (1 - pvalue) / pd, 2) * br;
  //    err_sqrd = err_sqrd + err;
  //  }

  
  //Left error
  double pd = _hbeta->Integral();
  double wl = _hbeta->Integral(1,obs_bin)/_hbetaUnscaled->Integral(1,obs_bin);
  double bl = _hbetaUnscaled->Integral(1,obs_bin);
  double err = TMath::Power(wl* (pvalue) / (pd), 2) * bl;
  err_sqrd = err_sqrd + err;
  
  //Right error
  double pd = _hbeta->Integral();
  double wr = _hbeta->Integral(obs_bin,_hbeta->GetXaxis()->GetNbins())/_hbetaUnscaled->Integral(obs_bin,_hbeta->GetXaxis()->GetNbins());
  double br = _hbetaUnscaled->Integral(obs_bin,_hbeta->GetXaxis()->GetNbins());
  double err = TMath::Power(wr* (1 - pvalue) / (pd), 2) * br;
  err_sqrd = err_sqrd + err;
  
  double final_error = TMath::Sqrt(err_sqrd);
  
  return final_error;
  
}

double get_pvalue_error_vector(vector<double> vec_dchisq,vector<double> vec_weight, double data_critval, double pvalue) // takes the p-value not CL!!!
{
  double error = 0;
  double total_right=0;
  double total_left=0;
  
  std::vector<double> _right_tail_weight;
  std::vector<double> _left_tail_weight;

  for(int i=0;i<vec_dchisq.size();i++){
    
    double weight = vec_weight[i]; 

    if(vec_dchisq[i]<data_critval){
      _left_tail_weight.push_back(weight);
    }
    else{
      _right_tail_weight.push_back(weight);
    }
  }
  
  for(int i=0;i<_right_tail_weight.size();i++){
    error+= TMath::Power(_right_tail_weight[i]*(1-pvalue),2);
    total_right+= _right_tail_weight[i];
  }
  
  for(int i=0;i<_left_tail_weight.size();i++){
    error+= TMath::Power(_left_tail_weight[i]*(pvalue),2);
    total_left+= _left_tail_weight[i];
  }
  
  error = error / (TMath::Power(total_left + total_right ,2));
  return TMath::Sqrt(error);  
}

