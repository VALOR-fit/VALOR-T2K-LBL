
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

void analyze_run_fc_marginalization_exe
( 
 string namepar0 = "#delta_{CP}",
 string namepar1 = "sin^{2} #theta_{13}",

 const bool UseFullChisqFit = true,
 const bool CheckSingleToy = false,
 const int  ToySel = 0,
 const bool PrintToy = false,
 
 //char infilename[300] = "fit_prova_fc_nexp_2_0.000_-999.000-job0-_0_1.root"
 //char infilename[300] = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/fit_fitmarg_fc_dcp_NH_nexp_20_-1.571.root"
 char infilename[300] = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/fit_fitmarg_fc_dcp_NH_nexp_20_0.000.root",
 //char infilename[300] = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/fit_fitmarg_fc_dcp_NH_nexp_20_1.571.root,

 const double mydeltachisq = 3.55154
 

  )
{
  char name[300];
  FitParamList * fparm  = util::NamedParameterList("nuosc3f"); 
  PhysParamList * iparm = new PhysParamList;
  vector<double> myvector;
  vector<double> myvector_nh;
  vector<double> myvector_ih;
  
  double deltachisq0BF_nh = conv::kBadResult; 
  double deltachisq0BF_ih = conv::kBadResult; 

  TFile *filein = new TFile(infilename,"READ");
  if(!filein->IsOpen()){
    cerr << endl;
    cerr << "File " << infilename << " is not open" << endl;
    cerr << endl;
  }

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
  TH1D *hErecoToy_sample[conv::kns] = {};
  TH1D *hErecoExp_1Rmu = new TH1D("hErecoExp_1Rmu","hErecoExp_1Rmu",1,0,1);
  TH1D *hErecoExp_1Re = new TH1D("hErecoExp_1Re","hErecoExp_1Re",1,0,1);
  TH1D *hErecoExp_RHC1Rmu = new TH1D("hErecoExp_RHC1Rmu","hErecoExp_RHC1Rmu",1,0,1);
  TH1D *hErecoExp_RHC1Re = new TH1D("hErecoExp_RHC1Re","hErecoExp_RHC1Re",1,0,1);
  TH1D *hErecoExp_sample[conv::kns] = {};
#else
  TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_sample[conv::kns] = {};
  TH2D *hErecoExp_1Rmu = new TH2D("hErecoExp2D_1Rmu","hErecoExp2D_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_1Re = new TH2D("hErecoExp2D_1Re","hErecoExp2D_1Re",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Rmu = new TH2D("hErecoExp2D_RHC1Rmu","hErecoExp2D_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Re = new TH2D("hErecoExp2D_RHC1Re","hErecoExp2D_RHC1Re",1,0,1,1,0,1); 
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

  // fit status
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

  // true variables
  treein->SetBranchAddress("truepar0", &truepar0 );
  treein->SetBranchAddress("truepar1", &truepar1 );
  treein->SetBranchAddress("trueMH",   &trueMH   );
  treein->SetBranchAddress("trueparname0",   &trueparname0   );
  treein->SetBranchAddress("trueparname1",   &trueparname1   );
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  
  // global fit output
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

  // normal hierarchy fit output
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

  // inverted hierarchy fit output
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
#else  
  treein->SetBranchAddress("hErecoToy_1Rmu",     &hErecoToy_1Rmu     );
  treein->SetBranchAddress("hErecoToy_1Re",      &hErecoToy_1Re      );
  treein->SetBranchAddress("hErecoToy_RHC1Rmu",  &hErecoToy_RHC1Rmu  );
  treein->SetBranchAddress("hErecoToy_RHC1Re",   &hErecoToy_RHC1Re   );
  treein->SetBranchAddress("hErecoExp_1Rmu",     &hErecoExp_1Rmu     );
  treein->SetBranchAddress("hErecoExp_1Re",      &hErecoExp_1Re      );
  treein->SetBranchAddress("hErecoExp_RHC1Rmu",  &hErecoExp_RHC1Rmu  );
  treein->SetBranchAddress("hErecoExp_RHC1Re",   &hErecoExp_RHC1Re   );
#endif


  ///////////// Declare histograms //////////////

  TH1D *hbestfitFullChisq_1D = new TH1D("hbestfitFullChisq_1D","hbestfitFullChisq_1D",201,-3.15,3.15);
  TH2D *hbestfitFullChisq_2D = new TH2D("hbestfitFullChisq_2D","hbestfitFullChisq_2D",201,0,0.3,201,-3.15,3.15);

  TH1D *hbestfit_1D = new TH1D("hbestfit_1D","hbestfit_1D",201,-3.15,3.15);
  TH1D *hbestfit_1D_nh = new TH1D("hbestfit_1D_nh","hbestfit_1D_nh",201,-3.15,3.15);
  TH1D *hbestfit_1D_ih = new TH1D("hbestfit_1D_ih","hbestfit_1D_ih",201,-3.15,3.15);

  TH2D *hbestfit_2D = new TH2D("hbestfit_2D","hbestfit_2D",201,0,0.3,201,-3.15,3.15);
  TH2D *hbestfit_2D_nh = new TH2D("hbestfit_2D_nh","hbestfit_2D_nh",201,0,0.3,201,-3.15,3.15);
  TH2D *hbestfit_2D_ih = new TH2D("hbestfit_2D_ih","hbestfit_2D_ih",201,0,0.3,201,-3.15,3.15);

  TH1D *hdeltachisq = new TH1D("hdeltachisq","hdeltachisq",100,-5,20);
  TH1D *hdeltachisq_nh = new TH1D("hdeltachisq_nh","hdeltachisq_nh",100,-5,20);
  TH1D *hdeltachisq_ih = new TH1D("hdeltachisq_ih","hdeltachisq_ih",100,-5,20);

  TH1D *hchisqmin = new TH1D("hchisqmin","hchisqmin",100,-10,100); 
  TH1D *hchisqminDiff = new TH1D("hchisqminDiff","hchisqminDiff",100,-0.3,3.3);
  TH1D *hchisqminFullChisq = new TH1D("hchisqminFullChisq","hchisqminFullChisq",100,-10,100);
  TH1D *hnchisqcalls = new TH1D("hnchisqcalls","hnchisqcalls",100,0,700);
  TH1D *hfitstatus = new TH1D("hfitstatus","hfitstatus",5,0,5);

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

  int Nin_1sigma_nh = 0;
  int Nin_90CL_nh   = 0;
  int Nin_95CL_nh   = 0;
  int Nin_2sigma_nh = 0;
  
  int Nin_1sigma_ih = 0;
  int Nin_90CL_ih   = 0;
  int Nin_95CL_ih   = 0;
  int Nin_2sigma_ih = 0;

  int NdiffChisq = 0;
  
  int Nentries = treein->GetEntries();
  Nentries = 10;

  cout << endl;
  cout << "Total # of toys = " << Nentries << endl;
  cout << endl;

  for(int ientry=0;ientry<Nentries;ientry++){

    if(CheckSingleToy){ 
      Nentries = 1;
      ientry = ToySel;
    }

    treein->GetEntry(ientry);
 

     // Print output
    
    if(PrintToy){

      cout << "/////////////////" << endl;
      cout << endl;
      cout << "Toy " << ientry << endl;
      cout << "ToyID " << ToyID << endl;
      cout << endl;
      
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
    
   



    // Find bestfit chisq and true chisq

    //cout << *fparm << endl;
        
    double deltachisq0BF_used = deltachisq0BF_glob;
    
    //if(!UseFullChisqFit) deltachisq0BF_used = chisqtrue0-chisqmin;
    //else                 deltachisq0BF_used = chisqtrue0-chisqminFullChisq;
    // // Print chisq
    // int npt_chisq_x = hchisq1D_out->GetNbinsX();
    // for(int ix=0;ix<npt_chisq_x;ix++){
    //   cout << hchisq1D_out->GetBinContent(ix+1) << " ";
    // }

    // Fill deltaChisq histogram
    
    if(*trueMH=="normal"){
      deltachisq0BF_nh = chisqtrue0 - chisqmin_nh;
      hdeltachisq_nh->Fill(deltachisq0BF_nh);
      myvector_nh.push_back(deltachisq0BF_nh);
      hbestfit_2D_nh->Fill(bestfit0_nh,bestfit1_nh);
      hbestfit_1D_nh->Fill(bestfit0_nh);
      // Fill variables for coverage study
      if(deltachisq0BF_nh <= GetConstChisq(IsFit2D,"1sigma") ) Nin_1sigma_nh++;
      if(deltachisq0BF_nh <= GetConstChisq(IsFit2D,"90CL")   ) Nin_90CL_nh++;
      if(deltachisq0BF_nh <= GetConstChisq(IsFit2D,"95CL")   ) Nin_95CL_nh++;
      if(deltachisq0BF_nh <= GetConstChisq(IsFit2D,"2sigma") ) Nin_2sigma_nh++;
    }
    else if(*trueMH=="inverted"){
      deltachisq0BF_ih = chisqtrue0 - chisqmin_ih;
      hdeltachisq_ih->Fill(deltachisq0BF_ih);
      myvector_ih.push_back(deltachisq0BF_ih);
      hbestfit_2D_ih->Fill(bestfit0_ih,bestfit1_ih);
      hbestfit_1D_ih->Fill(bestfit0_ih);
      // Fill variables for coverage study
      if(deltachisq0BF_ih <= GetConstChisq(IsFit2D,"1sigma") ) Nin_1sigma_ih++;
      if(deltachisq0BF_ih <= GetConstChisq(IsFit2D,"90CL")   ) Nin_90CL_ih++;
      if(deltachisq0BF_ih <= GetConstChisq(IsFit2D,"95CL")   ) Nin_95CL_ih++;
      if(deltachisq0BF_ih <= GetConstChisq(IsFit2D,"2sigma") ) Nin_2sigma_ih++;
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
    // Fill variables for coverage study
    if(deltachisq0BF_used <= GetConstChisq(IsFit2D,"1sigma") ) Nin_1sigma++;
    if(deltachisq0BF_used <= GetConstChisq(IsFit2D,"90CL")   ) Nin_90CL++;
    if(deltachisq0BF_used <= GetConstChisq(IsFit2D,"95CL")   ) Nin_95CL++;
    if(deltachisq0BF_used <= GetConstChisq(IsFit2D,"2sigma") ) Nin_2sigma++;

    //cout << deltachisq0BF_used << " " << deltachisq0BF_nh << endl;
    
    // Fill histograms for fit status
    //hfitstatus->Fill(fitstatus);
    //hnchisqcalls->Fill(nchisqcalls);
    
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
    
  } // end loop over tree entries

  // Calculate coverage
  double FracIn_1sigma = (double) Nin_1sigma / (double) Nentries;
  double FracIn_90CL   = (double) Nin_90CL   / (double) Nentries;
  double FracIn_95CL   = (double) Nin_95CL   / (double) Nentries;
  double FracIn_2sigma = (double) Nin_2sigma / (double) Nentries;  

  double FracIn_1sigma_nh = (double) Nin_1sigma_nh / (double) Nentries;
  double FracIn_90CL_nh   = (double) Nin_90CL_nh   / (double) Nentries;
  double FracIn_95CL_nh   = (double) Nin_95CL_nh   / (double) Nentries;
  double FracIn_2sigma_nh = (double) Nin_2sigma_nh / (double) Nentries;  

  double FracIn_1sigma_ih = (double) Nin_1sigma_ih / (double) Nentries;
  double FracIn_90CL_ih   = (double) Nin_90CL_ih   / (double) Nentries;
  double FracIn_95CL_ih   = (double) Nin_95CL_ih   / (double) Nentries;
  double FracIn_2sigma_ih = (double) Nin_2sigma_ih / (double) Nentries;  
  
  // Calculate p-value

  double pvalue = GetPvalue_rightside(myvector,mydeltachisq);
 
  // Calculate the quantile (p-value = gauss area of both tails --> divide by 2)
  
  TF1 f2("f2","ROOT::Math::gaussian_quantile(1.-(x/2.),1.) ",0.,1.);
  double zscore = f2.Eval(pvalue); 
  
  // Calculate critical value
  
  double FCCritVal_1sigma = CalcDLcrit(myvector,CLfrac_1sigma*100.,0.); 
  double FCCritVal_90CL   = CalcDLcrit(myvector,CLfrac_90CL*100.,0.); 
  double FCCritVal_95CL   = CalcDLcrit(myvector,CLfrac_95CL*100.,0.); 
  double FCCritVal_2sigma = CalcDLcrit(myvector,CLfrac_2sigma*100.,0.); 

  double FCCritVal_1sigma_nh = conv::kBadResult;
  double FCCritVal_90CL_nh = conv::kBadResult;
  double FCCritVal_95CL_nh = conv::kBadResult;
  double FCCritVal_2sigma_nh = conv::kBadResult;
  double FCCritVal_1sigma_ih = conv::kBadResult;
  double FCCritVal_90CL_ih = conv::kBadResult;
  double FCCritVal_95CL_ih = conv::kBadResult;
  double FCCritVal_2sigma_ih = conv::kBadResult;

  if(*trueMH=="normal"){
    FCCritVal_1sigma_nh = CalcDLcrit(myvector_nh,CLfrac_1sigma*100.,0.); 
    FCCritVal_90CL_nh   = CalcDLcrit(myvector_nh,CLfrac_90CL*100.,0.); 
    FCCritVal_95CL_nh   = CalcDLcrit(myvector_nh,CLfrac_95CL*100.,0.); 
    FCCritVal_2sigma_nh = CalcDLcrit(myvector_nh,CLfrac_2sigma*100.,0.); 
  }
  else if(*trueMH=="inverted"){
    FCCritVal_1sigma_ih = CalcDLcrit(myvector_ih,CLfrac_1sigma*100.,0.); 
    FCCritVal_90CL_ih   = CalcDLcrit(myvector_ih,CLfrac_90CL*100.,0.); 
    FCCritVal_95CL_ih   = CalcDLcrit(myvector_ih,CLfrac_95CL*100.,0.); 
    FCCritVal_2sigma_ih = CalcDLcrit(myvector_ih,CLfrac_2sigma*100.,0.); 
  }


  // Draw histograms

  // TCanvas *cnchisqcalls = new TCanvas("cnchisqcalls","cnchisqcalls");
  // hnchisqcalls->GetYaxis()->SetTitle("# entries");
  // hnchisqcalls->GetXaxis()->SetTitle("# of #chi^{2} computations");
  // hnchisqcalls->DrawCopy();

  // TCanvas *cfitstatus = new TCanvas("cfitstatus","cfitstatus");
  // hfitstatus->GetYaxis()->SetTitle("# entries");
  // hfitstatus->GetXaxis()->SetTitle("Fit Status");
  // hfitstatus->DrawCopy();

  TCanvas *cbestfit_nh = new TCanvas("hbestfit_nh","hbestfit_nh");
  if(IsFit2D){
    hbestfit_2D_nh->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_2D_nh->GetYaxis()->SetTitle(namepar1.c_str());
    hbestfit_2D_nh->DrawCopy("colz");
  }
  else{
    hbestfit_1D_nh->GetYaxis()->SetTitle("# entries");
    hbestfit_1D_nh->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_1D_nh->DrawCopy();
  }

  TCanvas *cbestfit_ih = new TCanvas("hbestfit_ih","hbestfit_ih");
  if(IsFit2D){
    hbestfit_2D_ih->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_2D_ih->GetYaxis()->SetTitle(namepar1.c_str());
    hbestfit_2D_ih->DrawCopy("colz");
  }
  else{
    hbestfit_1D_ih->GetYaxis()->SetTitle("# entries");
    hbestfit_1D_ih->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_1D_ih->DrawCopy();
  }

  TCanvas *cbestfit = new TCanvas("hbestfit","hbestfit");
  if(IsFit2D){
    hbestfit_2D->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_2D->GetYaxis()->SetTitle(namepar1.c_str());
    hbestfit_2D->DrawCopy("colz");
  }
  else{
    hbestfit_1D->GetYaxis()->SetTitle("# entries");
    hbestfit_1D->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfit_1D->DrawCopy();
  }

  TCanvas *cbestfitFullChisq = new TCanvas("hbestfitFullChisq","hbestfitFullChisq");
  if(IsFit2D){
    hbestfitFullChisq_2D->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfitFullChisq_2D->GetYaxis()->SetTitle(namepar1.c_str());
    hbestfitFullChisq_2D->DrawCopy("colz");
  }
  else{
    hbestfitFullChisq_1D->GetYaxis()->SetTitle("# entries");
    hbestfitFullChisq_1D->GetXaxis()->SetTitle(namepar0.c_str());
    hbestfitFullChisq_1D->DrawCopy();
  }

  TCanvas *cchisqminDiff = new TCanvas("hchisqminDiff","hchisqminDiff");
  hchisqminDiff->GetYaxis()->SetTitle("# entries");
  hchisqminDiff->GetXaxis()->SetTitle("#chi^{2}_{Descent} - #chi^{2}_{Global}");
  hchisqminDiff->DrawCopy();


  TCanvas *cdeltachisq = new TCanvas("cdeltachisq","cdeltachisq");
  hdeltachisq->GetYaxis()->SetTitle("# entries");
  hdeltachisq->GetXaxis()->SetTitle("#Delta #chi^{2}_{Descent} (no global)");
  hdeltachisq->DrawCopy();
  TLine *lcritval_90CL = new TLine(FCCritVal_90CL,0,FCCritVal_90CL,hdeltachisq->GetMaximum());
  TLine *lcritval_2sigma = new TLine(FCCritVal_2sigma,0,FCCritVal_2sigma,hdeltachisq->GetMaximum());
  lcritval_90CL->SetLineColor(kRed);
  lcritval_2sigma->SetLineColor(kBlue);
  lcritval_90CL->SetLineWidth(2);
  lcritval_2sigma->SetLineWidth(2);
  lcritval_90CL->Draw();
  lcritval_2sigma->Draw();
  TLine *ldeltachisq_data = new TLine(mydeltachisq,0,mydeltachisq,hdeltachisq->GetMaximum());
  ldeltachisq_data->SetLineColor(kGreen);
  ldeltachisq_data->SetLineWidth(2);
  if(mydeltachisq>0.) ldeltachisq_data->Draw();
  TLegend *ldeltachisq = new TLegend(0.6,0.45,0.8,0.75);
  ldeltachisq->SetFillColor(0);
  ldeltachisq->SetLineColor(0);
  ldeltachisq->SetTextFont(132);
  ldeltachisq->SetTextSize(0.07);
  ldeltachisq->AddEntry(lcritval_90CL,TString::Format("#Delta #chi^{2}_{crit}(90CL) = %1.2f",FCCritVal_90CL).Data(),"l");
  ldeltachisq->AddEntry(lcritval_2sigma,TString::Format("#Delta #chi^{2}_{crit} (2#sigma) = %1.2f",FCCritVal_2sigma).Data(),"l");
  if(mydeltachisq>0.) ldeltachisq->AddEntry(ldeltachisq_data,TString::Format("#Delta #chi^{2}_{data} = %1.3f",mydeltachisq).Data(),"l");
  ldeltachisq->Draw();
  TLatex *latPvalue = new TLatex(hdeltachisq->GetXaxis()->GetBinLowEdge(5),hdeltachisq->GetMaximum()*0.7, TString::Format("P-value = %1.3f",pvalue).Data() );
  latPvalue->Draw("same");
  TLatex *latZscore = new TLatex(hdeltachisq->GetXaxis()->GetBinLowEdge(5),hdeltachisq->GetMaximum()*0.7, TString::Format("Z = %1.2f#sigma",zscore).Data() );
  latZscore->Draw("same");


  TCanvas *cdeltachisq_nh = new TCanvas("cdeltachisq_nh","cdeltachisq_nh");
  hdeltachisq_nh->GetYaxis()->SetTitle("# entries");
  hdeltachisq_nh->GetXaxis()->SetTitle("#Delta #chi^{2}_{Descent} (no global)");
  hdeltachisq_nh->DrawCopy();
  TLine *lcritval_90CL_nh;
  TLine *lcritval_2sigma_nh;
  if(*trueMH=="normal"){
    lcritval_90CL_nh = new TLine(FCCritVal_90CL_nh,0,FCCritVal_90CL_nh,hdeltachisq_nh->GetMaximum());
    lcritval_2sigma_nh = new TLine(FCCritVal_2sigma_nh,0,FCCritVal_2sigma_nh,hdeltachisq_nh->GetMaximum());
    lcritval_90CL_nh->SetLineColor(kRed);
    lcritval_2sigma_nh->SetLineColor(kBlue);
    lcritval_90CL_nh->Draw();
    lcritval_2sigma_nh->Draw();
  }
  if(mydeltachisq>0.) ldeltachisq_data->Draw();
  TLegend *ldeltachisq_nh = new TLegend(0.6,0.45,0.8,0.75);
  ldeltachisq_nh->SetFillColor(0);
  ldeltachisq_nh->SetLineColor(0);
  ldeltachisq_nh->SetTextFont(132);
  ldeltachisq_nh->SetTextSize(0.07);
  ldeltachisq_nh->AddEntry(lcritval_90CL_nh,TString::Format("#Delta #chi^{2}_{crit}(90CL) = %1.2f",FCCritVal_90CL_nh).Data(),"l");
  ldeltachisq_nh->AddEntry(lcritval_2sigma_nh,TString::Format("#Delta #chi^{2}_{crit} (2#sigma) = %1.2f",FCCritVal_2sigma_nh).Data(),"l");
  if(mydeltachisq>0.) ldeltachisq_nh->AddEntry(ldeltachisq_data,TString::Format("#Delta #chi^{2}_{data} = %1.3f",mydeltachisq).Data(),"l");
  ldeltachisq_nh->Draw();


  TCanvas *cdeltachisq_ih = new TCanvas("cdeltachisq_ih","cdeltachisq_ih");
  hdeltachisq_ih->GetYaxis()->SetTitle("# entries");
  hdeltachisq_ih->GetXaxis()->SetTitle("#Delta #chi^{2}_{Descent} (no global)");
  hdeltachisq_ih->DrawCopy();
  TLine *lcritval_90CL_ih;
  TLine *lcritval_2sigma_ih;
  if(*trueMH=="inverted"){
    lcritval_90CL_ih = new TLine(FCCritVal_90CL_ih,0,FCCritVal_90CL_ih,hdeltachisq_ih->GetMaximum());
    lcritval_2sigma_ih = new TLine(FCCritVal_2sigma_ih,0,FCCritVal_2sigma_ih,hdeltachisq_ih->GetMaximum());
    lcritval_90CL_ih->SetLineColor(kRed);
    lcritval_2sigma_ih->SetLineColor(kBlue);
    lcritval_90CL_ih->Draw();
    lcritval_2sigma_ih->Draw();
  }
  if(mydeltachisq>0.) ldeltachisq_data->Draw();
  TLegend *ldeltachisq_ih = new TLegend(0.6,0.45,0.8,0.75);
  ldeltachisq_ih->SetFillColor(0);
  ldeltachisq_ih->SetLineColor(0);
  ldeltachisq_ih->SetTextFont(132);
  ldeltachisq_ih->SetTextSize(0.07);
  ldeltachisq_ih->AddEntry(lcritval_90CL_ih,TString::Format("#Delta #chi^{2}_{crit}(90CL) = %1.2f",FCCritVal_90CL_ih).Data(),"l");
  ldeltachisq_ih->AddEntry(lcritval_2sigma_ih,TString::Format("#Delta #chi^{2}_{crit} (2#sigma) = %1.2f",FCCritVal_2sigma_ih).Data(),"l");
  if(mydeltachisq>0.) ldeltachisq_ih->AddEntry(ldeltachisq_data,TString::Format("#Delta #chi^{2}_{data} = %1.3f",mydeltachisq).Data(),"l");
  ldeltachisq_ih->Draw();



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
      hchisq1D_out->GetXaxis()->SetTitle(namepar0.c_str());
      hchisq1D_out->Draw("histo");
      npt_chisq_x = hchisq1D_out->GetNbinsX();
      npt_chisq_y = hchisq1D_out->GetNbinsY();
    }
    TCanvas *cdeltachisq_out_GetEntry = new TCanvas("cdeltachisq_out_GetEntry","cdeltachisq_out_GetEntry");
    if(IsFit2D){
      hdeltachisq2D_out->GetXaxis()->SetTitle(namepar0.c_str());
      hdeltachisq2D_out->GetYaxis()->SetTitle(namepar1.c_str());
      hdeltachisq2D_out->Draw("colz");
    }    
    else{
      hdeltachisq1D_out->GetXaxis()->SetTitle(namepar0.c_str());
      hdeltachisq1D_out->Draw("histo");
    }
    cout << endl;
    
    cout << "npt_chisq_x = " << npt_chisq_x << endl;
    cout << "npt_chisq_y = " << npt_chisq_y << endl;
    
    for(int iy=0;iy<npt_chisq_y;iy++){
      for(int ix=0;ix<npt_chisq_x;ix++){
       	if(IsFit2D) cout << hchisq2D_out->GetBinContent(ix+1,iy+1) << " ";
        else        cout << hchisq1D_out->GetBinContent(ix+1) << " ";
      }
      cout << endl;
    }
    cout << endl;
    TCanvas *cErecoToy_1Rmu = new TCanvas("cErecoToy_1Rmu","cErecoToy_1Rmu");
    hErecoToy_1Rmu->Draw("colz");
    TCanvas *cErecoToy_1Re = new TCanvas("cErecoToy_1Re","cErecoToy_1Re");
    hErecoToy_1Re->Draw("colz");
    TCanvas *cErecoToy_RHC1Rmu = new TCanvas("cErecoToy_RHC1Rmu","cErecoToy_RHC1Rmu");
    hErecoToy_RHC1Rmu->Draw("colz");
    TCanvas *cErecoToy_RHC1Re = new TCanvas("cErecoToy_RHC1Re","cErecoToy_RHC1Re");
    hErecoToy_RHC1Re->Draw("colz");
   
    cout << endl;
    cout << "Events in fake data spectra:"<< endl;
    cout << " - # of 1Rmu events = "    << hErecoToy_1Rmu->Integral() << endl;
    cout << " - # of RHC1Rmu events = " << hErecoToy_RHC1Rmu->Integral() << endl;
    cout << " - # of 1Re events = "     << hErecoToy_1Re->Integral() << endl;
    cout << " - # of RHC1Re events = "  << hErecoToy_RHC1Re->Integral() << endl;
    cout << endl;
  }

  

  // Print plots
  
  string storedparname0 = *trueparname0;
  string storedparname1 = *trueparname1;
  if(storedparname1=="") storedparname1 = "1D";
  string storedmh       = *trueMH;
  double storetrue0      = truepar0;
  double storetrue1      = truepar1;
  
  cout << storetrue0 << " " << storetrue1 << endl;

  cdeltachisq->Print( TString::Format("deltchisq_glob_%s_%1.3f_%s_%1.3f_%s.pdf",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  cdeltachisq_nh->Print( TString::Format("deltchisq_nh_%s_%1.3f_%s_%1.3f_%s.pdf",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  cdeltachisq_ih->Print( TString::Format("deltchisq_ih_%s_%1.3f_%s_%1.3f_%s.pdf",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
 
  cchisqminDiff->Print( TString::Format("cchisqminDiff_%s_%1.3f_%s_%1.3f_%s.pdf",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  
  cbestfit->Print( TString::Format("cbestfit_%s_%1.3f_%s_%1.3f_%s.pdf",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  cbestfit_nh->Print( TString::Format("cbestfit_nh_%s_%1.3f_%s_%1.3f_%s.pdf",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()) );
  cbestfit_ih->Print( TString::Format("cbestfit_ih_%s_%1.3f_%s_%1.3f_%s.pdf",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str())  );

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
  cout << "/// Global Result ///" << endl;
  cout << " - P-value = " << pvalue << endl;
  cout << " - Z = " << zscore << " sigma" << endl;
  cout << " - Coverage:" << endl;
  cout << "     Fraction of toys inside 1sigma = " << FracIn_1sigma << endl;
  cout << "     Fraction of toys inside 90CL   = " << FracIn_90CL   << endl;
  cout << "     Fraction of toys inside 95CL   = " << FracIn_95CL   << endl;
  cout << "     Fraction of toys inside 2sigma = " << FracIn_2sigma << endl;
  cout << endl;
  cout << " - Feldman-Cousins critical values:" << endl;
  cout << "     Crit Value 1sigma = " << FCCritVal_1sigma << endl; 
  cout << "     Crit Value 90CL   = " << FCCritVal_90CL << endl; 
  cout << "     Crit Value 95CL   = " << FCCritVal_95CL << endl; 
  cout << "     Crit Value 2sigma = " << FCCritVal_2sigma << endl; 
  cout << endl;
  cout << "/// Normal Hierarchy Result ///" << endl;
  cout << " - Coverage:" << endl;
  cout << "     Fraction of toys inside 1sigma = " << FracIn_1sigma_nh << endl;
  cout << "     Fraction of toys inside 90CL   = " << FracIn_90CL_nh   << endl;
  cout << "     Fraction of toys inside 95CL   = " << FracIn_95CL_nh   << endl;
  cout << "     Fraction of toys inside 2sigma = " << FracIn_2sigma_nh << endl;
  cout << endl;
  cout << " - Feldman-Cousins critical values:" << endl;
  cout << "     Crit Value 1sigma = " << FCCritVal_1sigma_nh << endl; 
  cout << "     Crit Value 90CL = "   << FCCritVal_90CL_nh << endl; 
  cout << "     Crit Value 95CL = "   << FCCritVal_95CL_nh << endl; 
  cout << "     Crit Value 2sigma = " << FCCritVal_2sigma_nh << endl; 
  cout << endl;
  cout << "/// Inverted Hierarchy Result ///" << endl;
  cout << " - Coverage:" << endl;
  cout << "     Fraction of toys inside 1sigma = " << FracIn_1sigma_ih << endl;
  cout << "     Fraction of toys inside 90CL   = " << FracIn_90CL_ih   << endl;
  cout << "     Fraction of toys inside 95CL   = " << FracIn_95CL_ih   << endl;
  cout << "     Fraction of toys inside 2sigma = " << FracIn_2sigma_ih << endl;
  cout << endl;
  cout << " - Feldman-Cousins critical values:" << endl;
  cout << "     Crit Value 1sigma = " << FCCritVal_1sigma_ih << endl; 
  cout << "     Crit Value 90CL   = " << FCCritVal_90CL_ih << endl; 
  cout << "     Crit Value 95CL   = " << FCCritVal_95CL_ih << endl; 
  cout << "     Crit Value 2sigma = " << FCCritVal_2sigma_ih << endl; 
  cout << endl;

  
  ofstream * outfilep = 0;
  outfilep = new ofstream(TString::Format("critval_%s_%1.3f_%s_%1.3f_%s.pdf",storedparname0.c_str(),storetrue0,storedparname1.c_str(),storetrue1,storedmh.c_str()),ios_base::trunc);  
  ofstream & outfile = *outfilep;
  outfile << endl;
  outfile << endl;
  outfile << "/// True parameters ///" << endl;
  outfile << " - truepar0: " << *trueparname0 << " = " << truepar0 << endl;
  outfile << " - truepar1: " << *trueparname1 << " = " << truepar1 << endl;
  outfile << " - Mass Hierarchy: " << *trueMH << endl;
  outfile << endl;
  outfile << " - Difference between descent gradient and global minimum" << endl;
  outfile << " - # of toys with different minim and full chisq minimum: "
	  << NdiffChisq << " / " << Nentries << " toys" << endl;
  outfile << " - Mean difference = " << hchisqminDiff->GetMean() << endl;
  outfile << " - RMS difference = " << hchisqminDiff->GetRMS() << endl;
  outfile << endl;
  outfile << "/// Global Result ///" << endl;
  outfile << " - Coverage:" << endl;
  outfile << "     Fraction of toys inside 1sigma = " << FracIn_1sigma << endl;
  outfile << "     Fraction of toys inside 90CL   = " << FracIn_90CL   << endl;
  outfile << "     Fraction of toys inside 95CL   = " << FracIn_95CL   << endl;
  outfile << "     Fraction of toys inside 2sigma = " << FracIn_2sigma << endl;
  outfile << endl;
  outfile << " - Feldman-Cousins critical values:" << endl;
  outfile << "     Crit Value 1sigma = " << FCCritVal_1sigma << endl; 
  outfile << "     Crit Value 90CL   = " << FCCritVal_90CL << endl; 
  outfile << "     Crit Value 95CL   = " << FCCritVal_95CL << endl; 
  outfile << "     Crit Value 2sigma = " << FCCritVal_2sigma << endl; 
  outfile << endl;
  outfile << "/// Normal Hierarchy Result ///" << endl;
  outfile << " - Coverage:" << endl;
  outfile << "     Fraction of toys inside 1sigma = " << FracIn_1sigma_nh << endl;
  outfile << "     Fraction of toys inside 90CL   = " << FracIn_90CL_nh   << endl;
  outfile << "     Fraction of toys inside 95CL   = " << FracIn_95CL_nh   << endl;
  outfile << "     Fraction of toys inside 2sigma = " << FracIn_2sigma_nh << endl;
  outfile << endl;
  outfile << " - Feldman-Cousins critical values:" << endl;
  outfile << "     Crit Value 1sigma = " << FCCritVal_1sigma_nh << endl; 
  outfile << "     Crit Value 90CL = "   << FCCritVal_90CL_nh << endl; 
  outfile << "     Crit Value 95CL = "   << FCCritVal_95CL_nh << endl; 
  outfile << "     Crit Value 2sigma = " << FCCritVal_2sigma_nh << endl; 
  outfile << endl;
  outfile << "/// Inverted Hierarchy Result ///" << endl;
  outfile << " - Coverage:" << endl;
  outfile << "     Fraction of toys inside 1sigma = " << FracIn_1sigma_ih << endl;
  outfile << "     Fraction of toys inside 90CL   = " << FracIn_90CL_ih   << endl;
  outfile << "     Fraction of toys inside 95CL   = " << FracIn_95CL_ih   << endl;
  outfile << "     Fraction of toys inside 2sigma = " << FracIn_2sigma_ih << endl;
  outfile << endl;
  outfile << " - Feldman-Cousins critical values:" << endl;
  outfile << "     Crit Value 1sigma = " << FCCritVal_1sigma_ih << endl; 
  outfile << "     Crit Value 90CL   = " << FCCritVal_90CL_ih << endl; 
  outfile << "     Crit Value 95CL   = " << FCCritVal_95CL_ih << endl; 
  outfile << "     Crit Value 2sigma = " << FCCritVal_2sigma_ih << endl; 
  outfile << endl;
  outfile.close();
    
  return;
}


double GetConstChisq(bool Fit2D,string CL){
  
  const double CRITVAL_1SIGMA_1D = 1.;
  const double CRITVAL_90CL_1D =  2.71;
  const double CRITVAL_95CL_1D = 3.84;
  const double CRITVAL_2SIGMA_1D = 4.;

  const double CRITVAL_1SIGMA_2D = 2.3;
  const double CRITVAL_90CL_2D = 4.61;
  const double CRITVAL_95CL_2D = 5.99;
  const double CRITVAL_2SIGMA_2D = 6.18;

  const string CL1sigma = "1sigma";
  const string CL90 = "90CL";
  const string CL95 = "95CL";
  const string CL2sigma = "2sigma";
  
  if(!Fit2D){ // 1D fit
    if(CL==CL1sigma) return CRITVAL_1SIGMA_1D;
    if(CL==CL90)     return CRITVAL_90CL_1D;
    if(CL==CL95)     return CRITVAL_95CL_1D;
    if(CL==CL2sigma) return CRITVAL_2SIGMA_1D;
    else{
      cerr << endl;
      cerr << "The critical value for " << string << " is not available" << endl;
      cerr << "Possible choices are: " << endl;
      cerr << CL1sigma << endl;
      cerr << CL90 << endl;
      cerr << CL95 << endl;
      cerr << CL2sigma << endl;
      cerr << endl;
      exit(1);
    }
  }
  else{ // 2D fit
    if(CL==CL1sigma) return CRITVAL_1SIGMA_2D;
    if(CL==CL90)     return CRITVAL_90CL_2D;
    if(CL==CL95)     return CRITVAL_95CL_2D;
    if(CL==CL2sigma) return CRITVAL_2SIGMA_2D;
    else{
      cerr << endl;
      cerr << "The critical value for " << string << " is not available" << endl;
      cerr << "Possible choices are: " << endl;
      cerr << CL1sigma << endl;
      cerr << CL90 << endl;
      cerr << CL95 << endl;
      cerr << CL2sigma << endl;
      cerr << endl;
      exit(1);
    }
  }
  return -999;
}

double CalcDLcrit(vector<double> dchisq, double conflev, double sigmacl){
  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector   
  vector<double>::iterator dchisqIter;
  dchisqIter = dchisq.begin() + (int)(0.01*conflev*dchisq.size() + sigmacl);
  double critval = (*dchisqIter);
  return critval;
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
