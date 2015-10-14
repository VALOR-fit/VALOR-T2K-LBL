
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

using namespace std;
using namespace valor;
using namespace valor::conv;

void analyze_run_fit_marginalization
( 
 //char infilename[200] = "fit_submit_fitmarg_joint_dcp_nexp_2_-1.601_-999.000-job0-_0_1.root"
 //char infilename[200] = "fit_submit_fitmarg_joint_dcp_nexp_1_-1.601_-999.000-job0-_0_0.root"
 char infilename[200] = "fitmarg_joint_PROVA-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root"
)
{
  char name[200];
  FitParamList * fparm  = util::NamedParameterList("nuosc3f"); 
  PhysParamList * iparm = new PhysParamList;

  TFile *filein = new TFile(infilename,"READ");
  if(!filein->IsOpen()){
    cerr << endl;
    cerr << "File " << infilename << " is not open" << endl;
    cerr << endl;
  }

  TTree *treein = (TTree*)filein->Get("fitnt");
  if(!treein){
    cerr << endl;
    cerr << "Tree object is not empty" << endl;
    cerr << endl;
  }

  // Set Tree toy/fit
  int run=-999;
  string *tag;
  int nexpts=-999;
  int ToyID=-999;
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
  double deltachisq0BF=conv::kALargeNum;

  TH2D *hchisq2D_out = new TH2D("hchisq2D_out","hchisq2D_out",1,0,1,1,0,1);
  TH1D *hchisq1D_out = new TH1D("hchisq1D_out","hchisq1D_out",1,0,1);
  TH2D *hdeltachisq2D_out = new TH2D("hdeltachisq2D_out","hdeltachisq2D_out",1,0,1,1,0,1);
  TH1D *hdeltachisq1D_out = new TH1D("hdeltachisq1D_out","hdeltachisq1D_out",1,0,1);
  TH2D *hstep2D_out = new TH2D("hstep2D_out","hstep2D_out",1,0,1,1,0,1);
  TH1D *hstep1D_out = new TH1D("hstep1D_out","hstep1D_out",1,0,1);

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
  #ifndef __USE_2D_TEMPLATES__
    TH1D *hErecoToy_1Rmu = new TH1D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1);
    TH1D *hErecoToy_1Re = new TH1D("hErecoToy_1Re","hErecoToy_1Re",1,0,1);
    TH1D *hErecoToy_RHC1Rmu = new TH1D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1);
    TH1D *hErecoToy_RHC1Re = new TH1D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1);
  #else
    TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
    TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
    TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
    TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
  #endif
 
  // Create tree to store outputs
  treein->SetBranchAddress("best_fit_params",&fparm);
  treein->SetBranchAddress("input_params",&iparm);
  treein->SetBranchAddress("truepar0",&truepar0);
  treein->SetBranchAddress("truepar1",&truepar1);
  treein->SetBranchAddress("hchisq2D",&hchisq2D_out);
  treein->SetBranchAddress("hchisq1D",&hchisq1D_out);
  treein->SetBranchAddress("hdeltachisq2D",&hdeltachisq2D_out);
  treein->SetBranchAddress("hdeltachisq1D",&hdeltachisq1D_out);
  treein->SetBranchAddress("hstep2D",&hstep2D_out);
  treein->SetBranchAddress("hstep1D",&hstep1D_out);
  treein->SetBranchAddress("fitstatus",&fitstatus);
  treein->SetBranchAddress("chisqmin",&chisqmin);
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("deltachisqtrue",&deltachisqtrue);
  treein->SetBranchAddress("nfits",&nfits);
  treein->SetBranchAddress("nchisqcalls",&nchisqcalls);
  treein->SetBranchAddress("niterations",&niterations);
  treein->SetBranchAddress("bestfit0",&bestfit0);
  treein->SetBranchAddress("bestfit1",&bestfit1);
  treein->SetBranchAddress("bestfiterrNeg0",&bestfiterrNeg0);
  treein->SetBranchAddress("bestfiterrNeg1",&bestfiterrNeg1);
  treein->SetBranchAddress("bestfiterrPos0",&bestfiterrPos0);
  treein->SetBranchAddress("bestfiterrPos1",&bestfiterrPos1);
  treein->SetBranchAddress("bestfiterrSym0",&bestfiterrSym0);
  treein->SetBranchAddress("bestfiterrSym1",&bestfiterrSym1);
  treein->SetBranchAddress("bestfitFullChisq0", &bestfitFullChisq0 );
  treein->SetBranchAddress("bestfitFullChisq1", &bestfitFullChisq1 );
  treein->SetBranchAddress("chisqminFullChisq",&chisqminFullChisq);
  treein->SetBranchAddress("ToyID",&ToyID);
  treein->SetBranchAddress("Job",&run);
  treein->SetBranchAddress("Tag",&tag);
  treein->SetBranchAddress("nexpts",&nexpts);
  #ifndef __USE_2D_TEMPLATES__
    treein->SetBranchAddress("hErecoToy_1Rmu",&hErecoToy_1Rmu);
    treein->SetBranchAddress("hErecoToy_1Re",&hErecoToy_1Re);
    treein->SetBranchAddress("hErecoToy_RHC1Rmu",&hErecoToy_RHC1Rmu);
    treein->SetBranchAddress("hErecoToy_RHC1Re",&hErecoToy_RHC1Re);
  #else  
    treein->SetBranchAddress("hErecoToy_1Rmu",&hErecoToy_1Rmu);
    treein->SetBranchAddress("hErecoToy_1Re",&hErecoToy_1Re);
    treein->SetBranchAddress("hErecoToy_RHC1Rmu",&hErecoToy_RHC1Rmu);
    treein->SetBranchAddress("hErecoToy_RHC1Re",&hErecoToy_RHC1Re);
  #endif




  ///////////// Declare histograms //////////////
  ///////////////////////////////////////////////
  
  // Check if 1D or 2D fit
  bool IsFit2D;
  treein->GetEntry(1);
  if(bestfit1==-999) IsFit2D = false;
  else               IsFit2D = true;

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
  
  int Nentries = treein->GetEntries();
  for(int ientry=0;ientry<Nentries;ientry++){

    treein->GetEntry(ientry);
  
  } // end loop over tree entries

  



  /////////////// Get pulls //////////////
  ////////////////////////////////////////
  
  



  // Draw histograms



  //

  treein->GetEntry(1);
  
  TCanvas *cchisq_out_GetEntry = new TCanvas("cchisq_out_GetEntry","cchisq_out_GetEntry");
  //if(IsFit2D) hchisq2D_out->Draw("colz");
  //else        hchisq1D_out->Draw("histo");
  
  //hdeltachisq1D_out->Draw("histo");
  //hstep1D_out->Draw("histo");
  hchisq1D_out->Draw("histo");
  //hErecoToy_1Rmu->ProjectionY()->Draw("histo");
  //hErecoToy_1Re->ProjectionY()->Draw("histo");

  cout << endl;
  int npt_chisq_x = hchisq1D_out->GetNbinsX();
  for(int ix=0;ix<npt_chisq_x;ix++){
    cout << hchisq1D_out->GetBinContent(ix+1) << " ";
  }
  cout << endl;

  TCanvas *cErecoToy_1Rmu = new TCanvas("cErecoToy_1Rmu","cErecoToy_1Rmu");
  hErecoToy_1Rmu->Draw("colz");
  // TCanvas *cErecoToy_1Re = new TCanvas("cErecoToy_1Re","cErecoToy_1Re");
  // hErecoToy_1Re->Draw("colz");

  // Print results  
  cout << endl;
  cout << "truepar0 = " << truepar0 << endl;
  if(IsFit2D) cout << "truepar1 = " << truepar1 << endl;
  cout << endl;
  cout << "bestfit0 = " << bestfit0 << endl;
  if(IsFit2D) cout << "bestfit1 = " << bestfit1 << endl;
  cout << endl;
  
  
  return;
}
