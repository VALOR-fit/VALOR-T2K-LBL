
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

using namespace std;
using namespace valor;
using namespace valor::conv;

const double err_minos_par0[2] = {0,6.28};
const double err_minos_par1[2] = {0,6.28};
const double err_hesse_par0[2] = {0,20.};
const double err_hesse_par1[2] = {0,20.};

const double CLfrac_1sigma = 0.6827;  // 1sigma% CL
const double CLfrac_90CL   = 0.9;  // 90% CL
const double CLfrac_95CL   = 0.95;  // 95% CL


void analyze_run_fit_marginalization_exe
( 
 string namepar0 = "sin^{2} #theta_{13}",
 string namepar1 = "#delta_{CP}",

 bool UseFullChisqFit = true,
 bool CheckSingleToy = true,
 const int ToySel = 0, //488,

 //char infilename[300] = "../TECH_NOTE_14_3_16/BETABAR/fit_fitmarg_joint_betabar_react_margmh_10kToys_Run16_Asimov_A_nexp_1_0.500_-999.000-job0-_0_0.root"
 //char infilename[300] = "../TECH_NOTE_14_3_16/BETABAR/fit_fitmarg_joint_betabar_react_margmh_10kToys_Run16_Asimov_A_nexp_1_1.000_-999.000-job0-_0_0.root"
 //char infilename[300] = "../TECH_NOTE_14_3_16/BETABAR/fit_fitmarg_joint_betabar_react_margmh_10kToys_Run16_Asimov_A_nexp_1_1.500_-999.000-job0-_0_0.root"
 //char infilename[300] = "../TECH_NOTE_14_3_16/BETABAR/fit_fitmarg_joint_betabar_react_margmh_10kToys_Run16_Asimov_A_nexp_1_2.000_-999.000-job0-_0_0.root"
 //char infilename[300] = "../TECH_NOTE_14_3_16/BETABAR/fit_fitmarg_joint_betabar_react_margmh_10kToys_Run16_Asimov_A_nexp_1_2.500_-999.000-job0-_0_0.root"

 //char infilename[300] = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_nh_Data_Run16_nexp_1_-1.601_-999.000-job0-_0_0.root"
 //char infilename[300] = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth13_t2konly_nh_Data_Run16_nexp_1_0.022_-999.000-job0-_0_0.root"
 //char infilename[300] = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23_react_nh_Data_Run16_nexp_1_0.528_-999.000-job0-_0_0.root"
 //char infilename[300] = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dmsq_react_nh_Data_Run16_nexp_1_0.003_-999.000-job0-_0_0.root"
 
 char infilename[300] = "fit_fit_marg_PROVA_dcp_nexp_1_-1.601_-999.000-job0-_0_0.root"
  )
{
  char name[300];
  FitParamList * fparm  = util::NamedParameterList("nuosc3f"); 
  PhysParamList * iparm = new PhysParamList;
  vector<double> myvector;
  
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
  bool FlagFit2D=false;
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
  double deltachisq0BF=conv::kALargeNum;
  double deltachisq01=conv::kALargeNum;

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
  treein->SetBranchAddress("IsFit2D",&FlagFit2D);
  treein->SetBranchAddress("chisqmin",&chisqmin);
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("deltachisq0BF",&deltachisq0BF);
  treein->SetBranchAddress("deltachisq01",&deltachisq01);
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

  TH1D *hbestfitFullChisq_1D = new TH1D("hbestfitFullChisq_1D","hbestfitFullChisq_1D",201,-3.15,3.15);
  TH1D *hbestfit_1D = new TH1D("hbestfit_1D","hbestfit_1D",201,-3.15,3.15);
  TH2D *hbestfitFullChisq_2D = new TH2D("hbestfitFullChisq_2D","hbestfitFullChisq_2D",201,0,0.3,201,-3.15,3.15);
  TH2D *hbestfit_2D = new TH2D("hbestfit_2D","hbestfit_2D",201,0,0.3,201,-3.15,3.15);
  TH1D *hdeltachisqDiff = new TH1D("hdeltachisqDiff","hdeltachisqDiff",100,-5,20);
  TH1D *hchisqminDiff = new TH1D("hchisqminDiff","hchisqminDiff",100,-0.3,3.3);
  TH1D *hchisqmin = new TH1D("hchisqmin","hchisqmin",100,-10,100);
  TH1D *hchisqminFullChisq = new TH1D("hchisqminFullChisq","hchisqminFullChisq",100,-10,100);
  TH1D *hErrMINOSneg0 = new TH1D("hErrMINOSneg0","hErrMINOSneg0",1000,err_minos_par0[0],err_minos_par0[1]);
  TH1D *hErrMINOSneg1 = new TH1D("hErrMINOSneg1","hErrMINOSneg1",1000,err_minos_par1[0],err_minos_par1[1]);
  TH1D *hErrMINOSpos0 = new TH1D("hErrMINOSpos0","hErrMINOSpos0",1000,err_minos_par0[0],err_minos_par0[1]);
  TH1D *hErrMINOSpos1 = new TH1D("hErrMINOSpos1","hErrMINOSpos1",1000,err_minos_par1[0],err_minos_par1[1]);
  TH1D *hErrHESSEsym0 = new TH1D("hErrHESSEsym0","hErrHESSEsym0",1000,err_hesse_par0[0],err_hesse_par0[1]);
  TH1D *hErrHESSEsym1 = new TH1D("hErrHESSEsym1","hErrHESSEsym1",1000,err_hesse_par1[0],err_hesse_par1[1]);
  TH1D *hnchisqcalls = new TH1D("hnchisqcalls","hnchisqcalls",100,0,700);
  TH1D *hfitstatus = new TH1D("hfitstatus","hfitstatus",5,0,5);
  TH1D *hpullsHESSE0 = new TH1D("hpullsHESSE0","hpullsHESSE0",100,-5,5);
  TH1D *hpullsHESSE1 = new TH1D("hpullsHESSE1","hpullsHESSE1",100,-5,5);
  TH1D *hpullsMINOS0 = new TH1D("hpullsMINOS0","hpullsMINOS0",100,-5,5);
  TH1D *hpullsMINOS1 = new TH1D("hpullsMINOS1","hpullsMINOS1",100,-5,5);
  TH1D *hresidual0 = new TH1D("hresidual0","hresidual0",100,-5,5);
  TH1D *hresidual1 = new TH1D("hresidual1","hresidual1",100,-5,5);

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
  int Nin_90CL = 0;
  int Nin_95CL = 0;
  int Nin_2sigma = 0;
  
  int NdiffChisq = 0;
  
  int Nentries = treein->GetEntries();
  
  cout << "Nentries = " << Nentries << endl;

  //int NinPROVA = 0;


  
  for(int ientry=0;ientry<Nentries;ientry++){

    if(CheckSingleToy){ 
      Nentries = 1;
      ientry = ToySel;
    }

    treein->GetEntry(ientry);
 
    // Find bestfit chisq and true chisq
    
    //cout << *fparm << endl;
    
    int bin_bestfit0 = hdeltachisq1D_out->FindBin(bestfit0);
    double deltachisq_bestfit0 = hdeltachisq1D_out->GetBinContent(bin_bestfit0);

    int bin_truepar0 = hdeltachisq1D_out->FindBin(truepar0);
    double deltachisq_truepar0 = hdeltachisq1D_out->GetBinContent(bin_truepar0);        

    double deltachisq0BF_used = conv::kBadResult;
    if(!UseFullChisqFit) deltachisq0BF_used = chisqtrue0-chisqmin;
    else                 deltachisq0BF_used = chisqtrue0-chisqminFullChisq;
    
    // Print chisq
    int npt_chisq_x = hchisq1D_out->GetNbinsX();
    for(int ix=0;ix<npt_chisq_x;ix++){
      cout << hchisq1D_out->GetBinContent(ix+1) << " ";
    }
    cout << endl;
    cout << "truepar0 = " << truepar0 << endl;
    cout << "bestfit0 = " << bestfit0 << endl;
    cout << "deltachisq (fit) = " << chisqtrue0-chisqmin << endl;
    cout << "deltachisq (global) = " << chisqtrue0-chisqminFullChisq << endl;
    cout << "deltachisq0BF (tree) = " << deltachisq0BF << endl;
    cout << "chisq0 - chisq1 = " << chisqtrue0-chisqtrue1 << endl;
    cout << "deltachisq01 (tree) = " << deltachisq01 << endl;
    cout << endl;

    // Fill deltaChisq histogram
    
    hdeltachisqDiff->Fill(deltachisq0BF_used);

    myvector.push_back(deltachisq0BF_used);

    // Minimized chisq distribution
    hchisqmin->Fill(chisqmin);
    hbestfit_2D->Fill(bestfit0,bestfit1);
    hbestfit_1D->Fill(bestfit0);

    // Fill histograms for fit status
    hfitstatus->Fill(fitstatus);
    hnchisqcalls->Fill(nchisqcalls);
    
    // Full chisq distribution 
    hchisqminFullChisq->Fill(chisqminFullChisq);
    hbestfitFullChisq_2D->Fill(bestfitFullChisq0,bestfitFullChisq1);
    hbestfitFullChisq_1D->Fill(bestfitFullChisq0);

    // Fill variables for coverage study
    if(deltachisq0BF_used <= GetConstChisq(IsFit2D,"1sigma") ) Nin_1sigma++;
    if(deltachisq0BF_used <= GetConstChisq(IsFit2D,"90CL")   ) Nin_90CL++;
    if(deltachisq0BF_used <= GetConstChisq(IsFit2D,"95CL")   ) Nin_95CL++;
    if(deltachisq0BF_used <= GetConstChisq(IsFit2D,"2sigma") ) Nin_2sigma++;

    //if(deltachisq0BF_used <= 1.6212065) NinPROVA++;

    // Fill histograms for pull validation and errors
       
    // // cout << endl;
    // // cout << "deltaCP:" << endl;
    // // cout << "bestfit0 = " << bestfit0 << endl;
    // // cout << "truepar0 = " << truepar0 << endl;
    // // cout << "bestfiterrNeg0 = " << bestfiterrNeg0 << endl;
    // // cout << "bestfiterrPos0 = " << bestfiterrPos0 << endl;
 
    // bestfiterrNeg0 = TMath::Abs( TMath::Sin(bestfit0-bestfiterrNeg0) - TMath::Sin(bestfit0) );
    // //if( (TMath::Sin(bestfit0)-bestfiterrNeg0) < -1. ) bestfiterrNeg0 = TMath::Sin(bestfit0)+1.;
    
    // bestfiterrPos0 = TMath::Abs( TMath::Sin(bestfit0+bestfiterrPos0) - TMath::Sin(bestfit0) );
    // if( (TMath::Sin(bestfit0)+bestfiterrPos0) > +1. ) bestfiterrPos0 = 1.-TMath::Sin(bestfit0);
    
    // bestfit0 = TMath::Sin(bestfit0);
    // truepar0 = TMath::Sin(truepar0);
   
    // // cout << "sin(deltaCP):" << endl;
    // // cout << "bestfit0 = " << bestfit0 << endl;
    // // cout << "truepar0 = " << truepar0 << endl;
    // // cout << "bestfiterrNeg0 = " << bestfiterrNeg0 << endl;
    // // cout << "bestfiterrPos0 = " << bestfiterrPos0 << endl;
    
    hErrMINOSneg0->Fill(bestfiterrNeg0);
    hErrMINOSpos0->Fill(bestfiterrPos0);
    hErrMINOSneg1->Fill(bestfiterrNeg1);
    hErrMINOSpos1->Fill(bestfiterrPos1);
    hErrHESSEsym0->Fill(bestfiterrSym0);
    hErrHESSEsym1->Fill(bestfiterrSym1);

    double pullHESSE0 = -999;
    double pullHESSE1 = -999;
    pullHESSE0 = (truepar0 - bestfit0) / bestfiterrSym0;
    hpullsHESSE0->Fill(pullHESSE0);
    if(IsFit2D){
      pullHESSE1 = (truepar0 - bestfit0) / bestfiterrSym0;
      hpullsHESSE0->Fill(pullHESSE0);
    }

    double pullMINOS0 = -999;
    double pullMINOS1 = -999;
    if(bestfit0>truepar0) pullMINOS0 = (bestfit0 - truepar0) / (bestfiterrNeg0);
    else                  pullMINOS0 = (bestfit0 - truepar0) / (bestfiterrPos0);
    hpullsMINOS0->Fill(pullMINOS0);
    if(IsFit2D){
      if(bestfit1>truepar1) pullMINOS1 = (bestfit1 - truepar1) / (bestfiterrNeg1);
      else                  pullMINOS1 = (bestfit1 - truepar1) / (bestfiterrPos1);
    }

    double residual0 = -999;
    double residual1 = -999;
    residual0 = (bestfit0 - truepar0) / truepar0;
    hresidual0->Fill(residual0);
    if(IsFit2D){
      residual1 = (bestfit1 - truepar1) / truepar1;
      hresidual1->Fill(residual1);
    }

    // Check that full chisq min is equal to minimized best-fit
    if(chisqminFullChisq!=conv::kALargeNum){
      if( fabs(chisqmin-chisqminFullChisq)>1.E-6 ){
	double diff = chisqmin-chisqminFullChisq;
	hchisqminDiff->Fill(diff);
	cout << std::setprecision(8);
	cout << endl;
	cout << "Toy " << ientry << " --> ToyID = " << ToyID << endl;
	cout << "chisqmin = " << chisqmin
	     << " - bestfit0 = " << bestfit0
	     << " - bestfit1 = " << bestfit1 << endl;
	cout << "chisqminFullChisq = " << chisqminFullChisq 
	     << " - bestfitFullChisq0 = " << bestfitFullChisq0 
	     << " - bestfitFullChisq1 = " << bestfitFullChisq1 << endl;
	cout << "Difference = " << diff << endl;
	cout << "# of Chisq calls = " << nchisqcalls << endl;
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
  
  // Calculate critical value
  double FCCritVal_1sigma = CalcDLcrit(myvector,CLfrac_1sigma*100.,0.); 
  double FCCritVal_90CL   = CalcDLcrit(myvector,CLfrac_90CL*100.,0.); 
  double FCCritVal_95CL   = CalcDLcrit(myvector,CLfrac_95CL*100.,0.); 
  
  // Get pulls

  double FitMeanPar0_HESSE = -999;
  double FitRMSPar0_HESSE = -999;
  double FitMeanPar1_HESSE = -999;
  double FitRMSPar1_HESSE = -999;
  double FitMeanPar0_MINOS = -999;
  double FitRMSPar0_MINOS = -999;
  double FitMeanPar1_MINOS = -999;
  double FitRMSPar1_MINOS = -999;

  TF1 *fgaus = new TF1("fgaus","gaus",20,-20);
  hpullsHESSE0->Fit(fgaus,"Q");
  FitMeanPar0_HESSE = fgaus->GetParameter(1);
  FitRMSPar0_HESSE = fgaus->GetParameter(2); 
  if(IsFit2D){
    hpullsHESSE1->Fit(fgaus,"Q");
    FitMeanPar1_HESSE = fgaus->GetParameter(1);
    FitRMSPar1_HESSE = fgaus->GetParameter(2);
  }

  hpullsMINOS0->Fit(fgaus,"Q");
  FitMeanPar0_MINOS = fgaus->GetParameter(1);
  FitRMSPar0_MINOS = fgaus->GetParameter(2);
  if(IsFit2D){
    hpullsMINOS1->Fit(fgaus,"Q");
    FitMeanPar1_MINOS = fgaus->GetParameter(1);
    FitRMSPar1_MINOS = fgaus->GetParameter(2);
  }


  // Draw histograms

  TCanvas *cnchisqcalls = new TCanvas("cnchisqcalls","cnchisqcalls");
  hnchisqcalls->GetYaxis()->SetTitle("# entries");
  hnchisqcalls->GetXaxis()->SetTitle("# of #chi^{2} computations");
  hnchisqcalls->DrawCopy();

  TCanvas *cfitstatus = new TCanvas("cfitstatus","cfitstatus");
  hfitstatus->GetYaxis()->SetTitle("# entries");
  hfitstatus->GetXaxis()->SetTitle("Fit Status");
  hfitstatus->DrawCopy();

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

  TCanvas *cErrMINOSneg0 = new TCanvas("cErrMINOSneg0","cErrMINOSneg0");
  hErrMINOSneg0->GetYaxis()->SetTitle("# entries");
  
  if(IsFit2D){
    char name[300];
    sprintf(name,"#Delta %s (MINOS negative)",namepar0.c_str());
    hErrMINOSneg0->GetXaxis()->SetTitle(name);
  }
  else{
    char name[300];
    sprintf(name,"#Delta %s (MINOS negative)",namepar0.c_str());
    hErrMINOSneg0->GetXaxis()->SetTitle(name);
  }
  hErrMINOSneg0->DrawCopy();
  
  TCanvas *cErrMINOSpos0 = new TCanvas("cErrMINOSpos0","cErrMINOSpos0");
  hErrMINOSpos0->GetYaxis()->SetTitle("# entries");
  if(IsFit2D){
    char name[300];
    sprintf(name,"#Delta %s (MINOS positive)",namepar0.c_str());
    hErrMINOSpos0->GetXaxis()->SetTitle(name);
  }
  else{
    char name[300];
    sprintf(name,"#Delta %s (MINOS positive)",namepar0.c_str());
    hErrMINOSpos0->GetXaxis()->SetTitle(name);
  }
  hErrMINOSpos0->DrawCopy();
  
  TCanvas *cErrMINOSneg1;
  TCanvas *cErrMINOSpos1;
  if(IsFit2D){
    cErrMINOSneg1 = new TCanvas("cErrMINOSneg1","cErrMINOSneg1");
    hErrMINOSneg1->GetYaxis()->SetTitle("# entries");
    char name[300];
    sprintf(name,"#Delta %s (MINOS negative)",namepar1.c_str());
    hErrMINOSneg1->GetXaxis()->SetTitle(name);
    hErrMINOSneg1->DrawCopy();  
    cErrMINOSpos1 = new TCanvas("cErrMINOSpos1","cErrMINOSpos1");
    hErrMINOSpos1->GetYaxis()->SetTitle("# entries");
    char name[300];
    sprintf(name,"#Delta %s (MINOS positive)",namepar1.c_str());
    hErrMINOSpos1->GetXaxis()->SetTitle(name);
    hErrMINOSpos1->DrawCopy();
  }

  TCanvas *cErrHESSEsym0 = new TCanvas("cErrHESSEsym0","cErrHESSEsym0");
  hErrHESSEsym0->GetYaxis()->SetTitle("# entries");
  if(IsFit2D){
    char name[300];
    sprintf(name,"#Delta %s (HESSE)",namepar0.c_str());
    hErrHESSEsym0->GetXaxis()->SetTitle(name);
  }
  else{
    char name[300];
    sprintf(name,"#Delta %s (HESSE)",namepar0.c_str());
    hErrHESSEsym0->GetXaxis()->SetTitle(name);
  }
  hErrHESSEsym0->DrawCopy();

  TCanvas *cErrHESSEsym1;
  if(IsFit2D){  
    cErrHESSEsym1 = new TCanvas("cErrHESSEsym1","cErrHESSEsym1");
    hErrHESSEsym1->GetYaxis()->SetTitle("# entries");
    char name[300];
    sprintf(name,"#Delta %s (HESSE)",namepar1.c_str());
    hErrHESSEsym1->GetXaxis()->SetTitle(name);
    hErrHESSEsym1->DrawCopy();
  }
 
  TCanvas *cchisqminDiff = new TCanvas("hchisqminDiff","hchisqminDiff");
  hchisqminDiff->GetYaxis()->SetTitle("# entries");
  hchisqminDiff->GetXaxis()->SetTitle("#chi^{2}_{Descent} - #chi^{2}_{Global}");
  hchisqminDiff->DrawCopy();

  TCanvas *cdeltachisqDiff = new TCanvas("cdeltachisqDiff","cdeltachisqDiff");
  hdeltachisqDiff->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff->GetXaxis()->SetTitle("#Delta #chi^{2}_{Descent} (no global)");
  hdeltachisqDiff->DrawCopy();
  TLine *lcritval_90CL = new TLine(FCCritVal_90CL,0,FCCritVal_90CL,hdeltachisqDiff->GetMaximum());
  lcritval_90CL->SetLineColor(kRed);
  lcritval_90CL->Draw();

  TCanvas *cpullsHESSE0 = new TCanvas("cpullsHESSE0","cpullsHESSE0");
  TLatex *lpullsHESSE0;
  if(IsFit2D){
    char name[300];
    sprintf(name,"Pull %s (HESSE error)",namepar0.c_str());
    hpullsHESSE0->GetXaxis()->SetTitle(name);
  }
  else{
    char name[300];
    sprintf(name,"Pull %s (HESSE error)",namepar0.c_str());
    hpullsHESSE0->GetXaxis()->SetTitle(name);
   }
  hpullsHESSE0->GetYaxis()->SetTitle("# entries");
  hpullsHESSE0->DrawCopy();
  sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsHESSE0->GetMean(),hpullsHESSE0->GetRMS(),FitMeanPar0_HESSE,FitRMSPar0_HESSE); 
  lpullsHESSE0 = new TLatex(0.90,hpullsHESSE0->GetMaximum()*0.8,name);
  lpullsHESSE0->SetTextSize(0.07);
  lpullsHESSE0->Draw("same");

  TCanvas *cpullsHESSE1;
  TLatex *lpullsHESSE1;
  if(IsFit2D){
    cpullsHESSE1 = new TCanvas("cpullsHESSE1","cpullsHESSE1");
    char name[300];
    sprintf(name,"Pull %s (HESSE error)",namepar1.c_str());
    hpullsHESSE1->GetXaxis()->SetTitle(name);
    hpullsHESSE1->GetYaxis()->SetTitle("# entries");
    hpullsHESSE1->DrawCopy();
    sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsHESSE1->GetMean(),hpullsHESSE1->GetRMS(),FitMeanPar1_HESSE,FitRMSPar1_HESSE); 
    lpullsHESSE1 = new TLatex(0.90,hpullsHESSE1->GetMaximum()*0.8,name);
    lpullsHESSE1->SetTextSize(0.07);
    lpullsHESSE1->Draw();
  }

  TCanvas *cpullsMINOS0 = new TCanvas("cpullsMINOS0","cpullsMINOS0");
  if(IsFit2D){
    char name[300];
    sprintf(name,"Pull %s (MINOS error)",namepar0.c_str());
    hpullsMINOS0->GetXaxis()->SetTitle(name);
  }
  else{
    char name[300];
    sprintf(name,"Pull %s (MINOS error)",namepar0.c_str());
    hpullsMINOS0->GetXaxis()->SetTitle(name);
  }
  hpullsMINOS0->GetYaxis()->SetTitle("# entries");
  hpullsMINOS0->DrawCopy();
  sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsMINOS0->GetMean(),hpullsMINOS0->GetRMS(),FitMeanPar0_MINOS,FitRMSPar0_MINOS); 
  TLatex *lpullsMINOS0 = new TLatex(0.90,hpullsMINOS0->GetMaximum()*0.8,name);
  lpullsMINOS0->SetTextSize(0.07);
  lpullsMINOS0->Draw();

  TCanvas *cpullsMINOS1;
  TLatex *lpullsMINOS1;
  if(IsFit2D){
    cpullsMINOS1 = new TCanvas("cpullsMINOS1","cpullsMINOS1");
    char name[300];
    sprintf(name,"Pull %s (MINOS error)",namepar1.c_str());
    hpullsMINOS1->GetXaxis()->SetTitle(name);
    hpullsMINOS1->GetYaxis()->SetTitle("# entries");
    hpullsMINOS1->DrawCopy();
    sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsMINOS1->GetMean(),hpullsMINOS1->GetRMS(),FitMeanPar1_MINOS,FitRMSPar1_MINOS); 
    lpullsMINOS1 = new TLatex(0.90,hpullsMINOS1->GetMaximum()*0.8,name);
    lpullsMINOS1->SetTextSize(0.07);
    lpullsMINOS1->Draw();
  }

  TCanvas *cresidual0 = new TCanvas("cresidual0","cresidual0");

  if(IsFit2D){
    char name[300];
    sprintf(name,"Residual %s",namepar0.c_str());
    hresidual0->GetXaxis()->SetTitle(name);
  }
  else{
    char name[300];
    sprintf(name,"Residual %s",namepar0.c_str());
    hresidual0->GetXaxis()->SetTitle(name); 
  }
  hresidual0->GetYaxis()->SetTitle("# entries");
  hresidual0->DrawCopy();
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hresidual0->GetMean(),hresidual0->GetRMS()); 
  TLatex *lresidual0 = new TLatex(0.90,hresidual0->GetMaximum()*0.8,name);
  lresidual0->SetTextSize(0.07);
  lresidual0->Draw();

  TCanvas *cresidual1;
  TLatex *lresidual1;
  if(IsFit2D){
    cresidual1 = new TCanvas("cresidual1","cresidual1");
    char name[300];
    sprintf(name,"Residual %s",namepar1.c_str());
    hresidual1->GetXaxis()->SetTitle(name);
    hresidual1->GetYaxis()->SetTitle("# entries");
    hresidual1->DrawCopy();
    sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hresidual1->GetMean(),hresidual1->GetRMS()); 
    TLatex *lresidual1 = new TLatex(0.90,hresidual1->GetMaximum()*0.8,name);
    lresidual1->SetTextSize(0.07);
    lresidual1->Draw();
  }


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


  // Print results 
  
  cout << endl;
  cout << "truepar0 = " << truepar0 << endl;
  if(IsFit2D) cout << "truepar1 = " << truepar1 << endl;
  cout << endl;
  cout << "Difference between descent gradient and global minimum" << endl;
  cout << "# of toys with different minim and full chisq minimum: "
       << NdiffChisq << " / " << Nentries << " toys" << endl;
  cout << "Mean difference = " << hchisqminDiff->GetMean() << endl;
  cout << "RMS difference = " << hchisqminDiff->GetRMS() << endl;
  cout << endl;
  cout << "Informations on the fit errors" << endl;
  cout << "Param0: " << endl;
  cout << "  - average HESSE symmetric error = " 
       << hErrHESSEsym0->GetMean() 
       << " +- " << hErrHESSEsym0->GetRMS() << endl;
  cout << "  - average MINOS positive error = " 
       << hErrMINOSpos0->GetMean() 
       << " +- " << hErrMINOSpos0->GetRMS() << endl;
  cout << "  - average MINOS negative error = " 
       << hErrMINOSneg0->GetMean() 
       << " +- " << hErrMINOSneg0->GetRMS() << endl;
  if(IsFit2D){
    cout << "Param1: " << endl;
    cout << "  - average HESSE symmetric error = " 
	 << hErrHESSEsym1->GetMean() 
	 << " +- " << hErrHESSEsym1->GetRMS() << endl;
    cout << "  - average MINOS positive error = " 
	 << hErrMINOSpos1->GetMean() 
	 << " +- " << hErrMINOSpos1->GetRMS() << endl;
    cout << "  - average MINOS negative error = " 
	 << hErrMINOSneg1->GetMean() 
	 << " +- " << hErrMINOSneg1->GetRMS() << endl;
    cout << endl;
  }

  cout << endl;
  cout << "Pulls w/ HESSE errors:" << endl;
  cout << " - Mean Par0 = " << hpullsHESSE0->GetMean() << endl;
  cout << " - RMS Par0 = " << hpullsHESSE0->GetRMS() << endl;
  cout << " - Fit Mean Par0 = " << FitMeanPar0_HESSE << endl;
  cout << " - Fit RMS Par0 = " << FitRMSPar0_HESSE << endl;
  if(IsFit2D){
    cout << " - Mean Par1 = " << hpullsHESSE1->GetMean() << endl;
    cout << " - RMS Par1 = " << hpullsHESSE1->GetRMS() << endl;
    cout << " - Fit Mean Par1 = " << FitMeanPar1_HESSE << endl;
    cout << " - Fit RMS Par1 = " << FitRMSPar1_HESSE << endl;
  }
  cout << endl;
  cout << "Pulls w/ MINOS errors:" << endl;
  cout << " - Mean Par0 = " << hpullsMINOS0->GetMean() << endl;
  cout << " - RMS Par0 = " << hpullsMINOS0->GetRMS() << endl;
  cout << " - Fit Mean Par0 = " << FitMeanPar0_MINOS << endl;
  cout << " - Fit RMS Par0 = " << FitRMSPar0_MINOS << endl;
  if(IsFit2D){
    cout << " - Mean Par1 = " << hpullsMINOS1->GetMean() << endl;
    cout << " - RMS Par1 = " << hpullsMINOS1->GetRMS() << endl;
    cout << " - Fit Mean Par1 = " << FitMeanPar1_MINOS << endl;
    cout << " - Fit RMS Par1 = " << FitRMSPar1_MINOS << endl;
  }
  cout << endl;
  cout << "Coverage:" << endl;
  cout << "Fraction of toys inside 1sigma = " << FracIn_1sigma << endl;
  cout << "Fraction of toys inside 90CL   = " << FracIn_90CL   << endl;
  cout << "Fraction of toys inside 95CL   = " << FracIn_95CL   << endl;
  cout << "Fraction of toys inside 2sigma = " << FracIn_2sigma << endl;
  cout << endl;
  cout << "Feldman-Cousins critical values:" << endl;
  cout << "Crit Value 1sigma = " << FCCritVal_1sigma << endl; 
  cout << "Crit Value 90CL = " << FCCritVal_90CL << endl; 
  cout << "Crit Value 95CL = " << FCCritVal_95CL << endl; 
  cout << endl;

  //cout << "NinPROVA = " << NinPROVA << endl;
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
