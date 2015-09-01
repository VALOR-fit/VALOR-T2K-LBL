
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;

void analyze_run_fit_marginalization_readbinary_joint_2d
( 
 double trueval0 = 0.0243,
 double trueval1 = -1.609,
 //double trueval0 = 0.4,
 //double trueval1 = 2.509E-3,
 

 // NEW TOYS

 // sinsqth13 Vs dcp (Run 1-4)
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Kreact_Run1to4_Data_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Kreact_Run1to4_Data_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",

 // sinsqth13 Vs dcp (Run 1-4), no reactors
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Konly_Run1to4_Data_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Konly_Run1to4_Data_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 
 // sinsqth13 Vs dcp (Run 1-6)
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Kreact_Run1to6_DataSetA-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Kreact_Run1to6_DataSetB-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",

 char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Konly_Asimov_C_NH-sinsqtheta[12_0.304,13_0.042,23_0.528],dcp_1.910,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Konly_Asimov_C_IH-sinsqtheta[12_0.304,13_0.042,23_0.528],dcp_1.910,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 
 // sinsqth23 Vs dmsq (Run 1-6)
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth23Vsdmsq32_T2Kreact_Run1to6_DataSetA-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char inname[200]="../marginalization/joint_newtoys/fitmarg_sinsqth23Vsdmsq32_T2Kreact_Run1to6_DataSetD-sinsqtheta[12_0.304,13_0.024,23_0.400],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",

 const bool IsFit2D = true
  ) 
{
  char name[200];

  TFile *filein = new TFile(inname,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  
  PhysParamList *iparm;
  FitParamList *fparm; 
  vector<double> *chisqmin_subrange;
  vector<int> *niterations_subrange;
  vector<int> *nfits_subrange;
  vector<int> *fitstatus_subrange;
  vector<double> *bestfit0_subrange;
  vector<double> *bestfit1_subrange;  
  double truepar0=-999;
  double truepar1=-999;    
  int nfits=-999;
  int nchisqcalls=-999;
  int niterations=-999;
  int fitstatus=-999;
  double chisqmin=conv::kALargeNum;
  double chisqtrue=conv::kALargeNum;
  double chisqtrue0=conv::kALargeNum;
  double chisqtrue1=conv::kALargeNum;
  double deltachisqtrue=conv::kALargeNum;
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

  TBranch        *b_chisqmin_subrange;   //!

  treein->SetBranchAddress("best_fit_params",&fparm);
  treein->SetBranchAddress("input_params",&iparm);
  treein->SetBranchAddress("truepar0", &truepar0 );
  treein->SetBranchAddress("truepar1", &truepar1 );
  treein->SetBranchAddress("bestfit0", &bestfit0 );
  treein->SetBranchAddress("bestfit1", &bestfit1 );
  treein->SetBranchAddress("hchisq2D", &hchisq2D_out );
  treein->SetBranchAddress("hchisq1D", &hchisq1D_out );
  treein->SetBranchAddress("hpenaltyterm2D", &hpenaltyterm2D_out );
  treein->SetBranchAddress("hpenaltyterm1D", &hpenaltyterm1D_out );
  treein->SetBranchAddress("hdeltachisq2D", &hdeltachisq2D_out );
  treein->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D_out );
  treein->SetBranchAddress("hstep2D", &hstep2D_out );
  treein->SetBranchAddress("hstep1D", &hstep1D_out );
  treein->SetBranchAddress("fitstatus",&fitstatus);
  treein->SetBranchAddress("chisqmin",&chisqmin);
  treein->SetBranchAddress("chisqtrue",&chisqtrue);
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("deltachisqtrue",&deltachisqtrue);
  treein->SetBranchAddress("nfits",&nfits);
  treein->SetBranchAddress("nchisqcalls",&nchisqcalls);
  treein->SetBranchAddress("niterations",&niterations);
  treein->SetBranchAddress("chisqmin_subrange",&chisqmin_subrange,&b_chisqmin_subrange);
  treein->SetBranchAddress("niterations_subrange",&niterations_subrange);
  treein->SetBranchAddress("nfits_subrange",&nfits_subrange);
  treein->SetBranchAddress("fitstatus_subrange",&fitstatus_subrange);
  treein->SetBranchAddress("bestfit0_subrange",&bestfit0_subrange);
  treein->SetBranchAddress("bestfit1_subrange",&bestfit1_subrange);
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


  // Declare histograms

  TH1D *hbestfitFullChisq_1D = new TH1D("hbestfitFullChisq_1D","hbestfitFullChisq_1D",201,0,1);
  TH1D *hbestfit_1D = new TH1D("hbestfit_1D","hbestfit_1D",201,0,1);

  TH2D *hbestfitFullChisq_2D = new TH2D("hbestfitFullChisq_2D","hbestfitFullChisq_2D",201,0,0.3,201,-3.15,3.15);

  TH2D *hbestfit_2D = new TH2D("hbestfit_2D","hbestfit_2D",201,0,0.3,201,-3.15,3.15);

  TH1D *hdeltachisqDiff = new TH1D("hdeltachisqDiff","hdeltachisqDiff",100,0,20);
  TH1D *hchisqminDiff = new TH1D("hchisqminDiff","hchisqminDiff",100,0.,0.3);
  TH1D *hchisqmin = new TH1D("hchisqmin","hchisqmin",100,0,100);
  TH1D *hchisqminFullChisq = new TH1D("hchisqminFullChisq","hchisqminFullChisq",100,0,100);

  TH1D *hErrMINOSneg0 = new TH1D("hErrMINOSneg0","hErrMINOSneg0",1000,0,0.3);
  TH1D *hErrMINOSneg1 = new TH1D("hErrMINOSneg1","hErrMINOSneg1",1000,0.,6.28);
  TH1D *hErrMINOSpos0 = new TH1D("hErrMINOSpos0","hErrMINOSpos0",1000,0,0.3);
  TH1D *hErrMINOSpos1 = new TH1D("hErrMINOSpos1","hErrMINOSpos1",1000,0.,6.28);
  TH1D *hErrHESSEsym0 = new TH1D("hErrHESSEsym0","hErrHESSEsym0",1000,0,0.3);
  TH1D *hErrHESSEsym1 = new TH1D("hErrHESSEsym1","hErrHESSEsym1",1000,0.,20);

  TH1D *hnchisqcalls = new TH1D("hnchisqcalls","hnchisqcalls",100,0,700);
  TH1D *hfitstatus = new TH1D("hfitstatus","hfitstatus",5,0,5);
  int NdiffChisq = 0;

  TH1D *hpullsHESSE0 = new TH1D("hpullsHESSE0","hpullsHESSE0",100,-5,5);
  TH1D *hpullsHESSE1 = new TH1D("hpullsHESSE1","hpullsHESSE1",100,-5,5);
  TH1D *hpullsMINOS0 = new TH1D("hpullsMINOS0","hpullsMINOS0",100,-5,5);
  TH1D *hpullsMINOS1 = new TH1D("hpullsMINOS1","hpullsMINOS1",100,-5,5);

  TH1D *hresidual0 = new TH1D("hresidual0","hresidual0",100,-5,5);
  TH1D *hresidual1 = new TH1D("hresidual1","hresidual1",100,-5,5);

  // Check if 1D or 2D fit
  // bool IsFit2D;
  // treein->GetEntry(0);
  // if(bestfit1==-999) IsFit2D = false;
  // else               IsFit2D = true;
  
  cout << endl;
  if(IsFit2D) cout << "Fit 2D" << endl;
  else cout << "Fit 1D" << endl;
  cout << endl;

  // loop over the tree  
  
  int Nentries = treein->GetEntries();
  for(int ientry=0;ientry<Nentries;ientry++){

    treein->GetEntry(ientry);

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
    
    hfitstatus->Fill(fitstatus);

    hnchisqcalls->Fill(nchisqcalls);

    double diff = chisqmin-chisqminFullChisq;
    hchisqminDiff->Fill(diff);
   
    double diff_deltachisq = 0.;
    if(chisqmin!=conv::kBadResult) diff_deltachisq = chisqtrue0 - chisqmin;
    else diff_deltachisq = chisqtrue0 - chisqminFullChisq;
    hdeltachisqDiff->Fill(diff_deltachisq);

    if(chisqmin!=chisqminFullChisq){
      cout << std::setprecision(8);
      cout << endl;
      cout << "Toy " << ientry << " --> ToyID = " << ToyID << endl;
      cout << "chisqmin = " << chisqmin
    	   << " - bestfit0 = " << bestfit0
    	   << " - bestfit1 = " << bestfit1 << endl;
      cout << "chisqminFullChisq = " << chisqminFullChisq 
    	   << " - bestfitFullChisq0 = " << bestfitFullChisq0 
    	   << " - bestfitFullChisq1 = " << bestfitFullChisq1 << endl;
      cout << "Chisq_min - Chisq_true = " << diff_deltachisq << endl;
      cout << "Difference between full and minimized chisq= " << diff << endl;
      cout << "# of Chisq calls = " << nchisqcalls << endl;
      cout << endl; 
      NdiffChisq++;
    }
    else{
      continue;
    }

    hErrMINOSneg0->Fill(bestfiterrNeg0);
    hErrMINOSpos0->Fill(bestfiterrPos0);
    hErrMINOSneg1->Fill(bestfiterrNeg1);
    hErrMINOSpos1->Fill(bestfiterrPos1);
    hErrHESSEsym0->Fill(bestfiterrSym0);
    hErrHESSEsym1->Fill(bestfiterrSym1);


    // Minimized chisq distribution
    hchisqmin->Fill(chisqmin);
    hbestfit_2D->Fill(bestfit0,bestfit1);
    hbestfit_1D->Fill(bestfit0);

    // Full chisq distribution
    hchisqminFullChisq->Fill(chisqminFullChisq);
    hbestfitFullChisq_2D->Fill(bestfitFullChisq0,bestfitFullChisq1);
    hbestfitFullChisq_1D->Fill(bestfitFullChisq0);
  
  } // end loop tree

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
    hbestfit_2D->GetYaxis()->SetTitle("#delta_{CP}");
    hbestfit_2D->GetXaxis()->SetTitle("sin^{2} #bar{#theta}_{13}");
    //hbestfit_2D->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
    //hbestfit_2D->GetXaxis()->SetTitle("sin^{2} #bar{#theta}_{23}");
    hbestfit_2D->DrawCopy("colz");
  }
  else{
    hbestfit_1D->GetYaxis()->SetTitle("# entries");
    hbestfit_1D->GetXaxis()->SetTitle("sin^{2} #bar{#theta}_{13}");
    //hbestfit_1D->GetXaxis()->SetTitle("sin^{2} #bar{#theta}_{23}");
    hbestfit_1D->DrawCopy();
  }

  TCanvas *cbestfitFullChisq = new TCanvas("hbestfitFullChisq","hbestfitFullChisq");
  if(IsFit2D){
    hbestfitFullChisq_2D->GetYaxis()->SetTitle("#delta_{CP}");
    hbestfitFullChisq_2D->GetXaxis()->SetTitle("sin^{2} #bar{#theta}_{13}");
    //hbestfitFullChisq_2D->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
    //hbestfitFullChisq_2D->GetXaxis()->SetTitle("sin^{2} #bar{#theta}_{23}");
    hbestfitFullChisq_2D->DrawCopy("colz");
  }
  else{
    hbestfitFullChisq_1D->GetYaxis()->SetTitle("# entries");
    //hbestfitFullChisq_1D->GetXaxis()->SetTitle("sin^{2} #bar{#theta}_{23}");
    hbestfitFullChisq_1D->GetXaxis()->SetTitle("sin^{2} #bar{#theta}_{13}");
    hbestfitFullChisq_1D->DrawCopy();
  }

  TCanvas *cErrMINOSneg0 = new TCanvas("cErrMINOSneg0","cErrMINOSneg0");
  hErrMINOSneg0->GetYaxis()->SetTitle("# entries");
  //hErrMINOSneg0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{23} (MINOS negative)");
  hErrMINOSneg0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{13} (MINOS negative)");
  hErrMINOSneg0->DrawCopy();
  TCanvas *cErrMINOSpos0 = new TCanvas("cErrMINOSpos0","cErrMINOSpos0");
  hErrMINOSpos0->GetYaxis()->SetTitle("# entries");
  //hErrMINOSpos0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{23} (MINOS positive)");
  hErrMINOSpos0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{13} (MINOS positive)");
  hErrMINOSpos0->DrawCopy();

  TCanvas *cErrMINOSneg1 = new TCanvas("cErrMINOSneg1","cErrMINOSneg1");
  hErrMINOSneg1->GetYaxis()->SetTitle("# entries");
  //hErrMINOSneg1->GetXaxis()->SetTitle("#Delta #bar{m}^{2}_{32} (MINOS positive)");
  hErrMINOSneg1->GetXaxis()->SetTitle("#delta_{CP} (MINOS positive)");
  hErrMINOSneg1->DrawCopy();
  TCanvas *cErrMINOSpos1 = new TCanvas("cErrMINOSpos1","cErrMINOSpos1");
  hErrMINOSpos1->GetYaxis()->SetTitle("# entries");
  //hErrMINOSpos1->GetXaxis()->SetTitle("#Delta #bar{m}^{2}_{32} (MINOS negative)");
  hErrMINOSpos1->GetXaxis()->SetTitle("#delta_{CP} (MINOS negative)");
  hErrMINOSpos1->DrawCopy();

  TCanvas *cErrHESSEsym0 = new TCanvas("cErrHESSEsym0","cErrHESSEsym0");
  hErrHESSEsym0->GetYaxis()->SetTitle("# entries");
  //hErrHESSEsym0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{23} (HESSE)");
  hErrHESSEsym0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{13} (HESSE)");
  hErrHESSEsym0->DrawCopy();
  TCanvas *cErrHESSEsym1 = new TCanvas("cErrHESSEsym1","cErrHESSEsym1");
  hErrHESSEsym1->GetYaxis()->SetTitle("# entries");
  //hErrHESSEsym1->GetXaxis()->SetTitle("#Delta #bar{m}^{2}_{32} (HESSE)");
  hErrHESSEsym1->GetXaxis()->SetTitle("#delta_{CP} (HESSE)");
  hErrHESSEsym1->DrawCopy();

  TCanvas *cchisqminDiff = new TCanvas("hchisqminDiff","hchisqminDiff");
  hchisqminDiff->GetYaxis()->SetTitle("# entries");
  hchisqminDiff->GetXaxis()->SetTitle("#chi^{2}_{Descent} - #chi^{2}_{Global}");
  hchisqminDiff->DrawCopy();
  TCanvas *cdeltachisqDiff = new TCanvas("cdeltachisqDiff","cdeltachisqDiff");
  hdeltachisqDiff->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff->GetXaxis()->SetTitle("#Delta #chi^{2}_{Descent} (no global)");
  hdeltachisqDiff->DrawCopy();

  
  TCanvas *cpullsHESSE0 = new TCanvas("cpullsHESSE0","cpullsHESSE0");
  //hpullsHESSE0->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{23} (HESSE error)");
  hpullsHESSE0->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{13} (HESSE error)");
  hpullsHESSE0->GetYaxis()->SetTitle("# entries");
  hpullsHESSE0->DrawCopy();
  sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsHESSE0->GetMean(),hpullsHESSE0->GetRMS(),FitMeanPar0_HESSE,FitRMSPar0_HESSE); 
  TLatex *lpullsHESSE0 = new TLatex(0.90,hpullsHESSE0->GetMaximum()*0.8,name);
  lpullsHESSE0->SetTextSize(0.07);
  lpullsHESSE0->Draw("same");
  
  if(IsFit2D){
    TCanvas *cpullsHESSE1 = new TCanvas("cpullsHESSE1","cpullsHESSE1");
    //hpullsHESSE1->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{23} (HESSE error)");
    hpullsHESSE1->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{13} (HESSE error)");
    hpullsHESSE1->GetYaxis()->SetTitle("# entries");
    hpullsHESSE1->DrawCopy();
    sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsHESSE1->GetMean(),hpullsHESSE1->GetRMS(),FitMeanPar1_HESSE,FitRMSPar1_HESSE); 
    TLatex *lpullsHESSE1 = new TLatex(0.90,hpullsHESSE1->GetMaximum()*0.8,name);
    lpullsHESSE1->SetTextSize(0.07);
    lpullsHESSE1->Draw();
  }

  TCanvas *cpullsMINOS0 = new TCanvas("cpullsMINOS0","cpullsMINOS0");
  //hpullsMINOS0->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{23} (MINOS error)");
  hpullsMINOS0->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{13} (MINOS error)");
  hpullsMINOS0->GetYaxis()->SetTitle("# entries");
  hpullsMINOS0->DrawCopy();
  sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsMINOS0->GetMean(),hpullsMINOS0->GetRMS(),FitMeanPar0_MINOS,FitRMSPar0_MINOS); 
  TLatex *lpullsMINOS0 = new TLatex(0.90,hpullsMINOS0->GetMaximum()*0.8,name);
  lpullsMINOS0->SetTextSize(0.07);
  lpullsMINOS0->Draw();

  if(IsFit2D){
    TCanvas *cpullsMINOS1 = new TCanvas("cpullsMINOS1","cpullsMINOS1");
    //hpullsMINOS1->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{23} (MINOS error)");
    hpullsMINOS1->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{13} (MINOS error)");
    hpullsMINOS1->GetYaxis()->SetTitle("# entries");
    hpullsMINOS1->DrawCopy();
    sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsMINOS1->GetMean(),hpullsMINOS1->GetRMS(),FitMeanPar1_MINOS,FitRMSPar1_MINOS); 
    TLatex *lpullsMINOS1 = new TLatex(0.90,hpullsMINOS1->GetMaximum()*0.8,name);
    lpullsMINOS1->SetTextSize(0.07);
    lpullsMINOS1->Draw();
  }

  TCanvas *cresidual0 = new TCanvas("cresidual0","cresidual0");
  hresidual0->GetXaxis()->SetTitle("Residual sin^{2} #bar{#theta}_{13}");
  hresidual0->GetYaxis()->SetTitle("# entries");
  hresidual0->DrawCopy();
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hresidual0->GetMean(),hresidual0->GetRMS()); 
  TLatex *lresidual0 = new TLatex(0.90,hresidual0->GetMaximum()*0.8,name);
  lresidual0->SetTextSize(0.07);
  lresidual0->Draw();

  if(IsFit2D){
    TCanvas *cresidual1 = new TCanvas("cresidual1","cresidual1");
    hresidual1->GetXaxis()->SetTitle("Residual sin^{2} #bar{#theta}_{13}");
    hresidual1->GetYaxis()->SetTitle("# entries");
    hresidual1->DrawCopy();
    sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hresidual1->GetMean(),hresidual1->GetRMS()); 
    TLatex *lresidual1 = new TLatex(0.90,hresidual1->GetMaximum()*0.8,name);
    lresidual1->SetTextSize(0.07);
    lresidual1->Draw();
  }


  // Draw chisq distributions

  treein->GetEntry(0); // Get informations from the tree

  TCanvas *cchisq_out_GetEntry = new TCanvas("cchisq_out_GetEntry","cchisq_out_GetEntry");
  if(IsFit2D) hchisq2D_out->DrawCopy("colz");
  else        hchisq1D_out->DrawCopy();

  TCanvas *cdeltachisq_out_GetEntry = new TCanvas("cdeltachisq_out_GetEntry","cdeltachisq_out_GetEntry");
  if(IsFit2D) hdeltachisq2D_out->DrawCopy("colz");
  else        hdeltachisq1D_out->DrawCopy();

  TCanvas *cpenaltyterm_out_GetEntry = new TCanvas("cpenaltyterm_out_GetEntry","cpenaltyterm_out_GetEntry");
  if(IsFit2D) hpenaltyterm2D_out->DrawCopy("colz");
  else        hpenaltyterm1D_out->DrawCopy();


  // Draw contours

  TMarker * mPtBestFit = new TMarker(bestfitFullChisq0,bestfitFullChisq1,20);
  mPtBestFit->SetMarkerSize(1.5);
  mPtBestFit->SetMarkerColor(kBlue);

  //if(truepar0!=conv::kBadResult && truepar1!=conv::kBadResult){
  //trueval0 = truepar0;
  // trueval1 = truepar1;
  //}

  TMarker * mPtTrue = new TMarker(trueval0,trueval1,29);
  mPtTrue->SetMarkerSize(1.5);
  mPtTrue->SetMarkerColor(kBlack);

  TCanvas *ccontours = new TCanvas("ccontours","ccontours");
  TH2D *hcontour2D_68CL = (TH2D*)hdeltachisq2D_out->Clone("hcontour2D");
  hcontour2D_68CL->SetLineColor(kBlue);
  hcontour2D_68CL->SetLineWidth(2);
  hcontour2D_68CL->SetContour(1);
  hcontour2D_68CL->SetContourLevel(0,2.30);
  hcontour2D_68CL->Draw("cont3");
  TH2D *hcontour2D_90CL = (TH2D*)hdeltachisq2D_out->Clone("hcontour2D");
  hcontour2D_90CL->SetLineColor(kRed);
  hcontour2D_90CL->SetLineWidth(2);
  hcontour2D_90CL->SetContour(1);
  hcontour2D_90CL->SetContourLevel(0,4.61);
  hcontour2D_90CL->Draw("cont3 same");
  
  mPtBestFit->Draw("same");
  mPtTrue->Draw("same");

  TLegend *lcontours = new TLegend(0.5,0.6,0.8,0.8);
  lcontours->SetFillColor(0);
  lcontours->SetLineColor(0);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.05);
  lcontours->AddEntry(hcontour2D_68CL,"68% CL","l");
  lcontours->AddEntry(hcontour2D_90CL,"90% CL","l");
  lcontours->AddEntry(mPtBestFit,"VALOR best-fit","p");
  //lcontours->AddEntry(mPtTrue,"True value","p");
  lcontours->AddEntry(mPtTrue,"Old best-fit","p");
  lcontours->Draw();


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

  // cout << endl;
  // cout << std::setprecision(5);
  // cout << "Chisq grid:" << endl;
  // for(int ix=0;ix<hchisq2D_out->GetXaxis()->GetNbins();ix++){
  //   for(int iy=0;iy<hchisq2D_out->GetYaxis()->GetNbins();iy++){
  //     cout << hchisq2D_out->GetBinContent(ix+1,iy+1) << " ";
  //   }
  //   cout << endl;
  // }

  // double trueval0_oct = 0.62;
  // double trueval1_oct = 2.509E-3;
  // int bin_octant = hdeltachisq2D_out->FindBin(trueval0_oct,trueval1_oct);
  // //int bin_octant = hchisq2D_out->FindBin(trueval0_oct,trueval1_oct);  
  // double chisq_octant = hdeltachisq2D_out->GetBinContent(bin_octant);
  // //double chisq_octant = hchisq2D_out->GetBinContent(bin_octant);
  // cout << endl;
  // cout << "Check other octant: " << endl;
  // cout << " - sinsqth23 = " << trueval0_oct << endl;
  // cout << " - dmsq32 = " << trueval1_oct << endl;
  // cout << " chisq = " << chisq_octant << endl;
  // cout << endl;


  //// Save the contours  
  //TFile* f2 = new TFile(Form("%s_contour.root",filename),"recreate");
  //save_contour(f2,h1,h3);
  //h2->Write("chisq_map_2D");
  //bf->Write("best_fit_point_marker");
  //f2->Close();
  
    
  return;
}



void save_contour(TFile* f2, TH2D* h1=0, TH2D* h3=0)
{
    TCanvas* c2 = new TCanvas("c2","c",0,0,900,600);
    h3->Draw("CONT LIST");
    c2->Update();
    TObjArray *contours = 0;
    contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
    
    TList *list90 = (TList*)contours->At(0);
    TGraph *g;
    int it=0;
    TIter next(list90);
    while ((g = (TGraph*)next()))
    {
        char name[10];
        sprintf(name, "g90_valor_%d",it);
        g->SetLineWidth(2);
        g->SetLineStyle(1);
        g->SetLineColor(kRed);
        g->Write(name);
        it++;
    }
    
    c2->Clear();
    h1->Draw("CONT LIST");
    c2->Update();
    contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList *list68 = (TList*)contours->At(0);
    TIter next68(list68);
    it=0;
    while ((g = (TGraph*)next68()))
    {
        char name[10];
        sprintf(name, "g68_valor_%d",it);
        g->SetLineWidth(2);
        g->SetLineStyle(2);
        g->SetLineColor(kRed);
        g->Write(name);
        it++;
    }
    
    delete c2;
    delete contours;
    delete list68;
    delete list90;
    delete g;
    
}
