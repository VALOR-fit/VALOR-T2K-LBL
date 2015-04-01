
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

void analyze_run_fit_marginalization_readbinary
(
 //char infilename[200] = "fitmarg_sin2th23-3f+3fbar_10margtoys_3500toys.root",         // 1D (sin2th23)
 //char infilename[200] = "fitmarg_sin2th23Vsdmsq-3f+3fbar_10margtoys_3500toys.root", // 2D (sin2th23 Vs dmsq)
 //char infilename[200] = "fitmarg_sin2th23Vsdmsq_MINOSHESS-3f+3fbar_10margtoys_3500toys.root", // 2D (sin2th23 Vs dmsq)

 char infilename[200] = "fitmarg_sin2th13_MINOSHESSE-3f+3fbar_10margtoys_3500toys.root", // 1D (sin2th13)
 )
{
  TFile *filein = new TFile(infilename,"READ");
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
  TH2F *hchisq2D_out;
  TH1F *hchisq1D_out;
  TH2F *hstep2D_out;
  TH1F *hstep1D_out;
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

  treein->SetBranchAddress("best_fit_params",&fparm);
  treein->SetBranchAddress("input_params",&iparm);
  treein->SetBranchAddress("truepar0", &truepar0 );
  treein->SetBranchAddress("truepar1", &truepar1 );
  treein->SetBranchAddress("bestfit0", &bestfit0 );
  treein->SetBranchAddress("bestfit1", &bestfit1 );
  treein->SetBranchAddress("hchisq2D", &hchisq2D_out );
  treein->SetBranchAddress("hchisq1D", &hchisq1D_out );
  treein->SetBranchAddress("hstep2D", &hstep2D_out );
  treein->SetBranchAddress("hstep1D", &hstep1D_out );
  treein->SetBranchAddress("fitstatus",&fitstatus);
  treein->SetBranchAddress("chisqmin",&chisqmin);
  treein->SetBranchAddress("chisqtrue",&chisqtrue);
  treein->SetBranchAddress("nfits",&nfits);
  treein->SetBranchAddress("nchisqcalls",&nchisqcalls);
  treein->SetBranchAddress("niterations",&niterations);
  treein->SetBranchAddress("chisqmin_subrange",&chisqmin_subrange);
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


  // Declare histograms

  TH1F *hbestfitFullChisq_1D = new TH1F("hbestfitFullChisq_1D","hbestfitFullChisq_1D",201,0,1);
  TH1F *hbestfit_1D = new TH1F("hbestfit_1D","hbestfit_1D",201,0,1);

  TH2F *hbestfitFullChisq_2D = new TH2F("hbestfitFullChisq_2D","hbestfitFullChisq_2D",201,0,1,201,0.0002,0.02);
  TH2F *hbestfit_2D = new TH2F("hbestfit_2D","hbestfit_2D",201,0,1,201,0.0002,0.02);

  TH1F *hdeltachisqDiff = new TH1F("hdeltachisqDiff","hdeltachisqDiff",100,0,20);
  TH1F *hchisqminDiff = new TH1F("hchisqminDiff","hchisqminDiff",100,0.,0.3);
  TH1F *hchisqmin = new TH1F("hchisqmin","hchisqmin",100,0,100);
  TH1F *hchisqminFullChisq = new TH1F("hchisqminFullChisq","hchisqminFullChisq",100,0,100);

  TH1F *hErrMINOSneg0 = new TH1F("hErrMINOSneg0","hErrMINOSneg0",100,0,1.);
  TH1F *hErrMINOSneg1 = new TH1F("hErrMINOSneg1","hErrMINOSneg1",100,0,0.01);
  TH1F *hErrMINOSpos0 = new TH1F("hErrMINOSpos0","hErrMINOSpos0",100,0,1.);
  TH1F *hErrMINOSpos1 = new TH1F("hErrMINOSpos1","hErrMINOSpos1",100,0,0.01);
  TH1F *hErrHESSEsym0 = new TH1F("hErrHESSEsym0","hErrHESSEsym0",100,0,1.);
  TH1F *hErrHESSEsym1 = new TH1F("hErrHESSEsym1","hErrHESSEsym1",100,0,0.01);

  TH1F *hnchisqcalls = new TH1F("hnchisqcalls","hnchisqcalls",100,0,700);
  TH1F *hfitstatus = new TH1F("hfitstatus","hfitstatus",5,0,5);
  int NdiffChisq = 0;


  // Check if 1D or 2D fit
  bool IsFit2D;
  treein->GetEntry(0);
  if(bestfit1==-999) IsFit2D = false;
  else               IsFit2D = true;

  // loop over the tree  
  int Nentries = treein->GetEntries();
  for(int ientry=0;ientry<Nentries;ientry++){

    treein->GetEntry(ientry);

    hfitstatus->Fill(fitstatus);

    //if(nchisqcalls>700){
    //cout << "# of chisq calls = " << nchisqcalls << endl;
    //exit(1);
    //}
    hnchisqcalls->Fill(nchisqcalls);

    if(bestfitFullChisq1>0.006) continue; // Only first minimum

    double diff = chisqmin-chisqminFullChisq;
    //if(diff<0.04) continue;  // select large chisq difference

    if(chisqmin!=chisqminFullChisq){
      cout << endl;
      cout << "Toy " << ientry << endl;
      cout << "chisqmin = " << chisqmin
    	   << " - bestfit0 = " << bestfit0
    	   << " - bestfit1 = " << bestfit1 << endl;
      cout << "chisqminFullChisq = " << chisqminFullChisq 
    	   << " - bestfitFullChisq0 = " << bestfitFullChisq0 
    	   << " - bestfitFullChisq1 = " << bestfitFullChisq1 << endl;
      cout << "Difference = " << diff << endl;
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

    double diff_deltachisq = chisqtrue - chisqmin;
    if(diff_deltachisq>20){
      cout << "chisqtrue - chisqmin = " << diff_deltachisq << endl;
      exit(1);
    }
    hdeltachisqDiff->Fill(diff_deltachisq);
    
    hchisqminDiff->Fill(diff);

    // Minimized chisq distribution
    hchisqmin->Fill(chisqmin);
    hbestfit_2D->Fill(bestfit0,bestfit1);
    hbestfit_1D->Fill(bestfit0);

    if(bestfit0>0.6){
      cout << endl;
      cout << "Toy " << ientry << endl;
      cout << "bestfit0 = " << bestfit0 
	   << " + " << bestfiterrPos0 
	   << " - " << bestfiterrNeg0 
	   << endl;
      cout << "bestfit1 = " << bestfit1 
	   << " + " << bestfiterrPos1 
	   << " - " << bestfiterrNeg1 
	   << endl;
      cout << "chisqmin = " << chisqmin << endl;
      cout << endl;
    }
      
    // Full chisq distribution
    hchisqminFullChisq->Fill(chisqminFullChisq);
    hbestfitFullChisq_2D->Fill(bestfitFullChisq0,bestfitFullChisq1);
    hbestfitFullChisq_1D->Fill(bestfitFullChisq0);
  }
  

  // Print informations
  
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


  TCanvas *cnchisqcalls = new TCanvas("cnchisqcalls","cnchisqcalls");
  hnchisqcalls->GetYaxis()->SetTitle("# entries");
  hnchisqcalls->GetXaxis()->SetTitle("# of #chi^{2} computations");
  hnchisqcalls->DrawCopy();

  TCanvas *cfitstatus = new TCanvas("cfitstatus","cfitstatus");
  hfitstatus->GetYaxis()->SetTitle("# entries");
  hfitstatus->GetXaxis()->SetTitle("Fit Status");
  hfitstatus->DrawCopy();

  TCanvas *cbestfit = new TCanvas("hbestfit","hbestfit");
  if(IsFit2D) hbestfit_2D->DrawCopy("colz");
  else        hbestfit_1D->DrawCopy();

  TCanvas *cbestfitFullChisq = new TCanvas("hbestfitFullChisq","hbestfitFullChisq");
  if(IsFit2D) hbestfitFullChisq_2D->DrawCopy("colz");
  else        hbestfitFullChisq_1D->DrawCopy();

  TCanvas *cErrMINOSneg0 = new TCanvas("cErrMINOSneg0","cErrMINOSneg0");
  hErrMINOSneg0->GetYaxis()->SetTitle("# entries");
  hErrMINOSneg0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{23} (MINOS negative)");
  hErrMINOSneg0->DrawCopy();
  TCanvas *cErrMINOSpos0 = new TCanvas("cErrMINOSpos0","cErrMINOSpos0");
  hErrMINOSpos0->GetYaxis()->SetTitle("# entries");
  hErrMINOSpos0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{23} (MINOS positive)");
  hErrMINOSpos0->DrawCopy();

  TCanvas *cErrMINOSneg1 = new TCanvas("cErrMINOSneg1","cErrMINOSneg1");
  hErrMINOSneg1->GetYaxis()->SetTitle("# entries");
  hErrMINOSneg1->GetXaxis()->SetTitle("#Delta #bar{m}^{2}_{32} (MINOS positive)");
  hErrMINOSneg1->DrawCopy();
  TCanvas *cErrMINOSpos1 = new TCanvas("cErrMINOSpos1","cErrMINOSpos1");
  hErrMINOSpos1->GetYaxis()->SetTitle("# entries");
  hErrMINOSpos1->GetXaxis()->SetTitle("#Delta #bar{m}^{2}_{32} (MINOS negative)");
  hErrMINOSpos1->DrawCopy();

  TCanvas *cErrHESSEsym0 = new TCanvas("cErrHESSEsym0","cErrHESSEsym0");
  hErrHESSEsym0->GetYaxis()->SetTitle("# entries");
  hErrHESSEsym0->GetXaxis()->SetTitle("#Delta sin^{2} #bar{#theta}_{23} (HESSE)");
  hErrHESSEsym0->DrawCopy();
  TCanvas *cErrHESSEsym1 = new TCanvas("cErrHESSEsym1","cErrHESSEsym1");
  hErrHESSEsym1->GetYaxis()->SetTitle("# entries");
  hErrHESSEsym1->GetXaxis()->SetTitle("#Delta #bar{m}^{2}_{32} (HESSE)");
  hErrHESSEsym1->DrawCopy();

  TCanvas *cchisqminDiff = new TCanvas("hchisqminDiff","hchisqminDiff");
  hchisqminDiff->GetYaxis()->SetTitle("# entries");
  hchisqminDiff->GetXaxis()->SetTitle("#chi^{2}_{Descent} - #chi^{2}_{Global}");
  hchisqminDiff->DrawCopy();
  TCanvas *cdeltachisqDiff = new TCanvas("cdeltachisqDiff","cdeltachisqDiff");
  hdeltachisqDiff->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff->GetXaxis()->SetTitle("#Delta #chi^{2}_{Descent} (no global)");
  hdeltachisqDiff->DrawCopy();

  //treein->GetEntry(1707);
  //treein->GetEntry(1824);
  treein->GetEntry(3429);
  TCanvas *cchisq_out_GetEntry = new TCanvas("cchisq_out_GetEntry","cchisq_out_GetEntry");
  if(IsFit2D) hchisq2D_out->DrawCopy("colz");
  else        hchisq1D_out->DrawCopy();

  return;
}
