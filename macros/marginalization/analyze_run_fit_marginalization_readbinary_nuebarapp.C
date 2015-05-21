
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

void analyze_run_fit_marginalization_readbinary_nuebarapp
(
 //char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT_2-3f+normbar_10kToys.root",
 char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_4900Toys.root",
 )
{
  char name[200];

  TFile *filein = new TFile(infilename,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  
  vector<double> vec_chisq_01;
  vector<double> vec_chisq_0BF;
  vector<double> vec_evt_toys;

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
  double chisqtrue0=conv::kALargeNum;
  double chisqtrue1=conv::kALargeNum;
  double deltachisqtrue=conv::kALargeNum;
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
  int ToyID=-999;
  TH1D *hErecoToy;

  TBranch        *b_chisqmin_subrange;   //!

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
  treein->SetBranchAddress("hErecoToy",&hErecoToy);


  // Declare histograms

  TH1F *hbestfitFullChisq_1D = new TH1F("hbestfitFullChisq_1D","hbestfitFullChisq_1D",201,0,10);

  TH1F *hbestfit_1D = new TH1F("hbestfit_1D","hbestfit_1D",201,0-conv::kASmallNum,10);
  TH1F *hbestfit_1D_toyevt0 = new TH1F("hbestfit_1D_toyevt0","hbestfit_1D_toyevt0",201,0-conv::kASmallNum,10);
  TH1F *hbestfit_1D_toyevt1 = new TH1F("hbestfit_1D_toyevt1","hbestfit_1D_toyevt1",201,0-conv::kASmallNum,10);
  TH1F *hbestfit_1D_toyevt2 = new TH1F("hbestfit_1D_toyevt2","hbestfit_1D_toyevt2",201,0-conv::kASmallNum,10);
  THStack *sbestfit_1D = new THStack("sbestfit_1D","sbestfit_1D");
  
  TH1F *hdeltachisqDiff0BF = new TH1F("hdeltachisqDiff0BF","hdeltachisqDiff0BF",100,0,20);
  TH1F *hdeltachisqDiff0BF_toyevt0 = new TH1F("hdeltachisqDiff0BF_toyevt0","hdeltachisqDiff0BF_toyevt0",100,0,20);
  TH1F *hdeltachisqDiff0BF_toyevt1 = new TH1F("hdeltachisqDiff0BF_toyevt1","hdeltachisqDiff0BF_toyevt1",100,0,20);
  TH1F *hdeltachisqDiff0BF_toyevt2 = new TH1F("hdeltachisqDiff0BF_toyevt2","hdeltachisqDiff0BF_toyevt2",100,0,20);
  THStack *sdeltachisqDiff0BF = new THStack("sdeltachisqDiff0BF","sdeltachisqDiff0BF");

  TH1F *hdeltachisqDiff1BF = new TH1F("hdeltachisqDiff1BF","hdeltachisqDiff1BF",100,0,20);
  TH1F *hdeltachisqDiff1BF_toyevt0 = new TH1F("hdeltachisqDiff1BF_toyevt0","hdeltachisqDiff1BF_toyevt0",100,0,20);
  TH1F *hdeltachisqDiff1BF_toyevt1 = new TH1F("hdeltachisqDiff1BF_toyevt1","hdeltachisqDiff1BF_toyevt1",100,0,20);
  TH1F *hdeltachisqDiff1BF_toyevt2 = new TH1F("hdeltachisqDiff1BF_toyevt2","hdeltachisqDiff1BF_toyevt2",100,0,20);
  THStack *sdeltachisqDiff1BF = new THStack("sdeltachisqDiff1BF","sdeltachisqDiff1BF");

  TH1F *hdeltachisqDiff01 = new TH1F("hdeltachisqDiff01","hdeltachisqDiff01",200,-20,20);
  TH1F *hdeltachisqDiff01_toyevt0 = new TH1F("hdeltachisqDiff01_toyevt0","hdeltachisqDiff01_toyevt0",200,-20,20);
  TH1F *hdeltachisqDiff01_toyevt1 = new TH1F("hdeltachisqDiff01_toyevt1","hdeltachisqDiff01_toyevt1",200,-20,20);
  TH1F *hdeltachisqDiff01_toyevt2 = new TH1F("hdeltachisqDiff01_toyevt2","hdeltachisqDiff01_toyevt2",200,-20,20);
  THStack *sdeltachisqDiff01 = new THStack("sdeltachisqDiff01","sdeltachisqDiff01");

  TH1F *hchisqminDiff = new TH1F("hchisqminDiff","hchisqminDiff",100,0.,0.3); // difference between minim algorithm and global chisq

  TH1F *hchisqmin = new TH1F("hchisqmin","hchisqmin",100,0,100);
  TH1F *hchisqmin_toyevt0 = new TH1F("hchisqmin_toyevt0","hchisqmin_toyevt0",100,0,100);
  TH1F *hchisqmin_toyevt1 = new TH1F("hchisqmin_toyevt1","hchisqmin_toyevt1",100,0,100);
  TH1F *hchisqmin_toyevt2 = new TH1F("hchisqmin_toyevt2","hchisqmin_toyevt2",100,0,100);
  THStack *schisqmin = new THStack("schisqmin","schisqmin");

  TH1F *hchisq0 = new TH1F("hchisq0","hchisq0",100,0,100);
  TH1F *hchisq0_toyevt0 = new TH1F("hchisq0_toyevt0","hchisq0_toyevt0",100,0,100);
  TH1F *hchisq0_toyevt1 = new TH1F("hchisq0_toyevt1","hchisq0_toyevt1",100,0,100);
  TH1F *hchisq0_toyevt2 = new TH1F("hchisq0_toyevt2","hchisq0_toyevt2",100,0,100);
  THStack *schisq0 = new THStack("schisq0","schisq0");

  TH1F *hchisq1 = new TH1F("hchisq1","hchisq1",100,0,100);
  TH1F *hchisq1_toyevt0 = new TH1F("hchisq1_toyevt0","hchisq1_toyevt0",100,0,100);
  TH1F *hchisq1_toyevt1 = new TH1F("hchisq1_toyevt1","hchisq1_toyevt1",100,0,100);
  TH1F *hchisq1_toyevt2 = new TH1F("hchisq1_toyevt2","hchisq1_toyevt2",100,0,100);
  THStack *schisq1 = new THStack("schisq1","schisq1");

  TH1F *hchisq2 = new TH1F("hchisq2","hchisq2",100,0,100);
  TH1F *hchisq2_toyevt0 = new TH1F("hchisq2_toyevt0","hchisq2_toyevt0",100,0,100);
  TH1F *hchisq2_toyevt1 = new TH1F("hchisq2_toyevt1","hchisq2_toyevt1",100,0,100);
  TH1F *hchisq2_toyevt2 = new TH1F("hchisq2_toyevt2","hchisq0_toyevt2",100,0,100);
  THStack *schisq2 = new THStack("schisq2","schisq2");

  TH1F *hchisqminFullChisq = new TH1F("hchisqminFullChisq","hchisqminFullChisq",100,0,100);

  TH1F *hErrMINOSneg0 = new TH1F("hErrMINOSneg0","hErrMINOSneg0",1000,0,5.);
  TH1F *hErrMINOSneg1 = new TH1F("hErrMINOSneg1","hErrMINOSneg1",1000,0.,5.);
  TH1F *hErrMINOSpos0 = new TH1F("hErrMINOSpos0","hErrMINOSpos0",1000,0,5.);
  TH1F *hErrMINOSpos1 = new TH1F("hErrMINOSpos1","hErrMINOSpos1",1000,0.,5.);
  TH1F *hErrHESSEsym0 = new TH1F("hErrHESSEsym0","hErrHESSEsym0",1000,0,5.);
  TH1F *hErrHESSEsym1 = new TH1F("hErrHESSEsym1","hErrHESSEsym1",1000,0.,5.);

  TH1F *hnchisqcalls = new TH1F("hnchisqcalls","hnchisqcalls",100,0,700);
  TH1F *hfitstatus = new TH1F("hfitstatus","hfitstatus",5,0,5);
  int NdiffChisq = 0;

  TH1F *hpullsHESSE0 = new TH1F("hpullsHESSE0","hpullsHESSE0",100,-5,5);
  TH1F *hpullsHESSE1 = new TH1F("hpullsHESSE1","hpullsHESSE1",100,-5,5);
  TH1F *hpullsMINOS0 = new TH1F("hpullsMINOS0","hpullsMINOS0",100,-5,5);
  TH1F *hpullsMINOS1 = new TH1F("hpullsMINOS1","hpullsMINOS1",100,-5,5);

  TH1F *hresidual0 = new TH1F("hresidual0","hresidual0",100,-5,5);
  TH1F *hresidual1 = new TH1F("hresidual1","hresidual1",100,-5,5);

  // Check if 1D or 2D fit
  bool IsFit2D;
  treein->GetEntry(0);
  if(bestfit1==-999) IsFit2D = false;
  else               IsFit2D = true;

  // loop over the tree  
  
  int Nentries = treein->GetEntries();
  for(int ientry=0;ientry<Nentries;ientry++){

    treein->GetEntry(ientry);

    double nevents = (double)hErecoToy->Integral();    
    if(nevents<0){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events = " << nevents << " < 0" << endl;
      cerr << endl;
      exit(1);
    }
    
    //if(deltachisqtrue>-4) continue;
    //if(deltachisqtrue>-2 || deltachisqtrue<-2.2) continue;
    //if(deltachisqtrue<0. || deltachisqtrue>0.2) continue;
    //if( fabs(chisqmin-7.797102)>0.00001 ) continue;
    // if( (chisqtrue0 - chisqmin) > 0.4 ) continue;
    // else{
    //   cout << std::setprecision(8);
    //   cout << endl;
    //   cout << "Toy " << ientry << " --> ToyID = " << ToyID << endl;
    //   cout << "# events = " << nevents << endl;
    //   cout << "chisqmin = " << chisqmin
    // 	   << " - bestfit0 = " << bestfit0
    // 	   << " - bestfit1 = " << bestfit1 << endl;
    //   cout << "deltachisqtrue = " << deltachisqtrue << endl;
    // }


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
   
    hErrMINOSneg0->Fill(bestfiterrNeg0);
    hErrMINOSpos0->Fill(bestfiterrPos0);
    hErrMINOSneg1->Fill(bestfiterrNeg1);
    hErrMINOSpos1->Fill(bestfiterrPos1);
    hErrHESSEsym0->Fill(bestfiterrSym0);
    hErrHESSEsym1->Fill(bestfiterrSym1);

    double diff = chisqmin-chisqminFullChisq;
    hchisqminDiff->Fill(diff);

    double deltachisq_0BF = chisqtrue0 - chisqmin;
    hdeltachisqDiff0BF->Fill(deltachisq_0BF);  
    if(nevents==0)      hdeltachisqDiff0BF_toyevt0->Fill(deltachisq_0BF);
    else if(nevents==1) hdeltachisqDiff0BF_toyevt1->Fill(deltachisq_0BF);
    else if(nevents==2) hdeltachisqDiff0BF_toyevt2->Fill(deltachisq_0BF);

    double deltachisq_1BF = chisqtrue1 - chisqmin;
    hdeltachisqDiff1BF->Fill(deltachisq_1BF);
    if(nevents==0)      hdeltachisqDiff1BF_toyevt0->Fill(deltachisq_1BF);
    else if(nevents==1) hdeltachisqDiff1BF_toyevt1->Fill(deltachisq_1BF);
    else if(nevents==2) hdeltachisqDiff1BF_toyevt2->Fill(deltachisq_1BF);
    
    double deltachisq_01 = chisqtrue0 - chisqtrue1;
    hdeltachisqDiff01->Fill(deltachisq_01);
    if(nevents==0)      hdeltachisqDiff01_toyevt0->Fill(deltachisq_01);
    else if(nevents==1) hdeltachisqDiff01_toyevt1->Fill(deltachisq_01);
    else if(nevents==2) hdeltachisqDiff01_toyevt2->Fill(deltachisq_01);

    // Minimized chisq distribution
    hchisqmin->Fill(chisqmin);
    if(nevents==0)      hchisqmin_toyevt0->Fill(chisqmin);
    else if(nevents==1) hchisqmin_toyevt1->Fill(chisqmin);
    else if(nevents==2) hchisqmin_toyevt2->Fill(chisqmin);

    hchisq0->Fill(chisqtrue0);
    if(nevents==0)      hchisq0_toyevt0->Fill(chisqtrue0);
    else if(nevents==1) hchisq0_toyevt1->Fill(chisqtrue0);
    else if(nevents==2) hchisq0_toyevt2->Fill(chisqtrue0);

    hchisq1->Fill(chisqtrue1);
    if(nevents==0)      hchisq1_toyevt0->Fill(chisqtrue1);
    else if(nevents==1) hchisq1_toyevt1->Fill(chisqtrue1);
    else if(nevents==2) hchisq1_toyevt2->Fill(chisqtrue1);

    hbestfit_1D->Fill(bestfit0);
    if(nevents==0)      hbestfit_1D_toyevt0->Fill(bestfit0);
    else if(nevents==1) hbestfit_1D_toyevt1->Fill(bestfit0);
    else if(nevents==2) hbestfit_1D_toyevt2->Fill(bestfit0);

    // Full chisq distribution
    hchisqminFullChisq->Fill(chisqminFullChisq);
    hbestfitFullChisq_1D->Fill(bestfitFullChisq0);
  
  }

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

  // TCanvas *cnchisqcalls = new TCanvas("cnchisqcalls","cnchisqcalls");
  // hnchisqcalls->GetYaxis()->SetTitle("# entries");
  // hnchisqcalls->GetXaxis()->SetTitle("# of #chi^{2} computations");
  // hnchisqcalls->DrawCopy();

  // TCanvas *cfitstatus = new TCanvas("cfitstatus","cfitstatus");
  // hfitstatus->GetYaxis()->SetTitle("# entries");
  // hfitstatus->GetXaxis()->SetTitle("Fit Status");
  // hfitstatus->DrawCopy();

  // TCanvas *cbestfitFullChisq = new TCanvas("hbestfitFullChisq","hbestfitFullChisq");
  //   hbestfitFullChisq_1D->GetYaxis()->SetTitle("# entries");
  //   hbestfitFullChisq_1D->GetXaxis()->SetTitle("#bar{#beta}");
  //   hbestfitFullChisq_1D->DrawCopy();

  // TCanvas *cchisqminDiff = new TCanvas("hchisqminDiff","hchisqminDiff");
  // hchisqminDiff->GetYaxis()->SetTitle("# entries");
  // hchisqminDiff->GetXaxis()->SetTitle("#chi^{2}_{Descent} - #chi^{2}_{Global}");
  // hchisqminDiff->DrawCopy();

  TCanvas *cErrMINOSneg0 = new TCanvas("cErrMINOSneg0","cErrMINOSneg0");
  hErrMINOSneg0->GetYaxis()->SetTitle("# entries");
  hErrMINOSneg0->GetXaxis()->SetTitle("#Delta #bar{#beta} (MINOS negative)");
  hErrMINOSneg0->DrawCopy();
  TCanvas *cErrMINOSpos0 = new TCanvas("cErrMINOSpos0","cErrMINOSpos0");
  hErrMINOSpos0->GetYaxis()->SetTitle("# entries");
  hErrMINOSpos0->GetXaxis()->SetTitle("#Delta #bar{#beta} (MINOS positive)");
  hErrMINOSpos0->DrawCopy();

  TCanvas *cErrMINOSneg1 = new TCanvas("cErrMINOSneg1","cErrMINOSneg1");
  hErrMINOSneg1->GetYaxis()->SetTitle("# entries");
  hErrMINOSneg1->GetXaxis()->SetTitle("#delta_{CP} (MINOS positive)");
  hErrMINOSneg1->DrawCopy();
  TCanvas *cErrMINOSpos1 = new TCanvas("cErrMINOSpos1","cErrMINOSpos1");
  hErrMINOSpos1->GetYaxis()->SetTitle("# entries");
  hErrMINOSpos1->GetXaxis()->SetTitle("#delta_{CP} (MINOS negative)");
  hErrMINOSpos1->DrawCopy();

  TCanvas *cErrHESSEsym0 = new TCanvas("cErrHESSEsym0","cErrHESSEsym0");
  hErrHESSEsym0->GetYaxis()->SetTitle("# entries");
  hErrHESSEsym0->GetXaxis()->SetTitle("#Delta #bar{#beta} (HESSE)");
  hErrHESSEsym0->DrawCopy();
  TCanvas *cErrHESSEsym1 = new TCanvas("cErrHESSEsym1","cErrHESSEsym1");
  hErrHESSEsym1->GetYaxis()->SetTitle("# entries");
  hErrHESSEsym1->GetXaxis()->SetTitle("#delta_{CP} (HESSE)");
  hErrHESSEsym1->DrawCopy();

  TCanvas *cpullsHESSE0 = new TCanvas("cpullsHESSE0","cpullsHESSE0");
  hpullsHESSE0->GetXaxis()->SetTitle("Pull #bar{#beta} (HESSE error)");
  hpullsHESSE0->GetYaxis()->SetTitle("# entries");
  hpullsHESSE0->DrawCopy();
  sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsHESSE0->GetMean(),hpullsHESSE0->GetRMS(),FitMeanPar0_HESSE,FitRMSPar0_HESSE); 
  TLatex *lpullsHESSE0 = new TLatex(0.90,hpullsHESSE0->GetMaximum()*0.8,name);
  lpullsHESSE0->SetTextSize(0.07);
  lpullsHESSE0->Draw("same");
  
  if(IsFit2D){
    TCanvas *cpullsHESSE1 = new TCanvas("cpullsHESSE1","cpullsHESSE1");
    hpullsHESSE1->GetXaxis()->SetTitle("Pull sin^{2} #bar{#theta}_{13} (HESSE error)");
    hpullsHESSE1->GetYaxis()->SetTitle("# entries");
    hpullsHESSE1->DrawCopy();
    sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsHESSE1->GetMean(),hpullsHESSE1->GetRMS(),FitMeanPar1_HESSE,FitRMSPar1_HESSE); 
    TLatex *lpullsHESSE1 = new TLatex(0.90,hpullsHESSE1->GetMaximum()*0.8,name);
    lpullsHESSE1->SetTextSize(0.07);
    lpullsHESSE1->Draw();
  }

  TCanvas *cpullsMINOS0 = new TCanvas("cpullsMINOS0","cpullsMINOS0");
  hpullsMINOS0->GetXaxis()->SetTitle("Pull #bar{#beta} (MINOS error)");
  hpullsMINOS0->GetYaxis()->SetTitle("# entries");
  hpullsMINOS0->DrawCopy();
  sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsMINOS0->GetMean(),hpullsMINOS0->GetRMS(),FitMeanPar0_MINOS,FitRMSPar0_MINOS); 
  TLatex *lpullsMINOS0 = new TLatex(0.90,hpullsMINOS0->GetMaximum()*0.8,name);
  lpullsMINOS0->SetTextSize(0.07);
  lpullsMINOS0->Draw();

  if(IsFit2D){
    TCanvas *cpullsMINOS1 = new TCanvas("cpullsMINOS1","cpullsMINOS1");
    hpullsMINOS1->GetXaxis()->SetTitle("Pull #bar{#beta} (MINOS error)");
    hpullsMINOS1->GetYaxis()->SetTitle("# entries");
    hpullsMINOS1->DrawCopy();
    sprintf(name,"#splitline{#splitline{Mean=%1.2f}{RMS=%1.2f}}{#splitline{Fit Mean=%1.2f}{Fit RMS=%1.2f}}",hpullsMINOS1->GetMean(),hpullsMINOS1->GetRMS(),FitMeanPar1_MINOS,FitRMSPar1_MINOS); 
    TLatex *lpullsMINOS1 = new TLatex(0.90,hpullsMINOS1->GetMaximum()*0.8,name);
    lpullsMINOS1->SetTextSize(0.07);
    lpullsMINOS1->Draw();
  }

  TCanvas *cresidual0 = new TCanvas("cresidual0","cresidual0");
  hresidual0->GetXaxis()->SetTitle("Residual #bar{#beta}");
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

  TCanvas *cbestfit = new TCanvas("cbestfit","cbestfit");
  hbestfit_1D->GetYaxis()->SetTitle("# entries");
  hbestfit_1D->GetXaxis()->SetTitle("#bar{#beta}");
  hbestfit_1D->DrawCopy();
  hbestfit_1D_toyevt0->SetFillColor(kGreen);  
  hbestfit_1D_toyevt1->SetFillColor(kBlue);  
  hbestfit_1D_toyevt2->SetFillColor(kCyan);  
  sbestfit_1D->Add(hbestfit_1D_toyevt0);
  sbestfit_1D->Add(hbestfit_1D_toyevt1);
  sbestfit_1D->Add(hbestfit_1D_toyevt2);
  sbestfit_1D->Draw("same");
  cbestfit->SetLogy();

  TCanvas *cchisqmin = new TCanvas("cchisqmin","cchisqmin");
  hchisqmin->GetYaxis()->SetTitle("# entries");
  hchisqmin->GetXaxis()->SetTitle("#chi^{2} (#bar{#beta}_{BF})");
  hchisqmin->DrawCopy();
  hchisqmin_toyevt0->SetFillColor(kGreen);  
  hchisqmin_toyevt1->SetFillColor(kBlue);  
  hchisqmin_toyevt2->SetFillColor(kCyan);  
  schisqmin->Add(hchisqmin_toyevt0);
  schisqmin->Add(hchisqmin_toyevt1);
  schisqmin->Add(hchisqmin_toyevt2);
  schisqmin->Draw("same");
  cchisqmin->SetLogy();

  TCanvas *cchisq0 = new TCanvas("cchisq0","cchisq0");
  hchisq0->GetYaxis()->SetTitle("# entries");
  hchisq0->GetXaxis()->SetTitle("#chi^{2} (#bar{#beta=0})");
  hchisq0->DrawCopy();
  hchisq0_toyevt0->SetFillColor(kGreen);  
  hchisq0_toyevt1->SetFillColor(kBlue);  
  hchisq0_toyevt2->SetFillColor(kCyan);  
  schisq0->Add(hchisq0_toyevt0);
  schisq0->Add(hchisq0_toyevt1);
  schisq0->Add(hchisq0_toyevt2);
  schisq0->Draw("same");
  cchisq0->SetLogy();

  TCanvas *cchisq1 = new TCanvas("cchisq1","cchisq1");
  hchisq1->GetYaxis()->SetTitle("# entries");
  hchisq1->GetXaxis()->SetTitle("#chi^{2} (#bar{#beta=1})");
  hchisq1->DrawCopy();
  hchisq1_toyevt0->SetFillColor(kGreen);  
  hchisq1_toyevt1->SetFillColor(kBlue);  
  hchisq1_toyevt2->SetFillColor(kCyan);  
  schisq1->Add(hchisq1_toyevt0);
  schisq1->Add(hchisq1_toyevt1);
  schisq1->Add(hchisq1_toyevt2);
  schisq1->Draw("same");
  cchisq1->SetLogy();

  TCanvas *cdeltachisqDiff0BF = new TCanvas("cdeltachisqDiff0BF","cdeltachisqDiff0BF");
  hdeltachisqDiff0BF->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff0BF->GetXaxis()->SetTitle("#chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}_{BF})");
  hdeltachisqDiff0BF->DrawCopy();
  hdeltachisqDiff0BF_toyevt0->SetFillColor(kGreen);  
  hdeltachisqDiff0BF_toyevt1->SetFillColor(kBlue);  
  hdeltachisqDiff0BF_toyevt2->SetFillColor(kCyan);  
  sdeltachisqDiff0BF->Add(hdeltachisqDiff0BF_toyevt0);
  sdeltachisqDiff0BF->Add(hdeltachisqDiff0BF_toyevt1);
  sdeltachisqDiff0BF->Add(hdeltachisqDiff0BF_toyevt2);
  sdeltachisqDiff0BF->Draw("same");
  cdeltachisqDiff0BF->SetLogy();

  TCanvas *cdeltachisqDiff1BF = new TCanvas("cdeltachisqDiff1BF","cdeltachisqDiff1BF");
  hdeltachisqDiff1BF->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff1BF->GetXaxis()->SetTitle("#chi^{2}(#bar{#beta}=1) - #chi^{2}(#bar{#beta}_{BF})");
  hdeltachisqDiff1BF->DrawCopy();
  hdeltachisqDiff1BF_toyevt0->SetFillColor(kGreen);  
  hdeltachisqDiff1BF_toyevt1->SetFillColor(kBlue);  
  hdeltachisqDiff1BF_toyevt2->SetFillColor(kCyan);  
  sdeltachisqDiff1BF->Add(hdeltachisqDiff1BF_toyevt0);
  sdeltachisqDiff1BF->Add(hdeltachisqDiff1BF_toyevt1);
  sdeltachisqDiff1BF->Add(hdeltachisqDiff1BF_toyevt2);
  sdeltachisqDiff1BF->Draw("same");
  cdeltachisqDiff1BF->SetLogy();

  TCanvas *cdeltachisqDiff01 = new TCanvas("cdeltachisqDiff01","cdeltachisqDiff01");
  hdeltachisqDiff01->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff01->GetXaxis()->SetTitle("#chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)");
  hdeltachisqDiff01->DrawCopy();
  hdeltachisqDiff01_toyevt0->SetFillColor(kGreen);  
  hdeltachisqDiff01_toyevt1->SetFillColor(kBlue);  
  hdeltachisqDiff01_toyevt2->SetFillColor(kCyan);  
  sdeltachisqDiff01->Add(hdeltachisqDiff01_toyevt0);
  sdeltachisqDiff01->Add(hdeltachisqDiff01_toyevt1);
  sdeltachisqDiff01->Add(hdeltachisqDiff01_toyevt2);
  sdeltachisqDiff01->Draw("same");
  cdeltachisqDiff01->SetLogy();
  

  // 
  
  treein->GetEntry(493);
  //treein->GetEntry(4739);
  TCanvas *cchisq_out_GetEntry = new TCanvas("cchisq_out_GetEntry","cchisq_out_GetEntry");
  if(IsFit2D) hchisq2D_out->DrawCopy("colz");
  else        hchisq1D_out->DrawCopy();

  TCanvas *cErecoToy_GetEntry = new TCanvas("cErecoToy_GetEntry","cErecoToy_GetEntry");
  hErecoToy->DrawCopy();




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

  return;
}

