
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
// Compute p-values using deltaChisq_01, deltaChisq_0BF and # of events
// Store final results in a output root file to be used by other macros
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;

void CalcPvalue_nuebarapp
(
 // Mach3 - 4.5 x 10^20 POT

 // rate+shape (current stat), dcp toys
 TString tag = "pvalue_beta1_Mach3",
 double data = 4.15309, // # of data events
 double deltachisq_beta01_data = 1.9467, // data deltachisq 01
 double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Mach3_MargMH_4p5E20POT-3f+normbar_10kToys.root",


 // 4 x 10^20 POT

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta0",
 // double data = 1.75941, // # of data events
 // double deltachisq_beta01_data = -2.04059, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4E20POT-3f+normbar_10kToys.root",

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta0p5",
 // double data = 2.72552, // # of data events
 // double deltachisq_beta01_data = -0.145486, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4E20POT-3f+normbar_10kToys.root",

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta1",
 // double data = 3.69163, // # of data events
 // double deltachisq_beta01_data = 1.76424, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4E20POT-3f+normbar_10kToys.root",
 
 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta1p5",
 // double data = 4.65774, // # of data events
 // double deltachisq_beta01_data = 3.68919, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4E20POT-3f+normbar_10kToys.root",

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta2",
 // double data = 5.62385, // # of data events
 // double deltachisq_beta01_data = 5.62982, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4E20POT-3f+normbar_10kToys.root",

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta2p5",
 // double data = 6.58996, // # of data events
 // double deltachisq_beta01_data = 7.58645, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4E20POT-3f+normbar_10kToys.root",




 // 4.5 x 10^20 POT

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta0",
 // double data = 1.97934, // # of data events
 // double deltachisq_beta01_data = -2.2757, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // //char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_4900Toys.root",
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_10kToys_2.root",

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta0p5",
 // double data = 3.06621, // # of data events
 // double deltachisq_beta01_data = -0.156778, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // //char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_4900Toys.root",
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_10kToys_2.root",

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta1",
 // double data = 4.15309, // # of data events
 // double deltachisq_beta01_data = 1.98008, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // //char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_4900Toys.root",
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_10kToys_2.root",
 
 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta1p5",
 // double data = 5.23996, // # of data events
 // double deltachisq_beta01_data = 4.13569, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // //char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_4900Toys.root",
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_10kToys_2.root",

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta2",
 // double data = 6.32683, // # of data events
 // double deltachisq_beta01_data = 6.31073, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // //char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_4900Toys.root",
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_10kToys_2.root",

 // // rate+shape (current stat), dcp toys
 // TString tag = "pvalue_beta2p5",
 // double data = 7.41371, // # of data events
 // double deltachisq_beta01_data = 8.50564, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // //char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_4900Toys.root",
 // char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT-3f+normbar_10kToys_2.root",

 )
{

  // const double MaxHisto = 20.;
  // const double MinHisto = -20.;
  // const double NBinHisto = 200.;

  const double MaxHisto = 200.;
  const double MinHisto = -200.;
  const double NBinHisto = 2000.;
  
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
  
  TH1F *hnevents = new TH1F("hnevents","hnevents",NBinHisto,0.,MaxHisto);
  TH1F *hnevents_data = new TH1F("hnevents_data","hnevents_data",NBinHisto,0,MaxHisto);

  TH1F *hdeltachisqDiff0BF = new TH1F("hdeltachisqDiff0BF","hdeltachisqDiff0BF",NBinHisto,0,MaxHisto);
  TH1F *hdeltachisqDiff0BF_toyevt0 = new TH1F("hdeltachisqDiff0BF_toyevt0","hdeltachisqDiff0BF_toyevt0",NBinHisto,0,MaxHisto);
  TH1F *hdeltachisqDiff0BF_toyevt1 = new TH1F("hdeltachisqDiff0BF_toyevt1","hdeltachisqDiff0BF_toyevt1",NBinHisto,0,MaxHisto);
  TH1F *hdeltachisqDiff0BF_toyevt2 = new TH1F("hdeltachisqDiff0BF_toyevt2","hdeltachisqDiff0BF_toyevt2",NBinHisto,0,MaxHisto);
  THStack *sdeltachisqDiff0BF = new THStack("sdeltachisqDiff0BF","sdeltachisqDiff0BF");

  TH1F *hdeltachisqDiff1BF = new TH1F("hdeltachisqDiff1BF","hdeltachisqDiff1BF",NBinHisto,0,MaxHisto);
  TH1F *hdeltachisqDiff1BF_toyevt0 = new TH1F("hdeltachisqDiff1BF_toyevt0","hdeltachisqDiff1BF_toyevt0",NBinHisto,0,MaxHisto);
  TH1F *hdeltachisqDiff1BF_toyevt1 = new TH1F("hdeltachisqDiff1BF_toyevt1","hdeltachisqDiff1BF_toyevt1",NBinHisto,0,MaxHisto);
  TH1F *hdeltachisqDiff1BF_toyevt2 = new TH1F("hdeltachisqDiff1BF_toyevt2","hdeltachisqDiff1BF_toyevt2",NBinHisto,0,MaxHisto);
  THStack *sdeltachisqDiff1BF = new THStack("sdeltachisqDiff1BF","sdeltachisqDiff1BF");

  TH1F *hdeltachisqDiff01 = new TH1F("hdeltachisqDiff01","hdeltachisqDiff01",NBinHisto,MinHisto,MaxHisto);
  TH1F *hdeltachisqDiff01_toyevt0 = new TH1F("hdeltachisqDiff01_toyevt0","hdeltachisqDiff01_toyevt0",NBinHisto,MinHisto,MaxHisto);
  TH1F *hdeltachisqDiff01_toyevt1 = new TH1F("hdeltachisqDiff01_toyevt1","hdeltachisqDiff01_toyevt1",NBinHisto,MinHisto,MaxHisto);
  TH1F *hdeltachisqDiff01_toyevt2 = new TH1F("hdeltachisqDiff01_toyevt2","hdeltachisqDiff01_toyevt2",NBinHisto,MinHisto,MaxHisto);
  THStack *sdeltachisqDiff01 = new THStack("sdeltachisqDiff01","sdeltachisqDiff01");

  TH1F *hchisqmin = new TH1F("hchisqmin","hchisqmin",400,0,100);
  TH1F *hchisqmin_toyevt0 = new TH1F("hchisqmin_toyevt0","hchisqmin_toyevt0",400,0,100);
  TH1F *hchisqmin_toyevt1 = new TH1F("hchisqmin_toyevt1","hchisqmin_toyevt1",400,0,100);
  TH1F *hchisqmin_toyevt2 = new TH1F("hchisqmin_toyevt2","hchisqmin_toyevt2",400,0,100);
  THStack *schisqmin = new THStack("schisqmin","schisqmin");

  TH1F *hchisq0 = new TH1F("hchisq0","hchisq0",400,0,100);
  TH1F *hchisq0_toyevt0 = new TH1F("hchisq0_toyevt0","hchisq0_toyevt0",400,0,100);
  TH1F *hchisq0_toyevt1 = new TH1F("hchisq0_toyevt1","hchisq0_toyevt1",400,0,100);
  TH1F *hchisq0_toyevt2 = new TH1F("hchisq0_toyevt2","hchisq0_toyevt2",400,0,100);
  THStack *schisq0 = new THStack("schisq0","schisq0");

  TH1F *hchisq1 = new TH1F("hchisq1","hchisq1",400,0,100);
  TH1F *hchisq1_toyevt0 = new TH1F("hchisq1_toyevt0","hchisq1_toyevt0",400,0,100);
  TH1F *hchisq1_toyevt1 = new TH1F("hchisq1_toyevt1","hchisq1_toyevt1",400,0,100);
  TH1F *hchisq1_toyevt2 = new TH1F("hchisq1_toyevt2","hchisq1_toyevt2",400,0,100);
  THStack *schisq1 = new THStack("schisq1","schisq1");

  TH1F *hchisq2 = new TH1F("hchisq2","hchisq2",400,0,100);
  TH1F *hchisq2_toyevt0 = new TH1F("hchisq2_toyevt0","hchisq2_toyevt0",400,0,100);
  TH1F *hchisq2_toyevt1 = new TH1F("hchisq2_toyevt1","hchisq2_toyevt1",400,0,100);
  TH1F *hchisq2_toyevt2 = new TH1F("hchisq2_toyevt2","hchisq0_toyevt2",400,0,100);
  THStack *schisq2 = new THStack("schisq2","schisq2");

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
    hnevents->Fill(nevents);
    
    if(nevents<data) hnevents_data->Fill(nevents);
    
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

    // Fill vectors for pvalue calculations
    vec_evt_toys.push_back(nevents);
    vec_chisq_0BF.push_back(deltachisq_0BF);
    vec_chisq_01.push_back(deltachisq_01);
  }

  // Compute p-values

  double Ntoys = (double)Nentries;
  
  double pvalue_evt = GetPvalue(vec_evt_toys,Ntoys,data);
  double pvalue_01 = GetPvalue(vec_chisq_01,Ntoys,deltachisq_beta01_data);
  double pvalue_0BF = GetPvalue(vec_chisq_0BF,Ntoys,deltachisq_beta0BF_data);

 
  // Draw histograms
  
  // treein->GetEntry(493);
  // //treein->GetEntry(4739);
  // TCanvas *cchisq_out_GetEntry = new TCanvas("cchisq_out_GetEntry","cchisq_out_GetEntry");
  // if(IsFit2D) hchisq2D_out->DrawCopy("colz");
  // else        hchisq1D_out->DrawCopy();

  // TCanvas *cErecoToy_GetEntry = new TCanvas("cErecoToy_GetEntry","cErecoToy_GetEntry");
  // hErecoToy->DrawCopy();

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

  TCanvas *cPvalue_0BF = new TCanvas("cPvalue_0BF","cPvalue_0BF");
  cPvalue_0BF->SetLogy();
  hdeltachisqDiff0BF->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff0BF->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}_{BF})");
  hdeltachisqDiff0BF->DrawCopy();
  hdeltachisqDiff0BF_toyevt0->SetFillColor(kGreen);  
  hdeltachisqDiff0BF_toyevt1->SetFillColor(kBlue);  
  hdeltachisqDiff0BF_toyevt2->SetFillColor(kCyan);  
  sdeltachisqDiff0BF->Add(hdeltachisqDiff0BF_toyevt0);
  sdeltachisqDiff0BF->Add(hdeltachisqDiff0BF_toyevt1);
  sdeltachisqDiff0BF->Add(hdeltachisqDiff0BF_toyevt2);
  sdeltachisqDiff0BF->Draw("same"); 
  TLine *ldeltachisqDiff0BF = new TLine(deltachisq_beta0BF_data,0.,deltachisq_beta0BF_data,hdeltachisqDiff0BF->GetMaximum());
  ldeltachisqDiff0BF->SetLineColor(kRed);
  ldeltachisqDiff0BF->SetLineWidth(2);
  ldeltachisqDiff0BF->Draw();
  TLegend *lPvalue_0BF = new TLegend(0.5,0.6,0.9,0.9);
  lPvalue_0BF->SetFillColor(0);
  lPvalue_0BF->SetLineColor(0);
  lPvalue_0BF->SetTextFont(132);
  lPvalue_0BF->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_beta0BF_data); 
  lPvalue_0BF->AddEntry(ldeltachisqDiff0BF,name,"l");
  lPvalue_0BF->Draw();

  TCanvas *cPvalue_01 = new TCanvas("cPvalue_01","cPvalue_01");
  cPvalue_01->SetLogy();
  hdeltachisqDiff01->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff01->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)");
  hdeltachisqDiff01->DrawCopy();
  hdeltachisqDiff01_toyevt0->SetFillColor(kGreen);  
  hdeltachisqDiff01_toyevt1->SetFillColor(kBlue);  
  hdeltachisqDiff01_toyevt2->SetFillColor(kCyan);  
  sdeltachisqDiff01->Add(hdeltachisqDiff01_toyevt0);
  sdeltachisqDiff01->Add(hdeltachisqDiff01_toyevt1);
  sdeltachisqDiff01->Add(hdeltachisqDiff01_toyevt2);
  sdeltachisqDiff01->Draw("same"); 
  TLine *ldeltachisqDiff01 = new TLine(deltachisq_beta01_data,0.,deltachisq_beta01_data,hdeltachisqDiff01->GetMaximum());
  ldeltachisqDiff01->SetLineColor(kRed);
  ldeltachisqDiff01->SetLineWidth(2);
  ldeltachisqDiff01->Draw();
  TLegend *lPvalue_01 = new TLegend(0.5,0.6,0.9,0.9);
  lPvalue_01->SetFillColor(0);
  lPvalue_01->SetLineColor(0);
  lPvalue_01->SetTextFont(132);
  lPvalue_01->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_beta01_data); 
  lPvalue_01->AddEntry(ldeltachisqDiff01,name,"l");
  lPvalue_01->Draw();

  TCanvas *cPvalue_evt = new TCanvas("cPvalue_evt","cPvalue_evt");
  cPvalue_evt->SetLogy();
  hnevents->SetLineWidth(2);
  hnevents->GetXaxis()->SetTitle("# of observed events");
  hnevents->GetYaxis()->SetTitle("# of entries");
  hnevents_data->SetLineWidth(2);
  hnevents_data->SetFillColor(kBlue);
  hnevents->DrawCopy();
  hnevents_data->DrawCopy("same");
  TLine *ldata = new TLine(data,0.,data,hnevents->GetMaximum());
  ldata->SetLineColor(kRed);
  ldata->SetLineWidth(2);
  ldata->Draw();
  TLegend *lPvalues_evt = new TLegend(0.5,0.6,0.9,0.9);
  lPvalues_evt->SetFillColor(0);
  lPvalues_evt->SetLineColor(0);
  lPvalues_evt->SetTextFont(132);
  lPvalues_evt->SetTextSize(0.06);
  lPvalues_evt->AddEntry(hnevents,"# of toy events","l");
  sprintf(name,"# of data events = %1.3f",data); 
  lPvalues_evt->AddEntry(ldata,name,"l");
  lPvalues_evt->Draw();


  // Write output file
  
  TH1D *hNdata = new TH1D("hNdata","hNdata",1,0,1);
  hNdata->SetBinContent(1,data);
  TH1D *hdeltachisq_beta01 = new TH1D("hdeltachisq_beta01","hdeltachisq_beta01",1,0,1);
  hdeltachisq_beta01->SetBinContent(1,deltachisq_beta01_data);
  TH1D *hdeltachisq_beta0BF = new TH1D("hdeltachisq_beta0BF","hdeltachisq_beta0BF",1,0,1);
  hdeltachisq_beta0BF->SetBinContent(1,deltachisq_beta0BF_data);
  TH1D *hpvalue_rateonly = new TH1D("hpvalue_rateonly","hpvalue_rateonly",1,0,1);
  hpvalue_rateonly->SetBinContent(1,pvalue_evt);
  TH1D *hpvalue_rateshape01 = new TH1D("hpvalue_rateshape01","hpvalue_rateshape01",1,0,1);
  hpvalue_rateshape01->SetBinContent(1,pvalue_01);
  TH1D *hpvalue_rateshape0BF = new TH1D("hpvalue_rateshape0BF","hpvalue_rateshape0BF",1,0,1);
  hpvalue_rateshape0BF->SetBinContent(1,pvalue_0BF);

  char outname[200];
  sprintf(outname,"%s.root",tag.Data());
  TFile *outfile = new TFile(outname,"RECREATE");
  hNdata->Write();
  hdeltachisq_beta01->Write();
  hdeltachisq_beta0BF->Write();
  hpvalue_rateonly->Write();
  hpvalue_rateshape01->Write();
  hpvalue_rateshape0BF->Write();
  outfile->Close();

  // Print results
  
  cout << endl;
  cout << std::setprecision(10);
  cout << "truepar0 = " << truepar0 << endl;
  if(IsFit2D) cout << "truepar1 = " << truepar1 << endl;
  cout << endl;
  cout << "Calculate p-values:" << endl;
  cout << "p-value rate-only (# evt) = " << pvalue_evt << endl;
  cout << "p-value 01 (rate+shape) = " << pvalue_01 << endl;
  cout << "p-value 0BF (rate+shape) = " << pvalue_0BF << endl;
  cout << endl;

  return;
}



double GetPvalue(vector<double> vec_toys,int ntoys,double data){
  double ntot_toys = ntoys;
  if(ntot_toys != vec_toys.size()){    
    cout << endl;
    cout << "Warning:" << endl;
    cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
    cout << "Input total number of toys = " << ntot_toys << endl; 
    ntot_toys = vec_toys.size();
    cout << "Set total number of toys to vector size: " << ntot_toys << endl;
    cout << endl;
  }
  if(ntot_toys<=0.){
    cerr << endl;
    cerr << "ntot_toys = " << ntot_toys << endl;
    cerr << endl;
    exit(1);
  }
  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it >= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}
