
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

void DrawDeltaChisq_SigBkg_nuebar
(
 // Mach3 - 4.0108 x 10^20 POT (beta=1) --> DATA
 TString tag = "pvalue_beta1_Mach3_NH",
 double data = 3., // # of data events
 double deltachisq_beta01_data = 0.164213, // data deltachisq 01
 double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 char infilename_beta0[200] = "../../outputs/marginalization/nuebarapp/EPS_RELEASE_DEV_3/fitmarg_betabar_0_4.0108E20POT_EPS-3f+normbar_50kToys.root",
 char infilename_beta1[200] = "../../outputs/marginalization/nuebarapp/EPS_RELEASE_DEV_3/fitmarg_betabar_1_4.0108E20POT_EPS-3f+normbar_50kToys.root",
 )
{

  const double MaxHisto = 200.;
  const double MinHisto = -200.;
  const double NBinHisto = 2000.;
  
  // Declare histograms
  TH1F *hnevents_beta0 = new TH1F("hnevents_beta0","hnevents_beta0",NBinHisto,0.,MaxHisto);
  TH1F *hdeltachisqDiff01_beta0 = new TH1F("hdeltachisqDiff01_beta0","hdeltachisqDiff01_beta0",NBinHisto,MinHisto,MaxHisto);
  TH1F *hnevents_beta1 = new TH1F("hnevents_beta1","hnevents_beta1",NBinHisto,0.,MaxHisto);
  TH1F *hdeltachisqDiff01_beta1 = new TH1F("hdeltachisqDiff01_beta1","hdeltachisqDiff01_beta1",NBinHisto,MinHisto,MaxHisto);

  char name[200];
  double chisqtrue0=conv::kALargeNum;
  double chisqtrue1=conv::kALargeNum;
  TH1D *hErecoToy;
  
  // Beta=0

  TFile *filein = new TFile(infilename_beta0,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("hErecoToy",&hErecoToy);
  // loop over the tree  
  int Nentries = treein->GetEntries();
  //Nentries = 10000;
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
    hnevents_beta0->Fill(nevents);    
    double deltachisq_01 = chisqtrue0 - chisqtrue1;
    hdeltachisqDiff01_beta0->Fill(deltachisq_01);
  }


  // Beta=1

  filein = new TFile(infilename_beta1,"READ");
  treein = (TTree*)filein->Get("fitnt");
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("hErecoToy",&hErecoToy);
  // loop over the tree  
  int Nentries = treein->GetEntries();
  //Nentries = 10000;
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
    hnevents_beta1->Fill(nevents);    
    double deltachisq_01 = chisqtrue0 - chisqtrue1;
    hdeltachisqDiff01_beta1->Fill(deltachisq_01);
  }








  // Draw histograms

  TCanvas *cPvalue_01 = new TCanvas("cPvalue_01","cPvalue_01");
  cPvalue_01->SetLogy();
  hdeltachisqDiff01_beta0->Rebin(2);
  hdeltachisqDiff01_beta1->Rebin(2);
  hdeltachisqDiff01_beta0->GetYaxis()->SetTitle("# entries");
  hdeltachisqDiff01_beta0->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)");
  //hdeltachisqDiff01_beta0->SetFillColor(kGreen);
  hdeltachisqDiff01_beta0->SetLineColor(kGreen);
  hdeltachisqDiff01_beta1->SetLineColor(kRed);
  hdeltachisqDiff01_beta0->SetLineWidth(2);
  hdeltachisqDiff01_beta1->SetLineWidth(2);  
  hdeltachisqDiff01_beta0->DrawCopy();
  hdeltachisqDiff01_beta1->DrawCopy("same");  
  TLine *ldeltachisqDiff01 = new TLine(deltachisq_beta01_data,0.,deltachisq_beta01_data,hdeltachisqDiff01_beta0->GetMaximum());
  ldeltachisqDiff01->SetLineColor(kBlack);
  ldeltachisqDiff01->SetLineWidth(2);
  ldeltachisqDiff01->Draw();
  TLegend *lPvalue_01 = new TLegend(0.5,0.6,0.9,0.9);
  lPvalue_01->SetFillColor(0);
  lPvalue_01->SetLineColor(0);
  lPvalue_01->SetTextFont(132);
  lPvalue_01->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_beta01_data); 
  lPvalue_01->AddEntry(ldeltachisqDiff01,name,"l");
  lPvalue_01->AddEntry(hdeltachisqDiff01_beta0,"#Delta #chi^{2} (#bar{#beta}=0)","l");
  lPvalue_01->AddEntry(hdeltachisqDiff01_beta1,"#Delta #chi^{2} (#bar{#beta}=1)","l");
  lPvalue_01->Draw();

  TCanvas *cPvalue_rateonly = new TCanvas("cPvalue_rateonly","cPvalue_rateonly");
  cPvalue_rateonly->SetLogy();
  hnevents_beta0->Rebin(4);
  hnevents_beta1->Rebin(4);
  hnevents_beta0->GetYaxis()->SetTitle("# entries");
  hnevents_beta0->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)");
  hnevents_beta0->SetFillColor(kGreen);
  hnevents_beta0->SetLineColor(kGreen);
  hnevents_beta1->SetLineColor(kRed);
  hnevents_beta0->SetLineWidth(2);
  hnevents_beta1->SetLineWidth(2);  
  hnevents_beta0->DrawCopy();
  hnevents_beta1->DrawCopy("same");  
  TLine *lnevents = new TLine(data,0.,data,hnevents_beta0->GetMaximum());
  lnevents->SetLineColor(kBlack);
  lnevents->SetLineWidth(2);
  lnevents->Draw();
  TLegend *lPvalue_rateonly = new TLegend(0.5,0.6,0.9,0.9);
  lPvalue_rateonly->SetFillColor(0);
  lPvalue_rateonly->SetLineColor(0);
  lPvalue_rateonly->SetTextFont(132);
  lPvalue_rateonly->SetTextSize(0.06);
  sprintf(name,"# of data events",data); 
  lPvalue_rateonly->AddEntry(lnevents,name,"l");
  lPvalue_rateonly->AddEntry(hnevents_beta0,"\# of toy events (#bar{#beta}=0)","f");
  lPvalue_rateonly->AddEntry(hnevents_beta1,"\# of toy events (#bar{#beta}=1)","l");
  lPvalue_rateonly->Draw();

  return;
}

