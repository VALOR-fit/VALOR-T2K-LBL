
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

void analyze_run_fit_marginalization_readbinary_nuebarapp_rateonly
(
 //char infilename[200] = "../../outputs/marginalization/nuebarapp/fitmarg_T2Kreact_dcp_23_Valor_MargMH_4p5E20POT_2-3f+normbar_10kToys.root",
 char infilename[200] = "fit_beta1_nexp_10000_1.000_-999.000-job0-_0_9999.root",
 bool draw_predicted_event_distribution = false // Only for 1st fit
)
{
  char name[200];
  
  TFile *filein = new TFile(infilename,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  
  
  PhysParamList *iparm;
  FitParamList *fparm;
  TH1F *hchisq1D_out;
  int ToyID=-999;
  TH2F *hErecoToy=0;
  TH1D* hbeta0 = 0;
  TH1D* hbeta1 = 0;
  
  
  treein->SetBranchAddress("best_fit_params",&fparm);
  treein->SetBranchAddress("input_params",&iparm);
  treein->SetBranchAddress("hchisq1D", &hchisq1D_out );
  treein->SetBranchAddress("ToyID",&ToyID);
  treein->SetBranchAddress("hErecoToy_RHC1Re",&hErecoToy);
  treein->SetBranchAddress("hBetaBar0",&hbeta0);
  treein->SetBranchAddress("hBetaBar1",&hbeta1);
  
  
  // Declare histograms
  
  TH1F *hPValueBeta0 = new TH1F("hPValueBeta0","hPValueBeta0",200,0,1.5);
  TH1F *hPValueBeta0_toyevt0 = new TH1F("hPValueBeta0_toyevt0","hPValueBeta0_toyevt0",200,0,1.5);
  TH1F *hPValueBeta0_toyevt1 = new TH1F("hPValueBeta0_toyevt1","hPValueBeta0_toyevt1",200,0,1.5);
  TH1F *hPValueBeta0_toyevt2 = new TH1F("hPValueBeta0_toyevt2","hPValueBeta0_toyevt2",200,0,1.5);
  THStack *sPValueBeta0 = new THStack("sPValueBeta0","sPValueBeta0");
  
  TH1F *hPValueBeta1 = new TH1F("hPValueBeta1","hPValueBeta1",200,0,1.5);
  TH1F *hPValueBeta1_toyevt0 = new TH1F("hPValueBeta1_toyevt0","hPValueBeta1_toyevt0",200,0,1.5);
  TH1F *hPValueBeta1_toyevt1 = new TH1F("hPValueBeta1_toyevt1","hPValueBeta1_toyevt1",200,0,1.5);
  TH1F *hPValueBeta1_toyevt2 = new TH1F("hPValueBeta1_toyevt2","hPValueBeta1_toyevt2",200,0,1.5);
  THStack *sPValueBeta1 = new THStack("sPValueBeta1","sPValueBeta1");
  
  
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
    
    
    double pvalue_beta0 = hchisq1D_out->GetBinContent(1);
    
    hPValueBeta0->Fill(pvalue_beta0);
    
    if(nevents==0)      hPValueBeta0_toyevt0->Fill(pvalue_beta0);
    else if(nevents==1) hPValueBeta0_toyevt1->Fill(pvalue_beta0);
    else if(nevents==2) hPValueBeta0_toyevt2->Fill(pvalue_beta0);
    
    double pvalue_beta1 = hchisq1D_out->GetBinContent(2);
    
    hPValueBeta1->Fill(pvalue_beta1);
    if(nevents==0)      hPValueBeta1_toyevt0->Fill(pvalue_beta1);
    else if(nevents==1) hPValueBeta1_toyevt1->Fill(pvalue_beta1);
    else if(nevents==2) hPValueBeta1_toyevt2->Fill(pvalue_beta1);
    
    cout << "P-Values:    Beta=0: " << pvalue_beta0 << ", Beta=1: " << pvalue_beta1 << endl;
  }
  
  
  // Draw histograms
  
  if(draw_predicted_event_distribution){
    treein->GetEntry(0); // Do it for 1st fit in tree
    double nevents = (double)hErecoToy->Integral();

    TCanvas *cBeta0Distribution = new TCanvas("cBeta0Distribution","cBeta0Distribution");
    hbeta0->GetYaxis()->SetTitle("# entries");
    hbeta0->GetXaxis()->SetTitle("# events");
    hbeta0->SetLineColor(kBlue);
    hbeta0->SetLineWidth(5);
    hbeta0->GetYaxis()->SetTitleSize(0.06);
    hbeta0->GetXaxis()->SetTitleSize(0.06);
    hbeta0->GetYaxis()->SetTitleFont(132);
    hbeta0->GetXaxis()->SetTitleFont(132);
    hbeta0->GetYaxis()->SetLabelFont(132);
    hbeta0->GetXaxis()->SetLabelFont(132);
    hbeta0->GetYaxis()->SetLabelSize(0.05);
    hbeta0->GetXaxis()->SetLabelSize(0.05);
    hbeta0->DrawCopy();
    cBeta0Distribution->SetLogy();
    TLine *expEvents = new TLine(nevents,hbeta0->GetMinimum(),nevents,hbeta0->GetMaximum());
    expEvents->SetLineWidth(4);
    expEvents->SetLineColor(kRed);
    expEvents->SetLineStyle(2);
    expEvents->Draw();
    TLegend *leg = new TLegend(0.54,0.71,0.89,0.91);
    leg->SetBorderSize(1);
    leg->SetLineWidth(3);
    leg->SetTextFont(132);
    leg->SetTextSize(0.05);
    leg->AddEntry(expEvents,"# Events data","L");
    leg->AddEntry(hbeta0,"# Events prediction","L");
    leg->Draw();

    
    TCanvas *cBeta1Distribution = new TCanvas("cBeta1Distribution","cBeta1Distribution");
    hbeta1->GetYaxis()->SetTitle("# entries");
    hbeta1->GetXaxis()->SetTitle("# events");
    hbeta1->GetYaxis()->SetTitleSize(0.06);
    hbeta1->GetXaxis()->SetTitleSize(0.06);
    hbeta1->GetYaxis()->SetTitleFont(132);
    hbeta1->GetXaxis()->SetTitleFont(132);
    hbeta1->GetYaxis()->SetLabelFont(132);
    hbeta1->GetXaxis()->SetLabelFont(132);
    hbeta1->GetYaxis()->SetLabelSize(0.05);
    hbeta1->GetXaxis()->SetLabelSize(0.05);

    hbeta1->SetLineColor(kBlue);
    hbeta1->SetLineWidth(5);
    hbeta1->DrawCopy();
    cBeta1Distribution->SetLogy();
    expEvents->Draw();
    leg->Draw();

  }
  
  
  
  TCanvas *cPValueBeta0 = new TCanvas("cPValueBeta0","cPValueBeta0");
  hPValueBeta0->GetYaxis()->SetTitle("# entries");
  hPValueBeta0->GetXaxis()->SetTitle("p-value");
  hPValueBeta0->SetLineWidth(2);
  hPValueBeta0->DrawCopy();
  hPValueBeta0_toyevt0->SetFillColor(kBlack);
  hPValueBeta0_toyevt0->SetLineColor(kBlack);
  hPValueBeta0_toyevt1->SetFillColor(kBlue);
  hPValueBeta0_toyevt1->SetLineColor(kBlue);
  hPValueBeta0_toyevt2->SetFillColor(kGreen+2);
  hPValueBeta0_toyevt2->SetLineColor(kGreen+2);

  sPValueBeta0->Add(hPValueBeta0_toyevt0);
  sPValueBeta0->Add(hPValueBeta0_toyevt1);
  sPValueBeta0->Add(hPValueBeta0_toyevt2);
  sPValueBeta0->Draw("same");
  hPValueBeta0->DrawCopy("same");

  cPValueBeta0->SetLogy();
  
  TLegend *leg = new TLegend(0.66,0.70,0.89,0.94);
  leg->SetBorderSize(1);
  leg->SetLineWidth(2);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);
  leg->AddEntry(hPValueBeta0,"Total","L");
  leg->AddEntry(hPValueBeta0_toyevt0,"# Events = 0","F");
  leg->AddEntry(hPValueBeta0_toyevt1,"# Events = 1","F");
  leg->AddEntry(hPValueBeta0_toyevt2,"# Events = 2","F");
  leg->Draw();
  
  TCanvas *cPValueBeta1 = new TCanvas("cPValueBeta1","cPValueBeta1");
  hPValueBeta1->GetYaxis()->SetTitle("# entries");
  hPValueBeta1->GetXaxis()->SetTitle("p-value");
  hPValueBeta1->SetLineWidth(2);
  hPValueBeta1->DrawCopy();
  hPValueBeta1_toyevt0->SetFillColor(kBlack);
  hPValueBeta1_toyevt0->SetLineColor(kBlack);
  hPValueBeta1_toyevt1->SetFillColor(kBlue);
  hPValueBeta1_toyevt1->SetLineColor(kBlue);
  hPValueBeta1_toyevt2->SetFillColor(kGreen+2);
  hPValueBeta1_toyevt2->SetLineColor(kGreen+2);
  sPValueBeta1->Add(hPValueBeta1_toyevt0);
  sPValueBeta1->Add(hPValueBeta1_toyevt1);
  sPValueBeta1->Add(hPValueBeta1_toyevt2);
  sPValueBeta1->Draw("same");
  hPValueBeta1->DrawCopy("same");
  cPValueBeta1->SetLogy();
  leg->Draw();
  
  
  return;
}

