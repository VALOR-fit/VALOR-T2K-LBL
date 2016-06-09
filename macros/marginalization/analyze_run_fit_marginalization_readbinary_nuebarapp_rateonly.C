
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
 bool draw_predicted_event_distribution = false, // Only for 1st fit
 double exp_bkg = 1.718
 )
{
  char name[200];
  
  TFile *filein = new TFile(infilename,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  
  
  PhysParamList *iparm;
  FitParamList *fparm;
  TH1D *hchisq1D_out;
  int ToyID=-999;
  TH2D *hErecoToy=0;
  TH1D* hbeta0 = 0;
  TH1D* hbeta1 = 0;
  TH1D* hbeta0Unscaled = 0;
  TH1D* hbeta1Unscaled = 0;
  std::vector<double> right_tail_weight;
  std::vector<double> left_tail_weight;
  
  
  treein->SetBranchAddress("best_fit_params",&fparm);
  treein->SetBranchAddress("input_params",&iparm);
  treein->SetBranchAddress("hchisq1D", &hchisq1D_out );
  treein->SetBranchAddress("ToyID",&ToyID);
  treein->SetBranchAddress("hErecoToy_RHC1Re",&hErecoToy);
  treein->SetBranchAddress("hBetaBar0",&hbeta0);
  treein->SetBranchAddress("hBetaBar1",&hbeta1);
  treein->SetBranchAddress("hBetaBar0Unscaled",&hbeta0Unscaled);
  treein->SetBranchAddress("hBetaBar1Unscaled",&hbeta1Unscaled);
  
  
  // Declare histograms
  
  TH1D *hPValueBeta0 = new TH1D("hPValueBeta0","hPValueBeta0",200,0,1.5);
  TH1D *hPValueBeta0_toyevt0 = new TH1D("hPValueBeta0_toyevt0","hPValueBeta0_toyevt0",200,0,1.5);
  TH1D *hPValueBeta0_toyevt1 = new TH1D("hPValueBeta0_toyevt1","hPValueBeta0_toyevt1",200,0,1.5);
  TH1D *hPValueBeta0_toyevt2 = new TH1D("hPValueBeta0_toyevt2","hPValueBeta0_toyevt2",200,0,1.5);
  TH1D *hPValueBeta0_toyevt3 = new TH1D("hPValueBeta0_toyevt3","hPValueBeta0_toyevt3",200,0,1.5);
  TH1D *hPValueBeta0_toyevt4 = new TH1D("hPValueBeta0_toyevt4","hPValueBeta0_toyevt4",200,0,1.5);
  TH1D *hPValueBeta0_toyevt5 = new TH1D("hPValueBeta0_toyevt5","hPValueBeta0_toyevt5",200,0,1.5);
  TH1D *hPValueBeta0_toyevt6 = new TH1D("hPValueBeta0_toyevt6","hPValueBeta0_toyevt6",200,0,1.5);
  TH1F *hPValueBeta0_toyevt7 = new TH1F("hPValueBeta0_toyevt7","hPValueBeta0_toyevt7",200,0,1.5);
  TH1F *hPValueBeta0_toyevt8 = new TH1F("hPValueBeta0_toyevt8","hPValueBeta0_toyevt8",200,0,1.5);
  TH1F *hPValueBeta0_toyevt9 = new TH1F("hPValueBeta0_toyevt9","hPValueBeta0_toyevt9",200,0,1.5);
  
  TH1F *hPValueBeta0_toyevtH = new TH1F("hPValueBeta0_toyevtH","hPValueBeta0_toyevtH",200,0,1.5);
  TH1F *hPValueBeta0_bkg = new TH1F("hPValueBeta0bkg","hPValueBeta0bkg",200,0,1.5);
  
  THStack *sPValueBeta0 = new THStack("sPValueBeta0","sPValueBeta0");
  
  TH1D *hPValueBeta1 = new TH1D("hPValueBeta1","hPValueBeta1",200,0,1.5);
  TH1D *hPValueBeta1_toyevt0 = new TH1D("hPValueBeta1_toyevt0","hPValueBeta1_toyevt0",200,0,1.5);
  TH1D *hPValueBeta1_toyevt1 = new TH1D("hPValueBeta1_toyevt1","hPValueBeta1_toyevt1",200,0,1.5);
  TH1D *hPValueBeta1_toyevt2 = new TH1D("hPValueBeta1_toyevt2","hPValueBeta1_toyevt2",200,0,1.5);
  TH1D *hPValueBeta1_toyevt3 = new TH1D("hPValueBeta1_toyevt3","hPValueBeta1_toyevt3",200,0,1.5);
  TH1D *hPValueBeta1_toyevt4 = new TH1D("hPValueBeta1_toyevt4","hPValueBeta1_toyevt4",200,0,1.5);
  TH1D *hPValueBeta1_toyevt5 = new TH1D("hPValueBeta1_toyevt5","hPValueBeta1_toyevt5",200,0,1.5);
  TH1D *hPValueBeta1_toyevt6 = new TH1D("hPValueBeta1_toyevt6","hPValueBeta1_toyevt6",200,0,1.5);
  TH1F *hPValueBeta1_toyevt7 = new TH1F("hPValueBeta1_toyevt7","hPValueBeta1_toyevt7",200,0,1.5);
  TH1F *hPValueBeta1_toyevt8 = new TH1F("hPValueBeta1_toyevt8","hPValueBeta1_toyevt8",200,0,1.5);
  TH1F *hPValueBeta1_toyevt9 = new TH1F("hPValueBeta1_toyevt9","hPValueBeta1_toyevt9",200,0,1.5);
  TH1F *hPValueBeta1_toyevtH = new TH1F("hPValueBeta1_toyevtH","hPValueBeta1_toyevtH",200,0,1.5);
  
  THStack *sPValueBeta1 = new THStack("sPValueBeta1","sPValueBeta1");
  
  
  // loop over the tree
  
  int Nentries = treein->GetEntries();
  cout << "Nfits: " << Nentries << endl;
  for(int ientry=0;ientry<Nentries;ientry++){
    
    treein->GetEntry(ientry);
    
    double nevents = (double)hErecoToy->Integral();
    cout << "Nevents: " << nevents << endl;
    if(nevents<0){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events = " << nevents << " < 0" << endl;
      cerr << endl;
      exit(1);
    }
    
    
    double pvalue_beta0 = hchisq1D_out->GetBinContent(1);
    
    hPValueBeta0->Fill(pvalue_beta0);
    
    if(nevents > exp_bkg){
      hPValueBeta0_bkg->Fill(pvalue_beta0);
    }
    
    if(nevents==0)      hPValueBeta0_toyevt0->Fill(pvalue_beta0);
    else if(nevents==1) hPValueBeta0_toyevt1->Fill(pvalue_beta0);
    else if(nevents==2) hPValueBeta0_toyevt2->Fill(pvalue_beta0);
    else if(nevents==3) hPValueBeta0_toyevt3->Fill(pvalue_beta0);
    else if(nevents==4) hPValueBeta0_toyevt4->Fill(pvalue_beta0);
    else if(nevents==5) hPValueBeta0_toyevt5->Fill(pvalue_beta0);
    else if(nevents==6) hPValueBeta0_toyevt6->Fill(pvalue_beta0);
    else if(nevents==7) hPValueBeta0_toyevt7->Fill(pvalue_beta0);
    else if(nevents==8) hPValueBeta0_toyevt8->Fill(pvalue_beta0);
    else if(nevents==9) hPValueBeta0_toyevt9->Fill(pvalue_beta0);
    else     hPValueBeta0_toyevtH->Fill(pvalue_beta0);
    
    double pvalue_beta1 = hchisq1D_out->GetBinContent(2);
    
    hPValueBeta1->Fill(pvalue_beta1);
    if(nevents==0)      hPValueBeta1_toyevt0->Fill(pvalue_beta1);
    else if(nevents==1) hPValueBeta1_toyevt1->Fill(pvalue_beta1);
    else if(nevents==2) hPValueBeta1_toyevt2->Fill(pvalue_beta1);
    else if(nevents==3) hPValueBeta1_toyevt3->Fill(pvalue_beta1);
    else if(nevents==4) hPValueBeta1_toyevt4->Fill(pvalue_beta1);
    else if(nevents==5) hPValueBeta1_toyevt5->Fill(pvalue_beta1);
    else if(nevents==6) hPValueBeta1_toyevt6->Fill(pvalue_beta1);
    else if(nevents==7) hPValueBeta1_toyevt7->Fill(pvalue_beta1);
    else if(nevents==8) hPValueBeta1_toyevt8->Fill(pvalue_beta1);
    else if(nevents==9) hPValueBeta1_toyevt9->Fill(pvalue_beta1);
    else     hPValueBeta1_toyevtH->Fill(pvalue_beta1);
    
  }
  
  
  // Draw histograms
  
  if(draw_predicted_event_distribution){
    treein->GetEntry(0); // Do it for 1st fit in tree
    double nevents = (double)hErecoToy->Integral();
    
    TCanvas *cBeta0Distribution = new TCanvas("cBeta0Distribution","cBeta0Distribution");
    hbeta0->GetYaxis()->SetTitle("Number of entries");
    hbeta0->GetXaxis()->SetTitle("Number of events");
    hbeta0->SetLineColor(kBlue);
    hbeta0->SetLineWidth(5);
    hbeta0->GetYaxis()->SetRangeUser(0.0002,hbeta0->GetMaximum()*10);
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
    TLine *expEvents = new TLine(nevents,0,nevents,hbeta0->GetBinContent(hbeta0->GetMaximumBin()));
    expEvents->SetLineWidth(4);
    expEvents->SetLineColor(kRed);
    expEvents->SetLineStyle(2);
    expEvents->Draw();
    TLegend *leg = new TLegend(0.54,0.71,0.89,0.91);
    leg->SetBorderSize(1);
    leg->SetLineWidth(3);
    leg->SetTextFont(132);
    leg->SetTextSize(0.05);
    leg->AddEntry(expEvents,"Number of events data","L");
    leg->AddEntry(hbeta0,"Number of events prediction","L");
    char name[200];
    sprintf(name,"P-Value: %1.4f",hPValueBeta0->GetMean());
    leg->AddEntry(expEvents,name,"l");
    leg->Draw();
    
    double beta0_error = get_error(nevents,hbeta0,hbeta0Unscaled);
    cout << "Beta0 pvalue: " << pvalue_beta0 << "\t   Error: " <<  beta0_error <<  endl;
    
    
    TCanvas *cBeta1Distribution = new TCanvas("cBeta1Distribution","cBeta1Distribution");
    hbeta1->GetYaxis()->SetTitle("Number of entries");
    hbeta1->GetXaxis()->SetTitle("Number of events");
    hbeta1->GetYaxis()->SetRangeUser(0.02,hbeta1->GetMaximum()*10);
    hbeta1->GetYaxis()->SetTitleSize(0.06);
    hbeta1->GetXaxis()->SetTitleSize(0.06);
    hbeta1->GetYaxis()->SetTitleFont(132);
    hbeta1->GetXaxis()->SetTitleFont(132);
    hbeta1->GetYaxis()->SetLabelFont(132);
    hbeta1->GetXaxis()->SetLabelFont(132);
    hbeta1->GetYaxis()->SetLabelSize(0.05);
    hbeta1->GetXaxis()->SetLabelSize(0.05);
    
    TLine *expEvents2 = new TLine(nevents,hbeta1->GetMinimum(),nevents,hbeta1->GetMaximum());
    expEvents2->SetLineWidth(4);
    expEvents2->SetLineColor(kRed);
    expEvents2->SetLineStyle(2);
    
    
    hbeta1->SetLineColor(kBlue);
    hbeta1->SetLineWidth(5);
    hbeta1->DrawCopy();
    cBeta1Distribution->SetLogy();
    expEvents2->Draw();
    TLegend *leg2 = new TLegend(0.54,0.71,0.89,0.91);
    leg2->SetBorderSize(1);
    leg2->SetLineWidth(3);
    leg2->SetTextFont(132);
    leg2->SetTextSize(0.05);
    leg2->AddEntry(expEvents,"Number of events data","L");
    leg2->AddEntry(hbeta0,"Number of events prediction","L");
    leg2->AddEntry(expEvents,Form("P-Value: %lf",hPValueBeta1->GetMean()),"l");
    leg2->Draw();
    
    double beta1_error = get_left_error(nevents,hbeta1,hbeta1Unscaled);
    cout << "Beta1 pvalue: " << pvalue_beta1 << "\t   Error: " <<  beta1_error <<  endl;
    
    TFile fout("nuebarapp_output.root","recreate");
    hbeta0->Write("hbeta0");
    hbeta0Unscaled->Write("hbeta0Unscaled");
    hbeta1->Write("hbeta1");
    hbeta1Unscaled->Write("hbeta1Unscaled");
    fout.Close();
  }
  
  
  
  TCanvas *cPValueBeta0 = new TCanvas("cPValueBeta0","cPValueBeta0");
  hPValueBeta0->GetYaxis()->SetTitle("Number of entries");
  hPValueBeta0->GetXaxis()->SetTitle("p-value");
  hPValueBeta0->SetLineWidth(2);
  hPValueBeta0->DrawCopy();
  hPValueBeta0_toyevt0->SetFillColor(kBlack);
  hPValueBeta0_toyevt0->SetLineColor(kBlack);
  hPValueBeta0_toyevt1->SetFillColor(kBlue);
  hPValueBeta0_toyevt1->SetLineColor(kBlue);
  hPValueBeta0_toyevt2->SetFillColor(kGreen+2);
  hPValueBeta0_toyevt2->SetLineColor(kGreen+2);
  hPValueBeta0_toyevt3->SetFillColor(kBlue+2);
  hPValueBeta0_toyevt3->SetLineColor(kBlue+2);
  hPValueBeta0_toyevt4->SetFillStyle(3001);
  hPValueBeta0_toyevt4->SetFillColor(kBlue);
  hPValueBeta0_toyevt4->SetLineColor(kBlue);
  hPValueBeta0_toyevt5->SetFillStyle(3001);
  hPValueBeta0_toyevt5->SetFillColor(kGreen+2);
  hPValueBeta0_toyevt5->SetLineColor(kGreen+2);
  hPValueBeta0_toyevt6->SetFillStyle(3001);
  hPValueBeta0_toyevt6->SetFillColor(kBlue+2);
  hPValueBeta0_toyevt6->SetLineColor(kBlue+2);
  
  hPValueBeta0_toyevt7->SetFillStyle(3013);
  hPValueBeta0_toyevt7->SetFillColor(kBlue);
  hPValueBeta0_toyevt7->SetLineColor(kBlue);
  hPValueBeta0_toyevt8->SetFillStyle(3013);
  hPValueBeta0_toyevt8->SetFillColor(kGreen+2);
  hPValueBeta0_toyevt8->SetLineColor(kGreen+2);
  hPValueBeta0_toyevt9->SetFillStyle(3013);
  hPValueBeta0_toyevt9->SetFillColor(kBlue+2);
  hPValueBeta0_toyevt9->SetLineColor(kBlue+2);
  
  sPValueBeta0->Add(hPValueBeta0_toyevt0);
  sPValueBeta0->Add(hPValueBeta0_toyevt1);
  sPValueBeta0->Add(hPValueBeta0_toyevt2);
  sPValueBeta0->Add(hPValueBeta0_toyevt3);
  sPValueBeta0->Add(hPValueBeta0_toyevt4);
  sPValueBeta0->Add(hPValueBeta0_toyevt5);
  sPValueBeta0->Add(hPValueBeta0_toyevt6);
  sPValueBeta0->Add(hPValueBeta0_toyevt7);
  sPValueBeta0->Add(hPValueBeta0_toyevt8);
  sPValueBeta0->Add(hPValueBeta0_toyevt9);
  
  sPValueBeta0->Draw("same");
  hPValueBeta0->DrawCopy("same");
  
  cPValueBeta0->SetLogy();
  
  std::stringstream mean0;
  std::stringstream rms0;
  mean0 << fixed << setprecision( 3 ) << hPValueBeta0->GetMean();
  rms0 << fixed << setprecision( 3 ) << hPValueBeta0->GetRMS();
  
  cout << "PValue0Mean: " << hPValueBeta0->GetMean() << endl;
  TLegend *leg = new TLegend(0.66,0.40,0.89,0.94);
  leg->SetBorderSize(1);
  leg->SetLineWidth(2);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);
  leg->AddEntry(hPValueBeta0,"Total","L");
  leg->AddEntry(hPValueBeta0,Form("Mean: %s",mean0.str().c_str()),"L");
  leg->AddEntry(hPValueBeta0,Form("RMS: %s",rms0.str().c_str()),"L");
  leg->AddEntry(hPValueBeta0_toyevt0,"Number of events = 0","F");
  leg->AddEntry(hPValueBeta0_toyevt1,"Number of events = 1","F");
  leg->AddEntry(hPValueBeta0_toyevt2,"Number of events = 2","F");
  leg->AddEntry(hPValueBeta0_toyevt3,"Number of events = 3","F");
  leg->AddEntry(hPValueBeta0_toyevt4,"Number of events = 4","F");
  leg->AddEntry(hPValueBeta0_toyevt5,"Number of events = 5","F");
  leg->AddEntry(hPValueBeta0_toyevt6,"Number of events = 6","F");
  leg->AddEntry(hPValueBeta0_toyevt7,"Number of events = 7","F");
  leg->AddEntry(hPValueBeta0_toyevt8,"Number of events = 8","F");
  leg->AddEntry(hPValueBeta0_toyevt9,"Number of events = 9","F");
  leg->Draw();
  
  TCanvas *cPValueBeta1 = new TCanvas("cPValueBeta1","cPValueBeta1");
  hPValueBeta1->GetYaxis()->SetTitle("Number of entries");
  hPValueBeta1->GetXaxis()->SetTitle("p-value");
  hPValueBeta1->SetLineWidth(2);
  hPValueBeta1->DrawCopy();
  hPValueBeta1_toyevt0->SetFillColor(kBlack);
  hPValueBeta1_toyevt0->SetLineColor(kBlack);
  hPValueBeta1_toyevt1->SetFillColor(kBlue);
  hPValueBeta1_toyevt1->SetLineColor(kBlue);
  hPValueBeta1_toyevt2->SetFillColor(kGreen+2);
  hPValueBeta1_toyevt2->SetLineColor(kGreen+2);
  hPValueBeta1_toyevt3->SetFillColor(kBlue+2);
  hPValueBeta1_toyevt3->SetLineColor(kBlue+2);
  hPValueBeta1_toyevt4->SetFillStyle(3001);
  hPValueBeta1_toyevt4->SetFillColor(kBlue);
  hPValueBeta1_toyevt4->SetLineColor(kBlue);
  hPValueBeta1_toyevt5->SetFillStyle(3001);
  hPValueBeta1_toyevt5->SetFillColor(kGreen+2);
  hPValueBeta1_toyevt5->SetLineColor(kGreen+2);
  hPValueBeta1_toyevt6->SetFillStyle(3001);
  hPValueBeta1_toyevt6->SetFillColor(kBlue+2);
  hPValueBeta1_toyevt6->SetLineColor(kBlue+2);
  
  hPValueBeta1_toyevt7->SetFillStyle(3013);
  hPValueBeta1_toyevt7->SetFillColor(kBlue);
  hPValueBeta1_toyevt7->SetLineColor(kBlue);
  hPValueBeta1_toyevt8->SetFillStyle(3013);
  hPValueBeta1_toyevt8->SetFillColor(kGreen+2);
  hPValueBeta1_toyevt8->SetLineColor(kGreen+2);
  hPValueBeta1_toyevt9->SetFillStyle(3013);
  hPValueBeta1_toyevt9->SetFillColor(kBlue+2);
  hPValueBeta1_toyevt9->SetLineColor(kBlue+2);
  
  sPValueBeta1->Add(hPValueBeta1_toyevt0);
  sPValueBeta1->Add(hPValueBeta1_toyevt1);
  sPValueBeta1->Add(hPValueBeta1_toyevt2);
  sPValueBeta1->Add(hPValueBeta1_toyevt3);
  sPValueBeta1->Add(hPValueBeta1_toyevt4);
  sPValueBeta1->Add(hPValueBeta1_toyevt5);
  sPValueBeta1->Add(hPValueBeta1_toyevt6);
  sPValueBeta1->Add(hPValueBeta1_toyevt7);
  sPValueBeta1->Add(hPValueBeta1_toyevt8);
  sPValueBeta1->Add(hPValueBeta1_toyevt9);
  
  
  sPValueBeta1->Draw("same");
  hPValueBeta1->DrawCopy("same");
  cPValueBeta1->SetLogy();
  leg->Draw();
  
  std::stringstream mean1;
  std::stringstream rms1;
  mean1 << fixed << setprecision( 3 ) << hPValueBeta1->GetMean();
  rms1 << fixed << setprecision( 3 ) << hPValueBeta1->GetRMS();
  
  
  TLegend *leg2 = new TLegend(0.66,0.40,0.89,0.94);
  leg2->SetBorderSize(1);
  leg2->SetLineWidth(2);
  leg2->SetTextFont(132);
  leg2->SetTextSize(0.05);
  leg2->AddEntry(hPValueBeta1,"Total","L");
  leg2->AddEntry(hPValueBeta1,Form("Mean: %s",mean1.str().c_str()),"L");
  leg2->AddEntry(hPValueBeta1,Form("RMS: %s",rms1.str().c_str()),"L");
  leg2->AddEntry(hPValueBeta1_toyevt0,"Number of events = 0","F");
  leg2->AddEntry(hPValueBeta1_toyevt1,"Number of events = 1","F");
  leg2->AddEntry(hPValueBeta1_toyevt2,"Number of events = 2","F");
  leg2->AddEntry(hPValueBeta1_toyevt3,"Number of events = 3","F");
  leg2->AddEntry(hPValueBeta1_toyevt4,"Number of events = 4","F");
  leg2->AddEntry(hPValueBeta1_toyevt5,"Number of events = 5","F");
  leg2->AddEntry(hPValueBeta1_toyevt6,"Number of events = 6","F");
  leg2->AddEntry(hPValueBeta1_toyevt7,"Number of events = 7","F");
  leg2->AddEntry(hPValueBeta1_toyevt8,"Number of events = 8","F");
  leg2->AddEntry(hPValueBeta1_toyevt9,"Number of events = 9","F");
  
  leg2->Draw();
  
  cout << "\n\n Average p-values! \n" <<endl;
  cout << "Beta=0" << endl;
  cout << "0  events: " <<  hPValueBeta0_toyevt0->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt0->Integral() <<  endl;
  cout << "1  event : " <<  hPValueBeta0_toyevt1->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt1->Integral() << endl;
  cout << "2  events: " <<  hPValueBeta0_toyevt2->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt2->Integral() <<  endl;
  cout << "3  events: " <<  hPValueBeta0_toyevt3->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt3->Integral() <<  endl;
  cout << "4  events: " <<  hPValueBeta0_toyevt4->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt4->Integral() <<  endl;
  cout << "5  events: " <<  hPValueBeta0_toyevt5->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt5->Integral() <<  endl;
  cout << "6  events: " <<  hPValueBeta0_toyevt6->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt6->Integral() <<  endl;
  cout << "7  events: " <<  hPValueBeta0_toyevt7->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt7->Integral() <<  endl;
  cout << "8  events: " <<  hPValueBeta0_toyevt8->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt8->Integral() <<  endl;
  cout << "9  events: " <<  hPValueBeta0_toyevt9->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt9->Integral() <<  endl;
  cout << ">9 events: " <<  hPValueBeta0_toyevtH->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevtH->Integral() <<  endl;
  
  cout << "Beta=1" << endl;
  cout << "0  events: " <<  hPValueBeta1_toyevt0->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt0->Integral() <<  endl;
  cout << "1  event : " <<  hPValueBeta1_toyevt1->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt1->Integral() <<  endl;
  cout << "2  events: " <<  hPValueBeta1_toyevt2->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt2->Integral() <<  endl;
  cout << "3  events: " <<  hPValueBeta1_toyevt3->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt3->Integral() <<  endl;
  cout << "4  events: " <<  hPValueBeta1_toyevt4->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt4->Integral() <<  endl;
  cout << "5  events: " <<  hPValueBeta1_toyevt5->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt5->Integral() <<  endl;
  cout << "6  events: " <<  hPValueBeta1_toyevt6->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt6->Integral() <<  endl;
  cout << "7  events: " <<  hPValueBeta1_toyevt7->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt7->Integral() <<  endl;
  cout << "8  events: " <<  hPValueBeta1_toyevt8->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt8->Integral() <<  endl;
  cout << "9  events: " <<  hPValueBeta1_toyevt9->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevt9->Integral() <<  endl;
  cout << ">9 events: " <<  hPValueBeta1_toyevtH->GetMean() << "\t\t Toys: " << hPValueBeta1_toyevtH->Integral() <<  endl;
  
  
  cout << "NData > NBackground: " << hPValueBeta0_bkg->GetEntries() / Nentries << " Mean pvalue: " << hPValueBeta0_bkg->GetMean() << endl;
  cout << "Fraction below pvalue=0.1: " << hPValueBeta0->Integral(0,14)/ Nentries  << endl;
  cout << "Fraction below pvalue=0.05: " << hPValueBeta0->Integral(0,7)/ Nentries << endl;
  
  return;
}

double get_error(double obs, TH1D* hbeta, TH1D* hbetaUnscaled)
{
  int obs_bin = hbeta->FindBin(obs);
  double pvalue = hbeta->Integral(obs_bin,hbeta->GetXaxis()->GetNbins()) /  hbeta->Integral();
  double err_sqrd = 0;
  //  for(int i=1;i<obs_bin;i++)
  //  {
  //    double pd = hbeta->Integral();
  //    double wl = hbeta->GetBinContent(i)/hbetaUnscaled->GetBinContent(i);
  //    double bl = hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wl* (pvalue) / pd, 2) * bl;
  //    err_sqrd = err_sqrd + err;
  //  }
  
  //Left error
  double pd = hbeta->Integral();
  double wl = hbeta->Integral(1,obs_bin-1)/hbetaUnscaled->Integral(1,obs_bin-1);
  double bl = hbetaUnscaled->Integral(1,obs_bin-1);
  double err = TMath::Power(wl* (pvalue) / pd, 2) * bl;
  err_sqrd = err_sqrd + err;
  
  //Right error
  double pd = hbeta->Integral();
  double wr = hbeta->Integral(obs_bin,hbeta->GetXaxis()->GetNbins())/hbetaUnscaled->Integral(obs_bin,hbeta->GetXaxis()->GetNbins());
  double br = hbetaUnscaled->Integral(obs_bin,hbeta->GetXaxis()->GetNbins());
  double err = TMath::Power(wr* (1 - pvalue) / pd, 2) * br;
  err_sqrd = err_sqrd + err;
  
  //  for(int i=obs_bin;i<=hbeta->GetXaxis()->GetNbins();i++)
  //  {
  //    double pd = hbeta->Integral();
  //    double wr = hbeta->GetBinContent(i)/hbetaUnscaled->GetBinContent(i);
  //    double br = hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wr* (1 - pvalue) / pd, 2) * br;
  //    err_sqrd = err_sqrd + err;
  //  }
  double final_error = TMath::Sqrt(err_sqrd);
  return final_error;
}


double get_left_error(double obs, TH1D* hbeta, TH1D* hbetaUnscaled)
{
  int obs_bin = hbeta->FindBin(obs);
  double pvalue = hbeta->Integral(1,obs_bin) / hbeta->Integral() ;
  double err_sqrd = 0;
  //  for(int i=1;i<obs_bin;i++)
  //  {
  //    double pd = hbeta->Integral();
  //    double wl = hbeta->GetBinContent(i)/hbetaUnscaled->GetBinContent(i);
  //    double bl = hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wl* (pvalue) / pd, 2) * bl;
  //    err_sqrd = err_sqrd + err;
  //  }
  
  //Right error
  double pd = hbeta->Integral();
  double wr = hbeta->Integral(obs_bin+1,hbeta->GetXaxis()->GetNbins())/hbetaUnscaled->Integral(obs_bin+1,hbeta->GetXaxis()->GetNbins());
  double br = hbetaUnscaled->Integral(obs_bin +1,hbeta->GetXaxis()->GetNbins());
  double err = TMath::Power(wr* (pvalue) / pd, 2) * br;
  err_sqrd = err_sqrd + err;
  
  //Left error
  double pd = hbeta->Integral();
  double wl = hbeta->Integral(1,obs_bin)/hbetaUnscaled->Integral(1,obs_bin);
  double bl = hbetaUnscaled->Integral(1,obs_bin);
  double err = TMath::Power(wl* (1 - pvalue) / pd, 2) * bl;
  err_sqrd = err_sqrd + err;
  
  //  for(int i=obs_bin;i<=hbeta->GetXaxis()->GetNbins();i++)
  //  {
  //    double pd = hbeta->Integral();
  //    double wr = hbeta->GetBinContent(i)/hbetaUnscaled->GetBinContent(i);
  //    double br = hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wr* (1 - pvalue) / pd, 2) * br;
  //    err_sqrd = err_sqrd + err;
  //  }
  double final_error = TMath::Sqrt(err_sqrd);
  return final_error;
}

double get_error_vector( std::vector<double> _right_tail_weight, std::vector<double> _left_tail_weight, double _pvalue)
{
  double error = 0;
  double total_right=0;
  double total_left=0;
  for(int i=0;i<_right_tail_weight.size();i++)
  {
    error+= TMath::Power(_right_tail_weight[i]*(1-_pvalue),2);
    total_right+= _right_tail_weight[i];
  }
  
  for(int i=0;i<_left_tail_weight.size();i++)
  {
    error+= TMath::Power(_left_tail_weight[i]*(_pvalue),2);
    total_left+= _left_tail_weight[i];
  }
  
  error = error / (TMath::Power(total_left + total_right ,2));
  return TMath::Sqrt(error);
}

