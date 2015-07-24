void analyseFC
(
 double max_mixing = 1.0,
 double min_mixing = 0.9975,
 double max_dm     = 0.00245,
 double min_dm     = 0.00235
)
{
  TChain c("fcnt");
  c.Add("rootfiles/FC28*root");
  cout << c.GetEntries() << " entries" << endl;



  TString cut = Form("input_mixing23>%f&&input_mixing23<%f&&input_dm23sq>%f &&input_dm23sq<%f",min_mixing,max_mixing,min_dm,max_dm);
  cout << cut << endl;

  //c.Draw("input_mixing23:input_dm23sq","status!=3","colz");
  //c.Draw("input_mixing23:input_dm23sq","chisq_for_input_params-chisq<0","colz");


  TCanvas * can = new TCanvas();
  can->Divide(2,2);

  can->cd(1);
  c.Draw("chisq>>hbf", cut.Data());
  can->cd(2);
  c.Draw("chisq_for_input_params>>htr", cut.Data());

  TH1F *hbf = (TH1F*)gDirectory->Get("hbf");
  TH1F *htr = (TH1F*)gDirectory->Get("htr");
  TH1F *hbfcu = new TH1F(*hbf);
  TH1F *htrcu = new TH1F(*htr);

  hbfcu->ComputeIntegral();
  Double_t *integral = hbfcu->GetIntegral();
  hbfcu->SetContent(integral);
  can->cd(3);
  hbfcu->Draw();

  htrcu->ComputeIntegral();
  Double_t *integral = htrcu->GetIntegral();
  htrcu->SetContent(integral);
  can->cd(4);
  htrcu->Draw();

  can->SaveAs(Form("%f.pdf",max_mixing));


  
  TCanvas * c3 = new TCanvas();
  c3->Divide(1,2);

  c3->cd(1);
  c.Draw("mixing23", cut.Data());

  c3->cd(2);
  c.Draw("dm23sq", cut.Data());  

  c3->SaveAs(Form("bf%f.pdf",max_mixing));




  TCanvas * c2 = new TCanvas();
  c2->Divide(1,2);

  c2->cd(1);
  c.Draw("chisq_for_input_params-chisq>>hdc(10000,-0.1,50)", cut.Data());

  TH1F *hdc = (TH1F*)gDirectory->Get("hdc");
  TH1F *hdccu = new TH1F(*hdc);

  hdccu->ComputeIntegral();
  Double_t *integral = hdccu->GetIntegral();
  hdccu->SetContent(integral);
  c2->cd(2);
  hdccu->Draw();
 
  c2->SaveAs(Form("d%f.pdf",max_mixing));

 

  cout << endl << hdc->GetEntries() << endl;
  cout << hdccu->GetXaxis()->GetBinLowEdge(hdccu->FindFirstBinAbove(0.68,1)) << "\t";
  cout << hdccu->GetXaxis()->GetBinCenter(hdccu->FindFirstBinAbove(0.68,1)) << "\t";
  cout << hdccu->GetXaxis()->GetBinUpEdge(hdccu->FindFirstBinAbove(0.68,1)) << endl;



  TFile fc("delta_chisq_critical.root");
  TH2D * hfc = 0;
  fc.GetObject("hDeltaChi2Surface68",hfc);
  int ibin = hfc->GetXaxis()->FindBin((max_mixing+min_mixing)/2);
  int jbin = hfc->GetYaxis()->FindBin((max_dm+min_dm)/2);
  cout << hfc->GetBinContent(ibin,jbin) << endl;

  

  cut += "&&status==3&&chisq_for_input_params-chisq>=0";
  c.Draw("chisq_for_input_params-chisq>>hdc2(10000,-0.1,50)", cut.Data());

  TH1F *hdc2 = (TH1F*)gDirectory->Get("hdc2");
  TH1F *hdccu2 = new TH1F(*hdc2);

  hdccu2->ComputeIntegral();
  Double_t *integral = hdccu2->GetIntegral();
  hdccu2->SetContent(integral);
  cout << endl << hdc2->GetEntries() << endl;
  cout << hdccu2->GetXaxis()->GetBinLowEdge(hdccu2->FindFirstBinAbove(0.68,1)) << "\t";
  cout << hdccu2->GetXaxis()->GetBinCenter(hdccu2->FindFirstBinAbove(0.68,1)) << "\t";
  cout << hdccu2->GetXaxis()->GetBinUpEdge(hdccu2->FindFirstBinAbove(0.68,1)) << endl;


  cout << endl;
}
