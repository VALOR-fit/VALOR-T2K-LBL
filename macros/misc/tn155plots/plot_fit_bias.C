//make the fit residual plots for mixing23 & dm32sq, where 2 different POT are overlaid
void plot_fit_bias
(
 const char * file1 = "valor_fit_bias_study_nuis_sys_0.002400_0.370.root",
 const char * file2 = "valor_fit_bias_study_nuis_sys_0.002400_0.370_goalpot.root",
 const char * tag   = "0.00240_0.370",
 double mixlo = -0.1,
 double mixhi = +0.4,
 double dmlo  = -0.0005,
 double dmhi  = +0.0005
 )
{
  Color_t color[2] = {kBlack,kRed};
  gStyle->SetOptStat(0);

  TFile f1(file1);
  //get histograms
  TH1D * hm1 = 0;
  f1.GetObject("h_residuals_sinsqtheta23",hm1);
  TH1D * hd1 = 0;
  f1.GetObject("h_residuals_dm32sq",hd1);
  //scale
  hm1->Scale(1 / hm1->Integral());
  hd1->Scale(1 / hd1->Integral());
  //style
  hm1->SetLineColor(color[0]);
  hm1->SetLineWidth(1);
  hd1->SetLineColor(color[0]);
  hd1->SetLineWidth(1);
  //for "stat box"
  double hm1mean = hm1->GetMean();
  double hm1rms  = hm1->GetRMS();
  double hd1mean = hd1->GetMean();
  double hd1rms  = hd1->GetRMS();
  //range
  hm1->GetXaxis()->SetRangeUser(mixlo,mixhi);
  hd1->GetXaxis()->SetRangeUser(dmlo, dmhi);
  //title
  hm1->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
  hd1->GetXaxis()->SetTitle("|#Deltam^{2}_{32}|");
  hm1->GetYaxis()->SetTitle("Fraction of fits");
  hd1->GetYaxis()->SetTitle("Fraction of fits");
  hm1->GetYaxis()->SetTitleOffset(0.55);
  hd1->GetYaxis()->SetTitleOffset(0.55);
  //"stat box"
  TPaveText ptm1(0.85,0.8,1,1,"br NDC");
  ptm1.SetTextSize(0.06);
  ptm1.SetShadowColor(0);
  ptm1.SetFillColor(0);
  //ptm1.SetFillStyle(0);
  ptm1.AddText(Form("Mean \t %.4f", hm1mean));
  ptm1.AddText(Form("RMS  \t %.4f", hm1rms));
  ptm1.SetTextColor(color[0]);
  TPaveText ptd1(0.85,0.8,1,1,"br NDC");
  ptd1.SetTextSize(0.06);
  ptd1.SetShadowColor(0);
  ptd1.SetFillColor(0);
  //ptd1.SetFillStyle(0);
  ptd1.AddText(Form("Mean \t %.4g", hd1mean));
  ptd1.AddText(Form("RMS  \t %.4g", hd1rms));
  ptd1.SetTextColor(color[0]);
  //kill stats box
  hm1->SetStats(0);
  hd1->SetStats(0);

  TFile f2(file2);
  //get histograms
  TH1D * hm2 = 0;
  f2.GetObject("h_residuals_sinsqtheta23",hm2);
  TH1D * hd2 = 0;
  f2.GetObject("h_residuals_dm32sq",hd2);
  //scale
  hm2->Scale(1 / hm2->Integral());
  hd2->Scale(1 / hd2->Integral());
  //style
  hm2->SetLineColor(color[1]);
  hm2->SetLineWidth(1);
  hd2->SetLineColor(color[1]);
  hd2->SetLineWidth(1);
  //for "stat box"
  double hm2mean = hm2->GetMean();
  double hm2rms  = hm2->GetRMS();
  double hd2mean = hd2->GetMean();
  double hd2rms  = hd2->GetRMS();
  //range
  hm2->GetXaxis()->SetRangeUser(mixlo,mixhi);
  hd2->GetXaxis()->SetRangeUser(dmlo, dmhi);
  //title
  hm2->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
  hd2->GetXaxis()->SetTitle("|#Deltam^{2}_{32}|");
  hm2->GetYaxis()->SetTitle("Fraction of fits");
  hd2->GetYaxis()->SetTitle("Fraction of fits");
  hm2->GetYaxis()->SetTitleOffset(0.55);
  hd2->GetYaxis()->SetTitleOffset(0.55);
  //"stat box"
  TPaveText ptm2(0.85,0.8,1,0.6,"br NDC");
  ptm2.SetTextSize(0.06);
  ptm2.SetShadowColor(0);
  ptm2.SetFillColor(0);
  //ptm2.SetFillStyle(0);
  ptm2.AddText(Form("Mean \t %.4f", hm2mean));
  ptm2.AddText(Form("RMS  \t %.4f", hm2rms));
  ptm2.SetTextColor(color[1]);
  TPaveText ptd2(0.85,0.8,1,0.6,"br NDC");
  ptd2.SetTextSize(0.06);
  ptd2.SetShadowColor(0);
  ptd2.SetFillColor(0);
  //ptd2.SetFillStyle(0);
  ptd2.AddText(Form("Mean \t %.4g", hd2mean));
  ptd2.AddText(Form("RMS  \t %.4g", hd2rms));
  ptd2.SetTextColor(color[1]);
  //kill stats box
  hm2->SetStats(0);
  hd2->SetStats(0);


  //draw & save
  TCanvas * c = new TCanvas();
  c->Divide(1,2);
  c->GetPad(1)->SetLeftMargin(0.07);
  c->GetPad(2)->SetLeftMargin(0.07);
  c->cd(1);
  if(hm2->GetMaximum() > hm1->GetMaximum()) {
    hm2->Draw();
    hm1->Draw("same");
  }
  else {
    hm1->Draw();
    hm2->Draw("same");
  }
  ptm1.Draw();
  ptm2.Draw();
  c->cd(2);
  if(hd2->GetMaximum() > hd1->GetMaximum()) {
    hd2->Draw();
    hd1->Draw("same");
  }
  else {
    hd1->Draw();
    hd2->Draw("same");
  }
  ptd1.Draw();
  ptd2.Draw();
  c->SaveAs(Form("residuals_oscparams_%s.eps",tag));
  c->SaveAs(Form("residuals_oscparams_%s.pdf",tag));

}
