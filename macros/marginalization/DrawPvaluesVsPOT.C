
void DrawPvaluesVsPOT
(
)
{
  char name[200];
  
  const int npts = 5;
  double npot[npts];
  double dchisq_asimov[npts];
  double nobs_asimov[npts];
  double pvalue_rate_asimov[npts];
  double pvalue_shape_asimov[npts];
  double pvalue_rate_mean[npts];
  double pvalue_shape_mean[npts];

  // 4E20POT
  npot[0] = 4.0; 
  dchisq_asimov[0] = 1.72185;
  nobs_asimov[0] = 3.71812;
  pvalue_rate_asimov[0] = 0.113;
  pvalue_shape_asimov[0] = 0.0723;
  pvalue_rate_mean[0] = 0.20899822;
  pvalue_shape_mean[0] = 0.14327296;

  // 4.5E20POT
  npot[1] = 4.5; 
  dchisq_asimov[1] = 1.93242;
  nobs_asimov[1] = 4.18289;
  pvalue_rate_asimov[1] = 0.0538;
  pvalue_shape_asimov[1] = 0.0572;
  pvalue_rate_mean[1] = 0.18764391;
  pvalue_shape_mean[1] = 0.12671067;

  // 6.5E20POT
  npot[2] = 6.5; 
  dchisq_asimov[2] = 2.76351;
  nobs_asimov[2] = 6.04195;
  pvalue_rate_asimov[2] = 0.0319;
  pvalue_shape_asimov[2] = 0.0359;
  pvalue_rate_mean[2] = 0.12947445;
  pvalue_shape_mean[2] = 0.0863986;
  
  // 8E20POT
  npot[3] = 8.0; 
  dchisq_asimov[3] = 3.37486;
  nobs_asimov[3] = 7.43625;
  pvalue_rate_asimov[3] = 0.032;
  pvalue_shape_asimov[3] = 0.0237;
  pvalue_rate_mean[3] = 0.09954281;
  pvalue_shape_mean[3] = 0.066444;

  // 9E20POT
  npot[4] = 9.0; 
  dchisq_asimov[4] = 3.77671;
  nobs_asimov[4] = 8.36578;
  pvalue_rate_asimov[4] = 0.0239;
  pvalue_shape_asimov[4] = 0.0184;
  pvalue_rate_mean[4] = 0.08880659;
  pvalue_shape_mean[4] = 0.057858;
  
  //

  TGraph *gPvalue_nobs_asimov_Vs_pot = new TGraph(npts,npot,nobs_asimov);
  TGraph *gPvalue_rate_asimov_Vs_pot = new TGraph(npts,npot,pvalue_rate_asimov);
  TGraph *gPvalue_shape_asimov_Vs_pot = new TGraph(npts,npot,pvalue_shape_asimov);
  TGraph *gPvalue_rate_mean_Vs_pot = new TGraph(npts,npot,pvalue_rate_mean);
  TGraph *gPvalue_shape_mean_Vs_pot = new TGraph(npts,npot,pvalue_shape_mean);

  TCanvas *cnobs_Vs_pot_asimov = new TCanvas("cnobs_Vs_pot_asimov","cnobs_Vs_pot_asimov");
  gPvalue_nobs_asimov_Vs_pot->SetMinimum(0.);
  gPvalue_nobs_asimov_Vs_pot->SetMaximum(10.);
  gPvalue_nobs_asimov_Vs_pot->GetXaxis()->SetTitle("# of POT (#times 10^{20})");
  gPvalue_nobs_asimov_Vs_pot->GetYaxis()->SetTitle("# of expected events");
  gPvalue_nobs_asimov_Vs_pot->SetLineWidth(2);
  gPvalue_nobs_asimov_Vs_pot->Draw("apl");

  //

  TCanvas *crate_Vs_pot = new TCanvas("crate_Vs_pot","crate_Vs_pot");
  gPvalue_rate_mean_Vs_pot->SetMinimum(0.);
  gPvalue_rate_mean_Vs_pot->GetXaxis()->SetTitle("# of POT (#times 10^{20})");
  gPvalue_rate_mean_Vs_pot->GetYaxis()->SetTitle("p-value (rate-only)");
  gPvalue_rate_asimov_Vs_pot->SetLineColor(kBlue);
  gPvalue_rate_mean_Vs_pot->SetLineColor(kRed);
  gPvalue_rate_asimov_Vs_pot->SetLineWidth(2);
  gPvalue_rate_mean_Vs_pot->SetLineWidth(2);
  gPvalue_rate_mean_Vs_pot->Draw("apl");
  gPvalue_rate_asimov_Vs_pot->Draw("pl same");

  TLegend *lPvalue_rate_Vs_pot = new TLegend(0.6,0.6,0.9,0.9);
  lPvalue_rate_Vs_pot->SetFillColor(0);
  lPvalue_rate_Vs_pot->SetLineColor(0);
  lPvalue_rate_Vs_pot->SetTextFont(132);
  lPvalue_rate_Vs_pot->SetTextSize(0.06);
  lPvalue_rate_Vs_pot->AddEntry(gPvalue_rate_asimov_Vs_pot,"Asimov","l");
  lPvalue_rate_Vs_pot->AddEntry(gPvalue_rate_mean_Vs_pot,"Mean","l");
  lPvalue_rate_Vs_pot->Draw();

  //
  
  TCanvas *crate_Vs_pot_mean = new TCanvas("crate_Vs_pot_mean","crate_Vs_pot_mean");
  gPvalue_shape_mean_Vs_pot->SetMinimum(0.);
  gPvalue_shape_mean_Vs_pot->GetXaxis()->SetTitle("# of POT (#times 10^{20})");
  gPvalue_shape_mean_Vs_pot->GetYaxis()->SetTitle("p-value (rate+shape)");
  gPvalue_shape_asimov_Vs_pot->SetLineColor(kBlue);
  gPvalue_shape_mean_Vs_pot->SetLineColor(kRed);
  gPvalue_shape_asimov_Vs_pot->SetLineWidth(2);
  gPvalue_shape_mean_Vs_pot->SetLineWidth(2);
  gPvalue_shape_mean_Vs_pot->Draw("apl");
  gPvalue_shape_asimov_Vs_pot->Draw("pl same");

  TLegend *lPvalue_shape_Vs_pot = new TLegend(0.6,0.6,0.9,0.9);
  lPvalue_shape_Vs_pot->SetFillColor(0);
  lPvalue_shape_Vs_pot->SetLineColor(0);
  lPvalue_shape_Vs_pot->SetTextFont(132);
  lPvalue_shape_Vs_pot->SetTextSize(0.06);
  lPvalue_shape_Vs_pot->AddEntry(gPvalue_shape_asimov_Vs_pot,"Asimov","l");
  lPvalue_shape_Vs_pot->AddEntry(gPvalue_shape_mean_Vs_pot,"Mean","l");
  lPvalue_shape_Vs_pot->Draw();

  
}
