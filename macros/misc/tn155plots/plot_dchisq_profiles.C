void plot_dchisq_profiles
(
 const char * filename1 = "2013results/DATA/PROFILES/run1-4_nuis3fixdcp_sintheta_normal/mix/run1-4_nuis3fixdcp_sintheta_normal_mixprofile.root",
 const char * filename2 = "2013results/DATA/PROFILES/run1-4_nuis3fixdcp_sintheta_inverted/mix/run1-4_nuis3fixdcp_sintheta_inverted_mixprofile.root",
 const char * tag1      = "Normal mass hierarchy",
 const char * tag2      = "Inverted mass hierarchy",
 const char * fchisqbf1 = "2013results/DATA/run1-4_nuis3fixdcp_sintheta_normal/combination/run1-4_nuis3fixdcp_sintheta_normal-run_0-expt_6.root",
 const char * fchisqbf2 = "2013results/DATA/run1-4_nuis3fixdcp_sintheta_inverted/combination/run1-4_nuis3fixdcp_sintheta_inverted-run_0-expt_6.root",
 bool         plot_constdchisq = true,
 bool         is_mix23  = true
 )
{
  using namespace t2kana;

  // get chisq profiles
  TFile f1(filename1);
  TFile f2(filename2);
  TH1D * h1 = 0, * h2 = 0;
  if(is_mix23) {
    f1.GetObject("chisq_mixing23",h1);
    f2.GetObject("chisq_mixing23",h2);
  }
  else {
    f1.GetObject("chisq_dm32sq",h1);
    f2.GetObject("chisq_dm32sq",h2);
  }

  //subtract chisq_bf
  TFile fbf1(fchisqbf1);
  TFile fbf2(fchisqbf2);
  FitParamList * fparm1 = 0, * fparm2 = 0;
  fbf1.GetObject("fparm",fparm1);
  fbf2.GetObject("fparm",fparm2);
  double bf1 = fparm1->GetBestFitChisq();
  double bf2 = fparm2->GetBestFitChisq();
  for(int ix = 1; ix <= h1->GetNbinsX(); ix++) {
    double dchisq = h1->GetBinContent(ix) - bf1;
    if(dchisq < 0)
      cout << h1->GetBinCenter(ix) << " has a negative dchisq " << dchisq << " for histogram " << tag1 << endl;
    dchisq = TMath::Max(0., dchisq);
    h1->SetBinContent(ix, dchisq);
  }
  for(int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    double dchisq = TMath::Max(0., h2->GetBinContent(ix) - bf2);
    if(dchisq < 0)
      cout << h2->GetBinCenter(ix) << " has a negative dchisq " << dchisq << " for histogram " << tag2 << endl;
    h2->SetBinContent(ix, dchisq);
  }

  //format
  h2->SetLineColor(kRed);
  h2->SetLineStyle(2);
  h1->GetYaxis()->SetTitle("#Delta#chi^{2}");
  h2->GetYaxis()->SetTitle("#Delta#chi^{2}");
  h1->GetYaxis()->SetRangeUser(0,h1->GetMaximum()*1.05);
  if(!is_mix23) { // hide the dodgy peaks at ~ 0.00295)
    h1->GetXaxis()->SetRangeUser(0.0020,0.0029);
    h1->SetNdivisions(509);
  }

  //constdchisq lines
  double c68 = 1.00, c90 = 2.31;
  TLine * l68 = new TLine(h1->GetXaxis()->GetXmin(), c68, h1->GetXaxis()->GetXmax(), c68);
  l68->SetLineStyle(3);
  TLine * l90 = new TLine(h1->GetXaxis()->GetXmin(), c90, h1->GetXaxis()->GetXmax(), c90);
  l90->SetLineStyle(3);

  //legend
  TLegend * leg = new TLegend(0.3,0.6,0.7,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->AddEntry(h1,tag1,"l");
  leg->AddEntry(h2,tag2,"l");

  //draw
  TCanvas * c = new TCanvas();
  h1->Draw();
  h2->Draw("same");
  if(plot_constdchisq) {
    l68->Draw();
    l90->Draw();
  }
  leg->Draw();

  //save
  if(is_mix23) {
    c->SaveAs("dchisq_profiles_mix23.pdf");
    c->SaveAs("dchisq_profiles_mix23.eps");
  }
  else {
    c->SaveAs("dchisq_profiles_dm32sq.pdf");
    c->SaveAs("dchisq_profiles_dm32sq.eps");
  }

  //find the limits
  double h1_68_lo = -999;
  for(int ix = 1; ix <= h1->GetNbinsX(); ix++) {
    if(h1->GetBinContent(ix) <= c68)
      break;
    h1_68_lo = h1->GetBinCenter(ix);
  }
  double h1_68_hi = -999;
  for(int ix = h1->GetNbinsX(); ix >= 1; ix--) {
    if(h1->GetBinContent(ix) <= c68)
      break;
    h1_68_hi = h1->GetBinCenter(ix);
  }
  double h1_90_lo = -999;
  for(int ix = 1; ix <= h1->GetNbinsX(); ix++) {
    if(h1->GetBinContent(ix) <= c90)
      break;
    h1_90_lo = h1->GetBinCenter(ix);
  }
  double h1_90_hi = -999;
  for(int ix = h1->GetNbinsX(); ix >= 1; ix--) {
    if(h1->GetBinContent(ix) <= c90)
      break;
    h1_90_hi = h1->GetBinCenter(ix);
  }
  double h2_68_lo = -999;
  for(int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    if(h2->GetBinContent(ix) <= c68)
      break;
    h2_68_lo = h2->GetBinCenter(ix);
  }
  double h2_68_hi = -999;
  for(int ix = h2->GetNbinsX(); ix >= 1; ix--) {
    if(h2->GetBinContent(ix) <= c68)
      break;
    h2_68_hi = h2->GetBinCenter(ix);
  }
  double h2_90_lo = -999;
  for(int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    if(h2->GetBinContent(ix) <= c90)
      break;
    h2_90_lo = h2->GetBinCenter(ix);
  }
  double h2_90_hi = -999;
  for(int ix = h2->GetNbinsX(); ix >= 1; ix--) {
    if(h2->GetBinContent(ix) <= c90)
      break;
    h2_90_hi = h2->GetBinCenter(ix);
  }
  cout << endl << endl << "LIMITS" << endl
       << tag1 << " 68%:" << endl
       << h1_68_lo << "\t"
       << h1_68_hi << endl
       << tag1 << " 90%:" << endl
       << h1_90_lo << "\t"
       << h1_90_hi << endl
       << tag2 << " 68%:" << endl
       << h2_68_lo << "\t"
       << h2_68_hi << endl
       << tag2 << " 90%:" << endl
       << h2_90_lo << "\t"
       << h2_90_hi << endl
       << endl;
}

