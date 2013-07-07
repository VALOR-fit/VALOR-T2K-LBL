{

  TFile f0("$T2KOSCANA/2012results/NEWRESULTS/datafit/dchisq_2all.root");
  TFile f1("$T2KOSCANA/2012results/1DFC/dm23sq/1DFC_dm23sq-run_0-expts_3_3.root");
  TFile f2("$T2KOSCANA/2012results/1DFC/octant2/dm/octant2dm-run_0-expts_3_3-bins_-1_-1.root");

  f1.cd();
  double chisq_bf = 56.04;
  double dchisq_crit_68 = 1.00, dchisq_crit_90 = 2.71;
  chisqwithsyst_dm23sq_3->GetXaxis()->SetTitle("#Deltam^{2}_{32}");
  chisqwithsyst_dm23sq_3->GetYaxis()->SetTitle("#Delta#chi^{2}");
  TH1D * h68 = new TH1D(*chisqwithsyst_dm23sq_3);
  TH1D * h90 = new TH1D(*chisqwithsyst_dm23sq_3);
  TH1D * h2 = new TH1D(*chisqwithsyst_dm23sq_3);
  //h->GetXaxis()->SetTitle("#Delta#chi^{2}");
  //h2->GetXaxis()->SetTitle("#Delta#chi^{2}");
  for(int i = 1; i <= h68->GetNbinsX(); i++) {
    double chisq = h68->GetBinContent(i);
    h68->SetBinContent(i, chisq-chisq_bf-dchisq_crit_68);
    h90->SetBinContent(i, chisq-chisq_bf-dchisq_crit_90);
    h2->SetBinContent(i, chisq-chisq_bf);
  }
  f0.cd();
  TH2D * hdchisq_copy = new TH2D(*hdchisq);
  TH1D * hh = new TH1D("","",hdchisq_copy->GetNbinsY(), hdchisq_copy->GetYaxis()->GetXmin(), hdchisq_copy->GetYaxis()->GetXmax());
  for(int i = 1; i <= hh->GetNbinsX(); i++) {
    double dchisqmin = 999999;
    for(int j = 1; j <= hdchisq_copy->GetNbinsX(); j++) {
      double dchisq = hdchisq_copy->GetBinContent(j,i);
      if(dchisq < dchisqmin)
	dchisqmin = dchisq;
    }
    hh->SetBinContent(i, dchisqmin);
  }
  f2.cd();
  double chisq_bfo = 56.03;
  chisqwithsyst_dm23sq_3->GetYaxis()->SetTitle("#Delta#chi^{2}");
  TH1D * ho = new TH1D(*chisqwithsyst_dm23sq_3);
  for(int i = 1; i <= ho->GetNbinsX(); i++) {
    double chisq = ho->GetBinContent(i);
    ho->SetBinContent(i, chisq-chisq_bfo);
  }
  TCanvas * c = new TCanvas();
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.1);
  h68->Draw();
  TLine l(0.002,0,0.003,0);
  l->Draw();
  c->SaveAs("1D_contour_dm23sq.pdf");
  TGaxis::SetMaxDigits(3);
  h68->GetXaxis()->SetRangeUser(0.00225,0.00265);
  h68->GetXaxis()->SetNdivisions(510);
  h68->GetXaxis()->SetNoExponent(false);
  h68->Draw();
  TLine l2(0.00225,0,0.00265,0);
  l2->Draw();
  c->SaveAs("1D_contour_dm23sq_zoom.pdf");
  TGaxis::SetMaxDigits(5);
  h2->Draw("");
  c->SaveAs("1D_dchisq_dm23sq.pdf");
  c->SaveAs("1D_dchisq_dm23sq.eps");
  hh->Draw("same");
  c->SaveAs("1D_dchisq_dm23sq_profile.pdf");
  h68->GetXaxis()->SetRangeUser(0.00215,0.0028);
  h90->GetXaxis()->SetRangeUser(0.00215,0.0028);
  TLine l3(0.00215,0,0.0028,0);
  h68->SetLineColor(kRed);
  h90->Draw();
  h68->Draw("same");
  TLegend * leg = new TLegend(0.3,0.7,0.6,0.9);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  leg->AddEntry(h68,"const #Delta#chi^2 68% C.L.","l");
  leg->AddEntry(h90,"const #Delta#chi^2 90% C.L.","l");
  leg->Draw();
  l3->Draw();
  c->SaveAs("1D_contour_dm23sq_both.pdf");
  c->SaveAs("1D_contour_dm23sq_both.eps");
  //for TN
  h2->GetXaxis()->SetRangeUser(0.00215,0.0028);
  ho->GetXaxis()->SetRangeUser(0.00215,0.0028);
  ho->GetYaxis()->SetRangeUser(0, 4.8);
  ho->SetLineWidth(2);
  ho->SetLineColor(kRed);
  ho->Draw();
  h2->Draw("same");
  TLine l68(0.002145,1.00, 0.002805,1.00);
  TLine l90(0.002145,2.71, 0.002805,2.71);
  l68->SetLineStyle(3);
  l90->SetLineStyle(2);
  l68.Draw();
  l90.Draw();
  delete leg;
  leg = new TLegend(0.4,0.7,0.7,0.9);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  leg->AddEntry(h2,"Octant 1","l");
  leg->AddEntry(ho,"Octant 2","l");
  leg->Draw();
  c->SaveAs("1D_contour_dm23sq_lines.pdf");
  c->SaveAs("1D_contour_dm23sq_lines.eps"); //TN plot

  //find limits
  for(int i = 1; i <= h2->GetNbinsX(); i++) {
    if(h2->GetBinContent(i) > 2.71)
      cout << h2->GetBinCenter(i) << "\t" << h2->GetBinContent(i) << endl;
  }
  cout << "inversion point is 90% C.L. Boundary for octant 1" << endl;
  for(int i = 1; i <= h2->GetNbinsX(); i++) {
    if(h2->GetBinContent(i) > 1.00)
      cout << h2->GetBinCenter(i) << endl;
  }
  cout << "inversion point is 68% C.L. Boundary for octant 1" << endl;
  for(int i = 1; i <= ho->GetNbinsX(); i++) {
    if(ho->GetBinContent(i) > 2.71)
      cout << ho->GetBinCenter(i) << "\t" << ho->GetBinContent(i) << endl;
  }
  cout << "inversion point is 90% C.L. Boundary for octant 2" << endl;
  for(int i = 1; i <= ho->GetNbinsX(); i++) {
    if(ho->GetBinContent(i) > 1.00)
      cout << ho->GetBinCenter(i) << endl;
  }
  cout << "inversion point is 68% C.L. Boundary for octant 2" << endl;
}
