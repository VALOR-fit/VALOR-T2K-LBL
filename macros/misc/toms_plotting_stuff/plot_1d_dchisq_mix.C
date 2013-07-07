{
  TFile f1("$T2KOSCANA/2012results/1DFC/mixing23/1DFC-run_0-expts_3_3.root");
  TFile f2("$T2KOSCANA/2012results/1DFC/mixing23/fixdm/1D_mix_fixdm-run_0-expts_3_3.root");
  TFile f0("$T2KOSCANA/2012results/NEWRESULTS/datafit/dchisq_2all.root");
  TFile f3("$T2KOSCANA/2012results/1DFC/octant2/mix/octant2mix-run_0-expts_3_3-bins_-1_-1.root");

  f1.cd();
  double chisq_bf = 56.04;
  double dchisq_crit_68 = 1.00, dchisq_crit_90 = 2.71;
  chisqwithsyst_mixing23_3->GetYaxis()->SetTitle("#Delta#chi^{2}");
  TH1D * h68 = new TH1D(*chisqwithsyst_mixing23_3);
  TH1D * h90 = new TH1D(*chisqwithsyst_mixing23_3);
  TH1D * h2 = new TH1D(*chisqwithsyst_mixing23_3);
  //h->GetXaxis()->SetTitle("#Delta#chi^{2}");
  //h2->GetXaxis()->SetTitle("#Delta#chi^{2}");
  for(int i = 1; i <= h2->GetNbinsX(); i++) {
    double chisq = h2->GetBinContent(i);
    h68->SetBinContent(i, chisq-chisq_bf-dchisq_crit_68);
    h90->SetBinContent(i, chisq-chisq_bf-dchisq_crit_90);
    h2->SetBinContent(i, chisq-chisq_bf);
  }
  f2.cd();
  chisqwithsyst_mixing23_3->GetYaxis()->SetTitle("#Delta#chi^{2}");
  TH1D * hf = new TH1D(*chisqwithsyst_mixing23_3);
  for(int i = 1; i <= hf->GetNbinsX(); i++) {
    double chisq = hf->GetBinContent(i);
    hf->SetBinContent(i, chisq-chisq_bf);
  }
  f3.cd();
  double chisq_bfo = 56.03;
  chisqwithsyst_mixing23_3->GetYaxis()->SetTitle("#Delta#chi^{2}");
  TH1D * ho = new TH1D(*chisqwithsyst_mixing23_3);
  for(int i = 1; i <= ho->GetNbinsX(); i++) {
    double chisq = ho->GetBinContent(i);
    ho->SetBinContent(i, chisq-chisq_bfo);
  }
  f0.cd();
  TH2D * hdchisq_copy = new TH2D(*hdchisq);
  TH1D * hh = new TH1D("","",hdchisq_copy->GetNbinsX(), hdchisq_copy->GetXaxis()->GetXmin(), hdchisq_copy->GetXaxis()->GetXmax());
  TH1D * hhbf = new TH1D("","",hdchisq_copy->GetNbinsX(), hdchisq_copy->GetXaxis()->GetXmin(), hdchisq_copy->GetXaxis()->GetXmax());
  for(int i = 1; i <= hh->GetNbinsX(); i++) {
    double dchisqmin = 999999;
    for(int j = 1; j <= hdchisq_copy->GetNbinsY(); j++) {
      double dchisq = hdchisq_copy->GetBinContent(i,j);
      if(dchisq < dchisqmin)
	dchisqmin = dchisq;
    }
    hh->SetBinContent(i, dchisqmin);
    hhbf->SetBinContent(i, hdchisq_copy->GetBinContent(i, hdchisq_copy->GetXaxis()->FindBin(0.00244)));
  }
  TCanvas * c = new TCanvas();
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.1);
  h90->Draw();
  TLine l(0.9,0,1.0,0);
  l->Draw();
  c->SaveAs("1D_contour_mixing23.pdf");
  h90->GetXaxis()->SetRangeUser(0.94,0.97);
  h90->Draw();
  TLine l2(0.94,0,0.97,0);
  l2->Draw();
  c->SaveAs("1D_contour_mixing23_zoom.pdf");
  h2->Draw();
  c->SaveAs("1D_dchisq_mixing23.pdf");
  c->SaveAs("1D_dchisq_mixing23.eps");
  hh->Draw("same");
  c->SaveAs("1D_dchisq_mixing23_profile.pdf");
  h68->GetXaxis()->SetRangeUser(0.95,1.00);
  h90->GetXaxis()->SetRangeUser(0.95,1.00);
  TLine l3(0.95,0,1.00,0);
  h68->SetLineColor(kRed);
  h68->Draw();
  h90->Draw("same");
  TLegend * leg = new TLegend(0.5,0.7,0.8,0.9);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  leg->AddEntry(h68,"const #Delta#chi^2 68% C.L.","l");
  leg->AddEntry(h90,"const #Delta#chi^2 90% C.L.","l");
  leg->Draw();
  l3->Draw();
  c->SaveAs("1D_contour_mixing23_both.pdf");
  c->SaveAs("1D_contour_mixing23_both.eps");
  //for TN
  h2->GetXaxis()->UnZoom();
  //ho->GetXaxis()->SetRangeUser(0.95,1.00);
  ho->GetYaxis()->SetRangeUser(0,7);
  ho->SetLineColor(kRed);
  ho->SetLineWidth(2);
  ho->Draw();
  h2->Draw("same");
  TLine l68(0.8995,1.00,1.0005,1.00);
  TLine l90(0.8995,2.71,1.0005,2.71);
  l68->SetLineStyle(3);
  l90->SetLineStyle(2);
  l68.Draw();
  l90.Draw();
  delete leg;
  leg = new TLegend(0.5,0.7,0.8,0.9);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  leg->AddEntry(h2,"Octant 1","l");
  leg->AddEntry(ho,"Octant 2","l");
  leg->Draw();
  c->SaveAs("1D_contour_mixing23_lines.pdf");
  c->SaveAs("1D_contour_mixing23_lines.eps"); // the TN plot
  hf->SetLineColor(kBlue);
  h2->GetXaxis()->SetRangeUser(0.94,0.97);
  hf->Draw("same");
  c->SaveAs("1D_contour_mixing23_lines_withfix.pdf");
  hh->Draw();
  hhbf->Draw("same");
  c->SaveAs("temp.pdf");
  double xminhf = 0.94, xmaxhf = 0.97;
  hf->GetXaxis()->SetRangeUser(xminhf,xmaxhf);
  hf->Draw();
  TLine lhf90(xminhf, 2.71, xmaxhf, 2.71);

  TMarker m1(0.945,2.66,8);


  TMarker m2(0.950,2.56,8);


  TMarker m3(0.955,2.51,8);
  TMarker m4(0.960,2.49,8);
  TMarker m5(0.965,2.42,8);
  lhf90.Draw();
  m1.Draw();
  m2.Draw();
  m3.Draw();
  m4.Draw();
  m5.Draw();
  h2->Draw("same");
  c->SaveAs("compare_fc_constdchisq_fixdm.pdf");
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

