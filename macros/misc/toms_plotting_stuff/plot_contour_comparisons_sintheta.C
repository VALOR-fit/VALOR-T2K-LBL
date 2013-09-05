//plot valor sintheta contour vs world sintheta data
void plot_contour_comparisons_sintheta()
{
  //valor contour
  TFile fval("2012results/DATA/sintheta23/combine/combination/fitresults_3.root");
  TH2D * hval = 0;
  fval.GetObject("hdchisq",hval);
  TH2D * hval_68 = new TH2D(*hval);
  hval_68->SetContour(1);
  hval_68->SetContourLevel(0,2.30);
  hval_68->SetLineStyle(2);
  TH2D * hval_90 = new TH2D(*hval);
  hval_90->SetContour(1);
  hval_90->SetContourLevel(0,4.61);
  hval_90->SetLineStyle(1);

  //external contours
  TFile fext("inputs/contours/contours_nont2k.root");
  TGraph * gminos = 0;
  fext.GetObject("minos2013_combined_90CL_normal_sinsqtheta",gminos);
  gminos->SetName("gminos");
  TGraph * gsk = 0;
  fext.GetObject("sk2012_zenith3nu_90CL_normal_sinsqtheta",gsk);
  gsk->SetName("gsk");
  
  //valor best fit
  TMarker bfval(0.514,0.002441,29);
  bfval.SetMarkerSize(1.8);
  bfval.SetMarkerColor(kBlack);
 
  //canvas
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.16);

  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.3;
  double mixing23_max   = 0.75.;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1.8E-3;
  double dm23sq_max     = 4.2E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + TMath::CeilNint ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + TMath::CeilNint ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
  osc_param_space->GetXaxis()->SetTitle("sin^{2}(#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  //osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  //legend
  TLegend * leg = new TLegend(0.17,0.7,0.43,0.93);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);

  //style
  gminos->SetLineColor(kRed);
  gsk->SetLineColor(kBlue);
  hval_90->SetLineWidth(2);
  hval_68->SetLineWidth(2);
  gminos->SetLineWidth(1.5);
  gsk->SetLineWidth(1.5);

  //draw
  osc_param_space->Draw();
  gminos->Draw("l");
  gsk->Draw("l");
  hval_90->Draw("cont3same");
  bfval.Draw();
  //legend
  leg->AddEntry(hval_90,"T2K 90% C.L.","l");
  leg->AddEntry("gminos","MINOS 2013 90% C.L.","l");
  leg->AddEntry("gsk","SuperK 2012 zenith 90% C.L.","l");
  leg->Draw();
  //save
  c->SaveAs("compare_contours_external_sintheta.pdf");
  c->SaveAs("compare_contours_external_sintheta.eps");

}
