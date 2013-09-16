void plot_sensitivity_sinsq2theta(
		      const char * filename1 = "2013results/SENSITIVITY/sinsq2theta/run1-3bf/combination/final.root",
		      const char * filename2 = "2013results/SENSITIVITY/sinsq2theta/run1-3bfo2/combination/final.root",
		      double       in_mix    = 0.999,
		      double       in_dm     = 0.002441
)
{
  TFile fsenso1(filename1);
  TH2D * ho1_combchisq = 0;
  fsenso1.GetObject("chisq_surface",ho1_combchisq);

  //make the combined chisq surface a an average chisq surface (average)
  TH2D * ho1_avchisq = new TH2D(*ho1_combchisq);
  for(int ix = 1; ix <= ho1_combchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= ho1_combchisq->GetNbinsX(); iy++) {
      double chisq_combine = ho1_combchisq->GetBinContent(ix,iy);
      //hack for a missing job
      double chisq_average = chisq_combine / (double)3;
      if(chisq_combine < 1)
	chisq_average = 80;
      else if(chisq_combine < 100)
	chisq_average = chisq_combine / (double)1;
      else if(chisq_combine < 200)
	chisq_average = chisq_combine / (double)2;
      ho1_avchisq->SetBinContent(ix,iy, chisq_average);
    }//iy
  }//ix

  //find the minimum
  double chisq_min = ho1_avchisq->GetMinimum();
  int minbin = ho1_avchisq->GetMinimumBin();
  int minx, miny, minz;
  ho1_avchisq->GetBinXYZ(minbin, minx, miny, minz);
  cout << endl << endl << "min chisq " << chisq_min << " at " << ho1_avchisq->GetXaxis()->GetBinCenter(minx) << " " << ho1_avchisq->GetYaxis()->GetBinCenter(miny) << endl;

  //make the average chisq surface a dchisq surface (subtract min)
  TH2D * ho1_sens = new TH2D(*ho1_combchisq);
  for(int ix = 1; ix <= ho1_combchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= ho1_combchisq->GetNbinsX(); iy++) {
      double chisq_av = ho1_avchisq->GetBinContent(ix,iy);
      double dchisq   = chisq_av - chisq_min;
      ho1_sens->SetBinContent(ix, iy, dchisq);
    }//iy
  }//ix

  //set up the sensitivity histogram for contour drawing
  ho1_sens->SetContour(1);
  ho1_sens->SetContourLevel(0,4.61);
  ho1_sens->SetLineStyle(1);

  //////
  // do it all again for the second octant

  TFile fsenso2(filename2);
  TH2D * ho2_combchisq = 0;
  fsenso2.GetObject("chisq_surface",ho2_combchisq);

  //make the combined chisq surface a an average chisq surface (average)
  TH2D * ho2_avchisq = new TH2D(*ho2_combchisq);
  for(int ix = 1; ix <= ho2_combchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= ho2_combchisq->GetNbinsX(); iy++) {
      double chisq_combine = ho2_combchisq->GetBinContent(ix,iy);
      //hack for a missing job
      double chisq_average = chisq_combine / (double)3;
      if(chisq_combine < 1)
	chisq_average = 76;
      else if(chisq_combine < 100)
	chisq_average = chisq_combine / (double)1;
      else if(chisq_combine < 200)
	chisq_average = chisq_combine / (double)2;
      ho2_avchisq->SetBinContent(ix,iy, chisq_average);
    }//iy
  }//ix

  //find the minimum
  double chisq_min = ho2_avchisq->GetMinimum();
  int minbin = ho2_avchisq->GetMinimumBin();
  int minx, miny, minz;
  ho2_avchisq->GetBinXYZ(minbin, minx, miny, minz);
  cout << endl << endl << "min chisq " << chisq_min << " at " << ho2_avchisq->GetXaxis()->GetBinCenter(minx) << " " << ho2_avchisq->GetYaxis()->GetBinCenter(miny) << endl;

  //make the average chisq surface a dchisq surface (subtract min)
  TH2D * ho2_sens = new TH2D(*ho2_combchisq);
  for(int ix = 1; ix <= ho2_combchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= ho2_combchisq->GetNbinsX(); iy++) {
      double chisq_av = ho2_avchisq->GetBinContent(ix,iy);
      double dchisq   = chisq_av - chisq_min;
      ho2_sens->SetBinContent(ix, iy, dchisq);
    }//iy
  }//ix

  //set up the sensitivity histogram for contour drawing
  ho2_sens->SetContour(1);
  ho2_sens->SetContourLevel(0,4.61);
  ho2_sens->SetLineStyle(2);

  //Input marker
  TMarker * m_input = new TMarker(in_mix,in_dm,29);
  m_input->SetMarkerSize(1.8);

  //
  //Draw
  //
  // param space
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.16);
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.0;
  double mixing23_max   = 1.0;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 6E-3;
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
  osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  //get run1-3 result
  TFile f_t2k("inputs/contours/t2k_2013paper_likelihoodratio_contour_and_deltachisq.root");
  TGraph * g_t2k90o1 = 0;
  f_t2k.GetObject("g_contour_90_o1",g_t2k90o1);
  g_t2k90o1->SetLineStyle(1);
  //g_t2k90o1->SetLineColor(kRed);
  TGraph * g_t2k90o2 = 0;
  f_t2k.GetObject("g_contour_90_o2",g_t2k90o2);
  g_t2k90o2->SetLineStyle(2);
  g_t2k90o2->SetLineColor(kRed);

  //get externals
  //minos
  TFile fminos("$T2KOSCANA/inputs/contours/MINOS_NumuDisappearance_LogLikelihoodSurface_May2013.root");
  TH2D * h_minos = 0;
  fminos.GetObject("LogLikelihoodSurface",h_minos);
  h_minos->SetContour(1);
  h_minos->SetContourLevel(0,4.61/2);
  h_minos->SetLineStyle(1);
  h_minos->SetLineColor(kRed);
  //sk L/E
  TFile fskle("$T2KOSCANA/inputs/contours/contours_skloe.root");
  TGraph * g_skle = 0;
  fskle.GetObject("sk3f_loe_90cl",g_skle);
  g_skle->SetLineColor(kGreen+2);
  g_skle->SetLineStyle(1);
  g_skle->SetName("g_skle");
  //sk zenith
  TFile fskz("$T2KOSCANA/inputs/contours/contours_skzen.root");
  TGraph * g_skzen = 0;
  fskz.GetObject("sk3f_23vars_overlay_normal_q13_fixed",g_skzen);
  g_skzen->SetLineColor(kBlue);
  g_skzen->SetLineStyle(1);
  g_skzen->SetName("g_skzen");


  ho1_sens->SetLineStyle(2);
  ho1_sens->SetLineWidth(2);
  g_t2k90o1->SetLineWidth(2);
  h_minos->SetLineWidth(2);
  g_skle->SetLineWidth(2);
  g_skzen->SetLineWidth(2);

  //draw
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.0038);
  osc_param_space->Draw();
  h_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  g_t2k90o1->Draw("l");
  //g_t2k90o2->Draw("l");
  ho1_sens->Draw("cont3same");
  //ho2_sens->Draw("cont3same");
  m_input->Draw();
  //legend
  TLegend * leg = new TLegend(0.18, 0.65, 0.45, 0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);
  //leg->AddEntry(g_t2k90o1, "Run 1-3 (#theta_{23} #leq #pi/4) 90% C.L. result","l");
  //leg->AddEntry(g_t2k90o2, "Run 1-3 (#theta_{23} #geq #pi/4) 90% C.L. result","l");
  //leg->AddEntry(ho1_sens,  "Run 1-4 (#theta_{23} #leq #pi/4) 90% C.L. sensitivity","l");
  //leg->AddEntry(ho2_sens,  "Run 1-4 (#theta_{23} #geq #pi/4) 90% C.L. sensitivity","l");
  leg->AddEntry(g_t2k90o1, "T2K Run 1-3 3#nu 90% C.L. result","l");
  leg->AddEntry(ho1_sens,  "T2k Run 1-4 3#nu 90% C.L. sensitivity","l");
  leg->AddEntry("g_skzen", "SuperK zenith 2012 3#nu 90% C.L.","l");
  leg->AddEntry("g_skle",  "SuperK L/E 2012 2#nu 90% C.L.","l");
  leg->AddEntry(h_minos,   "MINOS 2013 2#nu 90% C.L.","l");
  leg->AddEntry(m_input,   "Run 1-3 best fit","p");
  leg->Draw();
  c->SaveAs("run123_with_run1234a_sensitivity.pdf");

}
