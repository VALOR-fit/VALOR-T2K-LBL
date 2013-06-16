//plot up to 4 contours

void contour_plots(
	   const char * f1,    // filename
	   const char * f2,
	   const char * f3,
	   const char * f4,
	   const char * h1,    // histogram name
	   const char * h2,
	   const char * h3,
	   const char * h4,
	   const char * leg1,  // legend
	   const char * leg2,
	   const char * leg3,
	   const char * leg4,
	   const char * tag,   // output tag
	   bool         plot_bestfit = false, //plot best fit?
	   double       m1 = -1,  // mixing angle point of best fit. If negative, find min from surface
	   double       d1 = -1,  // dm32sq point of best fit. If negative, find min from surface
	   double       m2 = -1,
	   double       d2 = -1,
	   double       m3 = -1,
	   double       d3 = -1,
	   double       m4 = -1,
	   double       d4 = -1
)
{
  // get histograms
  TFile f_1(f1);
  TFile f_2(f2);
  TFile f_3(f3);
  TFile f_4(f4);
  TH2D * h_1 = 0, * h_2 = 0, * h_3 = 0, * h_4 = 0;
  f_1.GetObject(h1,h_1);
  f_2.GetObject(h2,h_2);
  f_3.GetObject(h3,h_3);
  f_4.GetObject(h4,h_4);

  //get the 68% & 90% histograms
  TH2D * h1_90 = new TH2D(*h_1);
  TH2D * h1_68 = new TH2D(*h_1);
  TH2D * h2_90 = new TH2D(*h_2);
  TH2D * h2_68 = new TH2D(*h_2);
  TH2D * h3_90 = new TH2D(*h_3);
  TH2D * h3_68 = new TH2D(*h_3);
  TH2D * h4_90 = new TH2D(*h_4);
  TH2D * h4_68 = new TH2D(*h_4);
  //and set contour correctly
  h1_90->SetContour(1);
  h1_90->SetContourLevel(0, 4.61);
  h1_68->SetContour(1);
  h1_68->SetContourLevel(0, 2.30);
  h2_90->SetContour(1);
  h2_90->SetContourLevel(0, 4.61);
  h2_68->SetContour(1);
  h2_68->SetContourLevel(0, 2.30);
  h3_90->SetContour(1);
  h3_90->SetContourLevel(0, 4.61);
  h3_68->SetContour(1);
  h3_68->SetContourLevel(0, 2.30);
  h4_90->SetContour(1);
  h4_90->SetContourLevel(0, 4.61);
  h4_68->SetContour(1);
  h4_68->SetContourLevel(0, 2.30);

  //style
  //colour
  if(h1_90) h1_90->SetLineColor(kBlack);
  if(h1_68) h1_68->SetLineColor(kBlack);
  if(h2_90) h2_90->SetLineColor(kRed);
  if(h2_68) h2_68->SetLineColor(kRed);
  if(h3_90) h3_90->SetLineColor(kBlue);
  if(h3_68) h3_68->SetLineColor(kBlue);
  if(h4_90) h4_90->SetLineColor(kGreen+2);
  if(h4_68) h4_68->SetLineColor(kGreen+2);
  //width
  if(h1_90) h1_90->SetLineWidth(2);
  if(h1_68) h1_68->SetLineWidth(2);
  if(h2_90) h2_90->SetLineWidth(2);
  if(h2_68) h2_68->SetLineWidth(2);
  if(h3_90) h3_90->SetLineWidth(2);
  if(h3_68) h3_68->SetLineWidth(2);
  if(h4_90) h4_90->SetLineWidth(2);
  if(h4_68) h4_68->SetLineWidth(2);
  //style
  if(h1_68) h1_68->SetLineStyle(2);
  if(h2_68) h2_68->SetLineStyle(2);
  if(h3_68) h3_68->SetLineStyle(2);
  if(h4_68) h4_68->SetLineStyle(2);
  if(h1_90) h1_90->SetLineStyle(1);
  if(h2_90) h2_90->SetLineStyle(1);
  if(h3_90) h3_90->SetLineStyle(1);
  if(h4_90) h4_90->SetLineStyle(1);

  //make the markers
  //1
  if(m1 < 0 || d1 < 0) {
    double z1;
    if(h1_68) h1_68->GetMinimumBin(m1,d1,z1);
  }
  TMarker b1(m1,d1,29);
  b1.SetMarkerSize(1.8);
  b1.SetMarkerColor(kBlack);
  //2
  if(m2 < 0 || d2 < 0) {
    double z2;
    if(h2_68) h2_68->GetMinimumBin(m2,d2,z2);
  }
  TMarker b2(m2,d2,25);
  b2.SetMarkerSize(1.8);
  b2.SetMarkerColor(kRed);
  //3
  if(m3 < 0 || d3 < 0) {
    double z3;
    if(h3_68) h3_68->GetMinimumBin(m3,d3,z3);
  }
  TMarker b3(m3,d3,24);
  b3.SetMarkerSize(1.8);
  b3.SetMarkerColor(kBlue);
  //4
  if(m4 < 0 || d4 < 0) {
    double z4;
    if(h4_68) h4_68->GetMinimumBin(m4,d4,z4);
  }
  TMarker b4(m4,d4,27);
  b4.SetMarkerSize(1m.8);
  b4.SetMarkerColor(kGreen+2);

  //get externals
  TFile fext("$T2KOSCANA/inputs/contours/contours_nont2k.root");
  TGraph * g_minos = 0;
  fext.GetObject("minos2012_combined_90CL",g_minos);
  g_minos->SetLineColor(kBlack);
  TGraph * g_skle = 0;
  fext.GetObject("sk2012_le2nu_90CL",g_skle);
  g_skle->SetLineColor(kGreen+2);
  TGraph * g_skzen = 0;
  fext.GetObject("sk2012_zenith3nu_90CL",g_skzen);
  g_skzen->SetLineColor(kBlue);
  //names
  g_minos->SetName("g_minos");
  g_skle->SetName("g_skle");
  g_skzen->SetName("g_skzen");

  // param space
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.16);
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.9;
  double mixing23_max   = 1.;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 2E-3;
  double dm23sq_max     = 3E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + (int) ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
  osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  //osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  TLegend * leg = new TLegend(0.17,0.68,0.53,0.93);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);

  //Draw
  osc_param_space->Draw();
  if(h1_90) h1_90->Draw("cont3same");
  if(h1_68) h1_68->Draw("cont3same");
  if(h2_90) h2_90->Draw("cont3same");
  if(h2_68) h2_68->Draw("cont3same");
  if(h3_90) h3_90->Draw("cont3same");
  if(h3_68) h3_68->Draw("cont3same");
  if(h4_90) h4_90->Draw("cont3same");
  if(h4_68) h4_68->Draw("cont3same");
  if(h1_90) leg->AddEntry(h1_90, Form("%s 90%% C.L.", leg1), "l");
  if(h1_68) leg->AddEntry(h1_68, Form("%s 68%% C.L.", leg1), "l");
  if(h2_90) leg->AddEntry(h2_90, Form("%s 90%% C.L.", leg2), "l");
  if(h2_68) leg->AddEntry(h2_68, Form("%s 68%% C.L.", leg2), "l");
  if(h3_90) leg->AddEntry(h3_90, Form("%s 90%% C.L.", leg3), "l");
  if(h3_68) leg->AddEntry(h3_68, Form("%s 68%% C.L.", leg3), "l");
  if(h4_90) leg->AddEntry(h4_90, Form("%s 90%% C.L.", leg4), "l");
  if(h4_68) leg->AddEntry(h4_68, Form("%s 68%% C.L.", leg4), "l");
  leg->Draw();
  if(plot_bestfit) {
    b1.Draw();
    b2.Draw();
    b3.Draw();
    b4.Draw();
  }
  c->SaveAs(Form("compare_%s_only.pdf",tag));

  //Draw
  leg->Clear();
  osc_param_space->Draw();
  if(h1_90) h1_90->Draw("cont3same");
  //if(h1_68) h1_68->Draw("cont3same");
  if(h2_90) h2_90->Draw("cont3same");
  //if(h2_68) h2_68->Draw("cont3same");
  if(h3_90) h3_90->Draw("cont3same");
  //if(h3_68) h3_68->Draw("cont3same");
  if(h4_90) h4_90->Draw("cont3same");
  //if(h4_68) h4_68->Draw("cont3same");
  if(h1_90) leg->AddEntry(h1_90, Form("%s 90%%", leg1), "l");
  //if(h1_68) leg->AddEntry(h1_68, Form("%s 68%%", leg1), "l");
  if(h2_90) leg->AddEntry(h2_90, Form("%s 90%%", leg2), "l");
  //if(h2_68) leg->AddEntry(h2_68, Form("%s 68%%", leg2), "l");
  if(h3_90) leg->AddEntry(h3_90, Form("%s 90%%", leg3), "l");
  //if(h3_68) leg->AddEntry(h3_68, Form("%s 68%%", leg3), "l");
  if(h4_90) leg->AddEntry(h4_90, Form("%s 90%%", leg4), "l");
  //if(h4_68) leg->AddEntry(h4_68, Form("%s 68%%", leg4), "l");
  g_minos->Draw("l");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry("g_minos", "MINOS 2012 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->Draw();
  c->SaveAs(Form("compare_%s_ext.pdf",tag));
}
