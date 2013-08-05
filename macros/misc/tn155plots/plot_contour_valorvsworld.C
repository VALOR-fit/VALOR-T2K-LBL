void plot_contour_valorvsworld
(

)
{
  using std::cout;
  using namespace t2kana;

  double linewidth_val = 3; // for the new results
  double linewidth_ext = 1.5; // for old/external results
  double marker_size   = 1.8;

  //VaLOR data results
  TFile f_dc_normal("2013results/DATA/run1-4_nuis3fixdcp_sintheta_normal/combination/run1-4_nuis3fixdcp_sintheta_normal_sinsq2theta.root");
  TFile f_bf_normal("2013results/DATA/run1-4_nuis3fixdcp_sintheta_normal/combination/run1-4_nuis3fixdcp_sintheta_normal-run_0-expt_6.root");
  //TFile f_fc_normal("");

  //Fill in TMarker's with the best fit points
  FitParamList * fparm_normal = 0;
  f_bf_normal.GetObject("fparm",fparm_normal);
  double dm32sq_normal      = fparm_normal->BestFitValue(conv::kPhypDm23sq);
  double mix23_normal_s23sq = fparm_normal->BestFitValue(conv::kPhypMixing23);
  double mix23_normal_sinsq2theta23 = 1 - ((1 - (mix23_normal_s23sq * 2)) * (1 - (mix23_normal_s23sq * 2)));
  TMarker bf_normal(mix23_normal_sinsq2theta23, dm32sq_normal, 29);
  bf_normal.SetMarkerSize(marker_size);
  bf_normal.SetMarkerColor(kBlack);

  //normal const dchisq
  TH2D * h_dc_normal_90 = 0, * h_dc_normal_68 = 0;
  f_dc_normal.GetObject("hdchisq_90",h_dc_normal_90);
  f_dc_normal.GetObject("hdchisq_68",h_dc_normal_68);
  h_dc_normal_68->SetLineWidth(linewidth_val);
  h_dc_normal_90->SetLineWidth(linewidth_val);

  //externals + olds
  //minos
  TFile fminos("$T2KOSCANA/inputs/contours/MINOS_NumuDisappearance_LogLikelihoodSurface_May2013.root");
  TH2D * h_minos = 0;
  fminos.GetObject("LogLikelihoodSurface",h_minos);
  h_minos->SetContour(1);
  h_minos->SetContourLevel(0,4.61/2);
  h_minos->SetLineStyle(1);
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
  //icecube
  TFile fextern("$T2KOSCANA/inputs/contours/contours_nont2k.root");
  TGraph * g_ice = 0;
  fextern.GetObject("icecube20102011_90CL",g_ice);
  g_ice->SetLineColor(kBlue);
  g_ice->SetLineStyle(1);
  g_ice->SetName("g_ice");
  //t2k run1+2
  TFile fold("$T2KOSCANA/inputs/contours/contours_oldt2k.root");
  TGraph * g_run12_68_statsyst = 0;
  fold.GetObject("t2k2010_run1p2_68CL_statsyst",g_run12_68_statsyst);
  TGraph * g_run12_90_statsyst = 0;
  fold.GetObject("t2k2010_run1p2_90CL_statsyst",g_run12_90_statsyst);
  TGraph * g_run12_68_stat = 0;
  fold.GetObject("t2k2010_run1p2_68CL_stat",g_run12_68_stat);
  TGraph * g_run12_90_stat = 0;
  fold.GetObject("t2k2010_run1p2_90CL_stat",g_run12_90_stat);
  g_run12_68_stat->SetName("g_run12_68_stat");
  g_run12_90_stat->SetName("g_run12_90_stat");
  g_run12_68_statsyst->SetName("g_run12_68_statsyst");
  g_run12_90_statsyst->SetName("g_run12_90_statsyst");
  //t2k run1+2+3
  TFile frun123("$T2KOSCANA/inputs/contours/t2k_2013paper_likelihoodratio_contour_and_deltachisq.root");
  TGraph * g_run123_90_o1 = 0;
  frun123.GetObject("g_contour_90_o1",g_run123_90_o1);
  g_run123_90_o1->SetName("g_run123_90_o1");
  TGraph * g_run123_68_o1 = 0;
  frun123.GetObject("g_contour_68_o1",g_run123_68_o1);
  g_run123_68_o1->SetName("g_run123_68_o1");
  TGraph * g_run123_90_o2 = 0;
  frun123.GetObject("g_contour_90_o2",g_run123_90_o2);
  g_run123_90_o2->SetName("g_run123_90_o2");
  TGraph * g_run123_68_o2 = 0;
  frun123.GetObject("g_contour_68_o2",g_run123_68_o2);
  g_run123_68_o2->SetName("g_run123_68_o2");

  //width
  g_skle->SetLineWidth(linewidth_ext);
  g_skzen->SetLineWidth(linewidth_ext);
  g_ice->SetLineWidth(linewidth_ext);
  g_run12_68_stat->SetLineWidth(linewidth_ext);
  g_run12_90_stat->SetLineWidth(linewidth_ext);
  g_run12_68_statsyst->SetLineWidth(linewidth_ext);
  g_run12_90_statsyst->SetLineWidth(linewidth_ext);
  h_minos->SetLineWidth(linewidth_ext);
  g_run123_68_o1->SetLineWidth(linewidth_ext);
  g_run123_90_o1->SetLineWidth(linewidth_ext);
  g_run123_68_o2->SetLineWidth(linewidth_ext);
  g_run123_90_o2->SetLineWidth(linewidth_ext);

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
  //osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  TLegend * leg = new TLegend(0.17,0.65,0.43,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);

  //color
  g_run123_90_o1->SetLineStyle(3);
  g_run123_90_o2->SetLineStyle(2);
  g_run123_90_o1->SetLineColor(kBlack);
  g_run123_90_o2->SetLineColor(kBlack);
  h_minos->SetLineColor(kGreen+2);
  g_skle->SetLineColor(kBlue);
  g_skzen->SetLineStyle(2);
  g_skzen->SetLineColor(kBlue);


  //Draw Run 1-4 constant dchisq only with externals + old t2k
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0042);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  g_skzen->Draw("l");
  g_skle->Draw("l");
  h_minos->Draw("cont3same");
  g_run123_90_o1->Draw("l");
  g_run123_90_o2->Draw("l");
  h_dc_normal_90->Draw("cont3same");
  bf_normal.Draw();
  TLegend * leg2 = new TLegend(0.185,0.58,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(h_dc_normal_90, "T2K Run1-4a 3#nu (#Deltam^{2}_{32}>0) 90% CL","l");
  leg2->AddEntry("g_run123_90_o1", "T2K Run1-3 3#nu (#theta_{23} #leq #pi/4) 90% CL", "l");
  leg2->AddEntry("g_run123_90_o2", "T2K Run1-3 3#nu (#theta_{23} #geq #pi/4) 90% CL", "l");
  leg2->AddEntry(h_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg2->AddEntry("g_skzen", "SK zenith 2012 3#nu 90% CL", "l");
  leg2->AddEntry("g_skle", "SK L/E 2012 2#nu 90% CL", "l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(0.06);
  leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.17, 0.45, 0.30);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bf_normal, "T2K 3#nu (#Deltam^{2}_{32}>0) best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(0.06);
  legbf->Draw();
  c->SaveAs("contour_compare_external_sinsq2theta23.pdf");
  c->SaveAs("contour_compare_external_sinsq2theta23.png");
  c->SaveAs("contour_compare_external_sinsq2theta23.eps");
}

