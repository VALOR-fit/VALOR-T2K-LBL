void plot_contour_comparisons
(

)
{
  //VaLOR data results
  TFile f_dc_123("2012results/DATA/bigger_surface/28/fitresults_3.root");
  TFile f_dc_12("2012results/DATA/bigger_surface/28/fitresults_4.root");
  TFile f_dc_3("2012results/DATA/bigger_surface/28/fitresults_5.root");
  TFile f_dcall_123("2012results/DATA/2all/fitresults_3.root");
  TFile f_dcall_12("2012results/DATA/2all/fitresults_4.root");
  TFile f_dcall_3("2012results/DATA/2all/fitresults_5.root");
  TFile f_fc_123("2012results/FC/fc_contours.root");

  //28 const dchisq
  TH2D * h_dc_123 = 0;
  f_dc_123.GetObject("hdchisq",h_dc_123);
  TH2D * h_dc_123_68 = new TH2D(*h_dc_123);
  TH2D * h_dc_123_90 = new TH2D(*h_dc_123);
  h_dc_123_68->SetContour(1);
  h_dc_123_68->SetContourLevel(0,2.30);
  h_dc_123_68->SetLineStyle(2);
  h_dc_123_90->SetContour(1);
  h_dc_123_90->SetContourLevel(0,4.61);
  h_dc_123_90->SetLineStyle(1);

  TH2D * h_dc_12 = 0;
  f_dc_12.GetObject("hdchisq",h_dc_12);
  TH2D * h_dc_12_68 = new TH2D(*h_dc_12);
  TH2D * h_dc_12_90 = new TH2D(*h_dc_12);
  h_dc_12_68->SetContour(1);
  h_dc_12_68->SetContourLevel(0,2.30);
  h_dc_12_68->SetLineStyle(2);
  h_dc_12_90->SetContour(1);
  h_dc_12_90->SetContourLevel(0,4.61);
  h_dc_12_90->SetLineStyle(1);

  TH2D * h_dc_3 = 0;
  f_dc_3.GetObject("hdchisq",h_dc_3);
  TH2D * h_dc_3_68 = new TH2D(*h_dc_3);
  TH2D * h_dc_3_90 = new TH2D(*h_dc_3);
  h_dc_3_68->SetContour(1);
  h_dc_3_68->SetContourLevel(0,2.30);
  h_dc_3_68->SetLineStyle(2);
  h_dc_3_90->SetContour(1);
  h_dc_3_90->SetContourLevel(0,4.61);
  h_dc_3_90->SetLineStyle(1);

  //28 FC
  TH2D * h_fc_123_68 = 0;
  TH2D * h_fc_123_90 = 0;
  f_fc_123.GetObject("hCL68StatSyst",h_fc_123_68);
  f_fc_123.GetObject("hCL90StatSyst",h_fc_123_90);
  h_fc_123_68->SetContour(1);
  h_fc_123_68->SetContourLevel(0,0.0);
  h_fc_123_68->SetLineStyle(2);
  h_fc_123_90->SetContour(1);
  h_fc_123_90->SetContourLevel(0,0.0);
  h_fc_123_90->SetLineStyle(1);

  //2all const dchisq
  TH2D * h_dcall_123 = 0;
  f_dcall_123.GetObject("hdchisq",h_dcall_123);
  TH2D * h_dcall_123_68 = new TH2D(*h_dcall_123);
  TH2D * h_dcall_123_90 = new TH2D(*h_dcall_123);
  h_dcall_123_68->SetContour(1);
  h_dcall_123_68->SetContourLevel(0,2.30);
  h_dcall_123_68->SetLineStyle(2);
  h_dcall_123_90->SetContour(1);
  h_dcall_123_90->SetContourLevel(0,4.61);
  h_dcall_123_90->SetLineStyle(1);

  TH2D * h_dcall_12 = 0;
  f_dcall_12.GetObject("hdchisq",h_dcall_12);
  TH2D * h_dcall_12_68 = new TH2D(*h_dcall_12);
  TH2D * h_dcall_12_90 = new TH2D(*h_dcall_12);
  h_dcall_12_68->SetContour(1);
  h_dcall_12_68->SetContourLevel(0,2.30);
  h_dcall_12_68->SetLineStyle(2);
  h_dcall_12_90->SetContour(1);
  h_dcall_12_90->SetContourLevel(0,4.61);
  h_dcall_12_90->SetLineStyle(1);

  TH2D * h_dcall_3 = 0;
  f_dcall_3.GetObject("hdchisq",h_dcall_3);
  TH2D * h_dcall_3_68 = new TH2D(*h_dcall_3);
  TH2D * h_dcall_3_90 = new TH2D(*h_dcall_3);
  h_dcall_3_68->SetContour(1);
  h_dcall_3_68->SetContourLevel(0,2.30);
  h_dcall_3_68->SetLineStyle(2);
  h_dcall_3_90->SetContour(1);
  h_dcall_3_90->SetContourLevel(0,4.61);
  h_dcall_3_90->SetLineStyle(1);


  //externals + olds
  TFile fext("$T2KOSCANA/inputs/contours/contours_nont2k.root");
  TGraph * g_minos = 0;
  fext.GetObject("minos2012_combined_90CL",g_minos);
  TGraph * g_skle = 0;
  fext.GetObject("sk2012_le2nu_90CL",g_skle);
  TGraph * g_skzen = 0;
  fext.GetObject("sk2012_zenith3nu_90CL",g_skzen);

  TFile fold("$T2KOSCANA/inputs/contours/contours_oldt2k.root");
  TGraph * g_run12_68_statsyst = 0;
  fold.GetObject("t2k2010_run1p2_68CL_statsyst",g_run12_68_statsyst);
  TGraph * g_run12_90_statsyst = 0;
  fold.GetObject("t2k2010_run1p2_90CL_statsyst",g_run12_90_statsyst);
  TGraph * g_run12_68_stat = 0;
  fold.GetObject("t2k2010_run1p2_68CL_stat",g_run12_68_stat);
  TGraph * g_run12_90_stat = 0;
  fold.GetObject("t2k2010_run1p2_90CL_stat",g_run12_90_stat);

  //so we can put the graphs in the legends
  g_minos->SetName("g_minos");
  g_skle->SetName("g_skle");
  g_skzen->SetName("g_skzen");
  g_run12_68_stat->SetName("g_run12_68_stat");
  g_run12_90_stat->SetName("g_run12_90_stat");
  g_run12_68_statsyst->SetName("g_run12_68_statsyst");
  g_run12_90_statsyst->SetName("g_run12_90_statsyst");


  // param space
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.05);
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.6;
  double mixing23_max   = 1.;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 4E-3;
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
  osc_param_space->GetYaxis()->SetTitleOffset(1.65);

  TLegend * leg = new TLegend(0.17,0.72,0.40,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  //colors
  h_dc_123_68->SetLineColor(kRed);
  h_dc_123_90->SetLineColor(kRed);
  h_dc_12_68->SetLineColor(kBlue);
  h_dc_12_90->SetLineColor(kBlue);
  h_dc_3_68->SetLineColor(kGreen+2);
  h_dc_3_90->SetLineColor(kGreen+2);
  h_dcall_123_68->SetLineColor(kBlack);
  h_dcall_123_90->SetLineColor(kBlack);
  h_dcall_12_68->SetLineColor(kBlue+5);
  h_dcall_12_90->SetLineColor(kBlue+5);
  h_dcall_3_68->SetLineColor(kGreen+2+5);
  h_dcall_3_90->SetLineColor(kGreen+2+5);
  h_fc_123_68->SetLineColor(kBlack);
  h_fc_123_90->SetLineColor(kBlack);
  g_run12_68_statsyst->SetLineColor(kRed);
  g_run12_90_statsyst->SetLineColor(kRed);

  //get rid of the crap at the top of the FC contour
  h_fc_123_68->GetYaxis()->SetRangeUser(0,0.33);
  h_fc_123_90->GetYaxis()->SetRangeUser(0,0.33);


  //Run 1-3 const dchisq vs FC
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_dc_123_68->Draw("cont3same");
  h_dc_123_90->Draw("cont3same");
  h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_123_68, "const#Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_123_90, "const#Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_fc_123_68, "FC 68% CL","l");
  leg->AddEntry(h_fc_123_90, "FC 90% CL","l");
  leg->Draw();
  c->SaveAs("compare_FC_constchisq.eps");

  //Run 1-3 const dchisq vs FC + externals to this plot
  osc_param_space->GetXaxis()->SetRangeUser(0.8,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0035);
  osc_param_space->Draw();  
  h_dc_123_68->Draw("cont3same");
  h_dc_123_90->Draw("cont3same");
  h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_123_68, "const#Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_123_90, "const#Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_fc_123_68, "FC 68% CL","l");
  leg->AddEntry(h_fc_123_90, "FC 90% CL","l");
  g_minos->Draw("l");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry("g_minos", "MINOS 2012 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->Draw();
  c->SaveAs("compare_FC_constchisq_external.eps");

  //run1+2 2012 const dchisq vs 2010 FC
  osc_param_space->GetXaxis()->SetRangeUser(0.75,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0035);
  osc_param_space->Draw();
  h_dc_12_68->Draw("cont3same");
  h_dc_12_90->Draw("cont3same");
  g_run12_68_statsyst->Draw("l");
  g_run12_90_statsyst->Draw("l");
  leg->Clear();
  leg->AddEntry(h_dc_12_68, "2012 const#Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_12_90, "2012 const#Delta#chi^{2} 90% CL","l");
  leg->AddEntry("g_run12_68_statsyst", "2010 FC 68% CL", "l");
  leg->AddEntry("g_run12_90_statsyst", "2010 FC 90% CL", "l");
  leg->Draw();
  c->SaveAs("compare_run12_2012_2010.eps");

  //runs 123,12,3 constdchisq comparisons
  osc_param_space->GetXaxis()->SetRangeUser(0.75,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0035);
  osc_param_space->Draw();
  h_dc_123_68->Draw("cont3same");
  h_dc_123_90->Draw("cont3same");
  h_dc_12_68->Draw("cont3same");
  h_dc_12_90->Draw("cont3same");
  h_dc_3_68->Draw("cont3same");
  h_dc_3_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_123_68, "RunI-III const#Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_123_90, "RunI-III const#Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dc_12_68, "RunI+II const#Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_12_90, "RunI+II const#Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dc_3_68, "RunIII const#Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_3_90, "RunIII const#Delta#chi^{2} 90% CL","l");
  leg->Draw();
  c->SaveAs("compare_runs.eps");

  //Run 1-3 28, 2all const dchisq
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->Draw();  
  h_dc_123_68->Draw("cont3same");
  h_dc_123_90->Draw("cont3same");
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_123_68, "2+8 const#Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_123_90, "2+8 const#Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dcall_123_68, "2+N const#Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "2+N const#Delta#chi^{2} 90% CL","l");
  leg->Draw();
  c->SaveAs("compare_2all_28.eps");



}
