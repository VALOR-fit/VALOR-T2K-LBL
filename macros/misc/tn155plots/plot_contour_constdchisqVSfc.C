//plots an FC contour with a constant dchisq contour
//also produces an FC-only plot
void plot_contour_constdchisqVSfc
(
 const char * tag        = "normal",
 const char * dcfilename = "2013results/DATA/run1-4full_skdetmathfix_run1-4abanff_nuis4_sintheta_normal/combination/run1-4full_skdetmathfix_run1-4abanff_nuis4_sintheta_normal-run_0-expt_8.root", //"2013results/DATA/run1-4_nuis3fixdcp_sintheta_normal/combination/run1-4_nuis3fixdcp_sintheta_normal-run_0-expt_6.root",
 const char * fcfilename = "2013results/FC/fc_contour_normal_nuis3fixdcp.root",
 bool draw_fc            = false
 )
{
  using namespace t2kana;

  double linewidth = 2;
  double marker_size = 1.8;
  double legend_text_size = 0.05;

  //get dc contour
  TFile fdc(dcfilename);
  TH2D * hdc68 = 0, * hdc90 = 0;
  fdc.GetObject("dchisq_cont68",hdc68);
  fdc.GetObject("dchisq_cont90",hdc90);
  hdc68->SetLineWidth(2);
  hdc90->SetLineWidth(2);
  hdc68->SetLineStyle(2);
  hdc90->Smooth();
  hdc68->Smooth();

  //get fc contour
  TFile ffc(fcfilename);
  TH2D * hfc68 = 0, * hfc90 = 0;
  ffc.GetObject("hCL68StatSyst",hfc68);
  ffc.GetObject("hCL90StatSyst",hfc90);
  hfc68->SetLineWidth(2);
  hfc90->SetLineWidth(2);
  hfc68->SetLineStyle(2);
  //get rid of crap at the bottom of the plot
  hfc68->GetYaxis()->SetRangeUser(0.002,0.003);
  hfc90->GetYaxis()->SetRangeUser(0.002,0.003);

  //get bestfit point
  FitParamList * fparm = 0;
  fdc.GetObject("fparm",fparm);
  double dm32sq      = fparm->BestFitValue(conv::kPhypDm23sq);
  double mix23       = fparm->BestFitValue(conv::kPhypMixing23);
  TMarker bf(mix23, dm32sq, 29);
  bf.SetMarkerSize(marker_size);
  bf.SetMarkerColor(kBlack);

  //get externals
  TFile fext("inputs/contours/contours_nont2k.root");
  TGraph * g_minos = 0;
  fext.GetObject("minos2013_combined_90CL_normal_sinsqtheta",g_minos);
  g_minos->SetLineStyle(1);
  g_minos->SetName("g_minos");
  TGraph * g_sk = 0;
  fext.GetObject("sk2012_zenith3nu_90CL_normal_sinsqtheta",g_sk);
  g_sk->SetLineStyle(1);
  g_sk->SetName("g_sk");

  //get old t2k
  TFile fold("2012results/DATA/sintheta23/combine/combination/fitresults_3.root");
  TH2D * h_t2k123 = 0;
  fold.GetObject("hdchisq",h_t2k123);
  TH2D * h_t2k123_68 = new TH2D(*h_t2k123);
  TH2D * h_t2k123_90 = new TH2D(*h_t2k123);
  h_t2k123_68->SetContour(1);
  h_t2k123_68->SetContourLevel(0,2.30);
  h_t2k123_68->SetLineStyle(3);
  h_t2k123_90->SetContour(1);
  h_t2k123_90->SetContourLevel(0,4.61);
  h_t2k123_90->SetLineStyle(2);

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
  int mixing23_nbins = 1 + (int)( ((mixing23_max - mixing23_min) / mixing23_step) + 0.01);
  int dm23sq_nbins   = 1 + (int)( ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  ) + 0.01);
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

  //color
  hdc68->SetLineColor(kRed);
  hdc90->SetLineColor(kRed);
  hfc68->SetLineColor(kBlack);
  hfc90->SetLineColor(kBlack);
  h_t2k123_68->SetLineColor(kBlack);
  h_t2k123_90->SetLineColor(kBlack);
  g_minos->SetLineColor(kBlue);
  g_sk->SetLineColor(kGreen+2);

  //Draw constdchisq vs FC
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(0.4,0.65);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.0032);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  hdc68->Draw("cont3same");
  hdc90->Draw("cont3same");
  hfc68->Draw("cont3same");
  hfc90->Draw("cont3same");
  bf.Draw();
  TLegend * leg2 = new TLegend(0.185,0.7,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(hfc68, "Run 1-4a FC 68% CL","l");
  leg2->AddEntry(hfc90, "Run 1-4a FC 90% CL","l");
  leg2->AddEntry(hdc68, "Run 1-4a const-#Delta#chi^{2} 68% CL","l");
  leg2->AddEntry(hdc90, "Run 1-4a const-#Delta#chi^{2} 90% CL","l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(legend_text_size);
  leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.17, 0.45, 0.30);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bf, "Run 1-4a best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(legend_text_size);
  legbf->Draw();
  c->SaveAs(Form("contour_compare_%s_dcvsfc.pdf",tag));
  c->SaveAs(Form("contour_compare_%s_dcvsfc.png",tag));
  c->SaveAs(Form("contour_compare_%s_dcvsfc.eps",tag));

  //Draw constdchisq vs FC + externals
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(0.3,0.7);
  osc_param_space->GetYaxis()->SetRangeUser(0.0017,0.0045);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  g_sk->Draw("l");
  g_minos->Draw("l");
  h_t2k123_90->Draw("cont3same");
  //hdc68->Draw("cont3same");
  hdc90->Draw("cont3same");
  //hfc68->Draw("cont3same");
  if(draw_fc)
    hfc90->Draw("cont3same");
  bf.Draw();
  TLegend * leg2 = new TLegend(0.185,0.64,0.45,0.96);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  //leg2->AddEntry(hfc68, "T2K FC Run 1-4 (#Deltam^{2}_{32}>0) 68% CL","l");
  if(draw_fc)
    leg2->AddEntry(hfc90, "T2K FC Run 1-4 (#Deltam^{2}_{32}>0) 90% CL","l");
  //leg2->AddEntry(hdc68, "T2K const-#Delta#chi^{2} Run1-4 (#Deltam^{2}_{32}>0) 68% CL","l");
  if(draw_fc)
    leg2->AddEntry(hdc90, "T2K const-#Delta#chi^{2} Run1-4 (#Deltam^{2}_{32}>0) 90% CL","l");
  else
    leg2->AddEntry(hdc90, "T2K Run1-4 (#Deltam^{2}_{32}>0) 90% CL","l");
  leg2->AddEntry(h_t2k123_90, "T2K Run 1-3 (#Deltam^{2}_{32}>0) 90% CL","l");
  leg2->AddEntry(g_minos, "MINOS 2013 (#Deltam^{2}_{32}>0) 90% CL","l");
  leg2->AddEntry(g_sk, "SK 2012 (#Deltam^{2}_{32}>0) 90% CL","l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(legend_text_size);
  leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.14, 0.45, 0.27);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bf, "T2K 3#nu (#Deltam^{2}_{32}>0) best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(legend_text_size);
  legbf->Draw();
  c->SaveAs(Form("contour_compare_%s_dcvsfc_ext.pdf",tag));
  c->SaveAs(Form("contour_compare_%s_dcvsfc_ext.png",tag));
  c->SaveAs(Form("contour_compare_%s_dcvsfc_ext.eps",tag));

  //Draw FC contour only
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(0.4,0.65);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.0032);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  hfc68->Draw("cont3same");
  hfc90->Draw("cont3same");
  bf.Draw();
  TLegend * leg2 = new TLegend(0.185,0.8,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(hfc68, "Run 1-4a FC 68% CL","l");
  leg2->AddEntry(hfc90, "Run 1-4a FC 90% CL","l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(legend_text_size);
  leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.17, 0.45, 0.30);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bf, "Run 1-4a best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(legend_text_size);
  legbf->Draw();
  c->SaveAs(Form("contour_compare_%s_fconly.pdf",tag));
  c->SaveAs(Form("contour_compare_%s_fconly.png",tag));
  c->SaveAs(Form("contour_compare_%s_fconly.eps",tag));

}

