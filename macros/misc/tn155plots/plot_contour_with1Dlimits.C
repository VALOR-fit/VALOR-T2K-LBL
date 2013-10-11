// plot sinsqtheta contour, with & without external results
//also add the 1D limits as lines (request by Chang Kee, if I understood him correctly)
//

void plot_contour_with1Dlimits() {
  using namespace t2kana;

  const char * filename = "2013results/DATA/run1-4full_skdetmathfix_run1-4abanff_nuis4_sintheta_normal/combination/run1-4full_skdetmathfix_run1-4abanff_nuis4_sintheta_normal-run_0-expt_8.root";
  const char * tag = "normal";
  bool external = false;
  TFile f(filename);

  //style constants
  double marker_size = 1.8;
  double contour_width = 2;
  double limit_width = 1;
  double legend_text_size = 0.06;

  external = true;

  //edge of canvas
  double ymin = 0.002;
  double ymax = 0.003;
  double xmin = 0.39;
  double xmax = 0.64;
  if(external) {
    double ymin = 0.002;
    double ymax = 0.0035;
    double xmin = 0.33;
    double xmax = 0.7;
  }

  //get contours
  TH2D * h68 = 0, * h90 = 0;
  f.GetObject("dchisq_cont68",h68);
  h68->SetLineWidth(contour_width);
  h68->SetLineColor(kBlack);
  h68->SetLineStyle(2);
  f.GetObject("dchisq_cont90",h90);
  h90->SetLineWidth(contour_width);
  h90->SetLineColor(kBlack);
  h90->SetLineStyle(1);

  h68->Smooth();
  h90->Smooth();

  //get bestfit point
  FitParamList * fparm = 0;
  f.GetObject("fparm",fparm);
  double dm32sq      = fparm->BestFitValue(conv::kPhypDm23sq);
  double mix23       = fparm->BestFitValue(conv::kPhypMixing23);
  TMarker bf(mix23, dm32sq, 29);
  bf.SetMarkerSize(marker_size);
  bf.SetMarkerColor(kBlack);

  //get limits
  double mix_lo = 0.454;
  double mix_hi = 0.574;
  double dm_lo  = 0.00237;
  double dm_hi  = 0.00262;

  //get externals
  TFile fext("inputs/contours/contours_nont2k.root");
  TGraph * g_minos = 0;
  fext.GetObject("minos2013_combined_90CL_normal_sinsqtheta",g_minos);
  g_minos->SetLineStyle(1);
  g_minos->SetName("g_minos");
  g_minos->SetLineColor(kBlue);
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


  //param space
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


  //draw
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(xmin,xmax);
  osc_param_space->GetYaxis()->SetRangeUser(ymin,ymax);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  h68->Draw("cont3same");
  h90->Draw("cont3same");
  if(external) {
    g_sk->Draw("l");
    g_minos->Draw("l");
    h_t2k123_90->Draw("cont3same");
  }
  bf.Draw();
  TLine * l1 = new TLine(xmin,dm_lo,xmax,dm_lo);
  TLine * l2 = new TLine(xmin,dm_hi,xmax,dm_hi);
  TLine * l3 = new TLine(mix_lo,ymin,mix_lo,ymax);
  TLine * l4 = new TLine(mix_hi,ymin,mix_hi,ymax);
  l1->SetLineWidth(limit_width);
  l2->SetLineWidth(limit_width);
  l3->SetLineWidth(limit_width);
  l4->SetLineWidth(limit_width);
  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  l3->SetLineStyle(2);
  l4->SetLineStyle(2);
  l1->Draw();
  l2->Draw();
  l3->Draw();
  l4->Draw();
  TLegend * leg2 = new TLegend(0.185,0.7,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(h68, "Run 1-4 const-#Delta#chi^{2} 68% CL","l");
  leg2->AddEntry(h90, "Run 1-4 const-#Delta#chi^{2} 90% CL","l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(legend_text_size);
  //leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.17, 0.45, 0.30);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bf, "Run 1-4 best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(legend_text_size);
  //legbf->Draw();
  if(!external) {
    c->SaveAs(Form("contour_with_1dlimits_%s.pdf",tag));
    c->SaveAs(Form("contour_with_1dlimits_%s.eps",tag));
  }
  else {
   c->SaveAs(Form("contour_with_1dlimits_%s_ext.pdf",tag));
   c->SaveAs(Form("contour_with_1dlimits_%s_ext.eps",tag));
  }
}
