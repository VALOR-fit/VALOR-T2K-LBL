void plot_contours_with_profiles
(
 bool plot_fc = false,
 bool plot_constdchisq_lines = false,
 bool plot_kit = false
)
{
  //FC
  TFile f_val_cont("2012results/NEWRESULTS/FINALFCCONTOURS/fc_oldsurf_newres_5E5_nosub.root");
  TFile f_vo2_cont("2012results/FC/fc_contours_2all_octant2.root");
  TFile f_kit_cont("2012results/DATA/KIT/of_fc_kit.root");
  //const dchisq
  TFile f_val_cont_dchi("2012results/NEWRESULTS/datafit/dchisq_2all.root");
  TFile f_vo2_cont_dchi("2012results/DATA/octant2_23only/combination/fitresults_3.root");

  TCanvas * ctemp = new TCanvas();

  double chisq_bf = 56.04;


  TH2D * h_fcall_68 = 0;
  TH2D * h_fcall_90 = 0;
  TH2D * h_fco2_68 = 0;
  TH2D * h_fco2_90 = 0;

  //VaLOR FC
  if(plot_fc) {
    f_val_cont.GetObject("hCL68StatSyst",h_fcall_68);
    f_val_cont.GetObject("hCL90StatSyst",h_fcall_90);
    h_fcall_68->SetContour(1);
    h_fcall_68->SetContourLevel(0,0.0);
    h_fcall_68->SetLineStyle(2);
    h_fcall_90->SetContour(1);
    h_fcall_90->SetContourLevel(0,0.0);
    h_fcall_90->SetLineStyle(1);
    
    //VaLOR FC octant 2
    f_vo2_cont.GetObject("hCL68StatSyst",h_fco2_68);
    f_vo2_cont.GetObject("hCL90StatSyst",h_fco2_90);
    h_fco2_68->SetContour(1);
    h_fco2_68->SetContourLevel(0,0.0);
    h_fco2_68->SetLineStyle(2);
    h_fco2_90->SetContour(1);
    h_fco2_90->SetContourLevel(0,0.0);
    h_fco2_90->SetLineStyle(1);
  }
  else {
    //VaLOR const dchisq
    TH2D * h_fcall = 0;
    f_val_cont_dchi.GetObject("hdchisq",h_fcall);
    h_fcall_68 = new TH2D(*h_fcall);
    h_fcall_90 = new TH2D(*h_fcall);
    h_fcall_68->SetContour(1);
    h_fcall_68->SetContourLevel(0,2.30);
    h_fcall_68->SetLineStyle(2);
    h_fcall_90->SetContour(1);
    h_fcall_90->SetContourLevel(0,4.61);
    h_fcall_90->SetLineStyle(1);

    //VaLOR const dchisq octant 2
    TH2D * h_fco2 = 0;
    f_vo2_cont_dchi.GetObject("hdchisq",h_fco2);
    TH2D * h_fco2_68 = new TH2D(*h_fco2);
    TH2D * h_fco2_90 = new TH2D(*h_fco2);
    h_fco2_68->SetContour(1);
    h_fco2_68->SetContourLevel(0,2.30);
    h_fco2_68->SetLineStyle(2);
    h_fco2_90->SetContour(1);
    h_fco2_90->SetContourLevel(0,4.61);
    h_fco2_90->SetLineStyle(1);
  }

  h_fcall_68->SetLineWidth(2);
  h_fcall_90->SetLineWidth(2);
  h_fcall_68->SetLineColor(kBlack);
  h_fcall_90->SetLineColor(kBlack);
  h_fco2_68->SetLineWidth(2);
  h_fco2_90->SetLineWidth(2);
  h_fco2_68->SetLineColor(kRed);
  h_fco2_90->SetLineColor(kRed);

  //KIT FC
  TH2D * h_fc_kit_68 = 0;
  TH2D * h_fc_kit_90 = 0;
  f_kit_cont.GetObject("hfccontkit68",h_fc_kit_68);
  f_kit_cont.GetObject("hfccontkit90",h_fc_kit_90);
  h_fc_kit_68->SetLineStyle(2);
  h_fc_kit_90->SetLineStyle(1);

  //bestfits
  TMarker bestfit_2all(1.000,0.002443,29);
  bestfit_2all.SetMarkerSize(1.8);
  bestfit_2all.SetMarkerColor(kBlack);
  TMarker bestfit_2o2(0.999,0.002441,24);
  bestfit_2o2.SetMarkerSize(1.8);
  bestfit_2o2.SetMarkerColor(kRed);
  TMarker bestfit_kit(1.000,0.002444,29);
  bestfit_kit.SetMarkerSize(1.8);
  bestfit_kit.SetMarkerColor(kRed);

  //dmsq marg
  TFile f_val_dm("2012results/1DFC/dm23sq/1DFC_dm23sq-run_0-expts_3_3.root");
  TH1D * h_val_dm = 0;
  f_val_dm.GetObject("chisqwithsyst_dm23sq_3",h_val_dm);
  for(int i = 1; i <= h_val_dm->GetNbinsX(); i++) {
     double chisq = h_val_dm->GetBinContent(i);
     h_val_dm->SetBinContent(i, chisq-chisq_bf);
  }

  //mixing marg
  TFile f_val_mix("2012results/1DFC/mixing23/1DFC-run_0-expts_3_3.root");
  TH1D * h_val_mix = 0;
  f_val_mix.GetObject("chisqwithsyst_mixing23_3",h_val_mix);
  for(int i = 1; i <= h_val_mix->GetNbinsX(); i++) {
     double chisq = h_val_mix->GetBinContent(i);
     h_val_mix->SetBinContent(i, chisq-chisq_bf);
  }

  //dmsq marg octant 2
  TFile f_vo2_dm("2012results/1DFC/octant2/dm/octant2dm-run_0-expts_3_3-bins_-1_-1.root");
  TH1D * h_vo2_dm = 0;
  f_vo2_dm.GetObject("chisqwithsyst_dm23sq_3",h_vo2_dm);
  for(int i = 1; i <= h_vo2_dm->GetNbinsX(); i++) {
     double chisq = h_vo2_dm->GetBinContent(i);
     h_vo2_dm->SetBinContent(i, chisq-chisq_bf);
  }

  //mixing marg octant 2
  TFile f_vo2_mix("2012results/1DFC/octant2/mix/octant2mix-run_0-expts_3_3-bins_-1_-1.root");
  TH1D * h_vo2_mix = 0;
  f_vo2_mix.GetObject("chisqwithsyst_mixing23_3",h_vo2_mix);
  for(int i = 1; i <= h_vo2_mix->GetNbinsX(); i++) {
     double chisq = h_vo2_mix->GetBinContent(i);
     h_vo2_mix->SetBinContent(i, chisq-chisq_bf);
  }

  //canvas
  double top_fraction = 0.2;
  double right_fraction = 0.2;
  double left_fraction = 0.16;
  TCanvas * c = new TCanvas();
  c->SetLeftMargin(left_fraction);
  c->cd(1);
  c->SetTopMargin(top_fraction);
  c->SetRightMargin(right_fraction);
  
  // param space
  double mixing23_min   = 0.90;
  double mixing23_max   = 1..;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1.95E-3;
  double dm23sq_max     = 3.15E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int)( ((mixing23_max - mixing23_min) / mixing23_step) + 0.01);
  int dm23sq_nbins   = 1 + (int)( ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  ) + 0.01);
  TH2D * osc_param_space2 = new TH2D("osc_param_space2", "",
				     mixing23_nbins, mixing23_min, mixing23_max,
				     dm23sq_nbins, dm23sq_min, dm23sq_max);
  osc_param_space2->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space2->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  osc_param_space2->GetXaxis()->CenterTitle();
  osc_param_space2->GetYaxis()->CenterTitle();
  osc_param_space2->GetXaxis()->SetTitleOffset(1.25);
  osc_param_space2->GetYaxis()->SetTitleOffset(1.4);
  osc_param_space2->Draw();  

  //draw contour
  h_fcall_68->SetLineColor(kBlack);;
  h_fcall_90->SetLineColor(kBlack);;
  h_fc_kit_68->SetLineColor(kRed);;
  h_fc_kit_90->SetLineColor(kRed);;
  if(plot_kit) {
    h_fc_kit_68->Draw("cont3same");
    h_fc_kit_90->Draw("cont3same");
  }
  h_fcall_68->Draw("cont3same");
  h_fcall_90->Draw("cont3same");
  h_fco2_68 ->Draw("cont3same");
  h_fco2_90 ->Draw("cont3same");
  bestfit_2all.SetMarkerStyle(23);
  bestfit_2all.Draw();
  bestfit_2o2.Draw();
  bestfit_kit.SetMarkerStyle(24);
  bestfit_kit.SetMarkerColor(kRed);
  if(plot_kit)
    bestfit_kit.Draw();

  //legend
  TLegend * leg2 = new TLegend(0.18,0.59,0.41,0.77);
  leg2->SetFillColor(0);
  //leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(h_fcall_68, "T2K (#theta_{23} #leq #pi/4) 68% C.L.","l");
  leg2->AddEntry(h_fcall_90, "T2K (#theta_{23} #leq #pi/4) 90% C.L.","l");
  leg2->AddEntry(h_fco2_68,  "T2K (#theta_{23} #geq #pi/4) 68% C.L.","l");
  leg2->AddEntry(h_fco2_90,  "T2K (#theta_{23} #geq #pi/4) 90% C.L.","l");
  if(plot_kit) {
    leg2->AddEntry(h_fc_kit_68, "T2K 68% C.L. alternative analysis","l");
    leg2->AddEntry(h_fc_kit_90, "T2K 90% C.L. alternative analysis","l");
  }
  leg2->SetX2NDC(0.43);
  leg2->SetTextSize(osc_param_space2->GetXaxis()->GetTitleSize());
  leg2->SetTextSize(osc_param_space2->GetXaxis()->GetTitleSize()*2);
  leg2->SetTextSize(0.05);
  leg2->Draw();

  TLegend * leg3 = new TLegend(0.18, 0.185, 0.41, 0.305);
  leg3->SetFillColor(0);
  //leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextFont(132);
  //bestfit_2all->SetName("bestfit_2all");
  //bestfit_kit->SetName("bestfit_kit");
  leg3->AddEntry(&bestfit_2all, "T2K best fit (#theta_{23} #leq #pi/4)", "p");
  leg3->AddEntry(&bestfit_2o2,  "T2K best fit (#theta_{23} #geq #pi/4)", "p");
  if(plot_kit)
    leg3->AddEntry(&bestfit_kit, "T2K best fit alternative analysis", "p");
  leg3->SetTextSize(0.05);
  leg3->Draw();

  //draw mixing pad
  TPad * pad_top = new TPad("pad_top","pad_top",0,0,1,1);
  pad_top->SetBottomMargin(1-top_fraction);
  pad_top->SetRightMargin(right_fraction);
  pad_top->SetLeftMargin(left_fraction);
  pad_top->SetFillColor(0);
  pad_top->SetFillStyle(0);
  pad_top->Draw();
  pad_top->cd(0);

  //draw mixing frame
  Int_t mix_xnbins  = h_val_mix->GetNbinsX();
  Int_t mix_ynbins  = h_val_mix->GetNbinsY();
  Float_t mix_ymin = 0;
  Float_t mix_ymax = 6;//1.*h_val_mix->GetMaximum();
  Float_t mix_xmin = osc_param_space2->GetXaxis()->GetXmin();
  Float_t mix_xmax = osc_param_space2->GetXaxis()->GetXmax();
  TH2D * mix_frame = new TH2D("","",mix_xnbins,mix_xmin,mix_xmax,mix_ynbins,mix_ymin,mix_ymax);
  mix_frame->GetXaxis()->SetLabelSize(0);
  mix_frame->GetYaxis()->SetLabelFont(132);
  mix_frame->GetYaxis()->SetNdivisions(3,5,0);
  mix_frame->GetYaxis()->SetTitle("#Delta#chi^{2}");
  mix_frame->GetYaxis()->SetTitleOffset(osc_param_space2->GetYaxis()->GetTitleOffset());
  mix_frame->GetYaxis()->CenterTitle();
  mix_frame->GetYaxis()->SetLabelOffset(0.077);
  mix_frame->GetYaxis()->SetTickLength(0.12);
  mix_frame->Draw();

  //draw mixing data octant 2
  TLine *l1, *l2;
  int mix_nbins = h_vo2_mix->GetNbinsX();
  for (Int_t i=1;i<=mix_nbins;i++) {
    Float_t width = 0.5;
    if(i == 1 || i == mix_nbins)
      width = 0.25;
    Float_t mix_ybinmax = h_vo2_mix->GetBinContent(i);
    Float_t mix_ybinmax_previous = h_vo2_mix->GetBinContent(i-1);
    Float_t mix_xbin = h_vo2_mix->GetBinCenter(i);
    Float_t mix_xbinmin = mix_xbin -width*h_vo2_mix->GetBinWidth(i);
    Float_t mix_xbinmax = mix_xbin +width*h_vo2_mix->GetBinWidth(i);
    cout << mix_xbin << "\t" << mix_ybinmax << "\t" << mix_xbinmax << endl;
    if(i == 1) {
      cout << mix_ybinmax << endl;
      //l1 = new TLine(mix_xmin,mix_ybinmax,mix_xbinmax,mix_ybinmax); //bleeds into the y-axis line, and up x-edge is incorrect
      l1 = new TLine(0.900105,mix_ybinmax,0.9005,mix_ybinmax);
      l1->SetLineColor(kRed);
      l1->Draw();
      l1->Print();
    }
    if(mix_xbinmin < mix_xmin) continue;
    if(mix_ybinmax > mix_ymax) continue;
    l1 = new TLine(mix_xbinmin,mix_ybinmax,mix_xbinmax,mix_ybinmax);
    l1->SetLineColor(kRed);
    if(i == 2)
      l1->Print();
    l1->Draw(); //left-right
    double extra = 0.039;
    double extra_high = extra, extra_low = extra;
    if(i == 1)
      extra_high == 0;
    l2 = new TLine(mix_xbinmin,mix_ybinmax_previous+extra_high,mix_xbinmin,mix_ybinmax-extra_low);
    l2->SetLineColor(kRed);
    if(l2->GetY2() > mix_ymax || l2->GetY1() > mix_ymax) continue;
    if(i != mix_nbins)
      l2->Draw(); //up-down
    //cout << mix_xbin << endl;
  }

  //draw mixing data
  int mix_nbins = h_val_mix->GetNbinsX();
  for (Int_t i=1;i<=mix_nbins;i++) {
    Float_t width = 0.5;
    if(i == 1 || i == mix_nbins)
      width = 0.25;
    Float_t mix_ybinmax = h_val_mix->GetBinContent(i);
    Float_t mix_ybinmax_previous = h_val_mix->GetBinContent(i-1);
    Float_t mix_xbin = h_val_mix->GetBinCenter(i);
    Float_t mix_xbinmin = mix_xbin -width*h_val_mix->GetBinWidth(i);
    Float_t mix_xbinmax = mix_xbin +width*h_val_mix->GetBinWidth(i);
    if(mix_xbinmin < mix_xmin)
      continue;
    if(mix_ybinmax > mix_ymax)
      continue;
    l1 = new TLine(mix_xbinmin,mix_ybinmax,mix_xbinmax,mix_ybinmax);
    l1->Draw(); //left-right
    double extra = 0.039;
    double extra_high = extra, extra_low = extra;
    if(i == 1)
      extra_high == 0;
    l2 = new TLine(mix_xbinmin,mix_ybinmax_previous+extra_high,mix_xbinmin,mix_ybinmax-extra_low);
    //l2->SetLineColor(kBlue);
    if(l2->GetY2() > mix_ymax || l2->GetY1() > mix_ymax) continue;
    if(i != mix_nbins)
      l2->Draw(); //up-down
    //cout << mix_xbin << endl;
  }

  //draw mix const dchisq lines
  TLine * l68_mix = new TLine(mix_xmin,1.00,mix_xmax,1.00);
  TLine * l90_mix = new TLine(mix_xmin,2.71,mix_xmax,2.71);
  l68_mix->SetLineStyle(2);
  l90_mix->SetLineStyle(3);
  if(plot_constdchisq_lines) {
    l68_mix->Draw();
    l90_mix->Draw();
  }

  //draw the mix histogram itself, to debug the method
  h_val_mix->SetLabelSize(0);
  h_val_mix->GetXaxis()->SetTitle("");
  h_val_mix->GetXaxis()->SetRangeUser(mixing23_min,mixing23_max);
  //h_val_mix->GetYaxis()->SetRangeUser(0,4);
  h_val_mix->GetYaxis()->SetLabelFont(132);
  h_val_mix->SetLineStyle(2);
  h_val_mix->SetLineWidth(1);
  h_val_mix->SetLineColor(kBlue);
  //h_val_mix->Draw("same");
  //and the octant 2 histogram
  h_vo2_mix->SetLabelSize(0);
  h_vo2_mix->GetXaxis()->SetTitle("");
  //h_vo2_mix->GetXaxis()->SetRangeUser(mixing23_min,mixing23_max);
  //h_vo2_mix->GetYaxis()->SetRangeUser(0,4);
  h_vo2_mix->GetYaxis()->SetLabelFont(132);
  h_vo2_mix->SetLineStyle(1);
  h_vo2_mix->SetLineWidth(1);
  h_vo2_mix->SetLineColor(kBlue);
  //h_vo2_mix->Draw("same");

  //draw dm pad
  TPad * pad_right = new TPad("pad_right","pad_right",0,0,1,1);
  pad_right->SetTopMargin(top_fraction);
  pad_right->SetLeftMargin(1-right_fraction);
  pad_right->SetFillColor(0);
  pad_right->SetFillStyle(0);
  pad_right->Draw();
  pad_right->cd(0);

  //draw dm frame
  Int_t xnbins  = 100;//h_val_dm->GetNbinsX();
  Int_t ynbins  = h_val_dm->GetNbinsY();
  Float_t xmin = 0;
  Float_t xmax = 6;//10;//1.*h_val_dm->GetMaximum();
  Float_t ymin = osc_param_space2->GetYaxis()->GetXmin();
  Float_t ymax = osc_param_space2->GetYaxis()->GetXmax();
  TH2D * dm_frame = new TH2D("","",xnbins,xmin,xmax,ynbins,ymin,ymax);
  dm_frame->GetYaxis()->SetLabelSize(0);
  dm_frame->GetYaxis()->SetLabelFont(132);
  dm_frame->GetXaxis()->SetNdivisions(5,5,0);
  dm_frame->GetXaxis()->SetTitle("#Delta#chi^{2}");
  dm_frame->GetXaxis()->CenterTitle();
  dm_frame->GetXaxis()->SetTitleOffset(1.25);
  dm_frame->GetXaxis()->SetLabelOffset(0.05);
  /*
  for(int i = 1; i <= xnbins; i++) {
    double center = dm_frame->GetXaxis()->GetBinUpEdge(i);
    if(TMath::Abs(center-0) < 1E-8)
      dm_frame->GetXaxis()->SetBinLabel(i,"0");
    else if(TMath::Abs(center-2) < 1E-8)
      dm_frame->GetXaxis()->SetBinLabel(i,"2");
    else if(TMath::Abs(center-4) < 1E-8)
      dm_frame->GetXaxis()->SetBinLabel(i,"4");
    else if(TMath::Abs(center-6) < 1E-8)
      dm_frame->GetXaxis()->SetBinLabel(i,"6");
    else if(TMath::Abs(center-8) < 1E-8)
      dm_frame->GetXaxis()->SetBinLabel(i,"8");
    else if(TMath::Abs(center-10) < 1E-8)
      dm_frame->GetXaxis()->SetBinLabel(i,"10");
    else
      dm_frame->GetXaxis()->SetBinLabel(i,"");
  }
  dm_frame->GetXaxis()->SetNdivisions(5,5,5);
  dm_frame->GetXaxis()->LabelsOption("h");
  dm_frame->GetXaxis()->SetLabelSize(0.07);
  */
  dm_frame->GetXaxis()->SetTickLength(0.12);
  dm_frame->Draw("axis");

  //draw dm data octant 2
  Float_t xbinmin = ymin;
  int nbins = h_vo2_dm->GetNbinsX();
  for (Int_t i=1;i<=nbins;i++) {
    Float_t width = 0.5;
    if(i == 1 || i == nbins)
      width = 0.25;
    Float_t xbinmax = h_vo2_dm->GetBinContent(i);
    Float_t xbinmax_previous = h_vo2_dm->GetBinContent(i-1);
    Float_t ybin = h_vo2_dm->GetBinCenter(i);
    Float_t ybinmin = ybin -width*h_vo2_dm->GetBinWidth(i);
    Float_t ybinmax = ybin +width*h_vo2_dm->GetBinWidth(i);
    if(ybinmin < ymin && ybinmax < ymin) continue;
    if(ybinmin < ymin) ybinmin = ymin;
    l1 = new TLine(xbinmax,ybinmin,xbinmax,ybinmax);
    if(xbinmax > xmax && i > 20) {
      l2 = new TLine(xbinmax_previous-extra,ybinmin,xmax,ybinmin);
      l2->SetLineColor(kRed);
      //l2->Draw();
      break;
    }
    if(xbinmax > xmax) continue;
    l1->SetLineColor(kRed);
    //l1->Draw(); //up-down
    if(ybinmin < ymin) continue;
    double extra = 0.044;
    if(xbinmax_previous < xbinmax)
      l2 = new TLine(xbinmax_previous-extra,ybinmin,xbinmax+extra,ybinmin);
    else
      l2 = new TLine(xbinmax_previous+extra,ybinmin,xbinmax-extra,ybinmin);
    l2->SetLineColor(kRed);
    //l2->Draw(); //left-right
    //cout << ybin << endl;
  }

  //draw dm data
  Float_t xbinmin = ymin;
  int nbins = h_val_dm->GetNbinsX();
  for (Int_t i=1;i<=nbins;i++) {
    Float_t width = 0.5;
    if(i == 1 || i == nbins)
      width = 0.25;
    Float_t xbinmax = h_val_dm->GetBinContent(i);
    Float_t xbinmax_previous = h_val_dm->GetBinContent(i-1);
    Float_t ybin = h_val_dm->GetBinCenter(i);
    Float_t ybinmin = ybin -width*h_val_dm->GetBinWidth(i);
    Float_t ybinmax = ybin +width*h_val_dm->GetBinWidth(i);
    if(ybinmin < ymin && ybinmax < ymin) continue;
    if(ybinmin < ymin) ybinmin = ymin;
    l1 = new TLine(xbinmax,ybinmin,xbinmax,ybinmax);
    if(xbinmax > xmax && i > 20) {
      l2 = new TLine(xbinmax_previous-extra,ybinmin,xmax,ybinmin);
      l2->Draw();
      break;
    }
    if(xbinmax > xmax) continue;
    l1->Draw(); //up-down
    if(ybinmin < ymin) continue;
    double extra = 0.025;
    if(xbinmax_previous < xbinmax)
      l2 = new TLine(xbinmax_previous-extra,ybinmin,xbinmax+extra,ybinmin);
    else
      l2 = new TLine(xbinmax_previous+extra,ybinmin,xbinmax-extra,ybinmin);
    l2->Draw(); //left-right
    //cout << ybin << endl;
  }

  //draw dm const dchisq lines
  TLine * l68_dm = new TLine(1.00,ymin,1.00,ymax);
  TLine * l90_dm = new TLine(2.71,ymin,2.71,ymax);
  l68_dm->SetLineStyle(2);
  l90_dm->SetLineStyle(3);
  //l68_dm->Draw();
  //l90_dm->Draw();

  //draw the markers again (draws them over dm plot)
  c->cd();
  if(plot_kit)
    bestfit_kit.Draw();
  bestfit_2o2.Draw();
  bestfit_2all.Draw();

  if(plot_fc) {
    c->SaveAs("contour_compare_methods2.pdf");
    c->SaveAs("contour_compare_methods2.eps");
  }
  else {
    c->SaveAs("contour_compare_methods2_constdchisq.pdf");
    c->SaveAs("contour_compare_methods2_constdchisq.eps");
  }
  leg2->Clear();
  delete leg2;

  

}
