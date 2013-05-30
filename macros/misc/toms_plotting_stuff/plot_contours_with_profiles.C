void plot_contours_with_profiles()
{
  TFile f_val_cont("2012results/NEWRESULTS/FINALFCCONTOURS/fc_oldsurf_newres_5E5_nosub.root");
  TFile f_kit_cont("2012results/DATA/KIT/of_fc_kit.root");


  double chisq_bf = 56.04;

  //VaLOR FC
  TH2D * h_fcall_68 = 0;
  TH2D * h_fcall_90 = 0;
  f_val_cont.GetObject("hCL68StatSyst",h_fcall_68);
  f_val_cont.GetObject("hCL90StatSyst",h_fcall_90);
  h_fcall_68->SetContour(1);
  h_fcall_68->SetContourLevel(0,0.0);
  h_fcall_68->SetLineStyle(2);
  h_fcall_90->SetContour(1);
  h_fcall_90->SetContourLevel(0,0.0);
  h_fcall_90->SetLineStyle(1);

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
  double mixing23_min   = 0.92;
  double mixing23_max   = 1..;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 4E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + TMath::CeilNint ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + TMath::CeilNint ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  TH2D * osc_param_space2 = new TH2D("osc_param_space2", "",
				     mixing23_nbins, mixing23_min, mixing23_max,
				     dm23sq_nbins, 0.00205, 0.00305);
  osc_param_space2->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space2->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  //osc_param_space2->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space2->GetYaxis()->SetTitleOffset(1.4);
  osc_param_space2->Draw();  

  //draw contour
  h_fcall_68->SetLineColor(kBlack);;
  h_fcall_90->SetLineColor(kBlack);;
  h_fc_kit_68->SetLineColor(kRed);;
  h_fc_kit_90->SetLineColor(kRed);;
  h_fc_kit_68->Draw("cont3same");
  h_fc_kit_90->Draw("cont3same");
  h_fcall_68->Draw("cont3same");
  h_fcall_90->Draw("cont3same");
  bestfit_2all.Draw();
  bestfit_kit.SetMarkerStyle(24);
  bestfit_kit.SetMarkerColor(kRed);
  bestfit_kit.Draw();

  //legend
  TLegend * leg2 = new TLegend(0.17,0.7,0.5,0.92);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(h_fcall_68, "T2K 68% C.L.","l");
  leg2->AddEntry(h_fcall_90, "T2K 90% C.L.","l");
  leg2->AddEntry(h_fc_kit_68, "T2K 68% C.L. alternative analysis","l");
  leg2->AddEntry(h_fc_kit_90, "T2K 90% C.L. alternative analysis","l");
  leg2->SetX2NDC(0.43);
  leg2->SetTextSize(osc_param_space2->GetXaxis()->GetTitleSize());
  leg2->SetTextSize(osc_param_space2->GetXaxis()->GetTitleSize()*2);
  leg2->SetTextSize(0.05);
  leg2->Draw();

  //draw 1d mixing
  TPad * pad_top = new TPad("pad_top","pad_top",0,0,1,1);
  pad_top->SetBottomMargin(1-top_fraction);
  pad_top->SetRightMargin(right_fraction);
  pad_top->SetLeftMargin(left_fraction);
  pad_top->SetFillColor(0);
  pad_top->SetFillStyle(0);
  pad_top->Draw();
  pad_top->cd(0);
  h_val_mix->SetLabelSize(0);
  h_val_mix->GetXaxis()->SetTitle("");
  h_val_mix->GetXaxis()->SetRangeUser(mixing23_min,mixing23_max);
  //h_val_mix->GetYaxis()->SetRangeUser(0,4);
  h_val_mix->GetYaxis()->SetLabelFont(132);
  h_val_mix->Draw("");

  //draw 1D dm
  TPad * pad_right = new TPad("pad_right","pad_right",0,0,1,1);
  pad_right->SetTopMargin(top_fraction);
  pad_right->SetLeftMargin(1-right_fraction);
  pad_right->SetFillColor(0);
  pad_right->SetFillStyle(0);
  pad_right->Draw();
  pad_right->cd(0);
  h_val_dm->SetLabelSize(0);
  h_val_dm->GetXaxis()->SetTitle("");
  h_val_dm->GetXaxis()->SetRangeUser(dm23sq_min,dm23sq_max);

  //TCanvas * ctemp = new TCanvas();
  h_val_dm->Draw("");
  /*
  TImage *img = TImage::Create();
  img->FromPad(ctemp);
  img->Flip(90);
  pad_right->cd();
  img->Draw("x");
  pad_right->Update();
  */
  c->SaveAs("contour_compare_methods2.pdf");
  leg2->Clear();
  delete leg2;

  

}
