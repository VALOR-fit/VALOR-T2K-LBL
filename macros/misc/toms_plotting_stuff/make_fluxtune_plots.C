void make_fluxtune_plots() {

  TCanvas * c = new TCanvas();

  TFile f("inputs/official/2012/flux/11b_v3.2_run1-3c/sk_tuned11bv3.2_11anom_run1-run3c_fine.root");
  TH1D *hnumu = 0, *hnumub = 0, *hnue = 0, *hnueb = 0;
  f.GetObject("enu_sk_tuned11b_numu_ratio",hnumu);
  f.GetObject("enu_sk_tuned11b_numub_ratio",hnumub);
  f.GetObject("enu_sk_tuned11b_nue_ratio",hnue);
  f.GetObject("enu_sk_tuned11b_nueb_ratio",hnueb);

  hnumu->Draw("P");
  c->SaveAs("fluxtune_numu.eps");
  hnumub->Draw("P");
  c->SaveAs("fluxtune_numub.eps");
  hnue->Draw("P");
  c->SaveAs("fluxtune_nue.eps");
  hnueb->Draw("P");
  c->SaveAs("fluxtune_nueb.eps");

 
}
