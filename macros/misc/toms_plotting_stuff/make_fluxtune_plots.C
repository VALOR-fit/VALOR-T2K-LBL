void make_fluxtune_plots
(
 double erecomax = 5.0
)
{

  TCanvas * c = new TCanvas();

  TFile f("inputs/official/2012/flux/11b_v3.2_run1-3c/sk_tuned11bv3.2_11anom_run1-run3c_fine.root");
  TH1D *hnumu = 0, *hnumub = 0, *hnue = 0, *hnueb = 0;
  f.GetObject("enu_sk_tuned11b_numu_ratio",hnumu);
  f.GetObject("enu_sk_tuned11b_numub_ratio",hnumub);
  f.GetObject("enu_sk_tuned11b_nue_ratio",hnue);
  f.GetObject("enu_sk_tuned11b_nueb_ratio",hnueb);

  hnumu->GetXaxis()->SetRangeUser(0, erecomax - 0.001);
  hnumub->GetXaxis()->SetRangeUser(0, erecomax - 0.001);
  hnue->GetXaxis()->SetRangeUser(0, erecomax - 0.001);
  hnueb->GetXaxis()->SetRangeUser(0, erecomax - 0.001);

  hnumu->SetMarkerStyle(6);
  hnumub->SetMarkerStyle(6);
  hnue->SetMarkerStyle(6);
  hnueb->SetMarkerStyle(6);

  hnumu->Draw("P");
  c->SaveAs("fluxtune_numu.eps");
  hnumub->Draw("P");
  c->SaveAs("fluxtune_numub.eps");
  hnue->Draw("P");
  c->SaveAs("fluxtune_nue.eps");
  hnueb->Draw("P");
  c->SaveAs("fluxtune_nueb.eps");

 
}
