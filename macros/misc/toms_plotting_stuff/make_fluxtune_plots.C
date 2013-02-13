void make_fluxtune_plots
(
 double erecomax = 10.0,
 const char * filename = "inputs/official/2012/flux/11b_v3.2_run1-3c/sk_tuned11bv3.2_11anom_run1-run3c_fine.root",
 const char * prod     = "11b"
)
{

  TCanvas * c = new TCanvas();

  TFile f(filename);
  TH1D *hnumu = 0, *hnumub = 0, *hnue = 0, *hnueb = 0;
  f.GetObject(Form("enu_sk_tuned%s_numu_ratio",prod),hnumu);
  f.GetObject(Form("enu_sk_tuned%s_numub_ratio",prod),hnumub);
  f.GetObject(Form("enu_sk_tuned%s_nue_ratio",prod),hnue);
  f.GetObject(Form("enu_sk_tuned%s_nueb_ratio",prod),hnueb);

  hnumu->GetXaxis()->SetRangeUser(0, erecomax - 0.001);
  hnumub->GetXaxis()->SetRangeUser(0, erecomax - 0.001);
  hnue->GetXaxis()->SetRangeUser(0, erecomax - 0.001);
  hnueb->GetXaxis()->SetRangeUser(0, erecomax - 0.001);

  for(int i = 1; i <= hnumu->GetNbinsX(); i++) {
    hnumu->SetBinError(i,0);
    hnumub->SetBinError(i,0);
    hnue->SetBinError(i,0);
    hnueb->SetBinError(i,0);
  }

  hnumu->Draw("");
  c->SaveAs("fluxtune_numu.eps");
  hnumub->Draw("");
  c->SaveAs("fluxtune_numub.eps");
  hnue->Draw("");
  c->SaveAs("fluxtune_nue.eps");
  hnueb->Draw("");
  c->SaveAs("fluxtune_nueb.eps");

}
