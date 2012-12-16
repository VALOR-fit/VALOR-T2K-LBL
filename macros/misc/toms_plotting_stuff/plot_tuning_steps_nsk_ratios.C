void plot_tuning_steps_nsk_ratios
(

 )
{
  using namespace t2kana;

  TFile f_tun("alltune/nsk_vs_osc_params.root");
  TFile f_flux("fluxtune/nsk_vs_osc_params.root");
  TFile f_nom("notune/nsk_vs_osc_params.root");
  TFile f_inv("inverted/nsk_vs_osc_params.root");
  TFile f_2f("2flavour/nsk_vs_osc_params.root");

  //histograms
  TH2D * h_tun = 0;
  f_tun.GetObject("nsk_vs_osc_params",h_tun);
  TH2D * h_flux = 0;
  f_flux.GetObject("nsk_vs_osc_params",h_flux);
  TH2D * h_nom = 0;
  f_nom.GetObject("nsk_vs_osc_params",h_nom);
  TH2D * h_inv = 0;
  f_inv.GetObject("nsk_vs_osc_params",h_inv);
  TH2D * h_2f = 0;
  f_2f.GetObject("nsk_vs_osc_params",h_2f);

  //ratios
  TH2D * h_ratio_tunflux = new TH2D(*h_tun);
  TH2D * h_ratio_nomflux = new TH2D(*h_tun);
  TH2D * h_ratio_invtun = new TH2D(*h_tun);
  TH2D * h_ratio_2ftun = new TH2D(*h_tun);
  for(int i = 1; i <= h_tun->GetNbinsX(); i++) {
    for(int j = 1; j <= h_tun->GetNbinsY(); j++) {
      double n_tun = h_tun->GetBinContent(i,j);
      double n_flux = h_flux->GetBinContent(i,j);
      double n_nom = h_nom->GetBinContent(i,j);
      double n_inv = h_inv->GetBinContent(i,j);
      double n_2f = h_2f->GetBinContent(i,j);
      h_ratio_tunflux->SetBinContent(i,j, (n_tun-n_flux)/n_flux);
      h_ratio_nomflux->SetBinContent(i,j, (n_flux-n_nom)/n_nom);
      h_ratio_invtun->SetBinContent(i,j, (n_inv-n_tun)/n_tun);
      h_ratio_2ftun->SetBinContent(i,j, (n_2f-n_tun)/n_tun);
    }
  }


  TCanvas * c = new TCanvas();
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  h_ratio_tunflux->GetYaxis()->SetTitleOffset(1.3);
  h_ratio_nomflux->GetYaxis()->SetTitleOffset(1.3);
  h_tun->GetYaxis()->SetTitleOffset(1.3);
  h_ratio_invtun->GetYaxis()->SetTitleOffset(1.3);
  h_ratio_2ftun->GetYaxis()->SetTitleOffset(1.3);

  h_ratio_tunflux->Draw("colz");
  c->SaveAs("nsk_diff_tunflux.eps");
  h_ratio_nomflux->Draw("colz");
  c->SaveAs("nsk_diff_nomflux.eps");
  h_tun->Draw("colz");
  c->SaveAs("nsk_nominalbanff.eps");
  h_ratio_invtun->Draw("colz");
  c->SaveAs("nsk_diff_inv.eps");
  h_ratio_2ftun->Draw("colz");
  c->SaveAs("nsk_diff_2f.eps");

}
