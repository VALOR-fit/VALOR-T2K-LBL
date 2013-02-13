void oneovermaqe2
(
 int rebin = 2
)
{
  TFile fm("maqe/sk_ereco_spectrum-twk_sys_throw_sinsq2theta[12_0.857,13_0.098,23_1.000],dcp_0.000,dmsq[12_0.00007500,23_0.002400]eV2.root");
  TFile fo("oneovermaqe2/sk_ereco_spectrum-twk_sys_throw_sinsq2theta[12_0.857,13_0.098,23_1.000],dcp_0.000,dmsq[12_0.00007500,23_0.002400]eV2.root");

  TH1D * hm = 0;
  fm.GetObject("BANFF_correlated_nsk",hm);
  TH1D * ho = 0;
  fo.GetObject("BANFF_correlated_nsk",ho);

  TCanvas * c = new TCanvas();
  c->Divide(2,1);
  gStyle->SetOptFit(1);

  hm->Rebin(rebin);
  ho->Rebin(rebin);

  hm->GetXaxis()->SetRangeUser(45,70);
  ho->GetXaxis()->SetRangeUser(45,70);

  hm->GetXaxis()->SetTitle("N_{SK}");
  ho->GetXaxis()->SetTitle("N_{SK}");
  hm->GetYaxis()->SetTitle("");


  hm->SetLineColor(kRed);
  ho->SetLineColor(kRed);

  c->cd(1);
  hm->Draw();
  hm->Fit("gaus","+");
 
  c->cd(2);
  ho->Draw();
  ho->Fit("gaus","+");



  c->SaveAs("oneovermaqe2test.eps");

}
