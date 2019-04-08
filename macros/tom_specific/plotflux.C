{
  TCanvas * c = new TCanvas();

  TFile ft2k("inputs/official/2018/flux/tuned13av3/run1-9a_v1/sk_tuned13av3_13anom_run1-9a_numode_fine.root");
  TH1D * ht2k = 0;
  ft2k.GetObject("enu_sk_tuned13a_numu", ht2k);
  TH1D * ht2k2 = 0;
  ft2k.GetObject("enu_sk_13a_nom_numu", ht2k2);

  TFile fhk("inputs/hk/flux/hk_flux_and_weights_tochi.root");
  TH1D * hhk = 0;
  fhk.GetObject("enu_sk_nomkorea_numode_numu", hhk);
  TH1D * hhk2 = 0;
  fhk.GetObject("enu_sk_tunedkorea_numode_numu", hhk2);

  cout << "T2K tuned   " << ht2k->Integral() << endl;
  cout << "T2K nominal " << ht2k2->Integral() << endl;
  cout << "HK nominal  " << hhk->Integral() << endl;
  cout << "HK tuned    " << hhk2->Integral() << endl;

  double max = (TMath::Max(ht2k->GetMaximum(), hhk->GetMaximum()));
  max = 1.05 * TMath::Max(ht2k2->GetMaximum(), max);
  ht2k->GetXaxis()->SetRangeUser(0,1.5);
  ht2k->GetYaxis()->SetRangeUser(0, max);
  ht2k.Draw();
  ht2k2->SetLineColor(kGreen+2);
  ht2k2->SetMarkerColor(kGreen+2);
  ht2k2.Draw("SAME");
  hhk->SetLineColor(kRed);
  hhk->SetMarkerColor(kRed);
  hhk.Draw("SAME");
  hhk2->SetLineColor(kMagenta);
  hhk2->SetMarkerColor(kMagenta);
  hhk2->SetLineStyle(kSolid);
  hhk2.Draw("SAME LP");

  c->SaveAs("flux.pdf");  
}
