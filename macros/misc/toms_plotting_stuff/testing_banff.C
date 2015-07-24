//plot the 2012 BANFF-tuned nominal & error envelope, then compare with another nominal spectrum
void testing_banff()
{
  TFile fenv("$T2KOSCANA/2012results/THROWS_PAC/sk_ereco_spectrum-twk_sys_throw_sinsq2theta\[12_0.857,13_0.098,23_1.000\],dcp_0.000,dmsq\[12_0.00007500,23_0.002400\]eV2.root");
  TH1D * nominal = 0;
  fenv.GetObject("nominal",nominal);
  TH1D * banffplus = 0, * banffminus = 0;
  fenv.GetObject("BANFF_correlated_envpluserr",banffplus);
  fenv.GetObject("BANFF_correlated_envminuserr",banffminus);
  banffplus->SetLineStyle(2);
  banffminus->SetLineStyle(3);

  //convert +- err to nominal +- err
  for(int i = 1; i <= nominal->GetNbinsX(); i++) {
    double nev = nominal->GetBinContent(i);
    double errp = banffplus->GetBinContent(i);
    double errm = banffminus->GetBinContent(i);
  }

  TFile fnew("testBANFF/13v0/sk_ereco_spectrum-sinsq2theta[12_0.857,13_0.098,23_1.000],dcp_0.000,dmsq[12_0.00007500,23_0.002400]eV2.root");
  TH1D * hnew = 0;
  fnew.GetObject("sk_1Rmutotal",hnew);
  hnew->SetLineColor(kRed);

  TLegend * l = new TLegend(0.7,0.8,0.9,0.9);
  l->AddEntry(nominal,"BANFF12v8","l");
  l->AddEntry(hnew,"BANFF13v0","l");
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  TCanvas * c = new TCanvas();  
  nominal->GetXaxis()->SetRangeUser(0,6.999);
  banffplus->GetXaxis()->SetRangeUser(0,6.999);
  banffplus->Draw();
  nominal->Draw("same");
  banffminus->Draw("same");
  hnew->Draw("same");
  l->Draw();
  c->SaveAs("testbanffenv.pdf");
  banffplus->GetXaxis()->SetRangeUser(0,2.999);
  l->Draw();
  c->SaveAs("testbanffenvzoom.pdf");

  nominal->Draw();
  hnew->Draw("same");
  l->Draw();
  c->SaveAs("testbanff.pdf");

  TH1D * ratio =new TH1D(*nominal);
  ratio->Divide(hnew);
  ratio->GetYaxis()->SetTitle("BANFF12v8 / BANFF13v0");
  ratio->GetXaxis()->SetRangeUser(0.251,6.999);
  ratio->Draw();
  c->SaveAs("testbanffratio.pdf");
}
