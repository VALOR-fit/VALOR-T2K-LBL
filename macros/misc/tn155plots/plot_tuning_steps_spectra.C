void plot_tuning_steps_spectra
(
 double erecomax = 7.0,
 double templemax = 10.0
)
{
  using namespace t2kana;

  TFile f_tun_noosc ("2013results/SPECTRA/FULLRUN1-4FIXSPLINE/tuning_steps/banff_tuned/sk_ereco_spectrum-sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00000000,23_0.000000]eV2.root");
  TFile f_tun_osc ("2013results/SPECTRA/FULLRUN1-4FIXSPLINE/tuning_steps/banff_tuned/sk_ereco_spectrum-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.root");
  TFile f_flux_noosc ("2013results/SPECTRA/FULLRUN1-4FIXSPLINE/tuning_steps/flux_tuned/sk_ereco_spectrum-sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00000000,23_0.000000]eV2.root");
  TFile f_flux_osc ("2013results/SPECTRA/FULLRUN1-4FIXSPLINE/tuning_steps/flux_tuned/sk_ereco_spectrum-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.root");
  TFile f_nom_noosc ("2013results/SPECTRA/FULLRUN1-4FIXSPLINE/tuning_steps/nominal/sk_ereco_spectrum-sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00000000,23_0.000000]eV2.root");
  TFile f_nom_osc ("2013results/SPECTRA/FULLRUN1-4FIXSPLINE/tuning_steps/nominal/sk_ereco_spectrum-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.root");

  TH1D * h_ereco_tun_noosc = 0;
  f_tun_noosc.GetObject("sk_1Rmutotal",h_ereco_tun_noosc);
  TH1D * h_ereco_tun_osc = 0;
  f_tun_osc.GetObject("sk_1Rmutotal",h_ereco_tun_osc);
  TH1D * h_ereco_flux_noosc = 0;
  f_flux_noosc.GetObject("sk_1Rmutotal",h_ereco_flux_noosc);
  TH1D * h_ereco_flux_osc = 0;
  f_flux_osc.GetObject("sk_1Rmutotal",h_ereco_flux_osc);
  TH1D * h_ereco_nom_noosc = 0;
  f_nom_noosc.GetObject("sk_1Rmutotal",h_ereco_nom_noosc);
  TH1D * h_ereco_nom_osc = 0;
  f_nom_osc.GetObject("sk_1Rmutotal",h_ereco_nom_osc);

  //change to events per 50 MeV
  for(int ix = 1; ix <= h_ereco_tun_noosc->GetNbinsX(); ix++) {
    h_ereco_tun_noosc ->SetBinContent(ix, (h_ereco_tun_noosc ->GetBinContent(ix) / h_ereco_tun_noosc ->GetBinWidth(ix)) * 0.05);
    h_ereco_tun_osc   ->SetBinContent(ix, (h_ereco_tun_osc   ->GetBinContent(ix) / h_ereco_tun_osc   ->GetBinWidth(ix)) * 0.05);
    h_ereco_flux_noosc->SetBinContent(ix, (h_ereco_flux_noosc->GetBinContent(ix) / h_ereco_flux_noosc->GetBinWidth(ix)) * 0.05);
    h_ereco_flux_osc  ->SetBinContent(ix, (h_ereco_flux_osc  ->GetBinContent(ix) / h_ereco_flux_osc  ->GetBinWidth(ix)) * 0.05);
    h_ereco_nom_noosc ->SetBinContent(ix, (h_ereco_nom_noosc ->GetBinContent(ix) / h_ereco_nom_noosc ->GetBinWidth(ix)) * 0.05);
    h_ereco_nom_osc   ->SetBinContent(ix, (h_ereco_nom_osc   ->GetBinContent(ix) / h_ereco_nom_osc   ->GetBinWidth(ix)) * 0.05);
  }//ix

  //colors
  h_ereco_tun_noosc->SetLineColor(kBlack);
  h_ereco_tun_osc->SetLineColor(kBlack);
  h_ereco_flux_noosc->SetLineColor(kRed);
  h_ereco_flux_osc->SetLineColor(kRed);
  h_ereco_nom_noosc->SetLineColor(kBlue);
  h_ereco_nom_osc->SetLineColor(kBlue);

  //width
  h_ereco_tun_noosc->SetLineWidth(1);
  h_ereco_tun_osc->SetLineWidth(1);
  h_ereco_flux_noosc->SetLineWidth(1);
  h_ereco_flux_osc->SetLineWidth(1);
  h_ereco_nom_noosc->SetLineWidth(1);
  h_ereco_nom_osc->SetLineWidth(1);

  //range
  h_ereco_tun_noosc->GetXaxis()->SetRangeUser(0,erecomax-conv::kOneMeV);
  h_ereco_tun_osc->GetXaxis()->SetRangeUser(0,erecomax-conv::kOneMeV);
  h_ereco_flux_noosc->GetXaxis()->SetRangeUser(0,erecomax-conv::kOneMeV);
  h_ereco_flux_osc->GetXaxis()->SetRangeUser(0,erecomax-conv::kOneMeV);
  h_ereco_nom_noosc->GetXaxis()->SetRangeUser(0,erecomax-conv::kOneMeV);
  h_ereco_nom_osc->GetXaxis()->SetRangeUser(0,erecomax-conv::kOneMeV);

  //titles
  h_ereco_tun_noosc->GetXaxis()->SetTitle("E_{reco} / GeV");
  h_ereco_tun_osc->GetXaxis()->SetTitle("E_{reco} / GeV");
  h_ereco_flux_noosc->GetXaxis()->SetTitle("E_{reco} / GeV");
  h_ereco_flux_osc->GetXaxis()->SetTitle("E_{reco} / GeV");
  h_ereco_nom_noosc->GetXaxis()->SetTitle("E_{reco} / GeV");
  h_ereco_nom_osc->GetXaxis()->SetTitle("E_{reco} / GeV");
  h_ereco_tun_noosc->GetYaxis()->SetTitle("Events per 0.050 GeV");
  h_ereco_tun_osc->GetYaxis()->SetTitle("Events per 0.050 GeV");
  h_ereco_flux_noosc->GetYaxis()->SetTitle("Events per 0.050 GeV");
  h_ereco_flux_osc->GetYaxis()->SetTitle("Events per 0.050 GeV");
  h_ereco_nom_noosc->GetYaxis()->SetTitle("Events per 0.050 GeV");
  h_ereco_nom_osc->GetYaxis()->SetTitle("Events per 0.050 GeV");

  TCanvas * c = new TCanvas();
  c->Divide(2,2);
  TLegend * leg = new TLegend(0.25,0.9,0.88,0.7);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  //leg->SetTextFont(42);
  leg->SetTextSize(0.05);

  //nom vs flux tuned
  leg->Clear();
  leg->AddEntry(h_ereco_nom_osc,"Nominal","l");
  leg->AddEntry(h_ereco_flux_osc,"Flux tuning applied","l");
  c->cd(1);
  h_ereco_flux_noosc->Draw();
  h_ereco_nom_noosc->Draw("same");
  leg->Draw();
  c->cd(2);
  h_ereco_flux_osc->Draw();
  h_ereco_nom_osc->Draw("same");
  leg->Draw();
  c->cd(3);
  TH1D * h_ereco_nomflux_noosc_ratio = new TH1D(*h_ereco_flux_noosc);
  h_ereco_nomflux_noosc_ratio->Divide(h_ereco_flux_noosc,h_ereco_nom_noosc);
  h_ereco_nomflux_noosc_ratio->SetLineColor(kBlack);
  h_ereco_nomflux_noosc_ratio->GetYaxis()->SetTitle("Flux tuned / Nominal");
  h_ereco_nomflux_noosc_ratio->GetYaxis()->SetRangeUser(0.9,1.4);
  h_ereco_nomflux_noosc_ratio->Draw();
  c->cd(4);
  TH1D * h_ereco_nomflux_osc_ratio = new TH1D(*h_ereco_flux_osc);
  h_ereco_nomflux_osc_ratio->Divide(h_ereco_flux_osc,h_ereco_nom_osc);
  h_ereco_nomflux_osc_ratio->SetLineColor(kBlack);
  h_ereco_nomflux_osc_ratio->GetYaxis()->SetTitle("Flux tuned / Nominal");
  h_ereco_nomflux_osc_ratio->GetYaxis()->SetRangeUser(0.9,1.4);
  h_ereco_nomflux_osc_ratio->Draw();
  c->SaveAs("spectra_nom_fluxtuned.eps");
  c->SaveAs("spectra_nom_fluxtuned.pdf");

  //flux tuned vs banff
  leg->Clear();
  leg->AddEntry(h_ereco_flux_osc,"Flux tuning applied","l");
  leg->AddEntry(h_ereco_tun_osc,"BANFF tuning applied","l");
  c->cd(1);
  h_ereco_flux_noosc->Draw();
  h_ereco_tun_noosc->Draw("same");
  leg->Draw();
  c->cd(2);
  h_ereco_flux_osc->Draw();
  h_ereco_tun_osc->Draw("same");
  leg->Draw();
  c->cd(3);
  TH1D * h_ereco_fluxbanff_noosc_ratio = new TH1D(*h_ereco_tun_noosc);
  h_ereco_fluxbanff_noosc_ratio->Divide(h_ereco_tun_noosc,h_ereco_flux_noosc);
  h_ereco_fluxbanff_noosc_ratio->SetLineColor(kBlack);
  h_ereco_fluxbanff_noosc_ratio->GetYaxis()->SetTitle("Banff tuned / Flux tuned");
  h_ereco_fluxbanff_noosc_ratio->GetYaxis()->SetRangeUser(0.7,1.1);
  h_ereco_fluxbanff_noosc_ratio->Draw();
  c->cd(4);
  TH1D * h_ereco_fluxbanff_osc_ratio = new TH1D(*h_ereco_tun_osc);
  h_ereco_fluxbanff_osc_ratio->Divide(h_ereco_tun_osc,h_ereco_flux_osc);
  h_ereco_fluxbanff_osc_ratio->SetLineColor(kBlack);
  h_ereco_fluxbanff_osc_ratio->GetYaxis()->SetTitle("Banff tuned / Flux tuned");
  h_ereco_fluxbanff_osc_ratio->GetYaxis()->SetRangeUser(0.7,1.1);
  h_ereco_fluxbanff_osc_ratio->Draw();
  c->SaveAs("spectra_tun_banfftuned.eps");
  c->SaveAs("spectra_tun_banfftuned.pdf");

}
