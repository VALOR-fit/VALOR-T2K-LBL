void plot_tuning_steps_spectra
(
 double erecomax = 7.0,
 double templemax = 10.0
)
{
  using namespace t2kana;

  TFile f_tun_noosc ("alltune/tun_noosc.root");
  TFile f_tun_osc ("alltune/tun_osc.root");
  TFile f_flux_noosc ("fluxtune/flux_noosc.root");
  TFile f_flux_osc ("fluxtune/flux_osc.root");
  TFile f_nom_noosc ("notune/nom_noosc.root");
  TFile f_nom_osc ("notune/nom_osc.root");
  TFile f_nom_noosc_signue("notune_noosc_signuestillon/nom_noosc_withsignue.root");

  TH1D * h_ereco_tun_noosc = 0;
  f_tun_noosc.GetObject("ereco_spectrum",h_ereco_tun_noosc);
  TH1D * h_ereco_tun_osc = 0;
  f_tun_osc.GetObject("ereco_spectrum",h_ereco_tun_osc);
  TH1D * h_ereco_flux_noosc = 0;
  f_flux_noosc.GetObject("ereco_spectrum",h_ereco_flux_noosc);
  TH1D * h_ereco_flux_osc = 0;
  f_flux_osc.GetObject("ereco_spectrum",h_ereco_flux_osc);
  TH1D * h_ereco_nom_noosc = 0;
  f_nom_noosc.GetObject("ereco_spectrum",h_ereco_nom_noosc);
  TH1D * h_ereco_nom_osc = 0;
  f_nom_osc.GetObject("ereco_spectrum",h_ereco_nom_osc);

  //colors
  h_ereco_tun_noosc->SetLineColor(kBlack);
  h_ereco_tun_osc->SetLineColor(kBlack);
  h_ereco_flux_noosc->SetLineColor(kRed);
  h_ereco_flux_osc->SetLineColor(kRed);
  h_ereco_nom_noosc->SetLineColor(kBlue);
  h_ereco_nom_osc->SetLineColor(kBlue);

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

  TCanvas * c = new TCanvas();
  c->Divide(2,1);
  TLegend * leg = new TLegend(0.92,0.9,0.45,0.7);
  leg->SetFillColor(0);
  leg->SetLineColor(0);

  //nom vs flux tuned
  leg->Clear();
  c->cd(1);
  h_ereco_flux_noosc->Draw();
  h_ereco_nom_noosc->Draw("same");
  c->cd(2);
  h_ereco_flux_osc->Draw();
  h_ereco_nom_osc->Draw("same");
  leg->AddEntry(h_ereco_nom_osc,"Nominal","l");
  leg->AddEntry(h_ereco_flux_osc,"Flux tuning applied","l");
  leg->Draw();
  c->SaveAs("spectra_nom_fluxtuned.eps");

  //flux tuned vs banff
  leg->Clear();
  c->cd(1);
  h_ereco_flux_noosc->Draw();
  h_ereco_tun_noosc->Draw("same");
  c->cd(2);
  h_ereco_flux_osc->Draw();
  h_ereco_tun_osc->Draw("same");
  leg->AddEntry(h_ereco_flux_osc,"Flux tuning applied","l");
  leg->AddEntry(h_ereco_tun_osc,"BANFF tuning applied","l");
  leg->Draw();
  c->SaveAs("spectra_tun_banfftuned.eps");

  //Now get the templates
  const int ntempl = 28;
  TH2D * h_template[ntempl] = {0};
  for(int i = 0; i < ntempl; i++) {
    f_nom_noosc_signue.GetObject(Form("nsk_fin_0_%d",i),h_template[i]);
    //h_template[i]->SetTitle("");
    h_template[i]->GetXaxis()->SetRangeUser(0,templemax-conv::kOneMeV);
    h_template[i]->GetYaxis()->SetRangeUser(0,templemax-conv::kOneMeV);
  }
  TCanvas * ctempl = new TCanvas();
  ctempl->Divide(2, 3);
  for(int i = 0; i < ntempl; i++) {
    ctempl->cd((i%6) +1)->SetRightMargin(0.12);
    h_template[i]->Draw("colz");
    if(i%6 == 5) {
      ctempl->SaveAs(Form("template_%d.eps",i));
    }
  }
  //and the signue
  TCanvas * ctemplsig = new TCanvas();
  ctemplsig->Divide(2,2);
  for(int i = ntempl-4; i < ntempl; i++) {
    ctemplsig->cd((i%4) +1)->SetRightMargin(0.15);
    h_template[i]->Draw("colz");
  }
  ctemplsig->SaveAs("template_27.eps");

}
