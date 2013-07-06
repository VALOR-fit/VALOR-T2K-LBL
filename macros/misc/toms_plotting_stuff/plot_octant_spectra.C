void plot_octant_spectra
(
 double erecomax = 7.0,
 double templemax = 10.0
)
{
  using namespace t2kana;

  TFile f_o1("2012results/SPECTRA/0.95/octant1/sk_ereco_spectrum-sinsq2theta[12_0.857,13_0.098,23_0.950],dcp_0.000,dmsq[12_0.00007500,23_0.002440]eV2.root");
  TFile f_o2("2012results/SPECTRA/0.95/octant2/sk_ereco_spectrum-sinsq2theta[12_0.857,13_0.098,23_0.950],dcp_0.000,dmsq[12_0.00007500,23_0.002440]eV2.root");

  TH1D * h_o1 = 0, *h_o2 = 0;
  f_o1->GetObject("sk_1Rmutotal",h_o1);
  f_o2->GetObject("sk_1Rmutotal",h_o2);

  //colors
  h_o1->SetLineColor(kBlack);
  h_o2->SetLineColor(kRed);

  //range
  h_o1->GetXaxis()->SetRangeUser(0,erecomax-conv::kOneMeV);
  h_o2->GetXaxis()->SetRangeUser(0,erecomax-conv::kOneMeV);

  //titles
  h_o1->GetXaxis()->SetTitle("E_{reco} / GeV");
  h_o2->GetXaxis()->SetTitle("E_{reco} / GeV");
  h_o1->GetYaxis()->SetTitle("Events per bin");
  h_o2->GetYaxis()->SetTitle("Events per bin");

  TCanvas * c = new TCanvas();
  TLegend * leg = new TLegend(1.1,0.9,0.3,0.7);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  //draw
  leg->Clear();
  leg->AddEntry(h_o1,"Octant 1","l");
  leg->AddEntry(h_o2,"Octant 2","l");
  h_o1->Draw();
  h_o2->Draw("same");
  leg->Draw();
  c->SaveAs("spectra_octants.eps");

}
