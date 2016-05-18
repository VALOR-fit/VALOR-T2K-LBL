
void FormatCanvas(TCanvas *canvas)
{
  canvas->Divide(1);
  canvas->cd(1);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.15);
  //gPad->SetTopMargin(0.12);
  //gPad->SetRightMargin(0.16);
}

void FormatHistLong(TH1D *hist, int colour, int style)
{
  hist->SetTitle("");
  hist->SetLineColor(colour);
  hist->SetLineWidth(2);
  hist->SetLineStyle(style);
  hist->SetStats(0);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
}

void plot_flux_weights
(
 const bool doNu = true  // true --> Nu mode (run 1-6), false --> NuBar mode (run 5-6)
 )
{
  char *name_nu = "/Users/davidesgalaberna/Desktop/ETH_phd/T2K/ANA/VALOR/valor/softw/t2xk/branches/2016_v1.0.0_evolution/inputs/official/2015/flux/June2015/sk_tuned13av1.1_13anom_run1-6a_fine.root";
  char *name_nubar = "/Users/davidesgalaberna/Desktop/ETH_phd/T2K/ANA/VALOR/valor/softw/t2xk/branches/2016_v1.0.0_evolution/inputs/official/2015/flux/June2015/sk_tuned13av1.1_13anom_run5c-6e_antinumode_fine.root";

  
  TFile *weightfile;
  if(doNu) weightfile = TFile::Open(name_nu);
  else     weightfile = TFile::Open(name_nubar);

  weightfile->cd();
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  FormatCanvas(c1);
  enu_sk_tuned13a_numu_ratio->GetXaxis()->SetTitle("True E_{#nu} (GeV)");
  enu_sk_tuned13a_numu_ratio->GetYaxis()->SetTitle("13a tuned v1.1 / 13a nominal");
  enu_sk_tuned13a_numu_ratio->GetXaxis()->SetRange(0, 200);
  enu_sk_tuned13a_numu_ratio->GetYaxis()->SetRangeUser(0.2, 1.6);
  FormatHistLong(enu_sk_tuned13a_numu_ratio, 2, 1);
  enu_sk_tuned13a_numu_ratio->Draw("hist");

  if(doNu) c1->SaveAs("flux_ratio_13av1p1_to_13anom_numu_numode.eps");
  else     c1->SaveAs("flux_ratio_13av1p1_to_13anom_numu_nubarmode.eps");


  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  FormatCanvas(c2);
  enu_sk_tuned13a_numub_ratio->GetXaxis()->SetTitle("True E_{#nu} (GeV)");
  enu_sk_tuned13a_numub_ratio->GetYaxis()->SetTitle("13a tuned v1.1 / 13a nominal");
  enu_sk_tuned13a_numub_ratio->GetXaxis()->SetRange(0, 200);
  enu_sk_tuned13a_numub_ratio->GetYaxis()->SetRangeUser(0.2, 1.6);
  FormatHistLong(enu_sk_tuned13a_numub_ratio, 2, 1);
  enu_sk_tuned13a_numub_ratio->Draw("hist");

  if(doNu) c2->SaveAs("flux_ratio_13av1p1_to_13anom_numubar_numode.eps");
  else     c2->SaveAs("flux_ratio_13av1p1_to_13anom_numubar_nubarmode.eps");

  TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
  FormatCanvas(c3);
  enu_sk_tuned13a_nue_ratio->GetXaxis()->SetTitle("True E_{#nu} (GeV)");
  enu_sk_tuned13a_nue_ratio->GetYaxis()->SetTitle("13a tuned v1.1 / 13a nominal");
  enu_sk_tuned13a_nue_ratio->GetXaxis()->SetRange(0, 200);
  enu_sk_tuned13a_numu_ratio->GetYaxis()->SetRangeUser(0.2, 1.6);
  FormatHistLong(enu_sk_tuned13a_nue_ratio, 2, 1);
  enu_sk_tuned13a_nue_ratio->Draw("hist");

  if(doNu) c3->SaveAs("flux_ratio_13av1p1_to_13anom_nue_numode.eps");
  else     c3->SaveAs("flux_ratio_13av1p1_to_13anom_nue_nubarmode.eps");


  TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
  FormatCanvas(c4);
  enu_sk_tuned13a_nueb_ratio->GetXaxis()->SetTitle("True E_{#nu} (GeV)");
  enu_sk_tuned13a_nueb_ratio->GetYaxis()->SetTitle("13a tuned v1.1 / 13a nominal");
  enu_sk_tuned13a_nueb_ratio->GetXaxis()->SetRange(0, 200);
  enu_sk_tuned13a_nueb_ratio->GetYaxis()->SetRangeUser(0.2, 1.6);
  FormatHistLong(enu_sk_tuned13a_nueb_ratio, 2, 1);
  enu_sk_tuned13a_nueb_ratio->Draw("hist");

  if(doNu) c4->SaveAs("flux_ratio_13av1p1_to_13anom_nuebar_numode.eps");
  else     c4->SaveAs("flux_ratio_13av1p1_to_13anom_nuebar_nubarmode.eps");

}
