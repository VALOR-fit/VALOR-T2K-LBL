//plot mc statistical error, compared with current-pot stat error & systematic tweak error

void mcstaterr_spectrum()
{
  TCanvas * c = new TCanvas();
  double width = 0.5;//linewidth

  TFile f84  ("$T2KOSCANA/2013results/SPECTRA/mcstaterror/numu/Template_mcstaterr_84.root");
  TFile f84d ("$T2KOSCANA/2013results/SPECTRA/mcstaterror/numu/Template_mcstaterr_84_doublestat.root");
  TFile fsyst("$T2KOSCANA/2013results/SPECTRA/mcstaterror/numu/sk_ereco_spectrum-twk_sys_throw_sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00000000,23_0.000000]eV2.root");

  TH1D * h84 = 0, * h84d = 0, * h1 = 0, * hsyst = 0, * hnom = 0, * hstat = 0;

  //get nominal
  f84.GetObject("spe",hnom);
  int nbins = hnom->GetNbinsX();
  hnom->SetLineWidth(width);

  //all systematics tweaked
  fsyst.GetObject("all_pluserr",hsyst);
  hsyst->SetLineColor(kMagenta);
  hsyst->SetLineWidth(width);

  // current POT stat error
  fsyst.GetObject("nominal",hstat);
  hstat->SetLineColor(kBlack);
  hstat->SetLineWidth(width);
  for(int ix = 1; ix <= nbins; ix++) {
    double nev = hstat->GetBinContent(ix);
    if(nev < 1E-9) continue;
    hstat->SetBinContent(ix, 1 / TMath::Sqrt(nev));
    cout << hstat->GetBinContent(ix) << endl;
  }

  //84-bin 2D-error
  f84.GetObject("spe_with2derr",h84);
  h84->SetLineColor(kRed);
  h84->SetLineWidth(width);
  for(int ix = 1; ix <= nbins; ix++) {
    double env = h84->GetBinContent(ix);
    double nom = hnom->GetBinContent(ix);
    h84->SetBinContent(ix, ((nom > 0) ? (env-nom)/nom : 0));
  } 

  //1D-error
  f84.GetObject("spe_err",h1);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(width);
  for(int ix = 1; ix <= nbins; ix++) {
    double env = h1->GetBinContent(ix);
    double nom = hnom->GetBinContent(ix);
    h1->SetBinContent(ix, ((nom > 0) ? (env-nom)/nom : 0));
  }

  //84-bin 2D-error double MC stats
  f84d.GetObject("spe_with2derr",h84d);
  h84d->SetLineColor(kCyan);
  h84d->SetLineWidth(width);
  h84d->Scale(0.5);
  for(int ix = 1; ix <= nbins; ix++) {
    double env = h84d->GetBinContent(ix);
    double nom = hnom->GetBinContent(ix);
    h84d->SetBinContent(ix, ((nom > 0) ? (env-nom)/nom : 0));
  } 

  //draw
  //double ymax = TMath::Max(hsyst->GetMaximum() * 1.05, hstat->GetMaximum() * 1.05);
  double ymax = hsyst->GetMaximum() * 1.05;
  TH1F * hframe = (TH1F*) c->DrawFrame(0,0,30,ymax);
  hframe->GetYaxis()->SetTitle("Fractional error");
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->Draw();
  hsyst->Draw("same");
  //hstat->Draw("same");
  h84->Draw("same");
  h1->Draw("same");
  h84d->Draw("same");

  //legend
  TLegend * l = 0;
  //l = new TLegend(0.25,0.7,0.6,0.9);
  l = new TLegend(0.5,0.6,0.85,0.8);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.04);
  l->AddEntry(h84,"84 etrue bins * 32 modes","l");
  l->AddEntry(h1,      " 1 etrue bin  *  1 mode","l");
  l->AddEntry(h84d,"84 etrue bins * 32 modes","l");
  l->AddEntry((TObject*)0,"    + double MC stats","");
  l->AddEntry(hsyst,"All systematic tweaks","l");
  //l->AddEntry(hstat,"Run 1-4 stat error","l");
  l->Draw();

  c->Update();
  c->SaveAs("mcstaterr_spectrum.pdf");
  c->SaveAs("mcstaterr_spectrum.eps");
  hframe->GetXaxis()->SetRangeUser(0,2.9999);
  c->Update();
  c->SaveAs("mcstaterr_spectrum_zoom.pdf");
  c->SaveAs("mcstaterr_spectrum_zoom.eps");

  cout << hnom->Integral(-1,99999) << " total events in nominal" << endl;

  //ratio to systematic error
  TH1D * ratio = new TH1D(*h84);
  for(int ix = 1; ix <= nbins; ix++) {
    double syst = hsyst->GetBinContent(ix);
    double stat = h84->GetBinContent(ix);
    ratio->SetBinContent(ix,((syst > 0) ? stat/syst : 0));
    //ratio->SetBinContent(ix,((syst > 0) ? (syst-stat)/syst : 0));
    //cout << stat << "\t" << syst << "\t" << ((syst > 0) ? (syst-stat)/syst : 0) << "\t" << ((stat > 0) ? (syst-stat)/stat : 0) << endl;
  }
  //ratio->GetYaxis()->SetTitle("(all syst error - stat error)/(all syst error)");
  ratio->GetYaxis()->SetTitle("(stat error)/(all syst error)");
  ratio->Draw();
  c->SaveAs("mcstaterr_ratio.pdf");
  c->SaveAs("mcstaterr_ratio.eps");



  //ratio to stat error
  TH1D * ratio_mcstat_stat = new TH1D(*h84);
  TH1D * ratio_syst_stat = new TH1D(*h84);
  for(int ix = 1; ix <= nbins; ix++) {
    double syst   = hsyst->GetBinContent(ix);
    double mcstat = h84->GetBinContent(ix);
    double stat   = hstat->GetBinContent(ix);
    ratio_mcstat_stat->SetBinContent(ix, ((stat > 0) ? mcstat/stat : 0));
    ratio_syst_stat  ->SetBinContent(ix, ((stat > 0) ? syst  /stat : 0));
  }
  ratio_mcstat_stat->GetYaxis()->SetTitle("");
  ratio_syst_stat->GetYaxis()->SetTitle("");
  ratio_mcstat_stat->SetLineColor(kBlack);
  ratio_syst_stat->GetXaxis()->SetTitle("E_{reco} (GeV)");
  ratio_syst_stat->Draw("");
  ratio_mcstat_stat->Draw("same");
  TLegend * l2 = 0;
  l2 = new TLegend(0.4,0.65,0.75,0.8);
  l2->SetFillStyle(0);
  l2->SetBorderSize(0);
  l2->SetTextSize(0.04);
  l2->AddEntry(ratio_mcstat_stat,"(MC stat error) / (Run1-4 stat error)","l");
  l2->AddEntry(ratio_syst_stat,  "(all syst error) / (Run1-4 stat error)","l");
  l2->Draw();

  c->SaveAs("mcstaterr_ratio_stat.pdf");
  c->SaveAs("mcstaterr_ratio_stat.eps");

}
