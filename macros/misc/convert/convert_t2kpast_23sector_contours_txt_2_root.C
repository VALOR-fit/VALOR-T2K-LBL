// Read past-T2K anayses 23-sector contours stored in this package 
// in text format and save them in ROOT format.
// CA.

{
  TString basedir = gSystem->Getenv("T2KOSCANA");

  TGraph * gr_run1p2_90CL_statsyst = new TGraph( basedir + "/inputs/external/contours/t2k2010_run1+2_90CL_statsyst.data" );
  TGraph * gr_run1p2_90CL_stat     = new TGraph( basedir + "/inputs/external/contours/t2k2010_run1+2_90CL_stat.data" );
  TGraph * gr_run1p2_68CL_statsyst = new TGraph( basedir + "/inputs/external/contours/t2k2010_run1+2_68CL_statsyst.data" );
  TGraph * gr_run1p2_68CL_stat     = new TGraph( basedir + "/inputs/external/contours/t2k2010_run1+2_68CL_stat.data" );

  TCanvas * c = new TCanvas("c","",20,20,800,500);
  TH2F * hframe = (TH2F*) c->DrawFrame(0.82, 1.5E-3, 1.00, 4.5E-3);
  hframe->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  hframe->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
  hframe->Draw();

  gr_run1p2_90CL_statsyst -> SetLineColor(kBlue);
  gr_run1p2_90CL_stat     -> SetLineColor(kBlue);
  gr_run1p2_68CL_statsyst -> SetLineColor(kRed);
  gr_run1p2_68CL_stat     -> SetLineColor(kRed);

  gr_run1p2_90CL_statsyst -> SetLineStyle(kSolid);
  gr_run1p2_90CL_stat     -> SetLineStyle(kDashed);
  gr_run1p2_68CL_statsyst -> SetLineStyle(kSolid);
  gr_run1p2_68CL_stat     -> SetLineStyle(kDashed);

  gr_run1p2_90CL_statsyst -> Draw("LSAME");
  gr_run1p2_90CL_stat     -> Draw("LSAME");
  gr_run1p2_68CL_statsyst -> Draw("LSAME");
  gr_run1p2_68CL_stat     -> Draw("LSAME");

  TLegend * legend = new TLegend(0.2, 0.7, 0.8, 0.9);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(gr_run1p2_90CL_statsyst, "T2K Run1+2 (1.41E20 POT), 2010b analysis, FC 90% CL (stat+syst)", "L");
  legend->AddEntry(gr_run1p2_90CL_stat,     "T2K Run1+2 (1.41E20 POT), 2010b analysis, FC 90% CL (stat)",      "L");
  legend->AddEntry(gr_run1p2_68CL_statsyst, "T2K Run1+2 (1.41E20 POT), 2010b analysis, FC 68% CL (stat+syst)", "L");
  legend->AddEntry(gr_run1p2_68CL_stat,     "T2K Run1+2 (1.41E20 POT), 2010b analysis, FC 68% CL (stat)",      "L");
  legend->Draw();

  c->Update();

  TFile fout("t2kpast_contours.root","recreate");
  gr_run1p2_90CL_statsyst -> Write("t2k2010_run1p2_90CL_statsyst");
  gr_run1p2_90CL_stat     -> Write("t2k2010_run1p2_90CL_stat");
  gr_run1p2_68CL_statsyst -> Write("t2k2010_run1p2_68CL_statsyst");
  gr_run1p2_68CL_stat     -> Write("t2k2010_run1p2_68CL_stat");
  fout.Close();
}
