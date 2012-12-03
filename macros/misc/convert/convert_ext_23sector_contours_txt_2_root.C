// Read external 23-sector contours stored in this package 
// in text format and save them in ROOT format.
// CA.

{
  TString basedir = gSystem->Getenv("T2KOSCANA");
  
  TGraph * grMINOS90CL    = new TGraph( basedir + "/inputs/external/contours/minos2012_10p71E20POT+37p9ktyr_90CL.data" );
  TGraph * grSKZenith90CL = new TGraph( basedir + "/inputs/external/contours/sk2012_runs1to4_zenith3nu_90CL.data" );
  TGraph * grSKLE90CL     = new TGraph( basedir + "/inputs/external/contours/sk2012_runs1to4_LE2nu_90CL.data" );

  TCanvas * c = new TCanvas("c","",20,20,800,500);
  TH2F * hframe = (TH2F*) c->DrawFrame(0.82, 1.5E-3, 1.00, 4.5E-3);
  hframe->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  hframe->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
  hframe->Draw();

  grMINOS90CL    -> SetLineColor(kBlack);
  grSKZenith90CL -> SetLineColor(kRed);
  grSKLE90CL     -> SetLineColor(kBlue);

  grMINOS90CL    -> Draw("LSAME");
  grSKZenith90CL -> Draw("LSAME");
  grSKLE90CL     -> Draw("LSAME");

  TLegend * legend = new TLegend(0.2, 0.7, 0.8, 0.9);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(grMINOS90CL,    "MINOS 10.71E20 POT + 37.9 kton*yrs analysis (2 flavour), 90% CL", "L");
  legend->AddEntry(grSKZenith90CL, "SK1-4 zenith-angle analysis (3 flavour), 90% CL", "L");
  legend->AddEntry(grSKLE90CL,     "SK1-4 L/E analysis (2 flavour), 90% CL", "L");
  legend->Draw();

  c->Update();

  TFile fout("contours.root","recreate");
  grMINOS90CL    -> Write("minos2012_combined_90CL");
  grSKZenith90CL -> Write("sk2012_zenith3nu_90CL");
  grSKLE90CL     -> Write("sk2012_le2nu_90CL");
  fout.Close();
}
