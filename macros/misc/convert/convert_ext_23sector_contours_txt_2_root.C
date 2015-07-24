// Read external 23-sector contours stored in this package 
// in text format and save them in ROOT format.
// CA.

{
  TString basedir = gSystem->Getenv("T2KOSCANA");
  
  TGraph * grMINOS90CL    = new TGraph( basedir + "/inputs/contours/minos2012_10p71E20POT+37p9ktyr_90CL.data" );
  TGraph * grMINOS90CL_sinsqtheta = new TGraph( basedir + "/inputs/contours/minos2013_10p71E20POT+3p4E20POT+37p9ktyr_90CL_normal_sinsqtheta.data" );
  TGraph * grSKZenith90CL = new TGraph( basedir + "/inputs/contours/sk2012_runs1to4_zenith3nu_90CL.data" );
  TGraph * grSKZenith90CL_sinsqtheta = new TGraph( basedir + "/inputs/contours/sk2012_runs1to4_zenith3nu_90CL_normal_sinsqtheta.data" );
  TGraph * grSKLE90CL     = new TGraph( basedir + "/inputs/contours/sk2012_runs1to4_LE2nu_90CL.data" );
  TGraph * grICE90CL      = new TGraph( basedir + "/inputs/contours/icecube_20102011_90CL.data" );

  TCanvas * c = new TCanvas("c","",20,20,800,500);
  TH2F * hframe = (TH2F*) c->DrawFrame(0.7, 1E-3, 1.00, 4E-3);
  hframe->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  hframe->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
  hframe->Draw();

  grMINOS90CL    -> SetLineColor(kBlack);
  grSKZenith90CL -> SetLineColor(kRed);
  grSKLE90CL     -> SetLineColor(kBlue);
  grICE90CL      -> SetLineColor(kGreen+2);

  grMINOS90CL    -> Draw("LSAME");
  grSKZenith90CL -> Draw("LSAME");
  grSKLE90CL     -> Draw("LSAME");
  grICE90CL      -> Draw("LSAME");

  TLegend * legend = new TLegend(0.2, 0.7, 0.8, 0.9);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(grMINOS90CL,    "MINOS 10.71E20 POT + 37.9 kton*yrs analysis (2 flavour), 90% CL", "L");
  legend->AddEntry(grSKZenith90CL, "SK1-4 zenith-angle analysis (3 flavour), 90% CL", "L");
  legend->AddEntry(grSKLE90CL,     "SK1-4 L/E analysis (2 flavour), 90% CL", "L");
  legend->AddEntry(grICE90CL,      "ICECUBE 2010-2011 (2 flavour), 90% CL", "L");
  legend->Draw();

  c->Update();

  TCanvas * c2 = new TCanvas("c2","",20,20,800,500);
  TH2F * hframe2 = (TH2F*) c->DrawFrame(0.25, 1E-3, 0.75, 4E-3);
  hframe2->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  hframe2->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
  hframe2->Draw();
  
  grMINOS90CL_sinsqtheta    -> SetLineColor(kBlack);
  grSKZenith90CL_sinsqtheta -> SetLineColor(kRed);

  grMINOS90CL_sinsqtheta -> Draw("LSAME");
  grSKZenith90CL_sinsqtheta -> Draw("LSAME");

  TFile fout("contours.root","recreate");
  grMINOS90CL    -> Write("minos2012_combined_90CL");
  grMINOS90CL_sinsqtheta -> Write("minos2013_combined_90CL_normal_sinsqtheta");
  grSKZenith90CL -> Write("sk2012_zenith3nu_90CL");
  grSKLE90CL     -> Write("sk2012_le2nu_90CL");
  grSKZenith90CL_sinsqtheta -> Write("sk2012_zenith3nu_90CL_normal_sinsqtheta");
  grICE90CL      -> Write("icecube20102011_90CL");
  fout.Close();
}
