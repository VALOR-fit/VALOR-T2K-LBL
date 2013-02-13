void plot_dchisq_surfaces()
{
  TFile f_fine_28("$T2KOSCANA/2012results/DATA/fine_grain_28/fitresults_3.root");
  TFile f_fine_2all("$T2KOSCANA/2012results/DATA/fine_grain_2all/fitresults_3.root");

  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.1);
  c->SetLeftMargin(0.155);
		  
  //dchisq's & contours
  TH2D * dchisq_fine_28 = 0;
  f_fine_28.GetObject("hdchisq",dchisq_fine_28);
  dchisq_fine_28->SetLineColor(kRed);
  dchisq_fine_28->SetContour(255);
  dchisq_fine_28->GetYaxis()->SetTitleOffset(1.65);

  /*  TPaletteAxis *palette = (TPaletteAxis*)dchisq_fine_28->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  */

  dchisq_fine_28->GetZaxis()->SetRangeUser(0,20);
  dchisq_fine_28->Draw("colz");
  c->SaveAs("dchisq_fine_28.eps");
  cout << dchisq_fine_28->GetBinContent(dchisq_fine_28->GetMaximumBin()) << endl;

  
  //dchisq's & contours
  TH2D * dchisq_fine_2all = 0;
  f_fine_2all.GetObject("hdchisq",dchisq_fine_2all);
  dchisq_fine_2all->SetLineColor(kRed);
  dchisq_fine_2all->SetContour(255);
  dchisq_fine_2all->GetYaxis()->SetTitleOffset(1.65);

  /*  TPaletteAxis *palette = (TPaletteAxis*)dchisq_fine_2all->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  */

  dchisq_fine_2all->GetZaxis()->SetRangeUser(0,20);
  dchisq_fine_2all->Draw("colz");
  c->SaveAs("dchisq_fine_2all.eps");
  cout << dchisq_fine_2all->GetBinContent(dchisq_fine_2all->GetMaximumBin()) << endl;
  
}
