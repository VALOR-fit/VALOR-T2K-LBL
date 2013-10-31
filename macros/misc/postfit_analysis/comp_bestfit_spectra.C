//------------------------------------------------------------------------------------------------
// Compare best-fit spectra from 2 different analyses.
// Superimpose, plot difference and ratio.
//
// C.Andreopoulos
//------------------------------------------------------------------------------------------------

TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}

void comp_bestfit_spectra
(
  const char * file1,                  /* 1st best-fit spectrum */
  const char * label1 = "analysis a",  /* label for the legend */
  const char * file2,                  /* 2nd best-fit spectrum */
  const char * label2 = "analysis b",  /* label for the legend */
  const char * datafile  = 0,          /* if set, superimpose data */
  const char * nooscfile = 0,          /* if set, superimpose no oscillation prediction	 */
  bool rebin = false,
  double emax = 6,
  double nmax = 20,
  double maxdiff  = 0.05,              /* set y-axis of middle canvas to [0-maxdiff, 0+maxdiff] */
  double maxratio = 0.05,              /* set y-axis of bottom canvas to [1-maxratio,1+maxratio] */
  const char * tag = "",               /* tag for the output file */
  const char * label1short = "A",      /* label for the y-axis of ratio & subtraction canvases */
  const char * label2short = "B",      /* label for the y-axis of ratio & subtraction canvases */
  const char * labeldata   = "Run 1-4",/* label the data in the legend */
  bool different_pot = false,          /* if set, will rescale hist2 to hist1 pot */
  double pot_scale = 6.5700/6.3933
)

{
  const double emin = 0;

  // if rebin, use this binning
  const int kNErecoBins1Rmu = 16;
  const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = {
   0.000,  0.300,  0.400,  0.500,
   0.600,  0.700,  0.800,  0.900,  
   1.000,  1.100,  1.250,  1.500,  2.000, 3.000,  4.000, 6.000, 20.000
  };

  gStyle->SetOptStat(0);

  TFile * f1 = new TFile(file1,"read");
  f1->cd();
  TH1D * hfit1 = (TH1D *) f1->Get("pdf_1Rmu");
  assert(hfit1);

  TFile * f2 = new TFile(file2,"read");
  f2->cd();
  TH1D * hfit2 = (TH1D *) f2->Get("pdf_1Rmu");
  assert(hfit2);
  if(different_pot)
    hfit2->Scale(pot_scale);

  TH1D * hdata = 0;
  if(datafile){
    TFile * fd = new TFile(datafile,"read");
    TH1D * hdata = (TH1D *) fd->Get("obs_1Rmu");
    assert(hdata);
  }

  TH1D * hnull = 0;
  if(nooscfile) {
    TFile * fnull = new TFile(nooscfile,"read");
    hnull = (TH1D *) fnull->Get("sk_1Rmutotal");
    assert(hnull);
  }

  // create local copies, rebinned or not
  TH1D * histo_fit1 = 0;
  TH1D * histo_fit2 = 0;
  TH1D * histo_data = 0;
  TH1D * histo_null = 0;
  if(rebin) {
    histo_fit1 = (TH1D*) hfit1 -> Rebin(kNErecoBins1Rmu, "histo_fit1", kErecoBinEdges1Rmu);
    histo_fit2 = (TH1D*) hfit2 -> Rebin(kNErecoBins1Rmu, "histo_fit2", kErecoBinEdges1Rmu);
    if(datafile)
      histo_data = (TH1D*) hdata -> Rebin(kNErecoBins1Rmu, "histo_data", kErecoBinEdges1Rmu);
    if(nooscfile) 
      histo_null = (TH1D*) hnull -> Rebin(kNErecoBins1Rmu, "histo_null", kErecoBinEdges1Rmu);
  } else {
    histo_fit1 = new TH1D (*hfit1);
    histo_fit2 = new TH1D (*hfit2);
    if(datafile)
      histo_data = new TH1D (*hdata);
    if(nooscfile) 
      histo_null = new TH1D (*hnull);
  }

  TH1D * diff  = new TH1D(*histo_fit1);
  TH1D * ratio = new TH1D(*histo_fit1);
  for(int i=1; i < histo_fit1->GetXaxis()->GetNbins(); i++) {
	double n1 = histo_fit1->GetBinContent(i);
	double n2 = histo_fit2->GetBinContent(i);
	double vratio = (n1!=0 && n2!=0) ? n1/n2 : 0;
	double vdiff  = n1-n2;
	diff->SetBinContent(i,vdiff);
	ratio->SetBinContent(i,vratio);
  }

  TCanvas * c1 = new TCanvas("c1","",10,10,600,700);
  TPad * upp  = new TPad("upp", "",0.02,0.40,1.00,0.97);
  TPad * low1 = new TPad("low1","",0.02,0.20,1.00,0.40);
  TPad * low2 = new TPad("low2","",0.02,0.00,1.00,0.20);
  low1->SetBottomMargin(.25);
  low2->SetBottomMargin(.25);
  upp->Draw();
  low1->Draw();
  low2->Draw();
  upp->cd();

  TH1F * hframe = (TH1F*)upp->DrawFrame(emin,0,emax,nmax);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Events per 0.050 GeV");

  if(nooscfile) {
    histo_null->SetLineWidth(1);
    histo_null->SetLineColor(kBlack);
    histo_fit1->SetLineStyle(kDotted);
  }

  histo_fit1->SetLineWidth(2);
  histo_fit1->SetLineColor(kAzure-3);
  histo_fit1->SetLineStyle(kDashed);

  histo_fit2->SetLineWidth(2);
  histo_fit2->SetLineColor(kRed);
  histo_fit2->SetLineStyle(kSolid);

  if(datafile) {
    histo_data->SetMarkerSize(0.5);
    histo_data->SetMarkerStyle(8);
    histo_data->SetLineStyle(kDashed);
  }

  if(nooscfile) {
    TH1D * histoper50_null = spectrum_to_per_50mev(histo_null);
    histoper50_null->Draw("same");
  }
  TH1D * histoper50_fit2 = spectrum_to_per_50mev(histo_fit2);
  TH1D * histoper50_fit1 = spectrum_to_per_50mev(histo_fit1);
  histoper50_fit2->Draw("same");
  histoper50_fit1->Draw("same");
  if(datafile) {
    TH1D * histoper50_data = spectrum_to_per_50mev(histo_data);
    histoper50_data->Draw("PE03X0 same");
  }

  TLegend * legend = new TLegend(0.4, 0.6, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.05);
  legend->SetHeader("1 #mu-like ring spectrum");
  if(datafile)
    legend->AddEntry(histo_data, Form("%s dataset",labeldata), "P");
  if(nooscfile)
    legend->AddEntry(histo_null, "No oscillation expectation", "L");
  legend->AddEntry(histo_fit1, label1, "L");
  legend->AddEntry(histo_fit2, label2, "L");
  legend->Draw();

  low1->cd();
  TH1F * hframelow1 = (TH1F*)low1->DrawFrame(emin,-maxdiff,emax,maxdiff);
  hframelow1->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframelow1->GetYaxis()->SetTitle(Form("#splitline{%s-%s}{(Events per bin)}",label1short,label2short));
  hframelow1->GetYaxis()->SetLabelSize(0.09);
  hframelow1->GetYaxis()->SetTitleSize(0.11);
  hframelow1->GetYaxis()->SetTitleOffset(0.62);
  hframelow1->GetYaxis()->SetNdivisions(505);
  hframelow1->GetXaxis()->SetLabelSize(0.11);
  hframelow1->GetXaxis()->SetTitleSize(0.11);
  low1->SetGridy();
  diff->Draw("same");

  low2->cd();
  TH1F * hframelow2 = (TH1F*)low2->DrawFrame(emin,1-maxratio,emax,1+maxratio);
  hframelow2->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframelow2->GetYaxis()->SetTitle(Form("%s / %s", label1short, label2short));
  hframelow2->GetYaxis()->SetLabelSize(0.09);
  hframelow2->GetYaxis()->SetTitleSize(0.11);
  hframelow2->GetYaxis()->SetTitleOffset(0.62);
  hframelow2->GetYaxis()->SetNdivisions(505);
  hframelow2->GetXaxis()->SetLabelSize(0.11);
  hframelow2->GetXaxis()->SetTitleSize(0.11);
  low2->SetGridy();
  ratio->Draw("same");

  c1->SaveAs(Form("spectrum_comparison_%s.eps",tag));
  c1->SaveAs(Form("spectrum_comparison_%s.pdf",tag));

  TCanvas * c2 = new TCanvas();
  hframe->Draw();
  if(nooscfile) 
    histo_null->Draw("same");
  histo_fit2->Draw("same");
  histo_fit1->Draw("same");
  if(datafile)
    histo_data->Draw("PE03X0 same");
  legend->Draw();

  c2->SaveAs(Form("spectrum_comparison_%s_only.pdf",tag));
}
