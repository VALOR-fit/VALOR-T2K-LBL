//------------------------------------------------------------------------------------------------
// Show best-fit spectrum & ratio. 
//
// C.Andreopoulos
//------------------------------------------------------------------------------------------------

void show_bestfit_spectrum_and_ratio
(
  const char * bestfitfile,       /* contains best-fit spectrum    */
  const char * label = "VaLOR",   /* label for analysis/fitter     */
  const char * datafile,          /* contains the data */
  const char * nooscfile,         /* contains the no-oscillation prediction */
  bool rebin = false,             /* rebin? */
  double emax = 6,                /* max energy plotted */
  double nmax = 20                /* max oy y axis in spectrum plots */
)

{
  gStyle->SetOptStat(0);

  const double emin = 0;

  // binning used if spectra are rebinned
  const int kNErecoBins1Rmu = 16;
  const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = 
  {
    0.000,  0.300,  0.400,  0.500,
    0.600,  0.700,  0.800,  0.900,  
    1.000,  1.100,  1.250,  1.500,  2.000, 3.000,  4.000, 6.000, 20.000
  };

  TFile * ffit = new TFile(bestfitfile,"read");
  ffit->cd();
  TH1D * hfit = (TH1D *) ffit->Get("sk_1Rmutotal");
  assert(hfit);
  TFile * fd = new TFile(datafile,"read");
  TH1D * hdata = (TH1D *) fd->Get("obs_spectrum");
  assert(hdata);
  TFile * fnull = new TFile(nooscfile,"read");
  TH1D * hnull = (TH1D *) fnull->Get("sk_1Rmutotal");
  assert(hnull);

  // create local copies, rebinned or not
  TH1D * histo_fit  = 0;
  TH1D * histo_data = 0;
  TH1D * histo_null = 0;
  if(rebin) {
    histo_fit  = (TH1D*) hfit  -> Rebin(kNErecoBins1Rmu, "histo_fit ", kErecoBinEdges1Rmu);
    histo_data = (TH1D*) hdata -> Rebin(kNErecoBins1Rmu, "histo_data", kErecoBinEdges1Rmu);
    histo_null = (TH1D*) hnull -> Rebin(kNErecoBins1Rmu, "histo_null", kErecoBinEdges1Rmu);
  } else {
    histo_fit  = new TH1D (*hfit);
    histo_data = new TH1D (*hdata);
    histo_null = new TH1D (*hnull);
  }

  TH1D * ratio_fit  = new TH1D(*histo_fit);
  TH1D * ratio_data = new TH1D(*histo_fit);
  for(int i=1; i <= histo_fit->GetXaxis()->GetNbins(); i++) {
	double nf   = histo_fit  -> GetBinContent(i);
	double nd   = histo_data -> GetBinContent(i);
	double n0   = histo_null -> GetBinContent(i);

        double rf =  (n0>0) ? nf/n0 : 0;
        double rd =  (n0>0) ? nd/n0 : 0;

	double derr = TMath::Sqrt(nd);
	double rerr = (nd>0) ? rd*(derr/nd) : 0;

	ratio_fit ->SetBinContent(i,rf);
	ratio_data->SetBinContent(i,rd);
	ratio_data->SetBinError(i,rerr);

	cout << "E = " << histo_fit->GetBinCenter(i) << ", ratio(data) = " << rd << endl;
  }


  histo_null->SetLineWidth(1);
  histo_null->SetLineColor(kBlack);
  histo_fit ->SetLineStyle(kDotted);

  histo_fit->SetLineWidth(2);
  histo_fit->SetLineColor(kRed);
  histo_fit->SetLineStyle(kSolid);

  histo_data->SetMarkerSize(1);
  histo_data->SetMarkerStyle(8);
  histo_data->SetLineStyle(kDashed);

  ratio_fit->SetLineWidth(2);
  ratio_fit->SetLineColor(kRed);
  ratio_fit->SetLineStyle(kSolid);

  ratio_data->SetMarkerSize(1);
  ratio_data->SetMarkerStyle(8);
  ratio_data->SetLineStyle(kDashed);

  TCanvas * c1 = new TCanvas;//("c1","",10,10,600,700); // a squarer default canvas is used in the TN
  TPad * upp = new TPad("upp","",0.05,0.30,0.95,0.97);
  TPad * low = new TPad("low","",0.05,0.03,0.95,0.28);
  low->SetBottomMargin(.2);
  upp->Draw();
  low->Draw();
  upp->cd();

  TH1F * hframe = (TH1F*)upp->DrawFrame(emin,0,emax,nmax);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Events per bin");

  histo_null->Draw("same");
  histo_fit ->Draw("same");
  histo_data->Draw("PE03X0 same");

  TLegend * legend = new TLegend(0.5, 0.6, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetHeader("1 #mu-like ring spectrum");
  legend->AddEntry(histo_data, "Run1+2+3 dataset", "P");
  legend->AddEntry(histo_null, "No oscillation expectation", "L");
  legend->AddEntry(histo_fit, label, "L");
  legend->Draw();

  low->cd();
  TH1F * hframelow = (TH1F*)low->DrawFrame(emin,-1,emax,2);
  hframelow->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframelow->GetYaxis()->SetTitle("Ratio");
  hframelow->GetYaxis()->SetLabelSize(0.09);
  hframelow->GetYaxis()->SetTitleSize(0.11);
  hframelow->GetYaxis()->SetTitleOffset(0.5);
  hframelow->GetYaxis()->SetNdivisions(505);
  hframelow->GetXaxis()->SetLabelSize(0.11);
  hframelow->GetXaxis()->SetTitleSize(0.11);
  low->SetGridy();
  ratio_fit->Draw("same");
  ratio_data->Draw("PE03X0 same");

  c1->SaveAs("spectrum_and_ratio_coarse.eps");
}
