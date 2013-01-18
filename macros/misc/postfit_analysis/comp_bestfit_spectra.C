//
//  Compare best-fit spectra from different analyses
//  Runs on the output of the fitter.
//

void comp_bestfit_spectra
(
	const char * file1,
	const char * label1 = "analysis a",
	const char * file2,
	const char * label2 = "analysis b",
	const char * datafile = 0,
	const char * nooscfile = 0,
	bool rebin = false,
        double emax = 6,
	double nmax = 20
)

{
  double emin = 0;

  const int kNErecoBins1Rmu = 16;
  const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = {
   0.000,  0.300,  0.400,  0.500,
   0.600,  0.700,  0.800,  0.900,  
   1.000,  1.100,  1.250,  1.500,  2.000, 3.000,  4.000, 6.000, 20.000
  };

  gStyle->SetOptStat(0);

  TFile * f1 = new TFile(file1,"read");
  f1->cd();
  TH1D * hfit1 = (TH1D *) f1->Get("sk_1Rmutotal");
  assert(hfit1);

  TFile * f2 = new TFile(file2,"read");
  f2->cd();
  TH1D * hfit2 = (TH1D *) f2->Get("sk_1Rmutotal");
  assert(hfit2);

  TH1D * hdata = 0;
  if(datafile){
    TFile * fd = new TFile(datafile,"read");
    TH1D * hdata = (TH1D *) fd->Get("obs_spectrum");
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
    histo_data = (TH1D*) hdata -> Rebin(kNErecoBins1Rmu, "histo_data", kErecoBinEdges1Rmu);
    histo_null = (TH1D*) hnull -> Rebin(kNErecoBins1Rmu, "histo_null", kErecoBinEdges1Rmu);
  } else {
    histo_fit1 = new TH1D (*hfit1);
    histo_fit2 = new TH1D (*hfit2);
    histo_data = new TH1D (*hdata);
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
  TPad * upp  = new TPad("upp", "",0.05,0.40,0.95,0.97);
  TPad * low1 = new TPad("low1","",0.05,0.20,0.95,0.38);
  TPad * low2 = new TPad("low2","",0.05,0.02,0.95,0.18);
  low1->SetBottomMargin(.2);
  low2->SetBottomMargin(.2);
  upp->Draw();
  low1->Draw();
  low2->Draw();
  upp->cd();

  TH1F * hframe = (TH1F*)upp->DrawFrame(emin,0,emax,nmax);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Events per bin");

  histo_null->SetLineWidth(1);
  histo_null->SetLineColor(kBlack);
  histo_fit1->SetLineStyle(kDotted);

  histo_fit1->SetLineWidth(2);
  histo_fit1->SetLineColor(kAzure-3);
  histo_fit1->SetLineStyle(kDashed);

  histo_fit2->SetLineWidth(2);
  histo_fit2->SetLineColor(kRed);
  histo_fit2->SetLineStyle(kSolid);

  histo_data->SetMarkerSize(1);
  histo_data->SetMarkerStyle(8);
  histo_data->SetLineStyle(kDashed);

  histo_null->Draw("same");
  histo_fit1->Draw("same");
  histo_fit2->Draw("same");
  histo_data->Draw("PE03X0 same");

  TLegend * legend = new TLegend(0.5, 0.6, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetHeader("1 #mu-like ring spectrum");
  legend->AddEntry(histo_data, "Run1+2+3 dataset", "P");
  legend->AddEntry(histo_null, "No oscillation expectation", "L");
  legend->AddEntry(histo_fit1, label1, "L");
  legend->AddEntry(histo_fit2, label2, "L");
  legend->Draw();

  low1->cd();
  TH1F * hframelow1 = (TH1F*)low1->DrawFrame(emin,-0.2,emax,0.2);
  hframelow1->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframelow1->GetYaxis()->SetTitle("\'2+N\'-\'2+8\' (Events per bin)");
  hframelow1->GetYaxis()->SetLabelSize(0.09);
  hframelow1->GetYaxis()->SetTitleSize(0.11);
  hframelow1->GetYaxis()->SetTitleOffset(0.62);
  hframelow1->GetYaxis()->SetNdivisions(505);
  hframelow1->GetXaxis()->SetLabelSize(0.11);
  hframelow1->GetXaxis()->SetTitleSize(0.11);
  low1->SetGridy();
  diff->Draw("same");

  low2->cd();
  TH1F * hframelow2 = (TH1F*)low2->DrawFrame(emin,0,emax,2);
  hframelow2->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframelow2->GetYaxis()->SetTitle("\'2+N\' / \'2+8\'");
  hframelow2->GetYaxis()->SetLabelSize(0.09);
  hframelow2->GetYaxis()->SetTitleSize(0.11);
  hframelow2->GetYaxis()->SetTitleOffset(0.62);
  hframelow2->GetYaxis()->SetNdivisions(505);
  hframelow2->GetXaxis()->SetLabelSize(0.11);
  hframelow2->GetXaxis()->SetTitleSize(0.11);
  low2->SetGridy();
  ratio->Draw("same");

//  c1->SaveAs("marg_dchisq_comp.eps");
}
