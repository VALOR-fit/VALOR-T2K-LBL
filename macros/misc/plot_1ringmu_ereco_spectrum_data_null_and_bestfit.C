//
// ------------------------------------------------------------------------------------------------
// Plot 1-ring mu E_{reco} spectrum (data, best-fit oscillation hypothesis and null hypothesis).
// This is the macro used for Figs. 17,18,19,20 of T2K-TN-64
//
// C.Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void plot_1ringmu_ereco_spectrum_data_null_and_bestfit
(
  /* fitter output ROOT file (data, best-fit spectrum total + contributions from various modes) */
  const char * fitfile,                 

  /* ROOT file with null hypothesis spectrum (eg from test_superk_spectrum_prediction_contributions.C) */
  const char * nullpredictionfile = 0,  

  const char * label   = 0,        /* legend label */
  bool   rebin         = false,    /* re-bin E_{reco} spectra (fewer variable size bins) */
  bool   show_numuccqe = false,    /* show numu CCQE component of the best-fit spectrum?  */
  double emax          = 6.0,      /* max E_{reco} in plot  */
  double nmax          = 20.0      /* max of N_{events} axis in plot  */
)
{
  const int kNErecoBins1Rmu = 15;
  const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = {
   0.000,  0.300,  0.400,  0.500,  0.600,  0.700,  0.800, 0.900,  
   1.000,  1.200,  1.500,  2.000,  3.000,  4.000,  6.000, 10.000 
  };


  TFile * ffit  = new TFile(fitfile, "read");
  TH1D * hmcall      = 0;
  TH1D * hmcnumuccqe = 0;
  TH1D * hdata       = 0;
  ffit->cd();
  if(rebin) {
     hmcall      = (TH1D*) bestfit_1Rmu       -> Rebin(kNErecoBins1Rmu, "hmcall",      kErecoBinEdges1Rmu);
     hmcnumuccqe = (TH1D*) bestfit_1Rmu_numuccqe  -> Rebin(kNErecoBins1Rmu, "hmcnumuccqe", kErecoBinEdges1Rmu);
     hdata       = (TH1D*) data_1Rmu              -> Rebin(kNErecoBins1Rmu, "hdata",       kErecoBinEdges1Rmu);
  } else {
     hmcall      = new TH1D(*bestfit_1Rmu);
     hmcnumuccqe = new TH1D(*bestfit_1Rmu_numuccqe);
     hdata       = new TH1D(*data_1Rmu);
  }

  cout << "Read " << hdata->GetEntries() << " events " << endl;


  TFile * fnull   = 0;
  TH1D *  hmcnull = 0;
  if(nullpredictionfile) {
     fnull = new TFile(nullpredictionfile, "read");
     fnull->cd();
     if(rebin) {
        hmcnull = (TH1D*) total -> Rebin(kNErecoBins1Rmu, "hmcnull",  kErecoBinEdges1Rmu);
     } else {
        hmcnull = new TH1D(*total);
     }
  }

  int marker_size  = 1;
  int marker_style = 8;
  if(rebin) {
    marker_size  = 1.0;
    marker_style = 8;
  }
  int numuccqe_color = kGray+3;

  if(hmcnull) {
     hmcnull  -> SetLineWidth(2);
     hmcnull  -> SetLineColor(kBlack);
     hmcnull  -> SetLineStyle(kDashed);
  }
  hmcall      -> SetLineWidth(2);
  hmcall      -> SetLineColor(kBlack);
  hmcall      -> SetLineStyle(kSolid);
  hmcnumuccqe -> SetFillStyle(3354);
  hmcnumuccqe -> SetLineColor(numuccqe_color);
  hmcnumuccqe -> SetFillColor(numuccqe_color);
  hdata       -> SetMarkerStyle(marker_style);
  hdata       -> SetMarkerSize(marker_size);
  hdata       -> SetLineWidth(1);

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TH1F * hframe = (TH1F*)c->DrawFrame(0,0,emax,nmax);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Number of events");
  if(hmcnull) {
     hmcnull->Draw("same");
  }
  hmcall->Draw("same");
  if(show_numuccqe) {
     hmcnumuccqe->Draw("same");
  }
  hdata->Draw("perrsame");

  TLegend * legend = new TLegend(0.6,0.6,0.9,0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  if(hmcnull) {
     legend->AddEntry(hmcnull,"Null hypothesis (total)","L");
  }
  legend->AddEntry(hmcall,"Best-fit oscillation hypothesis (total)","L");
  if(show_numuccqe) {
     legend->AddEntry(hmcnumuccqe,"Best-fit oscillation hypothesis (#nu_{#mu} CCQE)","LF");
  }
  if(!label) {
    label = "data";
  }
  legend->AddEntry(hdata,label,"LP");
  legend->Draw();
}

