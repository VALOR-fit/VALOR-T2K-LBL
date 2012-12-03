// ------------------------------------------------------------------------------------------------
// Plot SK  1-ring mu E_{reco} spectrum error envelope.
// The input file is the one saved out by physin/test_superk_spectrum_prediction_err_envelope.C
// This is the macro used for Figs. 13,14 of T2K-TN-64.
//
// C.Andreopoulos
// ------------------------------------------------------------------------------------------------

void plot_1ringmu_ereco_spectrum_err_envelope(const char * filename)
{
  TFile * f = new TFile(filename,"read");

  TH1D * hErecoFracErrP = (TH1D*) f->Get("spectrum_frac_err_plus");	
  TH1D * hErecoFracErrM = (TH1D*) f->Get("spectrum_frac_err_minus");
  TH1D * hErecoMax      = (TH1D*) f->Get("spectrum_error_envelope_upper_boundary");	
  TH1D * hErecoMin      = (TH1D*) f->Get("spectrum_error_envelope_lower_boundary");
  TH1D * hErecoDef      = (TH1D*) f->Get("spectrum_nominal");

  hErecoDef->SetLineStyle(kDashed);
  hErecoMin->SetLineStyle(kSolid);
  hErecoMax->SetLineStyle(kSolid);
  hErecoMin->SetLineWidth(2);
  hErecoMax->SetLineWidth(2);
  hErecoFracErrP->SetLineWidth(2);
  hErecoFracErrM->SetLineWidth(2);

  double emin = 0;
  double emax = 4;
  double nmax = 2.4;
  double nfracerrmax = 1.3;

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->cd();
  TPad *pad1 = new TPad("pad1","pad1",.02,.30,.98,.98);
  pad1->SetFillColor(0);
  pad1->SetFrameFillColor(0);
  pad1->SetBorderSize(0);
  pad1->Draw();
  pad1->cd();
  TH2D * hf1 = (TH2D*) pad1->DrawFrame(emin,0,emax,nmax);
  hErecoMax  -> Draw("same");
  hErecoMin  -> Draw("same");
  hErecoDef  -> Draw("same");
  hf1->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hf1->GetYaxis()->SetTitle("Number of events");
  c->cd();
  TPad *pad2 = new TPad("pad1","pad1",.02,.02,.98,.28);
  pad2->SetFillColor(0);
  pad2->SetFrameFillColor(0);
  pad2->SetBorderSize(0);
  pad2->Draw();
  pad2->cd();
  TH2D * hf2 = (TH2D*) pad2->DrawFrame(emin,-1*nfracerrmax,emax,nfracerrmax);
  hErecoFracErrP -> Draw("same");
  hErecoFracErrM -> Draw("same");
  hf2->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hf2->GetYaxis()->SetTitle("Fractional error");
  c->Update();
}
