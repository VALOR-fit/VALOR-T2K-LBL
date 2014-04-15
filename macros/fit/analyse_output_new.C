//
// Macro to make a series of plots, using the result of the macro run_fit_xml_seeded.C
//
// 1. Bestfit spectrum, ratio, (comparison with other spectrum - TO IMPLEMENT), ...
// 2. Pull plot, best fit sys point plot
// 3. chisq
// 4. dchisq
// 5. contour (constant dchisq method)
//
// This macro has been used to make plots for TN155
//

void set_line_style(TGraph *contour, Int_t color, Double_t width, Int_t style)
{
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
}
void set_line_style(TH1 *hist, Int_t color, Double_t width, Int_t style)
{
  hist->SetLineWidth(width);
  hist->SetLineColor(color);
  hist->SetLineStyle(style);
}
void set_line_style(TLine *line, Int_t color, Double_t width, Int_t style)
{
  line->SetLineWidth(width);
  line->SetLineColor(color);
  line->SetLineStyle(style);
}
void set_fill_style(TH1 *hist, Int_t color, Int_t style)
{
  hist->SetFillColor(color);
  hist->SetFillStyle(style);
}
void set_marker_style(TH1 *hist, Int_t color, Double_t size, Int_t style)
{
  hist->SetMarkerSize(size);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(style);
}
void set_marker_style(TMarker *marker, Int_t color, Double_t size, Int_t style)
{
  marker->SetMarkerSize(size);
  marker->SetMarkerColor(color);
  marker->SetMarkerStyle(style);
}
void save_canvas(TCanvas * c, const char * name, const char * tag) {
  c->SaveAs(Form("%s_%s.pdf",tag,name));
  c->SaveAs(Form("%s_%s.eps",tag,name));
  c->SaveAs(Form("%s_%s.C",  tag,name));
}
void set_legend_style(TLegend *leg, Double_t size = 0.05, Int_t style = 0) {
  leg->SetFillColor(0);
  leg->SetFillStyle(style);
  leg->SetBorderSize(0);
  leg->SetTextSize(size);
}

TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}



void analyse_output_new
(
 const char * infile,
 const char * tag,
 const char * config_file          = "analysis_configuration.xml",   /* for generating no-osc spectrum */
 double       spectrum_max_xaxis   = 7,                              /* max of spectrum axes */
 double       mixing_min           = 0.93,                           /* chisq & contour boundaries */
 double       mixing_max           = 1.0,
 double       dm_min               = 0.002,
 double       dm_max               = 0.003,
 const char * mixparamchoice       = "sinsq2theta",                   /* sinsq2theta or sinsqtheta, for the contour x-axis */
 const char * other_spectrum_file  = "",
 const char * other_spectrum_name  = "",
 bool         smooth               = false,
 int          dmsqmode             = 0 //0 = |dm32sq|, 1 = |dm31sq|
)
{
  using namespace std;
  using namespace t2kana;

  //global options
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  TFile fin(infile,"READ");
  if(fin.IsZombie()) {
    cerr << "Cannot read input file: " << infile;
    exit(1);
  }

  //a regular canvas
  TCanvas * c = new TCanvas();
  TLegend * l;

  //////////////////////////////////////////////////////////////////////
  //SPECTRUM

  // 1 spectrum only
  // get no osc spectrum
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
  PhysHypothesis ph(conv::kHypOsc3F); //initialisation is to no oscillations
  skexp->Generate(ph);
  const int      n_ereco_bins    = util::NumErecoBins (conv::kSKsFCFV1RingMuLike);
  const double * ereco_bin_edges = util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike);
  TH1D * spectrum_noosc = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  for(int ibin = 1; ibin <= n_ereco_bins; ibin++)
    spectrum_noosc->SetBinContent(ibin, skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingMuLike)[ibin-1]);   
  // get bestfit spectrum & observation
  TH1D * spectrum_obs = 0, * spectrum_bf = 0;
  fin.GetObject("pdf_1Rmu",spectrum_bf);
  fin.GetObject("obs_1Rmu",spectrum_obs);
  //axis range
  spectrum_obs->  GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_bf->   GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_noosc->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  //axis title
  spectrum_obs->GetXaxis()->  SetTitle("E_{reco} / GeV");
  spectrum_obs->GetYaxis()->  SetTitle("Events per bin");
  spectrum_noosc->GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_noosc->GetYaxis()->SetTitle("Events per bin");
  spectrum_bf->GetXaxis()->   SetTitle("E_{reco} / GeV");
  spectrum_bf->GetYaxis()->   SetTitle("Events per bin");
  //styles
  set_line_style  (spectrum_noosc, kBlack, 2,   1);
  set_line_style  (spectrum_bf,    kRed,   2,   1);
  set_line_style  (spectrum_obs,   kBlack, 1.0, 1);
  set_marker_style(spectrum_obs,   kBlack, 1.0, 8);
  spectrum_obs->SetLineStyle(1);
  //legend
  delete l;
  l = new TLegend(0.43,0.7,0.9,0.9);
  set_legend_style(l);
  l->AddEntry(spectrum_obs,  "Run 1-4 data","P");
  l->AddEntry(spectrum_noosc,"No oscillation expectation","l");
  l->AddEntry(spectrum_bf,   "Best fit","l");
  //draw
  c->cd(1);
  spectrum_noosc->Draw();
  spectrum_obs->  Draw("PE0same");
  spectrum_bf->   Draw("same");
  l->Draw();
  save_canvas(c,"spectrum",tag);
  //make per 50 mev
  TH1D * spectrumper50_noosc = spectrum_to_per_50mev(spectrum_noosc);
  TH1D * spectrumper50_obs   = spectrum_to_per_50mev(spectrum_obs);
  TH1D * spectrumper50_bf    = spectrum_to_per_50mev(spectrum_bf);
  //draw
  c->cd(1);
  spectrumper50_noosc->Draw();
  spectrumper50_obs->  Draw("PE0same");
  spectrumper50_bf->   Draw("same");
  l->Draw();
  save_canvas(c,"spectrumper50",tag);

  // 2 ratio only
  // get ratio
  spectrum_obs->Sumw2();
  TH1D * spectrum_bf_ratio  = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  TH1D * spectrum_obs_ratio = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  spectrum_bf_ratio-> Divide(spectrum_bf, spectrum_noosc);
  spectrum_obs_ratio->Divide(spectrum_obs, spectrum_noosc);
  //make sure bf has no errors
  for(int ir = 0; ir <= n_ereco_bins+1; ir++) {
    spectrum_bf_ratio->SetBinError(ir,0);
  }
  TLine * spectrum_line = new TLine(0, 1, spectrum_bf_ratio->GetXaxis()->GetBinUpEdge(spectrum_bf_ratio->GetXaxis()->FindBin(spectrum_max_xaxis-conv::kOneMeV)),1);
  //axis range
  spectrum_obs_ratio->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_bf_ratio-> GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  //axis title
  spectrum_obs_ratio->GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_obs_ratio->GetYaxis()->SetTitle("Ratio to no oscillations");
  spectrum_bf_ratio-> GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_bf_ratio-> GetYaxis()->SetTitle("Ratio to no oscillations");
  //styles
  set_line_style  (spectrum_bf_ratio,    kRed,   2,   1);
  set_line_style  (spectrum_obs_ratio,   kBlack, 1.0, 1);
  set_marker_style(spectrum_obs_ratio,   kBlack, 0.5, 8);
  set_line_style  (spectrum_line,        kBlack, 2,   2);
  //legend
  delete l;
  l = new TLegend(0.65,0.8,0.9,0.9);
  set_legend_style(l);
  l->AddEntry(spectrum_obs_ratio,  "Run 1-4 data","P");
  l->AddEntry(spectrum_bf_ratio,   "Best fit","l");
  //draw
  c->cd(1);
  spectrum_obs_ratio->Draw("PE0");
  spectrum_bf_ratio-> Draw("same");
  spectrum_line->     Draw();
  l->Draw();
  save_canvas(c,"spectrum_ratio",tag);

  // 3 rebinned spectrum only
  const int nbins = 20;
  //double binedges[nbins+1] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,2.0,2.5,3.0,4.0,6.0,10.0,30.0};
  double binedges[nbins+1] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,2.0,2.5,3.0,4.0,7.0,10.0,30.0};
  TH1D * spectrum_noosc_rebin = spectrum_noosc->Rebin(nbins,"spectrum_noosc_rebin",binedges);
  TH1D * spectrum_obs_rebin   = spectrum_obs->  Rebin(nbins,"spectrum_obs_rebin",  binedges);
  TH1D * spectrum_bf_rebin    = spectrum_bf->   Rebin(nbins,"spectrum_bf_rebin",   binedges);
  //axis range
  spectrum_obs_rebin->  GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_bf_rebin->   GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_noosc_rebin->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  //axis title
  spectrum_obs_rebin->  GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_obs_rebin->  GetYaxis()->SetTitle("Events per bin");
  spectrum_noosc_rebin->GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_noosc_rebin->GetYaxis()->SetTitle("Events per bin");
  spectrum_bf_rebin->   GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_bf_rebin->   GetYaxis()->SetTitle("Events per bin");
  //styles
  set_line_style  (spectrum_noosc_rebin, kBlack, 2,   1);
  set_line_style  (spectrum_bf_rebin,    kRed,   2,   1);
  set_line_style  (spectrum_obs_rebin,   kBlack, 1.0, 1);
  set_marker_style(spectrum_obs_rebin,   kBlack, 1.0, 8);
  //legend
  l = new TLegend(0.43,0.7,0.9,0.9);
  set_legend_style(l);
  l->AddEntry(spectrum_obs_rebin,  "Run 1-4 data","P");
  l->AddEntry(spectrum_noosc_rebin,"No oscillation expectation","l");
  l->AddEntry(spectrum_bf_rebin,   "Best fit","l");
  //draw
  c->cd(1);
  spectrum_noosc_rebin->Draw();
  spectrum_obs_rebin->  Draw("PE0same");
  spectrum_bf_rebin->   Draw("same");
  l->Draw();
  save_canvas(c,"spectrum_rebin",tag);
  //make per 50 mev
  TH1D * spectrumper50_noosc_rebin = spectrum_to_per_50mev(spectrum_noosc_rebin);
  TH1D * spectrumper50_obs_rebin = spectrum_to_per_50mev(spectrum_obs_rebin);
  TH1D * spectrumper50_bf_rebin = spectrum_to_per_50mev(spectrum_bf_rebin);
  //draw
  c->cd(1);
  spectrumper50_noosc_rebin->Draw();
  spectrumper50_obs_rebin->  Draw("PE0same");
  spectrumper50_bf_rebin->   Draw("same");
  l->Draw();
  save_canvas(c,"spectrumper50_rebin",tag);

  // 4 rebinned ratio only
  // get ratio
  spectrum_obs->Sumw2();
  TH1D * spectrum_bf_ratio_rebin  = new TH1D("","", nbins, binedges);
  TH1D * spectrum_obs_ratio_rebin = new TH1D("","", nbins, binedges);
  spectrum_bf_ratio_rebin-> Divide(spectrum_bf_rebin, spectrum_noosc_rebin);
  spectrum_obs_ratio_rebin->Divide(spectrum_obs_rebin, spectrum_noosc_rebin);
  //make sure bf has no errors
  for(int ir = 0; ir <= nbins+1; ir++) {
    spectrum_bf_ratio_rebin->SetBinError(ir,0);
  }
  TLine * spectrum_line_rebin = new TLine(0, 1, spectrum_bf_ratio_rebin->GetXaxis()->GetBinUpEdge(spectrum_bf_ratio_rebin->GetXaxis()->FindBin(spectrum_max_xaxis-conv::kOneMeV)),1);
  //axis range
  spectrum_obs_ratio_rebin->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_bf_ratio_rebin-> GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  //axis title
  spectrum_obs_ratio_rebin->GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_obs_ratio_rebin->GetYaxis()->SetTitle("Ratio to no oscillations");
  spectrum_bf_ratio_rebin-> GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_bf_ratio_rebin-> GetYaxis()->SetTitle("Ratio to no oscillations");
  //styles
  set_line_style  (spectrum_bf_ratio_rebin,    kRed,   2,   1);
  set_line_style  (spectrum_obs_ratio_rebin,   kBlack, 1.0, 1);
  set_marker_style(spectrum_obs_ratio_rebin,   kBlack, 0.5, 8);
  set_line_style  (spectrum_line_rebin,        kBlack, 2,   2);
  //legend
  delete l;
  l = new TLegend(0.65,0.8,0.9,0.9);
  set_legend_style(l);
  l->AddEntry(spectrum_obs_ratio_rebin,  "Run 1-4 data","P");
  l->AddEntry(spectrum_bf_ratio_rebin,   "Best fit","l");
  //draw
  c->cd(1);
  spectrum_obs_ratio_rebin->Draw("PE0");
  spectrum_bf_ratio_rebin-> Draw("same");
  spectrum_line_rebin->     Draw();
  l->Draw();
  save_canvas(c,"spectrum_ratio_rebin",tag);

  // 5 spectrum broken into mode categories
  //get the indiviudal spectra
  TH1D * spectrum_modes[conv::knm] = {0};
  for(int im = 0; im < conv::knm; im++) {
    spectrum_modes[im] = 0;
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    fin.GetObject(Form("pdf_1Rmu_%s",util::AsStringCompact(m).c_str()),spectrum_modes[im]);
  }//im
  //fill the categories
  TH1D * spectrum_modes_numuccqe   = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  TH1D * spectrum_modes_numuccnqe  = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  TH1D * spectrum_modes_numubarcc  = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  TH1D * spectrum_modes_nuecc      = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  TH1D * spectrum_modes_nc         = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  double total_integral = 0;
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    if(!spectrum_modes[im]) {
      cerr << "Spectrum for mode " << util::AsStringCompact(m) << " not found" << endl;
      continue;
    }
    for(int ir = 1; ir <= n_ereco_bins; ir++) {
      if(m == conv::kSKmNuMuCCQE || m == conv::kSKmOscNuMuCCQE)
	spectrum_modes_numuccqe->SetBinContent (ir, spectrum_modes_numuccqe->GetBinContent(ir)  + spectrum_modes[im]->GetBinContent(ir));
      else if((util::IsNuMu(m) || util::IsOscNuMu(m)) && util::IsCC(m))
	spectrum_modes_numuccnqe->SetBinContent(ir, spectrum_modes_numuccnqe->GetBinContent(ir) + spectrum_modes[im]->GetBinContent(ir));
      else if(util::IsNC(m))
	spectrum_modes_nc->SetBinContent       (ir, spectrum_modes_nc->GetBinContent(ir)        + spectrum_modes[im]->GetBinContent(ir));
      else if(util::IsNuMuBar(m))
	spectrum_modes_numubarcc->SetBinContent(ir, spectrum_modes_numubarcc->GetBinContent(ir) + spectrum_modes[im]->GetBinContent(ir));
      else //it's nue/nuebar/oscnue & CC
	spectrum_modes_nuecc->SetBinContent    (ir, spectrum_modes_nuecc->GetBinContent(ir)     + spectrum_modes[im]->GetBinContent(ir));
    }//ir
    total_integral += spectrum_modes[im]->Integral();
  }//im

  cout << endl << spectrum_modes_nc->Integral() + spectrum_modes_nuecc->Integral()
    + spectrum_modes_numubarcc->Integral() + spectrum_modes_numuccnqe->Integral() + spectrum_modes_numuccqe->Integral()
       << endl << spectrum_bf->Integral() 
       << endl << total_integral << endl << endl;
    
  //style
  set_line_style  (spectrum_bf,    kBlack, 2,   1);
  set_line_style  (spectrum_obs,   kBlack, 1.0, 1);
  set_marker_style(spectrum_obs,   kBlack, 1.0, 8);
  set_line_style(spectrum_modes_numubarcc, kBlack, 1, 1);
  set_line_style(spectrum_modes_numuccnqe, kBlack, 1, 1);
  set_line_style(spectrum_modes_numuccqe,  kBlack, 1, 1);
  set_line_style(spectrum_modes_nuecc,     kBlack, 1, 1);
  set_line_style(spectrum_modes_nc,        kBlack, 1, 1);
  set_fill_style(spectrum_modes_numubarcc, kGreen,  1001);
  set_fill_style(spectrum_modes_numuccnqe, kBlack,  3004);
  set_fill_style(spectrum_modes_numuccqe,  kWhite,  1001);
  set_fill_style(spectrum_modes_nuecc,     kBlue,   1001);
  set_fill_style(spectrum_modes_nc,        kRed,    1001);
  //make stack
  THStack * stack = new THStack("","");
  stack->Add(spectrum_modes_nc);
  stack->Add(spectrum_modes_nuecc);
  stack->Add(spectrum_modes_numubarcc);
  stack->Add(spectrum_modes_numuccnqe);
  //stack->Add(spectrum_modes_numuccqe);
  //prelim draw
  stack->Draw();
  //axis range
  spectrum_obs->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  stack->       GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  //axis title
  stack->       GetXaxis()->SetTitle("E_{reco} / GeV");
  stack->       GetYaxis()->SetTitle("Events per bin");
  spectrum_obs->GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_obs->GetYaxis()->SetTitle("Events per bin");
  //legend
  delete l;
  l = new TLegend(0.5,0.58,0.9,0.93);
  set_legend_style(l);
  l->AddEntry(spectrum_obs,  "Run 1-4 data","P");
  l->AddEntry(spectrum_bf,   "Best fit (total)","l");
  l->AddEntry(spectrum_modes_numuccqe,  "Best fit (#nu_{#mu} CCQE)", "f");
  l->AddEntry(spectrum_modes_numuccnqe, "Best fit (#nu_{#mu} CCnonQE)", "f");
  l->AddEntry(spectrum_modes_numubarcc, "Best fit (#bar{#nu_{#mu}} CC)", "f");
  l->AddEntry(spectrum_modes_nuecc,     "Best fit (#nu_{e}, #bar{#nu_{e}} CC)", "f");
  l->AddEntry(spectrum_modes_nc,        "Best fit (NC)", "f");
  //draw
  spectrum_obs->Draw("PE0");
  stack->       Draw("same");
  spectrum_bf-> Draw("same"); //to make sure the total & indivual modes disagree (stupid e scale)
  spectrum_obs->Draw("PE0same");
  l->Draw();
  save_canvas(c,"spectrum_components",tag);
  //make per 50 mev
  TH1D * spectrumper50_modes_nc = spectrum_to_per_50mev(spectrum_modes_nc);
  TH1D * spectrumper50_modes_nuecc = spectrum_to_per_50mev(spectrum_modes_nuecc);
  TH1D * spectrumper50_modes_numubarcc = spectrum_to_per_50mev(spectrum_modes_numubarcc);
  TH1D * spectrumper50_modes_numuccnqe = spectrum_to_per_50mev(spectrum_modes_numuccnqe);
  TH1D * spectrumper50_modes_numuccqe = spectrum_to_per_50mev(spectrum_modes_numuccqe);
  THStack * stackper50 = new THStack("","");
  stackper50->Add(spectrumper50_modes_nc);
  stackper50->Add(spectrumper50_modes_nuecc);
  stackper50->Add(spectrumper50_modes_numubarcc);
  stackper50->Add(spectrumper50_modes_numuccnqe);
  //stackper50->Add(spectrumper50_modes_numuccqe);
  //draw
  stackper50->Draw();
  stackper50->       GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  stackper50->       GetXaxis()->SetTitle("E_{reco} / GeV");
  stackper50->       GetYaxis()->SetTitle("Events per 0.050 GeV");
  set_line_style  (spectrumper50_bf,    kBlack, 2,   1);
  set_line_style  (spectrumper50_obs,   kBlack, 1.0, 1);
  set_marker_style(spectrumper50_obs,   kBlack, 1.0, 8);
  spectrumper50_obs->Draw("PE0");
  stackper50->       Draw("same");
  spectrumper50_bf-> Draw("same"); //to make sure the total & indivual modes disagree (stupid e scale)
  spectrumper50_obs->Draw("PE0same");
  l->Draw();
  save_canvas(c,"spectrumper50_components",tag);


  //print out the ratios below, at, above the osc peak
  double nobs_below = spectrum_obs  ->Integral(1,                                          spectrum_obs->FindBin(0.5-conv::kOneMeV));
  double nobs_at    = spectrum_obs  ->Integral(spectrum_obs  ->FindBin(0.5+conv::kOneMeV), spectrum_obs->FindBin(1.0-conv::kOneMeV));
  double nobs_above = spectrum_obs  ->Integral(spectrum_obs  ->FindBin(1.0+conv::kOneMeV), spectrum_obs->GetNbinsX());
  double nobs_total = spectrum_obs  ->Integral(1,                                          spectrum_obs->GetNbinsX());
  double nexp_below = spectrum_noosc->Integral(1,                                          spectrum_noosc->FindBin(0.5-conv::kOneMeV));
  double nexp_at    = spectrum_noosc->Integral(spectrum_noosc->FindBin(0.5+conv::kOneMeV), spectrum_noosc->FindBin(1.0-conv::kOneMeV));
  double nexp_above = spectrum_noosc->Integral(spectrum_noosc->FindBin(1.0+conv::kOneMeV), spectrum_noosc->GetNbinsX());
  double nexp_total = spectrum_noosc->Integral(1,                                          spectrum_noosc->GetNbinsX());

  cout << "      \t" << "nobs"     << "\t" << "nexp(noosc)" << "\t"   << "ratio%" << endl;
  cout << "below \t" << nobs_below << "\t" << nexp_below    << "\t\t" << (nobs_below / nexp_below) * 100 << endl;
  cout << "at    \t" << nobs_at    << "\t" << nexp_at       << "\t\t" << (nobs_at    / nexp_at   ) * 100 << endl;
  cout << "above \t" << nobs_above << "\t" << nexp_above    << "\t\t" << (nobs_above / nexp_above) * 100 << endl;
  cout << "total \t" << nobs_total << "\t" << nexp_total    << "\t\t" << (nobs_total / nexp_total) * 100 << endl;

  //////////////////////////////////////////////////////////////////////
  //PULL PLOT

  FitParamList * fparm = 0;
  fin.GetObject("fparm",fparm);
  cout << fparm->GetParListName() << endl;

  //count free parameters
  int np = fparm->GetNumberParameters();
  int npfree = 0;
  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    if(fparm->IsNuisance(p))
      npfree++;
  }
  cout << npfree << " free parameters (" << np << " total)" << endl;
  //make the histogram
  TH1D * pulls    = new TH1D("","",npfree,0,npfree);
  TH1D * values   = new TH1D("","",npfree,0,npfree);
  TH1D * errratio = new TH1D("","",npfree,0,npfree); // postfit/prefit error
  //get deviation from best fit (in sigma)
  for(int ip = 0, ipfree = 1; ip < np; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    if(!fparm->IsNuisance(p)) continue;
    TString label;
    if(ipfree%2)
      label = Form("%s                         ",util::AsStringPaperRooTex(p).c_str());
    else
      label = Form("%s",util::AsStringPaperRooTex(p).c_str());
    pulls->   GetXaxis()->SetBinLabel(ipfree, label.Data());
    values->  GetXaxis()->SetBinLabel(ipfree, label.Data());
    errratio->GetXaxis()->SetBinLabel(ipfree, label.Data());
    //get best fit info
    double best = fparm->BestFitValue(p);
    double error = 0;
    if(fparm->HesseErrorsSet())
      error = fparm->HesseError(p);
    else
      error = fparm->MigradError(p);
    //osc params
    if(p == conv::kPhypMixing13) {
      if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
	pulls-> SetBinContent(ipfree, (best - conv::kGlobFitSinsq2Theta13) / error);
	values->SetBinContent(ipfree, (best - conv::kGlobFitSinsq2Theta13) / (conv::kGlobFitErrSinsq2Theta13 * conv::kGlobFitSinsq2Theta13));
	errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrSinsq2Theta13 * conv::kGlobFitSinsq2Theta13));
      }
      else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
	pulls-> SetBinContent(ipfree, (best - conv::kGlobFitSinsqTheta13) / error);
	values->SetBinContent(ipfree, (best - conv::kGlobFitSinsqTheta13) / (conv::kGlobFitErrSinsqTheta13 * conv::kGlobFitSinsqTheta13));
	errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrSinsqTheta13 * conv::kGlobFitSinsqTheta13));
	cout << util::AsString(p) << "\t" << error << "\t" << conv::kGlobFitErrSinsqTheta13 << "\t" << error / conv::kGlobFitErrSinsqTheta13 << endl;
      }
    }
    else if(p == conv::kPhypMixing12) {
      if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
	pulls-> SetBinContent(ipfree, (best - conv::kGlobFitSinsq2Theta12) / error);
	values->SetBinContent(ipfree, (best - conv::kGlobFitSinsq2Theta12) / (conv::kGlobFitErrSinsq2Theta12 * conv::kGlobFitSinsq2Theta12));
	errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrSinsq2Theta12 * conv::kGlobFitSinsq2Theta12));
      }
      else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
	pulls-> SetBinContent(ipfree, (best - conv::kGlobFitSinsqTheta12) / error);
	values->SetBinContent(ipfree, (best - conv::kGlobFitSinsqTheta12) / (conv::kGlobFitErrSinsqTheta12 * conv::kGlobFitSinsqTheta12));
	errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrSinsqTheta12 * conv::kGlobFitSinsqTheta12));
	cout << util::AsString(p) << "\t" << error << "\t" << conv::kGlobFitErrSinsqTheta12 << "\t" << error / conv::kGlobFitErrSinsqTheta12 << endl;
      }
    }
    else if(p == conv::kPhypCPV    ) {
      pulls-> SetBinContent(ipfree, (best - conv::kGlobFitDeltaCP) / error);
      values->SetBinContent(ipfree, (best - conv::kGlobFitDeltaCP) / conv::kpi);//conv::kGlobFitErrDeltaCP);
      errratio->SetBinContent(ipfree, error / conv::kpi);
      cout << best << endl;
      cout << error << endl;
      cout << 2*conv::kpi << endl;
    }
    else if(p == conv::kPhypDm12sq) {
      pulls-> SetBinContent(ipfree, (best - conv::kGlobFitDm12sq) / error);
      values->SetBinContent(ipfree, (best - conv::kGlobFitDm12sq) / (conv::kGlobFitErrDm12sq * conv::kGlobFitDm12sq));
      errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrDm12sq * conv::kGlobFitDm12sq));
	cout << util::AsString(p) << "\t" << error << "\t" << conv::kGlobFitErrDm12sq << "\t" << error / conv::kGlobFitErrDm12sq << endl;
    }
    else if(util::IsBANFF(p)) {
      pulls-> SetBinContent(ipfree, (best - skexp->GetBANFF()->BestFitValue(p)) / error);
      values->SetBinContent(ipfree, (best - skexp->GetBANFF()->BestFitValue(p)) / skexp->GetBANFF()->OneSigmaError(p));
      errratio->SetBinContent(ipfree, error / skexp->GetBANFF()->OneSigmaError(p));
    }
    else if(util::IsXSec(p)) {
      //"nominal best fit" is 0
      //"1 sigma error" is 1
      pulls-> SetBinContent(ipfree, (best - 0) / error);
      values->SetBinContent(ipfree, (best - 0) / 1);
      errratio->SetBinContent(ipfree, error / 1);
    }
    else if(p == conv::kNuisSuperkEnergyScale) {
      //"nominal best fit" is 0
      //"1 sigma error" is 1
      pulls-> SetBinContent(ipfree, (best - 0) / error);
      values->SetBinContent(ipfree, (best - 0) / 1);
      errratio->SetBinContent(ipfree, error / 1);
    }
    else if(util::IsSKDetFSI(p)) {
      //"nominal best fit" is 1
      pulls-> SetBinContent(ipfree, (best - 1) / error);
      values->SetBinContent(ipfree, (best - 1) / skexp->GetSuperkDetectorFSI()->OneSigmaError(p));
      errratio->SetBinContent(ipfree, error / skexp->GetSuperkDetectorFSI()->OneSigmaError(p));
    }
    else {
      cerr << "Unknown pull for parameter " << util::AsStringCompact(p) << endl;
    }
    cout << util::AsString(p) << "\t" << errratio->GetBinContent(ipfree) << endl;
    ipfree++;
  }//ip
  //style
  pulls->   GetYaxis()->SetRangeUser(-1,1);
  pulls->   GetXaxis()->LabelsOption("v");
  values->  GetYaxis()->SetRangeUser(-1,1);
  values->  GetXaxis()->LabelsOption("v");
  errratio->GetYaxis()->SetRangeUser( 0,2);
  errratio->GetXaxis()->LabelsOption("v");
  TLine * lpull = new TLine(0,0,pulls->GetXaxis()->GetLast(),0);
  TLine * lerr  = new TLine(0,1,pulls->GetXaxis()->GetLast(),1);
  //title
  pulls-> GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma_{bf}");
  values->GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma_{input}");
  errratio->GetYaxis()->SetTitle("#sigma_{bf} / #sigma_{input}");
  //draw
  TCanvas * cpull = new TCanvas();
  cpull->SetGridx(1);
  cpull->SetBottomMargin(0.4);
  cpull->cd();
  pulls->Draw("P");
  lpull->Draw();
  save_canvas(cpull,"pulls",tag);
  cpull->cd();
  values->Draw("P");
  lpull->Draw();
  save_canvas(cpull,"values",tag);
  cpull->cd();
  errratio->Draw("P");
  lerr->Draw();
  save_canvas(cpull,"err_ratio",tag);

  //save the pull histogram
  TFile fout_pull(Form("%s_pulls.root",tag),"RECREATE");
  fout_pull.cd();
  pulls->Write("pulls");
  values->Write("values");
  errratio->Write("pulls");
  fout_pull.Close();

  //count number of parameters with |pull| > 0.1
  int nbigpull = 0;
  for(int i = 1; i <= pulls->GetNbinsX(); i++) {
    if(TMath::Abs(pulls->GetBinContent(i)) > 0.1) {
      cout << pulls->GetXaxis()->GetBinLabel(i) << " has a large pull: " << pulls->GetBinContent(i) << endl;
      nbigpull++;
    }
  }
  cout << endl << nbigpull << " parameters have |pull| > 0.1" << endl << endl;

  for(int i = 1; i <= values->GetNbinsX(); i++) {
    if(TMath::Abs(values->GetBinContent(i)) > 0.1) {
      cout << values->GetXaxis()->GetBinLabel(i) << " has a large value: " << values->GetBinContent(i) << endl;
      nbigpull++;
    }
  }

  //draw again
  pulls-> GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma");
  pulls-> GetYaxis()->SetRangeUser(-1.1,1.1);
  set_marker_style(values, kRed, values->GetMarkerSize(),25);
  //set_marker_style(values,   kRed, 1.0, 8);
  cpull->cd();
  pulls->Draw("P");
  values->Draw("Psame");
  lpull->Draw();
  delete l;
  l = new TLegend(0.3,0.03,0.5,0.13);
  //l->SetFillColor(kWhite);
  //l->SetTextSize(0.05);
  set_legend_style(l);
  l->AddEntry(pulls,  "#sigma = #sigma_{bf}",  "P");
  l->AddEntry(values, "#sigma = #sigma_{input}","P");
  l->Draw();
  save_canvas(cpull,"pulls_and_values",tag);

  delete cpull;

  //////////////////////////////////////////////////////////////////////
  //CHISQ & DCHISQ

  //get the histograms
  TH2D * hchisq = 0, * hdchisq = 0;
  fin.GetObject("chisq",  hchisq);
  fin.GetObject("dchisq", hdchisq);
  if(smooth)
    hdchisq->Smooth();
  //axis range
  hchisq-> GetXaxis()->SetRangeUser(mixing_min, mixing_max);
  hchisq-> GetYaxis()->SetRangeUser(dm_min,     dm_max);
  hdchisq->GetXaxis()->SetRangeUser(mixing_min, mixing_max);
  hdchisq->GetYaxis()->SetRangeUser(dm_min,     dm_max);
  //axis title
  hchisq-> GetXaxis()->SetTitle(util::AsStringPaperRooTex(conv::kPhypMixing23).c_str());
  //hchisq-> GetYaxis()->SetTitle(util::AsStringPaperRooTex(conv::kPhypDm23sq).c_str());
  hchisq-> GetYaxis()->SetTitleOffset(1.50);
  hdchisq->GetXaxis()->SetTitle(util::AsStringPaperRooTex(conv::kPhypMixing23).c_str());
  //hdchisq->GetYaxis()->SetTitle(util::AsStringPaperRooTex(conv::kPhypDm23sq).c_str());
  hdchisq->GetYaxis()->SetTitleOffset(1.50);
  switch(dmsqmode) {
  case 0:
    hchisq ->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
    hdchisq->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
    break;
  case 1:
    hchisq ->GetYaxis()->SetTitle("|#Deltam^{2}_{31}| (eV^{2}/c^{4})");
    hdchisq->GetYaxis()->SetTitle("|#Deltam^{2}_{31}| (eV^{2}/c^{4})");
    break;
  default:
    hchisq ->GetYaxis()->SetTitle("");
    hdchisq->GetYaxis()->SetTitle("");
    break;
  }

  //draw
  TCanvas * cchisq = new TCanvas();
  cchisq->SetRightMargin(0.12);
  cchisq->SetLeftMargin(0.15);
  cchisq->cd();
  hchisq-> Draw("colz");
  save_canvas(cchisq,"chisq",tag);
  cchisq->cd();
  hdchisq->Draw("colz");
  save_canvas(cchisq,"dchisq",tag);
  hdchisq->GetZaxis()->SetRangeUser(0,60);
  hdchisq->Draw("colz");
  save_canvas(cchisq,"dchisq_maxz60",tag);


  //////////////////////////////////////////////////////////////////////
  //CONTOUR

  //param space
  double mixing23_min   = 0.0;
  double mixing23_max   = 1.5;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 6E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int)( ((mixing23_max - mixing23_min) / mixing23_step) + 0.01);
  int dm23sq_nbins   = 1 + (int)( ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  ) + 0.01);
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
  osc_param_space->GetXaxis()->SetTitle(util::AsStringPaperRooTex(conv::kPhypMixing23).c_str());
  //osc_param_space->GetYaxis()->SetTitle(util::AsStringPaperRooTex(conv::kPhypDm23sq).c_str());
  switch(dmsqmode) {
  case 0:
    osc_param_space->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
    break;
  case 1:
    osc_param_space->GetYaxis()->SetTitle("|#Deltam^{2}_{31}| (eV^{2}/c^{4})");
    break;
  default:
    osc_param_space->GetYaxis()->SetTitle("");
    break;
  }
  osc_param_space->GetYaxis()->SetTitleOffset(1.30);

  //get the histograms
  TH2D * hcont68 = 0, * hcont90 = 0;
  fin.GetObject("dchisq_cont68", hcont68);
  fin.GetObject("dchisq_cont90", hcont90);
  TH2D * hcont68smooth = new TH2D(*hcont68);
  hcont68smooth->Smooth();
  TH2D * hcont90smooth = new TH2D(*hcont90);
  hcont90smooth->Smooth();
  //axis range
  osc_param_space->GetXaxis()->SetRangeUser(mixing_min, mixing_max);
  osc_param_space->GetYaxis()->SetRangeUser(dm_min,     dm_max);
  //marker
  TMarker * markbf = new TMarker(fparm->BestFitValue(conv::kPhypMixing23), fparm->BestFitValue(conv::kPhypDm23sq), 29);
  //style
  set_marker_style(markbf,  kBlack, 1.8, 29);
  set_line_style  (hcont68, kBlack, 2,   2);
  set_line_style  (hcont90, kBlack, 2,   1);
  //legend
  delete l;
  l = new TLegend(0.2,0.7,0.6,0.9);
  set_legend_style(l);
  l->AddEntry(hcont68,  "Run 1-4 68% C.L.","l");
  l->AddEntry(hcont90,  "Run 1-4 90% C.L.","l");
  l->AddEntry(markbf,   "Best fit","P");
  //draw
  TCanvas * ccontour = new TCanvas();
  ccontour->SetLeftMargin(0.15);
  ccontour->cd();
  osc_param_space->Draw();
  if(!smooth) {
    hcont68->Draw("cont3same");
    hcont90->Draw("cont3same");
  }
  else {
    hcont68smooth->Draw("cont3same");
    hcont90smooth->Draw("cont3same");
  }
  markbf->Draw();
  l->Draw();
  save_canvas(ccontour,"contour",tag);

  //again, with externals & old t2k
  //read them in. External:
  TFile fminos("$T2KOSCANA/inputs/contours/MINOS_NumuDisappearance_LogLikelihoodSurface_May2013.root");
  TH2D * ll_minos = 0;
  fminos.GetObject("LogLikelihoodSurface",ll_minos);
  ll_minos->SetContour(1);
  ll_minos->SetContourLevel(0,4.61/2);
  TFile fskle("$T2KOSCANA/inputs/contours/contours_skloe.root");
  TGraph * g_skle = 0;
  fskle.GetObject("sk3f_loe_90cl",g_skle);
  g_skle->SetName("g_skle");
  TFile fskz("$T2KOSCANA/inputs/contours/contours_skzen.root");
  TGraph * g_skzen = 0;
  fskz.GetObject("sk3f_23vars_overlay_normal_q13_fixed",g_skzen);
  g_skzen->SetName("g_skzen");
  TFile fice("$T2KOSCANA/inputs/contours/contours_nont2k.root");
  TGraph * g_ice = 0;
  fice.GetObject("icecube20102011_90CL",g_ice);
  g_ice->SetName("g_ice");
  // Old t2k:
  TFile fold("$T2KOSCANA/inputs/contours/contours_oldt2k.root");
  TGraph * g_run12_68_statsyst = 0;
  fold.GetObject("t2k2010_run1p2_68CL_statsyst",g_run12_68_statsyst);
  g_run12_68_statsyst->SetName("g_run12_68_statsyst");
  TGraph * g_run12_90_statsyst = 0;
  fold.GetObject("t2k2010_run1p2_90CL_statsyst",g_run12_90_statsyst);
  g_run12_90_statsyst->SetName("g_run12_90_statsyst");
  TGraph * g_run12_68_stat = 0;
  fold.GetObject("t2k2010_run1p2_68CL_stat",g_run12_68_stat);
  g_run12_68_stat->SetName("g_run12_68_stat");
  TGraph * g_run12_90_stat = 0;
  fold.GetObject("t2k2010_run1p2_90CL_stat",g_run12_90_stat);
  g_run12_90_stat->SetName("g_run12_90_stat");
  TFile ft2krun123("$T2KOSCANA/2012results/NEWRESULTS/FINALFCCONTOURS/fc_oldsurf_newres_5E5_nosub.root");
  TH2D * h_run123_68 = 0;
  TH2D * h_run123_90 = 0;
  ft2krun123.GetObject("hCL68StatSyst",h_run123_68);
  ft2krun123.GetObject("hCL90StatSyst",h_run123_90);
  h_run123_68->SetContour(1);
  h_run123_68->SetContourLevel(0,0.0);
  h_run123_90->SetContour(1);
  h_run123_90->SetContourLevel(0,0.0);
  //style
  set_marker_style(markbf,  kBlack, 1.8, 29);
  //set_line_style(hcont68, kBlack, 2,   2);
  set_line_style(hcont90, kBlack, 2,   1);
  //set_line_style(h_run123_68, kBlack, 2,   2);
  set_line_style(h_run123_90, kBlack, 2,   2);
  set_line_style(ll_minos, kRed, 2, 1);
  set_line_style(g_skzen, kBlue, 2, 1);
  set_line_style(g_skle,  kBlue, 2, 2);
  //set_line_style(g_ice,  kBlue, 2, 2);
  //legend
  delete l;
  l = new TLegend(0.2,0.55,0.6,0.9);
  set_legend_style(l);
  l->AddEntry(hcont90,  "T2K Run 1-4 3#nu 90% C.L.","l");
  l->AddEntry(h_run123_90,  "T2K Run 1-3 3#nu 90% C.L.","l");
  l->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  l->AddEntry(g_skzen, "SK zenith 2012 3#nu 90% CL", "l");
  l->AddEntry(g_skle, "SK L/E 2012 2#nu 90% CL", "l");
  //l->AddEntry(g_ice, "IceCube-79 2013 2#nu 90% CL", "l");
  l->AddEntry(markbf,   "T2K Run 1-4 Best fit","P");
  //axis range
  osc_param_space->GetXaxis()->SetRangeUser(0.85, 1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002, 0.0035);
  ccontour->cd();
  osc_param_space->Draw();
  ll_minos->Draw("cont3same");
  g_skzen->Draw("l");
  g_skle->Draw("l");
  h_run123_90->Draw("cont3same");
  hcont90->Draw("cont3same");
  markbf->Draw();
  l->Draw();
  save_canvas(ccontour,"contour_external",tag);


  //////////////////////////////////////////////////////////////////////
  //PRINT BEST FIT

  cout << endl << endl;
  cout << "Best fit chisq:    " << fparm->GetBestFitChisq() << endl;
  cout << "Chisq surface min: " << hchisq->GetMinimum() << endl;
  cout << "Dchisq min:        " << hchisq->GetMinimum() - fparm->GetBestFitChisq() << endl;
  cout << "Best fit mixing23: " << fparm->BestFitValue(conv::kPhypMixing23) << endl;
  cout << "Best fit dm32sq:   " << fparm->BestFitValue(conv::kPhypDm23sq)   << endl;
  hdchisq->GetZaxis()->UnZoom();
  cout << "Minimum of the dchisq surface: " << hdchisq->GetMinimum() << endl;
  int tempx, tempy, tempz;
  hdchisq->GetMinimumBin(tempx,tempy,tempz);
  cout << " at mixing23:      " << hdchisq->GetXaxis()->GetBinCenter(tempx)
       << " (" << hdchisq->GetXaxis()->GetBinLowEdge(tempx) << " to " << hdchisq->GetXaxis()->GetBinUpEdge(tempx) << ")" << endl;
  cout << " at dm32sq:        " << hdchisq->GetYaxis()->GetBinCenter(tempy)
       << " (" << hdchisq->GetYaxis()->GetBinLowEdge(tempy) << " to " << hdchisq->GetYaxis()->GetBinUpEdge(tempy) << ")" << endl;
  cout << endl;
}
