void set_contour_style(TGraph *contour, Int_t color, Int_t width, Int_t style)
{
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
}

void analyse_output
(
 const char * chisqfile,                                           /* the file with 'chisqwithsyst' in it. E.g. tag-run_0-expts_0_0.root */
 int          id,                                                  /* dataset id, from the xml file. Appears in the filenames, and some histograms */
 const char * path               = "./",                           /* dir with the root files in. fitresults_id.root, mockfit_id.root and chisqfile should all be in the same place */
 const char * config_file        = "analysis_configuration.xml",   /* for generate no-osc spectrum */
 double       spectrum_max_xaxis = 7.0,                            /* max of spectrum axes */
 bool         external_contours  = false                           /* draw external contours? */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas();

  TFile f_dchisq(Form("%s/fitresults_%d.root",path,id));
  TFile f_chisq(Form("%s/%s",path,chisqfile));
  TFile f_mockfit(Form("%s/mockfit_%d.root",path,id));

  cout << Form("%s/fitresults_%d.root",path,id) << endl;

  //spectra
  TH1D * spectrum_obs = 0, * spectrum_bf = 0;
  f_dchisq.GetObject("bf_spectrum",spectrum_bf);
  f_dchisq.GetObject("obs_spectrum",spectrum_obs);
  //make the no-osc spectrum
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
  PhysHypothesis ph(conv::kHypOsc3F);
  skexp->Generate(ph);
  const int      n_ereco_bins    = util::NumErecoBins (conv::kSKsFCFV1RingMuLike);
  const double * ereco_bin_edges = util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike);
  TH1D * spectrum_noosc = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  for(int ibin = 1; ibin <= n_ereco_bins; ibin++)
    spectrum_noosc->SetBinContent(ibin, skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingMuLike)[ibin-1]);
  //now copy the data & no-osc, and divide
  TH1D * spectrum_ratio_bf = new TH1D("","", n_ereco_bins, ereco_bin_edges), * spectrum_ratio_obs = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  spectrum_ratio_bf->Divide(spectrum_bf, spectrum_noosc);
  spectrum_ratio_obs->Divide(spectrum_obs, spectrum_noosc);
  //set the axes
  spectrum_obs->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_bf->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_noosc->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_ratio_obs->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  spectrum_ratio_bf->GetXaxis()->SetRangeUser(0, spectrum_max_xaxis - conv::kOneMeV);
  TLine * spectrum_line = new TLine(0,1, 7,1);
  //styles
  spectrum_bf->SetLineColor(kRed);
  spectrum_obs->SetLineStyle(2);
  spectrum_ratio_bf->SetLineColor(kRed);
  spectrum_ratio_obs->SetLineStyle(2);
  //Now draw
  TCanvas * cspec = new TCanvas();
  cspec->Divide(1,2);
  cspec->cd(1);
  spectrum_noosc->GetYaxis()->SetTitle("Events");
  spectrum_noosc->Draw();
  spectrum_obs->Draw("PEsame");
  spectrum_bf->Draw("same");
  cspec->cd(2);
  spectrum_ratio_obs->GetXaxis()->SetTitle("E_{reco} / GeV");
  spectrum_ratio_obs->GetYaxis()->SetTitle("Ratio");
  spectrum_ratio_obs->Draw("PE");
  spectrum_ratio_bf->Draw("same");
  spectrum_line->Draw("same");
  cspec->SaveAs(Form("spectrum_%d.eps",id));
  delete cspec;

  //pulls
  TH1D * pulls = 0;
  f_dchisq.GetObject("hsystpull",pulls);
  c->cd();
  pulls->GetYaxis()->SetRangeUser(-2,2);
  pulls->Draw();
  c->SaveAs(Form("pulls_%d.eps",id));

  //chisqwithsyst
  c->SetLeftMargin(0.19);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetRightMargin(0.15);
  TH2D * chisq_with_syst = 0;
  f_chisq.GetObject(Form("chisqwithsyst_%d",id),chisq_with_syst);
  chisq_with_syst->GetYaxis()->SetTitleOffset(1.65);
  chisq_with_syst->Draw("colz");
  cout << chisq_with_syst->GetZaxis()->GetTitleOffset() << endl;
  c->SaveAs(Form("chisq_with_syst_%d.eps",id));

  //dchisq's & contours
  TH2D * dchisq = 0;
  f_dchisq.GetObject("hdchisq",dchisq);
  dchisq->SetLineColor(kRed);
  //dchisq
  dchisq->SetContour(255);
  dchisq->GetYaxis()->SetTitleOffset(1.65);
  TPaletteAxis *palette = (TPaletteAxis*)dchisq->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  dchisq->GetZaxis()->SetRangeUser(0,dchisq->GetMaximum());
  dchisq->Draw("colz");
  c->SaveAs(Form("dchisq_%d.eps",id));
  cout << dchisq->GetZaxis()->GetTitleOffset() << endl;

  //copy to be able to set the same color for the different (68%, 90%) contours
  TH2D * dchisq68 = new TH2D(*dchisq), * dchisq90 = new TH2D(*dchisq);
  double crit68[1] = {2.30}, crit90[1] = {4.61};
  dchisq68->SetContour(1, crit68);
  //  dchisq68->SetContourLevel(1,2.30);
  dchisq90->SetContour(1, crit90);
  //dchisq90->SetContourLevel(1,4.61);
  dchisq68->SetLineStyle(2);
  //Make the tmarker for the bestfit point
  TVectorD * vec_results = 0;
  f_dchisq.GetObject("oscparams",vec_results);
  TVectorD & vec = * vec_results;
  TMarker *bf = new TMarker(vec[0], vec[2], 29);
  bf->SetMarkerSize(1.8);
  bf->SetMarkerColor(1);
  //legend
  TLegend * leg = new TLegend(0.215517,0.190678,0.594828,0.483051);
  leg->SetFillColor(0);
  leg->SetLineStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(dchisq90, "Const #Delta#chi^2 90\% ", "l");
  leg->AddEntry(dchisq68, "Const #Delta#chi^2 68\% ", "l");
  leg->AddEntry(bf, "Best fit", "p");
  // param space
  c->SetRightMargin(0.05);
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.6;
  double mixing23_max   = 1.;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1.0E-3;
  double dm23sq_max     = 4.0E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + (int) ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
  osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  //osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.65);
  osc_param_space->Draw();
  //actually draw
  dchisq90->Draw("cont3same");
  dchisq68->Draw("cont3same");
  bf->Draw("Psame");
  //add externals
  TGraph * minos2012_90cl = 0;
  TGraph * skle2012_90cl  = 0;
  TGraph * skzen2012_90cl = 0;
  if(external_contours) {
    const char * extcontourfile = Form("%s/inputs/contours/contours_nont2k.root", gSystem->Getenv("T2KOSCANA"));
    TFile f(extcontourfile,"read");
    minos2012_90cl = new TGraph(* (TGraph*)f.Get("minos2012_combined_90CL") );
    skle2012_90cl  = new TGraph(* (TGraph*)f.Get("sk2012_le2nu_90CL")       );
    skzen2012_90cl = new TGraph(* (TGraph*)f.Get("sk2012_zenith3nu_90CL")   );
    set_contour_style(minos2012_90cl, kBlack,     1, kSolid);
    set_contour_style(skle2012_90cl,  kMagenta-3, 1, kSolid);
    set_contour_style(skzen2012_90cl, kGreen+2,   1, kSolid);
    minos2012_90cl->Draw("l");
    skle2012_90cl ->Draw("l");
    skzen2012_90cl->Draw("l");
    f.Close();
    minos2012_90cl->SetName("minos2012_90cl");
    skzen2012_90cl->SetName("skzen2012_90cl");
    skle2012_90cl ->SetName("skle2012_90cl");
    leg->AddEntry("minos2012_90cl", "MINOS 2012 (10.71#times10^{20} POT + 37.9 kton*yr), 90% CL","L");
    leg->AddEntry("skzen2012_90cl", "SK1-4 zenith angle, 90% CL","L");
    leg->AddEntry("skle2012_90cl",  "SK1-4 L/E, 90% CL","L");
  }
  leg->Draw();
  c->SaveAs(Form("contours_%d.eps",id));

  delete c;
}
