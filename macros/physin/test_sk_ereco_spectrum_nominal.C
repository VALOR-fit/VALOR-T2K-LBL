//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

void test_sk_ereco_spectrum_nominal
(
  const char * sk_sample            = "1Rmu",                         /* options: "1Rmu", "1Re" */
  const char * hypothesis           = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12       = -99,                            /* */
  double       input_mixing13       = -99,                            /* */
  double       input_mixing23       = -99,                            /* */
  double       input_deltacp        = -99,                            /* */
  double       input_dm12sq         = -99,                            /* */
  double       input_dm23sq         = -99,                            /* */
  bool         disable_banff_tuning = false,                          /* */
  const char * config_file          = "analysis_configuration.xml",   /* */
  bool         print_spectrum       = false                           /* */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Plotting spectra for sample: " << util::AsString(s) << endl;
  cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;

  // Histogram booking
  const int nbins = util::NumErecoBins(s);
  const double * ereco_bin_edges = util::ErecoBinEdges(s);
  TH1D * hErecoTot = new TH1D("","",nbins,ereco_bin_edges);
  TH1D * hErecoMode[conv::knm]; //< contribution from each mode  

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // Disable BANFF tuning if requested
  if(disable_banff_tuning) {
    cout << "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..." << endl;
    skexp.GetBANFF().RestoreNominalMC();
  }

  // Generate prediction
  cout << "- Generating predictions..." << endl;
  skexp.Generate(ph);

  // Get energy spectra
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) << endl;
    TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
    assert(h);
    hErecoMode[im] = new TH1D(*h);
  }
  hErecoTot = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));

  //
  // Draw
  //

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TLegend * legend = new TLegend(0.3, 0.4, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot->Draw();
  hErecoTot->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot->GetYaxis()->SetTitle("Number of events");
  THStack * stack = new THStack("stack","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack->Add(hErecoMode[im]);
  }
  stack->Draw("same");
  legend->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend->Draw();
  c->Update();
  cout << fixed;
  double ntot = hErecoTot->Integral();
  cout << "Total number of events: " 
       << setw(10) << setfill(' ') << setprecision(6) << ntot << endl;
  const int nslices=6;
  double slice_max[nslices] = { 0.50, 0.75, 1.00, 1.25, 1.50, 2.00 };
  for(int islice=0; islice<nslices; islice++) {
    double esmax = slice_max[islice];
    cout << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << esmax << " GeV): " 
         << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(esmax)) << endl;
  }
  for(int im = 0; im < conv::knm; im++) {
     conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
     if(m == conv::kSKmUndefined) continue;
     double nmode = hErecoMode[im]->Integral();
     cout << "Contrib from mode: " 
          << setw(40) << setfill(' ') << util::AsString(m) << " : " 
          << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
          << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  }

  //
  // Save histograms in a ROOT file
  //

  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-%s%s.root",oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : "")) ,"recreate");
  hErecoTot -> Write(util::AsStringCompact(s,conv::kSKmUndefined,"","total").c_str());
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    hErecoMode[im] -> Write(util::AsStringCompact(s,m).c_str());
  }
  f.Close();

  //
  // Print-out numbers requested during the analysis review
  //

  if(print_spectrum) {  
    int nprec=5;
    cout << fixed;
    for(int ibin = 1; ibin <= hErecoTot->GetXaxis()->GetNbins(); ibin++) {
      cout << setw(6) << setfill(' ') << ibin
           << " : "
           << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinLowEdge(ibin)
           << " - "
           << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinUpEdge(ibin)
           << " : "
           << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot->GetBinContent(ibin)
           << endl;
    }//ibin
  }

}
