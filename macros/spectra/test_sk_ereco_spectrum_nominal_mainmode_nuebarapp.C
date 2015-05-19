//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

using namespace std;
using namespace valor;
using namespace valor::conv;

void test_sk_ereco_spectrum_nominal_mainmode_nuebarapp
(
 const char * tag                   = "prova",
 
 const char * sk_sample            = "RHC1Re",                         /* options: "1Rmu", "1Re", "MultiRmu" */
 const char * hypothesis           = "3f_osc_normbar",                 /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_osc_normbar, 3f_oscbar_normbar" */
 
 // double       input_mixing12       = 0.,                          
 // double       input_mixing13       = 0.,                          
 // double       input_mixing23       = 0.,                          
 // double       input_deltacp        = 0.,                          
 // double       input_dm12sq         = 0.,                          
 // double       input_dm23sq         = 0.,                          
 // double       input_betabar        = 0.,

 // double       input_mixing12       = -99,                         
 // double       input_mixing13       = -99,                         
 // double       input_mixing23       = -99,                         
 // double       input_deltacp        = -99,                         
 // double       input_dm12sq         = -99,                         
 // double       input_dm23sq         = -99,                       
 // double       input_betabar        = -99,

 // // Normal Hierarchy
 // double       input_mixing12       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = conv::kT2KReacBestFitSinsqTheta13_NH,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = conv::kT2KReacBestFitSinsqTheta23_NH, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = conv::kT2KReacBestFitDeltaCP_NH, //-99,    /* deltacp, */
 // double       input_dm12sq         = -99,    /* Dm^2_{12},  */
 // double       input_dm23sq         = conv::kT2KReacBestFitDmsq23_NH, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar       = -99,

 // Inverted Hierarchy
 double       input_mixing12       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13       = conv::kT2KReacBestFitSinsqTheta13_IH,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23       = conv::kT2KReacBestFitSinsqTheta23_IH, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp        = conv::kT2KReacBestFitDeltaCP_IH, //-99,    /* deltacp, */
 double       input_dm12sq         = -99,    /* Dm^2_{12},  */
 double       input_dm23sq         = conv::kT2KReacBestFitDmsq23_IH, //-99,    /* |Dm^2_{23}| */
 double       input_betabar       = -99,
  
 bool         disable_banff_tuning = false,                          /* */

 const char * config_file      = "analysis_configuration_run1to6_test_app.xml",      
 //const char * config_file      = "analysis_configuration_run1to6_test_app_fullstat.xml",      
 //const char * config_file      = "analysis_configuration_run1to6_test_prefit_app.xml",      
   
 bool         print_spectrum       = false,                          /* */
 const char * dmsqscheme          = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy           = "inverted"   ,
 const char * mixparamchoice      = "sinsqtheta",
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

  bool plot_per_50MeV = true
)
{
  using namespace std;
  using namespace valor;
  
  const char * plot_format = "eps";

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

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
  if(input_betabar < -1) input_betabar = conv::kGlobFitBetaBar;
 
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FNormBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypBetaBar,    input_betabar);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FBarNormBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypBetaBar,    input_betabar);
  }
  
  
  
  string oscparams = util::AsString(ph.Hypothesis(),
				    input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq,input_betabar);

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

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
  TH1D * hErecoMode_draw[conv::knm]; //< contribution from each mode  
  const int nbins_coarse = util::NumVariableXsEtrueBins(s,conv::kSKmNuMuCCQE);
  const double * ereco_bin_edges_coarse = util::VariableXsEtrueBinEdges(s,conv::kSKmNuMuCCQE);
  TH1D * hErecoTot_coarse;
  TH1D * hErecoMode_coarse[conv::knm]; //< contribution from each mode
  TH1D * hErecoMode_draw_coarse[conv::knm]; //< contribution from each mode

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // Disable BANFF tuning if requested
  if(disable_banff_tuning) {
    cout << "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..." << endl;
    skexp.GetBANFF()->RestoreNominalMC();
  }

  // Generate prediction
  cout << "- Generating predictions..." << endl;
  skexp.GetSuperkMC()->SetXMEC(1);
  skexp.Generate(ph);


  // Get energy spectra
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) << endl;
    TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
    assert(h);
    hErecoMode[im] = new TH1D(*h);
    hErecoMode_coarse[im] = (TH1D)hErecoMode[im]->Rebin(nbins_coarse,Form("%s_coarse",hErecoMode[im]->GetName()),ereco_bin_edges_coarse);
    
    if (plot_per_50MeV) {
      hErecoMode_draw[im]        = spectrum_to_per_50mev(hErecoMode[im]);
      hErecoMode_draw_coarse[im] = spectrum_to_per_50mev(hErecoMode_coarse[im]);
    }
    else {
      hErecoMode_draw[im]        = hErecoMode[im];
      hErecoMode_draw_coarse[im] = hErecoMode_coarse[im];
    }
      
  }
  hErecoTot = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  hErecoTot_coarse = (TH1D)hErecoTot->Rebin(nbins_coarse,Form("%s_coarse",hErecoTot->GetName()),ereco_bin_edges_coarse);
  
  hErecoTot->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot->GetYaxis()->SetTitle("Number of events");
  hErecoTot_coarse->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot_coarse->GetYaxis()->SetTitle("Number of events");
  
  
  
  //can save this to file &/or print it
  stringstream streamer;

  streamer << fixed;
  double ntot = hErecoTot->Integral();
  double ntot_coarse = hErecoTot->Integral();
  if(TMath::Abs(ntot - ntot_coarse) > conv::kASmallNum) {
    cerr << "Integrals are different between differently binned histograms "
	 << "50MeV: " << ntot
	 << "BANFF binning: " << ntot_coarse << endl;
    exit(1);
  }

  streamer << "Total number of events: " 
       << setw(10) << setfill(' ') << setprecision(6) << ntot << endl;

  const int nslices=6;
  double slice_max[nslices] = { 0.50, 0.75, 1.00, 1.25, 1.50, 2.00 };
  for(int islice=0; islice<nslices; islice++) {
    double esmax = slice_max[islice];
    streamer << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << esmax << " GeV): " 
         << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(esmax)) << endl;
  }
  for(int im = 0; im < conv::knm; im++) {
     conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
     if(m == conv::kSKmUndefined) continue;
     double nmode = hErecoMode[im]->Integral();
     streamer << "Contrib from mode: " 
          << setw(40) << setfill(' ') << util::AsString(m) << " : " 
          << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
          << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  }


  // Build mode/flavor

  TH1D * hErecoMode_NC = new TH1D("hErecoMode_NC","hErecoMode_NC",nbins,ereco_bin_edges);
  TH1D * hErecoMode_BkgOther = new TH1D("hErecoMode_BkgOther","hErecoMode_BkgOther",nbins,ereco_bin_edges);
  TH1D * hErecoMode_NumuBarToNueBar = new TH1D("hErecoMode_NumuBarToNueBar","hErecoMode_NumuBarToNueBar",nbins,ereco_bin_edges);
  TH1D * hErecoMode_NumuToNue = new TH1D("hErecoMode_NumuToNue","hErecoMode_NumuToNue",nbins,ereco_bin_edges);

  // flavor only
  TH1D * hErecoFlav_Numu = new TH1D("hErecoFlav_Numu","hErecoFlav_Numu",nbins,ereco_bin_edges);
  TH1D * hErecoFlav_NumuBar = new TH1D("hErecoFlav_NumuBar","hErecoFlav_NumuBar",nbins,ereco_bin_edges);
  TH1D * hErecoFlav_Nue = new TH1D("hErecoFlav_Nue","hErecoFlav_Nue",nbins,ereco_bin_edges);
  TH1D * hErecoFlav_NueBar = new TH1D("hErecoFlav_NueBar","hErecoFlav_NueBar",nbins,ereco_bin_edges);

  hErecoFlav_Numu->SetFillStyle(0);
  hErecoFlav_NumuBar->SetFillStyle(0);
  hErecoFlav_Nue->SetFillStyle(0);
  hErecoFlav_NueBar->SetFillStyle(0);
  
  hErecoMode_NC->SetFillStyle(0);
  hErecoMode_BkgOther->SetFillStyle(0);
  hErecoMode_NumuBarToNueBar->SetFillStyle(0);
  hErecoMode_NumuToNue->SetFillStyle(0);

  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
    string modestr = util::AsStringCompact(m);
    
    // Fill by flavor
    if(util::IsNuE(m) || util::IsOscNuE(m))
      hErecoFlav_Nue->Add(skexp.CurrPredictedErecoSpectrum(s,m));
    else if(util::IsNuEBar(m) || util::IsOscNuEBar(m))
      hErecoFlav_NueBar->Add(skexp.CurrPredictedErecoSpectrum(s,m));
    else if(util::IsNuMuBar(m) || util::IsOscNuMuBar(m))
      hErecoFlav_NumuBar->Add(skexp.CurrPredictedErecoSpectrum(s,m));
    else if(util::IsNuMu(m) || util::IsOscNuMu(m))
      hErecoFlav_Numu->Add(skexp.CurrPredictedErecoSpectrum(s,m));
    else{
      cout << modestr << endl;
      cout << "No flavors are selected!!!" << endl;
      exit(1);
    }
    
    cout << m << " " << modestr << " --> ";
    if(util::IsCCQE(m)) cout << "IsCCQE" << endl;
    else if(util::IsMEC(m)) cout << "IsMEC" << endl;
    else if(util::IsCCnQE(m)) cout << "IsCCnQE" << endl;
    else if(util::IsNC(m)) cout << "IsNC" << endl;
    else{
      cout << "No modes!!!" << endl;
      exit(1);
    }
    cout << endl;
    
    // Fill by flavor/mode
    if(util::IsNC(m))
      hErecoMode_NC->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // NC
    else if(util::IsOscNuEBar(m))
      hErecoMode_NumuBarToNueBar->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // NumuBar->NueBar
    else if(util::IsOscNuE(m))
      hErecoMode_NumuToNue->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // Numu->Nue
    else
      hErecoMode_BkgOther->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // Bkg Other
    
    //else{
    //cout << modestr << endl;
    //cout << "No mode / flavors are selected!!!" << endl;
    //exit(1);
    //}

  } // loop over modes

  // Draw Histo by mode/flavor

  //gStyle->SetHatchesSpacing();
  //gStyle->SetHatchesLineWidth(3.);

  if(plot_per_50MeV){
    hErecoMode_NC = spectrum_to_per_50mev(hErecoMode_NC); 
    hErecoMode_BkgOther = spectrum_to_per_50mev(hErecoMode_BkgOther); 
    hErecoMode_NumuBarToNueBar = spectrum_to_per_50mev(hErecoMode_NumuBarToNueBar); 
    hErecoMode_NumuToNue = spectrum_to_per_50mev(hErecoMode_NumuToNue); 
  }

  hErecoMode_NC->SetFillStyle(1001);
  hErecoMode_NC->SetLineColor(kBlack);
  //hErecoMode_NC->SetLineColor(kRed+2);
  hErecoMode_NC->SetFillColor(kRed);
  hErecoMode_NC->SetLineWidth(2);
  
  hErecoMode_BkgOther->SetFillStyle(1001);
  hErecoMode_BkgOther->SetLineColor(kBlack);
  //hErecoMode_BkgOther->SetLineColor(kBlue+2);
  hErecoMode_BkgOther->SetFillColor(kBlue);
  hErecoMode_BkgOther->SetLineWidth(2);

  hErecoMode_NumuBarToNueBar->SetFillStyle(1001);
  hErecoMode_NumuBarToNueBar->SetLineColor(kBlack);
  hErecoMode_NumuBarToNueBar->SetFillColor(kGreen);
  hErecoMode_NumuBarToNueBar->SetLineWidth(2);

  hErecoMode_NumuToNue->SetFillStyle(3004);
  hErecoMode_NumuToNue->SetLineColor(kBlack);
  hErecoMode_NumuToNue->SetFillColor(kGray+3);
  hErecoMode_NumuToNue->SetLineWidth(2);
  //hErecoMode_NumuToNue->SetFillStyle(3004);

  hErecoMode_NC->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NC->GetYaxis()->SetTitle("Events per 0.050 GeV");
  hErecoMode_BkgOther->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_BkgOther->GetYaxis()->SetTitle("Events per 0.050 GeV");
  hErecoMode_NumuToNue->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NumuToNue->GetYaxis()->SetTitle("Events per 0.050 GeV");
  hErecoMode_NumuBarToNueBar->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NumuBarToNueBar->GetYaxis()->SetTitle("Events per 0.050 GeV");

  THStack *sErecoTotFlavMode = new THStack("sErecoTotFlavMode","sErecoTotFlavMode");
  sErecoTotFlavMode->Add(hErecoMode_NC);
  sErecoTotFlavMode->Add(hErecoMode_BkgOther);
  sErecoTotFlavMode->Add(hErecoMode_NumuToNue);
  sErecoTotFlavMode->Add(hErecoMode_NumuBarToNueBar);

  TCanvas *cErecoTotFlavMode = new TCanvas("cErecoTotFlavMode","cErecoTotFlavMode");
  sErecoTotFlavMode->Draw();
  sErecoTotFlavMode->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  sErecoTotFlavMode->GetYaxis()->SetTitle("Events per 0.050 GeV");
  cErecoTotFlavMode->Modified();

  gStyle->SetOptStat(0);
  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu} mode, e-like Total","L");
  leg->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","F");
  leg->AddEntry(hErecoMode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","F");
  leg->AddEntry(hErecoMode_NC,"NC","F");
  leg->AddEntry(hErecoMode_BkgOther,"Bkg other","F");
  leg->Draw();


  // Draw Histo by flavor

  hErecoFlav_Nue->SetFillStyle(1001);
  hErecoFlav_Nue->SetFillColor(kYellow);
  hErecoFlav_Nue->SetLineColor(kYellow);

  hErecoFlav_NueBar->SetFillStyle(1001);
  hErecoFlav_NueBar->SetFillColor(kRed);
  hErecoFlav_NueBar->SetLineColor(kRed);

  hErecoFlav_Numu->SetFillStyle(1001);
  hErecoFlav_Numu->SetFillColor(kBlue);
  hErecoFlav_Numu->SetLineColor(kBlue);

  hErecoFlav_NumuBar->SetFillStyle(1001);
  hErecoFlav_NumuBar->SetFillColor(kWhite);
  hErecoFlav_NumuBar->SetLineColor(kBlack);
  hErecoFlav_NumuBar->SetLineWidth(2);

  THStack *sErecoTotFlav = new THStack("sErecoTotFlav","sErecoTotFlav");
  sErecoTotFlav->Add(hErecoFlav_NueBar);
  sErecoTotFlav->Add(hErecoFlav_Nue);
  sErecoTotFlav->Add(hErecoFlav_Numu);
  sErecoTotFlav->Add(hErecoFlav_NumuBar);

  TCanvas *cErecoTotFlav = new TCanvas("cErecoTotFlav","cErecoTotFlav");
  sErecoTotFlav->Draw("histo");

 

  //combine modes
  double nccqe = 0, nccnqe = 0, nnuecc = 0, nnc = 0, nmec = 0;
  double nnc1pic = 0, n1pi;
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    double nmode = hErecoMode[im]->Integral();
    if(util::IsNC(m))
      nnc += nmode;
    else if(util::IsMEC(m))
      nmec += nmode;
    else if(util::IsNuE(m) || util::IsNuEBar(m) || util::IsOscNuE(m))
      nnuecc += nmode;
    else if(m == conv::kSKmNuMuCCQE || m == conv::kSKmNuMuBarCCQE
	    || m == conv::kSKmOscNuMuCCQE)// || m == conv::kSKmOscNuMuBarCCQE)
      nccqe += nmode;
    else
      nccnqe += nmode;
    if(m == conv::kSKmNuMuNC1piC
       || m == conv::kSKmNuMuBarNC1piC
       || m == conv::kSKmNuENC1piC
       || m == conv::kSKmNuEBarNC1piC
       )
      nnc1pic += nmode;
    


    if(m == conv::kSKmNuMuNCoth
       || m == conv::kSKmNuMuBarNCoth
       || m == conv::kSKmNuENCoth
       || m == conv::kSKmNuEBarNCoth
       || m == conv::kSKmNuMuNC1piC
       || m == conv::kSKmNuMuBarNC1piC
       || m == conv::kSKmNuENC1piC
       || m == conv::kSKmNuEBarNC1piC
       || m == conv::kSKmNuMuCC1pi
       || m == conv::kSKmNuMuBarCC1pi
       || m == conv::kSKmNuECC1pi
       || m == conv::kSKmNuEBarCC1pi
       || m == conv::kSKmOscNuMuCC1pi
       || m == conv::kSKmOscNuECC1pi
       )
      n1pi += nmode;


  }
  streamer << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "CCQE" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nccqe << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nccqe/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "CCnQE" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nccnqe << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nccnqe/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "NC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnc << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnc/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "$\\nu_e$CC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnuecc/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "MEC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nmec << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nmec/ntot << " %)" << endl;

  streamer << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "Non-MEC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc+nnc+nccnqe+nccqe << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*(nnuecc+nnc+nccnqe+nccqe)/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "Non-MEC or NC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc+nccnqe+nccqe << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*(nnuecc+nccnqe+nccqe)/ntot << " %)" << endl;
  streamer << nnuecc+nnc+nccnqe+nccqe+nmec << endl;

  streamer << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "NC1piC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnc1pic << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnc1pic/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "1pi" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << n1pi << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*n1pi/ntot << " %)" << endl;
     
  
  //
  // Draw
  //
  
  TH1D * hErecoTot_draw        = hErecoTot;
  TH1D * hErecoTot_draw_coarse = hErecoTot_coarse;
  if (plot_per_50MeV) {
    hErecoTot_draw        = spectrum_to_per_50mev(hErecoTot);
    hErecoTot_draw_coarse = spectrum_to_per_50mev(hErecoTot_coarse);
  }

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,1000,1000);
  TLegend * legend = new TLegend(0.6, 0.4, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot_draw->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
  hErecoTot_draw->Draw();
  THStack * stack = new THStack("stack","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack->Add(hErecoMode_draw[im]);
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
  c->SaveAs(Form("sk_ereco_spectrum_%s_%f_%f%s.%s",sk_sample,input_mixing23, input_dm23sq, (disable_banff_tuning ? "_nobanff" : ""),plot_format));
  hErecoTot_draw->GetXaxis()->UnZoom();

  //draw coarse binned histograms
  TCanvas * c_coarse = new TCanvas("c_coarse","",20,20,500,500);
  TLegend * legend_coarse = new TLegend(0.6, 0.4, 0.9, 0.9);
  legend_coarse->SetName("legend_coarse");
  legend_coarse->SetFillColor(0);
  legend_coarse->SetBorderSize(0);
  hErecoTot_draw_coarse->GetXaxis()->SetRangeUser(0,7-conv::kASmallNum);
  hErecoTot_draw_coarse->Draw();
  THStack * stack_coarse = new THStack("stack_coarse","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack_coarse->Add(hErecoMode_draw_coarse[im]);
  }
  stack_coarse->Draw("same");
  legend_coarse->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend_coarse->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend_coarse->Draw();
  c_coarse->Update();
  c_coarse->SaveAs(Form("sk_ereco_spectrum_%s_%f_%f%s_coarse.%s",sk_sample,input_mixing23, input_dm23sq, (disable_banff_tuning ? "_nobanff" : ""),plot_format));

  // get tweaks
  PhysParamList nuisance_params = skexp.Tweaks();

  //
  // Save histograms in a ROOT file
  //

  TFile f(Form("sk_ereco_spectrum_%s_%s-%s%s-%s.root",sk_sample,tag,oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : ""),sk_sample) ,"recreate");
  nuisance_params.Write("nuisance_params");
  hErecoTot -> Write(util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str());
  hErecoTot -> Write("erec_sk"); //KDI
  hErecoTot_coarse->Write(Form("%s_coarse",util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str()));
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    hErecoMode[im] -> Write(util::AsStringCompact(s,m,"sk_","").c_str());
    hErecoMode_coarse[im] -> Write(Form("sk_%s_coarse",util::AsStringCompact(s,m).c_str()));
  }
  stack->Write();
  legend->Write();
  stack_coarse->Write();
  legend_coarse->Write();
  f.Close();




















  //
  // Print-out numbers requested during the analysis review
  //

  if(print_spectrum) {  
    int nprec=5;
    streamer << fixed;

    streamer << endl;
    for(int ibin = 1; ibin <= hErecoTot_coarse->GetXaxis()->GetNbins(); ibin++) {
      streamer << setw(6) << setfill(' ') << ibin
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinLowEdge(ibin)
	       << " - "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinUpEdge(ibin)
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot_coarse->GetBinContent(ibin)
	       << endl;
    }//ibin

    streamer << endl;
    for(int ibin = 1; ibin <= hErecoTot->GetXaxis()->GetNbins(); ibin++) {
      streamer << setw(6) << setfill(' ') << ibin
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinLowEdge(ibin)
	       << " - "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinUpEdge(ibin)
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot->GetBinContent(ibin)
	       << endl;
    }//ibin
  }

  //actually print
  cout << streamer.str() << endl;
  
  // get size of streamer
  long size;
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);

  //write to file
  ofstream outfile (Form("sk_ereco_spectrum_%s_%s-%s%s-%s.txt",sk_sample,tag,oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : ""),sk_sample));
  outfile.write (streamer.str().c_str(),size);
  outfile.close();

}


TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}


