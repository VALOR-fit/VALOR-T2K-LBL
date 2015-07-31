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

void test_sk_ereco_spectrum_nominal_mainmode_joint_2D
(
 const char * tag                   = "prova",
 
 const char * sk_sample            = "1Rmu",                         /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
 const char * hypothesis           = "3f_osc",                 /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_osc_normbar, 3f_oscbar_normbar" */
 
 // Mach3, T2K Run1-4 + reactors best fit (A)
 double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp        = -1.601, //-99,    /* deltacp, */
 double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 double       input_betabar        = -99,
 
 // // Mach3, T2K Run1-4 + reactors best fit - dcp=0 (B)
 // double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 0., //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar        = -99,

 // // Mach3, T2K only Run1-4 best fit (C)
 // double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0422,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 1.91, //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar        = -99,

 // // Mach3, T2K Run1-4 + reactors best fit - th23 non maximal (D)
 // double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.4, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = -1.601, //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar        = -99,

 bool         disable_banff_tuning = false,                          /* */
 
 const char * config_file      = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml",
 //const char * config_file      = "analysis_configuration_run1to4_joint_Nu6.57E20POT.xml",
 
   
 bool         store_spectrum     = true,
 bool         print_spectrum       = false,                          /* */
 const char * mixparamchoice      = "sinsqtheta",
 const char * hierarchy           = "normal"   ,
 const char * dmsqscheme          = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

  bool plot_per_50MeV = false
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
  const int nbins_Q2 = util::NumQ2Bins(s);
  const double * Q2_bin_edges = util::Q2BinEdges(s);
  const int nbins = util::NumErecoBins(s);
  const double * ereco_bin_edges = util::ErecoBinEdges(s);
  TH1D * hErecoTot = new TH1D("hErecoTot","hErecoTot",nbins,ereco_bin_edges); 
  TH2D * hErecoVsQ2Tot = new TH2D("hErecoVsQ2Tot","hErecoVsQ2Tot",nbins,ereco_bin_edges,nbins_Q2,Q2_bin_edges);

  TH1D * hErecoMode[conv::knm]; //< contribution from each mode  
  TH1D * hErecoMode_draw[conv::knm]; //< contribution from each mode  

  TH2D * hErecoVsQ2Mode[conv::knm]; //< contribution from each mode  
  TH2D * hErecoVsQ2Mode_draw[conv::knm]; //< contribution from each mode  

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
    TH2D * h = skexp.CurrPredictedErecoVsQ2Spectrum(s, m);  
    assert(h);
    hErecoVsQ2Mode[im] = new TH2D(*h);
    if (plot_per_50MeV) {
      hErecoVsQ2Mode_draw[im]        = spectrum_to_per_50mev(hErecoVsQ2Mode[im]);
    }
    else {
      hErecoVsQ2Mode_draw[im]        = hErecoVsQ2Mode[im];
    }      
  }
  
  hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));
  
  hErecoVsQ2Tot->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hErecoVsQ2Tot->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hErecoVsQ2Tot->GetZaxis()->SetTitle("Number of events");
  
  hErecoVsQ2Tot->Draw("colz");

  cout << "# Tot events = " << hErecoVsQ2Tot->Integral() << endl;


  cout << "CIAO" << endl;

  return;



  // //can save this to file &/or print it
  // stringstream streamer;

  // streamer << fixed;
  // double ntot = hErecoTot->Integral();
  // double ntot_coarse = hErecoTot->Integral();
  // if(TMath::Abs(ntot - ntot_coarse) > conv::kASmallNum) {
  //   cerr << "Integrals are different between differently binned histograms "
  // 	 << "50MeV: " << ntot
  // 	 << "BANFF binning: " << ntot_coarse << endl;
  //   exit(1);
  // }

  // streamer << "Total number of events: " 
  //      << setw(10) << setfill(' ') << setprecision(6) << ntot << endl;

  // const int nslices=6;
  // double slice_max[nslices] = { 0.50, 0.75, 1.00, 1.25, 1.50, 2.00 };
  // for(int islice=0; islice<nslices; islice++) {
  //   double esmax = slice_max[islice];
  //   streamer << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << esmax << " GeV): " 
  //        << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(esmax)) << endl;
  // }
  // for(int im = 0; im < conv::knm; im++) {
  //    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
  //    if(m == conv::kSKmUndefined) continue;
  //    double nmode = hErecoMode[im]->Integral();
  //    streamer << "Contrib from mode: " 
  //         << setw(40) << setfill(' ') << util::AsString(m) << " : " 
  //         << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
  //         << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  // }


  // // Build mode/flavor

  // TH1D * hErecoMode_NC = new TH1D("hErecoMode_NC","hErecoMode_NC",nbins,ereco_bin_edges);
  // TH1D * hErecoMode_NueBeam = new TH1D("hErecoMode_NueBeam","hErecoMode_NueBeam",nbins,ereco_bin_edges);
  // TH1D * hErecoMode_NumuBeam = new TH1D("hErecoMode_NumuBeam","hErecoMode_NumuBeam",nbins,ereco_bin_edges);
  // TH1D * hErecoMode_NumuBarBeam = new TH1D("hErecoMode_NumuBarBeam","hErecoMode_NumuBarBeam",nbins,ereco_bin_edges);
  // TH1D * hErecoMode_NumuBarToNueBar = new TH1D("hErecoMode_NumuBarToNueBar","hErecoMode_NumuBarToNueBar",nbins,ereco_bin_edges);
  // TH1D * hErecoMode_NumuToNue = new TH1D("hErecoMode_NumuToNue","hErecoMode_NumuToNue",nbins,ereco_bin_edges);

  // // flavor only
  // TH1D * hErecoFlav_Numu = new TH1D("hErecoFlav_Numu","hErecoFlav_Numu",nbins,ereco_bin_edges);
  // TH1D * hErecoFlav_NumuBar = new TH1D("hErecoFlav_NumuBar","hErecoFlav_NumuBar",nbins,ereco_bin_edges);
  // TH1D * hErecoFlav_Nue = new TH1D("hErecoFlav_Nue","hErecoFlav_Nue",nbins,ereco_bin_edges);
  // TH1D * hErecoFlav_NueBar = new TH1D("hErecoFlav_NueBar","hErecoFlav_NueBar",nbins,ereco_bin_edges);

  // hErecoFlav_Numu->SetFillStyle(0);
  // hErecoFlav_NumuBar->SetFillStyle(0);
  // hErecoFlav_Nue->SetFillStyle(0);
  // hErecoFlav_NueBar->SetFillStyle(0);
  
  // hErecoMode_NC->SetFillStyle(0);
  // hErecoMode_NueBeam->SetFillStyle(0);
  // hErecoMode_NumuBeam->SetFillStyle(0);
  // hErecoMode_NumuBarBeam->SetFillStyle(0);
  // hErecoMode_NumuBarToNueBar->SetFillStyle(0);
  // hErecoMode_NumuToNue->SetFillStyle(0);

  // for (int im = 0; im < conv::knm; im++) {
  //   conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
  //   if(m == conv::kSKmUndefined) continue;
  //   if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
  //   string modestr = util::AsStringCompact(m);
    
  //   // Fill by flavor
  //   if(util::IsNuE(m) || util::IsOscNuE(m))
  //     hErecoFlav_Nue->Add(skexp.CurrPredictedErecoSpectrum(s,m));
  //   else if(util::IsNuEBar(m) || util::IsOscNuEBar(m))
  //     hErecoFlav_NueBar->Add(skexp.CurrPredictedErecoSpectrum(s,m));
  //   else if(util::IsNuMuBar(m) || util::IsOscNuMuBar(m))
  //     hErecoFlav_NumuBar->Add(skexp.CurrPredictedErecoSpectrum(s,m));
  //   else if(util::IsNuMu(m) || util::IsOscNuMu(m))
  //     hErecoFlav_Numu->Add(skexp.CurrPredictedErecoSpectrum(s,m));
  //   else{
  //     cout << modestr << endl;
  //     cout << "No flavors are selected!!!" << endl;
  //     exit(1);
  //   }
    
  //   cout << m << " " << modestr << " --> ";
  //   if(util::IsCCQE(m)) cout << "IsCCQE" << endl;
  //   else if(util::IsMEC(m)) cout << "IsMEC" << endl;
  //   else if(util::IsCCnQE(m)) cout << "IsCCnQE" << endl;
  //   else if(util::IsNC(m)) cout << "IsNC" << endl;
  //   else{
  //     cout << "No modes!!!" << endl;
  //     exit(1);
  //   }
  //   cout << endl;
    
  //   // Fill by flavor/mode
    
  //   if(util::IsNC(m))    
  //     hErecoMode_NC->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // NC
  //   else if(util::IsOscNuEBar(m))
  //     hErecoMode_NumuBarToNueBar->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // NumuBar->NueBar
  //   else if(util::IsOscNuE(m))
  //     hErecoMode_NumuToNue->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // Numu->Nue
  //   else if( util::IsNuMuBar(m))
  //     hErecoMode_NumuBarBeam->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // NumuBar Beam
  //   else if( util::IsNuMu(m))
  //     hErecoMode_NumuBeam->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // Numu Beam
  //   else if( util::IsNuE(m) || util::IsNuEBar(m) )
  //     hErecoMode_NueBeam->Add(skexp.CurrPredictedErecoSpectrum(s,m)); // Nue & NueBar Beam
  //   else{
  //     cout << modestr << endl;
  //     cout << "No mode / flavors are selected!!!" << endl;
  //     exit(1);
  //   }

  // } // loop over modes


  
  // // Draw Histo by mode/flavor

  // //gStyle->SetHatchesSpacing();
  // //gStyle->SetHatchesLineWidth(3.);
  
  // if(plot_per_50MeV){
  //   hErecoMode_NC = spectrum_to_per_50mev(hErecoMode_NC); 
  //   hErecoMode_NumuBarBeam = spectrum_to_per_50mev(hErecoMode_NumuBarBeam); 
  //   hErecoMode_NumuBeam = spectrum_to_per_50mev(hErecoMode_NumuBeam); 
  //   hErecoMode_NueBeam = spectrum_to_per_50mev(hErecoMode_NueBeam); 
  //   hErecoMode_NumuBarToNueBar = spectrum_to_per_50mev(hErecoMode_NumuBarToNueBar); 
  //   hErecoMode_NumuToNue = spectrum_to_per_50mev(hErecoMode_NumuToNue); 
  // }

  // hErecoMode_NC->SetFillStyle(1001);
  // hErecoMode_NC->SetLineColor(kBlack);
  // //hErecoMode_NC->SetLineColor(kRed+2);
  // hErecoMode_NC->SetFillColor(kRed);
  // hErecoMode_NC->SetLineWidth(2);
  
  // hErecoMode_NumuBarBeam->SetFillStyle(1001);
  // hErecoMode_NumuBarBeam->SetLineColor(kBlack);
  // hErecoMode_NumuBarBeam->SetFillColor(kBlue);
  // hErecoMode_NumuBarBeam->SetLineWidth(2);

  // hErecoMode_NumuBeam->SetFillStyle(1001);
  // hErecoMode_NumuBeam->SetLineColor(kBlack);
  // hErecoMode_NumuBeam->SetFillColor(kCyan);
  // hErecoMode_NumuBeam->SetLineWidth(2);

  // hErecoMode_NueBeam->SetFillStyle(1001);
  // hErecoMode_NueBeam->SetLineColor(kBlack);
  // hErecoMode_NueBeam->SetFillColor(kYellow);
  // hErecoMode_NueBeam->SetLineWidth(2);

  // hErecoMode_NumuBarToNueBar->SetFillStyle(1001);
  // hErecoMode_NumuBarToNueBar->SetLineColor(kBlack);
  // hErecoMode_NumuBarToNueBar->SetFillColor(kGreen);
  // hErecoMode_NumuBarToNueBar->SetLineWidth(2);

  // hErecoMode_NumuToNue->SetFillStyle(3004);
  // hErecoMode_NumuToNue->SetLineColor(kBlack);
  // hErecoMode_NumuToNue->SetFillColor(kGray+3);
  // hErecoMode_NumuToNue->SetLineWidth(2);
  // //hErecoMode_NumuToNue->SetFillStyle(3004);

  // hErecoMode_NC->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  // hErecoMode_NC->GetYaxis()->SetTitle("Events per 0.050 GeV");
  // hErecoMode_NumuBarBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  // hErecoMode_NumuBarBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
  // hErecoMode_NumuBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  // hErecoMode_NumuBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
  // hErecoMode_NueBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  // hErecoMode_NueBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
  // hErecoMode_NumuToNue->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  // hErecoMode_NumuToNue->GetYaxis()->SetTitle("Events per 0.050 GeV");
  // hErecoMode_NumuBarToNueBar->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  // hErecoMode_NumuBarToNueBar->GetYaxis()->SetTitle("Events per 0.050 GeV");
  
  // THStack *sErecoTotFlavMode = new THStack("sErecoTotFlavMode","sErecoTotFlavMode");
  // sErecoTotFlavMode->Add(hErecoMode_NC);
  // sErecoTotFlavMode->Add(hErecoMode_NueBeam);
  // sErecoTotFlavMode->Add(hErecoMode_NumuBeam);
  // sErecoTotFlavMode->Add(hErecoMode_NumuBarBeam);
  // sErecoTotFlavMode->Add(hErecoMode_NumuToNue);
  // sErecoTotFlavMode->Add(hErecoMode_NumuBarToNueBar);

  // TCanvas *cErecoTotFlavMode = new TCanvas("cErecoTotFlavMode","cErecoTotFlavMode");
  // sErecoTotFlavMode->Draw();
  // sErecoTotFlavMode->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  // sErecoTotFlavMode->GetYaxis()->SetTitle("Events per 0.050 GeV");
  // cErecoTotFlavMode->Modified();
  // gStyle->SetOptStat(0);

  // TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.9);
  // leg->SetFillColor(0);
  // leg->SetBorderSize(0);
  // leg->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu} mode, e-like Total","L");
  // leg->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","F");
  // leg->AddEntry(hErecoMode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","F");
  // leg->AddEntry(hErecoMode_NC,"NC","F");
  // leg->AddEntry(hErecoMode_NueBeam,"#nu_{e}/#bar{#nu}_{e} intrinsic","F");
  // leg->AddEntry(hErecoMode_NumuBeam,"#nu_{#mu} intrinsic","F");
  // leg->AddEntry(hErecoMode_NumuBarBeam,"#bar{#nu}_{#mu} intrinsic","F");
  // leg->Draw();

  // cout << "Tot = " << hErecoTot->Integral() << endl;
  // cout << "Tot 2 = " 
  //      << 
  //   hErecoMode_NumuBarToNueBar->Integral() + 
  //   hErecoMode_NC->Integral() + 
  //   hErecoMode_NumuBeam->Integral() + 
  //   hErecoMode_NumuBarBeam->Integral() + 
  //   hErecoMode_NueBeam->Integral() + 
  //   hErecoMode_NumuToNue->Integral()
  //      << endl; 
  // cout << "hErecoMode_NumuBarToNueBar = " << hErecoMode_NumuBarToNueBar->Integral() << endl; 
  // cout << "hErecoMode_NC = " << hErecoMode_NC->Integral() << endl; 
  // cout << "hErecoMode_NumuBeam = " << hErecoMode_NumuBeam->Integral() << endl; 
  // cout << "hErecoMode_NumuBarBeam = " << hErecoMode_NumuBarBeam->Integral() << endl; 
  // cout << "hErecoMode_NueBeam = " << hErecoMode_NueBeam->Integral() << endl; 
  // cout << "hErecoMode_NumuToNue = " << hErecoMode_NumuToNue->Integral() << endl;
  


  // TCanvas *cErecoTotFlavMode_norm = new TCanvas("cErecoTotFlavMode_norm","cErecoTotFlavMode_norm");
  // hErecoMode_NC_norm = (TH1D*) hErecoMode_NC->Clone("hErecoMode_NC_norm");
  // hErecoMode_NC_norm->Scale(1./hErecoMode_NC->Integral());
  // hErecoMode_NC_norm->SetFillColor(0);
  // hErecoMode_NC_norm->SetLineColor(kRed);

  // hErecoMode_NueBeam_norm = (TH1D*) hErecoMode_NueBeam->Clone("hErecoMode_NueBeam_norm");
  // hErecoMode_NueBeam_norm->Scale(1./hErecoMode_NueBeam->Integral());
  // hErecoMode_NueBeam_norm->SetFillColor(0);
  // hErecoMode_NueBeam_norm->SetLineColor(kYellow);

  // hErecoMode_NumuBeam_norm = (TH1D*) hErecoMode_NumuBeam->Clone("hErecoMode_NumuBeam_norm");
  // hErecoMode_NumuBeam_norm->Scale(1./hErecoMode_NumuBeam->Integral());
  // hErecoMode_NumuBeam_norm->SetFillColor(0);
  // hErecoMode_NumuBeam_norm->SetLineColor(kCyan);

  // hErecoMode_NumuBarBeam_norm = (TH1D*) hErecoMode_NumuBarBeam->Clone("hErecoMode_NumuBarBeam_norm");
  // hErecoMode_NumuBarBeam_norm->Scale(1./hErecoMode_NumuBarBeam->Integral());
  // hErecoMode_NumuBarBeam_norm->SetFillColor(0);
  // hErecoMode_NumuBarBeam_norm->SetLineColor(kBlue);

  // hErecoMode_NumuToNue_norm = (TH1D*) hErecoMode_NumuToNue->Clone("hErecoMode_NumuToNue_norm");
  // hErecoMode_NumuToNue_norm->Scale(1./hErecoMode_NumuToNue_norm->Integral());
  // hErecoMode_NumuToNue_norm->SetFillColor(0);
  // hErecoMode_NumuToNue_norm->SetLineColor(kGray+3);

  // hErecoMode_NumuBarToNueBar_norm = (TH1D*) hErecoMode_NumuBarToNueBar->Clone("hErecoMode_NumuBarToNueBar_norm");
  // hErecoMode_NumuBarToNueBar_norm->Scale(1./hErecoMode_NumuBarToNueBar->Integral());
  // hErecoMode_NumuBarToNueBar_norm->SetFillColor(0);
  // hErecoMode_NumuBarToNueBar_norm->SetLineColor(kGreen);
  
  // hErecoMode_NumuBarToNueBar_norm->DrawCopy("");
  // hErecoMode_NueBeam_norm->DrawCopy("same");
  // hErecoMode_NumuBeam_norm->DrawCopy("same");
  // hErecoMode_NumuBarBeam_norm->DrawCopy("same");
  // hErecoMode_NC_norm->DrawCopy("same");
  // hErecoMode_NumuToNue_norm->DrawCopy("same");
  
  // TLegend *leg_norm = new TLegend(0.6, 0.4, 0.9, 0.9);
  // leg_norm->SetFillColor(0);
  // leg_norm->SetBorderSize(0);
  // leg_norm->AddEntry(hErecoMode_NumuBarToNueBar_norm,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","L");
  // leg_norm->AddEntry(hErecoMode_NumuToNue_norm,"#nu_{#mu} #rightarrow #nu_{e}","L");
  // leg_norm->AddEntry(hErecoMode_NC_norm,"NC","L");
  // leg_norm->AddEntry(hErecoMode_NueBeam_norm,"#nu_{e}/#bar{#nu}_{e} intrinsic","L");
  // leg_norm->AddEntry(hErecoMode_NumuBeam_norm,"#nu_{#mu} intrinsic","L");
  // leg_norm->AddEntry(hErecoMode_NumuBarBeam_norm,"#bar{#nu}_{#mu} intrinsic","L");

  // leg_norm->Draw();
  
  
  // // Draw Histo by flavor

  // hErecoFlav_Nue->SetFillStyle(1001);
  // hErecoFlav_Nue->SetFillColor(kYellow);
  // hErecoFlav_Nue->SetLineColor(kYellow);

  // hErecoFlav_NueBar->SetFillStyle(1001);
  // hErecoFlav_NueBar->SetFillColor(kRed);
  // hErecoFlav_NueBar->SetLineColor(kRed);

  // hErecoFlav_Numu->SetFillStyle(1001);
  // hErecoFlav_Numu->SetFillColor(kBlue);
  // hErecoFlav_Numu->SetLineColor(kBlue);

  // hErecoFlav_NumuBar->SetFillStyle(1001);
  // hErecoFlav_NumuBar->SetFillColor(kWhite);
  // hErecoFlav_NumuBar->SetLineColor(kBlack);
  // hErecoFlav_NumuBar->SetLineWidth(2);

  // THStack *sErecoTotFlav = new THStack("sErecoTotFlav","sErecoTotFlav");
  // sErecoTotFlav->Add(hErecoFlav_NueBar);
  // sErecoTotFlav->Add(hErecoFlav_Nue);
  // sErecoTotFlav->Add(hErecoFlav_Numu);
  // sErecoTotFlav->Add(hErecoFlav_NumuBar);

  // TCanvas *cErecoTotFlav = new TCanvas("cErecoTotFlav","cErecoTotFlav");
  // sErecoTotFlav->Draw("histo");


  // //combine modes
  // double nccqe = 0, nccnqe = 0, nnuecc = 0, nnc = 0, nmec = 0;
  // double nnc1pic = 0, n1pi;
  // for(int im = 0; im < conv::knm; im++) {
  //   conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
  //   if(m == conv::kSKmUndefined) continue;
  //   double nmode = hErecoMode[im]->Integral();
  //   if(util::IsNC(m))
  //     nnc += nmode;
  //   else if(util::IsMEC(m))
  //     nmec += nmode;
  //   else if(util::IsNuE(m) || util::IsNuEBar(m) || util::IsOscNuE(m))
  //     nnuecc += nmode;
  //   else if(m == conv::kSKmNuMuCCQE || m == conv::kSKmNuMuBarCCQE
  // 	    || m == conv::kSKmOscNuMuCCQE)// || m == conv::kSKmOscNuMuBarCCQE)
  //     nccqe += nmode;
  //   else
  //     nccnqe += nmode;
  //   if(m == conv::kSKmNuMuNC1piC
  //      || m == conv::kSKmNuMuBarNC1piC
  //      || m == conv::kSKmNuENC1piC
  //      || m == conv::kSKmNuEBarNC1piC
  //      )
  //     nnc1pic += nmode;
    


  //   if(m == conv::kSKmNuMuNCoth
  //      || m == conv::kSKmNuMuBarNCoth
  //      || m == conv::kSKmNuENCoth
  //      || m == conv::kSKmNuEBarNCoth
  //      || m == conv::kSKmNuMuNC1piC
  //      || m == conv::kSKmNuMuBarNC1piC
  //      || m == conv::kSKmNuENC1piC
  //      || m == conv::kSKmNuEBarNC1piC
  //      || m == conv::kSKmNuMuCC1pi
  //      || m == conv::kSKmNuMuBarCC1pi
  //      || m == conv::kSKmNuECC1pi
  //      || m == conv::kSKmNuEBarCC1pi
  //      || m == conv::kSKmOscNuMuCC1pi
  //      || m == conv::kSKmOscNuECC1pi
  //      )
  //     n1pi += nmode;


  // }
  // streamer << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "CCQE" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << nccqe << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*nccqe/ntot << " %)" << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "CCnQE" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << nccnqe << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*nccnqe/ntot << " %)" << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "NC" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << nnc << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnc/ntot << " %)" << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "$\\nu_e$CC" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnuecc/ntot << " %)" << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "MEC" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << nmec << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*nmec/ntot << " %)" << endl;

  // streamer << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "Non-MEC" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc+nnc+nccnqe+nccqe << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*(nnuecc+nnc+nccnqe+nccqe)/ntot << " %)" << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "Non-MEC or NC" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc+nccnqe+nccqe << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*(nnuecc+nccnqe+nccqe)/ntot << " %)" << endl;
  // streamer << nnuecc+nnc+nccnqe+nccqe+nmec << endl;

  // streamer << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "NC1piC" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << nnc1pic << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnc1pic/ntot << " %)" << endl;
  // streamer << "Contrib from modes: " 
  // 	   << setw(40) << setfill(' ') << "1pi" << " : " 
  // 	   << setw(15) << setfill(' ') << setprecision(6) << n1pi << " events ("
  // 	   << setw(10) << setfill(' ') << setprecision(6) << 100*n1pi/ntot << " %)" << endl;
     
  // cout << endl;
  // cout << "Numu --> Nue = " << hErecoMode_NumuToNue->Integral() << endl;
  // cout << "NumuBar --> NueBar = " << hErecoMode_NumuBarToNueBar->Integral() << endl; 
  // cout << "NC = " << hErecoMode_NC->Integral() << endl; 
  // cout << "Numu beam = " << hErecoMode_NumuBeam->Integral() << endl; 
  // cout << "NumuBar beam = " << hErecoMode_NumuBarBeam->Integral() << endl; 
  // cout << "Nue beam = " << hErecoMode_NueBeam->Integral() << endl; 
  // cout << "Total = " << hErecoTot->Integral() << endl;
  // cout << endl;
  
  // //
  // // Draw
  // //
  
  // TH1D * hErecoTot_draw        = hErecoTot;
  // TH1D * hErecoTot_draw_coarse = hErecoTot_coarse;
  // if (plot_per_50MeV) {
  //   hErecoTot_draw        = spectrum_to_per_50mev(hErecoTot);
  //   hErecoTot_draw_coarse = spectrum_to_per_50mev(hErecoTot_coarse);
  // }

  // gStyle->SetOptStat(0);
  // TCanvas * c = new TCanvas("c","",20,20,1000,1000);
  // TLegend * legend = new TLegend(0.6, 0.4, 0.9, 0.9);
  // legend->SetFillColor(0);
  // legend->SetBorderSize(0);
  // hErecoTot_draw->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
  // hErecoTot_draw->Draw();
  // THStack * stack = new THStack("stack","spectrum contributions");
  // // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // // so as to enhance visibility in the stack of smaller components
  // for(int imstack = 0; imstack < conv::knm; imstack++) {
  //   conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
  //   if(m == conv::kSKmUndefined) continue;
  //   int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  //   stack->Add(hErecoMode_draw[im]);
  // }
  // stack->Draw("same");
  // legend->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  // for(int im = 0; im < conv::knm; im++) {
  //   conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
  //   if(m == conv::kSKmUndefined) continue;
  //   legend->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  // }
  // legend->Draw();
  // c->Update();
  // c->SaveAs(Form("sk_ereco_spectrum_%s_%f_%f%s.%s",sk_sample,input_mixing23, input_dm23sq, (disable_banff_tuning ? "_nobanff" : ""),plot_format));
  // hErecoTot_draw->GetXaxis()->UnZoom();

 
  // //draw coarse binned histograms
  // TCanvas * c_coarse = new TCanvas("c_coarse","",20,20,500,500);
  // TLegend * legend_coarse = new TLegend(0.6, 0.4, 0.9, 0.9);
  // legend_coarse->SetName("legend_coarse");
  // legend_coarse->SetFillColor(0);
  // legend_coarse->SetBorderSize(0);
  // hErecoTot_draw_coarse->GetXaxis()->SetRangeUser(0,7-conv::kASmallNum);
  // hErecoTot_draw_coarse->Draw();
  // THStack * stack_coarse = new THStack("stack_coarse","spectrum contributions");
  // // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // // so as to enhance visibility in the stack of smaller components
  // for(int imstack = 0; imstack < conv::knm; imstack++) {
  //   conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
  //   if(m == conv::kSKmUndefined) continue;
  //   int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  //   stack_coarse->Add(hErecoMode_draw_coarse[im]);
  // }
  // stack_coarse->Draw("same");
  // legend_coarse->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  // for(int im = 0; im < conv::knm; im++) {
  //   conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
  //   if(m == conv::kSKmUndefined) continue;
  //   legend_coarse->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  // }
  // legend_coarse->Draw();
  // c_coarse->Update();
  // c_coarse->SaveAs(Form("sk_ereco_spectrum_%s_%f_%f%s_coarse.%s",sk_sample,input_mixing23, input_dm23sq, (disable_banff_tuning ? "_nobanff" : ""),plot_format));

  // // get tweaks
  // PhysParamList nuisance_params = skexp.Tweaks();

  // //
  // // Save histograms in a ROOT file
  // //

  // TFile f(Form("sk_ereco_spectrum_%s_%s-%s%s-%s.root",sk_sample,tag,oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : ""),sk_sample) ,"recreate");
  // nuisance_params.Write("nuisance_params");
  // hErecoTot -> Write(util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str());
  // hErecoTot -> Write("erec_sk"); //KDI
  // hErecoTot_coarse->Write(Form("%s_coarse",util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str()));
  // for(int im = 0; im < conv::knm; im++) {
  //   conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
  //   if(m == conv::kSKmUndefined) continue;
  //   hErecoMode[im] -> Write(util::AsStringCompact(s,m,"sk_","").c_str());
  //   hErecoMode_coarse[im] -> Write(Form("sk_%s_coarse",util::AsStringCompact(s,m).c_str()));
  // }
  // stack->Write();
  // legend->Write();
  // stack_coarse->Write();
  // legend_coarse->Write();

  // hErecoMode_NC->Write();
  // hErecoMode_NumuBeam->Write();
  // hErecoMode_NumuBarBeam->Write();
  // hErecoMode_NueBeam->Write();
  // hErecoMode_NumuToNue->Write();
  // hErecoMode_NumuBarToNueBar->Write();

  // hErecoMode_NC_norm->Write();
  // hErecoMode_NumuBeam_norm->Write();
  // hErecoMode_NumuBarBeam_norm->Write();
  // hErecoMode_NueBeam_norm->Write();
  // hErecoMode_NumuToNue_norm->Write();
  // hErecoMode_NumuBarToNueBar_norm->Write();

  // f.Close();




  // //
  // // Print-out numbers requested during the analysis review
  // //

  // if(print_spectrum) {  
  //   int nprec=5;
  //   streamer << fixed;

  //   streamer << endl;
  //   for(int ibin = 1; ibin <= hErecoTot_coarse->GetXaxis()->GetNbins(); ibin++) {
  //     streamer << setw(6) << setfill(' ') << ibin
  // 	       << " : "
  // 	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinLowEdge(ibin)
  // 	       << " - "
  // 	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinUpEdge(ibin)
  // 	       << " : "
  // 	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot_coarse->GetBinContent(ibin)
  // 	       << endl;
  //   }//ibin

  //   streamer << endl;
  //   for(int ibin = 1; ibin <= hErecoTot->GetXaxis()->GetNbins(); ibin++) {
  //     streamer << setw(6) << setfill(' ') << ibin
  // 	       << " : "
  // 	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinLowEdge(ibin)
  // 	       << " - "
  // 	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinUpEdge(ibin)
  // 	       << " : "
  // 	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot->GetBinContent(ibin)
  // 	       << endl;
  //   }//ibin
  // }

  // //actually print
  // cout << streamer.str() << endl;
  
  // // get size of streamer
  // long size;
  // streamer.seekg(0,ifstream::end);
  // size=streamer.tellg();
  // streamer.seekg(0);

  // //write to file
  // ofstream outfile (Form("sk_ereco_spectrum_%s_%s-%s%s-%s.txt",sk_sample,tag,oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : ""),sk_sample));
  // outfile.write (streamer.str().c_str(),size);
  // outfile.close();

}


TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}


