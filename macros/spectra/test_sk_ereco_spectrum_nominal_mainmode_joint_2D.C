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
 
 const char * sk_sample            = "1Re",                         /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
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
 
   
 bool         store_spectrum      = true,
 bool         print_spectrum      = true,                          /* */
 const char * mixparamchoice      = "sinsqtheta",
 const char * hierarchy           = "normal"   ,
 const char * dmsqscheme          = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

 bool plot_per_50MeV = false,
 double ErecoRangeLow = 0.,
 double ErecoRangeUp = 1.5 
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
  const int nbins_ereco = util::NumErecoBins(s);
  const double * ereco_bin_edges = util::ErecoBinEdges(s);
  TH2D * hErecoVsQ2Tot = new TH2D("hErecoVsQ2Tot","hErecoVsQ2Tot",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);

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
  


  //Get energy spectra
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
  
  TCanvas *cErecoTot_Proj = new TCanvas("cErecoTot_Proj","cErecoTot_Proj");
  hErecoVsQ2Tot->ProjectionY()->Draw();

  TCanvas *cQ2Tot_Proj = new TCanvas("cQ2Tot_Proj","cQ2Tot_Proj");
  hErecoVsQ2Tot->ProjectionX()->Draw();

  TCanvas *cErecoVsQ2Tot = new TCanvas("cErecoVsQ2Tot","cErecoVsQ2Tot");
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Tot->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco
  hErecoVsQ2Tot->Draw("colz");


  // Build mode/flavor

  TH2D * hErecoVsQ2Mode_NC = new TH2D("hErecoVsQ2Mode_NC","hErecoVsQ2Mode_NC",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NueBeam = new TH2D("hErecoVsQ2Mode_NueBeam","hErecoVsQ2Mode_NueBeam",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NumuBeam = new TH2D("hErecoVsQ2Mode_NumuBeam","hErecoVsQ2Mode_NumuBeam",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NumuBarBeam = new TH2D("hErecoVsQ2Mode_NumuBarBeam","hErecoVsQ2Mode_NumuBarBeam",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NumuBarToNueBar = new TH2D("hErecoVsQ2Mode_NumuBarToNueBar","hErecoVsQ2Mode_NumuBarToNueBar",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NumuToNue = new TH2D("hErecoVsQ2Mode_NumuToNue","hErecoVsQ2Mode_NumuToNue",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);

  TH1D * hErecoMode_NC;
  TH1D * hErecoMode_NueBeam;
  TH1D * hErecoMode_NumuBeam;
  TH1D * hErecoMode_NumuBarBeam;
  TH1D * hErecoMode_NumuBarToNueBar;
  TH1D * hErecoMode_NumuToNue;

  TH1D * hQ2Mode_NC;
  TH1D * hQ2Mode_NueBeam;
  TH1D * hQ2Mode_NumuBeam;
  TH1D * hQ2Mode_NumuBarBeam;
  TH1D * hQ2Mode_NumuBarToNueBar;
  TH1D * hQ2Mode_NumuToNue;

  hErecoVsQ2Mode_NC->SetFillStyle(0);
  hErecoVsQ2Mode_NueBeam->SetFillStyle(0);
  hErecoVsQ2Mode_NumuBeam->SetFillStyle(0);
  hErecoVsQ2Mode_NumuBarBeam->SetFillStyle(0);
  hErecoVsQ2Mode_NumuBarToNueBar->SetFillStyle(0);
  hErecoVsQ2Mode_NumuToNue->SetFillStyle(0);

  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
    string modestr = util::AsStringCompact(m);
        
    //cout << m << " " << modestr << " --> ";
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
      hErecoVsQ2Mode_NC->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // NC
    else if(util::IsOscNuEBar(m))
      hErecoVsQ2Mode_NumuBarToNueBar->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // NumuBar->NueBar
    else if(util::IsOscNuE(m))
      hErecoVsQ2Mode_NumuToNue->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // Numu->Nue
    else if( util::IsNuMuBar(m))
      hErecoVsQ2Mode_NumuBarBeam->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // NumuBar Beam
    else if( util::IsNuMu(m))
      hErecoVsQ2Mode_NumuBeam->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // Numu Beam
    else if( util::IsNuE(m) || util::IsNuEBar(m) )
      hErecoVsQ2Mode_NueBeam->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // Nue & NueBar Beam
    else{
      cout << modestr << endl;
      cout << "No mode / flavors are selected!!!" << endl;
      exit(1);
    }
  } // loop over modes

  // Get 1D projections X: Q2
  hQ2Mode_NC              = (TH1D*) hErecoVsQ2Mode_NC             ->ProjectionX()->Clone("hQ2Mode_NC");
  hQ2Mode_NumuBarBeam     = (TH1D*) hErecoVsQ2Mode_NumuBarBeam    ->ProjectionX()->Clone("hQ2Mode_NumuBarBeam");
  hQ2Mode_NumuBeam        = (TH1D*) hErecoVsQ2Mode_NumuBeam       ->ProjectionX()->Clone("hQ2Mode_NumuBeam");
  hQ2Mode_NueBeam         = (TH1D*) hErecoVsQ2Mode_NueBeam        ->ProjectionX()->Clone("hQ2Mode_NueBeam");
  hQ2Mode_NumuBarToNueBar = (TH1D*) hErecoVsQ2Mode_NumuBarToNueBar->ProjectionX()->Clone("hQ2Mode_NumuBarToNueBar");
  hQ2Mode_NumuToNue       = (TH1D*) hErecoVsQ2Mode_NumuToNue      ->ProjectionX()->Clone("hQ2Mode_NumuToNue");
  
  // Get 1D projections Y: Ereco
  hErecoMode_NC              = (TH1D*) hErecoVsQ2Mode_NC             ->ProjectionY()->Clone("hErecoMode_NC");
  hErecoMode_NumuBarBeam     = (TH1D*) hErecoVsQ2Mode_NumuBarBeam    ->ProjectionY()->Clone("hErecoMode_NumuBarBeam");
  hErecoMode_NumuBeam        = (TH1D*) hErecoVsQ2Mode_NumuBeam       ->ProjectionY()->Clone("hErecoMode_NumuBeam");
  hErecoMode_NueBeam         = (TH1D*) hErecoVsQ2Mode_NueBeam        ->ProjectionY()->Clone("hErecoMode_NueBeam");
  hErecoMode_NumuBarToNueBar = (TH1D*) hErecoVsQ2Mode_NumuBarToNueBar->ProjectionY()->Clone("hErecoMode_NumuBarToNueBar");
  hErecoMode_NumuToNue       = (TH1D*) hErecoVsQ2Mode_NumuToNue      ->ProjectionY()->Clone("hErecoMode_NumuToNue");



  // Draw 2D templates

  TCanvas *cErecoVsQ2_NC = new TCanvas("cErecoVsQ2_NC","cErecoVsQ2_NC"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_NC->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco 
  hErecoVsQ2Mode_NC->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hErecoVsQ2Mode_NC->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hErecoVsQ2Mode_NC->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NC->Draw("colz");
  TLatex *lErecoVsQ2Mode_NC = new TLatex(1.,0.3,"NC");
  lErecoVsQ2Mode_NC->SetTextSize(0.07);
  lErecoVsQ2Mode_NC->Draw("same");

  TCanvas *cErecoVsQ2_NumuBarBeam = new TCanvas("cErecoVsQ2_NumuBarBeam","cErecoVsQ2_NumuBarBeam"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_NumuBarBeam->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco
  hErecoVsQ2Mode_NumuBarBeam->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hErecoVsQ2Mode_NumuBarBeam->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hErecoVsQ2Mode_NumuBarBeam->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NumuBarBeam->Draw("colz"); 
  TLatex *lErecoVsQ2Mode_NumuBarBeam = new TLatex(1.,0.3,"#bar{#nu}_{#mu} intrinsic");
  lErecoVsQ2Mode_NumuBarBeam->SetTextSize(0.07);
  lErecoVsQ2Mode_NumuBarBeam->Draw("same");

  TCanvas *cErecoVsQ2_NumuBeam = new TCanvas("cErecoVsQ2_NumuBeam","cErecoVsQ2_NumuBeam"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_NumuBeam->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco  
  hErecoVsQ2Mode_NumuBeam->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hErecoVsQ2Mode_NumuBeam->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hErecoVsQ2Mode_NumuBeam->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NumuBeam->Draw("colz");
  TLatex *lErecoVsQ2Mode_NumuBeam = new TLatex(1.,0.3,"#nu_{#mu} intrinsic");
  lErecoVsQ2Mode_NumuBeam->SetTextSize(0.07);
  lErecoVsQ2Mode_NumuBeam->Draw("same");

  TCanvas *cErecoVsQ2_NueBeam = new TCanvas("cErecoVsQ2_NueBeam","cErecoVsQ2_NueBeam"); 
  if(ErecoRangeLow!=ErecoRangeUp)   
    hErecoVsQ2Mode_NueBeam->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco  
  hErecoVsQ2Mode_NueBeam->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hErecoVsQ2Mode_NueBeam->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hErecoVsQ2Mode_NueBeam->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NueBeam->Draw("colz");
  TLatex *lErecoVsQ2Mode_NueBeam = new TLatex(1.,0.3,"#nu_{e}/#bar{#nu}_{e} intrinsic");
  lErecoVsQ2Mode_NueBeam->SetTextSize(0.07);
  lErecoVsQ2Mode_NueBeam->Draw("same");

  TCanvas *cErecoVsQ2_NumuBarToNueBar = new TCanvas("cErecoVsQ2_NumuBarToNueBar","cErecoVsQ2_NumuBarToNueBar"); 
  if(ErecoRangeLow!=ErecoRangeUp)   
    hErecoVsQ2Mode_NumuBarToNueBar->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco  
  hErecoVsQ2Mode_NumuBarToNueBar->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hErecoVsQ2Mode_NumuBarToNueBar->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hErecoVsQ2Mode_NumuBarToNueBar->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NumuBarToNueBar->Draw("colz");
  TLatex *lErecoVsQ2Mode_NumuBarToNueBar = new TLatex(1.,0.3,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}");
  lErecoVsQ2Mode_NumuBarToNueBar->SetTextSize(0.07);
  lErecoVsQ2Mode_NumuBarToNueBar->Draw("same");

  TCanvas *cErecoVsQ2_NumuToNue = new TCanvas("cErecoVsQ2_NumuToNue","cErecoVsQ2_NumuToNue"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_NumuToNue->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco  
  hErecoVsQ2Mode_NumuToNue->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hErecoVsQ2Mode_NumuToNue->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hErecoVsQ2Mode_NumuToNue->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NumuToNue->Draw("colz");
  TLatex *lErecoVsQ2Mode_NumuToNue = new TLatex(1.,0.3,"#nu_{#mu} #rightarrow #nu_{e}");
  lErecoVsQ2Mode_NumuToNue->SetTextSize(0.07);
  lErecoVsQ2Mode_NumuToNue->Draw("same");



  // Draw 1D templates: Ereco

  if(plot_per_50MeV){
    hErecoMode_NC = spectrum_to_per_50mev(hErecoMode_NC); 
    hErecoMode_NumuBarBeam = spectrum_to_per_50mev(hErecoMode_NumuBarBeam); 
    hErecoMode_NumuBeam = spectrum_to_per_50mev(hErecoMode_NumuBeam); 
    hErecoMode_NueBeam = spectrum_to_per_50mev(hErecoMode_NueBeam); 
    hErecoMode_NumuBarToNueBar = spectrum_to_per_50mev(hErecoMode_NumuBarToNueBar); 
    hErecoMode_NumuToNue = spectrum_to_per_50mev(hErecoMode_NumuToNue); 
  }

  hErecoMode_NC->SetFillStyle(1001);
  hErecoMode_NC->SetLineColor(kBlack);
  //hErecoMode_NC->SetLineColor(kRed+2);
  hErecoMode_NC->SetFillColor(kRed);
  hErecoMode_NC->SetLineWidth(2);
  
  hErecoMode_NumuBarBeam->SetFillStyle(1001);
  hErecoMode_NumuBarBeam->SetLineColor(kBlack);
  hErecoMode_NumuBarBeam->SetFillColor(kBlue);
  hErecoMode_NumuBarBeam->SetLineWidth(2);

  hErecoMode_NumuBeam->SetFillStyle(1001);
  hErecoMode_NumuBeam->SetLineColor(kBlack);
  hErecoMode_NumuBeam->SetFillColor(kCyan);
  hErecoMode_NumuBeam->SetLineWidth(2);

  hErecoMode_NueBeam->SetFillStyle(1001);
  hErecoMode_NueBeam->SetLineColor(kBlack);
  hErecoMode_NueBeam->SetFillColor(kYellow);
  hErecoMode_NueBeam->SetLineWidth(2);

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
  hErecoMode_NumuBarBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NumuBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NueBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NumuToNue->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NumuBarToNueBar->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  
  if(plot_per_50MeV){
    hErecoMode_NC->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NumuBarBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NumuBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NueBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NumuToNue->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NumuBarToNueBar->GetYaxis()->SetTitle("Events per 0.050 GeV");
  }
  else{
    hErecoMode_NC->GetYaxis()->SetTitle("Events");
    hErecoMode_NumuBarBeam->GetYaxis()->SetTitle("Events");
    hErecoMode_NumuBeam->GetYaxis()->SetTitle("Events");
    hErecoMode_NueBeam->GetYaxis()->SetTitle("Events");
    hErecoMode_NumuToNue->GetYaxis()->SetTitle("Events");
    hErecoMode_NumuBarToNueBar->GetYaxis()->SetTitle("Events");
  }

  THStack *sErecoTotFlavMode = new THStack("sErecoTotFlavMode","sErecoTotFlavMode");
  sErecoTotFlavMode->Add(hErecoMode_NC);
  sErecoTotFlavMode->Add(hErecoMode_NueBeam);
  sErecoTotFlavMode->Add(hErecoMode_NumuBeam);
  sErecoTotFlavMode->Add(hErecoMode_NumuBarBeam);
  sErecoTotFlavMode->Add(hErecoMode_NumuToNue);
  sErecoTotFlavMode->Add(hErecoMode_NumuBarToNueBar);
  
  TCanvas *cErecoTotFlavMode = new TCanvas("cErecoTotFlavMode","cErecoTotFlavMode");
  sErecoTotFlavMode->Draw();
  sErecoTotFlavMode->GetXaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco    
  sErecoTotFlavMode->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  if(plot_per_50MeV) sErecoTotFlavMode->GetYaxis()->SetTitle("Events per 0.050 GeV");
  else sErecoTotFlavMode->GetYaxis()->SetTitle("Events");
  cErecoTotFlavMode->Modified();
  gStyle->SetOptStat(0);

  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu} mode, e-like Total","L");
  leg->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","F");
  leg->AddEntry(hErecoMode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","F");
  leg->AddEntry(hErecoMode_NC,"NC","F");
  leg->AddEntry(hErecoMode_NueBeam,"#nu_{e}/#bar{#nu}_{e} intrinsic","F");
  leg->AddEntry(hErecoMode_NumuBeam,"#nu_{#mu} intrinsic","F");
  leg->AddEntry(hErecoMode_NumuBarBeam,"#bar{#nu}_{#mu} intrinsic","F");
  leg->Draw();



  // Draw 1D templates: Q2

  if(plot_per_50MeV){
    hQ2Mode_NC = spectrum_to_per_50mev(hQ2Mode_NC); 
    hQ2Mode_NumuBarBeam = spectrum_to_per_50mev(hQ2Mode_NumuBarBeam); 
    hQ2Mode_NumuBeam = spectrum_to_per_50mev(hQ2Mode_NumuBeam); 
    hQ2Mode_NueBeam = spectrum_to_per_50mev(hQ2Mode_NueBeam); 
    hQ2Mode_NumuBarToNueBar = spectrum_to_per_50mev(hQ2Mode_NumuBarToNueBar); 
    hQ2Mode_NumuToNue = spectrum_to_per_50mev(hQ2Mode_NumuToNue); 
  }

  hQ2Mode_NC->SetFillStyle(1001);
  hQ2Mode_NC->SetLineColor(kBlack);
  //hQ2Mode_NC->SetLineColor(kRed+2);
  hQ2Mode_NC->SetFillColor(kRed);
  hQ2Mode_NC->SetLineWidth(2);
  
  hQ2Mode_NumuBarBeam->SetFillStyle(1001);
  hQ2Mode_NumuBarBeam->SetLineColor(kBlack);
  hQ2Mode_NumuBarBeam->SetFillColor(kBlue);
  hQ2Mode_NumuBarBeam->SetLineWidth(2);

  hQ2Mode_NumuBeam->SetFillStyle(1001);
  hQ2Mode_NumuBeam->SetLineColor(kBlack);
  hQ2Mode_NumuBeam->SetFillColor(kCyan);
  hQ2Mode_NumuBeam->SetLineWidth(2);

  hQ2Mode_NueBeam->SetFillStyle(1001);
  hQ2Mode_NueBeam->SetLineColor(kBlack);
  hQ2Mode_NueBeam->SetFillColor(kYellow);
  hQ2Mode_NueBeam->SetLineWidth(2);

  hQ2Mode_NumuBarToNueBar->SetFillStyle(1001);
  hQ2Mode_NumuBarToNueBar->SetLineColor(kBlack);
  hQ2Mode_NumuBarToNueBar->SetFillColor(kGreen);
  hQ2Mode_NumuBarToNueBar->SetLineWidth(2);

  hQ2Mode_NumuToNue->SetFillStyle(3004);
  hQ2Mode_NumuToNue->SetLineColor(kBlack);
  hQ2Mode_NumuToNue->SetFillColor(kGray+3);
  hQ2Mode_NumuToNue->SetLineWidth(2);
  //hQ2Mode_NumuToNue->SetFillStyle(3004);
 
  if(plot_per_50MeV){
    hQ2Mode_NC->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NumuBarBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NumuBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NueBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NumuToNue->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NumuBarToNueBar->GetYaxis()->SetTitle("Events per 0.050 GeV");
  }
  else{
    hQ2Mode_NC->GetXaxis()->SetTitle("Q^{2} (GeV)");
    hQ2Mode_NumuBarBeam->GetXaxis()->SetTitle("Q^{2} (GeV)");
    hQ2Mode_NumuBeam->GetXaxis()->SetTitle("Q^{2} (GeV)");
    hQ2Mode_NueBeam->GetXaxis()->SetTitle("Q^{2} (GeV)");
    hQ2Mode_NumuToNue->GetXaxis()->SetTitle("Q^{2} (GeV)");
    hQ2Mode_NumuBarToNueBar->GetXaxis()->SetTitle("Q^{2} (GeV)");
  }

  THStack *sQ2TotFlavMode = new THStack("sQ2TotFlavMode","sQ2TotFlavMode");
  sQ2TotFlavMode->Add(hQ2Mode_NC);
  sQ2TotFlavMode->Add(hQ2Mode_NueBeam);
  sQ2TotFlavMode->Add(hQ2Mode_NumuBeam);
  sQ2TotFlavMode->Add(hQ2Mode_NumuBarBeam);
  sQ2TotFlavMode->Add(hQ2Mode_NumuToNue);
  sQ2TotFlavMode->Add(hQ2Mode_NumuBarToNueBar);

  TCanvas *cQ2TotFlavMode = new TCanvas("cQ2TotFlavMode","cQ2TotFlavMode");
  sQ2TotFlavMode->Draw();
  sQ2TotFlavMode->GetXaxis()->SetTitle("Q^{2} (GeV)");
  if(plot_per_50MeV) sQ2TotFlavMode->GetYaxis()->SetTitle("Events per 0.050 GeV");
  else sQ2TotFlavMode->GetYaxis()->SetTitle("Events");
  cQ2TotFlavMode->Modified();
  gStyle->SetOptStat(0);

  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hQ2Mode_NumuBarToNueBar,"#bar{#nu} mode, e-like Total","L");
  leg->AddEntry(hQ2Mode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","F");
  leg->AddEntry(hQ2Mode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","F");
  leg->AddEntry(hQ2Mode_NC,"NC","F");
  leg->AddEntry(hQ2Mode_NueBeam,"#nu_{e}/#bar{#nu}_{e} intrinsic","F");
  leg->AddEntry(hQ2Mode_NumuBeam,"#nu_{#mu} intrinsic","F");
  leg->AddEntry(hQ2Mode_NumuBarBeam,"#bar{#nu}_{#mu} intrinsic","F");
  leg->Draw();






  // Get total # of events

  double TotEvts = hErecoVsQ2Tot->Integral();

  double TotEvtsCheck = 
    hErecoVsQ2Mode_NC->Integral() + 
    hErecoVsQ2Mode_NumuBarBeam->Integral() + 
    hErecoVsQ2Mode_NumuBeam->Integral() + 
    hErecoVsQ2Mode_NueBeam->Integral() + 
    hErecoVsQ2Mode_NumuBarToNueBar->Integral() + 
    hErecoVsQ2Mode_NumuToNue->Integral();
  
  if( fabs(TotEvts-TotEvtsCheck)>0.0000001 ){
    cout << endl;
    cout << "TotEvts = " << TotEvts << endl;
    cout << "TotEvtsCheck = " << TotEvtsCheck << endl;
    cout << "There is an error!" << endl;
    cout << endl;
    exit(1);
  }

  cout << endl;
  cout << "Numu --> Nue = " << hErecoVsQ2Mode_NumuToNue->Integral() << endl;
  cout << "NumuBar --> NueBar = " << hErecoVsQ2Mode_NumuBarToNueBar->Integral() << endl; 
  cout << "NC = " << hErecoVsQ2Mode_NC->Integral() << endl; 
  cout << "Numu beam = " << hErecoVsQ2Mode_NumuBeam->Integral() << endl; 
  cout << "NumuBar beam = " << hErecoVsQ2Mode_NumuBarBeam->Integral() << endl; 
  cout << "Nue beam = " << hErecoVsQ2Mode_NueBeam->Integral() << endl; 
  cout << "Total = " << TotEvts << endl;
  cout << endl;




  // Print canvases

  if(print_spectrum){
    
    char name[200];

    sprintf(name,"%s_%s.pdf","cErecoVsQ2Tot",sk_sample);    
    cErecoVsQ2Tot ->Print(name);

    sprintf(name,"%s_%s.pdf","cErecoVsQ2_NC",sk_sample);    
    cErecoVsQ2_NC ->Print(name);
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_NumuBarBeam",sk_sample);    
    cErecoVsQ2_NumuBarBeam ->Print(name);
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_NumuBeam",sk_sample);    
    cErecoVsQ2_NumuBeam ->Print(name);
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_NueBeam",sk_sample);    
    cErecoVsQ2_NueBeam ->Print(name);
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_NumuBarToNueBar",sk_sample);    
    cErecoVsQ2_NumuBarToNueBar ->Print(name);
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_NumuToNue",sk_sample);    
    cErecoVsQ2_NumuToNue ->Print(name);

    sprintf(name,"%s_%s.pdf","cQ2TotFlavMode",sk_sample);    
    cQ2TotFlavMode ->Print(name);

    sprintf(name,"%s_%s.pdf","cErecoTotFlavMode",sk_sample);    
    cErecoTotFlavMode ->Print(name);
  }


  return;
}


TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}


