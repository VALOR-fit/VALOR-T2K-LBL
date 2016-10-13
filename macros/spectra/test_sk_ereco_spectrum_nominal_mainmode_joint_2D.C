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
 const char * tag                  = "tamesu",
 const char * sk_sample            = "MultiRe",                         /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu", "MultiRe" */
 const char * hypothesis           = "3f_osc",                 /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_osc_normbar, 3f_oscbar_normbar" */ 
 
 // // No Oscillations
 // double       input_mixing12       = 0.,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0., //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 0., //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar        = -99,

 // Mach3, T2K Run1-4 + reactors best fit (1)
 double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13       = 0.0217,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp        = -1.601, //-99,    /* deltacp, */
 double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 double       input_betabar        = 1., //0., //-99,

 // // Mach3, T2K Run1-4 + reactors best fit - dcp=0 + first octant (2)
 // double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0217,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.45, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 0., //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar        = 1., //0., //-99,

 // Mach3, T2K Run1-4 + reactors best fit (1)
 // double       input_mixing12       = 0.306,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0257,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.5, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 0., //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.6E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.49E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar        = -99,
 
 bool         disable_banff_tuning = true,                          /* */
 
 //
 // NEW CONFIGURATION FILES (4/3/16 BANFF RELEASE)
 //
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h.xml",
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_ExpSummer2016.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_FSIOnly.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_SKDetOnly.xml",
 
 //
 // OLD CONFIGURATION FILES
 //
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_NIWG2015c.xml",
 //const char * config_file      = "analysis_configuration_test_xsecterm_elike_RFGRPA_NIWG2015c.xml",
 //const char * config_file      = "analysis_configuration_test_q2_elike_RFGRPA_NIWG2015c.xml",
 //const char * config_file    = "analysis_configuration_test_theta_elike_RFGRPA_NIWG2015c_mrmc_test.xml",
 // const char * config_file    = "analysis_configuration.xml",
 const char * config_file     = "analysis_configuration_Run1-7cPOT_theta.xml",    //"analysis_configuration.xml",
 //const char * config_file    = "analysis_configuration_dummy_multiring.xml",
 //const char * config_file  = "analysis_configuration_test_theta_elike_RFGRPA.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_NIWG2015c_Nu6.57E20POT.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_NIWG2015c.xml",
 //const char * config_file      = "analysis_configuration_test_prefit_theta_1bin_NIWG2015c.xml",
 //const char * config_file      = "analysis_configuration_test_prefit_theta_elike_RFGRPA_NIWG2015c.xml",
 //const char * config_file      = "analysis_configuration_test_theta_1bin_NIWG2015c.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_NIWG2015c_Asimov.xml",
 //const char * config_file      = "analysis_configuration_test_theta_1bin_NIWG2015c_Asimov.xml",
 //const char * config_file      = "analysis_configuration_test_theta_1bin_NIWG2015c_SKDetonly_Asimov.xml",
 //const char * config_file      = "analysis_configuration_test_theta_1bin_NIWG2015c_FSISIonly_Asimov.xml",
 
 bool         store_spectrum     = true,
 bool         print_spectrum     = true,                          /* */
 const char * mixparamchoice     = "sinsqtheta",
 const char * hierarchy          = "normal"   ,
 const char * dmsqscheme         = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamoctant     = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

 bool plot_per_50MeV             = false,
 double ErecoRangeLow            = 0.,
 double ErecoRangeUp             = 30.,

 int          iexperiment         = 0, //4404, //4739,                       /* */
 //const char * rootfile = "data_run1to4/SKdata_Run14.root",
 // const char * rootfile = "data_run1to6/SKdata_Run16.root",
 const char * rootfile = "", 
 const bool draw_data_spectra = false
)
{
  using namespace std;
  using namespace valor;

  char name[200];
  
  const char * plot_format = "eps";

  double postlatex_X = 120.; // Theta
  double postlatex_Y = 0.05;
  double postlatex_reac_X = 120.; 
  double postlatex_reac_Y = 0.05;
  char titleX[200] = "#theta (degrees)";
  TLegend *legX = new TLegend(0.6, 0.4, 0.9, 0.9); 
  TLegend *legX_reac = new TLegend(0.6, 0.4, 0.9, 0.9); 

  // double postlatex_X = 1.0; // Q^2
  // double postlatex_Y = 0.3;
  // double postlatex_reac_X = 1.0; 
  // double postlatex_reac_Y = 0.3;
  // char titleX[200] = "Q^{2} (GeV^{2})";
  // TLegend *legX = new TLegend(0.6, 0.4, 0.9, 0.9); 
  // TLegend *legX_reac = new TLegend(0.6, 0.4, 0.9, 0.9); 

  // double postlatex_X = 0.3; // Qd^2
  // double postlatex_Y = 0.3;
  // char titleX[200] = "Q_{d}^{2} (GeV^{2})";
  // TLegend *legX = new TLegend(0.2, 0.4, 0.5, 0.9); 

  // double postlatex_X = 0.6; // energy hadronic system
  // double postlatex_Y = 0.3;
  // char titleX[200] = "#nu";
  // TLegend *legX = new TLegend(0.6, 0.4, 0.9, 0.9); 

  // double postlatex_X = 0.6; // inelasticity
  // double postlatex_Y = 0.3;
  // char titleX[200] = "y";
  // TLegend *legX = new TLegend(0.6, 0.4, 0.9, 0.9); 
 
  // double postlatex_X = 0.6; // +- term in xsec formula
  // double postlatex_Y = 0.3;
  // double postlatex_reac_X = 0.6; 
  // double postlatex_reac_Y = 0.3;
  // char titleX[200] = "XSec Term";
  // TLegend *legX = new TLegend(0.6, 0.4, 0.9, 0.9); 
  // TLegend *legX_reac = new TLegend(0.6, 0.4, 0.9, 0.9); 
  
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
    ph.SetCurr
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
				    input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq,input_betabar,-999,-999);

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  string sconfig(config);

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << sconfig << endl;
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

  TH1D * hErecoMode[conv::knm];      //< contribution from each mode  
  TH1D * hErecoMode_draw[conv::knm]; //< contribution from each mode  

  TH2D * hErecoVsQ2Modes[conv::knm];      //< contribution from each mode  
  TH2D * hErecoVsQ2Mode_draw[conv::knm]; //< contribution from each mode  

 
  // Data files

  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));
  SuperkObservation * skobs_fakedata = new SuperkObservation;
  cout << "Use root file" << endl;
  bool draw_data = true;
  TFile *infile = new TFile(rootfile2.Data(),"READ");
  if(!infile->IsOpen()){
    cout << endl;
    cout << "The data root file doesn't exist: don't draw data spectra" << endl;
    cout << endl;
    draw_data = false;
  }
  if(!draw_data_spectra){
    draw_data = false;
  }
  infile->Close();
  
  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(sconfig);

  // Disable BANFF tuning if requested
  if(disable_banff_tuning) {
    cout << "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..." << endl;
    //skexp.GetBANFF()->RestoreNominalMC();
    skexp.GetBANFF()->Tweak(0,0,"untuned");
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
    hErecoVsQ2Modes[im] = new TH2D(*h);
    if (plot_per_50MeV) {
      hErecoVsQ2Mode_draw[im] = spectrum_to_per_50mev(hErecoVsQ2Modes[im]);
    }
    else {
      hErecoVsQ2Mode_draw[im] = hErecoVsQ2Modes[im];
    }      
  }
  
  
  hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
  hErecoVsQ2Tot->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Tot->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Tot->GetZaxis()->SetTitle("Number of events");
  
  TH1D *hErecoTot = (TH1D*) hErecoVsQ2Tot->ProjectionY()->Clone("hErecoTot");
  TH1D *hQ2Tot = (TH1D*) hErecoVsQ2Tot->ProjectionX()->Clone("hQ2Tot");

  // Get data spectra and rescale MC to data

  TH2D *h_data;
  TH1D *h_data_ereco;
  double scale_ereco = 0.;
  TH1D *h_data_q2;
  double scale_q2 = 0.;
  TLatex *lScaleEreco;
  TLatex *lScaleQ2;
  if(draw_data){

    bool using_xml_file = false;
    skobs_fakedata->Read(rootfile2.Data(),iexperiment,using_xml_file);
    h_data = (TH2D*) skobs_fakedata->TwoDSpectrum(s)->Clone("h_data");
    h_data->SetMarkerSize(2);

    h_data_ereco = (TH1D*) h_data->ProjectionY()->Clone("h_data_ereco");
    scale_ereco = 
      (int) ( 2. * ( h_data_ereco->GetMaximum() / hErecoVsQ2Tot->ProjectionY()->GetMaximum() + 0.5 ) );
    
    h_data_q2 = (TH1D*) h_data->ProjectionX()->Clone("h_data_q2");
    scale_q2 = 
      (int) ( 2. * ( h_data_q2->GetMaximum() / hErecoVsQ2Tot->ProjectionX()->GetMaximum() + 0.5 ) );
    
    hErecoTot->Scale(scale_ereco);
    hQ2Tot->Scale(scale_q2);
    
    sprintf(name,"MCx%d",scale_ereco);
    double scale_ereco_x = hErecoTot->GetBinCenter(5);
    double scale_ereco_y = hErecoTot->GetMaximum() * 0.85;
    lScaleEreco = new TLatex(scale_ereco_x,scale_ereco_y,name);
    lScaleEreco->SetTextSize(0.07);
    
    sprintf(name,"MCx%d",scale_q2);
    double scale_q2_x = hQ2Tot->GetBinCenter(5);
    double scale_q2_y = hQ2Tot->GetMaximum() * 0.85;
    lScaleQ2 = new TLatex(scale_q2_x,scale_q2_y,name);
    lScaleQ2->SetTextSize(0.07);
  }
  
  //

  TCanvas *cErecoTot_Proj = new TCanvas("cErecoTot_Proj","cErecoTot_Proj");
  hErecoTot->Draw();
  if(draw_data){
    h_data_ereco->Draw("pe same");
    lScaleEreco->Draw("same");
  }
  
  TCanvas *cQ2Tot_Proj = new TCanvas("cQ2Tot_Proj","cQ2Tot_Proj");
  hQ2Tot->Draw();
  if(draw_data){
    h_data_q2->Draw("pe same");
    lScaleQ2->Draw("same");
  }

  TCanvas *cErecoVsQ2Tot = new TCanvas("cErecoVsQ2Tot","cErecoVsQ2Tot");
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Tot->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco
  hErecoVsQ2Tot->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }
  
  
  // Build mode/flavor

  TH2D * hErecoVsQ2Mode_NueBeam = new TH2D("hErecoVsQ2Mode_NueBeam","hErecoVsQ2Mode_NueBeam",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NumuBeam = new TH2D("hErecoVsQ2Mode_NumuBeam","hErecoVsQ2Mode_NumuBeam",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NumuBarBeam = new TH2D("hErecoVsQ2Mode_NumuBarBeam","hErecoVsQ2Mode_NumuBarBeam",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NumuBarToNueBar = new TH2D("hErecoVsQ2Mode_NumuBarToNueBar","hErecoVsQ2Mode_NumuBarToNueBar",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NumuToNue = new TH2D("hErecoVsQ2Mode_NumuToNue","hErecoVsQ2Mode_NumuToNue",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_NC = new TH2D("hErecoVsQ2Mode_NC","hErecoVsQ2Mode_NC",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_MEC = new TH2D("hErecoVsQ2Mode_MEC","hErecoVsQ2Mode_MEC",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_CCQE = new TH2D("hErecoVsQ2Mode_CCQE","hErecoVsQ2Mode_CCQE",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);
  TH2D * hErecoVsQ2Mode_CCnQE = new TH2D("hErecoVsQ2Mode_CCnQE","hErecoVsQ2Mode_CCnQE",nbins_Q2,Q2_bin_edges,nbins_ereco,ereco_bin_edges);

  TH1D * hErecoMode_NueBeam;
  TH1D * hErecoMode_NumuBeam;
  TH1D * hErecoMode_NumuBarBeam;
  TH1D * hErecoMode_NumuBarToNueBar;
  TH1D * hErecoMode_NumuToNue;
  TH1D * hErecoMode_NC;
  TH1D * hErecoMode_MEC;
  TH1D * hErecoMode_CCQE;
  TH1D * hErecoMode_CCnQE;

  TH1D * hQ2Mode_NueBeam;
  TH1D * hQ2Mode_NumuBeam;
  TH1D * hQ2Mode_NumuBarBeam;
  TH1D * hQ2Mode_NumuBarToNueBar;
  TH1D * hQ2Mode_NumuToNue;
  TH1D * hQ2Mode_NC;
  TH1D * hQ2Mode_MEC;
  TH1D * hQ2Mode_CCQE;
  TH1D * hQ2Mode_CCnQE;

  hErecoVsQ2Mode_NueBeam->SetFillStyle(0);
  hErecoVsQ2Mode_NumuBeam->SetFillStyle(0);
  hErecoVsQ2Mode_NumuBarBeam->SetFillStyle(0);
  hErecoVsQ2Mode_NumuBarToNueBar->SetFillStyle(0);
  hErecoVsQ2Mode_NumuToNue->SetFillStyle(0);
  hErecoVsQ2Mode_NC->SetFillStyle(0);
  hErecoVsQ2Mode_MEC->SetFillStyle(0);
  hErecoVsQ2Mode_CCQE->SetFillStyle(0);
  hErecoVsQ2Mode_CCnQE->SetFillStyle(0);

  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
    string modestr = util::AsStringCompact(m);
        
    //cout << m << " " << modestr << " --> ";
    if     (util::IsCCQE(m))  cout << "IsCCQE"  << endl;
    else if(util::IsMEC(m))   cout << "IsMEC"   << endl;
    else if(util::IsCCnQE(m)) cout << "IsCCnQE" << endl;
    else if(util::IsNC(m))    cout << "IsNC"    << endl;
    else{
      cout << "No modes!!!" << endl;
      exit(1);
    }
    
    // Fill by flavor/mode    
    if(util::IsNC(m))    
      //hErecoVsQ2Mode_NC->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // NC
      hErecoVsQ2Mode_NC->Add(hErecoVsQ2Modes[im]); // NC
    else if( util::IsOscNuEBar(m) )
      //hErecoVsQ2Mode_NumuBarToNueBar->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // NumuBar->NueBar
      hErecoVsQ2Mode_NumuBarToNueBar->Add(hErecoVsQ2Modes[im]); // NumuBar->NueBar
    else if( util::IsOscNuE(m) )
      //hErecoVsQ2Mode_NumuToNue->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // Numu->Nue
      hErecoVsQ2Mode_NumuToNue->Add(hErecoVsQ2Modes[im]); // Numu->Nue
    
    else if( util::IsNuMuBar(m) )
      //hErecoVsQ2Mode_NumuBarBeam->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // NumuBar Beam
      hErecoVsQ2Mode_NumuBarBeam->Add(hErecoVsQ2Modes[im]); // NumuBar Beam
    else if( util::IsNuMu(m))
      //hErecoVsQ2Mode_NumuBeam->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // Numu Beam
      hErecoVsQ2Mode_NumuBeam->Add(hErecoVsQ2Modes[im]); // Numu Beam
    else if( util::IsNuE(m) || util::IsNuEBar(m) )
      //hErecoVsQ2Mode_NueBeam->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // Nue & NueBar Beam
      hErecoVsQ2Mode_NueBeam->Add(hErecoVsQ2Modes[im]); // Nue & NueBar Beam
    else{
      cout << modestr << endl;
      cout << "No mode / flavors are selected!!!" << endl;
      exit(1);
    }

    if( util::IsCCQE(m) )    
      //hErecoVsQ2Mode_CCQE->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // CCQE
      hErecoVsQ2Mode_CCQE->Add(hErecoVsQ2Modes[im]); // CCQE
    else if( util::IsMEC(m) )    
      //hErecoVsQ2Mode_MEC->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // MEC {2p-2h}
      hErecoVsQ2Mode_MEC->Add(hErecoVsQ2Modes[im]); // MEC {2p-2h}
    else if( util::IsCCnQE(m) )    
      //hErecoVsQ2Mode_CCnQE->Add(skexp.CurrPredictedErecoVsQ2Spectrum(s,m)); // CCnQE
      hErecoVsQ2Mode_CCnQE->Add(hErecoVsQ2Modes[im]); // CCnQE
    else if( util::IsNC(m) ){} // already filled
    else{
      cout << modestr << endl;
      cout << "No mode / flavors are selected!!!" << endl;
      exit(1);
    }

  } // loop over modes
  
  // Get 1D projections X: Q2
  hQ2Mode_NumuBarBeam     = (TH1D*) hErecoVsQ2Mode_NumuBarBeam    ->ProjectionX()->Clone("hQ2Mode_NumuBarBeam");
  hQ2Mode_NumuBeam        = (TH1D*) hErecoVsQ2Mode_NumuBeam       ->ProjectionX()->Clone("hQ2Mode_NumuBeam");
  hQ2Mode_NueBeam         = (TH1D*) hErecoVsQ2Mode_NueBeam        ->ProjectionX()->Clone("hQ2Mode_NueBeam");
  hQ2Mode_NumuBarToNueBar = (TH1D*) hErecoVsQ2Mode_NumuBarToNueBar->ProjectionX()->Clone("hQ2Mode_NumuBarToNueBar");
  hQ2Mode_NumuToNue       = (TH1D*) hErecoVsQ2Mode_NumuToNue      ->ProjectionX()->Clone("hQ2Mode_NumuToNue");
  hQ2Mode_NC              = (TH1D*) hErecoVsQ2Mode_NC             ->ProjectionX()->Clone("hQ2Mode_NC");
  hQ2Mode_MEC             = (TH1D*) hErecoVsQ2Mode_MEC            ->ProjectionX()->Clone("hQ2Mode_MEC");
  hQ2Mode_CCQE            = (TH1D*) hErecoVsQ2Mode_CCQE           ->ProjectionX()->Clone("hQ2Mode_CCQE");
  hQ2Mode_CCnQE           = (TH1D*) hErecoVsQ2Mode_CCnQE          ->ProjectionX()->Clone("hQ2Mode_CCnQE");
  
  // Get 1D projections Y: Ereco
  hErecoMode_NumuBarBeam     = (TH1D*) hErecoVsQ2Mode_NumuBarBeam    ->ProjectionY()->Clone("hErecoMode_NumuBarBeam");
  hErecoMode_NumuBeam        = (TH1D*) hErecoVsQ2Mode_NumuBeam       ->ProjectionY()->Clone("hErecoMode_NumuBeam");
  hErecoMode_NueBeam         = (TH1D*) hErecoVsQ2Mode_NueBeam        ->ProjectionY()->Clone("hErecoMode_NueBeam");
  hErecoMode_NumuBarToNueBar = (TH1D*) hErecoVsQ2Mode_NumuBarToNueBar->ProjectionY()->Clone("hErecoMode_NumuBarToNueBar");
  hErecoMode_NumuToNue       = (TH1D*) hErecoVsQ2Mode_NumuToNue      ->ProjectionY()->Clone("hErecoMode_NumuToNue");
  hErecoMode_NC              = (TH1D*) hErecoVsQ2Mode_NC             ->ProjectionY()->Clone("hErecoMode_NC");
  hErecoMode_MEC             = (TH1D*) hErecoVsQ2Mode_MEC            ->ProjectionY()->Clone("hErecoMode_MEC");
  hErecoMode_CCQE            = (TH1D*) hErecoVsQ2Mode_CCQE           ->ProjectionY()->Clone("hErecoMode_CCQE");
  hErecoMode_CCnQE           = (TH1D*) hErecoVsQ2Mode_CCnQE          ->ProjectionY()->Clone("hErecoMode_CCnQE");
  

  // Draw 2D templates

  TCanvas *cErecoVsQ2_NumuBarBeam = new TCanvas("cErecoVsQ2_NumuBarBeam","cErecoVsQ2_NumuBarBeam"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_NumuBarBeam->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco
  hErecoVsQ2Mode_NumuBarBeam->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_NumuBarBeam->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_NumuBarBeam->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NumuBarBeam->Draw("colz"); 
  if(draw_data){
    h_data->Draw("same");
  }
  TLatex *lErecoVsQ2Mode_NumuBarBeam = new TLatex(postlatex_X,postlatex_Y,"#bar{#nu}_{#mu} intrinsic");
  lErecoVsQ2Mode_NumuBarBeam->SetTextSize(0.07);
  lErecoVsQ2Mode_NumuBarBeam->Draw("same");

  TCanvas *cErecoVsQ2_NumuBeam = new TCanvas("cErecoVsQ2_NumuBeam","cErecoVsQ2_NumuBeam"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_NumuBeam->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco  
  hErecoVsQ2Mode_NumuBeam->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_NumuBeam->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_NumuBeam->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NumuBeam->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }  
  TLatex *lErecoVsQ2Mode_NumuBeam = new TLatex(postlatex_X,postlatex_Y,"#nu_{#mu} intrinsic");
  lErecoVsQ2Mode_NumuBeam->SetTextSize(0.07);
  lErecoVsQ2Mode_NumuBeam->Draw("same");

  TCanvas *cErecoVsQ2_NueBeam = new TCanvas("cErecoVsQ2_NueBeam","cErecoVsQ2_NueBeam"); 
  if(ErecoRangeLow!=ErecoRangeUp)   
    hErecoVsQ2Mode_NueBeam->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco  
  hErecoVsQ2Mode_NueBeam->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_NueBeam->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_NueBeam->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NueBeam->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }
  TLatex *lErecoVsQ2Mode_NueBeam = new TLatex(postlatex_X,postlatex_Y,"#nu_{e}/#bar{#nu}_{e} intrinsic");
  lErecoVsQ2Mode_NueBeam->SetTextSize(0.07);
  lErecoVsQ2Mode_NueBeam->Draw("same");

  TCanvas *cErecoVsQ2_NumuBarToNueBar = new TCanvas("cErecoVsQ2_NumuBarToNueBar","cErecoVsQ2_NumuBarToNueBar"); 
  if(ErecoRangeLow!=ErecoRangeUp)   
    hErecoVsQ2Mode_NumuBarToNueBar->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco  
  hErecoVsQ2Mode_NumuBarToNueBar->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_NumuBarToNueBar->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_NumuBarToNueBar->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NumuBarToNueBar->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }
  TLatex *lErecoVsQ2Mode_NumuBarToNueBar = new TLatex(postlatex_X,postlatex_Y,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}");
  lErecoVsQ2Mode_NumuBarToNueBar->SetTextSize(0.07);
  lErecoVsQ2Mode_NumuBarToNueBar->Draw("same");

  TCanvas *cErecoVsQ2_NumuToNue = new TCanvas("cErecoVsQ2_NumuToNue","cErecoVsQ2_NumuToNue"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_NumuToNue->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco  
  hErecoVsQ2Mode_NumuToNue->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_NumuToNue->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_NumuToNue->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NumuToNue->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }
  TLatex *lErecoVsQ2Mode_NumuToNue = new TLatex(postlatex_X,postlatex_Y,"#nu_{#mu} #rightarrow #nu_{e}");
  lErecoVsQ2Mode_NumuToNue->SetTextSize(0.07);
  lErecoVsQ2Mode_NumuToNue->Draw("same");

  TCanvas *cErecoVsQ2_NC = new TCanvas("cErecoVsQ2_NC","cErecoVsQ2_NC"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_NC->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco 
  hErecoVsQ2Mode_NC->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_NC->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_NC->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_NC->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }
  TLatex *lErecoVsQ2Mode_NC = new TLatex(postlatex_X,postlatex_Y,"NC");
  lErecoVsQ2Mode_NC->SetTextSize(0.07);
  lErecoVsQ2Mode_NC->Draw("same");

  TCanvas *cErecoVsQ2_MEC = new TCanvas("cErecoVsQ2_MEC","cErecoVsQ2_MEC"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_MEC->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco 
  hErecoVsQ2Mode_MEC->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_MEC->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_MEC->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_MEC->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }
  TLatex *lErecoVsQ2Mode_MEC = new TLatex(postlatex_reac_X,postlatex_reac_Y,"MEC");
  lErecoVsQ2Mode_MEC->SetTextSize(0.07);
  lErecoVsQ2Mode_MEC->Draw("same");

  TCanvas *cErecoVsQ2_CCQE = new TCanvas("cErecoVsQ2_CCQE","cErecoVsQ2_CCQE"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_CCQE->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco 
  hErecoVsQ2Mode_CCQE->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_CCQE->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_CCQE->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_CCQE->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }
  TLatex *lErecoVsQ2Mode_CCQE = new TLatex(postlatex_reac_X,postlatex_reac_Y,"CCQE");
  lErecoVsQ2Mode_CCQE->SetTextSize(0.07);
  lErecoVsQ2Mode_CCQE->Draw("same");

  TCanvas *cErecoVsQ2_CCnQE = new TCanvas("cErecoVsQ2_CCnQE","cErecoVsQ2_CCnQE"); 
  if(ErecoRangeLow!=ErecoRangeUp) 
    hErecoVsQ2Mode_CCnQE->GetYaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco 
  hErecoVsQ2Mode_CCnQE->GetXaxis()->SetTitle(titleX);
  hErecoVsQ2Mode_CCnQE->GetYaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoVsQ2Mode_CCnQE->GetZaxis()->SetTitle("Number of events");
  hErecoVsQ2Mode_CCnQE->Draw("colz");
  if(draw_data){
    h_data->Draw("same");
  }  
  TLatex *lErecoVsQ2Mode_CCnQE = new TLatex(postlatex_reac_X,postlatex_reac_Y,"CCnQE");
  lErecoVsQ2Mode_CCnQE->SetTextSize(0.07);
  lErecoVsQ2Mode_CCnQE->Draw("same");



  // Draw 1D templates: Ereco

  if(plot_per_50MeV){
    hErecoMode_NumuBarBeam = spectrum_to_per_50mev(hErecoMode_NumuBarBeam); 
    hErecoMode_NumuBeam = spectrum_to_per_50mev(hErecoMode_NumuBeam); 
    hErecoMode_NueBeam = spectrum_to_per_50mev(hErecoMode_NueBeam); 
    hErecoMode_NumuBarToNueBar = spectrum_to_per_50mev(hErecoMode_NumuBarToNueBar); 
    hErecoMode_NumuToNue = spectrum_to_per_50mev(hErecoMode_NumuToNue); 
    hErecoMode_NC = spectrum_to_per_50mev(hErecoMode_NC); 
    hErecoMode_MEC = spectrum_to_per_50mev(hErecoMode_MEC); 
    hErecoMode_CCQE = spectrum_to_per_50mev(hErecoMode_CCQE); 
    hErecoMode_CCnQE = spectrum_to_per_50mev(hErecoMode_CCnQE); 
  }

  hErecoMode_NumuBarBeam->SetFillStyle(1001);
  hErecoMode_NumuBarBeam->SetLineColor(kBlack);
  hErecoMode_NumuBarBeam->SetFillColor(kBlue);
  hErecoMode_NumuBarBeam->SetLineWidth(2);
  if(draw_data) hErecoMode_NumuBarBeam->Scale(scale_ereco);
 
  hErecoMode_NumuBeam->SetFillStyle(1001);
  hErecoMode_NumuBeam->SetLineColor(kBlack);
  hErecoMode_NumuBeam->SetFillColor(kCyan);
  hErecoMode_NumuBeam->SetLineWidth(2);
  if(draw_data) hErecoMode_NumuBeam->Scale(scale_ereco);

  hErecoMode_NueBeam->SetFillStyle(1001);
  hErecoMode_NueBeam->SetLineColor(kBlack);
  hErecoMode_NueBeam->SetFillColor(kYellow);
  hErecoMode_NueBeam->SetLineWidth(2);
  if(draw_data) hErecoMode_NueBeam->Scale(scale_ereco);

  hErecoMode_NumuBarToNueBar->SetFillStyle(1001);
  hErecoMode_NumuBarToNueBar->SetLineColor(kBlack);
  hErecoMode_NumuBarToNueBar->SetFillColor(kGreen);
  hErecoMode_NumuBarToNueBar->SetLineWidth(2);
  if(draw_data) hErecoMode_NumuBarToNueBar->Scale(scale_ereco);

  hErecoMode_NumuToNue->SetFillStyle(3004);
  hErecoMode_NumuToNue->SetLineColor(kBlack);
  hErecoMode_NumuToNue->SetFillColor(kGray+3);
  hErecoMode_NumuToNue->SetLineWidth(2);
  if(draw_data) hErecoMode_NumuToNue->Scale(scale_ereco);

  hErecoMode_NC->SetFillStyle(1001);
  hErecoMode_NC->SetLineColor(kBlack);
  hErecoMode_NC->SetFillColor(kRed);
  hErecoMode_NC->SetLineWidth(2);
  if(draw_data) hErecoMode_NC->Scale(scale_ereco);

  hErecoMode_MEC->SetFillStyle(1001);
  hErecoMode_MEC->SetLineColor(kBlack);
  hErecoMode_MEC->SetFillColor(kBlue);
  hErecoMode_MEC->SetLineWidth(2);
  if(draw_data) hErecoMode_MEC->Scale(scale_ereco);
  
  hErecoMode_CCQE->SetFillStyle(1001);
  hErecoMode_CCQE->SetLineColor(kBlack);
  hErecoMode_CCQE->SetFillColor(kGreen);
  hErecoMode_CCQE->SetLineWidth(2);
  if(draw_data) hErecoMode_CCQE->Scale(scale_ereco);

  hErecoMode_CCnQE->SetFillStyle(1001);
  hErecoMode_CCnQE->SetLineColor(kBlack);
  hErecoMode_CCnQE->SetFillColor(kYellow);
  hErecoMode_CCnQE->SetLineWidth(2);
  if(draw_data) hErecoMode_CCnQE->Scale(scale_ereco);
  
  

  hErecoMode_NumuBarBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NumuBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NueBeam->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NumuToNue->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NumuBarToNueBar->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_NC->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_MEC->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_CCQE->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  hErecoMode_CCnQE->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  
  if(plot_per_50MeV){
    hErecoMode_NumuBarBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NumuBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NueBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NumuToNue->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NumuBarToNueBar->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_NC->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_MEC->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_CCQE->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hErecoMode_CCnQE->GetYaxis()->SetTitle("Events per 0.050 GeV");
  }
  else{
    hErecoMode_NumuBarBeam->GetYaxis()->SetTitle("Number of events");
    hErecoMode_NumuBeam->GetYaxis()->SetTitle("Number of events");
    hErecoMode_NueBeam->GetYaxis()->SetTitle("Number of events");
    hErecoMode_NumuToNue->GetYaxis()->SetTitle("Number of events");
    hErecoMode_NumuBarToNueBar->GetYaxis()->SetTitle("Number of events");
    hErecoMode_NC->GetYaxis()->SetTitle("Number of events");
    hErecoMode_MEC->GetYaxis()->SetTitle("Number of events");
    hErecoMode_CCQE->GetYaxis()->SetTitle("Number of events");
    hErecoMode_CCnQE->GetYaxis()->SetTitle("Number of events");
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
  if(draw_data){
    h_data_ereco->Draw("pe same");
    lScaleEreco->Draw("same");
  }
  sErecoTotFlavMode->GetXaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco    
  sErecoTotFlavMode->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  if(plot_per_50MeV) sErecoTotFlavMode->GetYaxis()->SetTitle("Events per 0.050 GeV");
  else sErecoTotFlavMode->GetYaxis()->SetTitle("Number of events");
  cErecoTotFlavMode->Modified();
  gStyle->SetOptStat(0);

  TLegend *legY = new TLegend(0.6, 0.4, 0.9, 0.9);
  //legY->SetFillColor(0);
  legY->SetFillStyle(0);
  legY->SetBorderSize(0);
  legY->AddEntry(hErecoMode_NumuBarToNueBar,"Total","L");
  legY->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","F");
  legY->AddEntry(hErecoMode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","F");
  legY->AddEntry(hErecoMode_NC,"NC","F");
  legY->AddEntry(hErecoMode_NueBeam,"#nu_{e}/#bar{#nu}_{e} intrinsic","F");
  legY->AddEntry(hErecoMode_NumuBeam,"#nu_{#mu} intrinsic","F");
  legY->AddEntry(hErecoMode_NumuBarBeam,"#bar{#nu}_{#mu} intrinsic","F");
  legY->Draw();


  THStack *sErecoTotReac = new THStack("sErecoTotReac","sErecoTotReac");
  sErecoTotReac->Add(hErecoMode_NC);
  sErecoTotReac->Add(hErecoMode_MEC);
  sErecoTotReac->Add(hErecoMode_CCQE);
  sErecoTotReac->Add(hErecoMode_CCnQE);
  
  TCanvas *cErecoTotReac = new TCanvas("cErecoTotReac","cErecoTotReac");
  sErecoTotReac->Draw();
  if(draw_data){
    h_data_ereco->Draw("pe same");
    lScaleEreco->Draw("same");
  }
  sErecoTotReac->GetXaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco    
  sErecoTotReac->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  if(plot_per_50MeV) sErecoTotReac->GetYaxis()->SetTitle("Events per 0.050 GeV");
  else sErecoTotReac->GetYaxis()->SetTitle("Number of events");
  cErecoTotReac->Modified();
  gStyle->SetOptStat(0);

  TLegend *legY_reac = new TLegend(0.6, 0.4, 0.9, 0.9);
  //legY_reac->SetFillColor(0);
  legY_reac->SetFillStyle(0);
  legY_reac->SetBorderSize(0);
  legY_reac->AddEntry(hErecoMode_NC,"NC","F");
  legY_reac->AddEntry(hErecoMode_MEC,"MEC","F");
  legY_reac->AddEntry(hErecoMode_CCQE,"CCQE","F");
  legY_reac->AddEntry(hErecoMode_CCnQE,"CCnQE","F");
  legY_reac->Draw();


  //Rebinning histograms to compare with Konosuke's presentation
  //Could also get the additional modes and stack them in the same 
  THStack *sErecoTotReac_rebin = new THStack("sErecoTotReac_rebin","sErecoTotReac_rebin");
  double xbins[11] = {0,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3};
  TH1D * hErecoMode_NC_rebin    = hErecoMode_NC   ->Rebin(10,"hErecoMode_NC_rebin"   , xbins);
  TH1D * hErecoMode_MEC_rebin   = hErecoMode_MEC  ->Rebin(10,"hErecoMode_MEC_rebin"  , xbins);
  TH1D * hErecoMode_CCQE_rebin  = hErecoMode_CCQE ->Rebin(10,"hErecoMode_CCQE_rebin" , xbins);
  TH1D * hErecoMode_CCnQE_rebin = hErecoMode_CCnQE->Rebin(10,"hErecoMode_CCnQE_rebin", xbins);
  
  sErecoTotReac_rebin->Add(hErecoMode_NC_rebin);
  sErecoTotReac_rebin->Add(hErecoMode_MEC_rebin);
  sErecoTotReac_rebin->Add(hErecoMode_CCQE_rebin);
  sErecoTotReac_rebin->Add(hErecoMode_CCnQE_rebin);

  TCanvas *cErecoTotReac_rebin = new TCanvas("cErecoTotReac_rebin","cErecoTotReac_rebin");
  sErecoTotReac_rebin->Draw();

  sErecoTotReac_rebin->GetXaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco    
  sErecoTotReac_rebin->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  if(plot_per_50MeV) sErecoTotReac_rebin->GetYaxis()->SetTitle("Events per 0.050 GeV");
  else sErecoTotReac_rebin->GetYaxis()->SetTitle("Number of events");
  cErecoTotReac_rebin->Modified();
  gStyle->SetOptStat(0);
  
  TLatex latex_reac_osc;

  // cout << "input_mixing12 " << input_mixing12 << "\n";
  // cout << "input_mixing13 " << input_mixing13 << "\n";
  // cout << "input_mixing23 " << input_mixing23 << "\n";
  // cout << "input_deltacp  " << input_deltacp  << "\n";
  // cout << "input_dm12sq   " << input_dm12sq   << "\n";
  // cout << "input_dm23sq   " << input_dm23sq   << "\n";
  
  // const char * input_mixing12_str = Form("%s%f","sin^{2} #theta_{12} = ",input_mixing12 );
  // const char * input_mixing13_str = Form("%s%f","sin^{2} #theta_{13} = ",input_mixing13 );
  // const char * input_mixing23_str = Form("%s%f","sin^{2} #theta_{23} = ",input_mixing23 );  
  // const char * input_deltacp_str  = Form("%s%f","#delta_{CP} = ",          input_deltacp);
  // const char * input_dm12sq_str   = Form("%s%f","#Delta m^{2}_{12} = ",    input_dm12sq ); 
  // const char * input_dm23sq_str   = Form("%s%f","#Delta m^{2}_{23} = ",    input_dm23sq ); 

  // latex_reac_osc.SetTextAlign(12);
  // latex_reac_osc.SetTextSize(0.05);
  // double oset = 0.05;
  // double pos_y   = 0.65;
  // latex_reac_osc.DrawLatex(0.1,pos_y,input_mixing12_str);
  // latex_reac_osc.DrawLatex(0.1,pos_y + 1*oset,input_mixing13_str);
  // latex_reac_osc.DrawLatex(0.1,pos_y + 2*oset,input_mixing23_str);
  // latex_reac_osc.DrawLatex(0.1,pos_y + 3*oset,input_deltacp_str );
  // latex_reac_osc.DrawLatex(0.1,pos_y + 4*oset,input_dm12sq_str  );  
  // latex_reac_osc.DrawLatex(0.1,pos_y + 5*oset,input_dm23sq_str  ); 


  TLegend *legY_reac = new TLegend(0.6, 0.4, 0.9, 0.9);
  //legY_reac->SetFillColor(0);
  legY_reac->SetFillStyle(0);
  legY_reac->SetBorderSize(0);
  legY_reac->AddEntry(hErecoMode_NC,"NC","F");
  legY_reac->AddEntry(hErecoMode_MEC,"MEC","F");
  legY_reac->AddEntry(hErecoMode_CCQE,"CCQE","F");
  legY_reac->AddEntry(hErecoMode_CCnQE,"CCnQE","F");
  legY_reac->Draw();


  THStack *sErecoTotFlavMode_rebin = new THStack("sErecoTotFlavMode_rebin","sErecoTotFlavMode_rebin");

  TH1D * hErecoMode_NueBeam_rebin         = hErecoMode_NueBeam         ->Rebin(10,"hErecoMode_NueBeam_rebin"         ,xbins);
  TH1D * hErecoMode_NumuBeam_rebin        = hErecoMode_NumuBeam        ->Rebin(10,"hErecoMode_NumuBeam_rebin"        ,xbins);
  TH1D * hErecoMode_NumuBarBeam_rebin     = hErecoMode_NumuBarBeam     ->Rebin(10,"hErecoMode_NumuBarBeam_rebin"     ,xbins);
  TH1D * hErecoMode_NumuToNue_rebin       = hErecoMode_NumuToNue       ->Rebin(10,"hErecoMode_NumuToNue_rebin"       ,xbins);
  TH1D * hErecoMode_NumuBarToNueBar_rebin = hErecoMode_NumuBarToNueBar ->Rebin(10,"hErecoMode_NumuBarToNueBar_rebin" ,xbins);

  sErecoTotFlavMode_rebin->Add(hErecoMode_NC_rebin);
  sErecoTotFlavMode_rebin->Add(hErecoMode_NueBeam_rebin);
  sErecoTotFlavMode_rebin->Add(hErecoMode_NumuBeam_rebin);
  sErecoTotFlavMode_rebin->Add(hErecoMode_NumuBarBeam_rebin);
  sErecoTotFlavMode_rebin->Add(hErecoMode_NumuToNue_rebin);
  sErecoTotFlavMode_rebin->Add(hErecoMode_NumuBarToNueBar_rebin);
  
  TCanvas *cErecoTotFlavMode_rebin = new TCanvas("cErecoTotFlavMode_rebin","cErecoTotFlavMode_rebin");
  sErecoTotFlavMode_rebin->Draw();
  if(draw_data){
    h_data_ereco->Draw("pe same");
    lScaleEreco->Draw("same");
  }
  sErecoTotFlavMode_rebin->GetXaxis()->SetRangeUser(ErecoRangeLow,ErecoRangeUp); // ereco    
  sErecoTotFlavMode_rebin->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  if(plot_per_50MeV) sErecoTotFlavMode_rebin->GetYaxis()->SetTitle("Events per 0.050 GeV");
  else sErecoTotFlavMode_rebin->GetYaxis()->SetTitle("Number of events");
  cErecoTotFlavMode_rebin->Modified();
  gStyle->SetOptStat(0);

  TLegend *legY = new TLegend(0.6, 0.4, 0.9, 0.9);
  //legY->SetFillColor(0);
  legY->SetFillStyle(0);
  legY->SetBorderSize(0);
  legY->AddEntry(hErecoMode_NumuBarToNueBar,"Total","L");
  legY->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","F");
  legY->AddEntry(hErecoMode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","F");
  legY->AddEntry(hErecoMode_NC,"NC","F");
  legY->AddEntry(hErecoMode_NueBeam,"#nu_{e}/#bar{#nu}_{e} intrinsic","F");
  legY->AddEntry(hErecoMode_NumuBeam,"#nu_{#mu} intrinsic","F");
  legY->AddEntry(hErecoMode_NumuBarBeam,"#bar{#nu}_{#mu} intrinsic","F");
  legY->Draw();

  
  // Draw 1D templates: Q2

  if(plot_per_50MeV){
    hQ2Mode_NumuBarBeam = spectrum_to_per_50mev(hQ2Mode_NumuBarBeam); 
    hQ2Mode_NumuBeam = spectrum_to_per_50mev(hQ2Mode_NumuBeam); 
    hQ2Mode_NueBeam = spectrum_to_per_50mev(hQ2Mode_NueBeam); 
    hQ2Mode_NumuBarToNueBar = spectrum_to_per_50mev(hQ2Mode_NumuBarToNueBar); 
    hQ2Mode_NumuToNue = spectrum_to_per_50mev(hQ2Mode_NumuToNue); 
    hQ2Mode_NC = spectrum_to_per_50mev(hQ2Mode_NC); 
    hQ2Mode_MEC = spectrum_to_per_50mev(hQ2Mode_MEC); 
    hQ2Mode_CCQE = spectrum_to_per_50mev(hQ2Mode_CCQE); 
    hQ2Mode_CCnQE = spectrum_to_per_50mev(hQ2Mode_CCnQE); 
  }

  hQ2Mode_NumuBarBeam->SetFillStyle(1001);
  hQ2Mode_NumuBarBeam->SetLineColor(kBlack);
  hQ2Mode_NumuBarBeam->SetFillColor(kBlue);
  hQ2Mode_NumuBarBeam->SetLineWidth(2);
  if(draw_data) hQ2Mode_NumuBarBeam->Scale(scale_q2);

  hQ2Mode_NumuBeam->SetFillStyle(1001);
  hQ2Mode_NumuBeam->SetLineColor(kBlack);
  hQ2Mode_NumuBeam->SetFillColor(kCyan);
  hQ2Mode_NumuBeam->SetLineWidth(2);
  if(draw_data) hQ2Mode_NumuBeam->Scale(scale_q2);

  hQ2Mode_NueBeam->SetFillStyle(1001);
  hQ2Mode_NueBeam->SetLineColor(kBlack);
  hQ2Mode_NueBeam->SetFillColor(kYellow);
  hQ2Mode_NueBeam->SetLineWidth(2);
  if(draw_data) hQ2Mode_NueBeam->Scale(scale_q2);

  hQ2Mode_NumuBarToNueBar->SetFillStyle(1001);
  hQ2Mode_NumuBarToNueBar->SetLineColor(kBlack);
  hQ2Mode_NumuBarToNueBar->SetFillColor(kGreen);
  hQ2Mode_NumuBarToNueBar->SetLineWidth(2);
  if(draw_data) hQ2Mode_NumuBarToNueBar->Scale(scale_q2);

  hQ2Mode_NumuToNue->SetFillStyle(3004);
  hQ2Mode_NumuToNue->SetLineColor(kBlack);
  hQ2Mode_NumuToNue->SetFillColor(kGray+3);
  hQ2Mode_NumuToNue->SetLineWidth(2);
  if(draw_data) hQ2Mode_NumuToNue->Scale(scale_q2);

  hQ2Mode_NC->SetFillStyle(1001);
  hQ2Mode_NC->SetLineColor(kBlack);
  hQ2Mode_NC->SetFillColor(kRed);
  hQ2Mode_NC->SetLineWidth(2);
  if(draw_data) hQ2Mode_NC->Scale(scale_q2);
  
  hQ2Mode_MEC->SetFillStyle(1001);
  hQ2Mode_MEC->SetLineColor(kBlack);
  hQ2Mode_MEC->SetFillColor(kBlue);
  hQ2Mode_MEC->SetLineWidth(2);
  if(draw_data) hQ2Mode_MEC->Scale(scale_q2);
  
  hQ2Mode_CCQE->SetFillStyle(1001);
  hQ2Mode_CCQE->SetLineColor(kBlack);
  hQ2Mode_CCQE->SetFillColor(kGreen);
  hQ2Mode_CCQE->SetLineWidth(2);
  if(draw_data) hQ2Mode_CCQE->Scale(scale_q2);
  
  hQ2Mode_CCnQE->SetFillStyle(1001);
  hQ2Mode_CCnQE->SetLineColor(kBlack);
  hQ2Mode_CCnQE->SetFillColor(kYellow);
  hQ2Mode_CCnQE->SetLineWidth(2);
  if(draw_data) hQ2Mode_CCnQE->Scale(scale_q2);
  
 
  if(plot_per_50MeV){
    hQ2Mode_NumuBarBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NumuBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NueBeam->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NumuToNue->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NumuBarToNueBar->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_NC->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_MEC->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_CCQE->GetYaxis()->SetTitle("Events per 0.050 GeV");
    hQ2Mode_CCnQE->GetYaxis()->SetTitle("Events per 0.050 GeV");
  }
  else{
    hQ2Mode_NumuBarBeam->GetXaxis()->SetTitle(titleX);
    hQ2Mode_NumuBeam->GetXaxis()->SetTitle(titleX);
    hQ2Mode_NueBeam->GetXaxis()->SetTitle(titleX);
    hQ2Mode_NumuToNue->GetXaxis()->SetTitle(titleX);
    hQ2Mode_NumuBarToNueBar->GetXaxis()->SetTitle(titleX);
    hQ2Mode_NC->GetXaxis()->SetTitle(titleX);
    hQ2Mode_MEC->GetXaxis()->SetTitle(titleX);
    hQ2Mode_CCQE->GetXaxis()->SetTitle(titleX);
    hQ2Mode_CCnQE->GetXaxis()->SetTitle(titleX);
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
  if(draw_data){
    h_data_q2->Draw("pe same");
    lScaleQ2->Draw("same");
  }
  sQ2TotFlavMode->GetXaxis()->SetTitle(titleX);
  if(plot_per_50MeV) sQ2TotFlavMode->GetYaxis()->SetTitle("Events per 0.050 GeV");
  else sQ2TotFlavMode->GetYaxis()->SetTitle("Number of events");
  cQ2TotFlavMode->Modified();
  gStyle->SetOptStat(0);

  //legX->SetFillColor(0);
  legX->SetFillStyle(0);
  legX->SetBorderSize(0);
  legX->AddEntry(hQ2Mode_NumuBarToNueBar,"Total","L");
  legX->AddEntry(hQ2Mode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","F");
  legX->AddEntry(hQ2Mode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","F");
  legX->AddEntry(hQ2Mode_NC,"NC","F");
  legX->AddEntry(hQ2Mode_NueBeam,"#nu_{e}/#bar{#nu}_{e} intrinsic","F");
  legX->AddEntry(hQ2Mode_NumuBeam,"#nu_{#mu} intrinsic","F");
  legX->AddEntry(hQ2Mode_NumuBarBeam,"#bar{#nu}_{#mu} intrinsic","F");
  legX->Draw();




  THStack *sQ2TotReac = new THStack("sQ2TotReac","sQ2TotReac");
  sQ2TotReac->Add(hQ2Mode_NC);
  sQ2TotReac->Add(hQ2Mode_MEC);
  sQ2TotReac->Add(hQ2Mode_CCQE);
  sQ2TotReac->Add(hQ2Mode_CCnQE);

  TCanvas *cQ2TotReac = new TCanvas("cQ2TotReac","cQ2TotReac");
  sQ2TotReac->Draw();
  if(draw_data){
    h_data_q2->Draw("pe same");
    lScaleQ2->Draw("same");
  }
  sQ2TotReac->GetXaxis()->SetTitle(titleX);
  if(plot_per_50MeV) sQ2TotReac->GetYaxis()->SetTitle("Events per 0.050 GeV");
  else sQ2TotReac->GetYaxis()->SetTitle("Number of events");
  cQ2TotReac->Modified();
  gStyle->SetOptStat(0);

  //legX_reac->SetFillColor(0);
  legX_reac->SetFillStyle(0);
  legX_reac->SetBorderSize(0);
  legX_reac->AddEntry(hQ2Mode_NC,"NC","F");
  legX_reac->AddEntry(hQ2Mode_MEC,"MEC","F");
  legX_reac->AddEntry(hQ2Mode_CCQE,"CCQE","F");
  legX_reac->AddEntry(hQ2Mode_CCnQE,"CCnQE","F");
  legX_reac->Draw();


 
  // // Draw Sig/Bkg ratio plots (likelihood)
  // // nexp = Bkg only, nobs = Sig + Bkg
  // TH2D *hSig;
  // TH2D *hSigpBkg;
  // if(s==kSKsFCFV1RingELike){
  //   hSig = (TH2D*) hErecoVsQ2Mode_NumuToNue->Clone("hSignal");
  // }
  // else if(s==kSKsFCFV1RingMuLike){
  //   hSig = (TH2D*) hErecoVsQ2Mode_NumuBeam->Clone("hSignal");
  // }
  // else if(s==kSKsRHCFCFV1RingELike){
  //   hSig = (TH2D*) hErecoVsQ2Mode_NumuBarToNueBar->Clone("hSignal");
  // }
  // else if(s==kSKsRHCFCFV1RingMuLike){
  //   hSig = (TH2D*) hErecoVsQ2Mode_NumuBarBeam->Clone("hSignal");
  // }  
  // TH2D *hLikeSigBkg = (TH2D*) hSig->Clone("hLikeSigBkg");
  // hLikeSigBkg->Reset();
  // int Nq2 = hErecoVsQ2Tot->GetXaxis()->GetNbins();
  // int Nereco = hErecoVsQ2Tot->GetYaxis()->GetNbins();  
  // for(int iq2=0;iq2<Nq2;iq2++){
  //   for(int iereco=0;iereco<Nereco;iereco++){
  //     double ntot = hErecoVsQ2Tot->GetBinContent(iq2+1,iereco+1);
  //     double nsig = hSig->GetBinContent(iq2+1,iereco+1);
  //     double nbkg = ntot - nsig;
  //     double nexp = nbkg;
  //     double nobs = ntot;
  //     double like = valor::util::ChisqPoisson(nobs,nexp);
  //     hLikeSigBkg->SetBinContent(iq2+1,iereco+1,like);
  //   }
  // }
  // TCanvas *cLikeSigBkg = new TCanvas("cLikeSigBkg","cLikeSigBkg");
  // hLikeSigBkg->Draw("colz");
  // if(draw_data){
  //   h_data->Draw("same");
  // }
  // hLikeSigBkg->GetZaxis()->SetTitle("#Delta #chi^{2}");
  // hLikeSigBkg_ProjEreco = (TH1D*) hLikeSigBkg->ProjectionY()->Clone("hLikeSigBkg_ProjEreco");
  // TCanvas *cLikeSigBkg_ProjEreco = new TCanvas("cLikeSigBkg_ProjEreco","cLikeSigBkg_ProjEreco");
  // hLikeSigBkg_ProjEreco->Draw();
  // if(draw_data){
  //   h_data_ereco->Draw("pe same");
  // }
  // hLikeSigBkg_ProjEreco->GetYaxis()->SetTitle("#Delta #chi^{2}");
  // hLikeSigBkg_ProjQ2 = (TH1D*) hLikeSigBkg->ProjectionX()->Clone("hLikeSigBkg_ProjQ2");
  // TCanvas *cLikeSigBkg_ProjQ2 = new TCanvas("cLikeSigBkg_ProjQ2","cLikeSigBkg_ProjQ2");
  // hLikeSigBkg_ProjQ2->Draw();
  // if(draw_data){
  //   h_data_q2->Draw("pe same");
  // }
  // hLikeSigBkg_ProjQ2->GetYaxis()->SetTitle("#Delta #chi^{2}");
  
  

  


  // Get total # of events

  double TotEvts = hErecoVsQ2Tot->Integral();

  double TotEvtsCheck = 
    hErecoVsQ2Mode_NC->Integral() + 
    hErecoVsQ2Mode_NumuBarBeam->Integral() + 
    hErecoVsQ2Mode_NumuBeam->Integral() + 
    hErecoVsQ2Mode_NueBeam->Integral() + 
    hErecoVsQ2Mode_NumuBarToNueBar->Integral() + 
    hErecoVsQ2Mode_NumuToNue->Integral();

  double TotEvtsCheck_reac = 
    hErecoVsQ2Mode_NC->Integral() + 
    hErecoVsQ2Mode_MEC->Integral() + 
    hErecoVsQ2Mode_CCQE->Integral() + 
    hErecoVsQ2Mode_CCnQE->Integral();
  
  if( fabs(TotEvts-TotEvtsCheck)>0.0000001 ){
    cout << endl;
    cout << "TotEvts = " << TotEvts << endl;
    cout << "TotEvtsCheck = " << TotEvtsCheck << endl;
    cout << "There is an error!" << endl;
    cout << endl;
    exit(1);
  }

  if( fabs(TotEvts-TotEvtsCheck_reac)>0.0000001 ){
    cout << endl;
    cout << "TotEvts = " << TotEvts << endl;
    cout << "TotEvtsCheck_reac = " << TotEvtsCheck_reac << endl;
    cout << "There is an error!" << endl;
    cout << endl;
    exit(1);
  }

  cout << endl;
  cout << "Numu --> Nue = " << hErecoVsQ2Mode_NumuToNue->Integral() << endl;
  cout << "NumuBar --> NueBar = " << hErecoVsQ2Mode_NumuBarToNueBar->Integral() << endl; 
  cout << "Numu beam = " << hErecoVsQ2Mode_NumuBeam->Integral() << endl; 
  cout << "NumuBar beam = " << hErecoVsQ2Mode_NumuBarBeam->Integral() << endl; 
  cout << "Nue beam = " << hErecoVsQ2Mode_NueBeam->Integral() << endl; 
  cout << "NC = " << hErecoVsQ2Mode_NC->Integral() << endl; 
  cout << "Total = " << TotEvts << endl;
  cout << endl;
  cout << "NC = " << hErecoVsQ2Mode_NC->Integral() << endl; 
  cout << "MEC = " << hErecoVsQ2Mode_MEC->Integral() << endl; 
  cout << "CCQE = " << hErecoVsQ2Mode_CCQE->Integral() << endl; 
  cout << "CCnQE = " << hErecoVsQ2Mode_CCnQE->Integral() << endl; 
  cout << endl;

  if(draw_data){
    cout << endl;
    cout << "Draw data spectra" << endl;
    cout << "Scale factor use for MC spectra is" << endl;
    cout << "Ereco --> " << scale_ereco << endl;
    cout << "scale_q2 --> " << scale_q2 << endl;
    cout << endl;
  }
  else {
    cout << endl;
    cout << "Don't draw data spectra" << endl;
    if(!draw_data_spectra){
      cout << "Chosen not to draw data spectra" << endl;
    }
    else{
      cout << "Input file " << rootfile << " doesn't exist" << endl;
      cout << endl;
    }
  }

   
  

  // Store spectra

  TFile *fout;
  if(store_spectrum){
    
    char name[200];
    sprintf(name,"%s_%s.root",tag,util::AsStringCompact(s).c_str());
    fout = new TFile(name,"RECREATE");
    
    hErecoVsQ2Tot->Write();
    hQ2Tot->Write();
    hErecoTot->Write();

    hErecoVsQ2Mode_NueBeam->Write();
    hErecoVsQ2Mode_NumuBeam->Write();
    hErecoVsQ2Mode_NumuBarBeam->Write();
    hErecoVsQ2Mode_NumuBarToNueBar->Write();
    hErecoVsQ2Mode_NumuToNue->Write();
    hErecoVsQ2Mode_NC->Write();
    hErecoVsQ2Mode_MEC->Write();
    hErecoVsQ2Mode_CCQE->Write();
    hErecoVsQ2Mode_CCnQE->Write();
    
    hErecoMode_NueBeam->Write();
    hErecoMode_NumuBeam->Write();
    hErecoMode_NumuBarBeam->Write();
    hErecoMode_NumuBarToNueBar->Write();
    hErecoMode_NumuToNue->Write();
    hErecoMode_NC->Write();
    hErecoMode_MEC->Write();
    hErecoMode_CCQE->Write();
    hErecoMode_CCnQE->Write();
    
    hQ2Mode_NueBeam->Write();
    hQ2Mode_NumuBeam->Write();
    hQ2Mode_NumuBarBeam->Write();
    hQ2Mode_NumuBarToNueBar->Write();
    hQ2Mode_NumuToNue->Write();
    hQ2Mode_NC->Write();
    hQ2Mode_MEC->Write();
    hQ2Mode_CCQE->Write();
    hQ2Mode_CCnQE->Write();

    fout->Close();
  }



  // Print canvases

  if(print_spectrum){
    
    char name[200];

    sprintf(name,"%s_%s.pdf","cErecoVsQ2Tot",sk_sample);    
    cErecoVsQ2Tot ->Print(name);

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
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_NC",sk_sample);    
    cErecoVsQ2_NC ->Print(name);
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_MEC",sk_sample);    
    cErecoVsQ2_MEC ->Print(name);
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_CCQE",sk_sample);    
    cErecoVsQ2_CCQE ->Print(name);
    sprintf(name,"%s_%s.pdf","cErecoVsQ2_CCnQE",sk_sample);    
    cErecoVsQ2_CCnQE ->Print(name);

    sprintf(name,"%s_%s.pdf","cQ2TotFlavMode",sk_sample);    
    cQ2TotFlavMode ->Print(name);

    sprintf(name,"%s_%s.pdf","cErecoTotFlavMode",sk_sample);    
    cErecoTotFlavMode ->Print(name);

    sprintf(name,"%s_%s.pdf","cErecoTotFlavMode_rebin",sk_sample);    
    cErecoTotFlavMode_rebin ->Print(name);
  }

  cout << "test_sk_ereco_spectrum_nominal_mainmode_joint_2D.C finished!" << endl;
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


