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

void event_rate_ellipse_fhc_Vs_rhc
(
 const string str_preliminary = "T2K Run 1-7c preliminary",

 const char * tag                  = "ErecoVsTheta_Run17_Asimov_A_dcpPi",
 const char * sk_sample            = "RHC1Re",                         /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
 const char * hypothesis           = "3f_osc",                 /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_osc_normbar, 3f_oscbar_normbar" */ 
 
 // // No Oscillations
 // double       input_mixing12       = 0.,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0., //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 0., //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar        = -99,

 // Asimov A
 double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13       = 0.0217,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 //double       input_deltacp        = 0., //-99,    /* deltacp, */
 //double       input_deltacp        = +TMath::Pi(), //-99,    /* deltacp, */
 double       input_deltacp        = -1.601, //-99,    /* deltacp, */
 //double       input_deltacp        = +TMath::Pi()/2, //-99,    /* deltacp, */
 double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 double       input_betabar        = 1., //0., //-99,

 // // Asimov B
 // double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13       = 0.0217,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23       = 0.45, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_deltacp        = 0., //-99,    /* deltacp, */
 // double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 // double       input_betabar        = 1., //0., //-99,

 bool         disable_banff_tuning = false,                          /* */

 //
 // NEW CONFIGURATION FILES (13/7/16 BANFF RELEASE W/ FIT BIAS AND ADDITIONAL UNCERTAINTY ON RUN 1-7c FLUX, FINAL)   
 //
 const char * config_file = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_Run17c_final.xml",


 //
 // NEW CONFIGURATION FILES (4/6/16 BANFF RELEASE W/ FIT BIAS AND ADDITIONAL UNCERTAINTY ON RUN 1-7 FLUX)
 //
 //const char * config_file      = "analysis_configuration_data_theta_1bin_NIWG2015c_1p1h_fitbias_FluxTunRun17.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_1bin_NIWG2015c_1p1h_fitbias_FluxTunRun17_prefit.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17_prefit.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_Run17.xml", // Run 1-6 flux tuning
 
 //
 // NEW CONFIGURATION FILES (14/3/16 BANFF RELEASE W/ FIT BIAS)
 // 
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias.xml",
 //const char * config_file      = "analysis_configuration_data_theta_elike_fhc_only_RFGRPA_NIWG2015c_1p1h_fitbias.xml",
 //const char * config_file      = "analysis_configuration_data_theta_elike_rhc_only_RFGRPA_NIWG2015c_1p1h_fitbias.xml",
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_ExpSummer2016.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_NIWG2015c_1p1h_fitbias.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FSIOnly.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_SKDetOnly.xml",
 
 bool         store_spectrum      = true,
 bool         print_spectrum      = true,                          /* */
 const char * mixparamchoice      = "sinsqtheta",
 const char * hierarchy           = "normal"   ,
 const char * dmsqscheme          = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

 bool plot_per_50MeV = false,
 double ErecoRangeLow = 0.,
 double ErecoRangeUp = 30.,



 int          iexperiment         = 0, //4404, //4739,                       /* */
 //const char * rootfile = "data_run1to4/SKdata_Run14.root",
 //const char * rootfile = "data_run1to6/SKdata_Run16.root", 
 //const char * rootfile = "data_run1to7b/SKdata_Run17.root", 
 const char * rootfile = "data_run1to7c/SKdata_Run17c.root", 
 const bool draw_data_spectra = true,
 const bool scale_mc = false,

 const bool change_banff_params = false


)
{
  using namespace std;
  using namespace valor;

  char name[300];

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
  cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;
 
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
    skexp.GetBANFF()->RestoreNominalMC();
  }

  //
  // - Open loop over values of deltaCP (each octant and mass hierarchy)
  //   and fill TGraph with a point of N_fhc_elike Vs N_rhc_elike
  // - Set lines of the TGraph
  // - Calculate N_fhc_elike Vs N_rhc_elike point from data and plot it on top of the MC one
  // 

    
  double deltaCP_first = -TMath::Pi();
  double deltaCP_last = +TMath::Pi();
  double dcp_npt = 9; //201;
  double deltaCP_width = (deltaCP_last - deltaCP_first) / (dcp_npt-1);
  
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p5 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p5 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p45 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p45 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p55 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p55 = new TGraph();

  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p5_4pt = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p5_4pt = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p45_4pt = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p45_4pt = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p55_4pt = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p55_4pt = new TGraph();

  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi = new TGraph();

  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2 = new TGraph();

  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_0 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_0 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_0 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0 = new TGraph();

  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_mpi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2 = new TGraph();
  TGraph *gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2 = new TGraph();


  for(int idcp=0;idcp<dcp_npt;idcp++){

    double Nfhc = -999;
    double Nrhc = -999;
 
    double dcp = deltaCP_first + idcp * deltaCP_width;
    
    cout << "dcp = " << dcp << endl;

    // 
    // sinsqth23 = 0.45 - Normal Hierarchy
    //
    double ss23 = 0.45;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhNormal);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23,ss23);
    // Generate prediction
    skexp.GetSuperkMC()->SetXMEC(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    gN_fhcVsrhc_elike_nh_sinsqth23_0p45->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p45->GetN(),Nfhc,Nrhc);

    if( fabs(dcp-(TMath::Pi())) <1.E-4 ||
  	fabs(dcp-(TMath::Pi()/2.)) <1.E-4 ||
  	fabs(dcp)<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()/2.))<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()))   <1.E-4 ){
      gN_fhcVsrhc_elike_nh_sinsqth23_0p45_4pt->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p45_4pt->GetN(),Nfhc,Nrhc);

      if(fabs(dcp-(-TMath::Pi())) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(-TMath::Pi()/2.)) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp)<1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_0->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_0->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(TMath::Pi()/2)) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2->GetN(),Nfhc,Nrhc);

      cout << "NH: sinsqth23 = " << ss23 << ", Nfhc = " << Nfhc << ", Nrhc = " << Nrhc << endl;
    }

    
    //
    // sinsqth23 = 0.5 - Normal Hierarchy
    //
    double ss23 = 0.5;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhNormal);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    //Generate prediction
    skexp.GetSuperkMC()->SetXMEC(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    gN_fhcVsrhc_elike_nh_sinsqth23_0p5->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p5->GetN(),Nfhc,Nrhc); 
    if( fabs(dcp-(TMath::Pi())) <1.E-4 ||
  	fabs(dcp-(TMath::Pi()/2.)) <1.E-4 ||
  	fabs(dcp)<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()/2.))<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()))   <1.E-4 ){
      
      gN_fhcVsrhc_elike_nh_sinsqth23_0p5_4pt->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_4pt->GetN(),Nfhc,Nrhc);
      
      if(fabs(dcp-(-TMath::Pi())) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(-TMath::Pi()/2.)) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp)<1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(TMath::Pi()/2)) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2->GetN(),Nfhc,Nrhc);

      cout << "NH: sinsqth23 = " << ss23 << ", Nfhc = " << Nfhc << ", Nrhc = " << Nrhc << endl;
    }
    
    // 
    // sinsqth23 = 0.55 - Normal Hierarchy
    //
    double ss23 = 0.55;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhNormal);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    // Generate prediction
    skexp.GetSuperkMC()->SetXMEC(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    gN_fhcVsrhc_elike_nh_sinsqth23_0p55->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p55->GetN(),Nfhc,Nrhc);
    if( fabs(dcp-(TMath::Pi())) <1.E-4 ||
  	fabs(dcp-(TMath::Pi()/2.)) <1.E-4 ||
  	fabs(dcp)<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()/2.))<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()))   <1.E-4 ){
      gN_fhcVsrhc_elike_nh_sinsqth23_0p55_4pt->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p55_4pt->GetN(),Nfhc,Nrhc);

      if(fabs(dcp-(-TMath::Pi())) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(-TMath::Pi()/2.)) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp)<1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_0->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_0->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(TMath::Pi()/2)) <1.E-4 ) 
	gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2->SetPoint(gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2->GetN(),Nfhc,Nrhc);

      cout << "NH: sinsqth23 = " << ss23 << ", Nfhc = " << Nfhc << ", Nrhc = " << Nrhc << endl; 
   }

    // 
    // sinsqth23 = 0.45 - Inverted Hierarchy
    //
    double ss23 = 0.45;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhInverted);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    // Generate prediction
    skexp.GetSuperkMC()->SetXMEC(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    gN_fhcVsrhc_elike_ih_sinsqth23_0p45->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p45->GetN(),Nfhc,Nrhc);
    if( fabs(dcp-(TMath::Pi())) <1.E-4 ||
  	fabs(dcp-(TMath::Pi()/2.)) <1.E-4 ||
  	fabs(dcp)<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()/2.))<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()))   <1.E-4 ){
      gN_fhcVsrhc_elike_ih_sinsqth23_0p45_4pt->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p45_4pt->GetN(),Nfhc,Nrhc);
 
      if(fabs(dcp-(-TMath::Pi())) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(-TMath::Pi()/2.)) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp)<1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(TMath::Pi()/2)) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->GetN(),Nfhc,Nrhc);

      cout << "IH: sinsqth23 = " << ss23 << ", Nfhc = " << Nfhc << ", Nrhc = " << Nrhc << endl;
    }

    //
    // sinsqth23 = 0.5 - Inverted Hierarchy
    //
    double ss23 = 0.5;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhInverted);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    //Generate prediction
    skexp.GetSuperkMC()->SetXMEC(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    gN_fhcVsrhc_elike_ih_sinsqth23_0p5->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p5->GetN(),Nfhc,Nrhc);   
    if( fabs(dcp-(TMath::Pi())) <1.E-4 ||
  	fabs(dcp-(TMath::Pi()/2.)) <1.E-4 ||
  	fabs(dcp)<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()/2.))<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()))   <1.E-4 ){
      gN_fhcVsrhc_elike_ih_sinsqth23_0p5_4pt->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p5_4pt->GetN(),Nfhc,Nrhc);

      if(fabs(dcp-(-TMath::Pi())) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(-TMath::Pi()/2.)) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_mpi2->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_mpi2->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp)<1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_0->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_0->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(TMath::Pi()/2)) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2->GetN(),Nfhc,Nrhc);

      cout << "IH: sinsqth23 = " << ss23 << ", Nfhc = " << Nfhc << ", Nrhc = " << Nrhc << endl;
    }
 
    // 
    // sinsqth23 = 0.55 - Inverted Hierarchy
    //
    double ss23 = 0.55;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhInverted);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    // Generate prediction
    skexp.GetSuperkMC()->SetXMEC(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    gN_fhcVsrhc_elike_ih_sinsqth23_0p55->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p55->GetN(),Nfhc,Nrhc);
    if( fabs(dcp-(TMath::Pi())) <1.E-4 ||
  	fabs(dcp-(TMath::Pi()/2.)) <1.E-4 ||
  	fabs(dcp)<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()/2.))<1.E-4 ||
  	fabs(dcp-(-TMath::Pi()))   <1.E-4 ){
      gN_fhcVsrhc_elike_ih_sinsqth23_0p55_4pt->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p55_4pt->GetN(),Nfhc,Nrhc);

      if(fabs(dcp-(-TMath::Pi())) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(-TMath::Pi()/2.)) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp)<1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0->GetN(),Nfhc,Nrhc);
      else if(fabs(dcp-(TMath::Pi()/2)) <1.E-4 ) 
	gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->SetPoint(gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->GetN(),Nfhc,Nrhc);

      cout << "IH: sinsqth23 = " << ss23 << ", Nfhc = " << Nfhc << ", Nrhc = " << Nrhc << endl;
    }

    // // 1Re # of events 
    // s = util::SKSampleFromString("1Re"); // requested sample
    // if(s == conv::kSKsUndefined) {
    //   cerr << "** Undefined sample!" << endl;
    //   return;
    // }
    // is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
    // hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
    // double Nfhc = hErecoVsQ2Tot->Integral();    
    // // RHC1Re # of events 
    // s = util::SKSampleFromString("RHC1Re"); // requested sample
    // if(s == conv::kSKsUndefined) {
    //   cerr << "** Undefined sample!" << endl;
    //   return;
    // }
    // is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
    // hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
    // double Nrhc = hErecoVsQ2Tot->Integral();

    cout << endl;
  }

  



  // Get data spectra and rescale MC to data

  TGraphErrors *gN_fhcVsrhc_elike_data = new TGraphErrors();
  if(draw_data){
    bool using_xml_file = false;
    skobs_fakedata->Read(rootfile2.Data(),iexperiment,using_xml_file);

    // 1Re # of events 
    s = util::SKSampleFromString("1Re"); // requested sample
    if(s == conv::kSKsUndefined) {
      cerr << "** Undefined sample!" << endl;
      return;
    }
    is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
    double nfhc = skobs_fakedata->TwoDSpectrum(s)->Integral();

    // RHC1Re # of events 
    s = util::SKSampleFromString("RHC1Re"); // requested sample
    if(s == conv::kSKsUndefined) {
      cerr << "** Undefined sample!" << endl;
      return;
    }
    is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
    double nrhc = skobs_fakedata->TwoDSpectrum(s)->Integral();
    
    nfhc = 32;
    nrhc = 4;
    
    double nfhc_err = sqrt(nfhc);
    double nrhc_err = sqrt(nrhc);

    cout << endl;
    cout << "Data: nfhc = " << nfhc << " +- " << nfhc_err << ", nrhc = " << nrhc << " +- " << nrhc_err << endl;
    cout << endl;

    int ipt = gN_fhcVsrhc_elike_data->GetN();
    gN_fhcVsrhc_elike_data->SetPoint(ipt,nfhc,nrhc);
    gN_fhcVsrhc_elike_data->SetPointError(ipt,nfhc_err,nrhc_err);
  }
 

  // Init histogram
  
  TH2D *hInit = new TH2D("hInit","hInit",500,0,50,500,0,50);
  hInit->GetXaxis()->SetRangeUser(10.,35.);
  hInit->GetYaxis()->SetRangeUser(3.,10.1);
  hInit->GetXaxis()->SetTitle("Number of #nu_{e} candidates");
  hInit->GetYaxis()->SetTitle("Number of #bar{#nu}_{e} candidates");


  // Draw 

  TCanvas *cEllipse_nh_sinsqth23 = new TCanvas("cEllipse_nh_sinsqth23","cEllipse_nh_sinsqth23");
  hInit->Draw();

  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi->SetMarkerStyle(24);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi->SetMarkerStyle(24);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi->SetMarkerStyle(24);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->SetMarkerStyle(20);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2->SetMarkerStyle(20);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->SetMarkerStyle(20);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_0->SetMarkerStyle(25);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0->SetMarkerStyle(25);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_0->SetMarkerStyle(25);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2->SetMarkerStyle(21);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2->SetMarkerStyle(21);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2->SetMarkerStyle(21);

  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi->SetMarkerStyle(24);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi->SetMarkerStyle(24);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi->SetMarkerStyle(24);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2->SetMarkerStyle(20);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_mpi2->SetMarkerStyle(20);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2->SetMarkerStyle(20);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0->SetMarkerStyle(25);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_0->SetMarkerStyle(25);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0->SetMarkerStyle(25);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->SetMarkerStyle(21);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2->SetMarkerStyle(21);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->SetMarkerStyle(21);

  // gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi->SetMarkerColor(kAzure+6);
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->SetMarkerColor(kAzure+6); 
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_0->SetMarkerColor(kAzure+6);
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2->SetMarkerColor(kAzure+6);
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi->SetMarkerColor(kAzure-6);
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->SetMarkerColor(kAzure-6); 
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_0->SetMarkerColor(kAzure-6);
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2->SetMarkerColor(kAzure-6);
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi->SetMarkerColor(kPink+6);
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2->SetMarkerColor(kPink+6); 
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0->SetMarkerColor(kPink+6);
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->SetMarkerColor(kPink+6);
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi->SetMarkerColor(kPink-6);
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2->SetMarkerColor(kPink-6); 
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0->SetMarkerColor(kPink-6);
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->SetMarkerColor(kPink-6);
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p45->SetLineColor(kAzure+6);
  // gN_fhcVsrhc_elike_nh_sinsqth23_0p55->SetLineColor(kAzure-6);
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p45->SetLineColor(kPink+6);
  // gN_fhcVsrhc_elike_ih_sinsqth23_0p55->SetLineColor(kPink-6);

  gN_fhcVsrhc_elike_nh_sinsqth23_0p45->SetLineColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55->SetLineColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->SetMarkerColor(kBlue); 
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_0->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->SetMarkerColor(kBlue); 
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_0->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45->SetLineColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55->SetLineColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2->SetMarkerColor(kRed); 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2->SetMarkerColor(kRed); 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->SetMarkerColor(kRed);
  
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45->SetLineStyle(1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55->SetLineStyle(1);

  gN_fhcVsrhc_elike_nh_sinsqth23_0p45->Draw("l");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55->Draw("l");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45->Draw("l");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55->Draw("l");

  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->Draw("p"); 
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_0->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->Draw("p"); 
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_0->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2->Draw("p");

  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2->Draw("p"); 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2->Draw("p"); 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->Draw("p");

  gN_fhcVsrhc_elike_data->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_data->SetMarkerSize(1);
  gN_fhcVsrhc_elike_data->SetMarkerStyle(4);
  gN_fhcVsrhc_elike_data->Draw("PE");
  
  double x1;
  double x2;
  double y1;
  double y2;

  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2->GetPoint(0,x1,y1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2->GetPoint(0,x2,y2);
  double Xfirst = (x1+x2) / 2. * 0.7; 
  double Yfirst = (y1+y2) / 2. * 0.7; 
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2->GetPoint(0,x1,y1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2->GetPoint(0,x2,y2);
  double Xlast = (x1+x2) / 2. * 1.2; 
  double Ylast = (y1+y2) / 2. * 1.2; 
  TLine *lMassHierarchy = new TLine(Xfirst,Yfirst,Xlast,Ylast);
  lMassHierarchy->SetLineColor(kBlack);
  lMassHierarchy->SetLineWidth(2);
  lMassHierarchy->SetLineStyle(2);
  //lMassHierarchy->Draw();
  
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->GetPoint(0,x1,y1);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->GetPoint(0,x2,y2);
  Xfirst = (x1+x2) / 2.; // * 0.7; 
  Yfirst = (y1+y2) / 2.; // * 0.7; 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->GetPoint(0,x1,y1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->GetPoint(0,x2,y2);
  Xlast = (x1+x2) / 2.; // * 1.2; 
  Ylast = (y1+y2) / 2.; // * 1.2;   

  
  TLine *lOctant = new TLine(Xfirst,Yfirst,Xlast,Ylast);
  lOctant->SetLineColor(kBlack);
  lOctant->SetLineWidth(2);
  lOctant->SetLineStyle(2);
  lOctant->Draw();

  // cout << "(" << Ylast << " - " << Yfirst << ") = " << Ylast - Yfirst
  //      << " / "
  //      << "(" << Xlast << " - " << Xfirst << ") = " << Xlast - Xfirst
  //      << " = " << (Ylast-Yfirst) / (Xlast-Xfirst)
  //      << endl;
  
  double angle = (lMassHierarchy->GetY2() - lMassHierarchy->GetY1()) / (lMassHierarchy->GetX2() - lMassHierarchy->GetX1());
  angle = TMath::ATan(angle) / TMath::Pi() * 180.;  
  //cout << "angle = " << angle << endl;    

  // TLatex *latNormal = new TLatex(Xfirst,Yfirst,"Normal");
  // latNormal->SetLineColor(kAzure);
  // latNormal->SetTextSize(0.05);
  // latNormal->SetTextAngle(angle);
  // latNormal->Draw("same");
  // TLatex *latInverted = new TLatex(Xfirst,Yfirst,"Inverted");
  // latInverted->SetLineColor(kAzure);
  // latInverted->SetTextSize(0.05);
  // latInverted->SetTextAngle(angle);
  // latInverted->Draw("same");

  TLatex *latOct1 = new TLatex(Xfirst,Yfirst,"sin^{2} #theta_{23} = 0.45");
  latOct1->SetLineColor(kBlack);
  latOct1->SetTextSize(0.05);
  //latOct1->SetTextAngle(angle);
  latOct1->Draw("same");
  TLatex *latOct2 = new TLatex(Xfirst,Yfirst,"sin^{2} #theta_{23} = 0.55");
  latOct2->SetLineColor(kBlack);
  latOct2->SetTextSize(0.05);
  //latOct2->SetTextAngle(angle);
  latOct2->Draw("same");

  TLegend *ldcp = new TLegend(0.6,0.45,0.8,0.65);
  ldcp->SetFillColor(0);
  ldcp->SetLineColor(0);
  ldcp->SetTextFont(132);
  ldcp->SetTextSize(0.05);
  ldcp->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p55,"Normal Hierarchy","l");
  ldcp->AddEntry(gN_fhcVsrhc_elike_ih_sinsqth23_0p55,"Inverted Hierarchy","l");
  ldcp->Draw();

  TLegend *lFhcVsRhc = new TLegend(0.6,0.45,0.8,0.65);
  lFhcVsRhc->SetFillColor(0);
  lFhcVsRhc->SetLineColor(0);
  lFhcVsRhc->SetTextFont(132);
  lFhcVsRhc->SetTextSize(0.05);
  lFhcVsRhc->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi,"#delta_{CP} = #pi","p");
  lFhcVsRhc->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2,"#delta_{CP} = +#pi/2","p");
  lFhcVsRhc->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0,"#delta_{CP} = 0","p");
  lFhcVsRhc->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2,"#delta_{CP} = -#pi/2","p");
  lFhcVsRhc->Draw();

  TLegend *ldata = new TLegend(0.6,0.45,0.8,0.65);
  ldata->SetFillColor(0);
  ldata->SetLineColor(0);
  ldata->SetTextFont(132);
  ldata->SetTextSize(0.05);
  ldata->AddEntry(gN_fhcVsrhc_elike_data,str_preliminary.c_str(),"pe");
  ldata->Draw();

  
  // // Get total # of events
  //double TotEvts = hErecoVsQ2Tot->Integral();
  // cout << "Total = " << TotEvts << endl;
  // if(draw_data){
  //   cout << "Draw data spectra" << endl;
  // }
  // else {
  //   cout << endl;
  //   cout << "Don't draw data spectra" << endl;
  //   if(!draw_data_spectra){
  //     cout << "Chosen not to draw data spectra" << endl;
  //   }
  //   else{
  //     cout << "Input file " << rootfile << " doesn't exist" << endl;
  //     cout << endl;
  //   }
  // }

   
  



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


void GetN_fhcVsrhc(PhysHypothesis ph,SuperkExpectation &skexp,double &Nfhc,double &Nrhc){

  // 1Re # of events 
  conv::SuperkSample_t s = util::SKSampleFromString("1Re"); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
  Nfhc = hErecoVsQ2Tot->Integral();
  
  // RHC1Re # of events 
  s = util::SKSampleFromString("RHC1Re"); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
  Nrhc = hErecoVsQ2Tot->Integral();

  return;
}
