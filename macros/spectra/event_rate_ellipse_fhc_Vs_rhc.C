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
 bool include_cc1pi                = false,
 const char * tag                  = "ellipse_Run18",
 const string str_preliminary = "T2K Run 1-8 preliminary",
 const char * sk_sample            = "1Re",                         /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
 const char * hypothesis           = "3f_osc",                 /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_osc_normbar, 3f_oscbar_normbar" */ 
 // Best fit values NH
 double       input_mixing12       = 0.304,
 double       input_mixing13       = 0.02235,
 double       input_mixing23       = 0.530,
 double       input_deltacp        = 0,
 double       input_dm12sq         = 7.53E-5,
 double       input_dm23sq         = 2.460E-3,
 double       input_betabar        = -99., //0., //-99,
 bool         disable_banff_tuning = false,                          /* */
 const char * config_file = "analysis_configuration_oa2017_run1-8_postfit.xml",
 bool         store_spectrum      = true,
 bool         print_spectrum      = true,                          /* */
 const char * mixparamchoice      = "sinsqtheta",
 const char * hierarchy           = "normal"   ,
 const char * dmsqscheme          = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */
 bool plot_per_50MeV = false,
 double ErecoRangeLow = 0.,
 double ErecoRangeUp = 30.,
 int          iexperiment         = 0,
 const char * rootfile = "data_run1to8/SKdata_Run18.root", 
 const bool draw_data_spectra = true,
 const bool scale_mc = false,
 const bool change_banff_params = false
)
{
  using namespace std;
  using namespace valor;

  cout << "include_cc1pi = " << include_cc1pi << endl;
  cout << "Tag: " << tag << endl;

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

  double ss12_nh = 0.304;
  double ss13_nh = 0.02235;
  double dm12_nh = 7.53E-5;
  double dm23_nh = 2.460E-3;

  double ss12_ih = 0.304;
  double ss13_ih = 0.0226125;
  double dm12_ih = 7.53E-5;
  double dm23_ih = 2.440E-3;

  for(int idcp=0;idcp<dcp_npt;idcp++){

    double Nfhc = -999;
    double Nrhc = -999;
    double Ncc1pi = -999;
 
    double dcp = deltaCP_first + idcp * deltaCP_width;
    
    cout << "dcp = " << dcp << endl;
    // Set NH oscillation parameter values
    ph.SetCurrValue(conv::kPhypMixing12, ss12_nh);
    ph.SetCurrValue(conv::kPhypMixing13, ss13_nh);
    ph.SetCurrValue(conv::kPhypDm12sq, dm12_nh);
    ph.SetCurrValue(conv::kPhypDm23sq, dm23_nh);

    // 
    // sinsqth23 = -1sigma - Normal Hierarchy
    //
    double ss23 = 0.404;//0.493;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhNormal);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23,ss23);
    // Generate prediction
    skexp.GetSuperkMC()->SetX2p2h(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    Ncc1pi = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    GetN_cc1pi(ph,skexp,Ncc1pi);
    if(include_cc1pi)
    {
       Nfhc += Ncc1pi;
    }
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
    // sinsqth23 = 0.530 - Best fit Normal Hierarchy
    //
    double ss23 = 0.530;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhNormal);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    //Generate prediction
    skexp.GetSuperkMC()->SetX2p2h(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    Ncc1pi = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    GetN_cc1pi(ph,skexp,Ncc1pi);
    if(include_cc1pi)
    {
       Nfhc += Ncc1pi;
    }

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
    // sinsqth23 = +1sigma - Normal Hierarchy
    //
    double ss23 = 0.624;//0.561;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhNormal);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    // Generate prediction
    skexp.GetSuperkMC()->SetX2p2h(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    Ncc1pi = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    GetN_cc1pi(ph,skexp,Ncc1pi);
    if(include_cc1pi)
    {
       Nfhc += Ncc1pi;
    }

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

    // Set IH oscillation parameter values
    ph.SetCurrValue(conv::kPhypMixing12, ss12_ih);
    ph.SetCurrValue(conv::kPhypMixing13, ss13_ih);
    ph.SetCurrValue(conv::kPhypDm12sq, dm12_ih);
    ph.SetCurrValue(conv::kPhypDm23sq, dm23_ih);

    // 
    // sinsqth23 = -1sigma - Inverted Hierarchy
    //
    double ss23 = 0.398;//0.496;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhInverted);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    // Generate prediction
    skexp.GetSuperkMC()->SetX2p2h(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    Ncc1pi = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    GetN_cc1pi(ph,skexp,Ncc1pi);
    if(include_cc1pi)
    {
       Nfhc += Ncc1pi;
    }

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
    // sinsqth23 = Best fit - Inverted Hierarchy
    //
    double ss23 = 0.530;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhInverted);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    //Generate prediction
    skexp.GetSuperkMC()->SetX2p2h(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    Ncc1pi = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    GetN_cc1pi(ph,skexp,Ncc1pi);
    if(include_cc1pi)
    {
       Nfhc += Ncc1pi;
    }

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
    // sinsqth23 = +1sigma - Inverted Hierarchy
    //
    double ss23 = 0.618;//0.559;
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
    skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhInverted);
    ph.SetCurrValue(conv::kPhypCPV,dcp);
    ph.SetCurrValue(conv::kPhypMixing23, ss23);
    // Generate prediction
    skexp.GetSuperkMC()->SetX2p2h(1);
    skexp.Generate(ph);
    Nfhc = -999;
    Nrhc = -999;
    Ncc1pi = -999;
    GetN_fhcVsrhc(ph,skexp,Nfhc,Nrhc);
    GetN_cc1pi(ph,skexp,Ncc1pi);
    if(include_cc1pi)
    {
       Nfhc += Ncc1pi;
    }

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
    if(include_cc1pi)
    {
       s = util::SKSampleFromString("MultiRe"); // requested sample
       if(s == conv::kSKsUndefined)
       {
         cerr << "** Undefined sample!" << endl;
         return;
       }
       is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
       nfhc += skobs_fakedata->TwoDSpectrum(s)->Integral();
    }

    // RHC1Re # of events 
    s = util::SKSampleFromString("RHC1Re"); // requested sample
    if(s == conv::kSKsUndefined) {
      cerr << "** Undefined sample!" << endl;
      return;
    }
    is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
    double nrhc = skobs_fakedata->TwoDSpectrum(s)->Integral();
    
    double nfhc_err = sqrt(nfhc);
    double nrhc_err = sqrt(nrhc);

    cout << endl;
    cout << "Data: nfhc = " << nfhc << " +- " << nfhc_err << ", nrhc = " << nrhc << " +- " << nrhc_err << endl;
    cout << endl;

    int ipt = gN_fhcVsrhc_elike_data->GetN();
    gN_fhcVsrhc_elike_data->SetPoint(ipt,nfhc,nrhc);
    gN_fhcVsrhc_elike_data->SetPointError(ipt,nfhc_err,nrhc_err);

    cout << "idcp = " << idcp << "" << endl;
  }

  // Init histogram
  
  TH2D *hInit = new TH2D("hInit","hInit",500,0,100,500,0,100);
  double xmax = 100.0;
  double xmin = 20.0;
  double ymin = 0.0;
  double ymax = 13.0;
  hInit->GetXaxis()->SetRangeUser(xmin,xmax);
  hInit->GetYaxis()->SetRangeUser(ymin,ymax);
  double ttx = (xmax + xmin) * 0.5;
  double ctx = xmin + (xmax + xmin) * 0.1;
  double ty = ymax;
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

  gN_fhcVsrhc_elike_nh_sinsqth23_0p5->SetLineColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2->SetMarkerColor(kBlue); 
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->SetLineColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_mpi2->SetMarkerColor(kRed); 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_0->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2->SetMarkerColor(kRed);

  gN_fhcVsrhc_elike_ih_sinsqth23_0p45->SetLineStyle(1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55->SetLineStyle(1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->SetLineStyle(1);

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


  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi->GetPoint(0,x1,y1);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi->GetPoint(0,x2,y2);
  Xfirst = (x1+x2) / 2.; // * 0.7; 
  Yfirst = (y1+y2) / 2.; // * 0.7; 

  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0->GetPoint(0,x1,y1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0->GetPoint(0,x2,y2);
  Xlast = (x1+x2) / 2.; // * 1.2; 
  Ylast = (y1+y2) / 2.; // * 1.2;  

  double xtag_m1 = x1;
  double ytag_m1 = y1;
  double xtag_p1 = x2;
  double ytag_p1 = y2;

  if(xtag_m1 > xtag_p1)
  {
     xtag_m1 *= 1.05;
     xtag_p1 *= 0.9;
  }
  else
  {
     xtag_m1 *= 0.9;
     xtag_p1 *= 1.05;
  }
  if(ytag_m1 > ytag_p1)
  {
     ytag_m1 *= 1.05;
     ytag_p1 *= 0.9;
  }
  else
  {
     ytag_m1 *= 0.9;
     ytag_p1 *= 1.05;
  }

  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->GetPoint(0,x1,y1);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->GetPoint(0,x2,y2);
  double Xfirst = (x1+x2) / 2.; 
  double Yfirst = (y1+y2) / 2.; 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->GetPoint(0,x1,y1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->GetPoint(0,x2,y2);
  double Xlast = (x1+x2) / 2.; 
  double Ylast = (y1+y2) / 2.; 


  TLine *lOctant = new TLine(Xfirst,Yfirst,Xlast,Ylast);
  lOctant->SetLineColor(kBlack);
  lOctant->SetLineWidth(2);
  lOctant->SetLineStyle(2);
  lOctant->Draw();

  double angle = (lMassHierarchy->GetY2() - lMassHierarchy->GetY1()) / (lMassHierarchy->GetX2() - lMassHierarchy->GetX1());
  angle = TMath::ATan(angle) / TMath::Pi() * 180.;  

//  TLatex *latOct1 = new TLatex(xtag_m1, ytag_m1,"sin^{2} #theta_{23,-1#sigma}=0.493 (0.496)");
  TLatex *latOct1 = new TLatex(xtag_m1, ytag_m1,"sin^{2} #theta_{23,-1#sigma}=0.404 (0.398)");
  latOct1->SetLineColor(kBlack);
  latOct1->SetTextSize(0.035);
  latOct1->Draw("same");
  TLatex *latOct2 = new TLatex(xtag_p1, ytag_p1,"sin^{2} #theta_{23,+1#sigma}=0.624 (0.618)");
  latOct2->SetLineColor(kBlack);
  latOct2->SetTextSize(0.035);
  latOct2->Draw("same");

  TLegend *ldcp = new TLegend(0.2,0.35,0.4,0.6);
  ldcp->SetFillColor(0);
  ldcp->SetLineColor(0);
  ldcp->SetTextFont(132);
  ldcp->SetTextSize(0.035);
  ldcp->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p55,"Normal Hierarchy","l");
  ldcp->AddEntry(gN_fhcVsrhc_elike_ih_sinsqth23_0p55,"Inverted Hierarchy","l");
  ldcp->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi,"#delta_{CP} = #pi","p");
  ldcp->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2,"#delta_{CP} = +#pi/2","p");
  ldcp->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0,"#delta_{CP} = 0","p");
  ldcp->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2,"#delta_{CP} = -#pi/2","p");
  ldcp->AddEntry(gN_fhcVsrhc_elike_data,"Data","pe");
  ldcp->Draw();

  TLatex* title = new TLatex(ttx, 1.01 * ty, "T2K Run1-8 Preliminary");
  title->SetTextSize(0.05);
  title->Draw();
  TLatex* caveat = new TLatex(ctx, 1.01 * ty, "Final systematics pending");
  caveat->SetTextSize(0.025);
  caveat->Draw();

  cEllipse_nh_sinsqth23->SaveAs(Form("%s_1sigma.pdf", tag));

  TCanvas *cEllipse_bestfit = new TCanvas("cEllipse_bestfit","cEllipse_bestfit");
  hInit->Draw();

  gN_fhcVsrhc_elike_nh_sinsqth23_0p5->Draw("l");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->Draw("l");

  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2->Draw("p"); 
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2->Draw("p");

  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_mpi2->Draw("p"); 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_0->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2->Draw("p");

  gN_fhcVsrhc_elike_data->Draw("PE");

  TLegend *ldcp_bf = new TLegend(0.2,0.35,0.4,0.6);
  ldcp_bf->SetFillColor(0);
  ldcp_bf->SetLineColor(0);
  ldcp_bf->SetTextFont(132);
  ldcp_bf->SetTextSize(0.035);
  ldcp_bf->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5,"Normal Hierarchy","l");
  ldcp_bf->AddEntry(gN_fhcVsrhc_elike_ih_sinsqth23_0p5,"Inverted Hierarchy","l");
  ldcp_bf->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi,"#delta_{CP} = #pi","p");
  ldcp_bf->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2,"#delta_{CP} = +#pi/2","p");
  ldcp_bf->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0,"#delta_{CP} = 0","p");
  ldcp_bf->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2,"#delta_{CP} = -#pi/2","p");
  ldcp_bf->AddEntry(gN_fhcVsrhc_elike_data,"Data","pe");
  ldcp_bf->Draw();

  title->Draw();
  caveat->Draw();

  cEllipse_bestfit->SaveAs(Form("%s_bestfit.pdf", tag));

  TCanvas *cEllipse_all = new TCanvas("cEllipse_all","cEllipse_all");
  hInit->Draw();

  gN_fhcVsrhc_elike_nh_sinsqth23_0p5->SetLineColor(kBlack);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5->SetLineStyle(1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->SetLineColor(kBlack);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->SetLineStyle(2);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_mpi2->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_0->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2->SetMarkerColor(kBlack);


  gN_fhcVsrhc_elike_nh_sinsqth23_0p5->Draw("l");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->Draw("l");

  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2->Draw("p"); 
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0->Draw("p");
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2->Draw("p");

  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_mpi2->Draw("p"); 
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_0->Draw("p");
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5_dcp_ppi2->Draw("p");

  gN_fhcVsrhc_elike_nh_sinsqth23_0p45->SetLineColor(kRed);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45->SetLineStyle(1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45->SetLineColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45->SetLineStyle(2);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55->SetLineColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55->SetLineStyle(1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55->SetLineColor(kBlue);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55->SetLineStyle(2);

  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_mpi2->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_0->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p45_dcp_ppi2->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_mpi2->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_0->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p55_dcp_ppi2->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_mpi2->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_0->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p45_dcp_ppi2->SetMarkerColor(kRed);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_mpi2->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_0->SetMarkerColor(kBlue);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p55_dcp_ppi2->SetMarkerColor(kBlue);

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

  gN_fhcVsrhc_elike_data->Draw("PE");

  TLegend *ldcp_all = new TLegend(0.15,0.19,0.35,0.59);
  ldcp_all->SetFillColor(0);
  ldcp_all->SetLineColor(0);
  ldcp_all->SetTextFont(132);
  ldcp_all->SetTextSize(0.035);
//  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5,"sin^{2} #theta_{23,bf} = 0.530 (0.530)","l");
//  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p45,"sin^{2} #theta_{23,-1#sigma} = 0.493 (0.496)","l");
//  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p55,"sin^{2} #theta_{23,+1#sigma} = 0.561 (0.559)","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5,"sin^{2} #theta_{23,T2K} = 0.530 (0.530)","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p45,"sin^{2} #theta_{23,NoVA_{1}} = 0.404 (0.398)","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p55,"sin^{2} #theta_{23,NoVA_{2}} = 0.624 (0.618)","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5,"Normal Hierarchy","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_ih_sinsqth23_0p5,"Inverted Hierarchy","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi,"#delta_{CP} = #pi","p");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2,"#delta_{CP} = +#pi/2","p");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0,"#delta_{CP} = 0","p");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2,"#delta_{CP} = -#pi/2","p");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_data,"Data","pe");
  ldcp_all->Draw();

  title->Draw();
  caveat->Draw();

  cEllipse_all->SaveAs(Form("%s_nova.pdf", tag));

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

void GetN_cc1pi(PhysHypothesis ph,SuperkExpectation &skexp,double &Nfhc){

  // MultiRe # of events 
  conv::SuperkSample_t s = util::SKSampleFromString("MultiRe"); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
  Nfhc = hErecoVsQ2Tot->Integral();
  
  return;
}
