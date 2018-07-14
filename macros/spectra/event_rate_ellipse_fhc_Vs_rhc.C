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
 bool include_cc1pi = false,
 const char * tag = "ellipse_Run18",
 const char* str_preliminary = "T2K Run 1-9c preliminary",
 const char * hypothesis = "3f_osc", /* 2f_numu_disap, 3f_osc, 3f_oscbar, 3f_osc_normbar, 3f_oscbar_normbar */ 
 bool         disable_banff_tuning = false,
 const char * config_file = "analysis_configuration_oa2018_run1-9c_v1_postfit.xml",
 const char * mixparamchoice = "sinsqtheta",
 const char * hierarchy = "normal",
 const char * dmsqscheme = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamoctant = "two", /* one==kNuOctantOne two==kNuOctantTwo */
 const char * rootfile = "official/2018/data/SK_data_Run19c_valor.root", 
 const bool draw_data_spectra = true
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

  PhysHypothesis ph(hypothesis);

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
  double dm23_nh = 2.440E-3;

  double ss12_ih = 0.304;
  double ss13_ih = 0.0226125;
  double dm12_ih = 7.53E-5;
  double dm23_ih = 2.410E-3;

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
    double ss23 = 0.450;
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
    double ss23 = 0.500;
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
    double ss23 = 0.550;
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
    double ss23 = 0.450;//0.496;
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
    double ss23 = 0.500;
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
    double ss23 = 0.550;//0.559;
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

    cout << endl;
  }

  // Get data spectra and rescale MC to data

  TGraphErrors *gN_fhcVsrhc_elike_data = new TGraphErrors();
  conv::SuperkSample_t s = conv::kSKsUndefined;
  if(draw_data){
    bool using_xml_file = false;
    skobs_fakedata->Read(rootfile2.Data(),0,using_xml_file);

    // 1Re # of events 
    s = util::SKSampleFromString("1Re"); // requested sample
    if(s == conv::kSKsUndefined) {
      cerr << "** Undefined sample!" << endl;
      return;
    }
    int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
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
  double xmin = 25.0;
  double ymin = 4.0;
  double ymax = 18.0;
  hInit->GetXaxis()->SetRangeUser(xmin,xmax);
  hInit->GetYaxis()->SetRangeUser(ymin,ymax);
  double ttx = (xmax + xmin) * 0.5;
  double ctx = xmin + (xmax + xmin) * 0.1;
  double ty = ymax;
  hInit->GetXaxis()->SetTitle("Neutrino mode 1Re candidates");
  hInit->GetYaxis()->SetTitle("Antineutrino mode 1Re candidates");

  // Draw 

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

  gN_fhcVsrhc_elike_data->SetMarkerColor(kBlack);
  gN_fhcVsrhc_elike_data->SetMarkerSize(1);
  gN_fhcVsrhc_elike_data->SetMarkerStyle(4);

  TLatex* title = new TLatex(ttx, 1.01 * ty, str_preliminary);
  title->SetTextSize(0.05);
  TLatex* caveat = new TLatex(ctx, 1.01 * ty, "Final systematics pending");
  caveat->SetTextSize(0.025);

  TCanvas *cEllipse_all = new TCanvas("cEllipse_all","cEllipse_all");
  hInit->Draw();

  gN_fhcVsrhc_elike_nh_sinsqth23_0p5->SetLineColor(kBlack);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5->SetLineStyle(1);
  gN_fhcVsrhc_elike_nh_sinsqth23_0p5->SetLineWidth(1);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->SetLineColor(kBlack);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->SetLineStyle(2);
  gN_fhcVsrhc_elike_ih_sinsqth23_0p5->SetLineWidth(1);
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

  TLegend *ldcp_all = new TLegend(0.15,0.19,0.35,0.64);
  ldcp_all->SetFillColor(0);
  ldcp_all->SetLineColor(0);
  ldcp_all->SetTextFont(132);
  ldcp_all->SetTextSize(0.035);
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5,"sin^{2} #theta_{23} = #color[1]{0.50}, #color[2]{0.45}, #color[4]{0.55}","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5,"#Deltam^{2}_{32} = 2.44#times10^{-3} eV^{2}/c^{4}","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_ih_sinsqth23_0p5,"#Deltam^{2}_{31} = -2.41#times10^{-3} eV^{2}/c^{4}","l");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi,"#delta_{CP} = #pi","p");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_ppi2,"#delta_{CP} = +#pi/2","p");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_0,"#delta_{CP} = 0","p");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_nh_sinsqth23_0p5_dcp_mpi2,"#delta_{CP} = -#pi/2","p");
  ldcp_all->AddEntry(gN_fhcVsrhc_elike_data,"Data (stat. errors only)","pe");
  ldcp_all->Draw();

  title->Draw();

  cEllipse_all->SaveAs(Form("%s_all.pdf", tag));
  cEllipse_all->SaveAs(Form("%s_all.png", tag));
  cEllipse_all->SaveAs(Form("%s_all.eps", tag));

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
