//
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Build.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

void generate_toymc_ensemble_joint_bar
(
 const char * tag                     = "toymc_nuebarapp_NuOscOnly_Toy23dcp_MH_betabar_0_rethrowunphysical_4p5E20POT",   

 int          run                     = 0,                   
 double       nexperiments            = 5000,               /* number of toy-MC expts to generate */
 const char * hypothesis              = "3f_osc_normbar", /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_osc_normbar, 3f_oscbar_normbar" */

 // double       input_mixing12_nh       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 // double       input_mixing13_nh       = -99,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 // double       input_mixing23_nh       = -99,   //conv::kT2KReacBestFitSinsqTheta23_NH, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 // double       input_deltacp_nh        = -99,   //conv::kT2KReacBestFitDeltaCP_NH, //-99,    /* deltacp, */
 // double       input_dm12sq_nh         = -99,    /* Dm^2_{12},  */
 // double       input_dm23sq_nh         = -99,   //conv::kT2KReacBestFitDmsq23_NH, //-99,    /* |Dm^2_{23}| */
 // double       input_mixing12_ih       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 // double       input_mixing13_ih       = -99,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 // double       input_mixing23_ih       = -99,   //conv::kT2KReacBestFitSinsqTheta23_IH, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 // double       input_deltacp_ih        = -99,   //conv::kT2KReacBestFitDeltaCP_IH, //-99,    /* deltacp, if -999 then randomize */
 // double       input_dm12sq_ih         = -99,    /* Dm^2_{12},  */
 // double       input_dm23sq_ih         = -99,   //conv::kT2KReacBestFitDmsq23_IH, //-99,    /* |Dm^2_{23}|, if -999 then randomize */
 // double       input_beta_bar          = -99,

 double       input_mixing12_nh       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13_nh       = conv::kT2KReacBestFitSinsqTheta13_NH,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23_nh       = -999,   //conv::kT2KReacBestFitSinsqTheta23_NH, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_deltacp_nh        = -999,   //conv::kT2KReacBestFitDeltaCP_NH, //-99,    /* deltacp, */
 double       input_dm12sq_nh         = -99,    /* Dm^2_{12},  */
 double       input_dm23sq_nh         = -999,   //conv::kT2KReacBestFitDmsq23_NH, //-99,    /* |Dm^2_{23}| */
 double       input_mixing12_ih       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13_ih       = conv::kT2KReacBestFitSinsqTheta13_IH,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23_ih       = -999,   //conv::kT2KReacBestFitSinsqTheta23_IH, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_deltacp_ih        = -999,   //conv::kT2KReacBestFitDeltaCP_IH, //-99,    /* deltacp, if -999 then randomize */
 double       input_dm12sq_ih         = -99,    /* Dm^2_{12},  */
 double       input_dm23sq_ih         = -999,   //conv::kT2KReacBestFitDmsq23_IH, //-99,    /* |Dm^2_{23}|, if -999 then randomize */
 double       input_beta_bar          = 0,

 const char * config_file      = "analysis_configuration_run1to6_test_app.xml",
 //const char * config_file      = "analysis_configuration_run1to6_test_app_fullstat.xml",

 bool rethrow_if_parameters_unphysical = true,                          // rethrow unphysical toys
 bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
 int          tweak_systematics       = 1,                              /* generate toy MC with randomized systematics? Numbers greater than 1 are special options */
 
 const bool   doMargMH                = true,      // marginalize Mass Hierarchy
 const double PriorProbNH             = 0.648861, // 0.5  // Prior for MH --> give the prob. of NH as input

 bool         save_in_xml_format      = true,     /* also save in txt format */
 bool         save_in_root_format     = true,                           /* */  
 const char * points_to_add_1Rmu      = 0,        /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_1Re       = 0,        /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_1Rmu_RHC  = 0,        /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_1Re_RHC   = 0,        /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_MultiRmu  = 0,        /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_MultiRe   = 0,        /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 //const char * xml_file_to_add         = "data_newstruct.xml",           /* xml with data in for when creating hybrid data/toy-MC experiments, eg data.xml. Always uses set 0 */
 const char * xml_file_to_add         = 0,        /* xml with data in for when creating hybrid data/toy-MC experiments, eg data.xml. Always uses set 0 */
 bool         save_fluct_in_xml       = true,     /* save the values of the variations in the xml file? */
 
 const char * dmsqscheme              = "auto",          /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy               = "normal",        /* normal or inverted */
 const char * mixparamchoice          = "sinsqtheta",    /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant          = "two",           /* one==kNuOctantOne two==kNuOctantTwo */

 TString file_priors = "../../inputs/priors/joint_T2Kreact_run14/valor/Priors_T2Kreact_run14_Valor.root" // Set to "" if no user priors are used (VALOR) 
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
  TH1D *hToy_dcp_nh;
  TH1D *hToy_dcp_ih;
  TH2D *hToy_23_nh;
  TH2D *hToy_23_ih;
  TH1D *hToy_MH = new TH1D("hToy_MH","0-->NH, 1-->IH",2,0,2);

  int n_points_per_line_xml = 5;
  
  // Use seed 0 for histo->GetRandom(), otherwise it extracts always the same distribution
  gRandom = new TRandom3(0);

  TRandom3 rand(0);

  TH1D * h_0 = new TH1D ("h_0", "h_0", 100, 0, 10);
  TH1D * h_1 = new TH1D ("h_1", "h_1", 100, 0, 10);
  TH1D * h_2 = new TH1D ("h_2","h_2", 100, 0, 10);
  TH1D * h_3 = new TH1D ("h_3","h_3", 100, 0, 10);
  TH1D * h_4 = new TH1D ("h_4","h_4", 100, 0, 10);
  TH2D * h_4_Escale = new TH2D ("h_4_Escale","h_4_Escale", 100, -5, 5, 100, 0, 10);

  // Analysis configuration
  const char * config = (new string(Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file)))->c_str();

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

  // Check whether some osc. params need to be randomized
  
  bool randomize_mixing12_nh   = (TMath::Abs(input_mixing12_nh + 999) < kASmallNum);
  bool randomize_mixing13_nh   = (TMath::Abs(input_mixing13_nh + 999) < kASmallNum);
  bool randomize_mixing23_nh   = (TMath::Abs(input_mixing23_nh + 999) < kASmallNum);
  bool randomize_deltacp_nh    = (TMath::Abs(input_deltacp_nh  + 999) < kASmallNum);
  bool randomize_dm12sq_nh     = (TMath::Abs(input_dm12sq_nh   + 999) < kASmallNum);
  bool randomize_dm23sq_nh     = (TMath::Abs(input_dm23sq_nh   + 999) < kASmallNum);

  bool randomize_mixing12_ih  = (TMath::Abs(input_mixing12_ih  + 999) < kASmallNum);
  bool randomize_mixing13_ih  = (TMath::Abs(input_mixing13_ih  + 999) < kASmallNum);
  bool randomize_mixing23_ih  = (TMath::Abs(input_mixing23_ih  + 999) < kASmallNum);
  bool randomize_deltacp_ih   = (TMath::Abs(input_deltacp_ih   + 999) < kASmallNum);
  bool randomize_dm12sq_ih    = (TMath::Abs(input_dm12sq_ih    + 999) < kASmallNum);
  bool randomize_dm23sq_ih    = (TMath::Abs(input_dm23sq_ih    + 999) < kASmallNum);

  bool randomize_beta_bar      = (TMath::Abs(input_beta_bar    + 999) < kASmallNum);

  //Set to defaults
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13_nh < -1 && input_mixing13_nh > -900) input_mixing13_nh = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_nh < -1 && input_mixing23_nh > -900) input_mixing23_nh = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_nh < -1 && input_mixing12_nh > -900) input_mixing12_nh = conv::kGlobFitSinsq2Theta12;

    if(input_mixing13_ih < -1 && input_mixing13_ih > -900) input_mixing13_ih = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_ih < -1 && input_mixing23_ih > -900) input_mixing23_ih = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_ih < -1 && input_mixing12_ih > -900) input_mixing12_ih = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13_nh < -1 && input_mixing13_nh > -900) input_mixing13_nh = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_nh < -1 && input_mixing23_nh > -900) input_mixing23_nh = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_nh < -1 && input_mixing12_nh > -900) input_mixing12_nh = conv::kGlobFitSinsqTheta12;
    
    if(input_mixing13_ih < -1 && input_mixing13_ih > -900) input_mixing13_ih = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_ih < -1 && input_mixing23_ih > -900) input_mixing23_ih = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_ih < -1 && input_mixing12_ih > -900) input_mixing12_ih = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq_nh < -1 && input_dm23sq_nh > -900)  input_dm23sq_nh = conv::kGlobFitDm23sq;
  if(input_dm12sq_nh < -1 && input_dm12sq_nh > -900)  input_dm12sq_nh = conv::kGlobFitDm12sq;
  if(input_deltacp_nh < -90 && input_deltacp_nh > -900) input_deltacp_nh = conv::kGlobFitDeltaCP;

  if(input_dm23sq_ih < -1 && input_dm23sq_ih > -900)  input_dm23sq_ih = conv::kGlobFitDm23sq;
  if(input_dm12sq_ih < -1 && input_dm12sq_ih > -900)  input_dm12sq_ih = conv::kGlobFitDm12sq;
  if(input_deltacp_ih < -90 && input_deltacp_ih > -900) input_deltacp_ih = conv::kGlobFitDeltaCP;
  
  if(input_beta_bar < -90 && input_beta_bar > -900) input_beta_bar = conv::kGlobFitBetaBar;


  // Take prior distributions
 
  TString name_prior_dcp_nh = "hPrior_Ldcp_nh";
  TString name_prior_dcp_ih = "hPrior_Ldcp_ih";
  TString name_prior_23_nh = "hPrior_L23_nh";
  TString name_prior_23_ih = "hPrior_L23_ih";
 
  TH1D *hPrior_dcp_nh = 0;
  TH2D *hPrior_23_nh = 0;
  TH1D *hPrior_dcp_ih = 0;
  TH2D *hPrior_23_ih = 0;
 
  bool UsePrior_dcp = false;
  bool UsePrior_23 = false;

  if(file_priors!=""){

    UsePrior_dcp = true;
    UsePrior_23 = true;
    
    cout << endl;
    cout << "Take priors from: " << file_priors << endl;
    TFile *fPrior = new TFile(file_priors.Data(),"READ");
    if(!fPrior->IsOpen() || !fPrior){
      cerr << endl;
      cerr << "Error: the file " << file_priors << " is not open!" << endl;
      cerr << endl;
      exit(1);
    }

    hPrior_dcp_nh = (TH1D*) fPrior->Get(name_prior_dcp_nh.Data());        
    if(hPrior_dcp_nh->IsZombie() || !hPrior_dcp_nh){
      cerr << endl;
      cerr << "Error: the histogram deltaCP normal hierarchy is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
    hPrior_dcp_ih = (TH1D*) fPrior->Get(name_prior_dcp_ih.Data());       
    if(hPrior_dcp_ih->IsZombie() || !hPrior_dcp_ih){
      cerr << endl;
      cerr << "Error: the histogram deltaCP inverted hierarchy is not taken!" << endl;
      cerr << endl;
      exit(1);
    }

    hPrior_23_nh = (TH2D*) fPrior->Get(name_prior_23_nh.Data()); 
    if(hPrior_23_nh->IsZombie() || !hPrior_23_nh){
      cerr << endl;
      cerr << "Error: the histogram par23_hist_normal is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
    hPrior_23_ih = (TH2D*) fPrior->Get(name_prior_23_ih.Data()); 
    if(hPrior_23_ih->IsZombie() || !hPrior_23_ih){
      cerr << endl;
      cerr << "Error: the histogram par23_hist_inverted is not taken!" << endl;
      cerr << endl;
      exit(1);
    }

    hToy_dcp_nh = (TH1D*) hPrior_dcp_nh->Clone("hToy_dcp_nh");
    hToy_dcp_nh->Reset();
    hToy_dcp_ih = (TH1D*) hPrior_dcp_ih->Clone("hToy_dcp_ih");
    hToy_dcp_ih->Reset();

    hToy_23_nh = (TH2D*) hPrior_23_nh->Clone("hToy_23_nh");
    hToy_23_nh->Reset();
    hToy_23_ih = (TH2D*) hPrior_23_ih->Clone("hToy_23_ih");
    hToy_23_ih->Reset();
        
  } // if(file_priors)



  //

  string tweak_sys_str = "None";
  switch(tweak_systematics) {
    case(0):
      tweak_sys_str = "None"; break;
    case(1):
      tweak_sys_str = "All"; break;
    case(2):
      tweak_sys_str = "Special low energy numu flux increase"; break;
    case(3):
      tweak_sys_str = "Tweak all but SK and FSI"; break;
    case(4):
      tweak_sys_str = "NC1pi norm = 0.5"; break;
    case(5):
      tweak_sys_str = "High energy >1.5 GeV flux increase"; break;
    case(6):
      tweak_sys_str = "BANFF"; break;
    case(7):
      tweak_sys_str = "XSec"; break;
    case(8):
      tweak_sys_str = "SK+FSI"; break;
    default:
      tweak_sys_str = "None";
  }
  
  string * full_points_to_add_1Rmu     = 0;
  string * full_points_to_add_1Re      = 0;
  string * full_points_to_add_1Rmu_RHC = 0;
  string * full_points_to_add_1Re_RHC  = 0;
  string * full_points_to_add_MultiRmu = 0;
  string * full_points_to_add_MultiRe  = 0;
  
  // Are we adding events from an XML file?
  if (xml_file_to_add) {
    TXMLEngine xml;
    const char * xml_path = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xml_file_to_add);
    XMLDocPointer_t xmldoc = xml.ParseFile(xml_path);
    if (xmldoc==0) {
      cerr << "Can not parse XML file to add: " << xml_path << endl;
      exit(1);
    }
    XMLNodePointer_t root_node = xml.DocGetRootElement(xmldoc);
    
    string s_1Rmu = util::XML2String(&xml, root_node, "superk/observation/data_set_0/one_mu_ring/ereco_values");
    full_points_to_add_1Rmu = new string( points_to_add_1Rmu ? Form("%s,%s",points_to_add_1Rmu,s_1Rmu) : s_1Rmu);

    string s_1Re = util::XML2String(&xml, root_node, "superk/observation/data_set_0/one_e_ring/ereco_values");
    full_points_to_add_1Re = new string( points_to_add_1Re ? Form("%s,%s",points_to_add_1Re,s_1Re) : s_1Re);
    
    string s_1Rmu_RHC = util::XML2String(&xml, root_node, "superk/observation/data_set_0/one_mu_ring_rhc/ereco_values");
    full_points_to_add_1Rmu_RHC = new string( points_to_add_1Rmu_RHC ? Form("%s,%s",points_to_add_1Rmu_RHC,s_1Rmu_RHC) : s_1Rmu_RHC);
    
    string s_1Re_RHC = util::XML2String(&xml, root_node, "superk/observation/data_set_0/one_e_ring_rhc/ereco_values");
    full_points_to_add_1Re_RHC = new string( points_to_add_1Re_RHC ? Form("%s,%s",points_to_add_1Re_RHC,s_1Re_RHC) : s_1Re_RHC);

    string s_MultiRmu = util::XML2String(&xml, root_node, "superk/observation/data_set_0/multi_mu_ring/ereco_values");
    full_points_to_add_MultiRmu = new string( points_to_add_MultiRmu ? Form("%s,%s",points_to_add_MultiRmu,s_MultiRmu) : s_MultiRmu);
    
    string s_MultiRe = util::XML2String(&xml, root_node, "superk/observation/data_set_0/multi_e_ring/ereco_values");
    full_points_to_add_MultiRe = new string( points_to_add_MultiRe ? Form("%s,%s",points_to_add_MultiRe,s_MultiRe) : s_MultiRe);
  }
  else {
    if (points_to_add_1Rmu    ) full_points_to_add_1Rmu     = new string(points_to_add_1Rmu    );
    if (points_to_add_1Re     ) full_points_to_add_1Re      = new string(points_to_add_1Re     );
    if (points_to_add_1Rmu_RHC) full_points_to_add_1Rmu_RHC = new string(points_to_add_1Rmu_RHC);
    if (points_to_add_1Re_RHC ) full_points_to_add_1Re_RHC  = new string(points_to_add_1Re_RHC );
    if (points_to_add_MultiRmu) full_points_to_add_MultiRmu = new string(points_to_add_MultiRmu);
    if (points_to_add_MultiRe ) full_points_to_add_MultiRe  = new string(points_to_add_MultiRe );
  }

  // Print out job summary
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Input physics hypothesis: " << hypothesis << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl;
  cout << "   - Randomize systematics : "      << tweak_sys_str << endl;
  cout << "   - Rethrow unphysical : " << ((rethrow_if_parameters_unphysical) ? "YES " : "NO") << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ? "YES " : "NO") << endl;

  cout << "   - True oscillation parameters (normal hierarchy): " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << ((randomize_mixing12_nh) ? "* randomized *" : input_mixing12_nh) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << ((randomize_mixing13_nh) ? "* randomized *" : input_mixing13_nh) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << ((randomize_mixing23_nh) ? "* randomized *" : input_mixing23_nh) << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = " << ((randomize_deltacp_nh) ? "* randomized *" : input_deltacp_nh) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq_nh)  ? "* randomized *" : input_dm12sq_nh ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << ((randomize_dm23sq_nh)  ? "* randomized *" : input_dm23sq_nh ) << " eV^2/c^4" << endl;

  cout << "   - True oscillation parameters (inverted hierarchy): " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << ((randomize_mixing12_ih) ? "* randomized *" : input_mixing12_ih) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << ((randomize_mixing13_ih) ? "* randomized *" : input_mixing13_ih) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << ((randomize_mixing23_ih) ? "* randomized *" : input_mixing23_ih) << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = " << ((randomize_deltacp_ih) ? "* randomized *" : input_deltacp_ih) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq_ih)  ? "* randomized *" : input_dm12sq_ih ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << ((randomize_dm23sq_ih)  ? "* randomized *" : input_dm23sq_ih ) << " eV^2/c^4" << endl;
  
  cout << endl;
  cout << "Same input parameters are used for Nu and NuBar modes are used" << endl;
  cout << endl;

  cout << "   - Adding the following data-points in each toy-MC expt: " << endl;
  cout << "    For 1R mu-like: "     << ((full_points_to_add_1Rmu)     ? *full_points_to_add_1Rmu     : "*None*") << endl;
  cout << "    For 1R e-like: "      << ((full_points_to_add_1Re)      ? *full_points_to_add_1Re      : "*None*") << endl;
  cout << "    For 1R mu-like RHC: " << ((full_points_to_add_1Rmu_RHC) ? *full_points_to_add_1Rmu_RHC : "*None*") << endl;
  cout << "    For 1R e-like RHC: "  << ((full_points_to_add_1Re_RHC)  ? *full_points_to_add_1Re_RHC  : "*None*") << endl;
  cout << "    For MultiR mu-like: " << ((full_points_to_add_MultiRmu) ? *full_points_to_add_MultiRmu : "*None*") << endl;
  cout << "    For MultiR e-like: "  << ((full_points_to_add_MultiRe)  ? *full_points_to_add_MultiRe  : "*None*") << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing toy-MC expts in XML format: "   << ((save_in_xml_format) ?   "YES " : "NO") << endl;
  cout << "     - Storing toy-MC expts in ROOT format:  " << ((save_in_root_format) ?  "YES " : "NO") << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  
  // Open files
  TFile * f = 0;
  if(save_in_root_format) {
    ostringstream root_filename;
    root_filename << tag << "-" << run << ".root";
    cout << root_filename.str();
    f = new TFile(root_filename.str().c_str(), "recreate");
  }

  ofstream fxml;
  ofstream ftxt;
  if(save_in_xml_format) {
    ostringstream xml_filename;
    xml_filename << tag << "-" << run << ".xml";
    fxml.open(xml_filename.str().c_str());
    fxml << "<?xml version=\"1.0\"?> " << endl;
    fxml << "<analysis_inputs>       " << endl;
    fxml << "  <superk>              " << endl;
    fxml << "    <observation>       " << endl;
    ostringstream txt_filename;
    txt_filename << tag << "-" << run << ".txt";
    ftxt.open(txt_filename.str().c_str());
  }

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  int weird = 0;
  for(int iexperiment=0; iexperiment < nexperiments; iexperiment++) {

    cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;
    
    // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs = new SuperkObservation;
     PhysHypothesis ph(hypothesis);

     double input_mixing12_used = conv::kBadResult;  
     double input_mixing13_used = conv::kBadResult;    
     double input_mixing23_used = conv::kBadResult;   
     double input_deltacp_used  = conv::kBadResult; 
     double input_dm12sq_used   = conv::kBadResult;   
     double input_dm23sq_used   = conv::kBadResult;   
     double input_beta_bar_used = conv::kBadResult;

     // Extract Mass Hierarchy

     conv::NuMassHierarchy_t hier;
     conv::NuMassHierarchy_t hierBar;
     
     int iMH = util::RandomizeMH(PriorProbNH); // 0 --> NH, 1 --> IH
     
     hToy_MH->Fill(iMH+conv::kASmallNum);
     
     if( iMH==0 ){ // Normal Hierarchy

       Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
       hier    = kNmhNormal;
       hierBar = kNmhNormal;
       skexp->GetSuperkMC()->SetMassHierarchy(hier);
       skexp->GetSuperkMC()->SetMassHierarchyBar(hierBar);
       
       // Randomize specified osc params for current toy MC expt
       if (randomize_mixing13_nh) { input_mixing13_nh = util::RandomizeMixing13(); }
       if (randomize_mixing12_nh) { input_mixing12_nh = util::RandomizeMixing12(); }
       if (randomize_dm12sq_nh  ) { input_dm12sq_nh   = util::RandomizeDm12sq();   }     
       if (randomize_dm23sq_nh && randomize_mixing23_nh && UsePrior_23){
	 double dm=0;
	 double mix23=0;
	 hPrior_23_nh->GetRandom2(mix23,dm);      
	 input_mixing23_nh = mix23;
	 input_dm23sq_nh = dm;
       }
       else{
	 if (randomize_dm23sq_nh  ) { input_dm23sq_nh   = util::RandomizeDm23sq();   }
	 if (randomize_mixing23_nh) { input_mixing23_nh = util::RandomizeMixing23(); }
       }
       if (randomize_deltacp_nh ) { 
	 if(UsePrior_dcp) input_deltacp_nh  = util::RandomizeUser(hPrior_dcp_nh);      
	 else input_deltacp_nh  = util::RandomizeCPV();   
       }

       input_mixing12_used = input_mixing12_nh;  
       input_mixing13_used = input_mixing13_nh;    
       input_mixing23_used = input_mixing23_nh;   
       input_deltacp_used  = input_deltacp_nh; 
       input_dm12sq_used   = input_dm12sq_nh;   
       input_dm23sq_used   = input_dm23sq_nh;   
     }
     
     else if(iMH==1){ // Inverted Hierarchy
     
       Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
       hier    = kNmhInverted;
       hierBar = kNmhInverted;
       skexp->GetSuperkMC()->SetMassHierarchy(hier);
       skexp->GetSuperkMC()->SetMassHierarchyBar(hierBar);
       
       // Randomize specified osc params for current toy MC expt
       if (randomize_mixing13_ih) { input_mixing13_ih = util::RandomizeMixing13(); }
       if (randomize_mixing12_ih) { input_mixing12_ih = util::RandomizeMixing12(); }
       if (randomize_dm12sq_ih  ) { input_dm12sq_ih   = util::RandomizeDm12sq();   }     
       if (randomize_dm23sq_ih && randomize_mixing23_ih && UsePrior_23){
	 double dm=0;
	 double mix23=0;
	 hPrior_23_ih->GetRandom2(mix23,dm);      
	 input_mixing23_ih = mix23;
	 input_dm23sq_ih = dm;
	 //hToy_23_ih->Fill(mix23,dm);
       }
       else{
	 if (randomize_dm23sq_ih  ) { input_dm23sq_ih   = util::RandomizeDm23sq();   }
	 if (randomize_mixing23_ih) { input_mixing23_ih = util::RandomizeMixing23(); }
       }
       if (randomize_deltacp_ih ) { 
	 if(UsePrior_dcp) input_deltacp_ih  = util::RandomizeUser(hPrior_dcp_ih);      
	 else input_deltacp_ih  = util::RandomizeCPV();   
	 //hToy_dcp_ih->Fill(input_deltacp_ih);
       }

       input_mixing12_used = input_mixing12_ih;  
       input_mixing13_used = input_mixing13_ih;    
       input_mixing23_used = input_mixing23_ih;   
       input_deltacp_used  = input_deltacp_ih; 
       input_dm12sq_used   = input_dm12sq_ih;   
       input_dm23sq_used   = input_dm23sq_ih;   
     }
     else{
       cerr << endl;
       cerr << "Mass Hierarchy is not defined" << endl;
       cerr << "iMH = " << iMH << endl;
       cerr << endl;
       exit(1);
     }

     if (randomize_beta_bar    ) { 
       input_beta_bar_used     = util::RandomizeBetaBar(); 
     }
     else{
       input_beta_bar_used     = input_beta_bar; 
     }
     
     // Set physics hypothesis
     if(ph.Hypothesis() == kHypNuMuDisapp2F) {
       ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
     }
     else if(ph.Hypothesis() == kHypOsc3F) {  
       ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
       ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
       ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
       ph.SetCurrValue(kPhypCPV    ,   input_deltacp_used);
       ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
     }
     else if(ph.Hypothesis() == kHypOsc3FBar) {
       ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
       ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
       ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
       ph.SetCurrValue(kPhypCPV    ,   input_deltacp_used);
       ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
       ph.SetCurrValue(kPhypMixing12Bar,  input_mixing12_used);
       ph.SetCurrValue(kPhypMixing13Bar,  input_mixing13_used);
       ph.SetCurrValue(kPhypMixing23Bar,  input_mixing23_used);
       ph.SetCurrValue(kPhypCPVBar    ,   input_deltacp_used);
       ph.SetCurrValue(kPhypDm12sqBar,    input_dm12sq_used);
       ph.SetCurrValue(kPhypDm23sqBar,    input_dm23sq_used);
     }
     else if(ph.Hypothesis() == kHypOsc3FNormBar) {
       ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
       ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
       ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
       ph.SetCurrValue(kPhypCPV    ,   input_deltacp_used);
       ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
       ph.SetCurrValue(kPhypBetaBar,    input_beta_bar_used);
     }
     else if(ph.Hypothesis() == kHypOsc3FBarNormBar) {
       ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
       ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
       ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
       ph.SetCurrValue(kPhypCPV    ,   input_deltacp_used);
       ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
       ph.SetCurrValue(kPhypMixing12Bar,  input_mixing12_used);
       ph.SetCurrValue(kPhypMixing13Bar,  input_mixing13_used);
       ph.SetCurrValue(kPhypMixing23Bar,  input_mixing23_used);
       ph.SetCurrValue(kPhypCPVBar    ,   input_deltacp_used);
       ph.SetCurrValue(kPhypDm12sqBar,    input_dm12sq_used);
       ph.SetCurrValue(kPhypDm23sqBar,    input_dm23sq_used);
       
       ph.SetCurrValue(kPhypBetaBar,    input_beta_bar_used);
     }
     
     cout << "- Physics hypothesis set: " << ph << endl;
     
     double mix12 = ph.CurrValue(kPhypMixing12);
     double mix13 = ph.CurrValue(kPhypMixing13);
     double mix23 = ph.CurrValue(kPhypMixing23);
     double dm12 = ph.CurrValue(kPhypDm12sq);
     double dm23 = ph.CurrValue(kPhypDm23sq);
     double dcp = ph.CurrValue(kPhypCPV);

     if(iMH==0){
       hToy_dcp_nh->Fill(dcp);
       hToy_23_nh->Fill(mix23,dm23);
     }
     else if(iMH==1){
       hToy_dcp_ih->Fill(dcp);
       hToy_23_ih->Fill(mix23,dm23);
     }
  
     //

     bool applied_negevents_condition = true;
     bool unphysical_tweak = true;

     //while(applied_negevents_condition) {
     while(applied_negevents_condition || unphysical_tweak) {
  
       if(tweak_systematics == 1) {
         // tweak all sources of systematic err (gaussian err)
         cout << "- Randomizing all systematic params in pdf..." << endl;
         skexp->TweakAllSystematics();
         //         skexp->TweakSystematic("sk_energy_scale",0.0,true);
         skexp->PrintSystematicVariations();
       }
       else if(tweak_systematics == 2) {
         // special tweak, increasing numu flux parameters <700MeV by 5%
         cout << "- Increasing numu flux parameters <700MeV by 5% in pdf..." << endl;
         TVectorD twkv(skexp->GetBANFF()->BestFitValues());
         twkv[0] += 0.05;
         twkv[1] += 0.05;
         twkv[2] += 0.05;
         twkv[3] += 0.05;
         skexp->GetBANFF()->Tweak(true,&twkv);
         skexp->PrintSystematicVariations();
       }
       else if(tweak_systematics == 3) {
         cout << "- Randomizing all but SK systematic params in pdf..." << endl;
         skexp->TweakAllBANFFSystematics();
         skexp->TweakAllUncorrelatedXSecSystematics();
         skexp->PrintSystematicVariations();
       }
       else if(tweak_systematics == 4) {
         skexp->TweakSystematic("xsec_norm_nc1piC",-0.5);
         skexp->PrintSystematicVariations();
       }
       else if(tweak_systematics == 5) {
         cout << "- Increasing flux parameters >1500MeV by 15% in pdf..." << endl;
         TVectorD twkv(skexp->GetBANFF()->BestFitValues());
         twkv[6]  += 0.15;
         twkv[7]  += 0.15;
         twkv[8]  += 0.15;
         twkv[9]  += 0.15;
         twkv[10] += 0.15;
         twkv[14] += 0.15;
         twkv[15] += 0.15;
         skexp->GetBANFF()->Tweak(true,&twkv);
         skexp->PrintSystematicVariations();
       }
       else if(tweak_systematics == 6) {
         cout << "- Randomizing all BANFF systematic params in pdf..." << endl;
         skexp->TweakAllBANFFSystematics();
         skexp->PrintSystematicVariations();
       }
       else if(tweak_systematics == 7) {
         cout << "- Randomizing all uncorrelated cross section systematic params in pdf..." << endl;
         skexp->TweakAllUncorrelatedXSecSystematics();
         skexp->PrintSystematicVariations();
       }
       else if(tweak_systematics == 8) {
         cout << "- Randomizing all SK+FSI systematic params in pdf..." << endl;
         skexp->TweakAllSKDetectorFSISystematics();
         skexp->PrintSystematicVariations();
       }

       // Rethrow unphysical parameters
       PhysParamList true_nuisance_params = skexp->Tweaks();
       if (rethrow_if_parameters_unphysical) {
	 unphysical_tweak = !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
	 if (unphysical_tweak) {cout<<"Parameters are unphysical - rethrowing"<<std::endl;}
       }
       else {
	 unphysical_tweak = false;
       }
       
       applied_negevents_condition = skobs->CreateMockData(*skexp, ph, stat_fluctuations);
       cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
     }//applied_negevents_condition
  
     double escale=0.;
#ifdef __FHC_RUNNING__
     escale = skexp->GetSuperkDetectorFSI()->CurrentValue(kNuisSuperkEnergyScale);
#endif     
#ifdef __FHC_RUNNING__
     // Get 1-ring mu-like Ereco spectrum FHC     
     TH1D* mock_data_1Rmu = new TH1D("","", util::NumErecoBins(kSKsFCFV1RingMuLike), util::ErecoBinEdges(kSKsFCFV1RingMuLike));
 
     if(util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingMuLike)>0){
       
       for(int ib = 1; ib <= util::NumErecoBins(kSKsFCFV1RingMuLike); ib++) {
	 mock_data_1Rmu->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsFCFV1RingMuLike)[ib-1]);
	 double events = skobs->ErecoSpectrumArray(kSKsFCFV1RingMuLike)[ib-1];
	 if(ib == 1) h_0->Fill(events);
	 if(ib == 2) h_1->Fill(events);
	 if(ib == 3) h_2->Fill(events);
	 if(ib == 4) h_3->Fill(events);
	 if(ib == 5) h_4->Fill(events);
	 if(ib == 5) h_4_Escale->Fill( escale, events);
       }
       mock_data_1Rmu->SetDirectory(0);     
     }
     

     TH1D* mock_data_1Re = new TH1D("","", util::NumErecoBins(kSKsFCFV1RingELike), util::ErecoBinEdges(kSKsFCFV1RingELike));
     
     if(util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingELike)>0){
       
       for(int ib = 1; ib <= util::NumErecoBins(kSKsFCFV1RingELike); ib++) {
	 mock_data_1Re->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsFCFV1RingELike)[ib-1]);
       }
       mock_data_1Re->SetDirectory(0);
     }
#endif

#ifdef __RHC_RUNNING__

     // Get 1-ring mu-like Ereco spectrum RHC
     TH1D* mock_data_1Rmu_RHC = new TH1D("","", util::NumErecoBins(kSKsRHCFCFV1RingMuLike), util::ErecoBinEdges(kSKsRHCFCFV1RingMuLike));
     
     if(util::ArrIdx::Instance()->SuperkSampleIdx(kSKsRHCFCFV1RingMuLike)>0){
       
       for(int ib = 1; ib <= util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++) {
	 mock_data_1Rmu_RHC->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib-1]);
       }
       mock_data_1Rmu_RHC->SetDirectory(0);     
     }
     
     // Get 1-ring e-like Ereco spectrum RHC
     TH1D* mock_data_1Re_RHC = new TH1D("","", util::NumErecoBins(kSKsRHCFCFV1RingELike), util::ErecoBinEdges(kSKsRHCFCFV1RingELike));

     if(util::ArrIdx::Instance()->SuperkSampleIdx(kSKsRHCFCFV1RingELike)>0){

       for(int ib = 1; ib <= util::NumErecoBins(kSKsRHCFCFV1RingELike); ib++) {
	 mock_data_1Re_RHC->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingELike)[ib-1]);
       }
       mock_data_1Re_RHC->SetDirectory(0);
     }

#endif


#ifdef __MULTIRING_FHC_RUNNING__
     // Get Multi-ring mu-like Ereco spectrum
     TH1D* mock_data_MultiRmu = new TH1D("","", util::NumErecoBins(kSKsFCFVMultiRingMuLike), util::ErecoBinEdges(kSKsFCFVMultiRingMuLike));

     if(util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFVMultiRingMuLike)>0){

       for(int ib = 1; ib <= util::NumErecoBins(kSKsFCFVMultiRingMuLike); ib++) {
	 mock_data_MultiRmu->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsFCFVMultiRingMuLike)[ib-1]);
       }
       mock_data_MultiRmu->SetDirectory(0);     
     }
     
     // Get Multi-ring e-like Ereco spectrum 
     TH1D* mock_data_MultiRe = new TH1D("","", util::NumErecoBins(kSKsFCFVMultiRingELike), util::ErecoBinEdges(kSKsFCFVMultiRingELike));

     if(util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFVMultiRingELike)>0){
       
       for(int ib = 1; ib <= util::NumErecoBins(kSKsFCFVMultiRingELike); ib++) {
	 mock_data_MultiRe->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsFCFVMultiRingELike)[ib-1]);
       }
       mock_data_MultiRe->SetDirectory(0);
     }
#endif
     

     
     // If an XML file was specified, add an entry
     
     if(save_in_xml_format) {

       // Dataset info
       fxml << "      " << Form("<data_set_%d>",iexperiment) << endl;
       if(save_fluct_in_xml) {
         fxml << "<!-- ";
         ph.Print(fxml);
         if(tweak_systematics)
           skexp->PrintSystematicVariations(fxml);
         fxml << " -->" << endl;
       }
     
       // FHC
#ifdef __FHC_RUNNING__
       // Now the 1Rmu part
       fxml << "        <one_mu_ring>    " << endl;
       fxml << "          <ereco_values> " << endl;
       ftxt << "1Rmu:";
       ftxt << endl;
       // Adding Ereco values as a comma separated list.
       // Assuming uniform distribution within each bin, randomize values within each bin.
       int ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_1Rmu->GetNbinsX(); ibin++) {
         double Nb   = mock_data_1Rmu->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_1Rmu->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_1Rmu->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
	 double Ecent = mock_data_1Rmu->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
         for(int iev = 0; iev < Nb; iev++) {
	   double E= Emin + rand.Uniform() * (Emax-Emin);
	   if(ievtot > 0) {
	     fxml << ", "; 
	     ftxt << ","; 
	     if ( ! (ievtot % n_points_per_line_xml) ) fxml << endl; 
	   }
	   fxml << E; 
	   ftxt << E; 
	   ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(full_points_to_add_1Rmu) {
         fxml << ",\n"<<endl; // should really have an XML comment here but TXMLEngine won't let us.
         fxml << *full_points_to_add_1Rmu;
         ftxt << "," << *full_points_to_add_1Rmu;
       }
       fxml << endl;
       ftxt << endl;
       fxml << "          </ereco_values> " << endl;
       fxml << "        </one_mu_ring>    " << endl;       
       
       //Now the 1Re part
       ftxt << endl;
       ftxt << "1Re:";
       ftxt << endl;
       fxml << "        <one_e_ring>    " << endl;
       fxml << "          <ereco_values> " << endl;       
       ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_1Re->GetNbinsX(); ibin++) {
         double Nb   = mock_data_1Re->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_1Re->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_1Re->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
 	 double Ecent = mock_data_1Re->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
        for(int iev = 0; iev < Nb; iev++) {
	   double E= Emin + rand.Uniform() * (Emax-Emin);
	   if(ievtot > 0) {
	     fxml << ", "; 
	     ftxt << ","; 
	     if ( ! (ievtot % n_points_per_line_xml) ) fxml << endl; 
	   }
	   fxml << E; 
	   ftxt << E; 
	   ievtot++;
	}//iev
       }//ibin
       // add extra data-points, if any
       if(full_points_to_add_1Re) {
         fxml << ",\n"<<endl;
         fxml << *full_points_to_add_1Re;
         ftxt << "," << *full_points_to_add_1Re;
       }
       fxml << endl;
       ftxt << endl;
       fxml << "          </ereco_values> " << endl;
       fxml << "        </one_e_ring>    " << endl;       
#endif
#ifdef __RHC_RUNNING__
       // Now the 1Rmu RHC part
       fxml << "        <one_mu_ring_rhc>    " << endl;
       fxml << "          <ereco_values> " << endl;       
       ftxt << "1Rmu RHC:";
       ftxt << endl;       
       // Adding Ereco values as a comma separated list.
       // Assuming uniform distribution within each bin, randomize values within each bin.
       int ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_1Rmu_RHC->GetNbinsX(); ibin++) {
         double Nb   = mock_data_1Rmu_RHC->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_1Rmu_RHC->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_1Rmu_RHC->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
	 double Ecent = mock_data_1Rmu_RHC->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
         for(int iev = 0; iev < Nb; iev++) {
	   double E= Emin + rand.Uniform() * (Emax-Emin);
            if(ievtot > 0) {
              fxml << ", "; 
              ftxt << ","; 
              if ( ! (ievtot % n_points_per_line_xml) ) fxml << endl; 
            }
            fxml << E; 
            ftxt << E; 
            ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(full_points_to_add_1Rmu_RHC) {
         fxml << ",\n"<<endl;
         fxml << *full_points_to_add_1Rmu_RHC;
         ftxt << "," << *full_points_to_add_1Rmu_RHC;
       }
       fxml << endl;
       ftxt << endl;
       fxml << "          </ereco_values> " << endl;
       fxml << "        </one_mu_ring_rhc>" << endl;
       
       //Now the 1Re RHC part
       ftxt << endl;
       ftxt << "1Re RHC:";
       ftxt << endl;
       fxml << "        <one_e_ring_rhc>"  << endl;
       fxml << "          <ereco_values> " << endl;       
       ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_1Re_RHC->GetNbinsX(); ibin++) {
         double Nb   = mock_data_1Re_RHC->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_1Re_RHC->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_1Re_RHC->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
 	 double Ecent = mock_data_1Re_RHC->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
	 for(int iev = 0; iev < Nb; iev++) {
	   double E= Emin + rand.Uniform() * (Emax-Emin);
	   if(ievtot > 0) {
	     fxml << ", "; 
	     ftxt << ","; 
              if ( ! (ievtot % n_points_per_line_xml) ) fxml << endl; 
            }
            fxml << E; 
            ftxt << E; 
            ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(full_points_to_add_1Re_RHC) {
         fxml << ",\n"<<endl;
         fxml << *full_points_to_add_1Re_RHC;
         ftxt << "," << *full_points_to_add_1Re_RHC;
       }
       fxml << endl;
       ftxt << endl;
       fxml << "          </ereco_values> " << endl;
       fxml << "        </one_e_ring_rhc>    " << endl;
#endif
#ifdef __MULTIRING_FHC_RUNNING__       
       // Now the MultiRmu part
       fxml << "        <multi_mu_ring>    " << endl;
       fxml << "          <ereco_values> " << endl;       
       ftxt << "MultiRmu :";
       ftxt << endl;       
       // Adding Ereco values as a comma separated list.
       // Assuming uniform distribution within each bin, randomize values within each bin.
       int ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_MultiRmu->GetNbinsX(); ibin++) {
         double Nb   = mock_data_MultiRmu->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_MultiRmu->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_MultiRmu->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
 	 double Ecent = mock_data_MultiRmu->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
	 for(int iev = 0; iev < Nb; iev++) {
 	   double E= Emin + rand.Uniform() * (Emax-Emin);
           if(ievtot > 0) {
              fxml << ", "; 
              ftxt << ","; 
              if ( ! (ievtot % n_points_per_line_xml) ) fxml << endl; 
            }
            fxml << E; 
            ftxt << E; 
            ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(full_points_to_add_MultiRmu) {
         fxml << ",\n"<<endl;
         fxml << *full_points_to_add_MultiRmu;
         ftxt << "," << *full_points_to_add_MultiRmu;
       }
       fxml << endl;
       ftxt << endl;
       fxml << "          </ereco_values> " << endl;
       fxml << "        </multi_mu_ring>" << endl;
       
       //Now the MultiRe part
       ftxt << endl;
       ftxt << "MultiRe:";
       ftxt << endl;
       fxml << "        <multi_e_ring>"  << endl;
       fxml << "          <ereco_values> " << endl;       
       ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_MultiRe->GetNbinsX(); ibin++) {
         double Nb   = mock_data_MultiRe->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_MultiRe->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_MultiRe->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
	 double Ecent = mock_data_MultiRe->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
         for(int iev = 0; iev < Nb; iev++) {
	   double E= Emin + rand.Uniform() * (Emax-Emin);
            if(ievtot > 0) {
              fxml << ", "; 
              ftxt << ","; 
              if ( ! (ievtot % n_points_per_line_xml) ) fxml << endl; 
            }
            fxml << E; 
            ftxt << E; 
            ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(full_points_to_add_MultiRe) {
         fxml << ",\n"<<endl;
         fxml << *full_points_to_add_MultiRe;
         ftxt << "," << *full_points_to_add_MultiRe;
       }
       fxml << endl;
       ftxt << endl;
       fxml << "          </ereco_values> " << endl;
       fxml << "        </multi_e_ring>    " << endl;
#endif
       fxml << "      " << Form("</data_set_%d>\n",iexperiment) << endl;       
     }//xml

     // If a ROOT filename was specified, save the ROOT histogram too
     if(f) {
       if(save_in_root_format) {
#ifdef __FHC_RUNNING__
           if(full_points_to_add_1Rmu) {
            vector<string> vecp = util::Tokens(*full_points_to_add_1Rmu,",");
            assert(vecp.size()>0);
            for(int ip=0; ip < vecp.size(); ip++) {
              mock_data_1Rmu->Fill(atof(vecp[ip].c_str()));
            }
          }
          if(full_points_to_add_1Re) {
            vector<string> vecp = util::Tokens(*full_points_to_add_1Re,",");
            assert(vecp.size()>0);
            for(int ip=0; ip < vecp.size(); ip++) {
              mock_data_1Re->Fill(atof(vecp[ip].c_str()));
            }
          }
#endif
#ifdef __RHC_RUNNING__	 
          if(full_points_to_add_1Rmu_RHC) {
            vector<string> vecp = util::Tokens(*full_points_to_add_1Rmu_RHC,",");
            assert(vecp.size()>0);
            for(int ip=0; ip < vecp.size(); ip++) {
              mock_data_1Rmu_RHC->Fill(atof(vecp[ip].c_str()));
            }
          }
          if(full_points_to_add_1Re_RHC) {
            vector<string> vecp = util::Tokens(*full_points_to_add_1Re_RHC,",");
            assert(vecp.size()>0);
            for(int ip=0; ip < vecp.size(); ip++) {
              mock_data_1Re_RHC->Fill(atof(vecp[ip].c_str()));
            }
          }
#endif
#ifdef __MULTIRING_FHC_RUNNING__
          if(full_points_to_add_MultiRmu) {
            vector<string> vecp = util::Tokens(*full_points_to_add_MultiRmu,",");
            assert(vecp.size()>0);
            for(int ip=0; ip < vecp.size(); ip++) {
              mock_data_MultiRmu->Fill(atof(vecp[ip].c_str()));
            }
          }
          if(full_points_to_add_MultiRe) {
            vector<string> vecp = util::Tokens(*full_points_to_add_MultiRe,",");
            assert(vecp.size()>0);
            for(int ip=0; ip < vecp.size(); ip++) {
              mock_data_MultiRe->Fill(atof(vecp[ip].c_str()));
            }
          }
#endif
       }
       f->cd();
       
#ifdef __FHC_RUNNING__
       mock_data_1Rmu->Write(Form("mock_data_1Rmu_%d",iexperiment));
       mock_data_1Re->Write(Form("mock_data_1Re_%d",iexperiment));
#endif
#ifdef __RHC_RUNNING__
       mock_data_1Rmu_RHC->Write(Form("mock_data_1Rmu_RHC_%d",iexperiment));
       mock_data_1Re_RHC->Write(Form("mock_data_1Re_RHC_%d",iexperiment));
#endif
#ifdef __MULTIRING_FHC_RUNNING__
       mock_data_MultiRmu->Write(Form("mock_data_MultiRmu_%d",iexperiment));
       mock_data_MultiRe->Write(Form("mock_data_MultiRe_%d",iexperiment));
#endif

       // Object for storing input oscillation and nuisance params in each toy-MC
       PhysParamList * iparm = new PhysParamList;
       iparm->Clear();
       // keep track of oscillation parameters
       PhysParamList true_oscillation_params = ph.AsPhysParamList();
       iparm->Append(true_oscillation_params);
       // keep track of all systematic params
#ifdef __FHC_RUNNING__
       PhysParamList true_nuisance_params = skexp->Tweaks();
       iparm->Append(true_nuisance_params);
#endif  
#ifdef __RHC_RUNNING__
       PhysParamList true_nuisance_params_rhc = skexp->Tweaks();
       iparm->Append(true_nuisance_params_rhc);
#endif  
#ifdef __MULTIRING_FHC_RUNNING__
       PhysParamList true_nuisance_params_multiring = skexp->Tweaks();
       iparm->Append(true_nuisance_params_multiring);
#endif  
     iparm->Write(Form("iparm_%d",iexperiment));
     
     delete iparm;
     }
     
  }//iexperiment

  // Draw histograms

  TCanvas *cToy_dcp_nh = new TCanvas("cToy_dcp_nh","cToy_dcp_nh");
  hToy_dcp_nh->Draw();

  TCanvas *cToy_dcp_ih = new TCanvas("cToy_dcp_ih","cToy_dcp_ih");
  hToy_dcp_ih->Draw();

  TCanvas *cToy_23_nh = new TCanvas("cToy_23_nh","cToy_23_nh");
  hToy_23_nh->Draw("colz");

  TCanvas *cToy_23_ih = new TCanvas("cToy_23_ih","cToy_23_ih");
  hToy_23_ih->Draw("colz");

  TCanvas *cToy_MH = new TCanvas("cToy_MH","cToy_MH");
  hToy_MH->Draw();
  
  // Write histograms to output file

  h_0->Write();
  h_1->Write();
  h_2->Write();
  h_3->Write();
  h_4->Write();
  h_4_Escale->Write();
  
  hPrior_dcp_nh->Write();  
  hPrior_dcp_ih->Write();  
  hToy_dcp_nh->Write();
  hToy_dcp_ih->Write();

  hPrior_23_nh->Write();
  hPrior_23_ih->Write();
  hToy_23_nh->Write();
  hToy_23_ih->Write();

  hToy_MH->Write();

  //
  // Close files
  //
  if(save_in_xml_format) {
    fxml << "    </observation>       " << endl;
    fxml << "  </superk>              " << endl;
    fxml << "</analysis_inputs>       ";// << endl;
    fxml.close();
  }
  if(f) {
    f->Close();
    delete f;
  }


}

