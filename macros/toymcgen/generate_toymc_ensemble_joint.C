//
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Build.h"

void generate_toymc_ensemble_joint
(
 const char * tag                     = "toymc_numudisap",                        /* */
 //const char * tag                     = "toymc_FHCpRHCpMultiRing",                        /* */
 //const char * tag                     = "toymc_FHCpRHC_DataRun1to4_MCRun5and6_ResultsRun1to4_IH",                        /* */
 //const char * tag                     = "toymc_FHCpRHC_DataRun1to4_MCRun5and6_Run1to4_IH",                        /* */
 //const char * tag                     = "toymc_FHCpRHC_MCRun1to6_th13DB_NH",                        /* */

 int          run                     = 0,                              /* */
 double       nexperiments            = 1,                             /* number of toy-MC expts to generate */
 const char * hypothesis              = "3f_oscbar",                       /* options: "2f_numu_disap", "3f_osc, 3f_oscbar" */

 double       input_mixing12          = -99, //0.306, //-999,                            /* 12-mixing, sinq2theta or sinsqtheta depending on config choice, if -999 then randomize  */
 double       input_mixing13          = -99, //0.0252, //-99, //0.02146 , // 0.0251,    /* 13-mixing, sinq2theta or sinsqtheta depending on config choice, if -999 then randomize  */
 double       input_mixing23          = -99, //0.533, //-99, //0.50,   // 0.513                         /* 23-mixing, sinq2theta or sinsqtheta depending on config choice */
 double       input_deltacp           = -99, //0.0,//-1.5707963,                            /* deltacp, if -999 then randomize */
 double       input_dm12sq            = -99, //7.58E-5,   //-999,                            /* Dm^2_{12}, if -999 then randomize */
 double       input_dm23sq            = -99, //2.481E-3, //-99, //2.32E-3,                            /* |Dm^2_{23}| */
 
 const char * config_file      = "analysis_configuration_1.7E20_prefit.xml",   /* analysis configuration file */

 bool rethrow_if_parameters_unphysical = true, // rethrow unphysical toys

  bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
  int          tweak_systematics       = 1,                              /* generate toy MC with randomized systematics? Numbers greater than 1 are special options */

  bool         save_in_xml_format      = true,                           /* also save in txt format */
  bool         save_in_root_format     = true,                           /* */  
  const char * points_to_add_1Rmu      = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  const char * points_to_add_1Re       = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  const char * points_to_add_1Rmu_RHC  = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  const char * points_to_add_1Re_RHC   = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  const char * points_to_add_MultiRmu  = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  const char * points_to_add_MultiRe   = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  //const char * xml_file_to_add         = "data_newstruct.xml",           /* xml with data in for when creating hybrid data/toy-MC experiments, eg data.xml. Always uses set 0 */
  const char * xml_file_to_add         = 0,           /* xml with data in for when creating hybrid data/toy-MC experiments, eg data.xml. Always uses set 0 */
  bool         save_fluct_in_xml       = true,                          /* save the values of the variations in the xml file? */
  const char * dmsqscheme              = "abs23",                        /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                       /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                   /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one"                           /* one==kNuOctantOne two==kNuOctantTwo */
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
  int n_points_per_line_xml = 5;

  TRandom rand(0);

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
  bool randomize_mixing12      = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13      = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_mixing23      = (TMath::Abs(input_mixing23 + 999) < kASmallNum);
  bool randomize_deltacp       = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq        = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);
  bool randomize_dm23sq        = (TMath::Abs(input_dm23sq   + 999) < kASmallNum);


  //Set to defaults
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1 && input_dm23sq > -900)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1 && input_dm12sq > -900)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90 && input_deltacp > -900) input_deltacp = conv::kGlobFitDeltaCP;

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
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << ((randomize_mixing12) ? "* randomized *" : input_mixing12) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << ((randomize_mixing13) ? "* randomized *" : input_mixing13) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << ((randomize_mixing23) ? "* randomized *" : input_mixing23) << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = " << ((randomize_deltacp) ? "* randomized *" : input_deltacp) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << ((randomize_dm23sq)  ? "* randomized *" : input_dm23sq ) << " eV^2/c^4" << endl;
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

     bool applied_negevents_condition = true;
     while(applied_negevents_condition) {

       // Randomize specified osc params for current toy MC expt
       if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
       if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
       if (randomize_mixing23) { input_mixing23 = util::RandomizeMixing23(); }
       if (randomize_deltacp ) { input_deltacp  = util::RandomizeCPV();       }
       if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();        }
       if (randomize_dm23sq  ) { input_dm23sq   = util::RandomizeDm23sq();        }

       // Set physics hypothesis
       if(ph.Hypothesis() == kHypNuMuDisapp2F) {
         ph.SetCurrValue(kPhypMixing23,  input_mixing23);
         ph.SetCurrValue(kPhypDm23sq,    input_dm23sq);
       }
       else if(ph.Hypothesis() == kHypOsc3F) {  
         ph.SetCurrValue(kPhypMixing12,  input_mixing12);
         ph.SetCurrValue(kPhypMixing13,  input_mixing13);
         ph.SetCurrValue(kPhypMixing23,  input_mixing23);
         ph.SetCurrValue(kPhypCPV    ,   input_deltacp);
         ph.SetCurrValue(kPhypDm12sq,    input_dm12sq);
         ph.SetCurrValue(kPhypDm23sq,    input_dm23sq);
       }





















       cout << "- Physics hypothesis set: " << ph << endl;

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
     TH1D* mock_data_1Re = new TH1D("","", util::NumErecoBins(kSKsFCFV1RingELike), util::ErecoBinEdges(kSKsFCFV1RingELike));
     for(int ib = 1; ib <= util::NumErecoBins(kSKsFCFV1RingELike); ib++) {
       mock_data_1Re->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsFCFV1RingELike)[ib-1]);
     }
     mock_data_1Re->SetDirectory(0);
#endif
#ifdef __RHC_RUNNING__
     // Get 1-ring mu-like Ereco spectrum RHC
     TH1D* mock_data_1Rmu_RHC = new TH1D("","", util::NumErecoBins(kSKsRHCFCFV1RingMuLike), util::ErecoBinEdges(kSKsRHCFCFV1RingMuLike));
     for(int ib = 1; ib <= util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++) {
       mock_data_1Rmu_RHC->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib-1]);
     }
     mock_data_1Rmu_RHC->SetDirectory(0);     
     // Get 1-ring e-like Ereco spectrum RHC
     TH1D* mock_data_1Re_RHC = new TH1D("","", util::NumErecoBins(kSKsRHCFCFV1RingELike), util::ErecoBinEdges(kSKsRHCFCFV1RingELike));
     for(int ib = 1; ib <= util::NumErecoBins(kSKsRHCFCFV1RingELike); ib++) {
       mock_data_1Re_RHC->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingELike)[ib-1]);
     }
     mock_data_1Re_RHC->SetDirectory(0);
#endif
#ifdef __MULTIRING_FHC_RUNNING__
     // Get Multi-ring mu-like Ereco spectrum
     TH1D* mock_data_MultiRmu = new TH1D("","", util::NumErecoBins(kSKsFCFVMultiRingMuLike), util::ErecoBinEdges(kSKsFCFVMultiRingMuLike));
     for(int ib = 1; ib <= util::NumErecoBins(kSKsFCFVMultiRingMuLike); ib++) {
       mock_data_MultiRmu->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsFCFVMultiRingMuLike)[ib-1]);
     }
     mock_data_MultiRmu->SetDirectory(0);     
     // Get Multi-ring e-like Ereco spectrum 
     TH1D* mock_data_MultiRe = new TH1D("","", util::NumErecoBins(kSKsFCFVMultiRingELike), util::ErecoBinEdges(kSKsFCFVMultiRingELike));
     for(int ib = 1; ib <= util::NumErecoBins(kSKsFCFVMultiRingELike); ib++) {
       mock_data_MultiRe->SetBinContent(ib, skobs->ErecoSpectrumArray(kSKsFCFVMultiRingELike)[ib-1]);
     }
     mock_data_MultiRe->SetDirectory(0);
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

  h_0->Write();
  h_1->Write();
  h_2->Write();
  h_3->Write();
  h_4->Write();
  h_4_Escale->Write();

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

