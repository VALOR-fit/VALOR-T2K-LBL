//
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void generate_toymc_ensemble_joint
(
  const char * tag                     = "toymc",                        /* */
  int          run                     = 0,                              /* */
  double       nexperiments            = 3,                             /* number of toy-MC expts to generate */
  const char * hypothesis              = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12          = -99,                            /* 12-mixing, sinq2theta or sinsqtheta depending on config choice, if -999 then randomize  */
  double       input_mixing13          = 0.0251,                            /* 13-mixing, sinq2theta or sinsqtheta depending on config choice, if -999 then randomize  */
  double       input_mixing23          = 0.513,                            /* 23-mixing, sinq2theta or sinsqtheta depending on config choice */
  double       input_deltacp           = 0.0,//-1.5707963,                            /* deltacp, if -999 then randomize */
  double       input_dm12sq            = -99,                            /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq            = 2.4E-3,                            /* |Dm^2_{23}| */
  const char * config_file             = "analysis_configuration_test.xml",   /* analysis configuration file */
  bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
  int          tweak_systematics       = 1,                              /* generate toy MC with randomized systematics? Numbers greater than 1 are special options */
  bool         save_in_xml_format      = true,                           /* also save in txt format */
  bool         save_in_root_format     = true,                           /* */
  const char * points_to_add_1Rmu      = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  const char * points_to_add_1Re       = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  const char * points_to_add_1Rmu_RHC      = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  const char * points_to_add_1Re_RHC       = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  bool         save_fluct_in_xml       = false,                           /* save the values of the variations in the xml file? */
  const char * dmsqscheme              = "abs23",                        /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                       /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                   /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one"                           /* one==kNuOctantOne two==kNuOctantTwo */
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  TRandom rand(0);

  TH1D * h_0 = new TH1D ("h_0", "h_0", 100, 0, 10);
  TH1D * h_1 = new TH1D ("h_1", "h_1", 100, 0, 10);
  TH1D * h_2 = new TH1D ("h_2","h_2", 100, 0, 10);
  TH1D * h_3 = new TH1D ("h_3","h_3", 100, 0, 10);
  TH1D * h_4 = new TH1D ("h_4","h_4", 100, 0, 10);
  TH2D * h_4_Escale = new TH2D ("h_4_Escale","h_4_Escale", 100, -5, 5, 100, 0, 10);

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

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
  bool randomize_deltacp       = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq        = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);
  bool randomized_osc_params   = (randomize_mixing12 || randomize_mixing13 || randomize_deltacp || randomize_dm12sq);

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

  // Print out job summary
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Input physics hypothesis: " << hypothesis << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl;
  cout << "   - Randomize systematics : "      << tweak_sys_str << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ? "YES " : "NO") << endl;
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << ((randomize_mixing12) ? "* randomized *" : input_mixing12) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << ((randomize_mixing13) ? "* randomized *" : input_mixing13) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << input_mixing23 << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = " << ((randomize_deltacp) ? "* randomized *" : input_deltacp) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << "   - Adding the following data-points in each toy-MC expt: " << endl;
  cout << "    For 1R mu-like: " << ((points_to_add_1Rmu) ? points_to_add_1Rmu : "* none *") << endl;
  cout << "    For 1R e-like: " << ((points_to_add_1Re) ? points_to_add_1Re : "* none *") << endl;
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
       if (randomize_deltacp ) { input_deltacp  = util::RandomizeCPV();       }
       if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();        }

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
	 //	 skexp->TweakSystematic("sk_energy_scale",0.0,true);
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
       applied_negevents_condition = skobs->CreateMockData(*skexp, ph, stat_fluctuations);
       cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
     }//applied_negevents_condition

     double escale = skexp->GetSuperkDetectorFSI()->CurrentValue(kNuisSuperkEnergyScale);
     cout << "escale = " << escale << endl;

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

     // If an XML file was specified, add an entry

     if(save_in_xml_format) {

       // Now the 1Rmu part
       fxml << "           " << Form("<data_set_%d>",iexperiment) << endl;
       fxml << "       <one_mu_ring>    " << endl;
       fxml << "         <ereco_values> " << endl;
       //fxml << "           " << Form("<data_set_%d>",iexperiment) << endl;
       
       if(save_fluct_in_xml) {
	 fxml << "\t\t" << "<!-- ";
	 ph.Print(fxml);
	 if(tweak_systematics)
	   skexp->PrintSystematicVariations(fxml);
	 fxml << " -->" << endl;
       }
       fxml << "\t\t";
       ftxt << "1Rmu:";
       ftxt << endl;
       
       // Adding Ereco values as a comma separated list.
       // Assuming uniform distribution within each bin, randomize values within each bin.
       int ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_1Rmu->GetNbinsX(); ibin++) {
         double Nb   = mock_data_1Rmu->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_1Rmu->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_1Rmu->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
         for(int iev = 0; iev < Nb; iev++) {
            double E = Emin + rand.Uniform() * (Emax-Emin);
            if(ievtot > 0) {
              fxml << ","; 
	      ftxt << ","; 
            }
            fxml << E; 
            ftxt << E; 
            ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(points_to_add_1Rmu) {
         fxml << "," << points_to_add_1Rmu;
         ftxt << "," << points_to_add_1Rmu;
       }
       fxml << endl;
       ftxt << endl;
       //fxml << "           " << Form("</data_set_%d>",iexperiment) << endl;
       fxml << "         </ereco_values> " << endl;
       fxml << "       </one_mu_ring>    " << endl;
       
              
       //Now the 1Re part
       ftxt << endl;
       ftxt << "1Re:";
       ftxt << endl;
       fxml << "       <one_e_ring>    " << endl;
       fxml << "         <ereco_values> " << endl;
       //fxml << "           " << Form("<data_set_%d>",iexperiment) << endl;
       
       ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_1Re->GetNbinsX(); ibin++) {
         double Nb   = mock_data_1Re->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_1Re->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_1Re->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
         for(int iev = 0; iev < Nb; iev++) {
            double E = Emin + rand.Uniform() * (Emax-Emin);
            if(ievtot > 0) {
              fxml << ","; 
       	      ftxt << ","; 
            }
            fxml << E; 
            ftxt << E; 
            ievtot++;
         }//iev
       }//ibin


       // add extra data-points, if any
       if(points_to_add_1Re) {
         fxml << "," << points_to_add_1Re;
         ftxt << "," << points_to_add_1Re;
       }

       fxml << endl;
       ftxt << endl;
       //fxml << "           " << Form("</data_set_%d>",iexperiment) << endl;
       fxml << "         </ereco_values> " << endl;
       fxml << "       </one_e_ring>    " << endl;
       

       // Now the 1Rmu RHC part
       fxml << "       <one_mu_ring_rhc>    " << endl;
       fxml << "         <ereco_values> " << endl;
       //fxml << "           " << Form("<data_set_%d>",iexperiment) << endl;
       fxml << "\t\t";
       ftxt << "1Rmu RHC:";
       ftxt << endl;
       
       // Adding Ereco values as a comma separated list.
       // Assuming uniform distribution within each bin, randomize values within each bin.
       int ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_1Rmu_RHC->GetNbinsX(); ibin++) {
         double Nb   = mock_data_1Rmu_RHC->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_1Rmu_RHC->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_1Rmu_RHC->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
         for(int iev = 0; iev < Nb; iev++) {
            double E = Emin + rand.Uniform() * (Emax-Emin);
            if(ievtot > 0) {
              fxml << ","; 
	      ftxt << ","; 
            }
            fxml << E; 
            ftxt << E; 
            ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(points_to_add_1Rmu_RHC) {
         fxml << "," << points_to_add_1Rmu_RHC;
         ftxt << "," << points_to_add_1Rmu_RHC;
       }
       fxml << endl;
       ftxt << endl;
       //fxml << "           " << Form("</data_set_%d>",iexperiment) << endl;
       fxml << "         </ereco_values> " << endl;
       fxml << "       </one_mu_ring_rhc>    " << endl;


       //Now the 1Re RHC part
       ftxt << endl;
       ftxt << "1Re RHC:";
       ftxt << endl;
       fxml << "       <one_e_ring_rhc>    " << endl;
       fxml << "         <ereco_values> " << endl;
       //fxml << "           " << Form("<data_set_%d>",iexperiment) << endl;
       
       ievtot = 0;
       for(int ibin = 1; ibin <= mock_data_1Re_RHC->GetNbinsX(); ibin++) {
         double Nb   = mock_data_1Re_RHC->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data_1Re_RHC->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data_1Re_RHC->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
         for(int iev = 0; iev < Nb; iev++) {
            double E = Emin + rand.Uniform() * (Emax-Emin);
            if(ievtot > 0) {
              fxml << ","; 
       	      ftxt << ","; 
            }
            fxml << E; 
            ftxt << E; 
            ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(points_to_add_1Re_RHC) {
         fxml << "," << points_to_add_1Re_RHC;
         ftxt << "," << points_to_add_1Re_RHC;
       }
       fxml << endl;
       ftxt << endl;
       //fxml << "           " << Form("</data_set_%d>",iexperiment) << endl;
       fxml << "         </ereco_values> " << endl;
       fxml << "       </one_e_ring_rhc>    " << endl;

       fxml << "           " << Form("</data_set_%d>",iexperiment) << endl;
       
     }//xml

     // If a ROOT filename was specified, save the ROOT histogram too
     if(f) {
       if(save_in_root_format) {
 	  if(points_to_add_1Rmu) {
	    vector<string> vecp = util::Tokens(points_to_add_1Rmu,",");
	    assert(vecp.size()>0);
	    for(int ip=0; ip < vecp.size(); ip++) {
	      mock_data_1Rmu->Fill(atof(vecp[ip].c_str()));
	    }
          }
	  if(points_to_add_1Re) {
	    vector<string> vecp = util::Tokens(points_to_add_1Re,",");
	    assert(vecp.size()>0);
	    for(int ip=0; ip < vecp.size(); ip++) {
	      mock_data_1Re->Fill(atof(vecp[ip].c_str()));
	    }
          }
	  
	  if(points_to_add_1Rmu_RHC) {
	    vector<string> vecp = util::Tokens(points_to_add_1Rmu_RHC,",");
	    assert(vecp.size()>0);
	    for(int ip=0; ip < vecp.size(); ip++) {
	      mock_data_1Rmu_RHC->Fill(atof(vecp[ip].c_str()));
	    }
          }
	  if(points_to_add_1Re_RHC) {
	    vector<string> vecp = util::Tokens(points_to_add_1Re_RHC,",");
	    assert(vecp.size()>0);
	    for(int ip=0; ip < vecp.size(); ip++) {
	      mock_data_1Re_RHC->Fill(atof(vecp[ip].c_str()));
	    }
          }
       }
       f->cd();
       mock_data_1Rmu->Write(Form("mock_data_1Rmu_%d",iexperiment));
       mock_data_1Rmu_RHC->Write(Form("mock_data_1Rmu_RHC_%d",iexperiment));
       mock_data_1Re->Write(Form("mock_data_1Re_%d",iexperiment));
       mock_data_1Re_RHC->Write(Form("mock_data_1Re_RHC_%d",iexperiment));

       // Object for storing input oscillation and nuisance params in each toy-MC
       PhysParamList * iparm = new PhysParamList;
       iparm->Clear();
       // keep track of oscillation parameters
       PhysParamList true_oscillation_params = ph.AsPhysParamList();
       iparm->Append(true_oscillation_params);
       // keep track of all systematic params
       PhysParamList true_nuisance_params = skexp->Tweaks();
       iparm->Append(true_nuisance_params);
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

