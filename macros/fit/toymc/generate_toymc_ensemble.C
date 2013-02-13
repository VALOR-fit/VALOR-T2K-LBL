//
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void generate_toymc_ensemble
(
  const char * tag                     = "toymc",                        /* */
  int          run                     = 0,                              /* */
  double       nexperiments            = 10,                             /* number of toy-MC expts to generate */
  const char * hypothesis              = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12          = -99,                            /* 12-mixing, sinq2theta or sinsqtheta depending on config choice, if -999 then randomize  */
  double       input_mixing13          = -99,                            /* 13-mixing, sinq2theta or sinsqtheta depending on config choice, if -999 then randomize  */
  double       input_mixing23          = -99,                            /* 23-mixing, sinq2theta or sinsqtheta depending on config choice */
  double       input_deltacp           = -99,                            /* deltacp, if -999 then randomize */
  double       input_dm12sq            = -99,                            /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq            = -99,                            /* |Dm^2_{23}| */
  const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
  bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
  int          tweak_systematics       = 1,                              /* generate toy MC with randomized systematics? Numbers greater than 1 are special options */
  bool         save_in_xml_format      = true,                           /* also save in txt format */
  bool         save_in_root_format     = true,                           /* */
  const char * points_to_add           = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
  bool         save_fluct_in_xml       = false                           /* save the values of the variations in the xml file? */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  TRandom rand(0);

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = kNmhNormal;

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
    case(1):
      tweak_sys_str = "All"; break;
    case(2):
      tweak_sys_str = "Special low energy numu flux increase"; break;
    case(3):
      tweak_sys_str = "Tweak all but SK and FSI"; break;
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
  cout << "       > " << util::AsString(kPhypDeltaCP)  << " = " << ((randomize_deltacp) ? "* randomized *" : input_deltacp) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << "   - Adding the following data-points in each toy-MC expt: " << ((points_to_add) ? points_to_add : "* none *") << endl;
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
    fxml << "       <one_mu_ring>    " << endl;
    fxml << "         <ereco_values> " << endl;
    ostringstream txt_filename;
    txt_filename << tag << "-" << run << ".txt";
    ftxt.open(txt_filename.str().c_str());
  }

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  for(int iexperiment=0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;

     // Randomize specified osc params for current toy MC expt
     if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
     if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
     if (randomize_deltacp ) { input_deltacp  = util::RandomizeDeltaCP();       }
     if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();        }

     // Set physics hypothesis
     PhysHypothesis ph(hypothesis);
     if(ph.Hypothesis() == kHypNuMuDisapp2F) {
       ph.SetCurrValue(kPhypMixing23,  input_mixing23);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq);
     } else
     if(ph.Hypothesis() == kHypOsc3F) {  
       ph.SetCurrValue(kPhypMixing12,  input_mixing12);
       ph.SetCurrValue(kPhypMixing13,  input_mixing13);
       ph.SetCurrValue(kPhypMixing23,  input_mixing23);
       ph.SetCurrValue(kPhypDeltaCP,   input_deltacp);
       ph.SetCurrValue(kPhypDm12sq,    input_dm12sq);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq);
     }
     //cout << "- Physics hypothesis set: " << ph << endl;

     if(tweak_systematics == 1) {
        // tweak all sources of systematic err (gaussian err)
        cout << "- Randomizing all systematic params in pdf..." << endl;
        skexp->TweakAllSystematics();
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
       skexp->TweakAllBANFFSystematics();
       skexp->TweakAllUncorrelatedXSecSystematics();
       skexp->PrintSystematicVariations();
     }

     // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs = new SuperkObservation;
     skobs->CreateMockData(*skexp, ph, stat_fluctuations);

     // Get 1-ring mu-like Ereco spectrum
     TH1D* mock_data = skobs->ErecoSpectrum(kSKsFCFV1RingMuLike);
     mock_data->SetDirectory(0);

     // If an XML file was specified, add an entry
     if(save_in_xml_format) {
       fxml << "           " << Form("<data_set_%d>",iexperiment) << endl;
       if(save_fluct_in_xml) {
	 fxml << "\t\t" << "<!-- ";
	 ph.Print(fxml);
	 if(tweak_systematics)
	   skexp->PrintSystematicVariations(fxml);
	 fxml << " -->" << endl;
       }
       fxml << "\t\t";
       // Adding Ereco values as a comma separated list.
       // Assuming uniform distribution within each bin, randomize values within each bin.
       int ievtot = 0;
       for(int ibin = 1; ibin <= mock_data->GetNbinsX(); ibin++) {
         double Nb   = mock_data->GetBinContent (ibin); // number of events in this Ereco bin
         double Emin = mock_data->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
         double Emax = mock_data->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
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
       if(points_to_add) {
         fxml << "," << points_to_add;
         ftxt << "," << points_to_add;
       }
       fxml << endl;
       ftxt << endl;
       fxml << "           " << Form("</data_set_%d>",iexperiment) << endl;
     }//xml

     // If a ROOT filename was specified, save the ROOT histogram too
     if(f) {
       if(save_in_root_format) {
 	  if(points_to_add) {
	    vector<string> vecp = util::Tokens(points_to_add,",");
	    assert(vecp.size()>0);
	    for(int ip=0; ip < vecp.size(); ip++) {
	      mock_data->Fill(atof(vecp[ip].c_str()));
	    }
          }
       }
       f->cd();
       mock_data->Write(Form("mock_data_%d",iexperiment));
     }

  }//iexperiment

  //
  // Close files
  //
  if(save_in_xml_format) {
    fxml << "         </ereco_values> " << endl;
    fxml << "       </one_mu_ring>    " << endl;
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

