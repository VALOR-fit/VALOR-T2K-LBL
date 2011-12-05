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
  double       input_mixing12          = 0.86,                           /* if -999 then randomize */
  double       input_mixing13          = 0.10,                           /* if -999 then randomize */
  double       input_mixing23          = 1.00,                           /* */
  double       input_deltacp           = 0.00,                           /* if -999 then randomize */
  double       input_dm12sq            = 7.58E-5,                        /* if -999 then randomize */
  double       input_dm23sq            = 2.32E-3,                        /* */
  const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
  bool         tweak_systematics       = true,                           /* generate toy MC with randomized systematics? */
  bool         save_in_xml_format      = true,                           /* */
  bool         save_in_root_format     = false,                          /* */
  const char * points_to_add           = 0                               /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  TRandom rand(0);
  bool include_statistical_fluctuations = true; // <- statistical fluctuations

  // Check whether some osc. params need to be randomized
  bool randomize_mixing12      = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13      = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_deltacp       = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq        = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);
  bool randomized_osc_params   = (randomize_mixing12 || randomize_mixing13 || randomize_deltacp || randomize_dm12sq);

  // Print out job summary
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Input physics hypothesis: " << hypothesis << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl;
  cout << "   - Randomize systematics : "      << ((tweak_systematics) ? "YES " : "NO") << endl;
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > mixing12 = " << ((randomize_mixing12) ? "* randomized *" : input_mixing12) << endl;
  cout << "       > mixing13 = " << ((randomize_mixing13) ? "* randomized *" : input_mixing13) << endl;
  cout << "       > mixing23 = " << input_mixing23 << endl;
  cout << "       > delta_cp      = " << ((randomize_deltacp) ? "* randomized *" : input_deltacp) << endl;
  cout << "       > Dm12sq        = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
  cout << "       > |Dm23sq|      = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << "   - Adding the following data-points in each toy-MC expt: " << ((points_to_add) ? points_to_add : "* none *") << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing toy-MC expts in XML format: "   << ((save_in_xml_format) ?   "YES " : "NO") << endl;
  cout << "     - Storing toy-MC expts in ROOT format:  " << ((save_in_root_format) ?  "YES " : "NO") << endl;
  cout << endl;

  // Open files
  TFile * f = 0;
  if(save_in_root_format) {
    ostringstream root_filename;
    root_filename << tag << "-" << run << ".root";
    f = new TFile(root_filename.str().c_str(), "recreate");
  }

  ofstream fxml;
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
     if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
       ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
       ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
     } else
     if(ph.Hypothesis() == conv::kHypOsc3F) {  
       ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
       ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
       ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
       ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
       ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
       ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
     }
     //cout << "- Physics hypothesis set: " << ph << endl;

     if(tweak_systematics) {
        // tweak all sources of systematic err (gaussian err)
        cout << "- Randomizing all systematic params in pdf..." << endl;
        skexp->TweakAllSystematics();
        skexp->PrintSystematicVariations();
     }

     // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs = new SuperkObservation;
     skobs->CreateMockData(*skexp, ph, include_statistical_fluctuations);

     // Get 1-ring mu-like Ereco spectrum
     TH1D* mock_data = skobs->ErecoSpectrum(conv::kSKsFCFV1RingMuLike);
     mock_data->SetDirectory(0);

     // If an XML file was specified, add an entry
     if(save_in_xml_format) {
       fxml << "           " << Form("<data_set_%d>",iexperiment) << endl;
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
            }
            fxml << E; 
            ievtot++;
         }//iev
       }//ibin
       // add extra data-points, if any
       if(points_to_add) {
         fxml << "," << points_to_add;
       }
       fxml << endl;
       fxml << "           " << Form("</data_set_%d>",iexperiment) << endl;
     }//xml

     // If a ROOT filename was specified, save the ROOT histogram too
     if(f) {
       if(save_in_root_format) {
          vector<string> vecp = util::Tokens(points_to_add,",");
          assert(vecp.size()>0);
          for(int ip=0; ip < vecp.size(); ip++) {
            mock_data->Fill(atof(vecp[ip].c_str()));
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

