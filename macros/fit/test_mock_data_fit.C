//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits for mock-data. 
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_mock_data_fit
(
   double       input_sinsq2theta12 = 0.86,                         /* */
   double       input_sinsq2theta13 = 0.10,                         /* */
   double       input_sinsq2theta23 = 1.00,                         /* */
   double       input_deltacp       = 0.00,                         /* */
   double       input_dm12sq        = 7.58E-5,                      /* */
   double       input_dm23sq        = 2.32E-3,                      /* */
   bool         fluctuate_stat      = true,                         /* include statistical fluctuations? */
   bool         include_systematics = false,                        /* include systematic variation in generated toy MC expt? */
   const char * config_file         = "analysis_configuration.xml", /* analysis configuration file */
   const char * fit_param_list      = "numudisapp2f",               /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
   const char * fit_option          = "minuit",                     /* oscillation fitter (see run_fitter.C) */
   const char * output_option       = "Observation+BestFitSpectra+ChisqSurface+Canvases",  /* what to write-out once the fit is done */
   bool         print_toy_mc        = false                         /* print-out toy-MC E_reco values? */
)
{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Systematic variations: "    << ((include_systematics) ? "YES " : "NO") << endl;
  cout << "   - Statistical fluctuations: " << ((fluctuate_stat) ? "YES " : "NO") << endl;
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > sinsq2theta12 = " << input_sinsq2theta12 << endl;
  cout << "       > sinsq2theta13 = " << input_sinsq2theta13 << endl;
  cout << "       > sinsq2theta23 = " << input_sinsq2theta23 << endl;
  cout << "       > delta_cp      = " << input_deltacp << endl;
  cout << "       > Dm12sq        = " << input_dm12sq << " eV^2/c^4" << endl;
  cout << "       > |Dm23sq|      = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit parameter list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "- Output option: " << output_option << endl;
  cout << endl;

  // instantiate fit params & set specified inputs as necessary
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;
  if(hyp == conv::kHypOsc3F) {
    fparm->Update(t2kana::conv::kPhypSinsq2Theta12, input_sinsq2theta12);
    fparm->Update(t2kana::conv::kPhypSinsq2Theta13, input_sinsq2theta13);
    fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
    fparm->Update(t2kana::conv::kPhypDeltaCP,       input_deltacp);
    fparm->Update(t2kana::conv::kPhypDm12sq,        input_dm12sq);
    fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);
  }
  else 
  if(hyp == conv::kHypNuMuDisapp2F) {
    fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
    fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);
  } 
  else {
    return;
  }

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  if(include_systematics) {
     // tweak all sources of systematic err (gaussian err)
     cout << "- Randomizing all systematic params in pdf..." << endl;
     skexp->TweakAllSystematics();
   //skexp->PrintSystematicVariations();
     // generate mock data using the current PDF (incl. systematic variations)
     // + statistical fluctuations (Poisson fluctuation in each bin)
     cout << "- Generating toy-MC..." << endl;
     skobs->CreateMockData(*skexp, *fparm, fluctuate_stat);
     // reset all systematics so that the PDF used in the fitter is the nominal one
     cout << "- Reseting all systematic params in pdf..." << endl;
     skexp->ResetSystematicVariations();
     //skexp->PrintSystematicVariations();
  } else {
     // generate mock data using the current PDF (without systematic variations)
     // with statistical fluctuations (Poisson fluctuation in each bin)
     cout << "- Generating toy-MC..." << endl;
     skobs->CreateMockData(*skexp, *fparm, fluctuate_stat);
  }

  // get fitter
  FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
  FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

  // do the fit
  fitter->DoTheFit();

  // write-out results
  fitter->Analyze(output_option);
  fitter->Save("mockfit.root");

  if(print_toy_mc) {
    // Print-out mock-data (Ereco values) as a comma separated list.
    // Assuming uniform distribution within each bin, randomize values within each bin.  
    cout << "\nToy-MC experiment, E_reco values: " << endl;
    TRandom3 rand(0);
    TH1D* mock_data = skobs->ErecoSpectrum(conv::kSKsFCFV1RingMuLike);
    bool is_first=true;
    for(int ibin = 1; ibin <= mock_data->GetNbinsX(); ibin++) {
      double Nb   = mock_data->GetBinContent (ibin); // number of events in this Ereco bin
      double Emin = mock_data->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
      double Emax = mock_data->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
      for(int iev = 0; iev < Nb; iev++) {
        double E = Emin + rand.Uniform() * (Emax-Emin);
        string delimeter = ",";
        if(is_first) {
          delimeter = "";
          is_first = false;
        }
        cout << delimeter << E;
      }//iev
    }//ibin
    cout << endl;
  }//print?

}
