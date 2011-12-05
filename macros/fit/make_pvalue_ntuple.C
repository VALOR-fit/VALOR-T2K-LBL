//
// ------------------------------------------------------------------------------------------------
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void make_pvalue_ntuple(
  const char * tag                     = "pv",                           /* */
  int          run                     = 0,                              /* */
  int          nexperiments            = 10,                             /* */
  const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
  bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics       = true,                           /* generate toy MC with randomized systematics? */
  const char * fit_param_list          = "numudisapp2f"                  /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  cout << "[Job inputs]" << endl;
  cout << "Job: " << tag << "/" << run << endl;
  cout << "Number of toy experiments : " << nexperiments << endl; 
  cout << "Systematic variations : "      << ((tweak_systematics) ?    "YES " : "NO") << endl;
  cout << "Statistical fluctuations : "   << ((stat_fluctuations) ?    "YES " : "NO") << endl;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // create file to store results
  TFile f(Form("%s-run_%d.root",tag,run), "recreate");

  // fit params - specify appropriate null hypothesis
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "Specified physics hypothesis: " << util::AsString(hyp) << endl;
  if(hyp == conv::kHypOsc3F) {
  //fparm->Update(t2kana::conv::kPhypMixing12, 0.);
    fparm->Update(t2kana::conv::kPhypMixing13, 0.);
    fparm->Update(t2kana::conv::kPhypMixing23, 0.);
    fparm->Update(t2kana::conv::kPhypDeltaCP,  0.);
  }
  else
  if(hyp == conv::kHypNuMuDisapp2F) {
    fparm->Update(t2kana::conv::kPhypMixing23, 0.);
  }
  else {
    return;
  }

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // store results
  TNtupleD * pvnt = new TNtupleD("pvnt","", "iexp:chisq");

  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     SuperkObservation * skobs = new SuperkObservation;
     if(tweak_systematics) {
        // tweak all sources of systematic err (gaussian err)
        cout << "- Randomizing all systematic params in pdf..." << endl;
        skexp->TweakAllSystematics();
        skexp->PrintSystematicVariations();
        // generate mock data using the current PDF (incl. systematic variations)
        // + statistical fluctuations if requested (Poisson fluctuation in each bin)
        cout << "- Generating toy-MC..." << endl;
        skobs->CreateMockData(*skexp, *fparm, stat_fluctuations);
        // reset all systematics so that the PDF used in the fitter is the nominal one
        cout << "- Reseting all systematic params in pdf..." << endl;
        skexp->ResetSystematicVariations();
        //skexp->PrintSystematicVariations();
     } else {
        // generate mock data using the current PDF (without systematic variations)
        // with statistical fluctuations if requested (Poisson fluctuation in each bin)
        cout << "- Generating toy-MC..." << endl;
        skobs->CreateMockData(*skexp, *fparm, stat_fluctuations);
     }

     // calc & store chisq
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     double chisq = fit_inputs->Chisquare();
     cout << "Chisq = " << chisq << endl;
     pvnt->Fill(iexperiment, chisq);

     delete skobs;
     delete fit_inputs;
  }

  // store ntuple with best-fit results
  f.cd();
  pvnt->Write("pvnt");

  f.Close();
}
