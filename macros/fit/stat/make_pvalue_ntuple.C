//
// ------------------------------------------------------------------------------------------------
// Make an ntuple used for estimating p-values for the null hypothesis
// Save two test-statistics: \chi^2 and Nsk
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void make_pvalue_ntuple
(
  const char * tag                     = "pv",                           /* */
  int          run                     = 0,                              /* */
  bool         exclude_noosc           = true,                           /* perform experiment to exclude no oscillations, or to exclude maximal mixing */
  int          nexperiments            = 1000,                           /* */
  const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
  bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics       = true,                           /* generate toy MC with randomized systematics? */
  const char * fit_param_list          = "nuosc3f_numudisapp_fix4"       /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = kNmhNormal;

  // Print settings
  cout << "[Job inputs]" << endl;
  cout << "Job: " << tag << "/" << run << endl;
  cout << "Number of toy experiments : " << nexperiments << endl; 
  cout << "Systematic variations : "      << ((tweak_systematics) ?    "YES " : "NO") << endl;
  cout << "Statistical fluctuations : "   << ((stat_fluctuations) ?    "YES " : "NO") << endl;
  cout << "Reading analysis configuration from: " << config << endl;
  cout << *Opt::Instance() << endl;

  // Create file & ntuple to store results
  TFile f(Form("%s-run_%d.root",tag,run), "recreate");

  // Fit params - specify appropriate null hypothesis
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "Specified physics hypothesis: " << util::AsString(hyp) << endl;
  if(exclude_noosc) {
    if(hyp == kHypOsc3F) {
      //fparm->Update(kPhypMixing12, 0.); // allow solar osc. as part of the null hypothesis
      fparm->Update(kPhypMixing13, 0.);
      fparm->Update(kPhypMixing23, 0.);
      fparm->Update(kPhypDeltaCP,  0.);
    }
    else
    if(hyp == kHypNuMuDisapp2F) {
      fparm->Update(kPhypMixing23, 0.);
    }
    else {
      return;
    }
  }
  else {
    if(hyp == kHypOsc3F) {
      //fparm->Update(kPhypMixing12, 0.); // allow solar osc. as part of the null hypothesis
      fparm->Update(kPhypMixing13, 0.);
      fparm->Update(kPhypMixing23, 1.);
      fparm->Update(kPhypDeltaCP,  0.);
    }
    else
    if(hyp == kHypNuMuDisapp2F) {
      fparm->Update(kPhypMixing23, 1.);
    }
    else {
      return;
    }
  }

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Output ntuple
  TNtupleD * pvnt = new TNtupleD("pvnt","", "chisq_1Rmu:nsk_1Rmu");

  // Toy MC expt loop
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

     // Calc chisq
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     double chisq = fit_inputs->Chisquare();

     // Calc Nsk
     double nsk = skobs->NumOfEvents(kSKsFCFV1RingMuLike);

     // Print test statistics
     cout << "1Rmu: chisq = " << chisq << ", Nsk = " << nsk << endl;

     pvnt->Fill(chisq,nsk);

     delete skobs;
     delete fit_inputs;
  }

  // store ntuple with best-fit results
  f.cd();
  pvnt->Write("pvnt");

  f.Close();
}
