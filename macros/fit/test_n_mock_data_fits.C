//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are generated on the fly.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_n_mock_data_fits(
  const char * tag                     = "toyexptfits",                  /* */
  int          run                     = 0,                              /* */
  int          nexperiments            = 1000,                           /* */
  double       input_sinsq2theta12     = 0.86,                           /* if -1 then randomize */
  double       input_sinsq2theta13     = 0.10,                           /* if -1 then randomize */
  double       input_sinsq2theta23     = 1.00,                           /* */
  double       input_deltacp           = 0.00,                           /* if -1 then randomize */
  double       input_dm12sq            = 7.58E-5,                        /* if -1 then randomize */
  double       input_dm23sq            = 2.32E-3,                        /* */
  const char * fit_param_list          = "numudisapp2f",                 /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fitter_name             = "minuit",                       /* fitter name (see run_fitter.C) */
  bool         fluctuate_stat          = true,                           /* include statistical fluctuations? */
  bool         include_systematics     = false,                          /* include systematic variation in generated toy MC expt? */
  bool         store_chisq_surf        = false,                          /* by default storing only best-fit values, enable this to store the chisq surfaces   */
  bool         store_toy_expts         = false,                          /* by default storing only best-fit values, enable this to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false,                          /* by default storing only best-fit values, enable this to store the best-fit spectra */
  const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
  bool         randomize_once          = false                            /* randomize non-23 parameters every experiment or every file? */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  cout << "[Job inputs]" << endl;
  cout << "Job: " << tag << "/" << run << endl;
  cout << "Fit Parameter List: " << fit_param_list << endl;
  cout << "Fitter: " << fitter_name << endl;
  cout << "Number of toy experiments: " << nexperiments << endl; 
  cout << "Systematics : "                << ((include_systematics) ?    "YES " : "NO") << endl;
  cout << "Statistical fluctuations: "    << ((fluctuate_stat) ?         "YES " : "NO") << endl;
  cout << "Storing chisq surfaces: "      << ((store_chisq_surf) ?       "YES " : "NO") << endl;
  cout << "Storing toy MC experiments:  " << ((store_toy_expts) ?        "YES " : "NO") << endl;
  cout << "Storing best-fit spectra: "    << ((store_bestfit_spectrum) ? "YES " : "NO") << endl;
  cout << "True oscillation parameters: " << endl;
  cout << " > sinsq2theta12 = " << input_sinsq2theta12 << endl;
  cout << " > sinsq2theta13 = " << input_sinsq2theta13 << endl;
  cout << " > sinsq2theta23 = " << input_sinsq2theta23 << endl;
  cout << " > delta_cp      = " << input_deltacp << endl;
  cout << " > Dm12sq        = " << input_dm12sq << " eV^2/c^4" << endl;
  cout << " > |Dm23sq|      = " << input_dm23sq << " eV^2/c^4" << endl;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // create file to store results
  TFile f(Form("%s-sinsq2theta23_%.2f-dm23sq_%.5feV2-run_%d.root",
          tag, input_sinsq2theta23,input_dm23sq,run), "recreate");

  // fit parameter list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "Specified physics hypothesis: " << util::AsString(hyp) << endl;

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // store oscillation fit results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","",
     "iexp:bestfit_sinsq2theta12:bestfit_sinsq2theta13:bestfit_sinsq2theta23:bestfit_deltacp:bestfit_dm12sq:bestfit_dm23sq:chisq:fitstatus:input_sinsq2theta12:input_sinsq2theta13:input_sinsq2theta23:input_deltacp:input_dm12sq:input_dm23sq");

  // init best-fit values
  double bestfit_sinsq2theta12 = -1;
  double bestfit_sinsq2theta13 = -1;
  double bestfit_sinsq2theta23 = -1;
  double bestfit_deltacp       = -1;
  double bestfit_dm12sq        = -1;
  double bestfit_dm23sq        = -1;
  double bestfit_chisq         = -1;

  //remember what the user started with
  double user_input_sinsq2theta12 = input_sinsq2theta12;
  double user_input_sinsq2theta13 = input_sinsq2theta13;
  double user_input_sinsq2theta23 = input_sinsq2theta23;
  double user_input_deltacp = input_deltacp;
  double user_input_dm12sq = input_dm12sq;
  double user_input_dm23sq = input_dm23sq;

  // loop over toy-MC experiments
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // randomize specified osc params for current toy MC expt
     if ( randomize_once ) {
       if (input_sinsq2theta12 < 0) { input_sinsq2theta12 = util::RandomizeSinsq2Theta12(); }
       if (input_sinsq2theta13 < 0) { input_sinsq2theta13 = util::RandomizeSinsq2Theta13(); }
       if (input_deltacp       < 0) { input_deltacp       = util::RandomizeDeltaCP();       }
       if (input_dm12sq        < 0) { input_dm12sq        = util::RandomizeDm12sq();        }
     }
     else {
       if (user_input_sinsq2theta12 < 0) { input_sinsq2theta12 = util::RandomizeSinsq2Theta12(); }
       if (user_input_sinsq2theta13 < 0) { input_sinsq2theta13 = util::RandomizeSinsq2Theta13(); }
       if (user_input_deltacp       < 0) { input_deltacp       = util::RandomizeDeltaCP();       }
       if (user_input_dm12sq        < 0) { input_dm12sq        = util::RandomizeDm12sq();        }
     }

     // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
     if(hyp == conv::kHypOsc3F) {
       fparm->Update(t2kana::conv::kPhypSinsq2Theta12, input_sinsq2theta12);
       fparm->Update(t2kana::conv::kPhypSinsq2Theta13, input_sinsq2theta13);
       fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
       fparm->Update(t2kana::conv::kPhypDeltaCP,       input_deltacp);
       fparm->Update(t2kana::conv::kPhypDm12sq,        input_dm12sq);
       fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);
       fparm->SetInitValue(t2kana::conv::kPhypSinsq2Theta12, input_sinsq2theta12);
       fparm->SetInitValue(t2kana::conv::kPhypSinsq2Theta13, input_sinsq2theta13);
       fparm->SetInitValue(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
       fparm->SetInitValue(t2kana::conv::kPhypDeltaCP,       input_deltacp);
       fparm->SetInitValue(t2kana::conv::kPhypDm12sq,        input_dm12sq);
       fparm->SetInitValue(t2kana::conv::kPhypDm23sq,        input_dm23sq);
     }
     else
     if(hyp == conv::kHypNuMuDisapp2F) {
       fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
       fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);
       fparm->SetInitValue(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
       fparm->SetInitValue(t2kana::conv::kPhypDm23sq,        input_dm23sq);
     }

     //
     // get superk spectrum observation (toy experiment)
     //
     SuperkObservation * skobs = new SuperkObservation;
     if(include_systematics) {
        // tweak all sources of systematic err (gaussian err)
        skexp->TweakAllSystematics();     
        // generate mock data using the current PDF (incl. systematic variations) 
        // + statistical fluctuations (Poisson fluctuation in each bin)
        skobs->CreateMockData(*skexp, *fparm, fluctuate_stat);
        // reset all systematics so that the PDF used in the fitter is the nominal one
        skexp->ResetSystematicVariations();
     } else {       
        // generate mock data using the current PDF (without systematic variations) 
        // with statistical fluctuations (Poisson fluctuation in each bin)
        skobs->CreateMockData(*skexp, *fparm, fluctuate_stat); 
     }

     // get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     NuMuDisappFitter * fitter = new NuMuDisappFitter(fit_inputs,fitter_name);
     fitter->ActualOscParams(input_sinsq2theta23, input_dm23sq);

     double true_chisq = fit_inputs->Chisquare();
     cout << "True-fit:"
	  << " sinsq2theta12 = "  << input_sinsq2theta12
	  << ", sinsq2theta13 = " << input_sinsq2theta13
	  << ", sinsq2theta23 = " << input_sinsq2theta23
	  << ", delta_cp = " << input_deltacp      
	  << ", Dm12^2 = "   << input_dm12sq << " eV^2/c^4"
	  << ", |Dm23^2| = " << input_dm23sq << " eV^2/c^4"
	  << " ==> chisq = " << true_chisq << endl;

     // do the fit
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();

     // get best-fit params
     bestfit_chisq = fit_inputs->Chisquare();
     if(hyp == conv::kHypOsc3F) {
        bestfit_sinsq2theta12 = fparm->CurrValue(t2kana::conv::kPhypSinsq2Theta12);
        bestfit_sinsq2theta13 = fparm->CurrValue(t2kana::conv::kPhypSinsq2Theta13);
        bestfit_sinsq2theta23 = fparm->CurrValue(t2kana::conv::kPhypSinsq2Theta23);
        bestfit_deltacp       = fparm->CurrValue(t2kana::conv::kPhypDeltaCP);
        bestfit_dm12sq        = fparm->CurrValue(t2kana::conv::kPhypDm12sq);
        bestfit_dm23sq        = fparm->CurrValue(t2kana::conv::kPhypDm23sq);
        cout << "Best-fit:"
             << " sinsq2theta12 = "  << bestfit_sinsq2theta12
             << ", sinsq2theta13 = " << bestfit_sinsq2theta13
             << ", sinsq2theta23 = " << bestfit_sinsq2theta23
             << ", delta_cp = " << bestfit_deltacp      
             << ", Dm12^2 = "   << bestfit_dm12sq << " eV^2/c^4"
             << ", |Dm23^2| = " << bestfit_dm23sq << " eV^2/c^4"
             << " ==> chisq = " << bestfit_chisq << endl;
     } else
     if(hyp == conv::kHypNuMuDisapp2F) {
        bestfit_sinsq2theta23 = fparm->CurrValue(t2kana::conv::kPhypSinsq2Theta23);
        bestfit_dm23sq        = fparm->CurrValue(t2kana::conv::kPhypDm23sq);
        cout << "Best-fit:"
             << " sinsq2theta23 = " << bestfit_sinsq2theta23
             << ", |Dm23^2| = " << bestfit_dm23sq << " eV^2/c^4"
             << " ==> chisq = " << bestfit_chisq << endl;
     }

     // store best-fit values
     nfitnt->Fill(iexperiment, 
                  bestfit_sinsq2theta12, bestfit_sinsq2theta13, bestfit_sinsq2theta23, 
                  bestfit_deltacp, bestfit_dm12sq, bestfit_dm23sq, bestfit_chisq, fitstatus,
		  input_sinsq2theta12, input_sinsq2theta13, input_sinsq2theta23, 
                  input_deltacp, input_dm12sq, input_dm23sq);

     // store toy-MC expt (if asked)
     if(store_toy_expts) {
        TH1D * toyexpspectrum = new TH1D(*skobs->ErecoSpectrum(kSKsFCFV1RingMuLike));
        f.cd();
        toyexpspectrum->SetName(Form("toyexp_ereco_spectrum_1Rmu_%d",iexperiment));
        toyexpspectrum->Write();
     }
     // store best-fit spectrum (if asked)
     if(store_bestfit_spectrum) {
        fparm->Update(t2kana::conv::kPhypSinsq2Theta23, bestfit_sinsq2theta23);
        fparm->Update(t2kana::conv::kPhypDm23sq,        bestfit_dm23sq);
        skexp->Generate(*fparm);
        TH1D * bestfitspectrum = new TH1D(*skexp->CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike));
        bestfitspectrum->SetName(Form("bestfit_ereco_spectrum_1Rmu_%d",iexperiment));
        f.cd();
        bestfitspectrum->Write();
     }
     // store chisq surface (if asked)
     if(store_chisq_surf) {
        //custom mode builds chisqsurf automatically, need to tell minuit mode to build it
        if(fitter_name=="minuit") {
	  cout << "Building chisqsurf" << endl;
	  //need to use Draw since BuildChisqSurface is a private method
	  fitter->Draw(new TCanvas(),"ChisqSurface");
	}
        TH2D * chisqsurf = new TH2D(*fitter->GetChisqSurface()); 
	if(chisqsurf) {
	  cout << "got the chisqsurf" << endl;
	  chisqsurf->SetName(Form("chisqsurf_%d",iexperiment)); 
	  f.cd();
	  chisqsurf->Write();
	}
	else
	  cout << "problem getting the chisqsurf" << endl;
     }

     delete skobs;
     delete fit_inputs;
     delete fitter;
  }

  // store ntuple with best-fit results
  f.cd();
  nfitnt->Write("nfitnt");

  f.Close();
}
