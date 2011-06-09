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
  double       input_sinsq2theta23     = 1.0,                            /* */
  double       input_dm23sq            = 2.32E-3,                        /* */
  const char * fitter_name             = "minuit",                       /* fitter name ("minuit", "custom" (2-D grid scan) */
  bool         fluctuate_stat          = true,                           /* include statistical fluctuations? */
  bool         include_systematics     = false,                          /* include systematic variation in generated toy MC expt? */
  bool         store_chisq_surf        = false,                          /* by default storing only best-fit values, enable this to store the chisq surfaces   */
  bool         store_toy_expts         = false,                          /* by default storing only best-fit values, enable this to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false,                          /* by default storing only best-fit values, enable this to store the best-fit spectra */
  const char * config_file             = "analysis_configuration.xml"    /* analysis configuration file */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  cout << "[Job inputs]" << endl;
  cout << "Job: " << tag << "/" << run << endl;
  cout << "Fitter : " << fitter_name << endl;
  cout << "Systematics : " << ((include_systematics) ? "YES " : "NO") << endl;
  cout << "Statistical fluctuations : " << ((fluctuate_stat) ? "YES " : "NO") << endl;
  cout << "Number of toy experiments : " << nexperiments << endl; 
  cout << "True oscillation parameters : (sinsq2theta23 = " 
       << input_sinsq2theta23 << ", dm23sq = " << input_dm23sq << " eV^2/c^4" << ")" << endl;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // create file to store results
  TFile f(Form("%s-sinsq2theta23_%.2f-dm23sq_%.5feV2-run_%d.root",
          tag, input_sinsq2theta23,input_dm23sq,run), "recreate");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // store results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","", "iexp:sinsq2theta23:dm23sq:chisq");

  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

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

     // do the fit
     fitter->DoTheFit();

     double bestfit_sinsq2theta23 = fparm->CurrValue(t2kana::conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(t2kana::conv::kPhypDm23sq);
     double bestfit_chisq         = fit_inputs->Chisquare();

     cout << "Best-fit:"
          << " sinsq2theta23 = " << bestfit_sinsq2theta23
          << ", dm23sq = " << bestfit_dm23sq << " eV^2/c^4"
          << " ==> chisq = " << bestfit_chisq << endl;

     // store best-fit values for the current toy-MC expt
     nfitnt->Fill(iexperiment, bestfit_sinsq2theta23, bestfit_dm23sq, bestfit_chisq);

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
        TH2D * chisqsurf = new TH2D(*fitter->GetChisqSurface()); 
        chisqsurf->SetName(Form("chisqsurf_%d",iexperiment)); 
        f.cd();
        chisqsurf->Write();
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
