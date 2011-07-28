//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are read from the input XML file.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_n_mock_data_fits_xml
(
  const char * tag,                                                /* */
  int          run,                                                /* */
  const char * xmlfile,                                            /* input XML file with toy MC experiments */
  int          first_experiment,                                   /* */
  int          last_experiment,                                    /* */
  const char * fitparamlist        = "numudisapp2f",               /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fitter_name         = "minuit",                     /* fitter name (see run_fitter.C) */
  const char * config_file         = "analysis_configuration.xml"  /* analysis configuration file */
)
{
  using namespace std;
  using namespace t2kana;

  cout << "[Job inputs]" << endl;
  cout << "Job: " << tag << "/" << run << endl;
  cout << "Fitter: " << fitter_name << endl;
  cout << "Reading experiments from: " << xmlfile << endl; 
  cout << "Fitting experiments: [" << first_experiment << ", " << last_experiment << "]" << endl; 

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // fit param list
  FitParamList * fparm = FitParamList::NamedParameterList(fitparamlist);

  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "Specified physics hypothesis: " << util::AsString(hyp) << endl;

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // store oscillation fit results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","", 
     "iexp:sinsq2theta12:sinsq2theta13:sinsq2theta23:deltacp:dm12sq:dm23sq:chisq:fitstatus");

  // init best-fit values
  double bestfit_sinsq2theta12 = -1;
  double bestfit_sinsq2theta13 = -1;
  double bestfit_sinsq2theta23 = -1;
  double bestfit_deltacp       = -1;
  double bestfit_dm12sq        = -1;
  double bestfit_dm23sq        = -1;
  double bestfit_chisq         = -1;

  // loop over toy-MC experiments
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // re-seed fitter (reset current values to specified init values)
     fparm->Reset();

     //
     // get superk spectrum observation (toy experiment)
     //
     SuperkObservation * skobs = new SuperkObservation;
     skobs->Read(xmlfile,iexperiment);

     // get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     NuMuDisappFitter * fitter = new NuMuDisappFitter(fit_inputs,fitter_name);

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
                  bestfit_deltacp, bestfit_dm12sq, bestfit_dm23sq, bestfit_chisq, fitstatus);

     delete skobs;
     delete fit_inputs;
     delete fitter;
  }

  //
  // store results
  //
  TFile f(Form("%s-run_%d-expts_%d_%d.root",
               tag, run, first_experiment, last_experiment), "recreate");
  nfitnt->Write("nfitnt");
  f.Close();
}
