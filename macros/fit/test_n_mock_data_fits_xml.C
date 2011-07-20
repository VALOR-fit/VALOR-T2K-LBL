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
  cout << "Fitter : " << fitter_name << endl;
  cout << "Reading experiments from : " << xmlfile << endl; 
  cout << "Fitting experiments : [" << first_experiment << ", " << last_experiment << "]" << endl; 

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // fit param list
  FitParamList * fparm = FitParamList::NamedParameterList(fitparamlist);
  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // store results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","", "iexp:sinsq2theta23:dm23sq:chisq");

  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // seed fitter
     if(hyp == conv::kHypOsc3F) {
       fparm->Update(conv::kPhypSinsq2Theta12, conv::kGlobFitSinsq2Theta12);
       fparm->Update(conv::kPhypSinsq2Theta13, conv::kGlobFitSinsq2Theta13);
       fparm->Update(conv::kPhypSinsq2Theta23, conv::kGlobFitSinsq2Theta23);  
       fparm->Update(conv::kPhypDeltaCP,       0.);
       fparm->Update(conv::kPhypDm12sq,        conv::kGlobFitDm12sq);
       fparm->Update(conv::kPhypDm23sq,        conv::kGlobFitDm23sq);
     }
     else
     if(hyp == conv::kHypNuMuDisapp2F) {
       fparm->Update(conv::kPhypSinsq2Theta23, conv::kGlobFitSinsq2Theta23);
       fparm->Update(conv::kPhypDm23sq,        conv::kGlobFitDm23sq);
     }

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

     double bestfit_sinsq2theta23 = fparm->CurrValue(t2kana::conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(t2kana::conv::kPhypDm23sq);
     double bestfit_chisq         = fit_inputs->Chisquare();

     nfitnt->Fill(iexperiment, bestfit_sinsq2theta23, bestfit_dm23sq, bestfit_chisq);

     cout << "Best-fit:"
          << " sinsq2theta23 = " << bestfit_sinsq2theta23
          << ", dm23sq = " << bestfit_dm23sq << " eV^2/c^4"
          << " ==> chisq = " << bestfit_chisq << endl;

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
