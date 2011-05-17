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
  double       input_sinsq2theta23 = 1.0,                          /* */
  double       input_dm23sq        = 2.32E-3,                      /* */
  const char * fitter_name         = "minuit",                     /* fitter name ("minuit", "custom" (2-D grid scan)) */
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
  cout << "True oscillation parameters : (sinsq2theta23 = " 
       << input_sinsq2theta23 << ", dm23sq = " << input_dm23sq << " eV^2/c^4" << ")" << endl;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // store results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","", "iexp:sinsq2theta23:dm23sq:chisq");

  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

     //
     // get superk spectrum observation (toy experiment)
     //
     SuperkObservation * skobs = new SuperkObservation;
     skobs->Read(xmlfile,iexperiment);

     // get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     NuMuDisappFitter * fitter = new NuMuDisappFitter(fit_inputs,fitter_name);
     //fitter->ActualOscParams(input_sinsq2theta23, input_dm23sq);

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
