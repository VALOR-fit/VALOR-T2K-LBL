//
// ------------------------------------------------------------------------------------------------
// Print chi^2 for the given dataset using the prediction with the input oscillation params
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
#include "Conventions.h"

void print_dataset_chisq
(
  double       input_mixing12      = 0.86,                         /* */
  double       input_mixing13      = 0.10,                         /* */
  double       input_mixing23      = 1.00,                         /* */
  double       input_deltacp       = 0.00,                         /* */
  double       input_dm12sq        = 7.58E-5,                      /* */
  double       input_dm23sq        = 2.32E-3,                      /* */
  const char * fit_param_list      = "numudisapp2f",               /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  int          dataset_id          = 0,                            /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file           = "inputs/data.xml",            /* */
  const char * config_file         = "analysis_configuration.xml"  /* analysis configuration file */
)
{
  using namespace std;
  using namespace t2kana;

  // get superk spectrum predictor & configure
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // set osc params
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "Specified physics hypothesis: " << util::AsString(hyp) << endl;
  if(hyp == conv::kHypOsc3F) {
    fparm->Update(t2kana::conv::kPhypMixing12, input_mixing12);
    fparm->Update(t2kana::conv::kPhypMixing13, input_mixing13);
    fparm->Update(t2kana::conv::kPhypMixing23, input_mixing23);
    fparm->Update(t2kana::conv::kPhypDeltaCP,  input_deltacp);
    fparm->Update(t2kana::conv::kPhypDm12sq,   input_dm12sq);
    fparm->Update(t2kana::conv::kPhypDm23sq,   input_dm23sq);
  }
  else
  if(hyp == conv::kHypNuMuDisapp2F) {
    fparm->Update(t2kana::conv::kPhypMixing23, input_mixing23);
    fparm->Update(t2kana::conv::kPhypDm23sq,   input_dm23sq);
  }
  else {
    return;
  }

  // generate prediction
  skexp->Generate(*fparm);
//TH1D * hExp = new TH1D(*skexp->CurrPredictedErecoSpectrum());

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);
//TH1D * hObs = new TH1D(*skobs->ErecoSpectrum());

  // get fit inputs object (chisq is calculated there)
  FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
  double chisq = fit_inputs->Chisquare();

  int ndof = 198;

  // report
  fparm->ShowOnlyCurrValuesInNextPrint();
  cout << *fparm;
  cout << "chisq      = " << chisq       << endl;
  cout << "ndof       = " << ndof        << endl;
  cout << "chisq/ndof = " << chisq/ndof  << endl;
}
