//
// ------------------------------------------------------------------------------------------------
// Print chi^2 for the given dataset using the prediction with the input oscillation params
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
#include "Conventions.h"

void print_dataset_chisq(
  double       input_sinsq2theta23 = 1.,                            /* */
  double       input_dm23sq        = 2.32E-3,                       /* */
  const char * config_file         = "analysis_configuration.xml",  /* analysis configuration file */
  int          dataset_id          = 0,                             /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file           = "inputs/data.xml")             /* */
{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
  skexp->Generate(*fparm);
  TH1D * hExp = new TH1D(*skexp->CurrPredictedErecoSpectrum());

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);
  TH1D * hObs = new TH1D(*skobs->ErecoSpectrum());

  // get fit inputs object (chisq is calculated there)
  FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

  double chisq = fit_inputs->Chisquare();
  cout << "Chisq(dm23sq = " << input_dm23sq << ", sinsq2theta23 = " << input_sinsq2theta23 << ") = " << chisq << endl;

}
