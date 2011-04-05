//
// ------------------------------------------------------------------------------------------------
// Calculate the probability of a dataset under the given oscillation parameter hypothesis
//
// Inputs: 
//	- sinsq2theta23
//	- dm23sq
//	- dataset_id (number in <data_set_%d> XML tag, see inputs/data.xml)
//      - data_file  
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
#include "Conventions.h"

void calculate_dataset_probability(
  double       input_sinsq2theta23 = 1.,
  double       input_dm23sq        = 2.4E-3,
  int          dataset_id          = 0, 
  const char * data_file           = "inputs/data.xml")
{
  using namespace std;
  using namespace t2kana;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");
  skexp->Generate(*fparm);
  TH1D * hExp = new TH1D(*skexp->CurrPredictedErecoSpectrum());

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);
  TH1D * hObs = new TH1D(*skobs->ErecoSpectrum());

  cout << fixed;
  double prob = 1.;

  for(int reco_bin = 1; reco_bin <= hExp->GetNbinsX(); reco_bin++) {
    double nexp = hExp->GetBinContent(reco_bin);
    double nobs = hObs->GetBinContent(reco_bin);
    double p    = TMath::Poisson(nobs, nexp);
    prob *= p;

    cout << "bin = " << setw(5)  << setfill(' ') << reco_bin 
         << ": nexp = " << setw(10) << setfill(' ') << setprecision(8) << nexp 
         << ", nobs = " << setw(6)  << setfill(' ') << setprecision(4) << nobs 
         << " ==> prob(bin) = " << setw(12) << setfill(' ') << setprecision(10) << p 
         << endl;
  }
  cout << scientific;
  cout << "prob(dataset) = " << setprecision(10) << prob << endl;
}
