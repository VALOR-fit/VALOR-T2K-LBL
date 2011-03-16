//
// ------------------------------------------------------------------------------------------------
// Print chi^2 for the given dataset using the prediction with the input oscillation params
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

void print_dataset_chisq(
  double       input_sinsq2theta23 = 1.,
  double       input_dm23sq        = 2.4E-3,
  int          dataset_id          = 0, 
  const char * data_file           = "inputs/data.xml")
{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("libT2KNuOsc");

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

  // get fit inputs object (chisq is calculated there)
  OscillationFitter::FitInputs * fit_inputs =
          new OscillationFitter::FitInputs(skobs,skexp,fparm);

  double chisq = fit_inputs->Chisquare();
  cout << "Chisq(dm23sq = " << input_dm23sq << ", sinsq2theta23 = " << input_sinsq2theta23 << ") = " << chisq << endl;

}
