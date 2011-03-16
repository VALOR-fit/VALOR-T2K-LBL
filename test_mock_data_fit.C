//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits for mock-data 
//
// Inputs:
//      - input_sinsq2theta23 :
//      - input_dm23sq        :
//      - fluctuate_stat      : include stat. fluctautions?
//      - fitter name         : oscillation fitter <"minuit", "custom" (2-D grid scan)>
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_mock_data_fit
(
   double       input_sinsq2theta23 = 0.9,
   double       input_dm23sq        = 2.4E-3,
   bool         fluctuate_stat      = true,
   const char * fitter_name         = "custom"
)
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

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->CreateMockData(*skexp, *fparm, fluctuate_stat);

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs = 
      new OscillationFitter::FitInputs(skobs,skexp,fparm);

  // get fitter
  OscillationFitter * fitter = new OscillationFitter(fit_inputs,fitter_name);
  fitter->ActualOscParams(input_sinsq2theta23, input_dm23sq);

  // do the fit
  fitter->DoTheFit();

  // draw results
  fitter->Draw(new TCanvas(), "StandardOscFitPlots+Minos90CL","mockfit");
  fitter->Draw(new TCanvas(), "ChisqSurface");
  fitter->Draw(new TCanvas(), "MarginalizedChisq");

  // Print-out mock-data (Ereco values) as a comma separated list.
  // Assuming uniform distribution within each bin, randomize values within each bin.
  TH1D* mock_data = skobs->ErecoSpectrum(conv::kSKsFCFV1RingMuLike);
  TRandom rand(0);
  bool is_first=true;
  for(int ibin = 1; ibin <= mock_data->GetNbinsX(); ibin++) {
    double Nb   = mock_data->GetBinContent (ibin); // number of events in this Ereco bin
    double Emin = mock_data->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
    double Emax = mock_data->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
    for(int iev = 0; iev < Nb; iev++) {
      double E = Emin + rand.Uniform() * (Emax-Emin);
      string delimeter = ",";
      if(is_first) {
        delimeter = "";
        is_first = false;
      }
      cout << delimeter << E;
    }//iev
  }//ibin
  cout << endl;

}
