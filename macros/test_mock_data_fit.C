//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits for mock-data. 
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_mock_data_fit
(
   double       input_sinsq2theta23 = 1.0,                          /* */
   double       input_dm23sq        = 2.32E-3,                      /* */
   bool         fluctuate_stat      = true,                         /* include statistical fluctuations? */
   bool         include_systematics = false,                        /* include systematic variation in generated toy MC expt? */
   const char * fitter_name         = "custom",                     /* oscillation fitter <"minuit", "custom" (2-D grid scan)> */
   const char * config_file         = "analysis_configuration.xml"  /* analysis configuration file */
)
{
  using namespace std;
  using namespace t2kana;

  cout << "[Job inputs]" << endl;
  cout << "Fitter : " << fitter_name << endl;
  cout << "Systematics : " << ((include_systematics) ? "YES " : "NO") << endl;
  cout << "Statistical fluctuations : " << ((fluctuate_stat) ? "YES " : "NO") << endl;
  cout << "True oscillation parameters : (sinsq2theta23 = "
       << input_sinsq2theta23 << ", dm23sq = " << input_dm23sq << " eV^2/c^4" << ")" << endl;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // get superk spectrum observation
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
