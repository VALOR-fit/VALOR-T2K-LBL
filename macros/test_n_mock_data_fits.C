//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are generated on the fly.
//
// Inputs:
//      - nexperiments, 
//      - input_sinsq2theta23,
//      - input_dm23sq,
//      - fitter name ("minuit", "custom" (2-D grid scan))
//      - include_systematics (include systematic variations in generated toy experiments)
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_n_mock_data_fits(
  const char * tag = "toyexptfits",
  int run=0,
  int nexperiments=1000, 
  double input_sinsq2theta23 = 1.0,
  double input_dm23sq = 2.4E-3,
  const char * fitter_name="minuit", 
  bool include_systematics=false)
{
  using namespace std;
  using namespace t2kana;

  cout << "[Job configuration]" << endl;
  cout << "Fitter : " << fitter_name << endl;
  cout << "Systematics : " << ((include_systematics) ? "YES " : "NO") << endl;
  cout << "Number of toy experiments : " << nexperiments << endl; 
  cout << "True oscillation parameters : (sinsq2theta23 = " 
       << input_sinsq2theta23 << ", dm23sq = " << input_dm23sq << " eV^2/c^4" << ")" << endl;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
//skexp->Read("inputs/analysis_configuration.xml");
  skexp->Read( Form("%s/inputs/analysis_configuration.xml", gSystem->GetEnv("T2KOSCANA")));

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
        bool fluctuate = true;
        skobs->CreateMockData(*skexp, *fparm, fluctuate);
        // reset all systematics so that the PDF used in the fitter is the nominal one
        skexp->ResetSystematicVariations();
     } else {       
        // generate mock data using the current PDF (without systematic variations) 
        // with statistical fluctuations (Poisson fluctuation in each bin)
        bool fluctuate = true;
        skobs->CreateMockData(*skexp, *fparm, fluctuate); 
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
  TFile f(Form("%s-sinsq2theta23_%.2f-dm23sq_%.5feV2-run_%d.root",
          tag, input_sinsq2theta23,input_dm23sq,run), "recreate");
  nfitnt->Write("nfitnt");
  f.Close();
}
