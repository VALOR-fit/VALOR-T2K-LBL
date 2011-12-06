//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance systematic study for given oscillation parameters.
// Note: Currently this macro works only for the effective 2-flavour disappearance analysis,
// not the 3-flavour disappearance or the 3-flavour global (disappearance + appearance) analysis.
//	
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void bestfit_systematic_study_fitting_job
(
   string       systematic,                                          /* tweaked systematic param */
   int          irun                = 0,                             /* run # */
   double       input_mixing23      = 0.90,                          /* */
   double       input_dm23sq        = 2.32E-3,                       /* */
   int          nexperiments        = 500,                           /* number of toy MC expt fits in this run */
   const char * config_file         = "analysis_configuration.xml"   /* analysis configuration */
)
{
  using namespace std;
  using namespace t2kana;

  TRandom3 rand(0);

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // Global options
  //Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  //Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  //cout << *Opt::Instance() << endl;

  // Create ntuple to store results
  TNtupleD * systnt = new TNtupleD("sysnt","", 
      "irun:iexp:input_mixing23:input_dm23sq:mixing23:dm23sq:chisq:status");

  // Get superk spectrum predictor & configure
  SuperkExpectation skexp;
  skexp.Read(config);

  // Fit params
  FitParamList * fparm = FitParamList::NumuDisapFitNoSyst();

  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << " ******* experiment = " << iexperiment << endl;

     // reset oscillation params following last fit
     fparm->Update(conv::kPhypMixing23, input_mixing23);
     fparm->Update(conv::kPhypDm23sq,   input_dm23sq);

     // tweak systematic in spectrum predictor
     skexp.TweakSystematic(systematic, rand.Gaus(0.,1.));

     // generate mock data using tweaked spectrum predictor - don't include statistical fluctuations
     SuperkObservation skobs;
     skobs.CreateMockData(skexp, *fparm, false);       

     // reset systematic variations so that the nominal spectrum prediction is used for the fit
     skexp.ResetSystematicVariations();

     // get fitter
     FitInputs fit_inputs(&skobs, &skexp, fparm);
     NuMuDisappFitter fitter(&fit_inputs,"minuit");

     // do the fit
     bool ok = fitter.DoTheFit();

     // store params
     double bestfit_mixing23 = fparm->CurrValue(conv::kPhypMixing23);
     double bestfit_dm23sq   = fparm->CurrValue(conv::kPhypDm23sq);
     double bestfit_chisq    = fit_inputs.Chisquare();
     double fit_status       = (ok) ? +1 : -1;

     systnt->Fill(irun, 
                  iexperiment, 
                  input_mixing23, 
                  input_dm23sq,
                  bestfit_mixing23, 
                  bestfit_dm23sq, 
                  bestfit_chisq, 
                  fit_status);
  }

  TFile f(Form("syst_%s-mixing23_%.2f-dm23sq_%.5feV2_%d.root",
          systematic.c_str(),input_mixing23,input_dm23sq,irun), "recreate"); 
  systnt->Write("sysnt");
  f.Close();
}
