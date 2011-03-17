//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance systematic study for given oscillation parameters
//
// Inputs
// - systematic:
//    Any of
//	- sk_energy_scale
//	- sk_ccqe_eff
//	- sk_ccnonqe_eff
//	- sk_nc_eff
//	- sk_nuecceff_eff
// 	- nd280_eff
//	- flux_norm
//	- flux_shape
//      - xsec_ccqe
//      - xsec_cc1pi_to_ccqe
//      - xsec_ccoth_to_ccqe
//      - xsec_nc_to_ccqe
//      - nuclear_effect
// - irun :
// - input_sinsq2theta23 :
// - input_dm23sq :
// - nexperiments :
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void bestfit_systematic_study_fitting_job(
   string systematic, 
   int    irun                = 0, 
   double input_sinsq2theta23 = 0.90,
   double input_dm23sq        = 2.4E-3,
   int    nexperiments        = 500
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

  TRandom3 rand(0);

  // create ntuple to store results
  TNtupleD * systnt = new TNtupleD("sysnt","", 
      "irun:iexp:input_sinsq2theta23:input_dm23sq:sinsq2theta23:dm23sq:chisq:status");

  // get superk spectrum predictor & configure
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << " ******* experiment = " << iexperiment << endl;

     // reset oscillation params following last fit
     fparm->Update(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(conv::kPhypDm23sq,        input_dm23sq);

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
     double bestfit_sinsq2theta23 = fparm->CurrValue(conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(conv::kPhypDm23sq);
     double bestfit_chisq         = fit_inputs.Chisquare();
     double fit_status            = (ok) ? +1 : -1;

     systnt->Fill(irun, 
                  iexperiment, 
                  input_sinsq2theta23, 
                  input_dm23sq,
                  bestfit_sinsq2theta23, 
                  bestfit_dm23sq, 
                  bestfit_chisq, 
                  fit_status);
  }

  TFile f(Form("syst_%s-sinsq2theta23_%.2f-dm23sq_%.5feV2_%d.root",
          systematic.c_str(),input_sinsq2theta23,input_dm23sq,irun), "recreate"); 
  systnt->Write("sysnt");
  f.Close();
}
