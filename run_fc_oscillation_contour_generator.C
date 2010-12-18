//
// ------------------------------------------------------------------------------------------------
// Run jobs needed to calculate the oscillation contour using the Feldman-Cousins method.
//
//
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void run_fc_oscillation_contour_generator(int run=0, string tag="fc")
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

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  const double input_sinsq2theta23 = 0.90;
  const double input_dm23sq        = 2.4E-3;

  // number of pseudo-experiment fits per job
  const int nexperiments = 100;

  // ntuple to store results
  TNtupleD * systnt = new TNtupleD("sysnt","", "iexp:sinsq2theta23:dm23sq:chisq");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->CreateMockData(*skexp, *fparm, true);  ///<-- generate mock data now (+ fluctuate)

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs =
     new OscillationFitter::FitInputs(skobs,skexp,fparm);
  
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
     cout << " ******* experiment = " << iexperiment << endl;

     // tweak sources of systematic err
     skexp->TweakAllSystematics();

     // reset oscillation params following last fit
     fparm->Update(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(conv::kPhypDm23sq,        input_dm23sq       );

     // get fitter
     OscillationFitter fitter(fit_inputs,"minuit");

     // do the fit
     fitter.DoTheFit();

     // store params
     double bestfit_sinsq2theta23 = fparm->CurrValue(conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(conv::kPhypDm23sq);
     double bestfit_chisq         = fit_inputs->Chisquare();
     systnt->Fill(iexperiment, bestfit_sinsq2theta23, bestfit_dm23sq, bestfit_chisq);
  }//iexperiment

  // write-out results
  TFile f(Form("%s-run_%d-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",
               tag.c_str(), run,input_sinsq2theta23,input_dm23sq), "recreate"); 

  systnt->Write("sysnt");
  f.Close();
}
