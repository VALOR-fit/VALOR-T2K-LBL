//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits for mock-data 
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void run_systematic_study()
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

  const double input_sinsq2theta23 = 0.97;
  const double input_dm23sq        = 2.5E-3;

  const int nexperiments = 100;
  double output_sinsq2theta23[nexperiments];
  double output_dm23sq       [nexperiments];

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->CreateMockData(*skexp, *fparm, false);

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs =
     new OscillationFitter::FitInputs(skobs,skexp,fparm);

  TRandom3 rand(0);

  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     skexp->TweakBeamNorm(rand.Gaus(0,1));

     fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

     // get fitter
     OscillationFitter * fitter = new OscillationFitter(fit_inputs,"minuit");

     // do the fit
     fitter->DoTheFit();

     output_sinsq2theta23[iexperiment] = fparm->CurrValue(t2kana::yr1numu::conv::kPhypSinsq2Theta23);
     output_dm23sq       [iexperiment] = fparm->CurrValue(t2kana::yr1numu::conv::kPhypDm23sq);

     delete fitter;
  }

  TH1D * hp1 = new TH1D("hs1","",100,0,1);
  TH1D * hp2 = new TH1D("hp2","",100,1E-4,1E-2);
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
    hp1->Fill( output_sinsq2theta23[iexperiment] );
    hp2->Fill( output_dm23sq[iexperiment]        );
  }

  TCanvas * c = new TCanvas("c","",20,20,700,700);
  c->Divide(1,2);
  c->cd(1);
  hp1->Draw();
  hp1->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  c->cd(2);
  hp2->Draw();
  hp2->GetXaxis()->SetTitle("|#Delta m^2_{23}|");

}
