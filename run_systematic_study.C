//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits for mock-data 
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void run_systematic_study(string systematic="superk/ccqe_efficiency")
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

  const int nexperiments = 10;
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

     double gtwk = rand.Gaus(0.,1.);

     if(systematic == "nd280/efficiency") {
	skexp->GetND280Inputs()->TweakEfficiency(gtwk);
     } 
     else 
     if(systematic == "superk/energy_scale") {
        skexp->GetSuperkInputs()->TweakEnergyScale(gtwk);
     } 
     else
     if(systematic == "superk/ccqe_efficiency") {
        skexp->GetSuperkInputs()->TweakEfficiency(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE, gtwk);
     } 
     else
     if(systematic == "superk/ccnonqe_efficiency") {
        skexp->GetSuperkInputs()->TweakEfficiency(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi, gtwk);
     } 
     else
     if(systematic == "superk/nc_efficiency") {
        skexp->GetSuperkInputs()->TweakEfficiency(conv::kSKsFCFV1RingMuLike, conv::kSKmNC, gtwk);
     }
     else
     if(systematic == "flux/norm") {   
        skexp->GetBeamInputs()->TweakNorm(gtwk);  
     }
     else
     if(systematic == "flux/shape") {   
	exit(1);
     }
     else
     if(systematic == "generator/ccqe_cross_section") {
        skexp->GetGeneratorInputs()->TweakSigmaModeSK(conv::kSKmNuMuCCQE, gtwk);
     }
     else
     if(systematic == "generator/ccnonqe_cross_section") {
        skexp->GetGeneratorInputs()->TweakSigmaModeSK(conv::kSKmNuMuCC1pi, gtwk);
     }
     else
     if(systematic == "generator/nc_cross_section") {
        skexp->GetGeneratorInputs()->TweakSigmaModeSK(conv::kSKmNC, gtwk);
     } 
     else {
	exit(1);
     }

     fparm->Update(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(conv::kPhypDm23sq,        input_dm23sq);

     // get fitter
     OscillationFitter fitter(fit_inputs,"minuit");

     // do the fit
     fitter.DoTheFit();

     output_sinsq2theta23[iexperiment] = fparm->CurrValue(conv::kPhypSinsq2Theta23);
     output_dm23sq       [iexperiment] = fparm->CurrValue(conv::kPhypDm23sq);
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
