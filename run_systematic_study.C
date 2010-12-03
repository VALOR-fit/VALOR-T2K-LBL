//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance systematic study
//
// The input (systematic) string can be any of:
// 	- nd280_eff
//	- sk_energy_scale
//	- sk_ccqe_eff
//	- sk_ccnonqe_eff
//	- sk_nc_eff
//	- flux_norm
//	- flux_shape
//      - xsec_ccqe
//      - xsec_cc1pi_to_ccqe
//      - xsec_ccoth_to_ccqe
//      - xsec_nc_to_ccqe
//      - nuclear_effect
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void run_systematic_study(string systematic, int run=0)
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

  const int nexperiments = 250;
  double output_sinsq2theta23[nexperiments];
  double output_dm23sq       [nexperiments];

  TRandom3 rand(0);

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->CreateMockData(*skexp, *fparm, false);     ///<-- generate mock data now
  //skobs->Read("inputs/analysis_configuration.xml"); ///<-- read mock data stored in file (do that if you run multiple jobs for the same systematic)

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs =
     new OscillationFitter::FitInputs(skobs,skexp,fparm);

  // ntuple to store results
  TNtupleD * systnt = new TNtupleD("sysnt","", "iexp:sinsq2theta23:dm23sq:chisq");

  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << " ******* experiment = " << iexperiment << endl;

     //     
     // ND280 detector systematics
     //     
     if(systematic == "nd280_eff") {
	skexp->GetND280Inputs()->TweakEfficiency(rand.Gaus(0.,1.));
     } 

     //     
     // SuperK detector systematics
     //     
     else 
     if(systematic == "sk_energy_scale") {
        skexp->GetSuperkInputs()->TweakEnergyScale(rand.Gaus(0.,1.));
     } 
     else
     if(systematic == "sk_ccqe_eff") {
        skexp->GetSuperkInputs()->TweakEfficiency(
           conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE, rand.Gaus(0.,1.));
     } 
     else
     if(systematic == "sk_ccnonqe_eff") {
        skexp->GetSuperkInputs()->TweakEfficiency(
           conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi, rand.Gaus(0.,1.));
     } 
     else
     if(systematic == "sk_nc_eff") {
        skexp->GetSuperkInputs()->TweakEfficiency(
           conv::kSKsFCFV1RingMuLike, conv::kSKmNC, rand.Gaus(0.,1.));
     }

     //     
     // flux systematics
     //     
     else
     if(systematic == "flux_norm") {   
        skexp->GetBeamInputs()->TweakNorm(rand.Gaus(0.,1.));  
     }
     else
     if(systematic == "flux_shape") {   
        skexp->GetBeamInputs()->TweakShape();  
     }

     //     
     // cross section systematics
     //     
     else
     if(systematic == "xsec_ccqe") {
        skexp->GetGeneratorInputs()->TweakCCQEXSec(rand.Gaus(0.,1.));
     }
     else
     if(systematic == "xsec_cc1pi_to_ccqe") {
        skexp->GetGeneratorInputs()->TweakCC1piToCCQEXSecRatio(rand.Gaus(0.,1.));
     }
     else
     if(systematic == "xsec_ccoth_to_ccqe") {
        skexp->GetGeneratorInputs()->TweakCCOthToCCQEXSecRatio(rand.Gaus(0.,1.));
     }
     else
     if(systematic == "xsec_nc_to_ccqe") {
        skexp->GetGeneratorInputs()->TweakNCToCCQEXSecRatio(rand.Gaus(0.,1.));
     }
     else
     if(systematic == "nuclear_effect") {
        skexp->GetGeneratorInputs()->TweakNuclEffect(rand.Gaus(0.,1.));
     }

     else {
        cerr << "Unknown systematic parameter!" << endl;
	exit(1);
     }

     fparm->Update(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(conv::kPhypDm23sq,        input_dm23sq);

     // get fitter
     OscillationFitter fitter(fit_inputs,"minuit");

     // do the fit
     fitter.DoTheFit();

     double bestfit_sinsq2theta23 = fparm->CurrValue(conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(conv::kPhypDm23sq);
     double bestfit_chisq         = fit_inputs->Chisquare();

     output_sinsq2theta23[iexperiment] = bestfit_sinsq2theta23;
     output_dm23sq       [iexperiment] = bestfit_dm23sq;

     systnt->Fill(iexperiment, bestfit_sinsq2theta23, bestfit_dm23sq, bestfit_chisq);
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

  TFile f(Form("syst_%s-sinsq2theta23_%.2f-dm23sq_%.5feV2_%d.root",
          systematic.c_str(),input_sinsq2theta23,input_dm23sq,run), "recreate"); 
  systnt->Write("sysnt");
  f.Close();
}
