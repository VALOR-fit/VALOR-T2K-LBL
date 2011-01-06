//
// ------------------------------------------------------------------------------------------------
// Run N numu disappearance fits and store best fit values
//
// Inputs:
//      - fitter name ("minuit", "custom" (2-D grid scan))
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_oscillation_contour_coverage(const char * fitter_name="minuit")
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

  gSystem->Load("libT2KDisapp2010a");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::numu10a;

  const double input_sinsq2theta23 = 1.0;
  const double input_dm23sq        = 2.4E-3;

  const int nexperiments = 5000;

  double output_sinsq2theta23[nexperiments];
  double output_dm23sq       [nexperiments];

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs = 
      new OscillationFitter::FitInputs(skobs,skexp,fparm);

  // store results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","", "iexp:sinsq2theta23:dm23sq:chisq");

  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     fparm->Update(t2kana::numu10a::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(t2kana::numu10a::conv::kPhypDm23sq,        input_dm23sq);

     skobs->CreateMockData(*skexp, *fparm);

     // get fitter
     OscillationFitter * fitter = new OscillationFitter(fit_inputs,fitter_name);
     fitter->ActualOscParams(input_sinsq2theta23, input_dm23sq);

     // do the fit
     fitter->DoTheFit();

     output_sinsq2theta23[iexperiment] = fparm->CurrValue(t2kana::numu10a::conv::kPhypSinsq2Theta23);
     output_dm23sq       [iexperiment] = fparm->CurrValue(t2kana::numu10a::conv::kPhypDm23sq);

     double bestfit_chisq = fit_inputs->Chisquare();

     nfitnt->Fill(iexperiment, 
         output_sinsq2theta23[iexperiment], output_dm23sq[iexperiment] , bestfit_chisq);

     cout << "Best-fit:"
          << " sinsq2theta23 = " << output_sinsq2theta23[iexperiment]
          << ", dm23sq = " << output_dm23sq[iexperiment] << " eV^2/c^4"
          << " ==> chisq = " << bestfit_chisq << endl;
  }

  // show results
  TH1D * hp1 = new TH1D("hs1","",100,0,1);
  TH1D * hp2 = new TH1D("hp2","",100,1E-4,1E-2);
  TH2D * hp3 = new TH2D("hp3","",100,0,1,100,1E-4,1E-2);
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
    hp1->Fill( output_sinsq2theta23[iexperiment] );
    hp2->Fill( output_dm23sq[iexperiment]        );
    hp3->Fill( output_sinsq2theta23[iexperiment], output_dm23sq[iexperiment] );
  }

  TCanvas * c = new TCanvas("c","",20,20,700,700);
  c->Divide(1,3);
  c->cd(1);
  hp1->Draw();
  hp1->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  c->cd(2);
  hp2->Draw();
  hp2->GetXaxis()->SetTitle("|#Delta m^2_{23}|");
  c->cd(3);
  hp3->SetMarkerStyle(20);
  hp3->SetMarkerSize(1);
  hp3->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  hp3->GetYaxis()->SetTitle("|#Delta m^2_{23}|");
  hp3->Draw();

  TFile f(Form("Nfits-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",
          input_sinsq2theta23,input_dm23sq), "recreate");
  nfitnt->Write("nfitnt");
  f.Close();
}
