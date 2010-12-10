//
// ------------------------------------------------------------------------------------------------
// Calc contour with systematics
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void run_syst_contour_calc(int run=0)
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

  //
  // ntuple & histograms to store results
  //
  const int nexperiments = 3;
  TH2D * chisq_surf[nexperiments];
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
  skobs->CreateMockData(*skexp, *fparm, false);     ///<-- generate mock data now
  //skobs->Read("inputs/analysis_configuration.xml"); ///<-- read mock data stored in file (do that if you run multiple jobs for the same systematic)

  // create fit inputs
  OscillationFitter::FitInputs * fit_inputs =
     new OscillationFitter::FitInputs(skobs,skexp,fparm);
  
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
     cout << " ******* experiment = " << iexperiment << endl;

     // tweak source of systematic err
     skexp->TweakAllSystematics( (iexperiment!=0) );

     // reset oscillation params following last fit
     fparm->Update(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(conv::kPhypDm23sq,        input_dm23sq);

     // get fitter
     OscillationFitter fitter(fit_inputs,"custom");

     // do the fit
     fitter.DoTheFit();

     // store params
     double bestfit_sinsq2theta23 = fparm->CurrValue(conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(conv::kPhypDm23sq);
     double bestfit_chisq         = fit_inputs->Chisquare();
     systnt->Fill(iexperiment, bestfit_sinsq2theta23, bestfit_dm23sq, bestfit_chisq);

     chisq_surf[iexperiment] = new TH2D(*fitter.GetChisqSurface());
     chisq_surf[iexperiment] -> SetName(Form("chisqsurf_%d",iexperiment));
  }

  //
  // average chisq over experiments in this run to make later
  // histogram processing easier
  //
  TH2D * chisq_surf_avg = new TH2D(*chisq_surf[0]);
  chisq_surf_avg->SetName("chisqsurf_avg");
  double chisq_avg_min = 99999999;
  for(int i = 1; i <= chisq_surf_avg->GetNbinsX(); i++) {
      for(int j = 1; j <= chisq_surf_avg->GetNbinsY(); j++) {
          double chisq_avg = 0;
          for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
              double chisq = chisq_surf[iexperiment]->GetBinContent(i,j);
              chisq_avg += chisq;
          }
          chisq_avg /= nexperiments;
          chisq_surf_avg->SetBinContent(i,j,chisq_avg);
          chisq_avg_min = TMath::Min(chisq_avg_min, chisq_avg);
      }//j
   }//i
  
   //
   // write-out results
   //
   TFile f(Form("contour-run_%d-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",
               run,input_sinsq2theta23,input_dm23sq), "recreate"); 

   chisq_surf_avg->Write();
   for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
	chisq_surf[iexperiment]->Write();
   }
   systnt->Write("sysnt");
   f.Close();
}
