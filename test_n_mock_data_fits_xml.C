//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are read from the input XML file.
//
// Inputs:
//      - fitter name ("minuit", "custom" (2-D grid scan))
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_n_mock_data_fits_xml(
  const char * xmlfile,
  int first_experiment, 
  int last_experiment, 
  double input_sinsq2theta23 = 1.0,
  double input_dm23sq = 2.4E-3,
  const char * fitter_name="minuit")
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

  cout << "Job configuration" << endl;
  cout << "Fitter : " << fitter_name << endl;
  cout << "Reading experiments from : " << xmlfile << endl; 
  cout << "Fitting experiments : [" << first_experiment << ", " << last_experiment << "]" << endl; 
  cout << "True oscillation parameters : (sinsq2theta23 = " 
       << input_sinsq2theta23 << ", dm23sq = " << input_dm23sq << " eV^2/c^4" << ")" << endl;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // store results
  TNtupleD * nfitnt = new TNtupleD("nfitnt","", "iexp:sinsq2theta23:dm23sq:chisq");

  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     fparm->Update(t2kana::numu10a::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(t2kana::numu10a::conv::kPhypDm23sq,        input_dm23sq);

     //
     // get superk spectrum observation (toy experiment)
     //
     SuperkObservation * skobs = new SuperkObservation;
     skobs->Read(xmlfile,iexperiment);

     // get fitter
     OscillationFitter::FitInputs * fit_inputs = new OscillationFitter::FitInputs(skobs,skexp,fparm);
     OscillationFitter * fitter = new OscillationFitter(fit_inputs,fitter_name);
     fitter->ActualOscParams(input_sinsq2theta23, input_dm23sq);

     // do the fit
     fitter->DoTheFit();

     double bestfit_sinsq2theta23 = fparm->CurrValue(t2kana::numu10a::conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(t2kana::numu10a::conv::kPhypDm23sq);
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
  TFile f(Form("Nfits-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",
          input_sinsq2theta23,input_dm23sq), "recreate");
  nfitnt->Write("nfitnt");
  f.Close();
}
