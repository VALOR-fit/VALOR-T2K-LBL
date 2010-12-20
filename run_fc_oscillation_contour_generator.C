//
// ------------------------------------------------------------------------------------------------
// Run jobs needed to calculate the oscillation contour using the Feldman-Cousins method.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void run_fc_oscillation_contour_generator(
  int input_dm23sq_bin, int input_sinsq2theta23_bin, int run=0, string tag="fc")
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

  // number of pseudo-experiment fits per job
  const int nexperiments = 200;

  //
  // The oscillation contour will be constructed using a (dm23sq, sinsq2theta23) 
  // 2-D histogram with the following binning.
  // Look-up the bin centre values for the input_dm23sq_bin and input_sinsq2theta23_bin.
  //

  // dm23sq, sinsq2theta23 range and step
  double sinsq2theta23_min   = 0.0;
  double sinsq2theta23_max   = 1.0;
  double sinsq2theta23_step  = 0.01;
  double dm23sq_min          = 1E-3;
  double dm23sq_max          = 6E-3;
  double dm23sq_step         = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
  int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );
  // param space
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       sinsq2theta23_nbins, 
       sinsq2theta23_min - 0.5*sinsq2theta23_step, 
       sinsq2theta23_max + 0.5*sinsq2theta23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  cout << "Oscillation param space: " << endl;
  cout << " - sinsq2theta23 : " << osc_param_space->GetXaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetXaxis()->GetXmin() << ", "
       << osc_param_space->GetXaxis()->GetXmax() << ")" << endl;
  cout << " - dm23sq        : " << osc_param_space->GetYaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetYaxis()->GetXmin() << ", "
       << osc_param_space->GetYaxis()->GetXmax() << ")" << endl;

  assert(input_sinsq2theta23_bin >= 1 && input_sinsq2theta23_bin <= sinsq2theta23_nbins);
  assert(input_dm23sq_bin        >= 1 && input_dm23sq_bin        <= dm23sq_nbins       );

  const double input_sinsq2theta23 = osc_param_space->GetXaxis()->GetBinCenter(input_sinsq2theta23_bin);
  const double input_dm23sq        = osc_param_space->GetYaxis()->GetBinCenter(input_dm23sq_bin);

  cout << "****** Will fit " << nexperiments << " pseudo-experiments with "
       << "Dm23sq_true = " << input_dm23sq << ", sinsq2theta23_true = " << input_sinsq2theta23 << endl;

  // ntuple to store results
  TNtupleD * systnt = new TNtupleD("sysnt","", 
     "iexp:run:input_sinsq2theta23_bin:input_dm23sq_bin:input_sinsq2theta23:input_dm23sq:sinsq2theta23:dm23sq:chisq");

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
  // good init value so that fitte converges fast to correct oscillation params
  fparm->SetInitValue(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->SetInitValue(conv::kPhypDm23sq,        input_dm23sq       );
  
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // tweak sources of systematic err
     skexp->TweakAllSystematics();

     // reset oscillation params following last fit
     fparm->Update(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(conv::kPhypDm23sq,        input_dm23sq       );
     
     // get fitter
     OscillationFitter fitter(fit_inputs,"minuit");

     // do the fit
     fitter.DoTheFit();

     // get best-fit params
     double bestfit_sinsq2theta23 = fparm->CurrValue(conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(conv::kPhypDm23sq);
     double bestfit_chisq         = fit_inputs->Chisquare();

     // store results for this pseudo-experiment
     systnt->Fill(iexperiment, 
                  run, 
                  input_sinsq2theta23_bin,
                  input_dm23sq_bin,
                  input_sinsq2theta23,
                  input_dm23sq,
                  bestfit_sinsq2theta23, 
                  bestfit_dm23sq, 
                  bestfit_chisq);
  }//iexperiment

  // write-out results
  TFile f(Form("%s-run_%d-sinsq2theta23bin_%d-dm23sqbin_%d.root",
               tag.c_str(), run,input_sinsq2theta23_bin,input_dm23sq_bin), "recreate"); 

  systnt->Write("sysnt");
  f.Close();
}
