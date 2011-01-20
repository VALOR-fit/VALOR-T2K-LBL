//
// ------------------------------------------------------------------------------------------------
// Run jobs needed to calculate the oscillation contour using the Feldman-Cousins method.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void run_feldman_cousins_deltachisq_calculator_rndbin(
  int irun=0, string tag="fc", bool include_systematics=true)
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

  // number of pseudo-experiment fits per job
  const int nexperiments = 3000;

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
       << osc_param_space->GetYaxis()->GetXmax() << ") eV^2/c^4" << endl;

  TRandom rand(0);

  // ntuple to store results
  TNtupleD * fcnt = new TNtupleD("fcnt","", 
     "iexp:irun:input_sinsq2theta23_bin:input_dm23sq_bin:input_sinsq2theta23:input_dm23sq:chisq_for_input_params:sinsq2theta23:dm23sq:chisq");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum PDF & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  cout << "****** Will fit " << nexperiments << " pseudo-experiments" << endl;
  
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     double rnd1 = rand.Rndm();
     double rnd2 = rand.Rndm();
     //double rnd1 = 0.2 + 0.8 * rand.Rndm();
     //double rnd2 = 0.8 * rand.Rndm();

     int input_sinsq2theta23_bin = 1 + TMath::Nint( (sinsq2theta23_nbins -1) * rnd1);
     int input_dm23sq_bin        = 1 + TMath::Nint( (dm23sq_nbins        -1) * rnd2);

     const double input_sinsq2theta23 = osc_param_space->GetXaxis()->GetBinCenter(input_sinsq2theta23_bin);
     const double input_dm23sq        = osc_param_space->GetYaxis()->GetBinCenter(input_dm23sq_bin);

     cout << "Dm23sq_true = " << input_dm23sq << ", sinsq2theta23_true = " << input_sinsq2theta23 << endl;

     // reset osc. params to desired input values (they now contain the best-(previous)fit values)
     fparm->Update(t2kana::numu10a::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(t2kana::numu10a::conv::kPhypDm23sq,        input_dm23sq);

     //
     // get superk spectrum observation (toy experiment)
     //
     SuperkObservation * skobs = new SuperkObservation;
     if(include_systematics) {
        // tweak all sources of systematic err (gaussian err)
        skexp->TweakAllSystematics();
        // generate mock data using the current PDF (incl. systematic variations)
        // + statistical fluctuations (Poisson fluctuation in each bin)
        bool fluctuate = true;
        skobs->CreateMockData(*skexp, *fparm, fluctuate);
        // reset all systematics so that the PDF used in the fitter is the nominal one
        skexp->ResetSystematicVariations();
     } else {
        // generate mock data using the current PDF (without systematic variations)
        // with statistical fluctuations (Poisson fluctuation in each bin)
        bool fluctuate = true;
        skobs->CreateMockData(*skexp, *fparm, fluctuate); 
     }
     
     // create fit inputs
     OscillationFitter::FitInputs * fit_inputs =
          new OscillationFitter::FitInputs(skobs,skexp,fparm);

     // get chisq for true (sinsq2theta,dmsq) point 
     double chisq_for_input_params = fit_inputs->Chisquare();
     cout << "Chisq(input/true oscillation params) = " << chisq_for_input_params << endl;

     // good init value so that fitte converges fast to correct oscillation params
     fparm->SetInitValue(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->SetInitValue(conv::kPhypDm23sq,        input_dm23sq       );

     // get fitter
     OscillationFitter fitter(fit_inputs,"minuit");

     // do the fit
     fitter.DoTheFit();

     // get best-fit params
     double bestfit_sinsq2theta23 = fparm->CurrValue(conv::kPhypSinsq2Theta23);
     double bestfit_dm23sq        = fparm->CurrValue(conv::kPhypDm23sq);
     double bestfit_chisq         = fit_inputs->Chisquare();

     cout << "Best-fit:" 
          << " sinsq2theta23 = " << bestfit_sinsq2theta23
          << ", dm23sq = " << bestfit_dm23sq << " eV^2/c^4"
          << " => chisq = " << bestfit_chisq << endl;

     // store results for this pseudo-experiment
     fcnt->Fill(iexperiment, 
                irun, 
                input_sinsq2theta23_bin,
                input_dm23sq_bin,
                input_sinsq2theta23,
                input_dm23sq,
                chisq_for_input_params,
                bestfit_sinsq2theta23, 
                bestfit_dm23sq, 
                bestfit_chisq);

     // clean-up
     delete fit_inputs;
     delete skobs;

  }//iexperiment

  // write-out results
  TFile f(Form("%s-run_%d.root",tag.c_str(), irun), "recreate"); 

  fcnt->Write("fcnt");
  f.Close();
}
