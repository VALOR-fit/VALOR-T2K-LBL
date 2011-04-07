//
// ------------------------------------------------------------------------------------------------
// Run a job with some of the toy expts fits needed to calculate the \Delta\chi^{2} values
// for constructing a Feldman-Cousins contour
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void feldman_cousins_fitting_job
(
  int          min_dm23sq_bin        =   1,                          /* */
  int          max_dm23sq_bin        =  51,                          /* */
  int          min_sinsq2theta23_bin =   1,                          /* */
  int          max_sinsq2theta23_bin = 101,                          /* */
  const char * tag                   = "fc",                         /* */ 
  int          irun                  = 0,                            /* */
  int          nexperiments          = 1000,                         /* */
  bool         include_systematics   = true,                         /* */
  const char * config_file           = "analysis_configuration.xml"  /* analysis configuration file */
)
{
  using namespace std;
  using namespace t2kana;

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
   
  //
  // Print out job summary
  //
  cout << "** Job: " << tag << "/" << irun << endl;
  cout << "** Number of toy MC experiments: " << nexperiments << endl;
  cout << "** Including systematics: " << ((include_systematics) ? "YES" : "NO") << endl;
  cout << "** Oscillation param space: " << endl;
  cout << " - sinsq2theta23 : " << osc_param_space->GetXaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetXaxis()->GetXmin() << ", "
       << osc_param_space->GetXaxis()->GetXmax() << ")" << endl;
  cout << " - dm23sq        : " << osc_param_space->GetYaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetYaxis()->GetXmin() << ", "
       << osc_param_space->GetYaxis()->GetXmax() << ") eV^2/c^4" << endl;
  cout << "** Input bin range: " << endl;
  cout << " - sinsq2theta23 : [" << min_sinsq2theta23_bin << ", " << max_sinsq2theta23_bin << "]" << endl;
  cout << " - dm23sq        : [" << min_dm23sq_bin        << ", " << max_dm23sq_bin        << "]" << endl;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // ntuple to store results
  TNtupleD * fcnt = new TNtupleD("fcnt","", 
     "iexp:irun:input_sinsq2theta23_bin:input_dm23sq_bin:input_sinsq2theta23:input_dm23sq:chisq_for_input_params:sinsq2theta23:dm23sq:chisq");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // get superk spectrum PDF & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  TRandom3 rand(0);
  
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     //
     // select a bin in the input range
     //
     int input_sinsq2theta23_bin = min_sinsq2theta23_bin + rand.Integer(1 + max_sinsq2theta23_bin - min_sinsq2theta23_bin);
     int input_dm23sq_bin        = min_dm23sq_bin        + rand.Integer(1 + max_dm23sq_bin        - min_dm23sq_bin       );

     input_sinsq2theta23_bin = TMath::Max(1,                   input_sinsq2theta23_bin);
     input_sinsq2theta23_bin = TMath::Min(sinsq2theta23_nbins, input_sinsq2theta23_bin);
     input_dm23sq_bin        = TMath::Max(1,                   input_dm23sq_bin       );
     input_dm23sq_bin        = TMath::Min(dm23sq_nbins,        input_dm23sq_bin       );

     //     
     // get dmsq and sinsq2theta at bin centre
     //     
     double input_sinsq2theta23 = osc_param_space->GetXaxis()->GetBinCenter(input_sinsq2theta23_bin);
     double input_dm23sq        = osc_param_space->GetYaxis()->GetBinCenter(input_dm23sq_bin       );

     input_sinsq2theta23 = TMath::Max(0.,input_sinsq2theta23); // protection from rnd-err
     input_sinsq2theta23 = TMath::Min(1.,input_sinsq2theta23); // protection from rnd-err

     cout << "* Dm23sq_{true} = " << input_dm23sq << ", sinsq2theta23_{true} = " << input_sinsq2theta23 << endl;

     //
     // reset osc. params to desired input values (they now contain the best-(previous)fit values)
     //
     fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

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
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

     // get chisq for true (sinsq2theta,dmsq) point 
     double chisq_for_input_params = fit_inputs->Chisquare();
     cout << "Chisq(input/true oscillation params) = " << chisq_for_input_params << endl;

     // good init value so that fitte converges fast to correct oscillation params
     fparm->SetInitValue(conv::kPhypSinsq2Theta23, input_sinsq2theta23);
     fparm->SetInitValue(conv::kPhypDm23sq,        input_dm23sq       );

     // get fitter
     NuMuDisappFitter fitter(fit_inputs,"minuit");

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
  TFile f(Form("%s-run_%d.root",tag, irun), "recreate"); 
  fcnt->Write("fcnt");
  f.Close();
}
