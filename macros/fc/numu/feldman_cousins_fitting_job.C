//
// ------------------------------------------------------------------------------------------------
// Run a job with some of the toy expts fits needed to calculate the \Delta\chi^{2} values
// for constructing a Feldman-Cousins contour
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void feldman_cousins_fitting_job
(
  const char * tag                   = "fc",                         /* */ 
  int          irun                  = 0,                            /* */
  int          nexperiments          = 1000,                         /* */
  int          min_dm23sq_bin        =   1,                          /* */
  int          max_dm23sq_bin        =  51,                          /* */
  int          min_mixing23_bin      =   1,                          /* */
  int          max_mixing23_bin      = 101,                          /* */
  const char * config_file           = "analysis_configuration.xml", /* analysis configuration file */
  bool         stat_fluctuations     = true,                         /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics     = true,                         /* generate toy MC with randomized systematics? */
  bool         tweak_non23oscparams  = true,                         /* generate toy MC with randomized dmsq12,theta12,theta13,dcp / only in 3-flavour fmwk */
  const char * fit_param_list        = "numudisapp2f",               /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option            = "minuit"                      /* fitter option (see run_fitter.C) */
)
{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  //
  // The oscillation contour will be constructed using a (dm23sq, mixing23) 
  // 2-D histogram with the following binning.
  // Look-up the bin centre values for the input_dm23sq_bin and input_mixing23_bin.
  //

  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.0;
  double mixing23_max   = 1.0;
  double mixing23_step  = 0.01;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 6E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + (int) ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  // param space
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
   
  //
  // Print out job summary
  //
  cout << "[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << irun << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy MC experiments: " << nexperiments << endl;
  cout << "   - Include statistical fluctuations: " << ((stat_fluctuations) ? "YES" : "NO") << endl;
  cout << "   - Randomizing systematics: " << ((tweak_systematics) ? "YES" : "NO") << endl;
  cout << "   - Randomizing non-'23'-sector oscillation params: " << ((tweak_non23oscparams) ? "YES" : "NO") << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "- Oscillation parameter space binning" << endl;
  cout << "   - mixing23 : " << osc_param_space->GetXaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetXaxis()->GetXmin() << ", "
       << osc_param_space->GetXaxis()->GetXmax() << ")" << endl;
  cout << "   - dm23sq   : " << osc_param_space->GetYaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetYaxis()->GetXmin() << ", "
       << osc_param_space->GetYaxis()->GetXmax() << ") eV^2/c^4" << endl;
  cout << "- Input bin range " << endl;
  cout << "   - mixing23 : [" << min_mixing23_bin << ", " << max_mixing23_bin << "]" << endl;
  cout << "   - dm23sq   : [" << min_dm23sq_bin   << ", " << max_dm23sq_bin   << "]" << endl;

  // get superk spectrum PDF & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  TRandom3 rand(0);

  // ntuple to store results
  TNtupleD * fcnt = new TNtupleD("fcnt","", 
     "iexp:irun:input_mixing23_bin:input_dm23sq_bin:input_mixing23:input_dm23sq:chisq_for_input_params:mixing23:dm23sq:chisq");
  
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
     cout << "\n\n** Pseudo-experiment = " << iexperiment << endl;

     //
     // Get the specified fit parameter list
     //
     FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
     conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
     if(iexperiment==0) {
        cout << "Physics hypothesis for current fitting job: " << util::AsString(hyp) << endl;
     }

     //     
     // get dm23sq and mixing23 at the centre of a bin randomly selected within the input range
     //     

     int input_mixing23_bin = min_mixing23_bin + rand.Integer(1 + max_mixing23_bin - min_mixing23_bin);
     int input_dm23sq_bin   = min_dm23sq_bin   + rand.Integer(1 + max_dm23sq_bin   - min_dm23sq_bin  );
     input_mixing23_bin = TMath::Max(1,              input_mixing23_bin);
     input_mixing23_bin = TMath::Min(mixing23_nbins, input_mixing23_bin);
     input_dm23sq_bin   = TMath::Max(1,              input_dm23sq_bin  );
     input_dm23sq_bin   = TMath::Min(dm23sq_nbins,   input_dm23sq_bin  );

     double input_mixing23 = osc_param_space->GetXaxis()->GetBinCenter(input_mixing23_bin);
     double input_dm23sq        = osc_param_space->GetYaxis()->GetBinCenter(input_dm23sq_bin       );
     input_mixing23 = TMath::Max(0.,input_mixing23); // protection from rnd-err
     input_mixing23 = TMath::Min(1.,input_mixing23); // protection from rnd-err

     cout << "- Toy-MC > Dm23sq_{true} = " 
          << input_dm23sq << " eV^2, mixing23_{true} = " << input_mixing23 << endl;

     //
     // if we are working in 3-flavour framework, randomize mixing12,mixing13,deltacp and dm12sq
     // (these parameters are marginalised in the 3-flavour numu disappearance analysis)
     //

     double input_mixing12 = -1;
     double input_mixing13 = -1;
     double input_deltacp  = -1;
     double input_dm12sq   = -1;
     if(hyp == conv::kHypOsc3F) {
         input_mixing12 = util::RandomizeMixing12(); 
         input_mixing13 = util::RandomizeMixing13(); 
         input_deltacp       = util::RandomizeDeltaCP();
         input_dm12sq        = util::RandomizeDm12sq(); 
         cout << "- Toy-MC > Randomized non-'23'-sector params for current toy-MC: " 
              << "Dm12sq_{true} = " << input_dm12sq  << " eV^2, mixing12_{true} = " << input_mixing12 
              << ", d_cp_{true} = " << input_deltacp      << ", mixing13_{true} = " << input_mixing13 
              << endl;
     }

     //
     // Set the generated oscillation parameter values
     //
     fparm->Update(t2kana::conv::kPhypMixing23, input_mixing23);
     fparm->Update(t2kana::conv::kPhypDm23sq,   input_dm23sq);
     if(hyp == conv::kHypOsc3F) {
       fparm->Update(t2kana::conv::kPhypMixing12, input_mixing12);
       fparm->Update(t2kana::conv::kPhypMixing13, input_mixing13);
       fparm->Update(t2kana::conv::kPhypDeltaCP,  input_deltacp);
       fparm->Update(t2kana::conv::kPhypDm12sq,   input_dm12sq);
     }

     //
     // get superk spectrum observation (toy experiment)
     //
     SuperkObservation * skobs = new SuperkObservation;
     if(tweak_systematics) {   
        // tweak all sources of systematic err (gaussian err)
        cout << "- Randomizing all systematic params in pdf..." << endl; 
        skexp->TweakAllSystematics();
        skexp->PrintSystematicVariations();
        // generate mock data using the current PDF (incl. systematic variations)
        cout << "- Generating toy-MC..." << endl; 
        skobs->CreateMockData(*skexp, *fparm, stat_fluctuations);
        // reset all systematics so that the PDF used in the fitter is the nominal one
        cout << "- Reseting all systematic params in pdf..." << endl; 
        skexp->ResetSystematicVariations();
        //skexp->PrintSystematicVariations();
     } else {
        // generate mock data using the current PDF (without systematic variations)
        cout << "- Generating toy-MC..." << endl; 
        skobs->CreateMockData(*skexp, *fparm, stat_fluctuations);
     }

     // create fit inputs
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

     // get chisq for true (input) "23"-sector point (sinsq2theta,dmsq) 
     // (marginalised non-"23"-sector params are treared as systematics and in the chisq_true calculation
     // they are set to nominal (constant) values used in the fit).
     fparm->ResetCurrToInitValues();
     fparm->Update(t2kana::conv::kPhypMixing23, input_mixing23);
     fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);
     double chisq_for_input_params = fit_inputs->Chisquare();
     cout << "- Chisq (input/true oscillation params) = " << chisq_for_input_params << endl;

     // get the specified fitter and perform the fit
     fparm->ResetCurrToInitValues();
     Superk1RmuErecoFitter fitter(fit_inputs,fit_option);
     fitter.DoTheFit();

     // get best-fit params
     double bestfit_mixing23 = fparm->BestFitValue(conv::kPhypMixing23);
     double bestfit_dm23sq   = fparm->BestFitValue(conv::kPhypDm23sq);
     double bestfit_chisq    = fit_inputs->BestFitChisquare();

     cout << "- Chisq (Best-fit params:" 
          << " mixing23 = " << bestfit_mixing23
          << ", dm23sq = " << bestfit_dm23sq << " eV^2/c^4"
          << ") = " << bestfit_chisq << endl;

     cout << "- DeltaChisq = " << chisq_for_input_params - bestfit_chisq << endl;

     // store results for this pseudo-experiment
     fcnt->Fill(iexperiment, 
                irun, 
                input_mixing23_bin,
                input_dm23sq_bin,
                input_mixing23,
                input_dm23sq,
                chisq_for_input_params,
                bestfit_mixing23, 
                bestfit_dm23sq, 
                bestfit_chisq);

     // clean-up
     delete fit_inputs;
     delete fparm;
     delete skobs;

  }//iexperiment

  // write-out results
  TFile f(Form("%s-run_%d.root",tag, irun), "recreate"); 
  fcnt->Write("fcnt");
  f.Close();
}
