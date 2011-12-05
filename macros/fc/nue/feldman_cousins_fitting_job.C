//
// ------------------------------------------------------------------------------------------------
// Run a job with some of the toy expts fits needed to calculate the \Delta\chi^{2} values
// for constructing a Feldman-Cousins contour.
// This is the nue appearance analysis version of the script.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void feldman_cousins_fitting_job
(
  const char * tag                   = "fc",                         /* */ 
  int          irun                  = 0,                            /* */
  int          nexperiments          = 1000,                         /* */
  int          min_deltacp_bin       =   1,                          /* */
  int          max_deltacp_bin       =  91,                          /* */
  int          min_mixing13_bin =   1,                          /* */
  int          max_mixing13_bin = 101,                          /* */
  const char * config_file           = "analysis_configuration.xml", /* analysis configuration file */
  bool         stat_fluctuations     = true,                         /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics     = true,                         /* generate toy MC with randomized systematics? */
  bool         tweak_oth_osc_params  = false,                        /* generate toy MC with randomized dm12sq,|dm23sq|,theta12,theta23 or fix them to global-best fit values?  */
  const char * fit_param_list        = "nuosc3f_nueapp_fix4",        /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option            = "minuit"                      /* fitter option (see run_fitter.C) */
)
{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  //
  // The oscillation contour will be constructed using a (deltacp, mixing13) 
  // 2-D histogram with the following binning.
  // Look-up the bin centre values for the input_deltacp_bin and input_mixing13_bin.
  //

  // deltacp, mixing13 range and step
  double mixing13_min   = 0.0;
  double mixing13_max   = 0.4;
  double mixing13_step  = 0.004;
  double deltacp_min         = -1 * conv::kpi;
  double deltacp_max         = conv::kpi;
  double deltacp_step        = 2*conv::kpi/90.;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing13_nbins = 1 + (int) ((mixing13_max - mixing13_min) / mixing13_step);
  int deltacp_nbins       = 1 + (int) ((deltacp_max       - deltacp_min      ) / deltacp_step      );
  // param space
  TH2D * osc_param_space = new TH2D("nue_osc_param_space", "",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);
   
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
  cout << "   - Randomizing other (than deltacp, theta13) oscillation params: " << ((tweak_oth_osc_params) ? "YES" : "NO") << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "- Oscillation parameter space binning" << endl;
  cout << "   - mixing13 : " << osc_param_space->GetXaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetXaxis()->GetXmin() << ", "
       << osc_param_space->GetXaxis()->GetXmax() << ")" << endl;
  cout << "   - deltacp       : " << osc_param_space->GetYaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetYaxis()->GetXmin() << ", "
       << osc_param_space->GetYaxis()->GetXmax() << ") eV^2/c^4" << endl;
  cout << "- Input bin range " << endl;
  cout << "   - mixing13 : [" << min_mixing13_bin << ", " << max_mixing13_bin << "]" << endl;
  cout << "   - deltacp       : [" << min_deltacp_bin       << ", " << max_deltacp_bin       << "]" << endl;

  // get superk spectrum PDF & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  TRandom3 rand(0);

  // ntuple to store results
  TNtupleD * fcnt = new TNtupleD("fcnt","", 
     "iexp:irun:input_mixing13_bin:input_deltacp_bin:input_mixing13:input_deltacp:chisq_for_input_params:mixing13:deltacp:chisq");
  
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
     cout << "\n\n** Pseudo-experiment = " << iexperiment << endl;

     //
     // Get the specified fit parameter list
     //
     FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
     conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
     assert(hyp == conv::kHypOsc3F);
     if(iexperiment==0) {
        cout << "Physics hypothesis for current fitting job: " << util::AsString(hyp) << endl;
     }

     //     
     // get deltacp and mixing13 at the centre of a bin randomly selected within the input range
     //     

     int input_mixing13_bin = min_mixing13_bin + rand.Integer(1 + max_mixing13_bin - min_mixing13_bin);
     int input_deltacp_bin       = min_deltacp_bin       + rand.Integer(1 + max_deltacp_bin       - min_deltacp_bin      );
     input_mixing13_bin = TMath::Max(1,                   input_mixing13_bin);
     input_mixing13_bin = TMath::Min(mixing13_nbins, input_mixing13_bin);
     input_deltacp_bin       = TMath::Max(1,                   input_deltacp_bin      );
     input_deltacp_bin       = TMath::Min(deltacp_nbins,       input_deltacp_bin      );

     double input_mixing13 = osc_param_space->GetXaxis()->GetBinCenter(input_mixing13_bin);
     double input_deltacp       = osc_param_space->GetYaxis()->GetBinCenter(input_deltacp_bin      );
     input_mixing13 = TMath::Max(0.,input_mixing13); // protection from rnd-err
     input_mixing13 = TMath::Min(1.,input_mixing13); // protection from rnd-err
     input_deltacp       = TMath::Max(-1*conv::kpi,input_deltacp); // protection from rnd-err
     input_deltacp       = TMath::Min(   conv::kpi,input_deltacp); // protection from rnd-err

     cout << "- Toy-MC > deltacp_{true} = " 
          << input_deltacp << ", mixing13_{true} = " << input_mixing13 << endl;

     fparm->Update(t2kana::conv::kPhypMixing13, input_mixing13);
     fparm->Update(t2kana::conv::kPhypDeltaCP,       input_deltacp);

     //
     // if requested, randomize mixing12,mixing23,|dm23sq| and dm12sq
     // (nue appearance analysis, these parameters are either marginalised or fixed)
     //

     if(tweak_oth_osc_params) 
     {
         double input_mixing12 = util::Randomizemixing12(); 
         double input_mixing23 = util::Randomizemixing23(); 
         double input_dm23sq        = util::RandomizeDm23sq();
         double input_dm12sq        = util::RandomizeDm12sq(); 

         cout << "- Toy-MC > Randomized oscillation params for current toy-MC: " 
              << "Dm12sq_{true} = "   << input_dm12sq  << " eV^2, mixing12_{true} = " << input_mixing12 
              << ", Dm23sq_{true} = " << input_dm23sq  << " eV^2, mixing23_{true} = " << input_mixing23 
              << endl;

         fparm->Update(t2kana::conv::kPhypMixing23, input_mixing23);
         fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);
         fparm->Update(t2kana::conv::kPhypMixing12, input_mixing12);
         fparm->Update(t2kana::conv::kPhypDm12sq,        input_dm12sq);
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

     // get chisq for true (input) oscillation point (deltacp, mixing13) 
     // (marginalised / other oscillation params are treared as systematics and in the chisq_true calculation
     // they are set to nominal (constant) values used in the fit).
     fparm->ResetCurrToInitValues();
     fparm->Update(t2kana::conv::kPhypMixing13, input_mixing13);
     fparm->Update(t2kana::conv::kPhypDeltaCP,       input_deltacp);
     double chisq_for_input_params = fit_inputs->Chisquare();
     cout << "- Chisq (input/true oscillation params) = " << chisq_for_input_params << endl;

     // get the specified fitter and perform the fit
     fparm->ResetCurrToInitValues();
     Superk1ReNFitter fitter(fit_inputs,fit_option);
     fitter.DoTheFit();

     // get best-fit params
     double bestfit_mixing13 = fparm->BestFitValue(conv::kPhypMixing13);
     double bestfit_deltacp       = fparm->BestFitValue(conv::kPhypDeltaCP);
     double bestfit_chisq         = fit_inputs->BestFitChisquare();

     cout << "- Chisq (Best-fit params:" 
          << " mixing13 = " << bestfit_mixing13
          << ", deltacp = " << bestfit_deltacp << ") = " << bestfit_chisq << endl;

     cout << "- DeltaChisq = " << chisq_for_input_params - bestfit_chisq << endl;

     // store results for this pseudo-experiment
     fcnt->Fill(iexperiment, 
                irun, 
                input_mixing13_bin,
                input_deltacp_bin,
                input_mixing13,
                input_deltacp,
                chisq_for_input_params,
                bestfit_mixing13, 
                bestfit_deltacp, 
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
