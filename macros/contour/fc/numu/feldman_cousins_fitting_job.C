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
  const char * tag                   = "fc",                         /* tag for output filename */ 
  int          irun                  = 0,                            /* run number */
  int          nexperiments          = 1000,                         /* number of toy MC expt fits in current run */
  int          min_dm23sq_bin        =   1,                          /* generate toy MC with this range of Dm^2_{23} */
  int          max_dm23sq_bin        =  51,                          /* -//- */
  int          min_mixing23_bin      =   1,                          /* generate toy MC with this range of 23 mixing (expressed in terms of either sinsq2theta, sinsqtheta, ...) */
  int          max_mixing23_bin      = 101,                          /* -//- */
  const char * bin_list              = "",                           /* if specified, run toy-MC expts only in the bins specified as long as they are in the range given above */
  const char * config_file           = "analysis_configuration.xml", /* analysis configuration file */
  bool         stat_fluctuations     = true,                         /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics     = true,                         /* generate toy MC with randomized systematics? */
  bool         tweak_non23oscparams  = true,                         /* generate toy MC with randomized dmsq12,theta12,theta13,dcp / only in 3-flavour fmwk */
  const char * fit_param_list        = "nuosc3f_numudisapp_fix4",    /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option            = "minuit"                      /* fitter option (see run_fitter.C) */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  cout << *Opt::Instance() << endl;

  //
  // The oscillation contour will be constructed using a (dm23sq, mixing23) 
  // 2-D histogram with the following binning.
  // Look-up the bin centre values for the input_dm23sq_bin and input_mixing23_bin.
  //

  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.5;
  double mixing23_max   = 1.0;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 6E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + TMath::CeilNint((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + TMath::CeilNint((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  // param space
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  //
  // If a bin list was specified, read the 2-D bin numbers
  //

  const int buffsize = 2048; 

  vector<int> mixing23_bin_list;
  vector<int> dm23sq_bin_list;

  bool use_bin_list        = false;
  bool bin_list_accessible = !gSystem->AccessPathName(bin_list);
  bool bin_list_input      = (strlen(bin_list) > 0);
  if(bin_list_input) {
    if(!bin_list_accessible) {
        cerr << "Input file: " << bin_list << " is not accessible." << endl;
        exit(1);
    }
    use_bin_list = true;
  }
  if(use_bin_list) {
     ifstream inpbins(bin_list, ios::in);
     if(!inpbins.is_open()) {
       cerr << "Can't open file: " << bin_list << endl;
       exit(1);
     } 
     // read all lines
     while(1) 
     {
        // ignore comment lines starting with #
        if(inpbins.peek() == '#') {
           inpbins.ignore(buffsize, '\n');
        } else {       
           int mixing23_bin = -1;
           int dm23sq_bin   = -1; 
           inpbins >> mixing23_bin >> dm23sq_bin;
           if(inpbins.eof()) break;
           mixing23_bin_list.push_back(mixing23_bin);
           dm23sq_bin_list.push_back(dm23sq_bin);
       }
     }//1 
  }//use bin list
  if(use_bin_list && mixing23_bin_list.size() == 0) {
     cerr << "Empty input list of bins... Discarding!" << endl; 
     use_bin_list = false;
  }
   
  //
  // Print out job summary
  //
  cout << "\n[Job inputs]" << endl;
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
  cout << "   - " << util::AsString(conv::kPhypMixing23) << " : " 
       << osc_param_space->GetXaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetXaxis()->GetXmin() << ", "
       << osc_param_space->GetXaxis()->GetXmax() << ")" << endl;
  cout << "   - " << util::AsString(conv::kPhypDm23sq) << " : " 
       << osc_param_space->GetYaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetYaxis()->GetXmin() << ", "
       << osc_param_space->GetYaxis()->GetXmax() << ") eV^2/c^4" << endl;
  cout << "- Input bin range " << endl;
  cout << "   - " << util::AsString(conv::kPhypMixing23) 
       << " : [" << min_mixing23_bin << ", " << max_mixing23_bin << "]" << endl;
  cout << "   - " << util::AsString(conv::kPhypDm23sq) 
       << " : [" << min_dm23sq_bin   << ", " << max_dm23sq_bin   << "]" << endl;
  cout << "- Input (" << util::AsString(conv::kPhypMixing23) << ", " << util::AsString(conv::kPhypDm23sq) << ") bin list " << endl; 
  if(!use_bin_list) {
    cout << "   - none" << endl;
  } else {
    const int nb = mixing23_bin_list.size();
    cout << "   - The following " << mixing23_bin_list.size() << " bins were specified: \n     ";
    for(int i = 0; i < nb; i++) {
        cout << "(" << mixing23_bin_list[i] << ", " << dm23sq_bin_list[i] << ")";
        if(i < nb-1) { 
          cout << ", "; 
        }
    }
    cout << endl;
  }

  // get superk spectrum PDF & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  TRandom3 rand(0);

  // Create n-tuple to store results
  TNtupleD * fcnt = new TNtupleD("fcnt","", 
     "iexp:irun:input_mixing23_bin:input_dm23sq_bin:input_mixing23:input_dm23sq:chisq_for_input_params:mixing23:dm23sq:chisq:status");
  
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
     cout << "\n\n** Pseudo-experiment = " << iexperiment << endl;

     // Get the specified fit parameter list
     FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
     conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
     if(iexperiment==0) {
        cout << "Physics hypothesis for current fitting job: " << util::AsString(hyp) << endl;
     }

     //     
     // Random selection of dm23sq and mixing23 bin
     //     
     
     int input_mixing23_bin = -1;
     int input_dm23sq_bin   = -1;
     
     if(!use_bin_list) {
       // If there was no input bin list, select a bin within the input range
       input_mixing23_bin = min_mixing23_bin + rand.Integer(1 + max_mixing23_bin - min_mixing23_bin);
       input_dm23sq_bin   = min_dm23sq_bin   + rand.Integer(1 + max_dm23sq_bin   - min_dm23sq_bin  );
     } else {
       // If a bin list was input, select a bin from that list and accept it if it is within the input range
       const int nb = mixing23_bin_list.size();
       while(1) {
           int i = rand.Integer(nb); // random in [0, nb-1] range
           input_mixing23_bin = mixing23_bin_list[i]; 
           input_dm23sq_bin   = dm23sq_bin_list[i];
           bool accept = false;
           if(input_mixing23_bin >= min_mixing23_bin &&
              input_mixing23_bin <= max_mixing23_bin &&
              input_dm23sq_bin   >= min_dm23sq_bin   &&
              input_dm23sq_bin   <= max_dm23sq_bin)
           {
             accept = true;
           }  
           if(accept) break; // exit inf loop
       }//1
     }
     // Sanity check on bin range
     input_mixing23_bin = TMath::Max(1,              input_mixing23_bin);
     input_mixing23_bin = TMath::Min(mixing23_nbins, input_mixing23_bin);
     input_dm23sq_bin   = TMath::Max(1,              input_dm23sq_bin  );
     input_dm23sq_bin   = TMath::Min(dm23sq_nbins,   input_dm23sq_bin  );

     cout << "- Toy-MC > Selected 23-sector params bin: " 
          << util::AsString(conv::kPhypMixing23) << " bin # = " << input_mixing23_bin << ", "
          << util::AsString(conv::kPhypDm23sq)   << " bin # = " << input_dm23sq_bin
          << endl;

     // Get dm23sq and mixing23 at the centre of the randomly selected bin
     double input_mixing23 = osc_param_space->GetXaxis()->GetBinCenter(input_mixing23_bin);
     double input_dm23sq   = osc_param_space->GetYaxis()->GetBinCenter(input_dm23sq_bin  );
     input_mixing23 = TMath::Max(0.,input_mixing23); // protection from rnd-err
     input_mixing23 = TMath::Min(1.,input_mixing23); // protection from rnd-err

     cout << "- Toy-MC > True 23-sector params: " 
          << util::AsString(conv::kPhypMixing23) << " = " << input_mixing23 << ", "
          << util::AsString(conv::kPhypDm23sq)   << " = " << input_dm23sq   << " eV^2"
          << endl;

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
         input_deltacp  = util::RandomizeDeltaCP();
         input_dm12sq   = util::RandomizeDm12sq(); 
         cout << "- Toy-MC > True (randomized) non-23-sector params: " 
              << util::AsString(conv::kPhypMixing12) << " = " << input_mixing12 << ", "
              << util::AsString(conv::kPhypDm12sq)   << " = " << input_dm12sq   << " eV^2, "
              << util::AsString(conv::kPhypMixing12) << " = " << input_mixing13 << ", "
              << util::AsString(conv::kPhypDeltaCP)  << " = " << input_deltacp  << " rad "
              << endl;
     }

     //
     // Set the generated oscillation parameter values
     //
     fparm->Update(conv::kPhypMixing23, input_mixing23);
     fparm->Update(conv::kPhypDm23sq,   input_dm23sq);
     if(hyp == conv::kHypOsc3F) {
       fparm->Update(conv::kPhypMixing12, input_mixing12);
       fparm->Update(conv::kPhypMixing13, input_mixing13);
       fparm->Update(conv::kPhypDeltaCP,  input_deltacp);
       fparm->Update(conv::kPhypDm12sq,   input_dm12sq);
     }

     //
     // get superk spectrum observation (toy experiment)
     //
     SuperkObservation * skobs = new SuperkObservation;
     //Use PhysHypothesis (instead of FitParamList) so that we don't need to Update the values of systematics with the values from TweakAllSystematics()
     PhysHypothesis physhyp(*fparm);
     if(tweak_systematics) {   
        // tweak all sources of systematic err (gaussian err)
        cout << "- Randomizing all systematic params in pdf..." << endl; 
        skexp->TweakAllSystematics();
        skexp->PrintSystematicVariations();
        // generate mock data using the current PDF (incl. systematic variations)
        cout << "- Generating toy-MC..." << endl; 
        skobs->CreateMockData(*skexp, physhyp, stat_fluctuations);
        // reset all systematics so that the PDF used in the fitter is the nominal one
        cout << "- Reseting all systematic params in pdf..." << endl; 
        skexp->ResetSystematicVariations();
        //skexp->PrintSystematicVariations();
     } else {
        // generate mock data using the current PDF (without systematic variations)
        cout << "- Generating toy-MC..." << endl; 
        skobs->CreateMockData(*skexp, physhyp, stat_fluctuations);
     }

     // get chisq for true (input) "23"-sector point (sinsq2theta,dmsq) 
     // (marginalised non-"23"-sector params are treared as systematics and in the chisq_true calculation
     // they are set to nominal (constant) values used in the fit).
     fparm->ResetCurrToInitValues();
     //help the fitter out by setting a new initial value
     fparm->SetInitValue(conv::kPhypMixing23, TMath::Min(input_mixing23,0.98));
     fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
     //spec func is special. it is nominally 0, but this is the fit boundary (so 'init' is 0.5)
     fparm->Update(conv::kNuisXSecSpecFunc, 0);

     // create fit inputs
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

     // get the specified fitter and perform the fit
     fparm->ResetCurrToInitValues();
     Superk1RmuErecoFitter fitter(fit_inputs,fit_option);
     fitter.DoTheFit();

     // get fitter status
     int migradfitstatus = fitter.FitMigradStatusCode();

     if(migradfitstatus != 4) { //converged
       double fitstatus    = (double)fitter.FitStatusCode();

       // get best-fit params
       double bestfit_mixing23 = fparm->BestFitValue(conv::kPhypMixing23);
       double bestfit_dm23sq   = fparm->BestFitValue(conv::kPhypDm23sq);
       double bestfit_chisq    = fit_inputs->BestFitChisquare();

       // Now do a 0+N parameter fit, to get the chisq true
       fparm->Update(conv::kPhypMixing23, input_mixing23);
       fparm->Update(conv::kPhypDm23sq,   input_dm23sq);
       fparm->FixParameter(conv::kPhypMixing23);
       fparm->FixParameter(conv::kPhypDm23sq);

       // Need to make a new FitInputs & fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm);
       Superk1RmuErecoFitter fitter2(fit_inputs2,fit_option);
       fitter2.DoTheFit();
       double chisq_for_input_params = fit_inputs2->Chisquare();

       //  Print out result
       cout << "- Chisq (Best-fit params:" 
	    << util::AsString(conv::kPhypMixing23) << " = " << bestfit_mixing23 << ", "
	    << util::AsString(conv::kPhypDm23sq)   << " = " << bestfit_dm23sq   << " eV^2"
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
		  bestfit_chisq,
		  fitstatus);

       delete fit_inputs2;
     }//converged

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
