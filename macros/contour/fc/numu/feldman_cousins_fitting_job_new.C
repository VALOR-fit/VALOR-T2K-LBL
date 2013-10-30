//
// ------------------------------------------------------------------------------------------------
// Run a job with some of the toy expts fits needed to calculate the \Delta\chi^{2} values
// for constructing a Feldman-Cousins contour 
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void feldman_cousins_fitting_job_new
(
  const char * tag                     = "fc",                                /* tag for output file */
  int          run                     = 0,                                   /* run number */
  int          nexperiments            = 1,                                   /* number of toy MC experiment fits */
  int          min_dm23sq_bin          =  25,                                 /* generate toy MC with this range of Dm^2_{23} */
  int          max_dm23sq_bin          =  36,                                 /* -//- */
  int          min_mixing23_bin        =  82,                                 /* generate toy MC with this range of 23 mixing (expressed in terms of either sinsq2theta, sinsqtheta, ...) */
  int          max_mixing23_bin        = 126,                                 /* -//- */
  const char * bin_list                = "",                                  /* if specified, run toy-MC expts only in the bins specified as long as they are in the range given above */
  const char * config_file             = "analysis_configuration.xml",        /* analysis configuration file */
  bool         stat_fluctuations       = true,                                /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics       = true,                                /* generate toy MC with randomized systematics? */
  bool         tweak_non23oscparams    = true,                                /* generate toy MC with randomized dmsq12,theta12,theta13,dcp / only in 3-flavour fmwk */
  const char * fit_param_list          = "nuosc3f_numudisapp_fix4+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option              = "minuit",                            /* fitter option */
  bool         store_toy_expts         = true,                                /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = true,                                /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one"                                /* one==kNuOctantOne two==kNuOctantTwo */
  //bool         apply_2nu_boundary      = true                                 /* if using the 2nu flavour osc prob, apply the physical boundary? */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

  //
  // The oscillation contour will be constructed using a (dm23sq, mixing23) 
  // 2-D histogram with the following binning.
  // Look-up the bin centre values for the input_dm23sq_bin and input_mixing23_bin.
  //

  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.25;
  double mixing23_max   = 0.75;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1.5E-3;
  double dm23sq_max     = 3.5E-3;
  double dm23sq_step    = 5E-5;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int)(((mixing23_max - mixing23_min) / mixing23_step) + 0.01);
  int dm23sq_nbins   = 1 + (int)(((dm23sq_max   - dm23sq_min  ) / dm23sq_step  ) + 0.01);
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

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - Randomize systematics : "      << ((tweak_systematics) ?      "YES " : "NO") << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ?      "YES " : "NO") << endl;
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
       << " : [" << min_mixing23_bin << ", " << max_mixing23_bin << "]"
       << " : [" << osc_param_space->GetXaxis()->GetBinCenter(min_mixing23_bin)
       << ", "   << osc_param_space->GetXaxis()->GetBinCenter(max_mixing23_bin)
       << "]" << endl;
  cout << "   - " << util::AsString(conv::kPhypDm23sq) 
       << " : [" << min_dm23sq_bin   << ", " << max_dm23sq_bin   << "]"
       << " : [" << osc_param_space->GetYaxis()->GetBinCenter(min_dm23sq_bin)
       << ", "   << osc_param_space->GetYaxis()->GetBinCenter(max_dm23sq_bin)
       << "]" << endl;
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
  cout << "Doing a 2D FC scanning dm32sq and mixing23" << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing toy MC experiments:  " << ((store_toy_expts) ?        "YES " : "NO") << endl;
  cout << "     - Storing best-fit spectra: "    << ((store_bestfit_spectrum) ? "YES " : "NO") << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;

  // Fit parameter list
  FitParamList * fparm      = FitParamList::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1     = FitParamList::NamedParameterList(fit_param_list); // for fit using the first seed
  FitParamList * fparm2     = FitParamList::NamedParameterList(fit_param_list); // for fit using the mirror of best fit seed
  FitParamList * fparm_true = FitParamList::NamedParameterList(fit_param_list); // for fit using the true 23 osc params (the 0+N)
  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each toy-MC.
  TObjArray * dat_array = new TObjArray;
  TObjArray * pdf_array = new TObjArray; //only used to initialise branch
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);
  TObjArray * pdf_array1 = new TObjArray;
  pdf_array1->SetOwner(false);
  TObjArray * pdf_array2 = new TObjArray;
  pdf_array2->SetOwner(false);
  TObjArray * pdf_array_true = new TObjArray;
  pdf_array_true->SetOwner(false);

  // Create file to store results
  TFile f(Form("%s-run_%d.root",tag, run), "recreate"); 

  //numbers to store results in
  int iexperiment, status_0N, status_2N, input_mixing23_bin, input_dm23sq_bin;
  double input_mixing23, input_dm23sq;
  double chisq_for_input_params, bestfit_mixing23, bestfit_dm23sq, bestfit_chisq;

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",        "t2kana::FitParamList",  &fparm,          1280000,0);
  fitnt->Branch("best_fit_params_fit1",   "t2kana::FitParamList",  &fparm1,         1280000,0);
  fitnt->Branch("best_fit_params_fit2",   "t2kana::FitParamList",  &fparm2,         1280000,0);
  fitnt->Branch("best_fit_params_true",   "t2kana::FitParamList",  &fparm_true,     1280000,0);
  fitnt->Branch("input_params",           "t2kana::PhysParamList", &iparm,          1280000,0);
  fitnt->Branch("data",                   "TObjArray",             &dat_array,      1280000,0);
  fitnt->Branch("best_fit_pdf",           "TObjArray",             &pdf_array,      1280000,0);
  fitnt->Branch("best_fit_pdf_fit1",      "TObjArray",             &pdf_array1,     1280000,0);
  fitnt->Branch("best_fit_pdf_fit2",      "TObjArray",             &pdf_array2,     1280000,0);
  fitnt->Branch("best_fit_pdf_true",      "TObjArray",             &pdf_array_true, 1280000,0);
  fitnt->Branch("iexp",                   &iexperiment);
  fitnt->Branch("irun",                   &run);
  fitnt->Branch("input_mixing23_bin",     &input_mixing23_bin);
  fitnt->Branch("input_dm23sq_bin",       &input_dm23sq_bin);
  fitnt->Branch("input_mixing23",         &input_mixing23);
  fitnt->Branch("input_dm23sq",           &input_dm23sq);
  fitnt->Branch("chisq_for_input_params", &chisq_for_input_params);
  fitnt->Branch("status_0N",              &status_0N);
  fitnt->Branch("mixing23",               &bestfit_mixing23);
  fitnt->Branch("dm23sq",                 &bestfit_dm23sq);
  fitnt->Branch("chisq",                  &bestfit_chisq);
  fitnt->Branch("status_2N",              &status_2N);
  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  // Loop over toy-MC experiments
  for(iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     //     
     // Random selection of dm23sq and mixing23 bin
     //     
     
     input_mixing23_bin = -1;
     input_dm23sq_bin   = -1;
     
     if(!use_bin_list) {
       // If there was no input bin list, select a bin within the input range
       input_mixing23_bin = min_mixing23_bin + Rnd::Instance()->RndOther().Integer(1 + max_mixing23_bin - min_mixing23_bin);
       input_dm23sq_bin   = min_dm23sq_bin   + Rnd::Instance()->RndOther().Integer(1 + max_dm23sq_bin   - min_dm23sq_bin  );
     } 
     else {
       // If a bin list was input, select a bin from that list and accept it if it is within the input range
       const int nb = mixing23_bin_list.size();
       while(1) {
	 int i = Rnd::Instance()->RndOther().Integer(nb); // random in [0, nb-1] range
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
       if(tweak_non23oscparams) {
         input_mixing12 = util::RandomizeMixing12(); 
         input_mixing13 = util::RandomizeMixing13(); 
         input_deltacp  = util::RandomizeDeltaCP();
         input_dm12sq   = util::RandomizeDm12sq(); 
       }
       else {
	 if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
	   input_mixing12 = conv::kGlobFitSinsq2Theta12;
	   input_mixing13 = conv::kGlobFitSinsq2Theta13;
	 }
	 else {
	   input_mixing12 = conv::kGlobFitSinsqTheta12;
	   input_mixing13 = conv::kGlobFitSinsqTheta13;
	 }
         input_deltacp  = conv::kGlobFitDeltaCP;
         input_dm12sq   = conv::kGlobFitDm12sq;
       }
       cout << "- Toy-MC > True (randomized) non-23-sector params: " 
	    << util::AsString(conv::kPhypMixing12) << " = " << input_mixing12 << ", "
	    << util::AsString(conv::kPhypDm12sq)   << " = " << input_dm12sq   << " eV^2, "
	    << util::AsString(conv::kPhypMixing12) << " = " << input_mixing13 << ", "
	    << util::AsString(conv::kPhypDeltaCP)  << " = " << input_deltacp  << " rad "
	    << endl;
     }

     //
     // Set the generated oscillation parameter values (used to generate toy-MC expt)
     //
     fparm1->Update(conv::kPhypMixing23, input_mixing23);
     fparm1->Update(conv::kPhypDm23sq,   input_dm23sq);
     if(hyp == conv::kHypOsc3F) {
       fparm1->Update(conv::kPhypMixing12, input_mixing12);
       fparm1->Update(conv::kPhypMixing13, input_mixing13);
       fparm1->Update(conv::kPhypDeltaCP,  input_deltacp);
       fparm1->Update(conv::kPhypDm12sq,   input_dm12sq);
     }

     // Reset all output spectra
     dat_array     ->Clear();
     pdf_array     ->Clear();
     pdf_array1    ->Clear();
     pdf_array2    ->Clear();
     pdf_array_true->Clear();

     // Reset input true params and add true oscillation params
     iparm->Clear();
     PhysParamList true_oscillation_params = fparm1->AsPhysParamList("current");
     iparm->Append(true_oscillation_params);

     //
     // get toy MC with/without systematics
     //
     SuperkObservation * skobs = new SuperkObservation;
     bool applied_negevents_condition = true;
     while(applied_negevents_condition) {
       if(tweak_systematics) {
	 // tweak all sources of systematic error
	 cout << "- Randomizing all systematic params in pdf..." << endl;
	 skexp->TweakAllSystematics();     
	 skexp->PrintSystematicVariations();
	 // keep track of all nuisance params set so as to add in output tree
	 PhysParamList true_nuisance_params = skexp->Tweaks();
	 iparm->Append(true_nuisance_params);
	 // generate mock data using the current PDF (incl. systematic variations) 
	 // + statistical fluctuations if requested (Poisson fluctuation in each bin)
	 cout << "- Generating toy-MC..." << endl;
	 //note that fparm does not contain the tweaks, so use only the PhysHyp here
	 applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), stat_fluctuations);
	 // reset all systematics so that the PDF used in the fitter is the nominal one
	 cout << "- Reseting all systematic params in pdf..." << endl;
	 skexp->ResetSystematicVariations();
	 //skexp->PrintSystematicVariations();
       }
       else {
	 //make sure systematics are reset
	 skexp->ResetSystematicVariations();
	 // generate mock data using the current PDF (without systematic variations) 
	 // with statistical fluctuations if requested (Poisson fluctuation in each bin)
	 cout << "- Generating toy-MC..." << endl;
	 applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), stat_fluctuations); 
       }
     }//applied_negevents_condition

     //
     // Get the 2+N fit results
     //

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

     // do the fit using the true osc param as a seed
     fparm1->ResetCurrToInitValues(); // reset osc params
     //help the fitter out by setting a new initial value
     fparm1->SetInitValue(conv::kPhypMixing23, TMath::Min(input_mixing23,0.98));
     fparm1->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
     fitter->DoTheFit();
     fparm1->ResetCurrToInitValues();
     int fitstatus1 = fitter->FitStatusCode("migrad_min");
     bestfit_chisq  = fit_inputs->BestFitChisquare();
     double bestfit_chisq2 = 0;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }

     //Do the fit again at the mirror point
     int fitstatus2 = -1;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       //get the correct intial values (best fit points from first fit)
       fparm2->CopyBestFitToInitValues(*fparm1);
       //set the mixing23 seed to be the mirror point of the initial best fit
       double sinsqtheta23_bestfit = fparm1->BestFitValue(conv::kPhypMixing23);
       double sinsqtheta13_bestfit = fparm1->BestFitValue(conv::kPhypMixing13);
       double centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       double reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
	 reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
	 reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm2->SetInitValue(conv::kPhypMixing23, reseed);
       // Get fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
       FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
       fparm2->ResetCurrToInitValues(); // reset osc params
       fitter2->DoTheFit();
       bestfit_chisq2 = fit_inputs2->BestFitChisquare();
       fitstatus2 = fitter2->FitStatusCode("migrad_min");
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
	 skexp->Write(pdf_array2, "pdf");
       }
     }

     // Print-out best-fit params, and set the "best of the two fits" branches correctly
     cout << "- Best-fit oscillation params (below) yield chisq = ";
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if(bestfit_chisq2 < bestfit_chisq) {
	 cout << bestfit_chisq2 << " from the second seed fit" << endl;
	 bestfit_chisq    = bestfit_chisq2;
	 bestfit_mixing23 = fparm2->BestFitValue(conv::kPhypMixing23);
	 bestfit_dm23sq   = fparm2->BestFitValue(conv::kPhypDm23sq);
	 status_2N        = fitstatus2;
	 fparm2->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2;
	 br_fparm->SetAddress(&fparm2);
	 br_bf_array->SetAddress(&pdf_array2);
       }
       else {
	 cout << bestfit_chisq << " from the first seed fit" << endl;
	 bestfit_mixing23 = fparm1->BestFitValue(conv::kPhypMixing23);
	 bestfit_dm23sq   = fparm1->BestFitValue(conv::kPhypDm23sq);
	 status_2N        = fitstatus1;
	 fparm1->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1;
	 br_fparm->SetAddress(&fparm1);
	 br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << bestfit_chisq << endl;
       bestfit_mixing23 = fparm1->BestFitValue(conv::kPhypMixing23);
       bestfit_dm23sq   = fparm1->BestFitValue(conv::kPhypDm23sq);
       status_2N        = fitstatus1;
       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1;
       br_fparm->SetAddress(&fparm1);
       br_bf_array->SetAddress(&pdf_array1);
     }
     
     //
     // Get the 0+N fit results
     //
     fparm_true->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm_true->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
     fparm_true->FixParameter(conv::kPhypMixing23);
     fparm_true->FixParameter(conv::kPhypDm23sq);

     // Need to make a new FitInputs & fitter
     FitInputs * fit_inputs_true = new FitInputs(skobs,skexp,fparm_true);
     FitterABC * fitter_true = util::GetFitter(fit_inputs_true,fit_option);
     fparm_true->ResetCurrToInitValues(); // reset osc params
     fitter_true->DoTheFit();
     status_0N = fitter_true->FitStatusCode("migrad_min");
     chisq_for_input_params = fit_inputs_true->BestFitChisquare();
     // Save spectra (best-fit prediction)
     if(store_bestfit_spectrum) {
       skexp->Write(pdf_array_true, "pdf");
     }
     cout << "- True oscillation params (below) yield chisq = "
	  << chisq_for_input_params << " from the second seed fit" << endl;
     fparm_true->ShowOnlyBestFitValuesInNextPrint();
     cout << *fparm_true;

     //  Print out result
     cout << "- Chisq (Best-fit params:" 
	  << util::AsString(conv::kPhypMixing23) << " = " << bestfit_mixing23 << ", "
	  << util::AsString(conv::kPhypDm23sq)   << " = " << bestfit_dm23sq   << " eV^2"
	  << ") = " << bestfit_chisq << endl;
     cout << "- Chisq (True) = " << chisq_for_input_params << endl;
     cout << "- DeltaChisq = " << chisq_for_input_params - bestfit_chisq << endl;

     // Add all results in the tree
     fitnt->Fill();

     // Clean-up
     delete skobs;
     delete fit_inputs;
     delete fitter;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       delete fit_inputs2;
       delete fitter2;
     }
     delete fit_inputs_true;
     delete fitter_true;

  }//iexperiment

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();
}
