//
// ------------------------------------------------------------------------------------------------
// Build chisq surface (incl. the effect of systematics) for the numu disappearance analysis.
// When all systematics are included this is a CPU-intensive task and it can be pallalelized,
// with each job calculating the chisq value for a single bin (or a row or a column of bins)
// in the osc parameter space. If bin numbers are left at -1, the task is not parallelized.
// The experiment (toy MC or data) is read from the input XML file.
// A range of experiments can also be read and chisq is averaged (for sensitivity studies).
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void build_chisq_surface_disapp_parallelized
(
  const char * tag,                                                         /* */
  const char * xmlfile,                                                     /* input XML file with toy MC experiments or actual T2K data */
  int          first_experiment    = 0,                                     /* */
  int          last_experiment     = 0,                                     /* */
  const char * config_file         = "analysis_configuration.xml",          /* analysis configuration file */
  const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisance",    /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option          = "minuit",                              /* fitter option */
  double       mixing23_min        = 0.9,
  double       mixing23_max        = 1.0,
  double       mixing23_step       = 0.001,
  double       dm32sq_min          = 0.002,
  double       dm32sq_max          = 0.003,
  double       dm32sq_step         = 0.00001,
  int          input_mixing23_bin  = -1,                                    /* if !=-1, will run only for the specified bin */
  int          input_dm32sq_bin    = -1,                                    /* if !=-1, will run only for the specified bin */
  const char * dmsqscheme          = "abs23",                               /* abs23==kNuDmsq12andAbs23 fogli==kNuDmsqFogliLisi */
  const char * hierarchy           = "normal",                              /* normal or inverted */
  const char * mixparamchoice      = "sinsq2theta",                         /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant      = "one",                                 /* one==kNuOctantOne two==kNuOctantTwo */
  bool         apply_2nu_boundary  = true                                   /* if using the 2nu flavour osc prob, apply the physical boundary? */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);

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

  // Calc number of bins & make sure input bin numbers are within range
  int mixing23_nbins = 1 + TMath::CeilNint((mixing23_max - mixing23_min) / mixing23_step);
  int dm32sq_nbins   = 1 + TMath::CeilNint((dm32sq_max -   dm32sq_min )  / dm32sq_step );
  if(input_mixing23_bin != -1) {
    assert(input_mixing23_bin >= 1 && input_mixing23_bin <= mixing23_nbins);
  }
  if(input_dm32sq_bin != -1) {
    assert(input_dm32sq_bin >= 1 && input_dm32sq_bin <= dm32sq_nbins);
  }

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC or data:" << endl;
  cout << "   - Reading experiments from: " << xmlfile2 << endl;
  cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit param list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "- Parameter ranges and step sizes: " << endl;
  cout << "   - mix23  : [" << mixing23_min << ", " << mixing23_max << "], step = " << mixing23_step 
       << ", nbins = " << mixing23_nbins << endl;
  cout << "   - dm32sq : [" << dm32sq_min << ", " << dm32sq_max << "], step = " << dm32sq_step 
       << ", nbins = " << dm32sq_nbins << endl;
  cout << "- Specified bins: " << endl;
  cout << "   - mix23  : " << input_mixing23_bin << endl;
  cout << "   - dm32sq : " << input_dm32sq_bin   << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit param list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Set whether the osc prob can be > 1 in 2 flavour world
  if(!apply_2nu_boundary && hyp == conv::kHypNuMuDisapp2F) {
    skexp->GetSuperkMC()->SetUnphysicalOscProb(true);
    fparm->SetMaxValue( conv::kPhypMixing23, 2);
  }

  // Output file
  ostringstream filename;
  filename << tag << "-";
  if(first_experiment == last_experiment) {
    filename << "expt_" << first_experiment;
  } else {
    filename = "expt_" << first_experiment << "_" << last_experiment;
  }
  if(input_mixing23_bin == -1 && input_dm32sq_bin == -1) {
    filename << "-all_bins";
  } else {
    if(input_mixing23_bin != -1) { filename << "-mix23bin_"  << input_mixing23_bin; }
    if(input_dm32sq_bin   != -1) { filename << "-dm32sqbin_" << input_dm32sq_bin;   }
  }
  filename << ".root";
  TFile f(filename.str().c_str(), "recreate");

  // Tree for storing best-fit params and each bin and for each expt
  int iexpt   = 0; // expt number
  int idm32sq = 0; // dm32sq bin
  int imix23  = 0; // mixing23 bin
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params", "t2kana::FitParamList",  &fparm,     1280000,0);
  fitnt->Branch("iexpt",   &iexpt,   "iexpt/I");
  fitnt->Branch("idm32sq", &idm32sq, "idm32sq/I");
  fitnt->Branch("imix23",  &imix23,  "imix23/I");

  // Define space. Min and max at the centre of first and last bin respectively.
  TH2D * chisq_surface = new TH2D("chisq_surface","VaLOR #chi^{2} surface for numu disappearance analysis",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm32sq_nbins,   dm32sq_min   - 0.5*dm32sq_step,   dm32sq_max   + 0.5*dm32sq_step  );
  chisq_surface->SetDirectory(0);
  chisq_surface->GetXaxis()->SetRangeUser(0.,1.);
  chisq_surface->Reset();

  // Loop over experiments
  int nexperiments = 0;
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {
     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // Re-seed fitter (reset current values to specified init values)
     fparm->ResetCurrToInitValues();

     // Read SuperK spectrum (toy MC or data)
     SuperkObservation * skobs = new SuperkObservation;
     skobs->Read(xmlfile2.Data(), iexperiment);

     // Loop over bins and calc the chisq fixing osc params to the true ones at the centre
     // of the bin and minimize with respect to all systematics (and, potentially, all other osc params).
     bool use_single_mixing23_bin = (input_mixing23_bin != -1);
     bool use_single_dm32sq_bin   = (input_dm32sq_bin   != -1);
     for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
        if(use_single_mixing23_bin) {
          bool bin_match = (input_mixing23_bin == i);
          if(!bin_match) continue;
        }
        for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
	   if(use_single_dm32sq_bin) {
	     bool bin_match = (j == input_dm32sq_bin);
	     if(!bin_match) continue;
	   }
           // init and fix params as needed
           fparm->ResetCurrToInitValues();
           double mixing23 = chisq_surface->GetXaxis()->GetBinCenter(i);
           double dm32sq   = chisq_surface->GetYaxis()->GetBinCenter(j);
           fparm->Update( conv::kPhypMixing23, mixing23 );
           fparm->Update( conv::kPhypDm23sq,   dm32sq   );
           fparm->SetInitValue( conv::kPhypMixing23, mixing23 );
           fparm->SetInitValue( conv::kPhypDm23sq,   dm32sq   );
           fparm->FixParameter(conv::kPhypMixing23);
           fparm->FixParameter(conv::kPhypDm23sq);
           cout << *fparm;
           // do the fit
           FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
           FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
         //fitter->SetHesseErrors(false);
           cout << "Computing chisq for bin centered at mix23 = " << mixing23 << ", Dm^2 = " << dm32sq << " eV^2" << endl;
           cout << "Fit has " << fparm->GetNumberFreeParameters() << " free parameters" << endl;
           fitter->DoTheFit();
           double chisq = fit_inputs->Chisquare();
           delete fitter;
           delete fit_inputs;               
           // store best-fit results in tree
           iexpt   = iexperiment;
           imix23  = i; 
           idm32sq = j;
           fitnt->Fill();
           // update chisq surface
	   chisq_surface->SetBinContent(i,j,chisq+chisq_surface->GetBinContent(i,j));
        }//j
     }//i
     nexperiments++;
     delete skobs;
  }//iexperiment

  // For sensitivity studies, chisq is averaged over N experiments.
  // For the chisq surface of a single dataset, the following has no effect.
  double scale_factor = 1./((double)nexperiments);
  chisq_surface->Sumw2();
  chisq_surface->Scale(scale_factor);

  // Set axes etc
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
    chisq_surface->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
  else
    chisq_surface->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  chisq_surface->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  chisq_surface->GetYaxis()->SetTitleOffset(1);
  chisq_surface->GetYaxis()->SetTitleSize(0.05);
  chisq_surface->GetXaxis()->SetTitleSize(0.05);

  f.cd();
  chisq_surface->Write("chisq_surface");
  fitnt->Write("fitnt");
  f.Close();
}
