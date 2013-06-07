//
// ------------------------------------------------------------------------------------------------
// Build chisq surface (incl. the effect of systematics) for the numu disappearance analysis.
// When all systematics are included this is a CPU-intensive task and it can be pallalelized,
// with each job calculating the chisq value for a single bin in the osc parameter space.
// The experiment (toy MC or data) is read from the input XML file.
// A range of experiments can also be read and chisq is averaged (for sensitivity studies)
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
  const char * fit_option          = "minuit",                              /* fitter option (see run_fitter.C) */
  int          input_mixing23_bin  = -1,                                    /* if !=-1, will run only for the specified bin */
  int          input_dm23sq_bin    = -1,                                    /* if !=-1, will run only for the specified bin */
  double       input_mixing23_min  = 0.9,
  double       input_mixing23_max  = 1.0,
  double       input_mixing23_step = 0.001,
  double       input_dm23sq_min    = 0.002,
  double       input_dm23sq_max    = 0.003,
  double       input_dm23sq_step   = 0.00001
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gOscDmsqScheme        = conv::kNuDmsq12andAbs23;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC or data:" << endl;
  cout << "   - Reading experiments from: " << xmlfile << endl;
  cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit param list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit param list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Output file
  ostringstream filename;
  filename << tag << "-";
  if(first_experiment == last_experiment) {
    filename = "expt_" << first_experiment;
  } else {
    filename = "expt_" << first_experiment << "_" << last_experiment;
  }
  if(input_mixing23_bin == -1 && input_dm23sq_bin == -1) {
    filename << "-all_bins";
  } else {
    if(input_mixing23_bin != -1) { filename << "-mix23bin_"  << input_mixing23_bin; }
    if(input_dm23sq_bin   != -1) { filename << "-dm32sqbin_" << input_dm23sq_bin;   }
  }
  filename << ".root";
  TFile f(Form("%s-expt_%d-bins_%d_%d.root", tag, experiment_id, input_mixing23_bin,input_dm23sq_bin), "recreate");

  // Define 2D osc param space range and binning
  double mixing23_min  = (input_mixing23_min  < 0) ? (conv::kMinChisqSurfaceSinsq2Theta23)  : (input_mixing23_min);
  double mixing23_max  = (input_mixing23_max  < 0) ? (conv::kMaxChisqSurfaceSinsq2Theta23)  : (input_mixing23_max);
  double mixing23_step = (input_mixing23_step < 0) ? (fparm->StepSize(conv::kPhypMixing23)) : (input_mixing23_step);
  double dm23sq_min    = (input_dm23sq_min    < 0) ? (conv::kMinChisqSurfaceDm23sq)         : (input_dm23sq_min);
  double dm23sq_max    = (input_dm23sq_max    < 0) ? (conv::kMaxChisqSurfaceDm23sq)         : (input_dm23sq_max);
  double dm23sq_step   = (input_dm23sq_step   < 0) ? (fparm->StepSize(conv::kPhypDm23sq))   : (input_dm23sq_step);

  // Make sure that min and max are at the bin centre
  int mixing23_nbins = 1 + TMath::CeilNint((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + TMath::CeilNint((dm23sq_max -   dm23sq_min )  / dm23sq_step );

  TH2D * chisq_surface = new TH2D("chisq_surface","VaLOR #chi^{2} surface for numu disappearance analysis",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm23sq_nbins,   dm23sq_min   - 0.5*dm23sq_step,   dm23sq_max   + 0.5*dm23sq_step  );
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
     skobs->Read(xmlfile, iexperiment);

     // Loop over bins and calc the chisq fixing osc params to the true ones at the centre
     // of the bin and minimize with respect to all systematics (and, potentially, all other osc params).
     for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
        for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
           double chisq = 0;
           if(i == input_mixing23_bin && j == input_dm23sq_bin) {
               fparm->ResetCurrToInitValues();
               double mixing23 = chisq_surface->GetXaxis()->GetBinCenter(i);
               double dm23sq   = chisq_surface->GetYaxis()->GetBinCenter(j);
               fparm->Update( conv::kPhypMixing23, mixing23 );
               fparm->Update( conv::kPhypDm23sq,   dm23sq   );
               fparm->SetInitValue( conv::kPhypMixing23, mixing23 );
               fparm->SetInitValue( conv::kPhypDm23sq,   dm23sq   );
               fparm->FixParameter(conv::kPhypMixing23);
               fparm->FixParameter(conv::kPhypDm23sq);
               cout << *fparm;
               FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
               FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
               fitter->SetHesseErrors(false);
               fitter->DoTheFit();
               chisq = fit_inputs->Chisquare();
               delete fitter;
               delete fit_inputs;
           }// = input bin?
	   chisq_surface->SetBinContent(i,j,chisq+chisq_surface->GetBinContent(i,j));
           nexperiments++;
        }//j
     }//i
     delete skobs;
  }//iexperiment

  // For sensitivity studies, chisq is averaged over N experiments.
  // For the chisq surface of a single dataset, the following has no effect.
  double scale_factor = 1./((double)nexperiments);
  chisq_surface->Sumw2();
  chisq_surface->Scale(scale_factor);

  // Set axes etc
  chisq_surface->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  chisq_surface->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  chisq_surface->GetYaxis()->SetTitleOffset(1);
  chisq_surface->GetYaxis()->SetTitleSize(0.05);
  chisq_surface->GetXaxis()->SetTitleSize(0.05);

  f.cd();
  chisq_surface_syst->Write("chisq_surface");
  f.Close();
}
