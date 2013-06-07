//
// ------------------------------------------------------------------------------------------------
// Build profiles of the chisq surface (incl. systematics) for the numu disappearance analysis.
// The experiment (toy MC or data) is read from the input XML file.
// A range of experiments can also be read and chisq is averaged (for sensitivity studies)
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void build_chisq_profiles_disapp
(
  const char * tag,                                                         /* */
  const char * xmlfile,                                                     /* input XML file with toy MC experiments or actual T2K data */
  int          first_experiment    = 0,                                     /* */
  int          last_experiment     = 0,                                     /* */
  const char * config_file         = "analysis_configuration.xml",          /* analysis configuration file */
  const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisance",    /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option          = "minuit",                              /* fitter option */
  double       input_mixing23_min  = 0.9,
  double       input_mixing23_max  = 1.0,
  double       input_mixing23_step = 0.001,
  double       input_dm32sq_min    = 0.002,
  double       input_dm32sq_max    = 0.003,
  double       input_dm32sq_step   = 0.00001
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
  cout << "- Job: " << tag << endl;
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
    filename << "expt_" << first_experiment;
  } else {
    filename << "expt_" << first_experiment << "_" << last_experiment;
  }
  filename << ".root";
  TFile f(filename.str().c_str(), "recreate");

  // Define osc param space range and binning
  double mixing23_min  = (input_mixing23_min  < 0) ? (conv::kMinChisqSurfaceSinsq2Theta23)  : (input_mixing23_min);
  double mixing23_max  = (input_mixing23_max  < 0) ? (conv::kMaxChisqSurfaceSinsq2Theta23)  : (input_mixing23_max);
  double mixing23_step = (input_mixing23_step < 0) ? (fparm->StepSize(conv::kPhypMixing23)) : (input_mixing23_step);
  double dm32sq_min    = (input_dm32sq_min    < 0) ? (conv::kMinChisqSurfaceDm23sq)         : (input_dm32sq_min);
  double dm32sq_max    = (input_dm32sq_max    < 0) ? (conv::kMaxChisqSurfaceDm23sq)         : (input_dm32sq_max);
  double dm32sq_step   = (input_dm32sq_step   < 0) ? (fparm->StepSize(conv::kPhypDm23sq))   : (input_dm32sq_step);

  // Make sure that min and max are at the bin centre
  int mixing23_nbins = 1 + TMath::CeilNint((mixing23_max - mixing23_min) / mixing23_step);
  int dm32sq_nbins   = 1 + TMath::CeilNint((dm32sq_max -   dm32sq_min )  / dm32sq_step );

  TH1D * chisq_mixing23 = new TH1D(
          "chisq_mixing23","VaLOR numu disappearance analysis #chi^{2} = f(mixing23)",
           mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step);
  chisq_mixing23->SetDirectory(0);
  chisq_mixing23->GetXaxis()->SetRangeUser(0.,1.);
  chisq_mixing23->Reset();

  TH1D * chisq_dm32sq = new TH1D(
          "chisq_dm32sq","VaLOR numu disappearance analysis #chi^{2} = f(Dm^2_{atm})",
           dm32sq_nbins, dm32sq_min - 0.5*dm32sq_step, dm32sq_max + 0.5*dm32sq_step);
  chisq_dm32sq->SetDirectory(0);
  chisq_dm32sq->Reset();

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
     for(int i = 1; i <= chisq_mixing23->GetNbinsX(); i++) {
        fparm->ResetCurrToInitValues();
        double mixing23 = chisq_mixing23->GetXaxis()->GetBinCenter(i);
        fparm->Update       ( conv::kPhypMixing23, mixing23 );
        fparm->SetInitValue ( conv::kPhypMixing23, mixing23 );
        fparm->SetInitValue ( conv::kPhypDm23sq,   conv::kInitFitDm23sq );
        fparm->FreeParameter( conv::kPhypDm23sq   );
        fparm->FixParameter ( conv::kPhypMixing23 );
        cout << *fparm;
        FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
        FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
      //fitter->SetHesseErrors(false);
        cout << "Computing chisq for bin centered at mix23 = " << mixing23 << endl;
        cout << "Fit has " << fparm->GetNumberFreeParameters() << " free parameters" << endl;
        fitter->DoTheFit();
        double chisq = fit_inputs->Chisquare();
        delete fitter;
        delete fit_inputs;
	chisq_mixing23->SetBinContent(i,chisq+chisq_mixing23->GetBinContent(i));
     }//i
     fparm->FreeParameter ( conv::kPhypMixing23 );

     // Loop over bins and calc the chisq fixing osc params to the true ones at the centre
     // of the bin and minimize with respect to all systematics (and, potentially, all other osc params).
     for(int i = 1; i <= chisq_dm32sq->GetNbinsX(); i++) {
        fparm->ResetCurrToInitValues();
        double dm32sq   = chisq_dm32sq->GetXaxis()->GetBinCenter(i);
        fparm->Update       ( conv::kPhypDm23sq,   dm32sq   );
        fparm->SetInitValue ( conv::kPhypDm23sq,   dm32sq   );
        fparm->SetInitValue ( conv::kPhypMixing23, conv::kInitFitSinsq2Theta23 );
        fparm->FreeParameter(conv::kPhypMixing23);
        fparm->FixParameter (conv::kPhypDm23sq  );
        cout << *fparm;
        FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
        FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
      //fitter->SetHesseErrors(false);
        cout << "Computing chisq for bin centered at Dm^2 = " << dm32sq << " eV^2" << endl;
        cout << "Fit has " << fparm->GetNumberFreeParameters() << " free parameters" << endl;
        fitter->DoTheFit();
        double chisq = fit_inputs->Chisquare();
        delete fitter;
        delete fit_inputs;
	chisq_dm32sq->SetBinContent(i,chisq+chisq_dm32sq->GetBinContent(i));
     }//i

     nexperiments++;
     delete skobs;

  }//iexperiment

  // For sensitivity studies, chisq is averaged over N experiments.
  // For the chisq surface of a single dataset, the following has no effect.
  double scale_factor = 1./((double)nexperiments);
  chisq_mixing23->Sumw2();
  chisq_dm32sq->Sumw2();
  chisq_mixing23->Scale(scale_factor);
  chisq_dm32sq->Scale(scale_factor);

  // Set axes etc
  chisq_mixing23->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  chisq_mixing23->GetXaxis()->SetTitleSize(0.05);
  chisq_dm32sq->GetXaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  chisq_dm32sq->GetXaxis()->SetTitleSize(0.05);

  // save
  f.cd();
  chisq_mixing23->Write();
  chisq_dm32sq  ->Write();
  f.Close();
}
