//
// ------------------------------------------------------------------------------------------------
// Build profiles of the chisq surface (incl. systematics) for the numu disappearance analysis.
// The experiment (toy MC or data) is read from the input XML file.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void build_chisq_profiles_disapp
(
 const char * tag,                                                         /* */
 const char * xmlfile,                                                     /* input XML file with toy MC experiment or actual T2K data */
 int          iexperiment,                                                 /* */
 const char * config_file         = "analysis_configuration.xml",          /* analysis configuration file */
 const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisance",    /* fit parameter list for the given physics hypothesis (see UtlFitParamListFactory.h) */
 const char * fit_option          = "minuit",                              /* fitter option */
 double       mixing23_min        = 0.3,
 double       mixing23_max        = 0.7,
 double       mixing23_step       = 0.001,
 double       dm32sq_min          = 0.002,
 double       dm32sq_max          = 0.003,
 double       dm32sq_step         = 0.00001,
 bool         do_mixing23         = true,
 bool         do_dm32sq           = false,
 int          input_mixing23_bin  = -1,                                    /* if !=-1, will run only for the specified bin */
 int          input_dm32sq_bin    = -1,                                    /* if !=-1, will run only for the specified bin */
 const char * dmsqscheme          = "abs23",                               /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy           = "normal",                              /* normal or inverted */
 const char * mixparamchoice      = "sinsqtheta",                          /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant      = "one"                                  /* one==kNuOctantOne two==kNuOctantTwo */
 //bool         apply_2nu_boundary  = true                                 /* if using the 2nu flavour osc prob, apply the physical boundary? */
 )
{
  using namespace std;
  using namespace valor;

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

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC or data:" << endl;
  cout << "   - Reading experiment from: " << xmlfile2 << endl;
  cout << "   - Using experiment iexperiment" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit param list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit param list
  FitParamList * fparm = util::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Output file
  ostringstream filename;
  filename << tag << "-";
  filename << "expt_" << iexperiment;
  filename << ".root";
  TFile f(filename.str().c_str(), "recreate");

  // Define osc param space, making sure that min and max are at the bin centre
  int mixing23_nbins = 1 + (int)(((mixing23_max - mixing23_min) / mixing23_step) + 0.01);
  int dm32sq_nbins   = 1 + (int)(((dm32sq_max -   dm32sq_min )  / dm32sq_step ) + 0.01);
  if(input_mixing23_bin != -1) {
    assert(input_mixing23_bin >= 1 && input_mixing23_bin <= mixing23_nbins);
  }
  if(input_dm32sq_bin != -1) {
    assert(input_dm32sq_bin >= 1 && input_dm32sq_bin <= dm32sq_nbins);
  }

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

  cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

  // Re-seed fitter (reset current values to specified init values)
  fparm->ResetCurrToInitValues();

  // Read SuperK spectrum (toy MC or data)
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(xmlfile2.Data(), iexperiment);

  // Loop over bins and calc the chisq fixing osc params to the true ones at the centre
  // of the bin and minimize with respect to all systematics (and, potentially, all other osc params).
  if(do_mixing23) {
    for(int i = 1; i <= chisq_mixing23->GetNbinsX(); i++) {
      if(input_mixing23_bin != -1 && input_mixing23_bin != i)
	continue;
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
      double chisq = fit_inputs->BestFitChisquare();
      delete fitter;
      delete fit_inputs;
      chisq_mixing23->SetBinContent(i,chisq+chisq_mixing23->GetBinContent(i));
    }//i
    fparm->FreeParameter ( conv::kPhypMixing23 );
  }//do_mixing23

  // Loop over bins and calc the chisq fixing osc params to the true ones at the centre
  // of the bin and minimize with respect to all systematics (and, potentially, all other osc params).
  if(do_dm32sq) {
    for(int i = 1; i <= chisq_dm32sq->GetNbinsX(); i++) {
      if(input_dm32sq_bin != -1 && input_dm32sq_bin != i)
	continue;
      fparm->ResetCurrToInitValues();
      double dm32sq   = chisq_dm32sq->GetXaxis()->GetBinCenter(i);
      fparm->Update       ( conv::kPhypDm23sq,   dm32sq   );
      fparm->SetInitValue ( conv::kPhypDm23sq,   dm32sq   );
      if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
	fparm->SetInitValue ( conv::kPhypMixing23, conv::kInitFitSinsq2Theta23 );
      else
	fparm->SetInitValue ( conv::kPhypMixing23, conv::kInitFitSinsqTheta23 );
      fparm->FreeParameter(conv::kPhypMixing23);
      fparm->FixParameter (conv::kPhypDm23sq  );
      cout << *fparm;
      FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
      FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
      //fitter->SetHesseErrors(false);
      cout << "Computing chisq for bin centered at Dm^2 = " << dm32sq << " eV^2" << endl;
      cout << "Fit has " << fparm->GetNumberFreeParameters() << " free parameters" << endl;
      fitter->DoTheFit();
      double chisq = fit_inputs->BestFitChisquare();
      delete fitter;
      delete fit_inputs;
      chisq_dm32sq->SetBinContent(i,chisq+chisq_dm32sq->GetBinContent(i));
    }//i
  }

  delete skobs;

  // Set axes etc
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
    chisq_mixing23->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  else
    chisq_mixing23->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
  chisq_mixing23->GetXaxis()->SetTitleSize(0.05);
  chisq_dm32sq->GetXaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  chisq_dm32sq->GetXaxis()->SetTitleSize(0.05);

  // save
  f.cd();
  chisq_mixing23->Write();
  chisq_dm32sq  ->Write();
  f.Close();
}
