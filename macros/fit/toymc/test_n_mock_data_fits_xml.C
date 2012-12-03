//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are read from the input XML file.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_n_mock_data_fits_xml
(
  const char * tag = "test",                                                /* */
  int          run = 0,                                                /* */
  const char * xmlfile = "/home/dealtry/spare_trunk/trunk/inputs/fakedata.xml",                                            /* input XML file with toy MC experiments */
  int          first_experiment = 1,                                   /* */
  int          last_experiment = 1,                                    /* */
  const char * config_file         = "analysis_configuration.xml", /* analysis configuration file */
  const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisanceselection",    /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option          = "minuit",                     /* fitter option (see run_fitter.C) */
  const char * output_option       = "Observation+BestFitSpectra+ChisqSurface+Canvases",  /* what to write-out once the fit is done */
  bool chisqwithsyst = true
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
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

  // Store oscillation fit results
  TFile f(Form("%s-run_%d-expts_%d_%d.root",
               tag, run, first_experiment, last_experiment), "recreate");
  TNtupleD * nfitnt = new TNtupleD("nfitnt","", 
     "iexp:mixing12:mixing13:mixing23:deltacp:dm12sq:dm23sq:chisq:fitstatus");

  // Loop over toy-MC experiments
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // Re-seed fitter (reset current values to specified init values)
     fparm->ResetCurrToInitValues();

     // Get superk spectrum observation (toy experiment)
     SuperkObservation * skobs = new SuperkObservation;
     skobs->Read(xmlfile,iexperiment);

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
     FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

     // Do the fit
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();

     // Write-out results
     fitter->Analyze(output_option);
     fitter->Save(Form("mockfit_%d.root",iexperiment));

     // Get best-fit chisquare
     double bestfit_chisq = fit_inputs->BestFitChisquare();

     // Get best-fit params
     double bestfit_mixing12 = -1;
     double bestfit_mixing13 = -1;
     double bestfit_mixing23 = -1;
     double bestfit_deltacp  = -1;
     double bestfit_dm12sq   = -1;
     double bestfit_dm23sq   = -1;
     if(hyp == conv::kHypOsc3F) {
        bestfit_mixing12 = fparm->BestFitValue(t2kana::conv::kPhypMixing12);
        bestfit_mixing13 = fparm->BestFitValue(t2kana::conv::kPhypMixing13);
        bestfit_mixing23 = fparm->BestFitValue(t2kana::conv::kPhypMixing23);
        bestfit_deltacp  = fparm->BestFitValue(t2kana::conv::kPhypDeltaCP);
        bestfit_dm12sq   = fparm->BestFitValue(t2kana::conv::kPhypDm12sq);
        bestfit_dm23sq   = fparm->BestFitValue(t2kana::conv::kPhypDm23sq);
     } 
     else
     if(hyp == conv::kHypNuMuDisapp2F) {
        bestfit_mixing23 = fparm->BestFitValue(t2kana::conv::kPhypMixing23);
        bestfit_dm23sq   = fparm->BestFitValue(t2kana::conv::kPhypDm23sq);
     }

     // Print-out best-fit params
     fparm->ShowOnlyBestFitValuesInNextPrint();
     cout << "- Best-fit oscillation params (below) yield chisq = " << bestfit_chisq 
          << *fparm;

     // Store best-fit values
     nfitnt->Fill(
           iexperiment, 
           bestfit_mixing12, bestfit_mixing13, bestfit_mixing23,
           bestfit_deltacp, bestfit_dm12sq, bestfit_dm23sq, 
           bestfit_chisq, fitstatus);

     if(chisqwithsyst) {
       //now get the with-syst chisq surface
       double mixing23_min  = conv::kMinChisqSurfaceSinsq2Theta23;
       double mixing23_max  = conv::kMaxChisqSurfaceSinsq2Theta23;
       double mixing23_step = fparm->StepSize(conv::kPhypMixing23);
       double dm23sq_min    = conv::kMinChisqSurfaceDm23sq;
       double dm23sq_max    = conv::kMaxChisqSurfaceDm23sq;
       double dm23sq_step   = fparm->StepSize(conv::kPhypDm23sq);

       // make sure that min and max are at the bin centre
       int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
       int dm23sq_nbins   = 1 + (int) ((dm23sq_max -   dm23sq_min )  / dm23sq_step );

       TH2D * chisq_surface_syst = new TH2D(
		"chisq_surface_syst","#chi^{2} surface_syst",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
       chisq_surface_syst->SetDirectory(0);
       chisq_surface_syst->GetXaxis()->SetRangeUser(0.,1.);

       for(int i = 1; i <= chisq_surface_syst->GetNbinsX(); i++) {
	 for(int j = 1; j <= chisq_surface_syst->GetNbinsY(); j++) {
	   fparm->ResetCurrToInitValues();
	   double mixing23 = chisq_surface_syst->GetXaxis()->GetBinCenter(i);
	   double dm23sq   = chisq_surface_syst->GetYaxis()->GetBinCenter(j);
	   fparm->Update( conv::kPhypMixing23, mixing23 );
	   fparm->Update( conv::kPhypDm23sq,   dm23sq   );
	   fparm->SetInitValue( conv::kPhypMixing23, mixing23 );
	   fparm->SetInitValue( conv::kPhypDm23sq,   dm23sq   );
	   fparm->FixParameter(conv::kPhypMixing23);
	   fparm->FixParameter(conv::kPhypDm23sq);
	   cout << *fparm;
	   //for some reason, need to make a new fitter when changing the FitParamList
	   FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm);
	   FitterABC * fitter2 = new Superk1RmuErecoFitter(fit_inputs,fit_option);
	   fitter2->DoTheFit();
	   double chisq = fit_inputs2->Chisquare();
	   chisq_surface_syst->SetBinContent(i,j,chisq);
	   delete fitter2;
	   delete fit_inputs2;
	 }//j
       }//i

       chisq_surface_syst->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
       chisq_surface_syst->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
       chisq_surface_syst->GetYaxis()->SetTitleOffset(1);
       chisq_surface_syst->GetYaxis()->SetTitleSize(0.05);
       chisq_surface_syst->GetXaxis()->SetTitleSize(0.05);

       f.cd();
       chisq_surface_syst->Write(Form("chisqwithsyst_%d",iexperiment));
     }


     delete skobs;
     delete fit_inputs;
     delete fitter;
  }

  //
  // Store results
  //
  nfitnt->Write("nfitnt");
  f.Close();
}
