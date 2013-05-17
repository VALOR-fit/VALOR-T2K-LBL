//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are read from the input XML file.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_n_mock_data_fits_xml_parallelized
(
  const char * tag,                                                /* */
  int          run,                                                /* */
  const char * xmlfile,                                            /* input XML file with toy MC experiments */
  int          first_experiment    = 0,                            /* */
  int          last_experiment     = 0,                            /* */
  const char * config_file         = "analysis_configuration.xml", /* analysis configuration file */
  const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisanceselection",    /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option          = "minuit",                     /* fitter option (see run_fitter.C) */
  const char * output_option       = "Observation+BestFitSpectra+ChisqSurface+Canvases",  /* what to write-out once the fit is done */
  int          input_mixing23_bin  = -1,
  int          input_dm23sq_bin    = -1,
  bool         chisqwithsyst       = true,
  double       chisqsyst_mixing23_min = 0.9,
  double       chisqsyst_mixing23_max = 1.0,
  double       chisqsyst_dm23sq_min   = 0.002,
  double       chisqsyst_dm23sq_max   = 0.003,
  bool         chisqwithsyst_mixing23 = true,
  double       chisqwithsyst_mixing23_step = 0.001,
  bool         chisqwithsyst_dm23sq = true,
  double       chisqwithsyst_dm23sq_step = 0.00001
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
     skexp->PrintSystematicVariations();
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();

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

     // Write-out results
     fitter->Analyze(output_option);
     fitter->Save(Form("mockfit_%d.root",iexperiment));

     if(chisqwithsyst) {
       //now get the with-syst chisq surface
       double mixing23_min  = (chisqsyst_mixing23_min < 0) ? (conv::kMinChisqSurfaceSinsq2Theta23) : (chisqsyst_mixing23_min);
       double mixing23_max  = (chisqsyst_mixing23_max < 0) ? (conv::kMaxChisqSurfaceSinsq2Theta23) : (chisqsyst_mixing23_max);
       double mixing23_step = (chisqwithsyst_mixing23_step < 0) ? (fparm->StepSize(conv::kPhypMixing23)) : (chisqwithsyst_mixing23_step);
       double dm23sq_min    = (chisqsyst_dm23sq_min   < 0) ? (conv::kMinChisqSurfaceDm23sq)        : (chisqsyst_dm23sq_min);
       double dm23sq_max    = (chisqsyst_dm23sq_max   < 0) ? (conv::kMaxChisqSurfaceDm23sq)        : (chisqsyst_dm23sq_max);
       double dm23sq_step   = (chisqwithsyst_dm23sq_step < 0) ? (fparm->StepSize(conv::kPhypDm23sq)) : (chisqwithsyst_dm23sq_step);

       // make sure that min and max are at the bin centre
       int mixing23_nbins = 1 + TMath::CeilNint((mixing23_max - mixing23_min) / mixing23_step);
       int dm23sq_nbins   = 1 + TMath::CeilNint((dm23sq_max -   dm23sq_min )  / dm23sq_step );

       TH2D * chisq_surface_syst = new TH2D(
		"chisq_surface_syst","#chi^{2} surface_syst",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
       chisq_surface_syst->SetDirectory(0);
       chisq_surface_syst->GetXaxis()->SetRangeUser(0.,1.);

       for(int i = 1; i <= chisq_surface_syst->GetNbinsX(); i++) {
	 for(int j = 1; j <= chisq_surface_syst->GetNbinsY(); j++) {
           double chisq = 0;
           if(i == input_mixing23_bin && j == input_dm23sq_bin) {
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
	       Superk1RmuErecoFitter * fitter2 = new Superk1RmuErecoFitter(fit_inputs,fit_option);
	       fitter2->SetHesseErrors(false);
	       fitter2->DoTheFit();
	       chisq = fit_inputs2->Chisquare();
	       delete fitter2;
	       delete fit_inputs2;
           }// = input bin?
	   chisq_surface_syst->SetBinContent(i,j,chisq);
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

     if(chisqwithsyst_mixing23) {
       //now get the with-syst chisq mixing23
       double mixing23_min  = (chisqsyst_mixing23_min < 0) ? (conv::kMinChisqSurfaceSinsq2Theta23) : (chisqsyst_mixing23_min);
       double mixing23_max  = (chisqsyst_mixing23_max < 0) ? (conv::kMaxChisqSurfaceSinsq2Theta23) : (chisqsyst_mixing23_max);
       double mixing23_step = (chisqwithsyst_mixing23_step < 0) ? (fparm->StepSize(conv::kPhypMixing23)) : (chisqwithsyst_mixing23_step);
       // make sure that min and max are at the bin centre
       int mixing23_nbins = 1 + TMath::CeilNint((mixing23_max - mixing23_min) / mixing23_step);

       TH1D * chisq_mixing23_syst = new TH1D(
					     "chisq_mixing23_syst","#chi^{2} mixing23_syst",
					     mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step);
       chisq_mixing23_syst->SetDirectory(0);
       chisq_mixing23_syst->GetXaxis()->SetRangeUser(0.,1.);

       fparm->SetInitValue (conv::kPhypDm23sq,   conv::kInitFitDm23sq);
       fparm->FreeParameter(conv::kPhypDm23sq);

       for(int i = 1; i <= chisq_mixing23_syst->GetNbinsX(); i++) {
	 fparm->ResetCurrToInitValues();
	 double mixing23 = chisq_mixing23_syst->GetXaxis()->GetBinCenter(i);
	 fparm->Update( conv::kPhypMixing23, mixing23 );
	 fparm->SetInitValue( conv::kPhypMixing23, mixing23 );
	 fparm->FixParameter(conv::kPhypMixing23);
	 cout << *fparm;
	 //for some reason, need to make a new fitter when changing the FitParamList
	 FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm);
	 Superk1RmuErecoFitter * fitter2 = new Superk1RmuErecoFitter(fit_inputs,fit_option);
	 fitter2->SetHesseErrors(false);
	 fitter2->DoTheFit();
	 double chisq = fit_inputs2->Chisquare();
	 chisq_mixing23_syst->SetBinContent(i,chisq);
	 delete fitter2;
	 delete fit_inputs2;
       }//i
       chisq_mixing23_syst->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
       chisq_mixing23_syst->GetXaxis()->SetTitleSize(0.05);
       f.cd();
       chisq_mixing23_syst->Write(Form("chisqwithsyst_mixing23_%d",iexperiment));
     }//chisqwithsyst_mixing23

     if(chisqwithsyst_dm23sq) {
       //now get the with-syst chisq dm23sq
       double dm23sq_min  = (chisqsyst_dm23sq_min < 0) ? (conv::kMinChisqSurfaceDm23sq) : (chisqsyst_dm23sq_min);
       double dm23sq_max  = (chisqsyst_dm23sq_max < 0) ? (conv::kMaxChisqSurfaceDm23sq) : (chisqsyst_dm23sq_max);
       double dm23sq_step = (chisqwithsyst_dm23sq_step < 0) ? (fparm->StepSize(conv::kPhypDm23sq)) : (chisqwithsyst_dm23sq_step);
       // make sure that min and max are at the bin centre
       int dm23sq_nbins = 1 + TMath::CeilNint((dm23sq_max - dm23sq_min) / dm23sq_step);

       TH1D * chisq_dm23sq_syst = new TH1D(
					     "chisq_dm23sq_syst","#chi^{2} dm23sq_syst",
					     dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
       chisq_dm23sq_syst->SetDirectory(0);
       chisq_dm23sq_syst->GetXaxis()->SetRangeUser(0.,1.);

       fparm->SetInitValue (conv::kPhypMixing23,   conv::kInitFitSinsq2Theta23);
       fparm->FreeParameter(conv::kPhypMixing23);

       for(int i = 1; i <= chisq_dm23sq_syst->GetNbinsX(); i++) {
	 fparm->ResetCurrToInitValues();
	 double dm23sq = chisq_dm23sq_syst->GetXaxis()->GetBinCenter(i);
	 fparm->Update( conv::kPhypDm23sq, dm23sq );
	 fparm->SetInitValue( conv::kPhypDm23sq, dm23sq );
	 fparm->FixParameter(conv::kPhypDm23sq);
	 cout << *fparm;
	 //for some reason, need to make a new fitter when changing the FitParamList
	 FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm);
	 Superk1RmuErecoFitter * fitter2 = new Superk1RmuErecoFitter(fit_inputs,fit_option);
	 fitter2->SetHesseErrors(false);
	 fitter2->DoTheFit();
	 double chisq = fit_inputs2->Chisquare();
	 chisq_dm23sq_syst->SetBinContent(i,chisq);
	 delete fitter2;
	 delete fit_inputs2;
       }//i
       chisq_dm23sq_syst->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
       chisq_dm23sq_syst->GetXaxis()->SetTitleSize(0.05);
       f.cd();
       chisq_dm23sq_syst->Write(Form("chisqwithsyst_dm23sq_%d",iexperiment));
     }//chisqwithsyst_dm23sq


     delete skobs;
     delete fit_inputs;
     delete fitter;
  }//iexperiment

  //
  // Store results
  //
  f.cd();
  nfitnt->Write("nfitnt");
  f.Close();
}
