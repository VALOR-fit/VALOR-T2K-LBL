//
// ------------------------------------------------------------------------------------------------
// Run fits to experiments read from an input XML file and store best-fit parameters and spectra.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_fit_xml
(
  const char * tag,                                                    /* */
  int          run,                                                    /* */
  const char * xmlfile,                                                /* input XML file with toy MC experiments */
  int          first_experiment = 0,                                   /* */
  int          last_experiment  = 0,                                   /* */
  const char * config_file      = "analysis_configuration.xml",        /* analysis configuration file */
  const char * fit_param_list   = "nuosc3f_numudisapp_fix4+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option       = "minuit",                            /* fitter option */
  bool         store_toy_expts         = true,                         /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = true,                         /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme          = "abs23",                               /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy           = "normal",                              /* normal or inverted */
  const char * mixparamchoice      = "sinsqtheta",                          /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant      = "one",                                 /* one==kNuOctantOne two==kNuOctantTwo */
  //bool         apply_2nu_boundary  = true                                   /* if using the 2nu flavour osc prob, apply the physical boundary? */
  double       dcp_value           = 0,                                      /* value of deltacp (in units of PI) to assume (when deltacp is not a free parameter) */
  bool         using_xml_file      = true                                   /* if false, skobs reads xmlfile as a root file */
)
{
  using namespace std;
  using namespace valor;

  // Analysis configuration
  TString config   = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
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
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiments from: " << xmlfile2 << endl; 
  cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl; 
  cout << "- Fit" << endl;
  cout << "   - Fit param list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit param list
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = util::NamedParameterList(fit_param_list);
  FitParamList * fparm2 = util::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  //make an assumption for deltacp
  if(fparm->IsFixed(conv::kPhypCPV    )) {
    fparm ->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm1->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm2->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
  }

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each expt.
  TObjArray * dat_array = new TObjArray; //only used to initialise branch
  TObjArray * pdf_array = new TObjArray; //only used to initialise branch
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);
  TObjArray * pdf_array1 = new TObjArray;
  pdf_array1->SetOwner(false);
  TObjArray * pdf_array2 = new TObjArray;
  pdf_array2->SetOwner(false);

  // Create file to store results
  TFile f(Form("%s-run_%d-expts_%d_%d.root", tag, run, first_experiment, last_experiment), "recreate");

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("best_fit_params_fit1", "valor::FitParamList",  &fparm1,     1280000,0);
  fitnt->Branch("best_fit_params_fit2", "valor::FitParamList",  &fparm2,     1280000,0);
  //fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("data",                 "TObjArray",             &dat_array,  1280000,0);
  fitnt->Branch("best_fit_pdf",         "TObjArray",             &pdf_array,  1280000,0);
  fitnt->Branch("best_fit_pdf_fit1",    "TObjArray",             &pdf_array1, 1280000,0);
  fitnt->Branch("best_fit_pdf_fit2",    "TObjArray",             &pdf_array2, 1280000,0);
  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  //Open the ROOT file that has the input PhysParamList (if applicable)
  TFile * f_sysparm = 0;
  if(!using_xml_file) {
    f_sysparm = new TFile(xmlfile2.Data(),"READ");
    if(f_sysparm->IsZombie()) {
      cerr << "ROOT file not found. Exiting..." << endl;
      exit(1);
    }
  }

  // Loop over toy-MC experiments
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     // Reset output spectra
     dat_array ->Clear();
     pdf_array ->Clear();
     pdf_array1->Clear();
     pdf_array2->Clear();

     // Re-seed fitter (reset current values to specified init values)
     fparm1->ResetCurrToInitValues();

     // Get superk spectrum observation (toy experiment)
     SuperkObservation * skobs = new SuperkObservation;
     skobs->Read(xmlfile2.Data(),iexperiment,using_xml_file);
     //Update the fitter initial values (osc & systematic)
     PhysParamList * iparm = 0;
     if(!using_xml_file) {
       f_sysparm->GetObject(Form("iparm_%d",iexperiment),iparm);
       fparm1->FromPhysParamList(*iparm, "initial", false);
     }
     
     cout << *fparm1 << endl;

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

     // Do the fit
     skexp->PrintSystematicVariations();
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();
     double bestfit_chisq = fparm->GetBestFitStatistic();

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }

     //Do the fit again
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
       bestfit_chisq2 = fparm2->GetBestFitStatistic();
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
	 fparm2->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2;
	 br_fparm->SetAddress(&fparm2);
	 br_bf_array->SetAddress(&pdf_array2);
       }
       else {
	 cout << bestfit_chisq << " from the first seed fit" << endl;
	 fparm1->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1;
	 br_fparm->SetAddress(&fparm1);
	 br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << bestfit_chisq << endl;
       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1;
       br_fparm->SetAddress(&fparm1);
       br_bf_array->SetAddress(&pdf_array1);
     }

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

  }//iexperiment
  if(f_sysparm) {
    if(f_sysparm->IsOpen())
      f_sysparm->Close();
    delete f_sysparm;
  }

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();
}
