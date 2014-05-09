//
// ------------------------------------------------------------------------------------------------
// Run fits to real or fake datasets in a general way, selecting if the fit is a single one or 
// a 1D fit (vs one of the oscillation parameters) or a 2D fit (vs two of the oscillation params)
// Initial choices for the tolerance and reactor constraint applied are made 
//
// Lorena Escudero
// ------------------------------------------------------------------------------------------------
//

void run_fit_data_joint
(
  const char * tag                     = "data_fit_test",                       /* tag for output file */
  const char * data_file               = "inputs/data.xml",                 /* data file to be read */
  int          data_id                 = 0,                                 /* data index inside the data file */
  const char * config_file             = "analysis_configuration.xml",      /* analysis configuration file */
  const char * fit_param_list          = "nuosc3f_12fix+nuisance",        /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option              = "minuit",                            /* fitter option */
  bool         store_bestfit_spectrum  = true,                                /* disable this if you don't want to store the best-fit spectra */
  bool         reactor                 = false,                               /* enable this if you want to apply the reactor constraint */
  const char * dmsqscheme              = "fogli",                         /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  double       tolerance               = 0.1,                                 /* value for the fitter tolerance */
  const char * param1                  = "none",                        /* first parameter to be fitting vs - default: none for a single fit */ 
  const char * param2                  = "none",                  /* second parameter to be fitting vs (for 2D) - default: none for a single fit or 1D */ 
  int          min_bin_param1          = 1 ,                                /* minimum bin for param1 */
  int          max_bin_param1          = 1 ,                                /* maximum bin for param1 */
  int          min_bin_param2          = 1 ,                                /* minimum bin for param2 */
  int          max_bin_param2          = 1                                  /* maximum bin for param2 */
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  const char * data_file_fullpath = Form("%s/%s", gSystem->Getenv("T2KOSCANA"), data_file);

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  
  Opt::Instance()->SetDefaultParameterValues();

  // Print settings
  cout << "\n[Data fit]" << endl;
  cout << "- File name: " << tag << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Reading data from: " << data_file << endl;
  cout << "            data id: " << data_id   << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "   - Fitter tolerance: " << tolerance << endl;
  cout << "  Fit vs " << param1 << " and " << param2 << endl; 
  cout << "  Applying reactor constraint: " << ((reactor)  ? "YES" : "NO" ) << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;

  // Fit parameter list
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = util::NamedParameterList(fit_param_list); //first fit (first octant)
  FitParamList * fparm2 = util::NamedParameterList(fit_param_list); //second fit (second octant)
  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Read the real (or fake) data set
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, data_id);

  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each toy-MC.
  TObjArray * dat_array = new TObjArray; //only used to initialise branch
  TObjArray * pdf_array = new TObjArray; //only used to initialise branch
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);
  TObjArray * pdf_array1 = new TObjArray;
  pdf_array1->SetOwner(false);
  TObjArray * pdf_array2 = new TObjArray;
  pdf_array2->SetOwner(false);
  double fittime1, fittime2;
  int fitstatus1, fitstatus2;
  double chisq1, chisq2;

  // Create file to store results
  TFile f(Form("%s.root", tag), "recreate");

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("fittime1",&fittime1,"fittime1/D");
  fitnt->Branch("fittime2",&fittime2,"fittime2/D");
  fitnt->Branch("fitstatus1",&fitstatus1,"fitstatus1/I");
  fitnt->Branch("fitstatus2",&fitstatus2,"fitstatus2/I");
  fitnt->Branch("chisq1",&chisq1,"chisq1/D");
  fitnt->Branch("chisq2",&chisq2,"chisq2/D");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("best_fit_params_fit1", "valor::FitParamList",  &fparm1,     1280000,0);
  fitnt->Branch("best_fit_params_fit2", "valor::FitParamList",  &fparm2,     1280000,0);
  fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("data",                 "TObjArray",             &dat_array,  1280000,0);
  fitnt->Branch("best_fit_pdf",         "TObjArray",             &pdf_array,  1280000,0);
  fitnt->Branch("best_fit_pdf_fit1",    "TObjArray",             &pdf_array1, 1280000,0);
  fitnt->Branch("best_fit_pdf_fit2",    "TObjArray",             &pdf_array2, 1280000,0);
  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  
  // Reset all output spectra
  dat_array ->Clear();
  pdf_array ->Clear();
  pdf_array1->Clear();
  pdf_array2->Clear();
  
  // Reset input true params and add true oscillation params
  iparm->Clear();
  PhysParamList true_oscillation_params = fparm1->AsPhysParamList("current");
  iparm->Append(true_oscillation_params);
  
  double mixing13_min = 0.0;
  double mixing13_max = 0.1;
  double mixing13_step = 0.005;

  double mixing23_min = 0.3;
  double mixing23_max = 0.7;
  double mixing23_step = 0.01;

  double dmsq_min = 2.0E-3;
  double dmsq_max = 3.0E-3;
  double dmsq_step = 0.05E-3;

  double deltacp_min = -conv::kpi;
  double deltacp_max = conv::kpi;
  double deltacp_step = (2*conv::kpi)/50.;

  //loop for the first parameter (if any)
  for (int i = min_bin_param1; i <= max_bin_param1; i++){
    if(param1 == "mixing13") {
      fparm1->SetInitValue(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
      fparm1->Update(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
      fparm1->FixParameter(conv::kPhypMixing13);
      fparm2->SetInitValue(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
      fparm2->Update(conv::kPhypMixing13,(mixing13_min + (i-1)*mixing13_step));
      fparm2->FixParameter(conv::kPhypMixing13);
    }
    if(param1 == "mixing23") {
      fparm1->SetInitValue(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
      fparm1->Update(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
      fparm1->FixParameter(conv::kPhypMixing23);
      fparm2->SetInitValue(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
      fparm2->Update(conv::kPhypMixing23,(mixing23_min + (i-1)*mixing23_step));
      fparm2->FixParameter(conv::kPhypMixing23);
    }
    if(param1 == "dmsq") {
      fparm1->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
      fparm1->Update(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
      fparm1->FixParameter(conv::kPhypDm23sq);
      fparm2->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
      fparm2->Update(conv::kPhypDm23sq,(dmsq_min + (i-1)*dmsq_step));
      fparm2->FixParameter(conv::kPhypDm23sq);
    }
    if(param1 == "deltacp") {
      fparm1->SetInitValue(conv::kPhypCPV    ,(deltacp_min + (i-1)*deltacp_step));
      fparm1->Update(conv::kPhypCPV    ,(deltacp_min + (i-1)*deltacp_step));
      fparm1->FixParameter(conv::kPhypCPV    );
      fparm2->SetInitValue(conv::kPhypCPV    ,(deltacp_min + (i-1)*deltacp_step));
      fparm2->Update(conv::kPhypCPV    ,(deltacp_min + (i-1)*deltacp_step));
      fparm2->FixParameter(conv::kPhypCPV    );
    }
    //loop for the second parameter (if any)
    for (int j = min_bin_param2; j <= max_bin_param2; j++){
      if(param2 == "mixing13") {
	fparm1->SetInitValue(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	fparm1->Update(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	fparm1->FixParameter(conv::kPhypMixing13);
	fparm2->SetInitValue(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	fparm2->Update(conv::kPhypMixing13,(mixing13_min + (j-1)*mixing13_step));
	fparm2->FixParameter(conv::kPhypMixing13);
      }
      if(param2 == "mixing23") {
	fparm1->SetInitValue(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	fparm1->Update(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	fparm1->FixParameter(conv::kPhypMixing23);
	fparm2->SetInitValue(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	fparm2->Update(conv::kPhypMixing23,(mixing23_min + (j-1)*mixing23_step));
	fparm2->FixParameter(conv::kPhypMixing23);
      }
      if(param2 == "dmsq") {
	fparm1->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	fparm1->Update(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	fparm1->FixParameter(conv::kPhypDm23sq);
	fparm2->SetInitValue(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	fparm2->Update(conv::kPhypDm23sq,(dmsq_min + (j-1)*dmsq_step));
	fparm2->FixParameter(conv::kPhypDm23sq);
      }
      if(param2 == "deltacp") {
	fparm1->SetInitValue(conv::kPhypCPV    ,(deltacp_min + (j-1)*deltacp_step));
	fparm1->Update(conv::kPhypCPV    ,(deltacp_min + (j-1)*deltacp_step));
	fparm1->FixParameter(conv::kPhypCPV    );
	fparm2->SetInitValue(conv::kPhypCPV    ,(deltacp_min + (j-1)*deltacp_step));
	fparm2->Update(conv::kPhypCPV    ,(deltacp_min + (j-1)*deltacp_step));
	fparm2->FixParameter(conv::kPhypCPV    );
      }
    
      //if theta23 is not fixed then initialize to 0.4 for the first fit (first octant) and 0.6 for the second one
      if(!(fparm1->IsFixed(conv::kPhypMixing23))) {     
	fparm1->SetInitValue(conv::kPhypMixing23,0.4);
	fparm1->Update(conv::kPhypMixing23,0.4);
	fparm2->SetInitValue(conv::kPhypMixing23,0.6);
	fparm2->Update(conv::kPhypMixing23,0.6);
      }
      //if reactor constraint enabled, set nuisance
      if(reactor) {
	fparm1->SetNuisance(conv::kPhypMixing13,true);
	fparm2->SetNuisance(conv::kPhypMixing13,true);
      }
      //reset systematic variations, otherwise fit is initialized at
      //best fit values of systematics from previous fit
      skexp->ResetSystematicVariations();
      
      // Get fitter
      FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
      FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
      
      // Chisq for true (input) params
      if(fparm1->Exists(conv::kNuisXSecSpecFunc))
	fparm1->Update(conv::kNuisXSecSpecFunc, 0); //init value of SpecFunc is 0.5
      double true_chisq = fit_inputs->Statistic();
      
      // Print-out true params
      //      fparm1->ShowOnlyCurrValuesInNextPrint();
      cout << "- True oscillation params (below) yield chisq = " << true_chisq
	   << *fparm1;
      
      // Do the fit
      fparm1->ResetCurrToInitValues(); // reset osc params
      fitter->SetTol(tolerance);
      /*      fitter->DoTheFit();
      fitstatus1 = fitter->FitStatusCode();
      fittime1 = fitter->FitStatusCode();
      chisq1  = fit_inputs->BestFitChisquare();*/
     // Save spectra for best-fit prediction
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }

     //now, only if theta23 is not fixed, repeat the fit at mirror point...
     if((Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) && (!(fparm1->IsFixed(conv::kPhypMixing23))) ) {
       skexp->ResetSystematicVariations();
       
       // Get fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
       FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
       
       // Print-out true params
       //       fparm2->ShowOnlyCurrValuesInNextPrint();
       cout << "- Second fit param list: " 
	    << *fparm2;
       
       // Do the fit
       fparm2->ResetCurrToInitValues(); // reset osc params
       fitter2->SetTol(tolerance);
       /*       fitter2->DoTheFit();
       fitstatus2 = fitter2->FitStatusCode();
       fittime2 = fitter2->FitStatusCode();
       chisq2  = fit_inputs2->BestFitChisquare();
       */
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
	 skexp->Write(pdf_array2, "pdf");
       }
     }
     
     // Print-out best-fit params, and set the "best of the two fits" branches correctly
     cout << "- Best-fit oscillation params (below) yield chisq = ";
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta  && (!(fparm1->IsFixed(conv::kPhypMixing23)))) {
       if(chisq2 < chisq1) {
	 cout << chisq2 << " from the second seed fit" << endl;
	 fparm2->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2;
	 br_fparm->SetAddress(&fparm2);
	 br_bf_array->SetAddress(&pdf_array2);
       }
       else {
	 cout << chisq1 << " from the first seed fit" << endl;
	 fparm1->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1;
	 br_fparm->SetAddress(&fparm1);
	 br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << chisq1 << endl;
       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1;
       br_fparm->SetAddress(&fparm1);
       br_bf_array->SetAddress(&pdf_array1);
     }

     // Add all results in the tree
     fitnt->Fill();

     // Clean-up
     delete fit_inputs;
     delete fitter;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta  && (!(fparm1->IsFixed(conv::kPhypMixing23)))) {
       delete fit_inputs2;
       delete fitter2;
     }

    } //j 
  } // i
  
  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();
}
