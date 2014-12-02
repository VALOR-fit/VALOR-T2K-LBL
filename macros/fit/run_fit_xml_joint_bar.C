//
// ------------------------------------------------------------------------------------------------
// Run fits to experiments read from an input XML file and store best-fit parameters and spectra.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//
 
void run_fit_xml_joint_bar
(
 const char * xmlfile = "fakedata/201411.fakedata/fakedata_201411_D.xml",                      /* input XML file Run5-6 fakedata */ 

 //const char * xmlfile = "data_newstruct.xml",                                                /* input XML file Run1-4 data */
 //const char * xmlfile = "data_newstruct_rhc.xml",                                         /* input XML file Run1-4 data (used by RHC analysis) */

 //const char * xmlfile = "data_newstruct_multiring.xml",                                                /* input XML file Run1-4 data (used by RHC analysis) */
 //const char * xmlfile = "toymc/toymc_FHCpRHC_DataRun1to4_MCRun5and6_th13DB_NH-0.xml",                 /* input XML file with toy MC experiments */
 //const char * xmlfile = "toymc_FHCpRHCpMultiRing-0.xml",
 //const char * xmlfile = "toymc_RHCpMultiRing_DataRun1to4_MCRun5and6_th13DB_NH-0.xml",

 const char * tag = "prova",                                                    /* */
 int          run = 0,                                                    /* */
 int          first_experiment = 0,                                   /* */
 int          last_experiment  = 0,                                   /* */
 
 const char * config_file      = "analysis_configuration_run1to5_test.xml",        /* analysis configuration file */
 //const char * config_file      = "analysis_configuration_test.xml",        /* analysis configuration file */
 //const char * config_file      = "analysis_configuration_test_rhc.xml",        /* analysis configuration file */
  //const char * config_file      = "analysis_configuration_test_multiring.xml",        /* analysis configuration file */
 //const char * config_file      = "analysis_configuration_test_multiring_mrmc.xml",        /* analysis configuration file */
 
 const char * fit_param_list          = "nuosc3fbar_fixnonbar_12fix",
 //const char * fit_param_list   = "nuosc3f_12fix+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
 //const char * fit_param_list   = "nuosc3f_12fix",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
 //const char * fit_param_list          = "nuosc3fbar_fixnonbar_12fix+nuisance",

 const char * fit_option       = "minuit",                            /* fitter option */
 //const char * fit_option       = "minuit+hesse",                            /* fitter option */
 
 bool         store_toy_expts         = false,                         /* disable this if you don't want to store the toy MC expts     */
 bool         store_bestfit_spectrum  = false,                         /* disable this if you don't want to store the best-fit spectra */
 //const char * dmsqscheme          = "fogli",                               /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * dmsqscheme          = "abs23",                               /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy           = "normal",                              /* normal or inverted */
  const char * mixparamchoice      = "sinsqtheta",                          /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant      = "one",                                 /* one==kNuOctantOne two==kNuOctantTwo */
  //bool         apply_2nu_boundary  = true                                   /* if using the 2nu flavour osc prob, apply the physical boundary? */
  double       dcp_value           = 0,                                      /* value of deltacp (in units of PI) to assume (when deltacp is not a free parameter) */ 
  bool         using_xml_file      = true,                                   /* if false, skobs reads xmlfile as a root file */
  bool         reactor             = true,                               /* enable this if you want to apply the reactor constraint */
  double       tolerance           = 0.1,                                 /* value for the fitter tolerance */

 TString optfixparams         = "mixing13, dmsq12"  // list of the params to be fixed in the fit   
 //TString optfixparams         = "mixing23, dmsq, dmsq12, mixing13, deltacp"  // list of the params to be fixed in the fit   
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  double bestfit_chisq,bestfit_chisq2;

  // Analysis configuration
  TString config   = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);
  
  FileStat_t dummy;
  
  // Check for config file, and XML data file.
  // If not found in inputs or config, see if it works as a normal path
  if (gSystem->GetPathInfo(config.Data(), dummy)) 
  {
    config = config_file;
    if (gSystem->GetPathInfo(config.Data(), dummy)) {
      cerr<<"Config file not found: "<<config_file<<endl;
      exit(-1);
    }
  }
  
  if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
  {
    xmlfile2 = xmlfile;
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
      cerr<<"Data input file not found: "<<xmlfile<<endl;
      exit(-1);
    }
  }

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
  
  Opt::Instance()->SetDefaultParameterValues();
  
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
  cout << "   - Fitter tolerance: " << tolerance << endl;
  cout << "  Applying reactor constraint: " << ((reactor)  ? "YES" : "NO" ) << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;
  

  // Fit param list
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = util::NamedParameterList(fit_param_list);
  FitParamList * fparm2 = util::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  //make an assumption for deltacp
  if(hyp==kHypOsc3F && fparm->IsFixed(conv::kPhypCPV)) {
    fparm ->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm1->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
    fparm2->SetInitValue(conv::kPhypCPV    , dcp_value * conv::kpi);
   }
  else if(hyp==kHypOsc3FBar && fparm->IsFixed(conv::kPhypCPVBar)) {
    fparm ->SetInitValue(conv::kPhypCPVBar    , dcp_value_bar * conv::kpi);
    fparm1->SetInitValue(conv::kPhypCPVBar    , dcp_value_bar * conv::kpi);
    fparm2->SetInitValue(conv::kPhypCPVBar    , dcp_value_bar * conv::kpi);
  }
  
  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

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
  TFile f(Form("fit-%s-run_%d-expts_%d_%d.root", tag, run, first_experiment, last_experiment), "recreate");
 
  double fittime1, fittime2;
  int fitstatus1, fitstatus2, fitstatus;
  int fitstatus_err1, fitstatus_err2, fitstatus_err;
  double chisq1, chisq2, chisq;

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("best_fit_params_fit1", "valor::FitParamList",  &fparm1,     1280000,0);
  fitnt->Branch("best_fit_params_fit2", "valor::FitParamList",  &fparm2,     1280000,0);
  fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("data",                 "TObjArray",             &dat_array,  1280000,0);
  fitnt->Branch("best_fit_pdf",         "TObjArray",             &pdf_array,  1280000,0);
  fitnt->Branch("best_fit_pdf_fit1",    "TObjArray",             &pdf_array1, 1280000,0);
  fitnt->Branch("best_fit_pdf_fit2",    "TObjArray",             &pdf_array2, 1280000,0);
  fitnt->Branch("fittime1",&fittime1,"fittime1/D");
  fitnt->Branch("fittime2",&fittime2,"fittime2/D");
  fitnt->Branch("fitstatus1",&fitstatus1,"fitstatus1/I");
  fitnt->Branch("fitstatus2",&fitstatus2,"fitstatus2/I");
  fitnt->Branch("fitstatus",&fitstatus,"fitstatus/I");
  fitnt->Branch("fitstatus_err1",&fitstatus_err1,"fitstatus_err1/I");
  fitnt->Branch("fitstatus_err2",&fitstatus_err2,"fitstatus_err2/I");
  fitnt->Branch("fitstatus_err",&fitstatus_err,"fitstatus_err/I");
  fitnt->Branch("chisq1",&chisq1,"chisq1/D");
  fitnt->Branch("chisq2",&chisq2,"chisq2/D");
  fitnt->Branch("chisq",&chisq,"chisq/D");

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




  ////////////////////////
  //                    //
  //   Fix parameters   //
  //                    //
  ////////////////////////
  
  if(hyp == kHypOsc3F){	
    if(optfixparams.Contains("mixing23")){ 
      fparm->FixParameter(kPhypMixing23);
      fparm1->FixParameter(kPhypMixing23);
      fparm2->FixParameter(kPhypMixing23);
    }
    if(optfixparams.Contains("mixing13")){ 
      fparm->FixParameter(kPhypMixing13);
      fparm1->FixParameter(kPhypMixing13);
      fparm2->FixParameter(kPhypMixing13);
    }
    if(optfixparams.Contains("dmsq")){ 
      fparm->FixParameter(kPhypDm23sq);
      fparm1->FixParameter(kPhypDm23sq);
      fparm2->FixParameter(kPhypDm23sq);
    }
    if(optfixparams.Contains("dmsq12")){ 
      fparm->FixParameter(kPhypDm12sq);
      fparm1->FixParameter(kPhypDm12sq);
      fparm2->FixParameter(kPhypDm12sq);
    }
    if(optfixparams.Contains("deltacp")){ 
      fparm->FixParameter(kPhypCPV);
      fparm1->FixParameter(kPhypCPV);
      fparm2->FixParameter(kPhypCPV);
    }
  }
  else if(hyp == kHypOsc3FBar){	
    if(optfixparams.Contains("mixing23")){ 
      fparm->FixParameter(kPhypMixing23Bar);
      fparm1->FixParameter(kPhypMixing23Bar);
      fparm2->FixParameter(kPhypMixing23Bar);
    }
    if(optfixparams.Contains("mixing13")){ 
      fparm->FixParameter(kPhypMixing13Bar);
      fparm1->FixParameter(kPhypMixing13Bar);
      fparm2->FixParameter(kPhypMixing13Bar);
    }
    if(optfixparams.Contains("dmsq")){ 
      fparm->FixParameter(kPhypDm23sqBar);
      fparm1->FixParameter(kPhypDm23sqBar);
      fparm2->FixParameter(kPhypDm23sqBar);
    }
    if(optfixparams.Contains("dmsq12")){ 
      fparm->FixParameter(kPhypDm12sqBar);
      fparm1->FixParameter(kPhypDm12sqBar);
      fparm2->FixParameter(kPhypDm12sqBar);
    }
    if(optfixparams.Contains("deltacp")){ 
      fparm->FixParameter(kPhypCPVBar);
      fparm1->FixParameter(kPhypCPVBar);
      fparm2->FixParameter(kPhypCPVBar);
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
 
     //if theta23 is not fixed then initialize to 0.4 for the first fit (first octant) and 0.6 for the second one
 
     if(hyp == kHypOsc3F) {
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
     }
     else if(hyp == kHypOsc3FBar){
       if(!(fparm1->IsFixed(conv::kPhypMixing23Bar))) {     
	 fparm1->SetInitValue(conv::kPhypMixing23Bar,0.4);
	 fparm1->Update(conv::kPhypMixing23Bar,0.4);
	 fparm2->SetInitValue(conv::kPhypMixing23Bar,0.6);
	 fparm2->Update(conv::kPhypMixing23Bar,0.6);
       }
       //if reactor constraint enabled, set nuisance
       if(reactor) {
	 fparm1->SetNuisance(conv::kPhypMixing13Bar,true);
	 fparm2->SetNuisance(conv::kPhypMixing13Bar,true);
       }
     }

     //reset systematic variations, otherwise fit is initialized at
     //best fit values of systematics from previous fit
     skexp->ResetSystematicVariations();
     
     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

     // Set SF Initial value to 0
     if(fparm1->Exists(conv::kNuisXSecSpecFunc)){
       fparm1->Update(conv::kNuisXSecSpecFunc, 0); //init value of SpecFunc is 0.5
     }

     // Chisq for true (input) params
     double true_chisq = fit_inputs->Statistic();
     
     // Print-out true params
     fparm1->ShowOnlyCurrValuesInNextPrint();
     
     cout << "- True oscillation params (below) yield chisq = " << true_chisq
	  << *fparm1;
     
     // Do the fit

     cout << endl;
     cout << endl;
     cout << "First Octant" << endl;
     cout << endl;
     cout << endl;

     fparm1->ResetCurrToInitValues(); // reset osc params
     
     skexp->PrintSystematicVariations();
     fitter->SetTol(tolerance);
     fitter->DoTheFit();
     fitstatus1 = fitter->FitStatusCode("migrad_min");
     fitstatus_err1 = fitter->FitStatusCode("hesse_err");     

     double bestfit_chisq  = fparm1->GetBestFitStatistic();
     chisq1 = bestfit_chisq;
     double bestfit_chisq2 = 0;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }


     //Do the fit again
     if( Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta &&
	 (!(fparm2->IsFixed(conv::kPhypMixing23)) && hyp == kHypOsc3F) || 
	 (!(fparm2->IsFixed(conv::kPhypMixing23Bar)) && hyp == kHypOsc3FBar)
	 ){
       
       cout << endl;
       cout << endl;
       cout << "Second Octant" << endl;
       cout << endl;
       cout << endl;

       skexp->ResetSystematicVariations();
       
       //get the correct intial values (best fit points from first fit)
       //fparm2->CopyBestFitToInitValues(*fparm1); // SF=0 not needed!!!
       fparm2->ResetCurrToInitValues(); // reset osc params
       if(fparm2->Exists(conv::kNuisXSecSpecFunc))
	 fparm2->Update(conv::kNuisXSecSpecFunc, 0); //init value of SpecFunc is 0.5

       cout << *fparm2 << endl;
       
       // Get fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
       FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
       fitter2->SetTol(tolerance);
       fitter2->DoTheFit();
       fitstatus2 = fitter->FitStatusCode("migrad_min");
       fitstatus_err2 = fitter->FitStatusCode("hesse_err");     
       bestfit_chisq2 = fparm2->GetBestFitStatistic();
       chisq2 = bestfit_chisq2;

       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
	 skexp->Write(pdf_array2, "pdf");
       }
     }

     // Print-out best-fit params, and set the "best of the two fits" branches correctly

     cout << "- Best-fit oscillation params (below) yield chisq = ";
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if(bestfit_chisq2 < bestfit_chisq) {

	 fitstatus = fitstatus2;
	 fitstatus_err = fitstatus_err2;
	 chisq = chisq2;

	 cout << bestfit_chisq2 << " from the second seed fit" << endl;
	 fparm2->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2;
	 br_fparm->SetAddress(&fparm2);
	 br_bf_array->SetAddress(&pdf_array2);
       }
       else {

	 fitstatus = fitstatus1;
	 fitstatus_err = fitstatus_err1;
	 chisq = chisq1;

	 cout << bestfit_chisq << " from the first seed fit" << endl;
	 fparm1->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1;
	 br_fparm->SetAddress(&fparm1);
	 br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << bestfit_chisq << endl;

       fitstatus = fitstatus1;
       fitstatus_err = fitstatus_err1;
       chisq = chisq1;

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
