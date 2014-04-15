//
// ------------------------------------------------------------------------------------------------
//
// Macro to create deltachisq 2D distribution for different toy experiments selected with certain seeds
// scanning all bins in theta23 for the joint analysis with systematics
// 
// ChangeLog:                                                                                                                              
// @ April 2013, Lorena
//
void chisq_vs_theta23_dm23sq
(
 const char * tag      = "sensitivity_point0_chisq_vs_theta23_dm23sq",                                  /* input SK data */
 const char * config_file    = "analysis_configuration.xml",                       /* analysis configuration file */
 const char * fit_param_list = "nuosc3f_12fix+nuisance",                          /* fit parameter list - options above  */
 const char * fit_option     = "minuit",                                           /* fitter option */
 //set initial values to different from default ones
 const char * hypothesis         = "3f_osc",                                       /* options: "2f_numu_disap", "3f_osc" */
 double       mixing12     = 0.311,                                         /* value from nue appearance TN52 */
 double       mixing13     = 0.025,                                           /* */
 double       mixing23     = 0.5,                                           /* value from nue appearance TN52 */
 double       deltacp      = 0.00,                                           /* */
 double       dm12sq       = 7.6E-5,                                         /* value from nue disapearance TN52 */
 double       dm23sq       = 2.36E-3 ,                                        /* value from nue disapearance TN52 */
 int          min_bin_theta23       = 1 ,                                        /* value from nue disapearance TN52 */
 int          max_bin_theta23       = 1 ,                                        /* value from nue disapearance TN52 */
 int          min_bin_dm23sq       = 1 ,                                        /* value from nue disapearance TN52 */
 int          max_bin_dm23sq       = 1 ,                                        /* value from nue disapearance TN52 */
 int          seed0       = 1,
 bool         normal       = true
)
{
  using namespace std;
  using namespace t2kana;
  
  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  if (normal) 
    Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  else 
    Opt::Instance()->gMassHierarchyAssumed = conv::kNmhInverted;
  Opt::Instance()->gOscDmsqScheme        = conv::kNuDmsqFogliLisi;

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  
  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);


  const char * outfile0 = Form("%s.root",tag);
  TFile f0(outfile0,"recreate");
  TTree * fit_tree = new TTree("fit_tree","fit_tree");
  
  //define variables for the tree
  int iexp;
  int fitstatus,fitstatus_disap, migradfitstatus;
  double fittime, fittime_disap;
  double bestfit_chisq,bestfit_disap_chisq;
  double true_chisq;
  double nsk_nom_1Rmu, nsk_bf_1Rmu, nsk_bf_disap_1Rmu;
  double nsk_nom_1Re, nsk_bf_1Re, nsk_bf_disap_1Re;

  double input_mixing12, bestfit_mixing12, bestfit_disap_mixing12 ;
  double input_mixing13, bestfit_mixing13, bestfit_disap_mixing13 ;
  double input_mixing23, bestfit_mixing23, bestfit_disap_mixing23 ;
  double input_deltacp, bestfit_deltacp, bestfit_disap_deltacp ;
  double input_dm12sq, bestfit_dm12sq, bestfit_disap_dm12sq;
  double input_dm23sq, bestfit_dm23sq, bestfit_disap_dm23sq;
 
  
 //set the branches for the tree
 fit_tree->Branch("iexp",&iexp,"iexp/I");
 fit_tree->Branch("fitstatus",&fitstatus,"fitstatus/I");
 fit_tree->Branch("fitstatus_disap",&fitstatus_disap,"fitstatus_disap/I");
 fit_tree->Branch("migradfitstatus",&migradfitstatus,"migradfitstatus/I");
 fit_tree->Branch("fittime",&fittime,"fittime/D");
 fit_tree->Branch("fittime_disap",&fittime_disap,"fittime_disap/D");
 fit_tree->Branch("true_chisq",&true_chisq,"true_chisq/D");
 fit_tree->Branch("bestfit_chisq",&bestfit_chisq,"bestfit_chisq/D");
 fit_tree->Branch("bestfit_disap_chisq",&bestfit_disap_chisq,"bestfit_disap_chisq/D");
 fit_tree->Branch("nsk_nom_1Rmu",&nsk_nom_1Rmu,"nsk_nom_1Rmu/D");
 fit_tree->Branch("nsk_nom_1Re",&nsk_nom_1Re,"nsk_nom_1Re/D");
 fit_tree->Branch("nsk_bf_1Rmu",&nsk_bf_1Rmu,"nsk_bf_1Rmu/D");
 fit_tree->Branch("nsk_bf_1Re",&nsk_bf_1Re,"nsk_bf_1Re/D");
 fit_tree->Branch("nsk_bf_disap_1Rmu",&nsk_bf_disap_1Rmu,"nsk_bf_disap_1Rmu/D");
 fit_tree->Branch("nsk_bf_disap_1Re",&nsk_bf_disap_1Re,"nsk_bf_disap_1Re/D");
  //inputs                                                                                                                                      
 fit_tree->Branch("input_mixing12",&input_mixing12,"input_mixing12/D");
 fit_tree->Branch("input_mixing13",&input_mixing13,"input_mixing13/D");
 fit_tree->Branch("input_mixing23",&input_mixing23,"input_mixing23/D");
 fit_tree->Branch("input_deltacp",&input_deltacp,"input_deltacp/D");
 fit_tree->Branch("input_dm12sq",&input_dm12sq,"input_dm12sq/D");
 fit_tree->Branch("input_dm23sq",&input_dm23sq,"input_dm23sq/D");
   //best fit                                                                                                                         
 fit_tree->Branch("bestfit_mixing12",&bestfit_mixing12,"bestfit_mixing12/D");
 fit_tree->Branch("bestfit_mixing13",&bestfit_mixing13,"bestfit_mixing13/D");
 fit_tree->Branch("bestfit_mixing23",&bestfit_mixing23,"bestfit_mixing23/D");
 fit_tree->Branch("bestfit_deltacp",&bestfit_deltacp,"bestfit_deltacp/D");
 fit_tree->Branch("bestfit_dm12sq",&bestfit_dm12sq,"bestfit_dm12sq/D");
 fit_tree->Branch("bestfit_dm23sq",&bestfit_dm23sq,"bestfit_dm23sq/D");

   //best fit                                                                                                                         
 fit_tree->Branch("bestfit_disap_mixing12",&bestfit_disap_mixing12,"bestfit_disap_mixing12/D");
 fit_tree->Branch("bestfit_disap_mixing13",&bestfit_disap_mixing13,"bestfit_disap_mixing13/D");
 fit_tree->Branch("bestfit_disap_mixing23",&bestfit_disap_mixing23,"bestfit_disap_mixing23/D");
 fit_tree->Branch("bestfit_disap_deltacp",&bestfit_disap_deltacp,"bestfit_disap_deltacp/D");
 fit_tree->Branch("bestfit_disap_dm12sq",&bestfit_disap_dm12sq,"bestfit_disap_dm12sq/D");
 fit_tree->Branch("bestfit_disap_dm23sq",&bestfit_disap_dm23sq,"bestfit_disap_dm23sq/D");

 // 
 // This macro is running only 1 toy experiment at once, scanning all bins 
 //
 const char * fit_param_list_disap = "nuosc3f_numudisapp_fix4+nuisance";

 int seed = seed0;

 for (int exp = 1; exp <21; i++){

   // Get specified fit parameter list & build fit inputs
   FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
   fparm->Update(conv::kPhypMixing12, mixing12);
   fparm->Update(conv::kPhypMixing13, mixing13);
   fparm->Update(conv::kPhypMixing23, mixing23);
   fparm->Update(conv::kPhypCPV    ,  deltacp);
   fparm->Update(conv::kPhypDm12sq,   dm12sq);
   fparm->Update(conv::kPhypDm23sq,   dm23sq);
   
   Rnd::Instance()->SetSeeds(seed++,seed++,seed++,seed++,seed++);
 
   // Get superk spectrum observation
   SuperkObservation * skobs = new SuperkObservation;
   skexp->TweakAllSystematics();
   skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), true); //with statistical fluct
   //   skobs->Read(data_file, dataset_id);

   //input values
   input_mixing12 = mixing12;
   input_mixing13 = mixing13;
   input_mixing23 = mixing23;
   input_dm12sq = dm12sq;
   input_dm23sq = dm23sq;
   input_deltacp = deltacp;
 
   double deltacp_init = deltacp;
   double mixing13_init = mixing13;
   
   // reset all systematics so that the PDF used in the fitter is the nominal one                                                                  
   cout << "- Reseting all systematic params in pdf..." << endl;
   skexp->ResetSystematicVariations();
   
   double mixing23_min = 0.3;
   double mixing23_max = 0.7;
   double step_mixing23 = 0.01;

   double dm23sq_min = 2.0E-3;
   double dm23sq_max = 3.0E-3;
   double step_dm23sq = 0.05E-3;
   
   
   for (int i = min_bin_theta23; i<= max_bin_theta23; i++){
     double mixing23_bin = mixing23_min + (i-1)*step_mixing23;
     
     for (int j = min_bin_dm23sq; j<= max_bin_dm23sq; j++){
       double dm23sq_bin = dm23sq_min + (j-1)*step_dm23sq;
       

       // cout << "mixing23_bin = " << mixing23_bin << endl;
       Opt::Instance()->gOscAnaTypeChoice     = conv::kAnaCombined;
       // FIRST FIT - JOINT
       FitParamList * fparm0 = FitParamList::NamedParameterList(fit_param_list);
       fparm0->SetInitValue(conv::kPhypMixing23, mixing23_bin);
       fparm0->Update(conv::kPhypMixing23, mixing23_bin);
       fparm0->FixParameter(conv::kPhypMixing23);
       fparm0->SetInitValue(conv::kPhypDm23sq, dm23sq_bin);
       fparm0->Update(conv::kPhypDm23sq, dm23sq_bin);
       fparm0->FixParameter(conv::kPhypDm23sq);
       // fparm0->ResetCurrToInitValues(); // reset osc params     
       fparm0->Update(conv::kPhypCPV    , deltacp_init);
       fparm0->SetInitValue(conv::kPhypCPV    , deltacp_init);
       fparm0->SetInitValue(conv::kPhypMixing13, mixing13_init);
       fparm0->Update(conv::kPhypMixing13, mixing13_init);
       fparm0->Print();
       
       FitInputs * fit_inputs0 = new FitInputs(skobs,skexp,fparm0);
       true_chisq = fit_inputs0->Chisquare();
       nsk_nom_1Rmu = 0;
       for(int ir = 0; ir < conv::kNErecoBins1Rmu; ir++)
	 nsk_nom_1Rmu += skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingMuLike)[ir];
       nsk_nom_1Re = 0;
       for(int ir = 0; ir < conv::kNErecoBins1Re; ir++)
	 nsk_nom_1Re += skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingELike)[ir];
       
       FitterABC * fitter0 = new FitterMinuit(fit_inputs0,fit_option);
       fitter0->SetStr(1);
       fitter0->SetTol(1.0);
       fitter0->DoTheFit();
       
       fitstatus       = fitter0->FitStatusCode();
       migradfitstatus = fitter0->FitMigradStatusCode();
       fittime = fitter0->FitTime();
       bestfit_chisq = fit_inputs0->BestFitChisquare();
       //   cout << "Chisq = " << bestfit_chisq << endl;
       nsk_bf_1Rmu = 0;
       for(int ir = 0; ir < conv::kNErecoBins1Rmu; ir++)
	 nsk_bf_1Rmu += skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingMuLike)[ir];
       nsk_bf_1Re = 0;
       for(int ir = 0; ir < conv::kNErecoBins1Re; ir++)
	 nsk_bf_1Re += skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingELike)[ir];
       
       bestfit_mixing12 = fparm0->BestFitValue(conv::kPhypMixing12);
       bestfit_mixing13 = fparm0->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing23 = fparm0->BestFitValue(conv::kPhypMixing23);
       bestfit_deltacp  = fparm0->BestFitValue(conv::kPhypCPV    );
       bestfit_dm12sq   = fparm0->BestFitValue(conv::kPhypDm12sq);
       bestfit_dm23sq   = fparm0->BestFitValue(conv::kPhypDm23sq);
   
       deltacp_init = fparm0->BestFitValue(conv::kPhypCPV    );
       mixing13_init = fparm0->BestFitValue(conv::kPhypMixing13);
       
       //Disappearance only:
       Opt::Instance()->gOscAnaTypeChoice     = conv::kAnaNuMuDisappOnly;
       // SECOND FIT - DISAP
       FitParamList * fparm_disap = FitParamList::NamedParameterList(fit_param_list_disap);
       fparm_disap->SetInitValue(conv::kPhypMixing23, mixing23_bin);
       fparm_disap->Update(conv::kPhypMixing23, mixing23_bin);
       fparm_disap->FixParameter(conv::kPhypMixing23);
       fparm_disap->SetInitValue(conv::kPhypDm23sq, dm23sq_bin);
       fparm_disap->Update(conv::kPhypDm23sq, dm23sq_bin);
       fparm_disap->FixParameter(conv::kPhypDm23sq);
       //       fparm_disap->ResetCurrToInitValues(); // reset osc params  
       fparm_disap->Print();                                                                                                    
       
       FitInputs * fit_inputs_disap = new FitInputs(skobs,skexp,fparm_disap);
       FitterABC * fitter_disap = new FitterMinuit(fit_inputs_disap,fit_option);
       fitter_disap->SetStr(1);
       fitter_disap->SetTol(1.0);
       fitter_disap->DoTheFit();
       
       fitstatus_disap       = fitter_disap->FitStatusCode();
       fittime_disap = fitter_disap->FitTime();
       bestfit_disap_chisq = fit_inputs_disap->BestFitChisquare();
       //   cout << "Chisq = " << bestfit_chisq << endl;
       nsk_bf_disap_1Rmu = 0;
       for(int ir = 0; ir < conv::kNErecoBins1Rmu; ir++)
	 nsk_bf_disap_1Rmu += skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingMuLike)[ir];
       
       bestfit_disap_mixing12 = fparm_disap->BestFitValue(conv::kPhypMixing12);
       bestfit_disap_mixing13 = fparm_disap->BestFitValue(conv::kPhypMixing13);
       bestfit_disap_mixing23 = fparm_disap->BestFitValue(conv::kPhypMixing23);
       bestfit_disap_deltacp  = fparm_disap->BestFitValue(conv::kPhypCPV    );
       bestfit_disap_dm12sq   = fparm_disap->BestFitValue(conv::kPhypDm12sq);
       bestfit_disap_dm23sq   = fparm_disap->BestFitValue(conv::kPhypDm23sq);
       
       fit_tree->Fill();
       
     }
   }
 }
 f0.cd();
 fit_tree->Write("fit_tree");
 f0.Close();
 cout << "Done!" << endl;
}
