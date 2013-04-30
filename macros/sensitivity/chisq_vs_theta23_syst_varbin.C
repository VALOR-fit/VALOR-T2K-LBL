//
// ------------------------------------------------------------------------------------------------
//
// Macro to create deltachisq 1D distribution for one toy experiment, scanning all bins in theta23
// for the joint analysis with systematics and perform the fit as well for the standalone disappearance in the same toy
// 
// ChangeLog:                                                                                                                              
// @ 12/04/2013, Lorena
//
void chisq_vs_theta23_syst_varbin
(
 const char * tag      = "sensitivity_point0_0",                                  /* input SK data */
 const char * config_file    = "analysis_configuration.xml",                       /* analysis configuration file */
 const char * fit_param_list = "nuosc3f_12fix+nuisance",                          /* fit parameter list - options above  */
 const char * fit_option     = "minuit",                                           /* fitter option */
 //set initial values to different from default ones
 const char * hypothesis         = "3f_osc",                                       /* options: "2f_numu_disap", "3f_osc" */
 double       mixing12     = 0.311,                                         /* value from nue appearance TN52 */
 double       mixing13     = 0.025,                                           /* */
 double       mixing23     = 0.38,                                           /* value from nue appearance TN52 */
 double       deltacp      = 0.00,                                           /* */
 double       dm12sq       = 7.6E-5,                                         /* value from nue appearance TN52 */
 double       dm23sq       = 2.4E-3 ,                                        /* value from nue appearance TN52 */
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

 // Get specified fit parameter list & build fit inputs
 FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
 fparm->Update(conv::kPhypMixing12, mixing12);
 fparm->Update(conv::kPhypMixing13, mixing13);
 fparm->Update(conv::kPhypMixing23, mixing23);
 fparm->Update(conv::kPhypDeltaCP,  deltacp);
 fparm->Update(conv::kPhypDm12sq,   dm12sq);
 fparm->Update(conv::kPhypDm23sq,   dm23sq);
 
 const char * fit_param_list_disap = "nuosc3f_numudisapp_fix4+nuisance";

 // Get superk spectrum observation
 SuperkObservation * skobs = new SuperkObservation;
 //   skobs->Read(data_file, dataset_id);
 bool applied_negevents_condition = true;
 while(applied_negevents_condition) {
   skexp->TweakAllSystematics();
   applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), true); //with statistical fluct
 }
 
 skexp->PrintSystematicVariations();
 
 //input values
 input_mixing12 = mixing12;
 input_mixing13 = mixing13;
 input_mixing23 = mixing23;
 input_dm12sq = dm12sq;
 input_dm23sq = dm23sq;
 input_deltacp = deltacp;

 
 // reset all systematics so that the PDF used in the fitter is the nominal one                                                                  
 cout << "- Reseting all systematic params in pdf..." << endl;
 skexp->ResetSystematicVariations();
 
 double mixing23_min = 0.59;
 double mixing23_max = 0.61;
 double step = 0.01;
 int nbins = 41;
 double mixing23_bin = mixing23_min;
 bool first = true;
 double mixing23_oct = 1.0 - mixing23;

 double mixing13_init = mixing13;
 double dm23sq_init = dm23sq;
 double dm23sq_disap_init = dm23sq;
 double deltacp_init = deltacp;

 while(mixing23_bin < mixing23_max){
 // for (int i = 1; i<= nbins; i++){
   if (!first) {
     if ((TMath::Abs( mixing23_bin - mixing23) < 0.04) || (TMath::Abs( mixing23_bin - mixing23_oct) < 0.04)) 
       mixing23_bin += 0.005; 
     else if ((TMath::Abs( mixing23_bin - mixing23) < 0.08) || (TMath::Abs( mixing23_bin - mixing23_oct) < 0.08)) 
       mixing23_bin += 0.01; 
     else mixing23_bin += 0.02;
   }
   //cout << "mixing23_bin = " << mixing23_bin << endl;
   
   Opt::Instance()->gOscAnaTypeChoice     = conv::kAnaCombined;
   // FIRST FIT - JOINT
   FitParamList * fparm0 = FitParamList::NamedParameterList(fit_param_list);
   fparm0->SetInitValue(conv::kPhypMixing23, mixing23_bin);
   fparm0->Update(conv::kPhypMixing23, mixing23_bin);
   fparm0->FixParameter(conv::kPhypMixing23);
   //   fparm0->ResetCurrToInitValues(); // reset osc params 
   fparm0->Update(conv::kPhypMixing13, mixing13_init);
   fparm0->Update(conv::kPhypDm23sq, dm23sq_init);
   fparm0->Update(conv::kPhypDeltaCP, deltacp_init);
   fparm0->SetInitValue(conv::kPhypMixing13, mixing13_init);
   fparm0->SetInitValue(conv::kPhypDm23sq, dm23sq_init);
   fparm0->SetInitValue(conv::kPhypDeltaCP, deltacp_init);
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
   bestfit_deltacp  = fparm0->BestFitValue(conv::kPhypDeltaCP);
   bestfit_dm12sq   = fparm0->BestFitValue(conv::kPhypDm12sq);
   bestfit_dm23sq   = fparm0->BestFitValue(conv::kPhypDm23sq);

   mixing13_init = fparm0->BestFitValue(conv::kPhypMixing13);
   dm23sq_init = fparm0->BestFitValue(conv::kPhypDm23sq);
   deltacp_init  = fparm0->BestFitValue(conv::kPhypDeltaCP);

   
   //Disappearance only:
   Opt::Instance()->gOscAnaTypeChoice     = conv::kAnaNuMuDisappOnly;
   
   FitParamList * fparm_disap = FitParamList::NamedParameterList(fit_param_list_disap);
   fparm_disap->SetInitValue(conv::kPhypMixing23, mixing23_bin);
   fparm_disap->Update(conv::kPhypMixing23, mixing23_bin);
   fparm_disap->FixParameter(conv::kPhypMixing23);
   //fparm_disap->ResetCurrToInitValues(); // reset osc params
   fparm_disap->Update(conv::kPhypDm23sq, dm23sq_disap_init);
   fparm_disap->SetInitValue(conv::kPhypDm23sq, dm23sq_disap_init);
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
   bestfit_disap_deltacp  = fparm_disap->BestFitValue(conv::kPhypDeltaCP);
   bestfit_disap_dm12sq   = fparm_disap->BestFitValue(conv::kPhypDm12sq);
   bestfit_disap_dm23sq   = fparm_disap->BestFitValue(conv::kPhypDm23sq);
   
   dm23sq_disap_init = fparm_disap->BestFitValue(conv::kPhypDm23sq);
       
   first = false;
   fit_tree->Fill();
 }
 
 f0.cd();
 fit_tree->Write("fit_tree");
 f0.Close();
 cout << "Done!" << endl;
}
