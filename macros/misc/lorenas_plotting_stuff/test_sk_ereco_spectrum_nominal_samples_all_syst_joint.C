//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

void test_sk_ereco_spectrum_nominal_samples_all_syst_joint
(
  const char * hypothesis           = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12       = -99,                            /* */
  double       input_mixing13       = -99,                            /* */
  double       input_mixing23       = -99,                            /* */
  double       input_deltacp        = -99,                            /* */
  double       input_dm12sq         = -99,                            /* */
  double       input_dm23sq         = -99,                            /* */
  bool         disable_banff_tuning = false,                          /* */
  const char * config_file          = "analysis_configuration.xml"   /* */
 )

{
  using namespace std;
  using namespace t2kana;

 // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);  

  conv::SuperkSample_t s_mu = conv::kSKsFCFV1RingMuLike;    
  conv::SuperkSample_t s_e = conv::kSKsFCFV1RingELike;                                                                                                                        
  int idx_mu = util::ArrIdx::Instance()->SuperkSampleIdx(s_mu);                                                                                                                
  int idx_e = util::ArrIdx::Instance()->SuperkSampleIdx(s_e);                                                                                                                  
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  Opt::Instance()->gOscAnaTypeChoice     = conv::kAnaCombined;
  Opt::Instance()->gOscDmsqScheme        = conv::kNuDmsqFogliLisi;

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  {
    conv::NuOscDmsqScheme_t odmsq = Opt::Instance()->gOscDmsqScheme;                                                                                                                                    
    switch (odmsq) {                                                                                                                                                                                        
    case (conv::kNuDmsqFogliLisi): 
      input_dm23sq  =  conv::kGlobFitDmsqFL;
      break;                                                                                                                                                                                                
    case (conv::kNuDmsq12andAbs23):            
      input_dm23sq = conv::kGlobFitDm23sq;
      break;                                                                                                                                                                                                
    }                       
  }
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
  PhysHypothesis ph(hypothesis);
  
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }

  //  for(int im = 0; im < conv::knm; im++) {
  // conv::SuperkMode_t m = conv::kSKReactionModes[idx_mu][im];
  //if(m == conv::kSKmUndefined) continue;  
  //cout << "idx_mu =" << idx_mu << endl;
  //cout << "im =" << im << endl;
  //cout << "nskmode_0[idx_mu][im] = " <<  nskmode_0[idx_mu][im]  << endl;
  //}

  const char * fit_param_list      = "nuosc3f_12fix";

  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  fparm->SetInitValue(conv::kPhypDeltaCP, input_deltacp);
  fparm->Update(conv::kPhypDeltaCP, input_deltacp);
  fparm->SetInitValue(conv::kPhypMixing13, input_mixing13);
  fparm->Update(conv::kPhypMixing13, input_mixing13);

  fparm->SetInitValue(conv::kPhypMixing12, input_mixing12);
  fparm->SetInitValue(conv::kPhypDm12sq,   input_dm12sq);
  fparm->Update(conv::kPhypMixing12, input_mixing12);
  fparm->Update(conv::kPhypDm12sq,   input_dm12sq);
  
  fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
  fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
  fparm->Update(conv::kPhypMixing23, input_mixing23);
  fparm->Update(conv::kPhypDm23sq,   input_dm23sq);      
  
    
 // Create & config superk spectrum predictor
  SuperkExpectation * skexp  = new SuperkExpectation;
  skexp->Read(config);
  
  // Disable BANFF tuning if requested
  if(disable_banff_tuning) {
    cout << "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..." << endl;
    skexp.GetBANFF()->RestoreNominalMC();
  }

  // Generate prediction
  cout << "- Generating predictions..." << endl;
  skexp.Generate(ph);
  SuperkObservation * skobs = new SuperkObservation;
  bool applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), false);

  //Book histograms
  const int ntweaks = 4;
  int * tweaks[ntweaks] = {+3, +1, -1, -3};
  // int * tweaks[ntweaks] = {+1, -1};

  TH1D * hErecoTotNominal[conv::kns];
  TH1D * hErecoModeNominal[conv::kns][conv::knm];
  TH1D * hErecoTot[conv::kns][ntweaks];
  TH1D * hErecoRatio[conv::kns][ntweaks];
  TH1D * hErecoMode[conv::kns][conv::knm][ntweaks];
  TH1D * hErecoModeRatio[conv::kns][conv::knm][ntweaks];
  
  double nsk_0[conv::kns] = {0}, nsk_p[conv::kns] = {0}, nsk_m[conv::kns] = {0};
  double frac_nsk_p[conv::kns] = {0}, frac_nsk_m[conv::kns] = {0};

  double nskmode_0[conv::kns][conv::knm] = {0}, nskmode_p[conv::kns][conv::knm] = {0}, nskmode_m[conv::kns][conv::knm] = {0};
  double frac_nskmode_p[conv::kns][conv::knm] = {0}, frac_nskmode_m[conv::kns][conv::knm] = {0};
  
  double nsk_0_numu[conv::kns]={0}, nsk_p_numu[conv::kns]={0}, nsk_m_numu[conv::kns]={0};
  double nsk_0_oscnumu[conv::kns]={0}, nsk_p_oscnumu[conv::kns]={0}, nsk_m_oscnumu[conv::kns]={0};
  double nsk_0_numubar[conv::kns]={0}, nsk_p_numubar[conv::kns]={0}, nsk_m_numubar[conv::kns]={0};
  double nsk_0_nue[conv::kns]={0}, nsk_p_nue[conv::kns]={0}, nsk_m_nue[conv::kns]={0};
  double nsk_0_nuebar[conv::kns]={0}, nsk_p_nuebar[conv::kns]={0}, nsk_m_nuebar[conv::kns]={0};
  double nsk_0_oscnuebar[conv::kns]={0}, nsk_p_oscnuebar[conv::kns]={0}, nsk_m_oscnuebar[conv::kns]={0};
  double nsk_0_nuesig[conv::kns]={0}, nsk_p_nuesig[conv::kns]={0}, nsk_m_nuesig[conv::kns]={0};


  double frac_p_numu[conv::kns]={0}, frac_m_numu[conv::kns]={0};
  double frac_p_oscnumu[conv::kns]={0}, frac_m_oscnumu[conv::kns]={0};
  double frac_p_numubar[conv::kns]={0}, frac_m_numubar[conv::kns]={0};
  double frac_p_nue[conv::kns]={0}, frac_m_nue[conv::kns]={0};
  double frac_p_nuebar[conv::kns]={0}, frac_m_nuebar[conv::kns]={0};
  double frac_p_oscnuebar[conv::kns]={0}, frac_m_oscnuebar[conv::kns]={0};
  double frac_p_nuesig[conv::kns]={0}, frac_m_nuesig[conv::kns]={0};

  //loop over samples
  for(int is=0; is<conv::kns; is++) {
    cout << "is = " << is << endl;
    conv::SuperkSample_t s = conv::kSKSamples[is];
    if(s==conv::kSKmUndefined) continue;
    
    // Print options
    cout << "\n[Job inputs]" << endl;
    cout << "- Reading analysis configuration from: " << config << endl;
    cout << "- Plotting spectra for sample: " << util::AsString(s) << endl;
    cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
    cout << "- Physics hypothesis: " << ph << endl;
    cout << *Opt::Instance() << endl;
    
    // Histogram booking
    const int nbins = util::NumErecoBins(s);
    const double * ereco_bin_edges = util::ErecoBinEdges(s);
    cout << "nbins=" << nbins << endl;
    hErecoTotNominal[is] = new TH1D("","",nbins,ereco_bin_edges);
    
    //cout << "knm=" << conv::knm << endl;
    // Get energy spectra
    for(int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
      if(m == conv::kSKmUndefined) continue;
      cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) ;
      TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
      assert(h);
      cout << " with entries =" << h->GetEntries() << endl;
      hErecoModeNominal[is][im] = new TH1D(*h); 
      nskmode_0[is][im] = hErecoModeNominal[is][im]->Integral(); 
      if(util::IsNuMu(m)) nsk_0_numu[is] + = hErecoModeNominal[is][im]->Integral(); 
      if(util::IsOscNuMu(m)) nsk_0_oscnumu[is] + = hErecoModeNominal[is][im]->Integral(); 
      if(util::IsNuMuBar(m)) nsk_0_numubar[is] + = hErecoModeNominal[is][im]->Integral(); 
      if(util::IsNuE(m)) nsk_0_nue[is] + = hErecoModeNominal[is][im]->Integral(); 
      if(util::IsNuEBar(m)) nsk_0_nuebar[is] + = hErecoModeNominal[is][im]->Integral(); 
      if(util::IsOscNuEBar(m)) nsk_0_oscnuebar[is] + = hErecoModeNominal[is][im]->Integral(); 
      if(util::IsSignalNuE(m)) nsk_0_nuesig[is] + = hErecoModeNominal[is][im]->Integral(); 
    }
    hErecoTotNominal[is] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
    cout << "Total =" << hErecoTotNominal[is]->Integral() << " events. " << endl;
    //save here the numbers for the nominal
  } //end of first loop over samples, to save the nominal histograms
    nsk_0[0] = hErecoTotNominal[0]->Integral(); 
    nsk_0[1] = hErecoTotNominal[1]->Integral(); 
    cout << "Nominal done!" << endl;

  //TWEAKING ... 
  int nprecisionnev = 4;
  stringstream streamer;

  //Here it goes the case...     

  string modetitle ; 

  int i_fsi = 9;
  modetitle = Form("skdetfsi_%03d",i_fsi);
  //  int i_banff = 32;
  //modetitle = "sk_energy_scale";//"banff_nc1pi0";//Form("banff_%02d",i_banff); //
     
  conv::PhysParamId_t  p = util::PhysParamIdFromString(modetitle);
  /*
  TVectorD v(conv::kNumOfParamsInBANFFfit);
  v = skexp.GetBANFF()->BestFitValues();
  */

  TVectorD v(conv::kNumOfSKDetFSIParams);
  v.Zero();
  for(int ib = 0; ib < conv::kNumOfSKDetFSIParams; ib++) {
    v[ib] = 1;
  }

  FitParamList * fparm0 = FitParamList::NamedParameterList(fit_param_list);
  
  fparm0->SetInitValue(conv::kPhypMixing12, input_mixing12);
  fparm0->SetInitValue(conv::kPhypDm12sq,   input_dm12sq);
  fparm0->Update(conv::kPhypMixing12, input_mixing12);
  fparm0->Update(conv::kPhypDm12sq,   input_dm12sq);
  
  fparm0->SetInitValue(conv::kPhypMixing23, input_mixing23);
  fparm0->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
  fparm0->Update(conv::kPhypMixing23, input_mixing23);
  fparm0->Update(conv::kPhypDm23sq,   input_dm23sq);   

  // LOOP OVER TWEAKS 
  for(int itweak = 0; itweak < ntweaks; itweak++) {
    int tweak = tweaks[itweak];
    cout << endl << "Tweaking by " << tweak << " sigma for systematics: " << modetitle << endl;                                                                                 
      
      //set histogram color                                                                                                                                                        
      Color_t histcolor = kBlack;
      switch(itweak) {
	case 0:
	  histcolor = kRed;
	  break;
      case 1:
	histcolor = kBlue;
	break;
      case 2:
	histcolor = kGreen+3;
	  break;
      case 3:
	histcolor = kOrange+7;
	break;
      default:
	histcolor = kBlack;
      }
      
      skexp.TweakSystematic("xsec_sf",0); //by default the spectral function is at 0.5
      // skexp.TweakSystematic(modetitle.c_str(),tweak,true);
      /*
      v[i_banff] = skexp.GetBANFF()->BestFitValues()[i_banff] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[i_banff][i_banff]));
      skexp.GetBANFF()->Tweak(true,&v);
      */


      /*
      TVectorD v(conv::kNumOfFSIParams);
       v.Zero();
      for(int ib = 0; ib < conv::kNumOfFSIParams; ib++) {
	//if (ib == (i_fsi+6)) v[ib] = 1 + (tweak * sqrt(skexp.GetFSI()->CovMatrix()[ib][ib]));
	if (ib == i_fsi) v[ib] = 1 + (tweak * sqrt(skexp.GetFSI()->CovMatrix()[ib][ib]));
	else v[ib] = 1;
      }
      v.Print();
      skexp.GetFSI()->Tweak(true,&v);
      */

      v[i_fsi] = 1 + (tweak * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()[i_fsi][i_fsi]));
      //  cout << "penalty = " << ((v[i_fsi]-1)/skexp.GetSuperkDetectorFSI()->CovMatrix()[i_fsi][i_fsi])*((v[i_fsi]-1)/skexp.GetSuperkDetectorFSI()->CovMatrix()[i_fsi][i_fsi]);

      skexp.GetSuperkDetectorFSI()->Tweak(true,&v);

      //skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak);
      // skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak);
      //skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak);


      //skexp.TweakSystematic("sk_1Rmu_nuecc_eff",tweak);

      /*    
      TVectorD v(conv::kNumOfSKDetParams1Re);
      v.Zero();
      for(int ib = 0; ib < conv::kNumOfSKDetParams1Re; ib++) {
	if (ib == i_skdet) v[ib] = 1 + (tweak * sqrt(skexp.GetSuperkDetector()->Eff1ReCovarRelErr()[ib][ib]));
	//if (ib == i_fsi) v[ib] = 1 + (tweak * sqrt(skexp.GetFSI()->CovMatrix()[ib][ib]));
	else v[ib] = 1;
      }
      v.Print();
      skexp.GetSuperkDetector()->TweakEfficiencies1Re(true,&v);
*/
      
      skexp.PrintSystematicVariations();
      skexp.Generate(ph);    

      FitInputs * fit_inputs0 = new FitInputs(skobs,skexp,fparm0);
      double true_chisq = fit_inputs0->Chisquare();//("poisson",conv::kSKsFCFV1RingELike,-1);
      cout << "true_chisq  = " << true_chisq  << endl;
      double penalty = skexp->GetSuperkDetectorFSI()->ChisqPenalty();
      cout << "penalty = " << penalty << endl;

      // Print options
      //cout << "\n[Job inputs]" << endl;
      //cout << "- Reading analysis configuration from: " << config << endl;
      //cout << "- Plotting spectra for sample: " << util::AsString(s) << endl;
      //cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
      //cout << "- Physics hypothesis: " << ph << endl;
      //cout << *Opt::Instance() << endl;
      

      //loop over samples
      for(int is=0; is<conv::kns; is++) {
	cout << "sample = " << is << endl;
	conv::SuperkSample_t s = conv::kSKSamples[is];
	if(s==conv::kSKmUndefined) continue;
	// Histogram booking
	const int nbins = util::NumErecoBins(s);
	const double * ereco_bin_edges = util::ErecoBinEdges(s);
	//	cout << "nbins=" << nbins << endl;
	hErecoTot[is][itweak] = new TH1D("","",nbins,ereco_bin_edges);
	
	cout << "knm=" << conv::knm << endl;
	// Get energy spectra
	for(int im = 0; im < conv::knm; im++) {
	  // cout << "itweak = " << itweak << endl;
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  if(m == conv::kSKmUndefined) continue;
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) ;
	  TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  assert(h);
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  cout << " with entries =" << h->GetEntries() << endl;
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  hErecoMode[is][im][itweak] = new TH1D(*h);
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  hErecoModeRatio[is][im][itweak] = new TH1D(*h);
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  hErecoModeRatio[is][im][itweak]->Divide(hErecoModeNominal[is][im]);
	  //cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	  /*
	  if (itweak == 1) {
	    cout << "teak1" << endl;
	    nskmode_p[is][im] =     hErecoMode[is][im][itweak]->Integral();
	    frac_nskmode_p[is][im] = (nskmode_0[is][im] > 0) ? 100.*(nskmode_p[is][im] - nskmode_0[is][im]) / nskmode_0[is][im] : 0.;
	    if (util::IsNuMu(m))        nsk_p_numu[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsOscNuMu(m))     nsk_p_oscnumu[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsNuMuBar(m))     nsk_p_numubar[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsNuE(m))         nsk_p_nue[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsNuEBar(m))      nsk_p_nuebar[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsOscNuEBar(m))   nsk_p_oscnuebar[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsSignalNuE(m))   nsk_p_nuesig[is] + = hErecoMode[is][im][itweak]->Integral();
	  }
	  if  (itweak == 2) {
	    cout << "teak2" << endl;
	    nskmode_m[is][im] =     hErecoMode[is][im][itweak]->Integral();
	    frac_nskmode_m[is][im] = (nskmode_0[is][im] > 0) ? 100.*(nskmode_m[is][im] - nskmode_0[is][im]) / nskmode_0[is][im] : 0.;
	    if (util::IsNuMu(m))        nsk_m_numu[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsOscNuMu(m))     nsk_m_oscnumu[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsNuMuBar(m))     nsk_m_numubar[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsNuE(m))         nsk_m_nue[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsNuEBar(m))      nsk_m_nuebar[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsNuOscEBar(m))   nsk_m_oscnuebar[is] + = hErecoMode[is][im][itweak]->Integral();
	    if (util::IsSignalNuE(m))   nsk_m_nuesig[is] + = hErecoMode[is][im][itweak]->Integral();
	  }
	  */
	  //	  cout << "im=" << im << endl;
	  //cout << "itweak = " << itweak << endl;
	}//mode
	cout << "finished!" << endl;

	hErecoTot[is][itweak] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
	hErecoRatio[is][itweak] = new TH1D(*hErecoTot[is][itweak]);
	//	cout << "hErecoTot[is][itweak="<< itweak << "]=" << hErecoTot[is][itweak]->Integral() << endl;
	hErecoRatio[is][itweak]->Divide(hErecoTotNominal[is]);
	hErecoTot[is][itweak]->SetLineColor(histcolor);
	hErecoRatio[is][itweak]->SetLineColor(histcolor);
	
      }//sample

      
  }//tweak
  




  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m0 = conv::kSKReactionModes[0][im];
    conv::SuperkMode_t m1 = conv::kSKReactionModes[1][im];
    if(m0 != conv::kSKmUndefined){
      // cout << "teak1" << endl;
      nskmode_p[0][im] =     hErecoMode[0][im][1]->Integral();
      frac_nskmode_p[0][im] = (nskmode_0[0][im] > 0) ? 100.*(nskmode_p[0][im] - nskmode_0[0][im]) / nskmode_0[0][im] : 0.;
      if (util::IsNuMu(m0))        nsk_p_numu[0] + = hErecoMode[0][im][1]->Integral();
      if (util::IsOscNuMu(m0))     nsk_p_oscnumu[0] + = hErecoMode[0][im][1]->Integral();
      if (util::IsNuMuBar(m0))     nsk_p_numubar[0] + = hErecoMode[0][im][1]->Integral();
      if (util::IsNuE(m0))         nsk_p_nue[0] + = hErecoMode[0][im][1]->Integral();
      if (util::IsNuEBar(m0))      nsk_p_nuebar[0] + = hErecoMode[0][im][1]->Integral();
      if (util::IsOscNuEBar(m0))   nsk_p_oscnuebar[0] + = hErecoMode[0][im][1]->Integral();
      if (util::IsSignalNuE(m0))   nsk_p_nuesig[0] + = hErecoMode[0][im][1]->Integral();
    }
    if(m1 != conv::kSKmUndefined){
      //  cout << "teak1" << endl;
      nskmode_p[1][im] =     hErecoMode[1][im][1]->Integral();
      frac_nskmode_p[1][im] = (nskmode_0[1][im] > 0) ? 100.*(nskmode_p[1][im] - nskmode_0[1][im]) / nskmode_0[1][im] : 0.;
      if (util::IsNuMu(m1))        nsk_p_numu[1] + = hErecoMode[1][im][1]->Integral();
      if (util::IsOscNuMu(m1))     nsk_p_oscnumu[1] + = hErecoMode[1][im][1]->Integral();
      if (util::IsNuMuBar(m1))     nsk_p_numubar[1] + = hErecoMode[1][im][1]->Integral();
      if (util::IsNuE(m1))         nsk_p_nue[1] + = hErecoMode[1][im][1]->Integral();
      if (util::IsNuEBar(m1))      nsk_p_nuebar[1] + = hErecoMode[1][im][1]->Integral();
      if (util::IsOscNuEBar(m1))   nsk_p_oscnuebar[1] + = hErecoMode[1][im][1]->Integral();
      if (util::IsSignalNuE(m1))   nsk_p_nuesig[1] + = hErecoMode[1][im][1]->Integral();
    }
    if(m0 != conv::kSKmUndefined){
      // cout << "teak1" << endl;
      nskmode_m[0][im] =     hErecoMode[0][im][2]->Integral();
      frac_nskmode_m[0][im] = (nskmode_0[0][im] > 0) ? 100.*(nskmode_m[0][im] - nskmode_0[0][im]) / nskmode_0[0][im] : 0.;
      if (util::IsNuMu(m0))        nsk_m_numu[0] + = hErecoMode[0][im][2]->Integral();
      if (util::IsOscNuMu(m0))     nsk_m_oscnumu[0] + = hErecoMode[0][im][2]->Integral();
      if (util::IsNuMuBar(m0))     nsk_m_numubar[0] + = hErecoMode[0][im][2]->Integral();
      if (util::IsNuE(m0))         nsk_m_nue[0] + = hErecoMode[0][im][2]->Integral();
      if (util::IsNuEBar(m0))      nsk_m_nuebar[0] + = hErecoMode[0][im][2]->Integral();
      if (util::IsOscNuEBar(m0))   nsk_m_oscnuebar[0] + = hErecoMode[0][im][2]->Integral();
      if (util::IsSignalNuE(m0))   nsk_m_nuesig[0] + = hErecoMode[0][im][2]->Integral();
    }
    if(m1 != conv::kSKmUndefined){
      // cout << "teak1" << endl;
      nskmode_m[1][im] =     hErecoMode[1][im][2]->Integral();
      frac_nskmode_m[1][im] = (nskmode_0[1][im] > 0) ? 100.*(nskmode_m[1][im] - nskmode_0[1][im]) / nskmode_0[1][im] : 0.;
      if (util::IsNuMu(m1))        nsk_m_numu[1] + = hErecoMode[1][im][2]->Integral();
      if (util::IsOscNuMu(m1))     nsk_m_oscnumu[1] + = hErecoMode[1][im][2]->Integral();
      if (util::IsNuMuBar(m1))     nsk_m_numubar[1] + = hErecoMode[1][im][2]->Integral();
      if (util::IsNuE(m1))         nsk_m_nue[1] + = hErecoMode[1][im][2]->Integral();
      if (util::IsNuEBar(m1))      nsk_m_nuebar[1] + = hErecoMode[1][im][2]->Integral();
      if (util::IsOscNuEBar(m1))   nsk_m_oscnuebar[1] + = hErecoMode[1][im][2]->Integral();
      if (util::IsSignalNuE(m1))   nsk_m_nuesig[1] + = hErecoMode[1][im][2]->Integral();
    }
  }

  cout << "finished!" << endl;
  
  nsk_p[0] = hErecoTot[0][1]->Integral(); 
  nsk_p[1] = hErecoTot[1][1]->Integral(); 
  nsk_m[0] = hErecoTot[0][2]->Integral(); 
  nsk_m[1] = hErecoTot[1][2]->Integral(); 

      frac_nsk_p[0] = (nsk_0[0] > 0) ? 100.*(nsk_p[0] - nsk_0[0]) / nsk_0[0] : 0.;
      frac_nsk_p[1] = (nsk_0[1] > 0) ? 100.*(nsk_p[1] - nsk_0[1]) / nsk_0[1] : 0.;
      frac_nsk_m[0] = (nsk_0[0] > 0) ? 100.*(nsk_m[0] - nsk_0[0]) / nsk_0[0] : 0.;
      frac_nsk_m[1] = (nsk_0[1] > 0) ? 100.*(nsk_m[1] - nsk_0[1]) / nsk_0[1] : 0.;


    frac_p_numu[0] = (nsk_0_numu[0] > 0) ? 100.*(nsk_p_numu[0] - nsk_0_numu[0]) / nsk_0_numu[0] : 0.;
    frac_m_numu[0] = (nsk_0_numu[0] > 0) ? 100.*(nsk_m_numu[0] - nsk_0_numu[0]) / nsk_0_numu[0] : 0.;
    frac_p_numubar[0] = (nsk_0_numubar[0] > 0) ? 100.*(nsk_p_numubar[0] - nsk_0_numubar[0]) / nsk_0_numubar[0] : 0.;
    frac_m_numubar[0] = (nsk_0_numubar[0] > 0) ? 100.*(nsk_m_numubar[0] - nsk_0_numubar[0]) / nsk_0_numubar[0] : 0.;
    frac_p_oscnumu[0] = (nsk_0_oscnumu[0] > 0) ? 100.*(nsk_p_oscnumu[0] - nsk_0_oscnumu[0]) / nsk_0_oscnumu[0] : 0.;
    frac_m_oscnumu[0] = (nsk_0_oscnumu[0] > 0) ? 100.*(nsk_m_oscnumu[0] - nsk_0_oscnumu[0]) / nsk_0_oscnumu[0] : 0.;
    frac_p_nue[0] = (nsk_0_nue[0] > 0) ? 100.*(nsk_p_nue[0] - nsk_0_nue[0]) / nsk_0_nue[0] : 0.;
    frac_m_nue[0] = (nsk_0_nue[0] > 0) ? 100.*(nsk_m_nue[0] - nsk_0_nue[0]) / nsk_0_nue[0] : 0.;
    frac_p_nuebar[0] = (nsk_0_nuebar[0] > 0) ? 100.*(nsk_p_nuebar[0] - nsk_0_nuebar[0]) / nsk_0_nuebar[0] : 0.;
    frac_m_nuebar[0] = (nsk_0_nuebar[0] > 0) ? 100.*(nsk_m_nuebar[0] - nsk_0_nuebar[0]) / nsk_0_nuebar[0] : 0.;
    frac_p_oscnuebar[0] = (nsk_0_oscnuebar[0] > 0) ? 100.*(nsk_p_oscnuebar[0] - nsk_0_oscnuebar[0]) / nsk_0_oscnuebar[0] : 0.;
    frac_m_oscnuebar[0] = (nsk_0_oscnuebar[0] > 0) ? 100.*(nsk_m_oscnuebar[0] - nsk_0_oscnuebar[0]) / nsk_0_oscnuebar[0] : 0.;
    frac_p_nuesig[0] = (nsk_0_nuesig[0] > 0) ? 100.*(nsk_p_nuesig[0] - nsk_0_nuesig[0]) / nsk_0_nuesig[0] : 0.;
    frac_m_nuesig[0] = (nsk_0_nuesig[0] > 0) ? 100.*(nsk_m_nuesig[0] - nsk_0_nuesig[0]) / nsk_0_nuesig[0] : 0.;

    frac_p_numu[1] = (nsk_0_numu[1] > 0) ? 100.*(nsk_p_numu[1] - nsk_0_numu[1]) / nsk_0_numu[1] : 0.;
    frac_m_numu[1] = (nsk_0_numu[1] > 0) ? 100.*(nsk_m_numu[1] - nsk_0_numu[1]) / nsk_0_numu[1] : 0.;
    frac_p_oscnumu[1] = (nsk_0_oscnumu[1] > 0) ? 100.*(nsk_p_oscnumu[1] - nsk_0_oscnumu[1]) / nsk_0_oscnumu[1] : 0.;
    frac_m_oscnumu[1] = (nsk_0_oscnumu[1] > 0) ? 100.*(nsk_m_oscnumu[1] - nsk_0_oscnumu[1]) / nsk_0_oscnumu[1] : 0.;
    frac_p_numubar[1] = (nsk_0_numubar[1] > 0) ? 100.*(nsk_p_numubar[1] - nsk_0_numubar[1]) / nsk_0_numubar[1] : 0.;
    frac_m_numubar[1] = (nsk_0_numubar[1] > 0) ? 100.*(nsk_m_numubar[1] - nsk_0_numubar[1]) / nsk_0_numubar[1] : 0.;
    frac_p_nue[1] = (nsk_0_nue[1] > 0) ? 100.*(nsk_p_nue[1] - nsk_0_nue[1]) / nsk_0_nue[1] : 0.;
    frac_m_nue[1] = (nsk_0_nue[1] > 0) ? 100.*(nsk_m_nue[1] - nsk_0_nue[1]) / nsk_0_nue[1] : 0.;
    frac_p_oscnuebar[1] = (nsk_0_oscnuebar[1] > 0) ? 100.*(nsk_p_oscnuebar[1] - nsk_0_oscnuebar[1]) / nsk_0_oscnuebar[1] : 0.;
    frac_m_oscnuebar[1] = (nsk_0_oscnuebar[1] > 0) ? 100.*(nsk_m_oscnuebar[1] - nsk_0_oscnuebar[1]) / nsk_0_oscnuebar[1] : 0.;
    frac_p_nuebar[1] = (nsk_0_nuebar[1] > 0) ? 100.*(nsk_p_nuebar[1] - nsk_0_nuebar[1]) / nsk_0_nuebar[1] : 0.;
    frac_m_nuebar[1] = (nsk_0_nuebar[1] > 0) ? 100.*(nsk_m_nuebar[1] - nsk_0_nuebar[1]) / nsk_0_nuebar[1] : 0.;
    frac_p_nuesig[1] = (nsk_0_nuesig[1] > 0) ? 100.*(nsk_p_nuesig[1] - nsk_0_nuesig[1]) / nsk_0_nuesig[1] : 0.;
    frac_m_nuesig[1] = (nsk_0_nuesig[1] > 0) ? 100.*(nsk_m_nuesig[1] - nsk_0_nuesig[1]) / nsk_0_nuesig[1] : 0.;

    
    //write out latex: 
    //write down table - flavours
    //    cout << "Write first table" << endl;    

    streamer << "\\begin{table}"   << endl;
    streamer << "\\scriptsize" << endl;
    streamer << "\\begin{tabular} {|p{2cm}||p{1.5cm}|p{1.0cm}|p{1.0cm}|p{1.0cm}|p{1.0cm}||p{1.5cm}|p{1.0cm}|p{1.0cm}|p{1.0cm}|p{1.0cm}|}" << endl;
    streamer << "\\hline" << endl;
    streamer << "   & \\multicolumn {5}{c||} {$1R_{\\mu}$ }  &  \\multicolumn {5}{c|} {$1R_{e}$ } \\\\ " << endl;  
    streamer << "\\hline \\hline" << endl;
    streamer << "{\\bf Component} &  "
	     << "{\\bf BANFF tuned} &  "
	     << "\\multicolumn {2}{ c|} { {\\bf BANFF tuned +1$\\sigma$} } &  "
	     << "\\multicolumn {2}{ c||} { {\\bf BANFF tuned -1$\\sigma$} } &"
	     << "{\\bf BANFF tuned} &  "
	     << "\\multicolumn {2}{ c|} { {\\bf BANFF tuned +1$\\sigma$} } &  "
	     << "\\multicolumn {2}{ c|} { {\\bf BANFF tuned -1$\\sigma$} } \\\\"
	     << endl;
    streamer << "\\hline" << endl;
    streamer << "   &  $N_{SK}$ & $N_{SK}$ & change (\\%) & $N_{SK}$ & change (\\%) & $N_{SK}$ & $N_{SK}$ & change (\\%) & $N_{SK}$ & change (\\%)" << "\\\\" << endl;
    
    streamer << "\\hline" << endl;
    streamer << "{\\bf Total} & {\\bf" << setprecision(nprecisionnev) << nsk_0[idx_mu] << "} & {\\bf " << setprecision(nprecisionnev) <<  nsk_p[idx_mu]<< "}  & {\\bf" 
	     << setprecision(nprecisionnev) << frac_nsk_p[idx_mu] << "} & {\\bf" << setprecision(nprecisionnev) << nsk_m[idx_mu]<< "}  & {\\bf" 
	     << setprecision(nprecisionnev) << frac_nsk_m[idx_mu] << "} & {\\bf" 
	     << setprecision(nprecisionnev) << nsk_0[idx_e] << "} & {\\bf " << setprecision(nprecisionnev) <<  nsk_p[idx_e]<< " } & {\\bf" 
	     << setprecision(nprecisionnev) << frac_nsk_p[idx_e] << " } & {\\bf" << setprecision(nprecisionnev) << nsk_m[idx_e]<< " } & {\\bf" 
	     << setprecision(nprecisionnev) << frac_nsk_m[idx_e] << "} \\\\" << endl;
    streamer << "\\hline" << endl;
    streamer << "$\\nu_{\\mu}$ & " << setprecision(nprecisionnev) << nsk_0_numu[idx_mu] << " & " << setprecision(nprecisionnev) <<  nsk_p_numu[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_numu[idx_mu] << " & " << setprecision(nprecisionnev) << nsk_m_numu[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_numu[idx_mu] << " & " 
	     << setprecision(nprecisionnev) << nsk_0_numu[idx_e] << " & " << setprecision(nprecisionnev) <<  nsk_p_numu[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_numu[idx_e] << "  & " << setprecision(nprecisionnev) << nsk_m_numu[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_numu[idx_e] << " \\\\" << endl;
    streamer << "\\hline" << endl;
 streamer << "$Osc. \\nu_{\\mu}$ & " << setprecision(nprecisionnev) << nsk_0_oscnumu[idx_mu] << " & " << setprecision(nprecisionnev) <<  nsk_p_oscnumu[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_oscnumu[idx_mu] << " & " << setprecision(nprecisionnev) << nsk_m_oscnumu[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_oscnumu[idx_mu] << " & " 
	     << setprecision(nprecisionnev) << nsk_0_oscnumu[idx_e] << " & " << setprecision(nprecisionnev) <<  nsk_p_oscnumu[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_oscnumu[idx_e] << "  & " << setprecision(nprecisionnev) << nsk_m_oscnumu[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_oscnumu[idx_e] << " \\\\" << endl;
    streamer << "\\hline" << endl;
    streamer << "$\\bar{\\nu_{\\mu}}$ & " << setprecision(nprecisionnev) << nsk_0_numubar[idx_mu] << " & " << setprecision(nprecisionnev) <<  nsk_p_numubar[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_numubar[idx_mu] << " & " << setprecision(nprecisionnev) << nsk_m_numubar[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_numubar[idx_mu] << " & " 
	     << setprecision(nprecisionnev) << nsk_0_numubar[idx_e] << " & " << setprecision(nprecisionnev) <<  nsk_p_numubar[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_numubar[idx_e] << "  & " << setprecision(nprecisionnev) << nsk_m_numubar[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_numubar[idx_e] << " \\\\" << endl;
    streamer << "\\hline" << endl;
    streamer << "$\\nu_{e}$ & " << setprecision(nprecisionnev) << nsk_0_nue[idx_mu] << " & " << setprecision(nprecisionnev) <<  nsk_p_nue[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_nue[idx_mu] << " & " << setprecision(nprecisionnev) << nsk_m_nue[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_nue[idx_mu] << " & " 
	     << setprecision(nprecisionnev) << nsk_0_nue[idx_e] << " & " << setprecision(nprecisionnev) <<  nsk_p_nue[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_nue[idx_e] << "  & " << setprecision(nprecisionnev) << nsk_m_nue[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_nue[idx_e] << " \\\\" << endl;
    streamer << "\\hline" << endl;
    streamer << "$\\bar{\\nu_{e}}$ & " << setprecision(nprecisionnev) << nsk_0_nuebar[idx_mu] << " & " << setprecision(nprecisionnev) <<  nsk_p_nuebar[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_nuebar[idx_mu] << " & " << setprecision(nprecisionnev) << nsk_m_nuebar[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_nuebar[idx_mu] << " & " 
	     << setprecision(nprecisionnev) << nsk_0_nuebar[idx_e] << " & " << setprecision(nprecisionnev) <<  nsk_p_nuebar[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_nuebar[idx_e] << "  & " << setprecision(nprecisionnev) << nsk_m_nuebar[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_nuebar[idx_e] << " \\\\" << endl;
    streamer << "\\hline" << endl;
    streamer << "Osc.$\\bar{\\nu_{e}}$ & " << setprecision(nprecisionnev) << nsk_0_oscnuebar[idx_mu] << " & " << setprecision(nprecisionnev) <<  nsk_p_oscnuebar[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_oscnuebar[idx_mu] << " & " << setprecision(nprecisionnev) << nsk_m_oscnuebar[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_oscnuebar[idx_mu] << " & " 
	     << setprecision(nprecisionnev) << nsk_0_oscnuebar[idx_e] << " & " << setprecision(nprecisionnev) <<  nsk_p_oscnuebar[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_oscnuebar[idx_e] << "  & " << setprecision(nprecisionnev) << nsk_m_oscnuebar[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_oscnuebar[idx_e] << " \\\\" << endl;
    streamer << "\\hline" << endl;
    streamer << "$Osc. \\nu_{e}$ & " << setprecision(nprecisionnev) << nsk_0_nuesig[idx_mu] << " & " << setprecision(nprecisionnev) <<  nsk_p_nuesig[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_nuesig[idx_mu] << " & " << setprecision(nprecisionnev) << nsk_m_nuesig[idx_mu]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_nuesig[idx_mu] << " & " 
	     << setprecision(nprecisionnev) << nsk_0_nuesig[idx_e] << " & " << setprecision(nprecisionnev) <<  nsk_p_nuesig[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_p_nuesig[idx_e] << "  & " << setprecision(nprecisionnev) << nsk_m_nuesig[idx_e]<< "  & " 
	     << setprecision(nprecisionnev) << frac_m_nuesig[idx_e] << " \\\\" << endl;
    streamer << "\\hline" << endl;
    streamer << "\\end{tabular}" << endl;
    //  streamer << "\\caption{Effect of " << util::AsStringTex(p) << "  with no oscillation (all mixing angles set to 0) for each flavour.}"<< endl;
    streamer << "\\caption{Effect of " << util::AsStringTex(p) << "  with $\\Delta m^{2}$ = " << input_dm23sq << " , $\\delta m^{2}$ = " << input_dm12sq << " , $sin^{2}\\theta_{12}$ = " << input_mixing12 << " , $sin^{2}\\theta_{13}$ = " << input_mixing13 << " , $sin^{2}\\theta_{23}$ = " << input_mixing23 << " , $\\delta_{CP}$ = " << input_deltacp << " and normal hierarchy.}" << endl;
    streamer << "\\end{table}" << endl;
    
    //  cout << "Write second table" << endl;    
    
    
    //write down table - the big one
    
    streamer << "\\begin{table}"   << endl;
    streamer << "\\scriptsize" << endl;
    streamer << "\\begin{tabular} {|p{2.5cm}||p{1.2cm}|p{1.0cm}|p{1.0cm}|p{1.0cm}|p{1.0cm}||p{1.2cm}|p{1.0cm}|p{1.0cm}|p{1.0cm}|p{1.0cm}|}" << endl;
    streamer << "\\hline" << endl;
    streamer << "   & \\multicolumn {5}{c||} {$1R_{\\mu}$ }  &  \\multicolumn {5}{c|} {$1R_{e}$ } \\\\ " << endl;  
    streamer << "\\hline \\hline" << endl;
    streamer << "{\\bf Component} &  "
	     << "{\\bf BANFF tuned} &  "
	     << "\\multicolumn {2}{c|} { {\\bf BANFF tuned +1$\\sigma$} } &  "
	     << "\\multicolumn {2}{c||} { {\\bf BANFF tuned -1$\\sigma$} } &"
	     << "{\\bf BANFF tuned} &  "
	     << "\\multicolumn {2}{c|} { {\\bf BANFF tuned +1$\\sigma$} } &  "
	     << "\\multicolumn {2}{c|} { {\\bf BANFF tuned -1$\\sigma$} } \\\\"
	     << endl;
    streamer << "\\hline" << endl;
    streamer << "   &  $N_{SK}$ & $N_{SK}$ & change (\\%) & $N_{SK}$ & change (\\%) & $N_{SK}$ & $N_{SK}$ & change (\\%) & $N_{SK}$ & change (\\%)" << "\\\\" << endl;
    
    streamer << "\\hline" << endl;
    streamer << "{\\bf Total} & {\\bf" << setprecision(nprecisionnev) << nsk_0[idx_mu] << "} & {\\bf " << setprecision(nprecisionnev) <<  nsk_p[idx_mu]<< "}  & {\\bf" 
	     << setprecision(nprecisionnev) << frac_nsk_p[idx_mu] << "} & {\\bf" << setprecision(nprecisionnev) << nsk_m[idx_mu]<< "}  & {\\bf" 
	     << setprecision(nprecisionnev) << frac_nsk_m[idx_mu] << "} & {\\bf" 
	     << setprecision(nprecisionnev) << nsk_0[idx_e] << "} & {\\bf " << setprecision(nprecisionnev) <<  nsk_p[idx_e]<< " } & {\\bf" 
	     << setprecision(nprecisionnev) << frac_nsk_p[idx_e] << " } & {\\bf" << setprecision(nprecisionnev) << nsk_m[idx_e]<< " } & {\\bf" 
	     << setprecision(nprecisionnev) << frac_nsk_m[idx_e] << "} \\\\" << endl;
    streamer << "\\hline" << endl;
    
    //  cout << "Writing in second table" << endl;    
    
    //For each mode ... number and order are different for the two samples
    //for now, write each by hand
    //find more elegant way in the future
    
    //cout << " nskmode_0[idx_mu][2] =" << nskmode_0[idx_mu][2] << endl;
    //cout << " nskmode_p[idx_mu][2] =" << nskmode_p[idx_mu][2] << endl;
    
    streamer << util::AsStringTex(conv::kSKmNuMuCCQE)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][0] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][0]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][0] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][0]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][0] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][0] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][0]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][0] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][0]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][0] << " \\\\" << endl;
   
 streamer << util::AsStringTex(conv::kSKmNuMuCC1pi)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][1] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][1] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][1] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][1] << " \\\\" << endl;
    streamer << util::AsStringTex(conv::kSKmNuMuCCcoh)  << " &" // << setprecision(nprecisionnev) << 
	     << nskmode_0[idx_mu][2] << " &  " 
      //<< setprecision(nprecisionnev) 
	     <<  nskmode_p[idx_mu][2]<< "  & " 
      //  << setprecision(nprecisionnev)
	     << frac_nskmode_p[idx_mu][2] << " & " //<< setprecision(nprecisionnev) 
	     << nskmode_m[idx_mu][2]<< "  & " 
      //<< setprecision(nprecisionnev) 
	     << frac_nskmode_m[idx_mu][2] << " & " 
      //<< setprecision(nprecisionnev) 
	     << nskmode_0[idx_e][2] << " &  " 
      //<< setprecision(nprecisionnev) 
	     <<  nskmode_p[idx_e][2]<< "  & " 
      //<< setprecision(nprecisionnev) 
	     << frac_nskmode_p[idx_e][2] << "  & "
      //<< setprecision(nprecisionnev) 
	     << nskmode_m[idx_e][2]<< "  & " 
      //<< setprecision(nprecisionnev) 
	     << frac_nskmode_m[idx_e][2] << " \\\\" << endl;
    streamer << util::AsStringTex(conv::kSKmNuMuCCoth)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][3] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][3] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][3] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][3] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][3] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][3] << " \\\\" << endl;
    
    streamer  << util::AsStringTex(conv::kSKmNuMuNC1piC)  << " & " << setprecision(nprecisionnev) << nskmode_0[idx_mu][4] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4] << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<"   \\\\" << endl;  
    
    streamer  << util::AsStringTex(conv::kSKmNuMuNC1pi0)  << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<" & "
	      << setprecision(nprecisionnev) << nskmode_0[idx_e][4] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4] << " \\\\ "  <<endl;
    streamer  << util::AsStringTex(conv::kSKmNuMuNCcoh)  << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<" & "
	      << setprecision(nprecisionnev) << nskmode_0[idx_e][5] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_e][5]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][5] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_e][5]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][5] << " \\\\ "  <<endl;
    
    streamer << util::AsStringTex(conv::kSKmNuMuNCoth_e)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][5] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][5]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][5] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][5]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][5] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][5+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][5+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][5+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][5+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][5+1] << " \\\\" << endl;
    
    streamer << "\\hline" << endl; // END OF NUMU 
    
    //OSCNUMU 

    streamer << util::AsStringTex(conv::kSKmOscNuMuCCQE)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][6] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][6]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][6] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][6]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][6] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][6+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][6+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][6+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][6+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][6+1] << " \\\\" << endl;

 streamer << util::AsStringTex(conv::kSKmOscNuMuCC1pi)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][7] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][7]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][7] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][7]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][7] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][7+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][7+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][7+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][7+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][7+1] << " \\\\" << endl;

 streamer << util::AsStringTex(conv::kSKmOscNuMuCCcoh)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][8] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][8]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][8] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][8]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][8] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][8+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][8+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][8+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][8+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][8+1] << " \\\\" << endl;

 streamer << util::AsStringTex(conv::kSKmOscNuMuCCoth)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][9] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][9]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][9] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][9]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][9] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][9+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][9+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][9+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][9+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][9+1] << " \\\\" << endl;
    streamer << "\\hline" << endl;

    //NUMUBAR


    streamer << util::AsStringTex(conv::kSKmNuMuBarCCQE)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+6] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+6]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+6] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+6]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+6] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+6+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+6+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+6+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+6+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+6+1] << " \\\\" << endl;
    streamer << util::AsStringTex(conv::kSKmNuMuBarCC1pi)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+7] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+7]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+7] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+7]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+7] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+7+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+7+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+7+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+7+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+7+1] << " \\\\" << endl;
    streamer << util::AsStringTex(conv::kSKmNuMuBarCCcoh)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+8] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+8]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+8] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+8]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+8] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+8+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+8+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+8+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+8+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+8+1] << " \\\\" << endl;
    streamer << util::AsStringTex(conv::kSKmNuMuBarCCoth)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+9] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+9]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+9] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+9]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+9] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+9+1] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+9+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+9+1] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+9+1]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+9+1] << " \\\\" << endl;
    
    streamer  << util::AsStringTex(conv::kSKmNuMuBarNC1piC)  << " & " << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+10] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+10]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+10] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+10]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+10] << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<"   \\\\" << endl;  
    
    streamer  << util::AsStringTex(conv::kSKmNuMuBarNC1pi0)  << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<" & "
	      << setprecision(nprecisionnev) << nskmode_0[idx_e][4+11] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+11]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+11] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+11]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+11] << " \\\\ "  <<endl;
    streamer  << util::AsStringTex(conv::kSKmNuMuBarNCcoh)  << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<" & "
	      << setprecision(nprecisionnev) << nskmode_0[idx_e][4+12] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+12]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+12] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+12]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+12] << " \\\\ "  <<endl;
    
    streamer << util::AsStringTex(conv::kSKmNuMuBarNCoth_e)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+11] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+11]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+11] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+11]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+11] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+11+2] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+11+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+11+2] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+11+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+11+2] << " \\\\" << endl;
    
    streamer << "\\hline" << endl; // END OF NUMUBAR 
    streamer << util::AsStringTex(conv::kSKmNuECCQE)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+12] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+12]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+12] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+12]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+12] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+12+2] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+12+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+12+2] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+12+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+12+2] << " \\\\" << endl;
    streamer << util::AsStringTex(conv::kSKmNuECC1pi)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+13] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+13]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+13] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+13]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+13] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+13+2] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+13+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+13+2] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+13+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+13+2] << " \\\\" << endl;
    streamer << util::AsStringTex(conv::kSKmNuECCcoh)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+14] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+14]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+14] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+14]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+14] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+14+2] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+14+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+14+2] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+14+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+14+2] << " \\\\" << endl;
    streamer << util::AsStringTex(conv::kSKmNuECCoth)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+15] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+15]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+15] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+15]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+15] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+15+2] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+15+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+15+2] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+15+2]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+15+2] << " \\\\" << endl;
    
    streamer  << util::AsStringTex(conv::kSKmNuENC1piC)  << " & " << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+16] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+16]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+16] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+16]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+16] << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<"   \\\\" << endl;   
    streamer  << util::AsStringTex(conv::kSKmNuENC1pi0)  << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<" & "
	      << setprecision(nprecisionnev) << nskmode_0[idx_e][4+18] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+18]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+18] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+18]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+18] << " \\\\ "  <<endl;
    streamer  << util::AsStringTex(conv::kSKmNuENCcoh)  << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<" & "
	      << setprecision(nprecisionnev) << nskmode_0[idx_e][4+19] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+19]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+19] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+19]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+19] << " \\\\ "  <<endl; 
    streamer << util::AsStringTex(conv::kSKmNuENCoth_e)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+17] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+17]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+17] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+17]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+17] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+17+3] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+17+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+17+3] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+17+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+17+3] << " \\\\" << endl; 
    streamer << "\\hline" << endl; // END OF NUE 
    
    streamer << util::AsStringTex(conv::kSKmNuEBarCCQE)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+18] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+18]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+18] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+18]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+18] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+18+3] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+18+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+18+3] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+18+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+18+3] << " \\\\" << endl; 
    streamer << util::AsStringTex(conv::kSKmNuEBarCC1pi)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+19] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+19]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+19] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+19]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+19] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+19+3] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+19+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+19+3] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+19+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+19+3] << " \\\\" << endl; 
    streamer << util::AsStringTex(conv::kSKmNuEBarCCcoh)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+20] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+20]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+20] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+20]<< "  & " 
	   << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+20] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+20+3] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+20+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+20+3] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+20+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+20+3] << " \\\\" << endl; 
    streamer << util::AsStringTex(conv::kSKmNuEBarCCoth)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+21] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+21]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+21] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+21]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+21] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+21+3] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+21+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+21+3] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+21+3]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+21+3] << " \\\\" << endl; 
    
    streamer  << util::AsStringTex(conv::kSKmNuEBarNC1piC)  << " & " << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+22] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+22]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+22] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+22]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+22] << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<"   \\\\" << endl;   
    streamer  << util::AsStringTex(conv::kSKmNuEBarNC1pi0)  << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<" & "
	      << setprecision(nprecisionnev) << nskmode_0[idx_e][4+25] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+25]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+25] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+25]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+25] << " \\\\ "  <<endl;
    streamer  << util::AsStringTex(conv::kSKmNuEBarNCcoh)  << " & " << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "& "
	      << setw(3)  << setfill('-') <<" & " << setw(3)  << setfill('-') << "&" << setw(3)  << setfill('-') <<" & "
	      << setprecision(nprecisionnev) << nskmode_0[idx_e][4+26] << " & " 
	      << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+26]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+26] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+26]<< "  & " 
	      << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+26] << " \\\\ "  <<endl; 
    streamer << util::AsStringTex(conv::kSKmNuEBarNCoth_e)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][4+23] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][4+23]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][4+23] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][4+23]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][4+23] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][4+23+4] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][4+23+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][4+23+4] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][4+23+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][4+23+4] << " \\\\" << endl; 
    streamer << "\\hline" << endl; // END OF NUEBAR

    //OSC NUE BAR
    streamer << util::AsStringTex(conv::kSKmOscNuEBarCCQE)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][28] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][28]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][28] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][28]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][28] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][32] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][32]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][32] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][32]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][32] << " \\\\" << endl; 

    streamer << util::AsStringTex(conv::kSKmOscNuEBarCC1pi)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][29] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][29]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][29] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][29]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][29] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][33] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][33]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][33] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][33]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][33] << " \\\\" << endl; 

    streamer << util::AsStringTex(conv::kSKmOscNuEBarCCcoh)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][30] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][30]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][30] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][30]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][30] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][34] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][34]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][34] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][34]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][34] << " \\\\" << endl; 

  streamer << util::AsStringTex(conv::kSKmOscNuEBarCCoth)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][31] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][31]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][31] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][31]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][31] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][35] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][35]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][35] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][35]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][35] << " \\\\" << endl; 

    streamer << "\\hline" << endl;


    streamer << util::AsStringTex(conv::kSKmSigNuECCQE)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][8+24] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][8+24]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][8+24] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][8+24]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][8+24] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][8+24+4] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][8+24+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][8+24+4] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][8+24+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][8+24+4] << " \\\\" << endl; 
    streamer << util::AsStringTex(conv::kSKmSigNuECC1pi)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][8+25] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][8+25]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][8+25] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][8+25]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][8+25] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][8+25+4] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][8+25+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][8+25+4] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][8+25+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][8+25+4] << " \\\\" << endl; 
    streamer << util::AsStringTex(conv::kSKmSigNuECCcoh)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][8+26] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][8+26]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][8+26] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][8+26]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][8+26] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][8+26+4] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][8+26+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][8+26+4] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][8+26+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][8+26+4] << " \\\\" << endl; 
    streamer << util::AsStringTex(conv::kSKmSigNuECCoth)  << " &" << setprecision(nprecisionnev) << nskmode_0[idx_mu][8+27] << " &  " 
	     << setprecision(nprecisionnev) <<  nskmode_p[idx_mu][8+27]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_mu][8+27] << " & " << setprecision(nprecisionnev) << nskmode_m[idx_mu][8+27]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_mu][8+27] << " & " 
	     << setprecision(nprecisionnev) << nskmode_0[idx_e][8+27+4] << " &  " << setprecision(nprecisionnev) <<  nskmode_p[idx_e][8+27+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_p[idx_e][8+27+4] << "  & " << setprecision(nprecisionnev) << nskmode_m[idx_e][8+27+4]<< "  & " 
	     << setprecision(nprecisionnev) << frac_nskmode_m[idx_e][8+27+4] << " \\\\" << endl; 
    
    
    streamer << "\\hline" << endl;
    streamer << "\\end{tabular}" << endl;
    //    streamer << "\\caption{Effect of " << util::AsStringTex(p) << "  with no oscillation (all mixing angles set to 0)}"<< endl;
    streamer << "\\caption{Effect of " << util::AsStringTex(p) << "  with $\\Delta m^{2}$ = " << input_dm23sq << " , $\\delta m^{2}$ = " << input_dm12sq << " , $sin^{2}\\theta_{12}$ = " << input_mixing12 << " , $sin^{2}\\theta_{13}$ = " << input_mixing13 << " , $sin^{2}\\theta_{23}$ = " << input_mixing23 << " , $\\delta_{CP}$ = " << input_deltacp << " and normal hierarchy.}" << endl;
//<< "  with oscillation parameters set to PDG 2012 values}"<< endl;
    streamer << "\\end{table}" << endl;
    
    
    //write down first table - the big one
    
    
    //DRAW
    double ereco_max_plot = 7.0;
    
    
    TCanvas * c = new TCanvas("c","",20,20,1000,1000);                                                                                                                           
    gStyle->SetOptStat(0);                                                                                                                                                       
    c->Divide(2,2);                                                                                                                                                              
    //1Rmu ereco spectra                                                                                                                                                         
    c->cd(1);                                                                                                                                                                    
    TLegend * legend = new TLegend(0.7, 0.6, 0.9, 0.9);                                                                            
    legend->SetFillColor(0);  
    legend->SetBorderSize(0);                                                                                                         
    legend->AddEntry(hErecoTotNominal[idx_mu], "nominal", "L" );                                                                                                                    
    hErecoTotNominal[idx_mu]->SetTitle(Form("%s",modetitle.c_str()));      
    hErecoTotNominal[idx_mu]->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hErecoTotNominal[idx_mu]->GetYaxis()->SetTitle("Number of events");
    hErecoTotNominal[idx_mu]->GetXaxis()->SetRangeUser(0,ereco_max_plot-conv::kOneMeV);                                                                                             
    hErecoTotNominal[idx_mu]->GetYaxis()->SetRangeUser(0,8);
    hErecoTotNominal[idx_mu]->Draw();                                                                                                                                               
    for(int itweak = 0; itweak < ntweaks; itweak++) {                                                                                 
      //for(int itweak = 1; itweak < 2; itweak++) {                                                                                 
      hErecoTot[idx_mu][itweak]->Draw("same");                                                                                        
      legend->AddEntry(hErecoTot[idx_mu][itweak], Form("%d sigma",tweaks[itweak]), "L" );                                             
    }                                                                                                                                                                            
    legend->Draw();                                                                                                                                                              
    //1Re ereco spectra                                                                                                                                                          
    c->cd(2);                                                                                                                                                                    
    TLegend * legend = new TLegend(0.7, 0.6, 0.9, 0.9);                                                                              
    legend->SetFillColor(0);                                                                                                      
    legend->SetBorderSize(0);                                                                                                      
    legend->AddEntry(hErecoTotNominal[idx_e], "nominal", "L" );                                                                                                                     
    hErecoTotNominal[idx_e]->SetTitle(Form("%s",modetitle.c_str()));   
    hErecoTotNominal[idx_e]->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hErecoTotNominal[idx_e]->GetYaxis()->SetTitle("Number of events");                                                                                                             
    hErecoTotNominal[idx_e]->GetYaxis()->SetRangeUser(0,2.5);
    hErecoTotNominal[idx_e]->Draw();                                                                                                                                                
    for(int itweak = 0; itweak < ntweaks; itweak++) {                                                                                    
      //for(int itweak = 1; itweak < 2; itweak++) { 
      hErecoTot[idx_e][itweak]->Draw("same");                                                                                     
      legend->AddEntry(hErecoTot[idx_e][itweak], Form("%d sigma",tweaks[itweak]), "L" );                                            
    }                                                                                                                                                                            
    legend->Draw();                                                                                                                                                              
    //1Rmu ratio                                                                                                                                                                 
    c->cd(3);                                                                                                                                                                    
    for(int itweak = 0; itweak < ntweaks; itweak++) {
      //for(int itweak = 1; itweak < 2; itweak++) {                          
   if(itweak==0){                                                                                                                  
	hErecoRatio[idx_mu][itweak]->GetXaxis()->SetRangeUser(0,ereco_max_plot-conv::kOneMeV);                                        
	hErecoRatio[idx_mu][itweak]->GetYaxis()->SetRangeUser(0.3,1.7);                                                              
	hErecoRatio[idx_mu][itweak]->GetXaxis()->SetTitle("E_{reco} (GeV)");                                                          
	hErecoRatio[idx_mu][itweak]->GetYaxis()->SetTitle("Ratio to nominal");                                                        
	hErecoRatio[idx_mu][itweak]->SetTitle("Ratio");                                                                               
	hErecoRatio[idx_mu][itweak]->Draw();                                                                                          
	continue;                                                                                                                     
	}                                                                                                                                                                          
      hErecoRatio[idx_mu][itweak]->Draw("same");                                                                                                                                 
    }                                         
    //1Rmu ratio                               
    c->cd(4);                                  
    for(int itweak = 0; itweak < ntweaks; itweak++) { 
      // for(int itweak = 1; itweak < 2; itweak++) { 
      if(itweak==0){                                  
	hErecoRatio[idx_e][itweak]->GetXaxis()->SetTitle("E_{reco} (GeV)");                                                                                                        
	hErecoRatio[idx_e][itweak]->GetYaxis()->SetTitle("Ratio to nominal");
	hErecoRatio[idx_e][itweak]->GetYaxis()->SetRangeUser(0.3,1.7);                                                                                   
	hErecoRatio[idx_e][itweak]->SetTitle("Ratio");                                                                                                                             
	hErecoRatio[idx_e][itweak]->Draw();                                                                                                                                        
	continue;                                                                                                                                                                  
      }                                                                                                                                                                            
    hErecoRatio[idx_e][itweak]->Draw("same");                                                                                                                                    
    }                                                                                                                                                                              
    c->SaveAs(Form("ereco_twk_%s.eps",modetitle.c_str()));    

    //
    // include the .eps plot
    
    string eps_file (Form("ereco_twk_%s.eps",modetitle.c_str()));
    streamer << "\\begin{figure}[htbp]" << endl;
    streamer << "\\begin{center}" << endl;
    streamer << "\\includegraphics[width=0.9\\textwidth]{" << eps_file.c_str() <<"}" << endl;
    //streamer << "\\caption{Effect of " << util::AsStringTex(p) << "  with no oscillation (all mixing angles set to 0)}"<< endl;
    streamer << "\\caption{Effect of " << util::AsStringTex(p) << "  with $\\Delta m^{2}$ = " << input_dm23sq << " , $\\delta m^{2}$ = " << input_dm12sq << " , $sin^{2}\\theta_{12}$ = " << input_mixing12 << " , $sin^{2}\\theta_{13}$ = " << input_mixing13 << " , $sin^{2}\\theta_{23}$ = " << input_mixing23 << " , $\\delta_{CP}$ = " << input_deltacp << " and normal hierarchy.}" << endl; //<< endl;//"  with oscillation parameters set to PDG 2012 values}"<< endl;
    streamer << "\\end{center}" << endl;
    streamer << "\\end{figure}" << endl;
    streamer << "\\clearpage" << endl;
    
   // write in the output file
    // }


  //  cout << streamer.str() << endl;
  long size;
  // get size of file                                                                                                                                                             
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);
 
  
  string output_file ("systematic_tables_and_plots_2013_skfsi.tex");
  ofstream outfile (output_file.c_str(),ios_base::app);

  // write document header to outfile                  
  string f;
  f += "\\documentclass[10pt,a4paper]{article}\n";
  f += ("\\usepackage{amsmath}\n");
  f += ("\\usepackage{graphicx}\n");
  f += ("\\usepackage{dcolumn}\n");
  f += ("\\usepackage[cm]{fullpage}\n");
  f += ("\\begin{document}\n");
  outfile.write (f.c_str(),f.size());
   //outfile.close();
  
  // write to outfile                                  
  outfile.write (streamer.str().c_str(),size);
  
  //after loop over modes 
  string f;
  f += ("\\end{document}\n");
  outfile.write (f.c_str(),f.size());
  outfile.close();
  
}
 
