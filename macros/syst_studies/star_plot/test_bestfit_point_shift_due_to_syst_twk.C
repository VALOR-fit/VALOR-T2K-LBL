//
// ------------------------------------------------------------------------------------------------
// Test best-fit point shift as a result of tweaking a single systematic parameter by +/- 1 sigma.
// Adapted from test_systematic_error_propagation.C
// The best-fit parameter shift is written out in a text file.
// One can assemble all such text files and produce the `star plot'.
//
// The input (systematic) string can be:
//
//      - osc_mixing12
//      - osc_mixing13
//      - osc_dm12sq
//      - osc_deltacp
//
//	- sk_energy_scale 
//	- sk_1Rmu_ccqe_other_eff
//	- sk_1Rmu_ccqe_ring_counting_eff : 
//             This energy-dependent correlated error can not be used in this test macro. Instead,
//             you can tweak specific components of the efficiency tweaking vector (neglecting,
//             correlation which, in the fitter, are handled using the Cholesky method). Use:
//	       - sk_1Rmu_ccqe_ring_counting_eff_0 
//	       - sk_1Rmu_ccqe_ring_counting_eff_1
//	       - sk_1Rmu_ccqe_ring_counting_eff_2
//	- sk_1Rmu_ccnonqe_eff
//	- sk_1Rmu_nc_eff
//	- sk_1Rmu_nuecc_eff
//      - banff
//             This correlated error can not be used in this test macro. Instead,
//             you can tweak specific components of the BANFF tweaking vector (neglecting,
//             correlation which, in the fitter, are handled using the Cholesky method). Use:
//             - banff_00
//             - ...
//             - banff_30
//             - banff_maqe
//             - banff_mares
//             - banff_ccqeE0
//             - banff_ccqeE1
//             - banff_ccqeE2
//             - banff_cc1piE0
//             - banff_cc1piE1
//             - banff_nc1pi0
//      - fsi
//             This correlated error can not be used in this test macro. Instead,
//             you can tweak specific components of the BANFF tweaking vector (neglecting,
//             correlation which, in the fitter, are handled using the Cholesky method). Use:
//             - fsi_00
//             - ...
//             - fsi_14
//      - xsec_pf
//      - xsec_sf
//      - xsec_wshape
//      - xsec_shape_ccoth
//      - xsec_pilessdelta
//      - xsec_enushape
//      - xsec_norm_ccqeE1
//      - xsec_norm_ccqeE2
//      - xsec_norm_cc1piE1
//      - xsec_norm_cccoh
//      - xsec_norm_nccoh
//      - xsec_norm_ncoth
//      - xsec_norm_nuetonumu
//      - xsec_norm_antinutonu
//      
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_bestfit_point_shift_due_to_syst_twk
(
 const char * systematic,                                                  /* options: see above */
 const char * sksample                   = "1Rmu",                         /* options: "1Rmu", "1Re" */
 const char * fit_param_list             = "nuosc3f_numudisapp_fix4",      /* */
 const char * fit_option                 = "minuit",                       /* oscillation fitter (see run_fitter.C) */
 double       input_mixing12             = -99,                            /* */
 double       input_mixing13             = -99,                            /* */
 double       input_mixing23             = -99,                            /* */
 double       input_deltacp              = -99,                            /* */
 double       input_dm12sq               = -99,                            /* */
 double       input_dm23sq               = -99,                            /* */
 const char * config_file                = "analysis_configuration.xml",   /* */
 bool         correlated_parameters      = false,                          /* tweak matrix parameters respecting correlations? */
 const int    throw_n_toys               = 0                               /* if <=0, will use +1 and -1 sigma tweaks. Otherwise, throw n toys & find standard deviation */
 )
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Get requested sample
  conv::SuperkSample_t s = util::SKSampleFromString(sksample);
  assert(s != conv::kSKsUndefined);

  //string for convieniance
  string sys(systematic);

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
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;

  // Instantiate fit params & set specified inputs as necessary
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;
  if(hyp == conv::kHypOsc3F) {
    fparm->SetInitValue(conv::kPhypMixing12, input_mixing12);
    fparm->SetInitValue(conv::kPhypMixing13, input_mixing13);
    fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
    fparm->SetInitValue(conv::kPhypDeltaCP,  input_deltacp);
    fparm->SetInitValue(conv::kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
  }
  else if(hyp == conv::kHypNuMuDisapp2F) {
    fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
    fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
  }  
  else {
    return;
  }

  bool fluctuate_stat = false; // <- don't fluctuate, best-fit point shift is purely due to systematic tweak

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << " - Reading analysis configuration from: " << config << endl;
  cout << " - SuperK sample: " << util::AsString(s) << endl;
  cout << " - Toy MC expt info: " << endl;
  cout << "    - Varied systematic parameter (by +/- 1 sigma): " << systematic << endl;
  cout << "    - Including statistical fluctuations? : " << fluctuate_stat << endl;
  cout << "    - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(conv::kPhypMixing12) << " = " << input_mixing12 << endl;
  cout << "       > " << util::AsString(conv::kPhypMixing13) << " = " << input_mixing13 << endl;
  cout << "       > " << util::AsString(conv::kPhypMixing23) << " = " << input_mixing23 << endl;
  cout << "       > " << util::AsString(conv::kPhypDeltaCP)  << " = " << input_deltacp << endl;
  cout << "       > " << util::AsString(conv::kPhypDm12sq)   << " = " << input_dm12sq << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(conv::kPhypDm23sq)   << " = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << " - Fit info: " << endl;
  cout << "    - Fit parameter list: " << fit_param_list << endl;
  cout << "    - Fitter option: " << fit_option << endl;
  cout << *Opt::Instance() << endl;

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  SuperkObservation * skobs = 0;
  FitInputs * fit_inputs = 0;
  FitterABC * fitter = 0;

  //initialise output
  double bestfit_mixing12_p1sig = -1;
  double bestfit_mixing13_p1sig = -1;
  double bestfit_mixing23_p1sig = -1;
  double bestfit_deltacp_p1sig  = -1;
  double bestfit_dm12sq_p1sig   = -1;
  double bestfit_dm23sq_p1sig   = -1;
  double bestfit_mixing12_m1sig = -1;
  double bestfit_mixing13_m1sig = -1;
  double bestfit_mixing23_m1sig = -1;
  double bestfit_deltacp_m1sig  = -1;
  double bestfit_dm12sq_m1sig   = -1;
  double bestfit_dm23sq_m1sig   = -1;
  
  if(throw_n_toys <= 0) {
    //
    //
    // +1 sigma (relative to tuned BANFF)
    //
    //
    
    cout << "\n\nApplying +1 sigma tweak" << endl;
    if(sys.find("banff_") == 0) {
      //find parameter to tweak
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
      //initialise tweak vector
      const int nbanff = conv::kNumOfParamsInBANFFfit;
      TVectorD twkv(nbanff);
      const TMatrixD & cov = skexp->GetBANFF()->CovMatrix();
      //tweak the parameter to tweak
      twkv[ip] = skexp->GetBANFF()->BestFitValue(p) + (1 * sqrt(cov[ip][ip]));
      //if we're tweaking the other parameters in a correlated way too
      double sigma = (twkv[ip] - skexp->GetBANFF()->BestFitValue(p)) / cov[ip][ip];
      for(int ipp = 0; ipp < nbanff; ipp++) {
	int position = (int)conv::kNuisBANFF00 + ipp;
	conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	if(pp == p)
	  continue;
	//set non-tweaked values to best fit
	twkv[ipp] = skexp->GetBANFF()->BestFitValue(pp);
	//now, tweak the other parameters according to the correlations between them and the fitted parameters
	if(correlated_parameters)
	  twkv[ipp] += sigma * cov[ip][ipp];
      }
      skexp->GetBANFF()->Tweak(true,&twkv);
    }
    else if(sys.find("skdetfsi_") == 0) {
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->SKDetFSISystIdx(p);
      //initialise tweak vector
      const int nskdetfsi = conv::kNumOfSKDetFSIParams;
      TVectorD twkv(nskdetfsi);
      const TMatrixD & cov = skexp->GetSuperkDetectorFSI()->CovMatrix();
      //tweak the parameter to tweak
      twkv[ip] = 1 + (1 * sqrt(cov[ip][ip]));
      //if we're tweaking the other parameters in a correlated way too
      double sigma = (twkv[ip] - 1) / cov[ip][ip];
      for(int ipp = 0; ipp < nskdetfsi; ipp++) {
	int position = (int)conv::kNuisSKDetFSI000 + ipp;
	conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	if(pp == p)
	  continue;
	//set non-tweaked values to best fit
	twkv[ipp] = 1;
	//now, tweak the other parameters according to the correlations between them and the fitted parameters
	if(correlated_parameters)
	  twkv[ipp] += sigma * cov[ip][ipp];
      }
      skexp->GetSuperkDetectorFSI()->Tweak(true,&twkv);
    }
    else if(sys.find("fsi_") == 0) {
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->FSISystIdx(p);
      //initialise tweak vector
      const int nfsi = conv::kNumOfFSIParams;
      TVectorD twkv(nfsi);
      const TMatrixD & cov = skexp->GetFSI()->CovMatrix();
      //tweak the parameter to tweak
      twkv[ip] = 1 + (1 * sqrt(cov[ip][ip]));
      //if we're tweaking the other parameters in a correlated way too
      double sigma = (twkv[ip] - 1) / cov[ip][ip];
      for(int ipp = 0; ipp < nfsi; ipp++) {
	int position = (int)conv::kNuisFSI000 + ipp;
	conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	if(pp == p)
	  continue;
	//set non-tweaked values to best fit
	twkv[ipp] = 1;
	//now, tweak the other parameters according to the correlations between them and the fitted parameters
	if(correlated_parameters)
	  twkv[ipp] += sigma * cov[ip][ipp];
      }
      skexp->GetFSI()->Tweak(true,&twkv);
    }
    else if(sys.find("osc_mixing12") == 0) {
      double tweak_mixing12 = input_mixing12 * (1 + conv::kGlobFitErrSinsq2Theta12);
      if(tweak_mixing12 >  fparm->MaxValue(conv::kPhypMixing12)) tweak_mixing12 = fparm->MaxValue(conv::kPhypMixing12);
      if(tweak_mixing12 <  fparm->MinValue(conv::kPhypMixing12)) tweak_mixing12 = fparm->MinValue(conv::kPhypMixing12);
      fparm->Update(conv::kPhypMixing12, tweak_mixing12);
      fparm->SetInitValue(conv::kPhypMixing12, tweak_mixing12);
    }
    else if(sys.find("osc_mixing13") == 0) {
      double tweak_mixing13 = input_mixing13 * (1 + conv::kGlobFitErrSinsq2Theta13);
      if(tweak_mixing13 >  fparm->MaxValue(conv::kPhypMixing13)) tweak_mixing13 = fparm->MaxValue(conv::kPhypMixing13);
      if(tweak_mixing13 <  fparm->MinValue(conv::kPhypMixing13)) tweak_mixing13 = fparm->MinValue(conv::kPhypMixing13);
      fparm->Update(conv::kPhypMixing13, tweak_mixing13);
      fparm->SetInitValue(conv::kPhypMixing13, tweak_mixing13);
    }
    else if(sys.find("osc_dm12sq") == 0) {
      double tweak_dm12sq = input_dm12sq * (1 + conv::kGlobFitErrDm12sq);
      if(tweak_dm12sq >  fparm->MaxValue(conv::kPhypDm12sq)) tweak_dm12sq = fparm->MaxValue(conv::kPhypDm12sq);
      if(tweak_dm12sq <  fparm->MinValue(conv::kPhypDm12sq)) tweak_dm12sq = fparm->MinValue(conv::kPhypDm12sq);
      fparm->Update(conv::kPhypDm12sq, tweak_dm12sq);
      fparm->SetInitValue(conv::kPhypDm12sq, tweak_dm12sq);
    }
    else if(sys.find("osc_deltacp") == 0) {
      double tweak_deltacp = input_deltacp + (TMath::Pi()/2);
      if(tweak_deltacp >  fparm->MaxValue(conv::kPhypDeltaCP)) tweak_deltacp = fparm->MaxValue(conv::kPhypDeltaCP);
      if(tweak_deltacp <  fparm->MinValue(conv::kPhypDeltaCP)) tweak_deltacp = fparm->MinValue(conv::kPhypDeltaCP);
      fparm->Update(conv::kPhypDeltaCP, tweak_deltacp);
      fparm->SetInitValue(conv::kPhypDeltaCP, tweak_deltacp);
    }
    else if(sys == "nominal") {
      //do nothing
    }
    else {
      skexp->TweakSystematic(systematic,+1.);
    }

    // Generate toy MC
    cout << "Generating toy-MC..." << endl;
    cout << "p.d.f parameters: " << endl;
    skexp->PrintSystematicVariations();
    skobs = new SuperkObservation;
    skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), fluctuate_stat);

    // Reset systematic variation
    skexp->ResetSystematicVariations();

    // Get fitter
    fit_inputs = new FitInputs(skobs,skexp,fparm);
    fitter = new FitterMinuit(fit_inputs,fit_option);

    // Do the fit
    fitter->DoTheFit();

    // Store best-fit params
    if(hyp == conv::kHypOsc3F) {
      bestfit_mixing12_p1sig = fparm->BestFitValue(conv::kPhypMixing12);
      bestfit_mixing13_p1sig = fparm->BestFitValue(conv::kPhypMixing13);
      bestfit_mixing23_p1sig = fparm->BestFitValue(conv::kPhypMixing23);
      bestfit_deltacp_p1sig  = fparm->BestFitValue(conv::kPhypDeltaCP);
      bestfit_dm12sq_p1sig   = fparm->BestFitValue(conv::kPhypDm12sq);
      bestfit_dm23sq_p1sig   = fparm->BestFitValue(conv::kPhypDm23sq);
    }
    else if(hyp == conv::kHypNuMuDisapp2F) {
      bestfit_mixing23_p1sig = fparm->BestFitValue(conv::kPhypMixing23);
      bestfit_dm23sq_p1sig   = fparm->BestFitValue(conv::kPhypDm23sq);
    }

    // Clean-up in preparation for the -1 sigma fit
    delete skobs;
    delete fit_inputs;
    delete fitter;
    skexp->ResetSystematicVariations();
    fparm->ResetCurrToInitValues();
    //need to do this, just incase InitValues were changed inside FitterKernel::DoTheFit()
    if(hyp == conv::kHypOsc3F) {
      fparm->SetInitValue(conv::kPhypMixing12, input_mixing12);
      fparm->SetInitValue(conv::kPhypMixing13, input_mixing13);
      fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
      fparm->SetInitValue(conv::kPhypDeltaCP,  input_deltacp);
      fparm->SetInitValue(conv::kPhypDm12sq,   input_dm12sq);
      fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
    }
    else if(hyp == conv::kHypNuMuDisapp2F) {
      fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
      fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
    }  
  
    //
    //
    // -1 sigma (relative to tuned BANFF)
    //
    //
    cout << "\n\nApplying -1 sigma tweak" << endl;
    if(sys.find("banff_") == 0) {
      //find parameter to tweak
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
      //initialise tweak vector
      const int nbanff = conv::kNumOfParamsInBANFFfit;
      TVectorD twkv(nbanff);
      const TMatrixD & cov = skexp->GetBANFF()->CovMatrix();
      //tweak the parameter to tweak
      twkv[ip] = skexp->GetBANFF()->BestFitValue(p) + (-1 * sqrt(cov[ip][ip]));
      //if we're tweaking the other parameters in a correlated way too
      double sigma = (twkv[ip] - skexp->GetBANFF()->BestFitValue(p)) / cov[ip][ip];
      for(int ipp = 0; ipp < nbanff; ipp++) {
	int position = (int)conv::kNuisBANFF00 + ipp;
	conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	if(pp == p)
	  continue;
	//set non-tweaked values to best fit
	twkv[ipp] = skexp->GetBANFF()->BestFitValue(pp);
	//now, tweak the other parameters according to the correlations between them and the fitted parameters
	if(correlated_parameters)
	  twkv[ipp] += sigma * cov[ip][ipp];
      }
      skexp->GetBANFF()->Tweak(true,&twkv);
    }
    else if(sys.find("skdetfsi_") == 0) {
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->SKDetFSISystIdx(p);
      //initialise tweak vector
      const int nskdetfsi = conv::kNumOfSKDetFSIParams;
      TVectorD twkv(nskdetfsi);
      const TMatrixD & cov = skexp->GetSuperkDetectorFSI()->CovMatrix();
      //tweak the parameter to tweak
      twkv[ip] = 1 (-1 * sqrt(cov[ip][ip]));
      //if we're tweaking the other parameters in a correlated way too
      double sigma = (twkv[ip] - 1) / cov[ip][ip];
      for(int ipp = 0; ipp < nskdetfsi; ipp++) {
	int position = (int)conv::kNuisSKDetFSI000 + ipp;
	conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	if(pp == p)
	  continue;
	//set non-tweaked values to best fit
	twkv[ipp] = 1;
	//now, tweak the other parameters according to the correlations between them and the fitted parameters
	if(correlated_parameters)
	  twkv[ipp] += sigma * cov[ip][ipp];
      }
      skexp->GetSuperkDetectorFSI()->Tweak(true,&twkv);
    }
    else if(sys.find("fsi_") == 0) {
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->FSISystIdx(p);
      //initialise tweak vector
      const int nfsi = conv::kNumOfFSIParams;
      TVectorD twkv(nfsi);
      const TMatrixD & cov = skexp->GetFSI()->CovMatrix();
      //tweak the parameter to tweak
      twkv[ip] = 1 + (-1 * sqrt(cov[ip][ip]));
      //if we're tweaking the other parameters in a correlated way too
      double sigma = (twkv[ip] - 1) / cov[ip][ip];
      for(int ipp = 0; ipp < nfsi; ipp++) {
	int position = (int)conv::kNuisFSI000 + ipp;
	conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	if(pp == p)
	  continue;
	//set non-tweaked values to best fit
	twkv[ipp] = 1;
	//now, tweak the other parameters according to the correlations between them and the fitted parameters
	if(correlated_parameters)
	  twkv[ipp] += sigma * cov[ip][ipp];
      }
      skexp->GetFSI()->Tweak(true,&twkv);
    }
    else if(sys.find("osc_mixing12") == 0) {
      double tweak_mixing12 = input_mixing12 * (1 - conv::kGlobFitErrSinsq2Theta12);
      if(tweak_mixing12 >  fparm->MaxValue(conv::kPhypMixing12)) tweak_mixing12 = fparm->MaxValue(conv::kPhypMixing12);
      if(tweak_mixing12 <  fparm->MinValue(conv::kPhypMixing12)) tweak_mixing12 = fparm->MinValue(conv::kPhypMixing12);
      fparm->Update(conv::kPhypMixing12, tweak_mixing12);
      fparm->SetInitValue(conv::kPhypMixing12, tweak_mixing12);
    }
    else if(sys.find("osc_mixing13") == 0) {
      double tweak_mixing13 = input_mixing13 * (1 - conv::kGlobFitErrSinsq2Theta13);
      if(tweak_mixing13 >  fparm->MaxValue(conv::kPhypMixing13)) tweak_mixing13 = fparm->MaxValue(conv::kPhypMixing13);
      if(tweak_mixing13 <  fparm->MinValue(conv::kPhypMixing13)) tweak_mixing13 = fparm->MinValue(conv::kPhypMixing13);
      fparm->Update(conv::kPhypMixing13, tweak_mixing13);
      fparm->SetInitValue(conv::kPhypMixing13, tweak_mixing13);
    }
    else if(sys.find("osc_dm12sq") == 0) {
      double tweak_dm12sq = input_dm12sq * (1 - conv::kGlobFitErrDm12sq);
      if(tweak_dm12sq >  fparm->MaxValue(conv::kPhypDm12sq)) tweak_dm12sq = fparm->MaxValue(conv::kPhypDm12sq);
      if(tweak_dm12sq <  fparm->MinValue(conv::kPhypDm12sq)) tweak_dm12sq = fparm->MinValue(conv::kPhypDm12sq);
      fparm->Update(conv::kPhypDm12sq, tweak_dm12sq);
      fparm->SetInitValue(conv::kPhypDm12sq, tweak_dm12sq);
    }
    else if(sys.find("osc_deltacp") == 0) {
      double tweak_deltacp = input_deltacp - (TMath::Pi()/2);
      if(tweak_deltacp >  fparm->MaxValue(conv::kPhypDeltaCP)) tweak_deltacp = fparm->MaxValue(conv::kPhypDeltaCP);
      if(tweak_deltacp <  fparm->MinValue(conv::kPhypDeltaCP)) tweak_deltacp = fparm->MinValue(conv::kPhypDeltaCP);
      fparm->Update(conv::kPhypDeltaCP, tweak_deltacp);
      fparm->SetInitValue(conv::kPhypDeltaCP, tweak_deltacp);
    }
    else if(sys == "nominal") {
      //do nothing
    }
    else {
      skexp->TweakSystematic(systematic,-1.);
    }

    // Generate toy MC
    cout << "Generating toy-MC..." << endl;
    cout << "p.d.f parameters: " << endl;
    skexp->PrintSystematicVariations();
    skobs = new SuperkObservation;
    skobs->CreateMockData(*skexp, *fparm, fluctuate_stat);

    // Reset systematic variation
    skexp->ResetSystematicVariations();

    // Get fitter
    fit_inputs = new FitInputs(skobs,skexp,fparm);
    fitter = new FitterMinuit(fit_inputs,fit_option);

    // Do the fit
    fitter->DoTheFit();

    // Store best-fit params
    if(hyp == conv::kHypOsc3F) {
      bestfit_mixing12_m1sig = fparm->BestFitValue(conv::kPhypMixing12);
      bestfit_mixing13_m1sig = fparm->BestFitValue(conv::kPhypMixing13);
      bestfit_mixing23_m1sig = fparm->BestFitValue(conv::kPhypMixing23);
      bestfit_deltacp_m1sig  = fparm->BestFitValue(conv::kPhypDeltaCP);
      bestfit_dm12sq_m1sig   = fparm->BestFitValue(conv::kPhypDm12sq);
      bestfit_dm23sq_m1sig   = fparm->BestFitValue(conv::kPhypDm23sq);
    }
    else if(hyp == conv::kHypNuMuDisapp2F) {
      bestfit_mixing23_m1sig = fparm->BestFitValue(conv::kPhypMixing23);
      bestfit_dm23sq_m1sig   = fparm->BestFitValue(conv::kPhypDm23sq);
    }
    if(sys == "xsec_sf") {
      //SF has no -1 sigma tweak!
      if(hyp == conv::kHypOsc3F) {
	bestfit_mixing12_m1sig = fparm->InitValue(conv::kPhypMixing12);
	bestfit_mixing13_m1sig = fparm->InitValue(conv::kPhypMixing13);
	bestfit_mixing23_m1sig = fparm->InitValue(conv::kPhypMixing23);
	bestfit_deltacp_m1sig  = fparm->InitValue(conv::kPhypDeltaCP);
	bestfit_dm12sq_m1sig   = fparm->InitValue(conv::kPhypDm12sq);
	bestfit_dm23sq_m1sig   = fparm->InitValue(conv::kPhypDm23sq);
      }
      else if(hyp == conv::kHypNuMuDisapp2F) {
	bestfit_mixing23_m1sig = fparm->InitValue(conv::kPhypMixing23);
	bestfit_dm23sq_m1sig   = fparm->InitValue(conv::kPhypDm23sq);
      }
    }

    // Clean-up in preparation for exiting
    delete skobs;
    delete fit_inputs;
    delete fitter;
  }//throw_n_toys <= 0
  else {

    //initialise
    double bestfit_mixing12_toys = 0;
    double bestfit_mixing13_toys = 0;
    double bestfit_mixing23_toys = 0;
    double bestfit_deltacp_toys  = 0;
    double bestfit_dm12sq_toys   = 0;
    double bestfit_dm23sq_toys   = 0;

    for(int itoy = 0; itoy < throw_n_toys; itoy++) {
      cout << "\n\Toy " << itoy << endl;

      // Reset systematic variation
      skexp->ResetSystematicVariations();

      //apply the tweak
      if(sys.find("banff_") == 0) {
	//find parameter to tweak
	conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
	int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
	//initialise tweak vector
	const int nbanff = conv::kNumOfParamsInBANFFfit;
	TVectorD twkv(nbanff);
	const TMatrixD & cov = skexp->GetBANFF()->CovMatrix();
	//tweak the parameter to tweak
	cout << *Rnd::Instance() << endl;
	int seed = -1;
	twkv[ip] = Rnd::Instance()->RndOther().Gaus(skexp->GetBANFF()->BestFitValue(p),sqrt(cov[ip][ip]));
	//if we're tweaking the other parameters in a correlated way too
	double sigma = (twkv[ip] - skexp->GetBANFF()->BestFitValue(p)) / cov[ip][ip];
	for(int ipp = 0; ipp < nbanff; ipp++) {
	  int position = (int)conv::kNuisBANFF00 + ipp;
	  conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	  if(pp == p)
	    continue;
	  //set non-tweaked values to best fit
	  twkv[ipp] = skexp->GetBANFF()->BestFitValue(pp);
	  //now, tweak the other parameters according to the correlations between them and the fitted parameters
	  if(correlated_parameters)
	    twkv[ipp] += sigma * cov[ip][ipp];
	}
	skexp->GetBANFF()->Tweak(true,&twkv);
      }
      else if(sys.find("skdetfsi_") == 0) {
	conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
	int ip = util::ArrIdx::Instance()->SKDetFSISystIdx(p);
	//initialise tweak vector
	const int nskdetfsi = conv::kNumOfSKDetFSIParams;
	TVectorD twkv(nskdetfsi);
	const TMatrixD & cov = skexp->GetSuperkDetectorFSI()->CovMatrix();
	//tweak the parameter to tweak
	twkv[ip] = Rnd::Instance()->RndOther().Gaus(1,sqrt(cov[ip][ip]));
	//if we're tweaking the other parameters in a correlated way too
	double sigma = (twkv[ip] - 1) / cov[ip][ip];
	for(int ipp = 0; ipp < nskdetfsi; ipp++) {
	  int position = (int)conv::kNuisSKDetFSI000 + ipp;
	  conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	  if(pp == p)
	    continue;
	  //set non-tweaked values to best fit
	  twkv[ipp] = 1;
	  //now, tweak the other parameters according to the correlations between them and the fitted parameters
	  if(correlated_parameters)
	    twkv[ipp] += sigma * cov[ip][ipp];
	}
	skexp->GetSuperkDetectorFSI()->Tweak(true,&twkv);
      }
      else if(sys.find("fsi_") == 0) {
	conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
	int ip = util::ArrIdx::Instance()->FSISystIdx(p);
	//initialise tweak vector
	const int nfsi = conv::kNumOfFSIParams;
	TVectorD twkv(nfsi);
	const TMatrixD & cov = skexp->GetFSI()->CovMatrix();
	//tweak the parameter to tweak
	twkv[ip] = Rnd::Instance()->RndOther().Gaus(1,sqrt(cov[ip][ip]));
	//if we're tweaking the other parameters in a correlated way too
	double sigma = (twkv[ip] - 1) / cov[ip][ip];
	for(int ipp = 0; ipp < nfsi; ipp++) {
	  int position = (int)conv::kNuisFSI000 + ipp;
	  conv::PhysParamId_t pp = (conv::PhysParamId_t)position;
	  if(pp == p)
	    continue;
	  //set non-tweaked values to best fit
	  twkv[ipp] = 1;
	  //now, tweak the other parameters according to the correlations between them and the fitted parameters
	  if(correlated_parameters)
	    twkv[ipp] += sigma * cov[ip][ipp];
	}
	skexp->GetFSI()->Tweak(true,&twkv);
      }
      else if(sys.find("osc_mixing12") == 0) {
	double tweak_mixing12 = util::RandomizeMixing12();
	if(tweak_mixing12 >  fparm->MaxValue(conv::kPhypMixing12)) tweak_mixing12 = fparm->MaxValue(conv::kPhypMixing12);
	if(tweak_mixing12 <  fparm->MinValue(conv::kPhypMixing12)) tweak_mixing12 = fparm->MinValue(conv::kPhypMixing12);
	fparm->Update(conv::kPhypMixing12, tweak_mixing12);
	fparm->SetInitValue(conv::kPhypMixing12, tweak_mixing12);
      }
      else if(sys.find("osc_mixing13") == 0) {
	double tweak_mixing13 = util::RandomizeMixing13();
	if(tweak_mixing13 >  fparm->MaxValue(conv::kPhypMixing13)) tweak_mixing13 = fparm->MaxValue(conv::kPhypMixing13);
	if(tweak_mixing13 <  fparm->MinValue(conv::kPhypMixing13)) tweak_mixing13 = fparm->MinValue(conv::kPhypMixing13);
	fparm->Update(conv::kPhypMixing13, tweak_mixing13);
	fparm->SetInitValue(conv::kPhypMixing13, tweak_mixing13);
      }
      else if(sys.find("osc_dm12sq") == 0) {
	double tweak_dm12sq = util::RandomizeDm12sq();
	if(tweak_dm12sq >  fparm->MaxValue(conv::kPhypDm12sq)) tweak_dm12sq = fparm->MaxValue(conv::kPhypDm12sq);
	if(tweak_dm12sq <  fparm->MinValue(conv::kPhypDm12sq)) tweak_dm12sq = fparm->MinValue(conv::kPhypDm12sq);
	fparm->Update(conv::kPhypDm12sq, tweak_dm12sq);
	fparm->SetInitValue(conv::kPhypDm12sq, tweak_dm12sq);
      }
      else if(sys.find("osc_deltacp") == 0) {
	double tweak_deltacp = util::RandomizeDeltaCP();
	if(tweak_deltacp >  fparm->MaxValue(conv::kPhypDeltaCP)) tweak_deltacp = fparm->MaxValue(conv::kPhypDeltaCP);
	if(tweak_deltacp <  fparm->MinValue(conv::kPhypDeltaCP)) tweak_deltacp = fparm->MinValue(conv::kPhypDeltaCP);
	fparm->Update(conv::kPhypDeltaCP, tweak_deltacp);
	fparm->SetInitValue(conv::kPhypDeltaCP, tweak_deltacp);
      }
      else if(sys.find("xsec_sf") == 0) {
	double tweak_sf = Rnd::Instance()->RndOther().Uniform(0.,1.);
	skexp->TweakSystematic(systematic, tweak_sf);
      }
      else if(sys == "nominal") {
	//do nothing
      }
      else if(sys == "mec") {
	double tweak_mec = Rnd::Instance()->RndOther().Gaus(-0.7,0.15);
	skexp->GetSuperkMC()->SetXMEC(tweak_mec);
      }
      else {
	double tweak = Rnd::Instance()->RndOther().Gaus(0.,1.);
	skexp->TweakSystematic(systematic, tweak);
      }

      // Generate toy MC
      cout << "Generating toy-MC..." << endl;
      cout << "p.d.f parameters: " << endl;
      skexp->PrintSystematicVariations();
      skobs = new SuperkObservation;
      skobs->CreateMockData(*skexp, *fparm, fluctuate_stat);

      // Reset systematic variation
      skexp->ResetSystematicVariations();
      skexp->GetSuperkMC()->SetXMEC(-1);

      // Get fitter
      fit_inputs = new FitInputs(skobs,skexp,fparm);
      fitter = new FitterMinuit(fit_inputs,fit_option);

      // Do the fit
      fitter->DoTheFit();

      // Store best-fit params
      if(hyp == conv::kHypOsc3F) {
	bestfit_mixing12_toys += fparm->BestFitValue(conv::kPhypMixing12);
	bestfit_mixing13_toys += fparm->BestFitValue(conv::kPhypMixing13);
	bestfit_mixing23_toys += fparm->BestFitValue(conv::kPhypMixing23);
	bestfit_deltacp_toys  += fparm->BestFitValue(conv::kPhypDeltaCP);
	bestfit_dm12sq_toys   += fparm->BestFitValue(conv::kPhypDm12sq);
	bestfit_dm23sq_toys   += fparm->BestFitValue(conv::kPhypDm23sq);
      }
      else if(hyp == conv::kHypNuMuDisapp2F) {
	bestfit_mixing23_toys += fparm->BestFitValue(conv::kPhypMixing23);
	bestfit_dm23sq_toys   += fparm->BestFitValue(conv::kPhypDm23sq);
      }

      //clean
      delete skobs;
      delete fit_inputs;
      delete fitter;

      //reset osc params to initial
      fparm->ResetCurrToInitValues();
      //need to do this, just incase InitValues were changed inside FitterKernel::DoTheFit()
      if(hyp == conv::kHypOsc3F) {
	fparm->SetInitValue(conv::kPhypMixing12, input_mixing12);
	fparm->SetInitValue(conv::kPhypMixing13, input_mixing13);
	fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
	fparm->SetInitValue(conv::kPhypDeltaCP,  input_deltacp);
	fparm->SetInitValue(conv::kPhypDm12sq,   input_dm12sq);
	fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
      }
      else if(hyp == conv::kHypNuMuDisapp2F) {
	fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
	fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
      }  
    }//end loop over itoy
    
    // Store best-fit params
    if(hyp == conv::kHypOsc3F) {
      bestfit_mixing12_p1sig = bestfit_mixing12_toys / throw_n_toys;
      bestfit_mixing13_p1sig = bestfit_mixing13_toys / throw_n_toys;
      bestfit_mixing23_p1sig = bestfit_mixing23_toys / throw_n_toys;
      bestfit_deltacp_p1sig  = bestfit_deltacp_toys  / throw_n_toys;
      bestfit_dm12sq_p1sig   = bestfit_dm12sq_toys   / throw_n_toys;
      bestfit_dm23sq_p1sig   = bestfit_dm23sq_toys   / throw_n_toys;
      bestfit_mixing12_m1sig = bestfit_mixing12_toys / throw_n_toys;
      bestfit_mixing13_m1sig = bestfit_mixing13_toys / throw_n_toys;
      bestfit_mixing23_m1sig = bestfit_mixing23_toys / throw_n_toys;
      bestfit_deltacp_m1sig  = bestfit_deltacp_toys  / throw_n_toys;
      bestfit_dm12sq_m1sig   = bestfit_dm12sq_toys   / throw_n_toys;
      bestfit_dm23sq_m1sig   = bestfit_dm23sq_toys   / throw_n_toys;
    }
    else if(hyp == conv::kHypNuMuDisapp2F) {
      bestfit_mixing23_p1sig = bestfit_mixing23_toys / throw_n_toys;
      bestfit_dm23sq_p1sig   = bestfit_dm23sq_toys   / throw_n_toys;
      bestfit_mixing23_m1sig = bestfit_mixing23_toys / throw_n_toys;
      bestfit_dm23sq_m1sig   = bestfit_dm23sq_toys   / throw_n_toys;
    }

  }//throw_n_toys > 0

  delete skexp;

  // Write-out results
  ostringstream outfilename;
  string cor = "";
  if(correlated_parameters)
    cor = "_correlated";
  outfilename << "star_plot_line-" << cor << systematic << ".data";
  ofstream outfile;
  outfile.open(outfilename.str().c_str());
  outfile << systematic 
	  << "\t" << bestfit_mixing12_p1sig 
	  << "\t" << bestfit_mixing13_p1sig 
	  << "\t" << bestfit_mixing23_p1sig 
	  << "\t" << bestfit_deltacp_p1sig  
	  << "\t" << bestfit_dm12sq_p1sig   
	  << "\t" << bestfit_dm23sq_p1sig   
	  << "\t" << bestfit_mixing12_m1sig 
	  << "\t" << bestfit_mixing13_m1sig 
	  << "\t" << bestfit_mixing23_m1sig 
	  << "\t" << bestfit_deltacp_m1sig  
	  << "\t" << bestfit_dm12sq_m1sig   
	  << "\t" << bestfit_dm23sq_m1sig   
	  << endl;
  outfile.close();

  cout << bestfit_mixing23_p1sig  << "\t" << bestfit_mixing23_m1sig << endl
       <<  bestfit_dm23sq_p1sig   << "\t" << bestfit_dm23sq_m1sig   << endl;
  cout << "Done!" << endl;
}
