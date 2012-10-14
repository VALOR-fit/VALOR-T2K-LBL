//
// ------------------------------------------------------------------------------------------------
// Calculate SuperK Ereco spectrum error envelope due to systematic parameter variations,
// for +3, +1, -1 sigma variations.
// Correlated parameter variations generated using diagonal of covariance matrix
//
// There are many run modes:
//
//   0 tweak all parameters (BANFF correlated)
//   1 tweak all BANFF (correlated)
//   2 tweak all BANFF flux (correlated)
//   3 tweak all BANFF xsec (correlated)
//   4 tweak all unconstrained xsec
//   5 tweak all xsec (BANFF (correlated) + unconstrained)
//   6 tweak all sk detector (without energy scale)
//   7 tweak all sk ccqe ring counting efficiency
//   8 tweak all sk detector
//   9 tweak all fsi
//  10 tweak all (except energy scale!) (BANFF correlated)
//  11 tweak all BANFF uncorrelated
//
// 100 tweak banff00
// 101 tweak banff01
// ...
// 100+kNumOfParamsInBANFFfit tweak last banff parameter
//
// 200 tweak xsec_pf
// 201 tweak xsec_sf
// 202 tweak xsec_wshape
// 203 tweak xsec_shape_ccoth
// 204 tweak xsec_pilessdelta
// 205 tweak xsec_enushape
// 206 tweak xsec_norm_ccqeE1
// 207 tweak xsec_norm_ccqeE2
// 208 tweak xsec_norm_cc1piE1
// 209 tweak xsec_norm_cccoh
// 210 tweak xsec_norm_nc1piC
// 211 tweak xsec_norm_nccoh
// 212 tweak xsec_norm_ncoth
// 213 tweak xsec_norm_nuetonumu
// 214 tweak xsec_norm_antinutonumu
// 215 tweak xsec_bindE
//
// 300 tweak sk_energy_scale
// 301 tweak sk_1Rmu_ccqe_other_eff
// 302 tweak sk_1Rmu_ccnonqe_eff
// 303 tweak sk_1Rmu_nc_eff
// 304 tweak sk_1Rmu_nuecc_eff
// 305 tweak sk_1Rmu_ccqe_ring_counting_eff_0
// 306 tweak sk_1Rmu_ccqe_ring_counting_eff_1
// 307 tweak sk_1Rmu_ccqe_ring_counting_eff_2
//
// 400 tweak fsi00
// 401 tweak fsi01
// ...
// 400+kNumOfFSIParams tweak last fsi parameter
//
// SPECIAL COMPARISON TWEAKS
// 500 ring counting 0 + other
// 501 ring counting 1 + other
// 502 ring counting 2 + other
// 503 all ring counting + other
// 504 tweak xsec_pilessdelta, put name it fsiabs
//
// 600 tweak unbanff00
// 601 tweak unbanff01
// ...
// 600+kNumOfParamsInBANFFfit tweak last uncorrelated banff parameter
//
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//
#define NTWEAKS 4 //2 or 4!
void test_sk_ereco_spectrum_twk_sys
(
  int          mode,                                                /* options: see above */
  bool         update_rootfile    = true,                           /* "update" ROOT output file? If false, will "recreate" */
  const char * hypothesis         = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12     = -99,                            /* */
  double       input_mixing13     = -99,                            /* */
  double       input_mixing23     = -99,                            /* */
  double       input_deltacp      = -99,                            /* */
  double       input_dm12sq       = -99,                            /* */
  double       input_dm23sq       = -99,                            /* */
  double       ereco_max_plot     = 7,                              /* maximum of the ereco axis (in GeV) to plot */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;
  int is = (int)s;

  //Check if the mode is valid
  string modetitle = "";
  switch(mode) {
    case 0:
      modetitle = "all"; break;
    case 1:
      modetitle = "BANFF"; break;
    case 2:
      modetitle = "BANFF-flux"; break;
    case 3:
      modetitle = "BANFF-xsec"; break;
    case 4:
      modetitle = "unconstrained-xsec"; break;
    case 5:
      modetitle = "xsec-all(BANFF+unconstrained)"; break;
    case 6:
      modetitle = "sk-detector(without-energy-scale)"; break;
    case 7:
      modetitle = "sk-ccqe-ring-counting-efficiency"; break;
    case 8:
      modetitle = "sk-detector"; break;
    case 9:
      modetitle = "fsi"; break;
    case 10:
      modetitle = "all(without-energy-scale)"; break;
    case 11:
      modetitle = "BANFF-uncorrelated"; break;

    case 200:
      modetitle = "xsec_pf"; break;
    case 201:
      modetitle = "xsec_sf"; break;
    case 202:
      modetitle = "xsec_wshape"; break;
    case 203:
      modetitle = "xsec_shape_ccoth"; break;
    case 204:
      modetitle = "xsec_pilessdelta"; break;
    case 205:
      modetitle = "xsec_enushape"; break;
    case 206:
      modetitle = "xsec_norm_ccqeE1"; break;
    case 207:
      modetitle = "xsec_norm_ccqeE2"; break;
    case 208:
      modetitle = "xsec_norm_cc1piE1"; break;
    case 209:
      modetitle = "xsec_norm_cccoh"; break;
    case 210:
      modetitle = "xsec_norm_nc1piC"; break;
    case 211:
      modetitle = "xsec_norm_nccoh"; break;
    case 212:
      modetitle = "xsec_norm_ncoth"; break;
    case 213:
      modetitle = "xsec_norm_nuetonumu"; break;
    case 214:
      modetitle = "xsec_norm_antinutonu"; break;
    case 215:
      modetitle = "xsec_bindE"; break;

    case 300:
      modetitle = "sk_energy_scale"; break;
    case 301:
      modetitle = "sk_1Rmu_ccqe_other_eff"; break;
    case 302:
      modetitle = "sk_1Rmu_ccnonqe_eff"; break;
    case 303:
      modetitle = "sk_1Rmu_nc_eff"; break;
    case 304:
      modetitle = "sk_1Rmu_nuecc_eff"; break;
    case 305:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_0"; break;
    case 306:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_1"; break;
    case 307:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_2"; break;

    case 500:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_0+other"; break;
    case 501:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_1+other"; break;
    case 502:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_2+other"; break;
    case 503:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_all+other"; break;
    case 504:
      modetitle = "xsec_fsiabs"; break;

    default:
      if(mode >= 100 && mode < 100+conv::kNumOfParamsInBANFFfit) {
	modetitle = Form("banff_%02d",mode-100); break;
      }
      else if(mode >= 400 && mode < 400+conv::kNumOfFSIParams) {
	modetitle = Form("fsi_%03d",mode-400); break;
      }
      if(mode >= 600 && mode < 600+conv::kNumOfParamsInBANFFfit) {
	modetitle = Form("unbanff_%02d",mode-100); break;
      }
      else {
	cout << "Unknown mode " << mode << endl; exit(1); break;
      }
  }//switch mode

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

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

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Physics hypothesis: \n" << ph << endl;
  cout << *Opt::Instance() << endl;

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);
  
  // Output histograms and graphs
#if NTWEAKS == 4
  const int ntweaks = 4;
  int * tweaks[ntweaks] = {+3, +1, -1, -3};
#elif NTWEAKS == 2
  const int ntweaks = 2;
  int * tweaks[ntweaks] = {+1, -1};
#else
  cerr << "Unknown NTWEAKS" << endl;
  exit(1);
#endif
  TH1D * hErecoNominal = 0;
  TH1D * hEreco[ntweaks];
  TH1D * hErecoRatio[ntweaks];
  TH1D * hErecoModeNominal[conv::knm];
  TH1D * hErecoMode[conv::knm][ntweaks];
  TH1D * hErecoModeRatio[conv::knm][ntweaks];

  // Open output file
  string rootmode = "update";
  if(!update_rootfile)
    rootmode = "recreate";
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_sys_%s.root",oscparams.c_str()),rootmode.c_str());

  // Get the nominal
  skexp.Generate(ph);
  hErecoNominal = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  hErecoNominal->SetDirectory(0);
  hErecoNominal->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoNominal->GetYaxis()->SetTitle("Number of events");
  hErecoNominal->SetLineWidth(2);
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) << endl;
    TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
    assert(h);
    hErecoModeNominal[im] = new TH1D(*h);
    hErecoModeNominal[im]->SetDirectory(0);
    hErecoModeNominal[im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hErecoModeNominal[im]->GetYaxis()->SetTitle("Number of events");
    hErecoModeNominal[im]->SetFillStyle(0);
  }
  if(!update_rootfile) {
    hErecoNominal->Write("nominal");
    for(int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
      if(m == conv::kSKmUndefined) continue;
      hErecoModeNominal[im]->Write(Form("nominal_%s",util::AsStringCompact(m).c_str()));
    }
  }

  // Tweak!
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

    if(mode < 100) {

      //all BANFF
      if(mode == 0 || mode == 1 || mode == 10) {
	TVectorD v(conv::kNumOfParamsInBANFFfit);
	v = skexp.GetBANFF()->BestFitValues();
	for(int ib = 0; ib < conv::kNumOfParamsInBANFFfit; ib++) {
	  v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	}
	skexp.GetBANFF()->Tweak(true,&v);
      }

      //BANFF flux only
      if(mode == 2) {
	TVectorD v(conv::kNumOfParamsInBANFFfit);
	v = skexp.GetBANFF()->BestFitValues();
	for(int ib = 0; ib < conv::kNumOfFluxParamsInBANFFfit; ib++) {
	  v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	}
	skexp.GetBANFF()->Tweak(true,&v);
      }

      //BANFF xsec only
      if(mode == 3 || mode == 5) {
	TVectorD v(conv::kNumOfParamsInBANFFfit);
	v = skexp.GetBANFF()->BestFitValues();
	for(int ib = conv::kNumOfFluxParamsInBANFFfit; ib < conv::kNumOfParamsInBANFFfit; ib++) {
	  v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	}
	skexp.GetBANFF()->Tweak(true,&v);
      }

      if(mode == 0 || mode == 4 || mode == 5 || mode == 10) {
	skexp.TweakSystematic("xsec_pf",tweak);
	skexp.TweakSystematic("xsec_sf",tweak);
	skexp.TweakSystematic("xsec_wshape",tweak);
	skexp.TweakSystematic("xsec_shape_ccoth",tweak);
	skexp.TweakSystematic("xsec_pilessdelta",tweak);
	skexp.TweakSystematic("xsec_enushape",tweak);
	/*
	  skexp.TweakSystematic("xsec_norm_ccqeE1",tweak);
	  skexp.TweakSystematic("xsec_norm_ccqeE2",tweak);
	  skexp.TweakSystematic("xsec_norm_cc1piE1",tweak);
	*/
	skexp.TweakSystematic("xsec_norm_cccoh",tweak);
	skexp.TweakSystematic("xsec_norm_nc1piC",tweak);
	skexp.TweakSystematic("xsec_norm_nccoh",tweak);
	skexp.TweakSystematic("xsec_norm_ncoth",tweak);
	skexp.TweakSystematic("xsec_norm_nuetonumu",tweak);
	skexp.TweakSystematic("xsec_norm_antinutonu",tweak);
	skexp.TweakSystematic("xsec_bindE",tweak);
      }

      if(mode == 0 || mode == 6 || mode == 8 || mode == 10) {
	if(mode == 0 || mode == 8)
	  skexp.TweakSystematic("sk_energy_scale",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_other_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccnonqe_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_nc_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_nuecc_eff",tweak);
      }

      if(mode == 0 || mode == 6 || mode == 7 || mode == 8 || mode == 10) {
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak);
      }

      if(mode == 0 || mode == 9 || mode == 10) {
	TVectorD v(conv::kNumOfFSIParams);
	v.Zero();
	for(int ib = 0; ib < conv::kNumOfFSIParams; ib++) {
	  v[ib] = 1 + (tweak * sqrt(skexp.GetFSI()->CovMatrix()[ib][ib]));
	}
	skexp.GetFSI()->Tweak(true,&v);
      }

      //all BANFF uncorrelated
      if(mode == 0 || mode == 1 || mode == 10) {
	TVectorD v(conv::kNumOfParamsInBANFFfit);
	v = skexp.GetBANFF()->CurrentUncorrelatedValues();
	for(int ib = 0; ib < conv::kNumOfParamsInBANFFfit; ib++) {
	  v[ib] = 0 + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	}
	skexp.GetBANFF()->Tweak(true,&v,false);
      }

    }//mode<100

    else {

      switch(mode) {
        case 200:
	  skexp.TweakSystematic("xsec_pf",tweak); break;
        case 201:
	  skexp.TweakSystematic("xsec_sf",tweak); break;
        case 202:
	  skexp.TweakSystematic("xsec_wshape",tweak); break;
        case 203:
	  skexp.TweakSystematic("xsec_shape_ccoth",tweak); break;
        case 204:
	  skexp.TweakSystematic("xsec_pilessdelta",tweak); break;
        case 205:
	  skexp.TweakSystematic("xsec_enushape",tweak); break;
        case 206:
	  skexp.TweakSystematic("xsec_norm_ccqeE1",tweak); break;
        case 207:
	  skexp.TweakSystematic("xsec_norm_ccqeE2",tweak); break;
        case 208:
	  skexp.TweakSystematic("xsec_norm_cc1piE1",tweak); break;
        case 209:
	  skexp.TweakSystematic("xsec_norm_cccoh",tweak); break;
        case 210:
	  skexp.TweakSystematic("xsec_norm_nc1piC",tweak); break;
        case 211:
	  skexp.TweakSystematic("xsec_norm_nccoh",tweak); break;
        case 212:
	  skexp.TweakSystematic("xsec_norm_ncoth",tweak); break;
        case 213:
	  skexp.TweakSystematic("xsec_norm_nuetonumu",tweak); break;
        case 214:
	  skexp.TweakSystematic("xsec_norm_antinutonu",tweak); break;
        case 215:
	  skexp.TweakSystematic("xsec_bindE",tweak); break;

        case 300:
	  skexp.TweakSystematic("sk_energy_scale",tweak); break;
        case 301:
	  skexp.TweakSystematic("sk_1Rmu_ccqe_other_eff",tweak); break;
        case 302:
	  skexp.TweakSystematic("sk_1Rmu_ccnonqe_eff",tweak); break;
        case 303:
	  skexp.TweakSystematic("sk_1Rmu_nc_eff",tweak); break;
        case 304:
	  skexp.TweakSystematic("sk_1Rmu_nuecc_eff",tweak); break;
        case 305:
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak); break;
        case 306:
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak); break;
        case 307:
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak); break;

        case 500: 
	  skexp.TweakSystematic("sk_energy_scale",tweak);
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak); break;
        case 501: 
	  skexp.TweakSystematic("sk_energy_scale",tweak);
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak); break;
        case 502: 
	  skexp.TweakSystematic("sk_energy_scale",tweak);
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak); break;
        case 503: 
	  skexp.TweakSystematic("sk_energy_scale",tweak);
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak);
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak);
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak); break;
        case 504:
	  skexp.TweakSystematic("xsec_pilessdelta",tweak); break;

        default:
	  if(mode >= 100 && mode < 100+conv::kNumOfParamsInBANFFfit) {
	    TVectorD v(conv::kNumOfParamsInBANFFfit);
	    v = skexp.GetBANFF()->BestFitValues();
	    int ib = mode - 100;
	    v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	    skexp.GetBANFF()->Tweak(true,&v);
	    break;
	  }
	  else if(mode >= 400 && mode < 400+conv::kNumOfFSIParams) {
	    TVectorD v(conv::kNumOfFSIParams);
	    for (int ipf = 0; ipf < conv::kNumOfFSIParams; ipf++) {
	      v[ipf] = 1;
	    }
	    int ib = mode - 400;
	    v[ib] = 1 + (tweak * sqrt(skexp.GetFSI()->CovMatrix()[ib][ib]));
	    skexp.GetFSI()->Tweak(true,&v);
	    break;
	  }
	  if(mode >= 600 && mode < 600+conv::kNumOfParamsInBANFFfit) {
	    TVectorD v(conv::kNumOfParamsInBANFFfit);
	    v = skexp.GetBANFF()->CurrentUncorrelatedValues();
	    int ib = mode - 600;
	    v[ib] = 0 + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	    skexp.GetBANFF()->Tweak(true,&v,false);
	    break;
	  }
	  else {
	    cout << "Unknown mode " << mode << endl; exit(1); break;
	  }
      }//switch(mode)
    }//mode>=100

    // generate predictions
    skexp.Generate(ph);
    //store spectrum
    hEreco[itweak] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
    hEreco[itweak]->SetDirectory(0);
    hEreco[itweak]->SetTitle(Form("%s_%d_sigma",modetitle.c_str(),tweak));
    hEreco[itweak]->SetLineColor(histcolor);
    hEreco[itweak]->Write(Form("%s_%d_sigma",modetitle.c_str(),tweak));		    
    hErecoRatio[itweak] = new TH1D(*hEreco[itweak]);
    hErecoRatio[itweak]->Divide(hErecoNominal);
    hErecoRatio[itweak]->Write(Form("%s_%d_sigma_ratio",modetitle.c_str(),tweak));
    for(int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
      if(m == conv::kSKmUndefined) continue;
      hErecoMode[im][itweak] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s,m));
      hErecoMode[im][itweak]->SetDirectory(0);
      hErecoMode[im][itweak]->SetTitle(Form("%s_%d_sigma_%s",modetitle.c_str(),tweak,util::AsStringCompact(m).c_str()));
      hErecoMode[im][itweak]->SetLineColor(histcolor);
      hErecoMode[im][itweak]->SetFillStyle(0);
      hErecoMode[im][itweak]->Write(Form("%s_%d_sigma_%s",modetitle.c_str(),tweak,util::AsStringCompact(m).c_str()));
      hErecoModeRatio[im][itweak] = new TH1D(*hErecoMode[im][itweak]);
      hErecoModeRatio[im][itweak]->Divide(hErecoModeNominal[im]);
      hErecoModeRatio[im][itweak]->Write(Form("%s_%d_sigma_ratio_%s",modetitle.c_str(),tweak,util::AsStringCompact(m).c_str()));
    }
  }//end loop over tweaks

  //
  // Draw
  //
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->Divide(1,2);
  c->cd(1);
  TLegend * legend = new TLegend(0.7, 0.4, 0.9, 0.85);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->AddEntry(hErecoNominal, "nominal", "L" );
  hErecoNominal->SetTitle(Form("%s",modetitle.c_str()));
  hErecoNominal->GetXaxis()->SetRangeUser(0,ereco_max_plot-conv::kOneMeV);
  hErecoNominal->Draw();
  for(int itweak = 0; itweak < ntweaks; itweak++) {
    hEreco[itweak]->Draw("same");
    legend->AddEntry(hEreco[itweak], Form("%d sigma",tweaks[itweak]), "L" );
  }
  legend->Draw();
  c->cd(2);
  for(int itweak = 0; itweak < ntweaks; itweak++) {
    if(itweak==0){
      hErecoRatio[itweak]->GetXaxis()->SetRangeUser(0,ereco_max_plot-conv::kOneMeV);
      hErecoRatio[itweak]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hErecoRatio[itweak]->GetYaxis()->SetTitle("Relative error");
      hErecoRatio[itweak]->SetTitle("Ratio");
      hErecoRatio[itweak]->Draw();
      continue;
    }
    hErecoRatio[itweak]->Draw("same");
  }
  c->SaveAs(Form("ereco_twk_%s.pdf",modetitle.c_str()));
  hErecoNominal->GetXaxis()->UnZoom();
  hErecoRatio[0]->GetXaxis()->UnZoom();
  c->SaveAs(Form("ereco_twk_%s_fullrange.pdf",modetitle.c_str()));
  // Draw mode by mode
  TCanvas * cm = new TCanvas("cm","",20,20,500,500);
  gStyle->SetOptStat(0);
  cm->Divide(1,2);
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cm->cd(1);
    TLegend * legend = new TLegend(0.7, 0.4, 0.9, 0.85);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->AddEntry(hErecoModeNominal[im], "nominal", "L" );
    hErecoModeNominal[im]->SetTitle(Form("%s_%s",modetitle.c_str(),util::AsStringCompact(m).c_str()));
    hErecoModeNominal[im]->GetXaxis()->SetRangeUser(0,ereco_max_plot-conv::kOneMeV);
    hErecoModeNominal[im]->Draw();
    for(int itweak = 0; itweak < ntweaks; itweak++) {
      hErecoMode[im][itweak]->Draw("same");
      legend->AddEntry(hErecoMode[im][itweak], Form("%d sigma",tweaks[itweak]), "L" );
    }
    legend->Draw();
    cm->cd(2);
    for(int itweak = 0; itweak < ntweaks; itweak++) {
      if(itweak==0){
	hErecoModeRatio[im][itweak]->GetXaxis()->SetRangeUser(0,ereco_max_plot-conv::kOneMeV);
	hErecoModeRatio[im][itweak]->GetXaxis()->SetTitle("E_{reco} (GeV)");
	hErecoModeRatio[im][itweak]->GetYaxis()->SetTitle("Relative error");
	hErecoModeRatio[im][itweak]->SetTitle("Ratio");
	hErecoModeRatio[im][itweak]->Draw();
	continue;
      }
      hErecoModeRatio[im][itweak]->Draw("same");
    }
    cm->SaveAs(Form("ereco_twk_%s_%s.pdf",modetitle.c_str(),util::AsStringCompact(m).c_str()));
    hErecoModeNominal[im]->GetXaxis()->UnZoom();
    hErecoModeRatio[im][0]->GetXaxis()->UnZoom();
    cm->SaveAs(Form("ereco_twk_%s_%s_fullrange.pdf",modetitle.c_str(),util::AsStringCompact(m).c_str()));
  }//loop over modes

  //
  // Save the change in nsk
  //
  double nsk_0 = hErecoNominal->Integral(1,hErecoNominal->GetNbinsX());
  double nsk_p  = hEreco[ntweaks/2 - 1]->Integral(1,hEreco[ntweaks/2 - 1]->GetNbinsX());
  double nsk_m  = hEreco[ntweaks/2    ]->Integral(1,hEreco[ntweaks/2    ]->GetNbinsX());
  double frac_nsk_p = (nsk_0 > 0) ? 100.*(nsk_p - nsk_0) / nsk_0 : 0.;
  double frac_nsk_m = (nsk_0 > 0) ? 100.*(nsk_m - nsk_0) / nsk_0 : 0.;
  double nskmode_p[conv::knm] = {0}, nskmode_m[conv::knm] = {0};
  double frac_nskmode_p[conv::knm] = {0}, frac_nskmode_m[conv::knm] = {0};
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    nskmode_p[im] = hErecoMode[im][ntweaks/2 - 1]->Integral(1,hErecoMode[im][ntweaks/2 - 1]->GetNbinsX());
    nskmode_m[im] = hErecoMode[im][ntweaks/2    ]->Integral(1,hErecoMode[im][ntweaks/2    ]->GetNbinsX());
    frac_nskmode_p[im] = (nsk_0 > 0) ? 100.*(nskmode_p[im] - nsk_0) / nsk_0 : 0.;
    frac_nskmode_m[im] = (nsk_0 > 0) ? 100.*(nskmode_m[im] - nsk_0) / nsk_0 : 0.;
  }

  cout << endl << "Nominal number of events: " << nsk_0 << endl
       << "   +1 sigma:              " << nsk_p << "\t" << frac_nsk_p << "%" << endl
       << "   -1 sigma:              " << nsk_m << "\t" << frac_nsk_m << "%" << endl;
#line 0
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << util::AsString(m) << "\t" 
	 << nskmode_p[im] << "   " << frac_nskmode_p[im] << "%" << "\t"
	 << nskmode_m[im] << "   " << frac_nskmode_m[im] << "%" << endl;
  }

  TH1D * summary_p = 0, * summary_m = 0, * categories_p = 0, * categories_m = 0;
  TH1D * summarymode_p[conv::knm], * summarymode_m[conv::knm], * categoriesmode_p[conv::knm], * categoriesmode_m[conv::knm];
  f.GetObject("summary_p",summary_p);
  f.GetObject("summary_m",summary_m);
  f.GetObject("categories_p",categories_p);
  f.GetObject("categories_m",categories_m);
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    summarymode_p[im] = 0;
    summarymode_m[im] = 0;
    categoriesmode_p[im] = 0;
    categoriesmode_m[im] = 0;
    f.GetObject(Form("summary_p_%s",util::AsStringCompact(m).c_str()),summarymode_p[im]);
    f.GetObject(Form("summary_m_%s",util::AsStringCompact(m).c_str()),summarymode_m[im]);
    f.GetObject(Form("categories_p_%s",util::AsStringCompact(m).c_str()),categoriesmode_p[im]);
    f.GetObject(Form("categories_m_%s",util::AsStringCompact(m).c_str()),categoriesmode_m[im]);
  }

  if(!summary_p) {
    int nbins = (2*conv::kNumOfParamsInBANFFfit) + 19 + 8 + conv::kNumOfFSIParams;
    summary_p = new TH1D("summary_p","summary_p",nbins,0,nbins);
    summary_p->SetLineColor(kRed);
    summary_p->GetYaxis()->SetTitle("% change");
    for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
      summary_p->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+1,"xsec_pf");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+2,"xsec_sf");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+3,"xsec_wshape");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+4,"xsec_shape_ccoth");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+5,"xsec_pilessdelta");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+6,"xsec_enushape");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+7,"xsec_norm_ccqeE1");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+8,"xsec_norm_ccqeE2");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+9,"xsec_norm_cc1piE1");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+10,"xsec_norm_cccoh");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+11,"xsec_norm_nc1piC");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+12,"xsec_norm_nccoh");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+13,"xsec_norm_ncoth");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+14,"xsec_norm_nuetonumu");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+15,"xsec_norm_antinutonu");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+16,"sk_energy_scale");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+17,"sk_1Rmu_ccqe_other_eff");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+18,"sk_1Rmu_ccnonqe_eff");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+19,"sk_1Rmu_nc_eff");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+20,"sk_1Rmu_nuecc_eff");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+21,"sk_1Rmu_ccqe_ring_counting_eff_0");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+22,"sk_1Rmu_ccqe_ring_counting_eff_1");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+23,"sk_1Rmu_ccqe_ring_counting_eff_2");
    summary_p->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+24,"xsec_bindE");
    for(int i = 0; i < conv::kNumOfFSIParams; i++)
      summary_p->GetXaxis()->SetBinLabel(i+conv::kNumOfParamsInBANFFfit+23,Form("fsi_%03d",i));
    for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
      summary_p->GetXaxis()->SetBinLabel(i+conv::kNumOfParamsInBANFFfit+23+conv::kNumOfFSIParams,Form("unbanff_%02d",i));
  }
  if(!summary_m) {
    int nbins =  (2*conv::kNumOfParamsInBANFFfit) + 19 + 8 + conv::kNumOfFSIParams;
    summary_m = new TH1D("summary_m","summary_m",nbins,0,nbins);
    summary_m->SetLineColor(kBlue);
    summary_m->GetYaxis()->SetTitle("% change");
    for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
      summary_m->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+1,"xsec_pf");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+2,"xsec_sf");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+3,"xsec_wshape");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+4,"xsec_shape_ccoth");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+5,"xsec_pilessdelta");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+6,"xsec_enushape");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+7,"xsec_norm_ccqeE1");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+8,"xsec_norm_ccqeE2");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+9,"xsec_norm_cc1piE1");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+10,"xsec_norm_cccoh");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+11,"xsec_norm_nc1piC");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+12,"xsec_norm_nccoh");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+13,"xsec_norm_ncoth");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+14,"xsec_norm_nuetonumu");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+15,"xsec_norm_antinutonu");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+16,"sk_energy_scale");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+17,"sk_1Rmu_ccqe_other_eff");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+18,"sk_1Rmu_ccnonqe_eff");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+19,"sk_1Rmu_nc_eff");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+20,"sk_1Rmu_nuecc_eff");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+21,"sk_1Rmu_ccqe_ring_counting_eff_0");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+22,"sk_1Rmu_ccqe_ring_counting_eff_1");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+23,"sk_1Rmu_ccqe_ring_counting_eff_2");
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+24,"xsec_bindE");
    for(int i = 0; i < conv::kNumOfFSIParams; i++)
      summary_m->GetXaxis()->SetBinLabel(i+conv::kNumOfParamsInBANFFfit+23,Form("fsi_%03d",i));
    for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
      summary_m->GetXaxis()->SetBinLabel(i+conv::kNumOfParamsInBANFFfit+23+conv::kNumOfFSIParams,Form("unbanff_%02d",i));
  }
  if(!categories_p) {
    categories_p = new TH1D("categories_p","categories_p",12,0,12);
    categories_p->SetLineColor(kRed);
    categories_p->GetYaxis()->SetTitle("% change");
    categories_p->GetXaxis()->SetBinLabel(1,"all");
    categories_p->GetXaxis()->SetBinLabel(2,"all(without-energy-scale)");
    categories_p->GetXaxis()->SetBinLabel(3,"BANFF");
    categories_p->GetXaxis()->SetBinLabel(4,"BANFF-flux");
    categories_p->GetXaxis()->SetBinLabel(5,"BANFF-xsec");
    categories_p->GetXaxis()->SetBinLabel(6,"unconstrained-xsec");
    categories_p->GetXaxis()->SetBinLabel(7,"xsec-all(BANFF+unconstrained)");
    categories_p->GetXaxis()->SetBinLabel(8,"sk-detector");
    categories_p->GetXaxis()->SetBinLabel(9,"sk-detector(without-energy-scale)");
    categories_p->GetXaxis()->SetBinLabel(10,"sk-ccqe-ring-counting-efficiency");
    categories_p->GetXaxis()->SetBinLabel(11,"fsi");
    categories_p->GetXaxis()->SetBinLabel(12,"BANFF-uncorrelated");
  }
  if(!categories_m) {
    categories_m = new TH1D("categories_m","categories_m",12,0,12);
    categories_m->SetLineColor(kRed);
    categories_m->GetYaxis()->SetTitle("% change");
    categories_m->GetXaxis()->SetBinLabel(1,"all");
    categories_m->GetXaxis()->SetBinLabel(2,"all(without-energy-scale)");
    categories_m->GetXaxis()->SetBinLabel(3,"BANFF");
    categories_m->GetXaxis()->SetBinLabel(4,"BANFF-flux");
    categories_m->GetXaxis()->SetBinLabel(5,"BANFF-xsec");
    categories_m->GetXaxis()->SetBinLabel(6,"unconstrained-xsec");
    categories_m->GetXaxis()->SetBinLabel(7,"xsec-all(BANFF+unconstrained)");
    categories_m->GetXaxis()->SetBinLabel(8,"sk-detector");
    categories_m->GetXaxis()->SetBinLabel(9,"sk-detector(without-energy-scale)");
    categories_m->GetXaxis()->SetBinLabel(10,"sk-ccqe-ring-counting-efficiency");
    categories_m->GetXaxis()->SetBinLabel(11,"fsi");
    categories_m->GetXaxis()->SetBinLabel(12,"BANFF-uncorrelated");
  }
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    if(!summarymode_p[im]) {
      int nbins = (2*conv::kNumOfParamsInBANFFfit) + 19 + 8 + conv::kNumOfFSIParams;
      summarymode_p[im] = new TH1D(Form("summarymode_p_%s",util::AsStringCompact(m).c_str()),Form("summarymode_p_%s",util::AsStringCompact(m).c_str()),nbins,0,nbins);
      summarymode_p[im]->SetLineColor(kRed);
      summarymode_p[im]->GetYaxis()->SetTitle("% change");
      for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
	summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+1,"xsec_pf");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+2,"xsec_sf");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+3,"xsec_wshape");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+4,"xsec_shape_ccoth");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+5,"xsec_pilessdelta");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+6,"xsec_enushape");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+7,"xsec_norm_ccqeE1");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+8,"xsec_norm_ccqeE2");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+9,"xsec_norm_cc1piE1");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+10,"xsec_norm_cccoh");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+11,"xsec_norm_nc1piC");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+12,"xsec_norm_nccoh");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+13,"xsec_norm_ncoth");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+14,"xsec_norm_nuetonumu");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+15,"xsec_norm_antinutonu");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+16,"sk_energy_scale");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+17,"sk_1Rmu_ccqe_other_eff");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+18,"sk_1Rmu_ccnonqe_eff");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+19,"sk_1Rmu_nc_eff");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+20,"sk_1Rmu_nuecc_eff");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+21,"sk_1Rmu_ccqe_ring_counting_eff_0");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+22,"sk_1Rmu_ccqe_ring_counting_eff_1");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+23,"sk_1Rmu_ccqe_ring_counting_eff_2");
      summarymode_p[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+24,"xsec_bindE");
      for(int i = 0; i < conv::kNumOfFSIParams; i++)
	summarymode_p[im]->GetXaxis()->SetBinLabel(i+conv::kNumOfParamsInBANFFfit+23,Form("fsi_%03d",i));
      for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
	summarymode_p[im]->GetXaxis()->SetBinLabel(i+conv::kNumOfParamsInBANFFfit+23+conv::kNumOfFSIParams,Form("unbanff_%02d",i));
    }
    if(!summarymode_m[im]) {
      int nbins = (2*conv::kNumOfParamsInBANFFfit) + 19 + 8 + conv::kNumOfFSIParams;
      summarymode_m[im] = new TH1D(Form("summarymode_m_%s",util::AsStringCompact(m).c_str()),Form("summarymode_m_%s",util::AsStringCompact(m).c_str()),nbins,0,nbins);
      summarymode_m[im]->SetLineColor(kBlue);
      summarymode_m[im]->GetYaxis()->SetTitle("% change");
      for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
	summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+1,"xsec_pf");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+2,"xsec_sf");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+3,"xsec_wshape");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+4,"xsec_shape_ccoth");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+5,"xsec_pilessdelta");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+6,"xsec_enushape");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+7,"xsec_norm_ccqeE1");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+8,"xsec_norm_ccqeE2");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+9,"xsec_norm_cc1piE1");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+10,"xsec_norm_cccoh");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+11,"xsec_norm_nc1piC");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+12,"xsec_norm_nccoh");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+13,"xsec_norm_ncoth");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+14,"xsec_norm_nuetonumu");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+15,"xsec_norm_antinutonu");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+16,"sk_energy_scale");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+17,"sk_1Rmu_ccqe_other_eff");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+18,"sk_1Rmu_ccnonqe_eff");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+19,"sk_1Rmu_nc_eff");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+20,"sk_1Rmu_nuecc_eff");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+21,"sk_1Rmu_ccqe_ring_counting_eff_0");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+22,"sk_1Rmu_ccqe_ring_counting_eff_1");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+23,"sk_1Rmu_ccqe_ring_counting_eff_2");
      summarymode_m[im]->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+24,"xsec_bindE");
      for(int i = 0; i < conv::kNumOfFSIParams; i++)
	summarymode_m[im]->GetXaxis()->SetBinLabel(i+conv::kNumOfParamsInBANFFfit+23,Form("fsi_%03d",i));
      for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
	summarymode_m[im]->GetXaxis()->SetBinLabel(i+conv::kNumOfParamsInBANFFfit+23+conv::kNumOfFSIParams,Form("unbanff_%02d",i));
    }
    if(!categoriesmode_p[im]) {
      categoriesmode_p[im] = new TH1D(Form("categoriesmode_p_%s",util::AsStringCompact(m).c_str()),Form("categoriesmode_p_%s",util::AsStringCompact(m).c_str()),12,0,12);
      categoriesmode_p[im]->SetLineColor(kRed);
      categoriesmode_p[im]->GetYaxis()->SetTitle("% change");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(1,"all");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(2,"all(without-energy-scale)");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(3,"BANFF");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(4,"BANFF-flux");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(5,"BANFF-xsec");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(6,"unconstrained-xsec");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(7,"xsec-all(BANFF+unconstrained)");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(8,"sk-detector");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(9,"sk-detector(without-energy-scale)");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(10,"sk-ccqe-ring-counting-efficiency");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(11,"fsi");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(12,"BANFF-uncorrelated");
    }
    if(!categoriesmode_m[im]) {
      categoriesmode_m[im] = new TH1D(Form("categoriesmode_m_%s",util::AsStringCompact(m).c_str()),Form("categoriesmode_m_%s",util::AsStringCompact(m).c_str()),12,0,12);
      categoriesmode_m[im]->SetLineColor(kRed);
      categoriesmode_m[im]->GetYaxis()->SetTitle("% change");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(1,"all");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(2,"all(without-energy-scale)");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(3,"BANFF");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(4,"BANFF-flux");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(5,"BANFF-xsec");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(6,"unconstrained-xsec");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(7,"xsec-all(BANFF+unconstrained)");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(8,"sk-detector");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(9,"sk-detector(without-energy-scale)");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(10,"sk-ccqe-ring-counting-efficiency");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(11,"fsi");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(12,"BANFF-uncorrelated");
    }
  }


  if(mode >= 100) {
    summary_p->Fill(modetitle.c_str(),frac_nsk_p);
    summary_m->Fill(modetitle.c_str(),frac_nsk_m);
    summary_p->Write();
    summary_m->Write();
    for(int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
      if(m == conv::kSKmUndefined) continue;
      summarymode_p[im]->Fill(modetitle.c_str(),frac_nskmode_p[im]);
      summarymode_m[im]->Fill(modetitle.c_str(),frac_nskmode_m[im]);
      summarymode_p[im]->Write();
      summarymode_m[im]->Write();

    }
  }
  else{
    categories_p->Fill(modetitle.c_str(),frac_nsk_p);
    categories_m->Fill(modetitle.c_str(),frac_nsk_m);
    categories_p->Write();
    categories_m->Write();
    for(int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
      if(m == conv::kSKmUndefined) continue;
      categoriesmode_p[im]->Fill(modetitle.c_str(),frac_nskmode_p[im]);
      categoriesmode_m[im]->Fill(modetitle.c_str(),frac_nskmode_m[im]);
      categoriesmode_p[im]->Write();
      categoriesmode_m[im]->Write();

    }
  }
  f.Close();
}
  
