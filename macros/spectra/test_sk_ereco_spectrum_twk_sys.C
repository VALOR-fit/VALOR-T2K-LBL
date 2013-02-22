//
// ------------------------------------------------------------------------------------------------
// Calculate SuperK Ereco spectrum error envelope due to systematic parameter variations,
// for +3, +1, -1 sigma variations.
// Correlated parameter variations generated using diagonal of covariance matrix
//
// There are many run modes:
//
//   0 tweak all parameters
//   1 tweak all BANFF
//   2 tweak all BANFF flux
//   3 tweak all BANFF xsec
//   4 tweak all unconstrained xsec
//   5 tweak all xsec (BANFF + unconstrained)
//   6 tweak all sk detector (without energy scale)
//   7 tweak all sk ccqe ring counting efficiency
//   8 tweak all sk detector
//   9 tweak all fsi
//  10 tweak all (except energy scale!)
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
// 209 tweak xsec_norm_cccoh
// 210 tweak xsec_norm_nc1piC
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
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//
#define NTWEAKS 4 //2 or 4 or 6!
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
  const char * config_file        = "analysis_configuration.xml",   /* */
  bool         save_all_canvases  = false                           /* if false, will not save the individual mode plots, or the 'fullrange' plots */
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
      else {
	cout << "Unknown mode " << mode << endl; exit(1); break;
      }
  }//switch mode
  //name for the output plots
  conv::PhysParamId_t p = util::PhysParamIdFromString(modetitle);
  string modetitleplot = modetitle;
  if(p != conv::kPhypUndefined)
    modetitleplot = util::AsStringCompact(p);

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
  int tweaks[ntweaks] = {+3, +1, -1, -3};
#elif NTWEAKS == 2
  const int ntweaks = 2;
  int tweaks[ntweaks] = {+1, -1};
#elif NTWEAKS == 6
  const int ntweaks = 6;
  double tweaks[ntweaks] = {+1., +0.4545, +0.2273, -0.2273, -0.4545, -1.};
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
  hErecoNominal->GetYaxis()->SetTitle("Events per bin");
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
    hErecoModeNominal[im]->GetYaxis()->SetTitle("Events per bin");
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
#if NTWEAKS == 6
    double tweak = tweaks[itweak];
#else
    int tweak = tweaks[itweak];
#endif

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
      case 4:
	histcolor = kYellow-4;
	break;
      case 5:
	histcolor = kMagenta+1;
	break;
      default:
	histcolor = kBlack;
    }

    if(mode < 100) {

      //all BANFF correlated
      if(mode == 0 || mode == 1 || mode == 10) {
	TVectorD v(conv::kNumOfParamsInBANFFfit);
	v = skexp.GetBANFF()->BestFitValues();
	for(int ib = 0; ib < conv::kNumOfParamsInBANFFfit; ib++) {
	  v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	}
	skexp.GetBANFF()->Tweak(true,&v);
      }

      //BANFF flux only correlated
      if(mode == 2) {
	TVectorD v(conv::kNumOfParamsInBANFFfit);
	v = skexp.GetBANFF()->BestFitValues();
	for(int ib = 0; ib < conv::kNumOfFluxParamsInBANFFfit; ib++) {
	  v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	}
	skexp.GetBANFF()->Tweak(true,&v);
      }

      //BANFF xsec only correlated
      if(mode == 3 || mode == 5) {
	TVectorD v(conv::kNumOfParamsInBANFFfit);
	v = skexp.GetBANFF()->BestFitValues();
	for(int ib = conv::kNumOfFluxParamsInBANFFfit; ib < conv::kNumOfParamsInBANFFfit; ib++) {
	  v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	}
	skexp.GetBANFF()->Tweak(true,&v);
      }

      //Uncorrelated xsec
      if(mode == 0 || mode == 4 || mode == 5 || mode == 10) {
	skexp.TweakSystematic("xsec_pf",tweak);
	skexp.TweakSystematic("xsec_sf",tweak);
	skexp.TweakSystematic("xsec_wshape",tweak);
	skexp.TweakSystematic("xsec_shape_ccoth",tweak);
	skexp.TweakSystematic("xsec_pilessdelta",tweak);
	skexp.TweakSystematic("xsec_bindE",tweak);
	skexp.TweakSystematic("xsec_norm_cccoh",tweak);
	skexp.TweakSystematic("xsec_norm_nc1piC",tweak);
	skexp.TweakSystematic("xsec_norm_ncoth",tweak);
	skexp.TweakSystematic("xsec_norm_nuetonumu",tweak);
	skexp.TweakSystematic("xsec_norm_antinutonu",tweak);
      }

      //SK det uncorrelated
      if(mode == 0 || mode == 6 || mode == 8 || mode == 10) {
	if(mode == 0 || mode == 8)
	  skexp.TweakSystematic("sk_energy_scale",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_other_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccnonqe_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_nc_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_nuecc_eff",tweak);
      }

      //SK det ring counting
      if(mode == 0 || mode == 6 || mode == 7 || mode == 8 || mode == 10) {
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak);
      }

      //FSI
      if(mode == 0 || mode == 9 || mode == 10) {
	TVectorD v(conv::kNumOfFSIParams);
	v.Zero();
	for(int ib = 0; ib < conv::kNumOfFSIParams; ib++) {
	  v[ib] = 1 + (tweak * sqrt(skexp.GetFSI()->CovMatrix()[ib][ib]));
	}
	skexp.GetFSI()->Tweak(true,&v);
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
        case 209:
	  skexp.TweakSystematic("xsec_norm_cccoh",tweak); break;
        case 210:
	  skexp.TweakSystematic("xsec_norm_nc1piC",tweak); break;
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
	  else {
	    cout << "Unknown mode " << mode << endl; exit(1); break;
	  }
      }//switch(mode)
    }//mode>=100

    // generate predictions
    skexp.PrintSystematicVariations();
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
  TLegend * legend = new TLegend(0.75, 0.53, 0.95, 0.9);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->AddEntry(hErecoNominal, "nominal", "L" );
  hErecoNominal->SetTitle(Form("%s",modetitle.c_str()));
  hErecoNominal->GetXaxis()->SetRangeUser(0,ereco_max_plot-conv::kOneMeV);
  hErecoNominal->Draw();
  double maxy = hErecoNominal->GetMaximum();
  for(int itweak = 0; itweak < ntweaks; itweak++) {
    if(mode == 201 && itweak != 1) continue; //SF is only +1 sigma!
    hEreco[itweak]->Draw("same");
    maxy = TMath::Max(maxy, hEreco[itweak]->GetMaximum());
    legend->AddEntry(hEreco[itweak], Form("%d sigma",tweaks[itweak]), "L" );
  }
  maxy *= 1.1;
  hErecoNominal->GetYaxis()->SetRangeUser(0,maxy);
  legend->Draw();
  c->cd(2);
  //find the relevant y axis range for the ratio
  double yminratio = 0.75, ymaxratio = 1.25;
  bool foundbig = false;
  for(int itweak = 0; itweak < ntweaks; itweak++) {
    if(mode == 201 && itweak != 1) continue; //SF is only +1 sigma!
    if(foundbig) break;
    for(int ibin = 1; ibin < hErecoRatio[itweak]->GetNbinsX(); ibin++) {
      if(hErecoRatio[itweak]->GetBinContent(ibin) > ymaxratio) {
	ymaxratio = 1.75;
	yminratio = 0.25;
	foundbig = true;
	break;
      }
    }
  }
  for(int itweak = 0; itweak < ntweaks; itweak++) {
    if(itweak==0){
      if(mode == 201) itweak++; //SF is only +1 sigma!
      hErecoRatio[itweak]->GetXaxis()->SetRangeUser(0,ereco_max_plot-conv::kOneMeV);
      hErecoRatio[itweak]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hErecoRatio[itweak]->GetYaxis()->SetTitle("Ratio to nominal");
      hErecoRatio[itweak]->GetYaxis()->SetRangeUser(yminratio,ymaxratio);
      hErecoRatio[itweak]->SetTitle("Ratio");
      hErecoRatio[itweak]->Draw();
      if(mode == 201) break; //SF is only +1 sigma!
      continue;
    }
    hErecoRatio[itweak]->Draw("same");
  }
  c->SaveAs(Form("ereco_twk_%s.eps",modetitleplot.c_str()));
  hErecoNominal->GetXaxis()->UnZoom();
  hErecoRatio[0]->GetXaxis()->UnZoom();
  if(save_all_canvases)
    c->SaveAs(Form("ereco_twk_%s_fullrange.eps",modetitleplot.c_str()));
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
	hErecoModeRatio[im][itweak]->GetYaxis()->SetTitle("Ratio to nominal");
	hErecoModeRatio[im][itweak]->SetTitle("Ratio");
	hErecoModeRatio[im][itweak]->Draw();
	continue;
      }
      hErecoModeRatio[im][itweak]->Draw("same");
    }
    if(save_all_canvases)
      cm->SaveAs(Form("ereco_twk_%s_%s.eps",modetitleplot.c_str(),util::AsStringCompact(m).c_str()));
    hErecoModeNominal[im]->GetXaxis()->UnZoom();
    hErecoModeRatio[im][0]->GetXaxis()->UnZoom();
    if(save_all_canvases)
      cm->SaveAs(Form("ereco_twk_%s_%s_fullrange.eps",modetitleplot.c_str(),util::AsStringCompact(m).c_str()));
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
    int nbins = conv::kNumOfNuisanceParams1Rmu;
    summary_p = new TH1D("summary_p","summary_p",nbins,0,nbins);
    summary_p->SetLineColor(kRed);
    summary_p->GetYaxis()->SetTitle("% change");
    int i = 0;
    for(; i < conv::kNumOfParamsInBANFFfit; i++)
      summary_p->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_pf"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_sf"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_wshape"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_shape_ccoth"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_pilessdelta"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_enushape"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_ccqeE1"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_ccqeE2"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_cc1piE1"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_cccoh"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_nc1piC"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_nccoh"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_ncoth"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_nuetonumu"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_norm_antinutonu"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"sk_energy_scale"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"sk_1Rmu_ccqe_other_eff"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"sk_1Rmu_ccnonqe_eff"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"sk_1Rmu_nc_eff"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"sk_1Rmu_nuecc_eff"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"sk_1Rmu_ccqe_ring_counting_eff_0"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"sk_1Rmu_ccqe_ring_counting_eff_1"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"sk_1Rmu_ccqe_ring_counting_eff_2"); i++;
    summary_p->GetXaxis()->SetBinLabel(i+1"xsec_bindE"); i++;
    for(j = 0; j < conv::kNumOfFSIParams; j++, i++)
      summary_p->GetXaxis()->SetBinLabel(i+1Form("fsi_%03d",j));
  }
  if(!summary_m) {
    int nbins =  conv::kNumOfNuisanceParams1Rmu;
    summary_m = new TH1D("summary_m","summary_m",nbins,0,nbins);
    summary_m->SetLineColor(kBlue);
    summary_m->GetYaxis()->SetTitle("% change");
    int i = 0;
    for(; i < conv::kNumOfParamsInBANFFfit; i++)
      summary_m->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
    summary_m->GetXaxis()->SetBinLabel(conv::kNumOfParamsInBANFFfit+1,"xsec_pf");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_sf");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_wshape");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_shape_ccoth");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_pilessdelta");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_enushape");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ccqeE1");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ccqeE2");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_cc1piE1");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_cccoh");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nc1piC");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nccoh");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ncoth");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nuetonumu");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_norm_antinutonu");
    summary_m->GetXaxis()->SetBinLabel(i+1,"sk_energy_scale");
    summary_m->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_other_eff");
    summary_m->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccnonqe_eff");
    summary_m->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_nc_eff");
    summary_m->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_nuecc_eff");
    summary_m->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_0");
    summary_m->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_1");
    summary_m->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_2");
    summary_m->GetXaxis()->SetBinLabel(i+1,"xsec_bindE");
    for(j = 0; j < conv::kNumOfFSIParams; j++, i++)
      summary_m->GetXaxis()->SetBinLabel(i+1,Form("fsi_%03d",j));
  }
  if(!categories_p) {
    categories_p = new TH1D("categories_p","categories_p",11,0,11);
    categories_p->SetLineColor(kRed);
    categories_p->GetYaxis()->SetTitle("% change");
    categories_p->GetXaxis()->SetBinLabel(i+1,"all");
    categories_p->GetXaxis()->SetBinLabel(i+1,"all(without-energy-scale)");
    categories_p->GetXaxis()->SetBinLabel(i+1,"BANFF");
    categories_p->GetXaxis()->SetBinLabel(i+1,"BANFF-flux");
    categories_p->GetXaxis()->SetBinLabel(i+1,"BANFF-xsec");
    categories_p->GetXaxis()->SetBinLabel(i+1,"unconstrained-xsec");
    categories_p->GetXaxis()->SetBinLabel(i+1,"xsec-all(BANFF+unconstrained)");
    categories_p->GetXaxis()->SetBinLabel(i+1,"sk-detector");
    categories_p->GetXaxis()->SetBinLabel(i+1,"sk-detector(without-energy-scale)");
    categories_p->GetXaxis()->SetBinLabel(i+1,"sk-ccqe-ring-counting-efficiency");
    categories_p->GetXaxis()->SetBinLabel(i+1,"fsi");
  }
  if(!categories_m) {
    categories_m = new TH1D("categories_m","categories_m",11,0,11);
    categories_m->SetLineColor(kRed);
    categories_m->GetYaxis()->SetTitle("% change");
    categories_m->GetXaxis()->SetBinLabel(i+1,"all");
    categories_m->GetXaxis()->SetBinLabel(i+1,"all(without-energy-scale)");
    categories_m->GetXaxis()->SetBinLabel(i+1,"BANFF");
    categories_m->GetXaxis()->SetBinLabel(i+1,"BANFF-flux");
    categories_m->GetXaxis()->SetBinLabel(i+1,"BANFF-xsec");
    categories_m->GetXaxis()->SetBinLabel(i+1,"unconstrained-xsec");
    categories_m->GetXaxis()->SetBinLabel(i+1,"xsec-all(BANFF+unconstrained)");
    categories_m->GetXaxis()->SetBinLabel(i+1,"sk-detector");
    categories_m->GetXaxis()->SetBinLabel(i+1,"sk-detector(without-energy-scale)");
    categories_m->GetXaxis()->SetBinLabel(i+1,"sk-ccqe-ring-counting-efficiency");
    categories_m->GetXaxis()->SetBinLabel(i+1,"fsi");
  }
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    if(!summarymode_p[im]) {
      int nbins = conv::kNumOfNuisanceParams1Rmu;
      summarymode_p[im] = new TH1D(Form("summarymode_p_%s",util::AsStringCompact(m).c_str()),Form("summarymode_p_%s",util::AsStringCompact(m).c_str()),nbins,0,nbins);
      summarymode_p[im]->SetLineColor(kRed);
      summarymode_p[im]->GetYaxis()->SetTitle("% change");
      int i = 0;
      for(; i < conv::kNumOfParamsInBANFFfit; i++)
	summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_pf");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_sf");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_wshape");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_shape_ccoth");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_pilessdelta");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_enushape");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ccqeE1");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ccqeE2");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_cc1piE1");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_cccoh");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nc1piC");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nccoh");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ncoth");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nuetonumu");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_antinutonu");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk_energy_scale");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_other_eff");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccnonqe_eff");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_nc_eff");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_nuecc_eff");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_0");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_1");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_2");
      summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec_bindE");
      for(j = 0; j < conv::kNumOfFSIParams; j++, i++)
	summarymode_p[im]->GetXaxis()->SetBinLabel(i+1,Form("fsi_%03d",j));
    }
    if(!summarymode_m[im]) {
      int nbins = conv::kNumOfNuisanceParams1Rmu;
      summarymode_m[im] = new TH1D(Form("summarymode_m_%s",util::AsStringCompact(m).c_str()),Form("summarymode_m_%s",util::AsStringCompact(m).c_str()),nbins,0,nbins);
      summarymode_m[im]->SetLineColor(kBlue);
      summarymode_m[im]->GetYaxis()->SetTitle("% change");
      int i = 0;
      for(; i < conv::kNumOfParamsInBANFFfit; i++)
	summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_pf");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_sf");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_wshape");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_shape_ccoth");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_pilessdelta");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_enushape");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ccqeE1");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ccqeE2");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_cc1piE1");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_cccoh");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nc1piC");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nccoh");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_ncoth");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_nuetonumu");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_norm_antinutonu");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk_energy_scale");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_other_eff");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccnonqe_eff");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_nc_eff");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_nuecc_eff");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_0");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_1");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk_1Rmu_ccqe_ring_counting_eff_2");
      summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec_bindE");
      for(j = 0; j < conv::kNumOfFSIParams; j++, i++)
	summarymode_m[im]->GetXaxis()->SetBinLabel(i+1,Form("fsi_%03d",j));
    }
    if(!categoriesmode_p[im]) {
      categoriesmode_p[im] = new TH1D(Form("categoriesmode_p_%s",util::AsStringCompact(m).c_str()),Form("categoriesmode_p_%s",util::AsStringCompact(m).c_str()),11,0,11);
      categoriesmode_p[im]->SetLineColor(kRed);
      categoriesmode_p[im]->GetYaxis()->SetTitle("% change");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"all");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"all(without-energy-scale)");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"BANFF");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"BANFF-flux");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"BANFF-xsec");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"unconstrained-xsec");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"xsec-all(BANFF+unconstrained)");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk-detector");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk-detector(without-energy-scale)");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"sk-ccqe-ring-counting-efficiency");
      categoriesmode_p[im]->GetXaxis()->SetBinLabel(i+1,"fsi");
    }
    if(!categoriesmode_m[im]) {
      categoriesmode_m[im] = new TH1D(Form("categoriesmode_m_%s",util::AsStringCompact(m).c_str()),Form("categoriesmode_m_%s",util::AsStringCompact(m).c_str()),11,0,11);
      categoriesmode_m[im]->SetLineColor(kRed);
      categoriesmode_m[im]->GetYaxis()->SetTitle("% change");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"all");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"all(without-energy-scale)");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"BANFF");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"BANFF-flux");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"BANFF-xsec");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"unconstrained-xsec");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"xsec-all(BANFF+unconstrained)");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk-detector");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk-detector(without-energy-scale)");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"sk-ccqe-ring-counting-efficiency");
      categoriesmode_m[im]->GetXaxis()->SetBinLabel(i+1,"fsi");
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
  
