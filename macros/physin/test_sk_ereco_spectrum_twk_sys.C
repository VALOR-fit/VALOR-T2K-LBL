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
//   6 tweak all sk detector
//   7 tweak all sk ccqe ring counting efficiency
//   8 tweak all fsi
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
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

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
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

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
  const int ntweaks = 4;
  int * tweaks[ntweaks] = {+3, +1, -1, -3};
  TH1D * hErecoNominal = 0;
  TH1D * hEreco[ntweaks];
  TH1D * hErecoRatio[ntweaks];

  // Get the nominal
  skexp.Generate(ph);
  hErecoNominal = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  hErecoNominal->SetDirectory(0);
  hErecoNominal->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoNominal->GetYaxis()->SetTitle("Number of events");
  hErecoNominal->SetLineWidth(2);

  // Open output file
  string rootmode = "update";
  if(!update_rootfile)
    rootmode = "recreate";
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_sys_%s.root",oscparams.c_str()),rootmode.c_str());

  string modetitle = "";
  switch(mode) {
    case: 0
      moodetitle = "all"; break;
    case: 1
      moodetitle = "BANFF"; break;
    case: 2
      moodetitle = "BANFF flux"; break;
    case: 3
      moodetitle = "BANFF xsec"; break;
    case: 4
      moodetitle = "unconstrained xsec"; break;
    case: 5
      moodetitle = "all xsec (BANFF + unconstrained)"; break;
    case: 6
      moodetitle = "sk detector"; break;
    case: 7
      moodetitle = "sk ccqe ring counting efficiency"; break;
    case: 8
      moodetitle = "fsi"; break;
    default: 8
      moodetitle = ""; break;
  }

  // Tweak!
  for(int itweak = 0; itweak < ntweaks; itweak++;) {
    int tweak = tweaks[itweak];

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

    if(mode<0 || (mode>8)// && mode<100) || (mode>150 && mode<200) || (mode>213 && mode<300) || (mode>307 && mode<400) || mode>415) {
      cout << "Invalid mode " << mode << endl;
      exit(1);
    }

    if(mode == 0 || mode == 1) {
      TVectorD v(conv::kNumOfParamsInBANFFfit);
      v = skexp.GetBANFF()->BestFitValues();
      for(int ib = 0; ib < conv::kNumOfParamsInBANFFfit; ib++) {
	v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * skexp.GetBANFF()->CovMatrix()[ib][ib]);
      }
      skexp.GetBANFF()->Tweak(true,&v);
    }

    if(mode == 2) {
      TVectorD v(conv::kNumOfParamsInBANFFfit);
      v = skexp.GetBANFF()->BestFitValues();
      for(int ib = 0; ib < conv::kNumOfFluxParamsInBANFFfit; ib++) {
	v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * skexp.GetBANFF()->CovMatrix()[ib][ib]);
      }
      skexp.GetBANFF()->Tweak(true,&v);
    }

    if(mode == 3 || mode == 5) {
      TVectorD v(conv::kNumOfParamsInBANFFfit);
      v = skexp.GetBANFF()->BestFitValues();
      for(int ib = conv::kNumOfFluxParamsInBANFFfit; ib < conv::kNumOfParamsInBANFFfit; ib++) {
	v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * skexp.GetBANFF()->CovMatrix()[ib][ib]);
      }
      skexp.GetBANFF()->Tweak(true,&v);
    }

    if(mode == 0 || mode == 4 || mode == 5) {
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
    }

    if(mode == 0 || mode == 6) {
      skexp.TweakSystematic("sk_energy_scale",tweak);
      skexp.TweakSystematic("sk_1Rmu_ccqe_other_eff",tweak);
      skexp.TweakSystematic("sk_1Rmu_ccnonqe_eff",tweak);
      skexp.TweakSystematic("sk_1Rmu_nc_eff",tweak);
      skexp.TweakSystematic("sk_1Rmu_nuecc_eff",tweak);
    }

    if(mode == 0 || mode == 6 || mode == 7) {
      skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak);
      skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak);
      skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak);
    }

    /*
    if(mode == 0 || mode == 8) {
      TVectorD v(conv::kNumOfFSIParams);
      v.Zero();
      for(int ib = 0; ib < conv::kNumOfFSIParams; ib++) {
	v[ib] = tweak * skexp.GetFSI()->CovMatrix()[ib][ib];
      }
      skexp.GetFSI()->Tweak(true,&v);
    }
    */

    // generate predictions
    skexp.Generate(ph);
    //store spectrum
    hEreco[itweak] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
    hEreco[itweak]->SetDirectory(0);
    hEreco[itweak]->SetTitle(Form("%s_%d_sigma",modetitle,tweak));
    hEreco[itweak]->SetLineColor(histcolor);
    hEreco[itweak]->Write(Form("%s_%d_sigma",modetitle,tweak));		    
    hErecoRatio[itweak] = new TH1D(*hEreco[itweak]);
    hErecoRatio[itweak]->Divide(hErecoNominal);
    hErecoRatio[itweak]->Write(Form("%s_%d_sigma_ratio",modetitle,tweak));
	
  }//end loop over tweaks

  //
  // Draw
  //
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->Divide(1,2);
  c->cd(1);
  hErecoNominal->Draw();
  for(int itweak = 0; itweak < ntweaks; itweak++) {
    hEreco[itweak]->Draw("same");
  }
  c->cd(2);
  for(int itweak = 0; itweak < ntweaks; itweak++) {
    if(itweak==0){
      hErecoRatio[itweak]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hErecoRatio[itweak]->GetYaxis()->SetTitle("Relative error");
      hErecoRatio[itweak]->Draw();
      continue;
    }
    hErecoRatio[itweak]->Draw("same");
  }
  c->SaveAs(Form("ereco_twk_%s.png",modetitle));

  f.Close();
}
