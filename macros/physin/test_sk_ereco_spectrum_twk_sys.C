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
  double       ereco_max_plot     = 7,                              /* maximum of the ereco axis (in GeV) to plot */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

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

    default:
      if(mode >= 100 && mode < 100+conv::kNumOfParamsInBANFFfit) {
	modetitle = Form("banff_%02d",mode-100); break;
      }
      else if(mode >= 400 && mode < 400+conv::kNumOfFSIParams) {
	modetitle = Form("fsi_%02d",mode-400); break;
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
  const int ntweaks = 4;
  int * tweaks[ntweaks] = {+3, +1, -1, -3};
  TH1D * hErecoNominal = 0;
  TH1D * hEreco[ntweaks];
  TH1D * hErecoRatio[ntweaks];

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
  if(!update_rootfile)
    hErecoNominal->Write("nominal");

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
      if(mode == 0 || mode == 1) {
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

      if(mode == 0 || mode == 6 || mode == 8) {
	if(mode == 0 || mode == 8)
	  skexp.TweakSystematic("sk_energy_scale",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_other_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccnonqe_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_nc_eff",tweak);
	skexp.TweakSystematic("sk_1Rmu_nuecc_eff",tweak);
      }

      if(mode == 0 || mode == 6 || mode == 7 || mode == 8) {
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak);
	skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak);
      }

      if(mode == 0 || mode == 9) {
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
  c->SaveAs(Form("ereco_twk_%s.png",modetitle.c_str()));

  //
  // Save the change in nsk
  //
  double nsk_0 = hErecoNominal->Integral(1,hErecoNominal->GetNbinsX());
  double nsk_p  = hEreco[1]->Integral(1,hEreco[1]->GetNbinsX());
  double nsk_m  = hEreco[2]->Integral(1,hEreco[2]->GetNbinsX());
  double frac_nsk_p = (nsk_0 > 0) ? 100.*(nsk_p - nsk_0) / nsk_0 : 0.;
  double frac_nsk_m = (nsk_0 > 0) ? 100.*(nsk_m - nsk_0) / nsk_0 : 0.;

  cout << endl << "Nominal number of events: " << nsk_0 << endl
       << "   +1 sigma:              " << nsk_p << "\t" << frac_nsk_p << "%" << endl
       << "   -1 sigma:              " << nsk_m << "\t" << frac_nsk_m << "%" << endl;
  if(mode>=100) {
    TH1D * summary_p = 0, * summary_m = 0;
    summary_p = f.GetObject("summary_p",summary_p);
    summary_m = f.GetObject("summary_m",summary_m);

    if(!summary_p) {
      summary_p = new TH1D("summary_p","summary_p",77,0,77);
      summary_p->SetLineColor(kRed);
      for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
	summary_p->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
      summary_p->GetXaxis()->SetBinLabel(40,"xsec_pf");
      summary_p->GetXaxis()->SetBinLabel(41,"xsec_sf");
      summary_p->GetXaxis()->SetBinLabel(42,"xsec_wshape");
      summary_p->GetXaxis()->SetBinLabel(43,"xsec_shape_ccoth");
      summary_p->GetXaxis()->SetBinLabel(44,"xsec_pilessdelta");
      summary_p->GetXaxis()->SetBinLabel(45,"xsec_enushape");
      summary_p->GetXaxis()->SetBinLabel(46,"xsec_norm_ccqeE1");
      summary_p->GetXaxis()->SetBinLabel(47,"xsec_norm_ccqeE2");
      summary_p->GetXaxis()->SetBinLabel(48,"xsec_norm_cc1piE1");
      summary_p->GetXaxis()->SetBinLabel(49,"xsec_norm_cccoh");
      summary_p->GetXaxis()->SetBinLabel(50,"xsec_norm_nc1piC");
      summary_p->GetXaxis()->SetBinLabel(51,"xsec_norm_nccoh");
      summary_p->GetXaxis()->SetBinLabel(52,"xsec_norm_ncoth");
      summary_p->GetXaxis()->SetBinLabel(53,"xsec_norm_nuetonumu");
      summary_p->GetXaxis()->SetBinLabel(54,"sk_energy_scale");
      summary_p->GetXaxis()->SetBinLabel(55,"sk_1Rmu_ccqe_other_eff");
      summary_p->GetXaxis()->SetBinLabel(56,"sk_1Rmu_ccnonqe_eff");
      summary_p->GetXaxis()->SetBinLabel(57,"sk_1Rmu_nc_eff");
      summary_p->GetXaxis()->SetBinLabel(58,"sk_1Rmu_nuecc_eff");
      summary_p->GetXaxis()->SetBinLabel(59,"sk_1Rmu_ccqe_ring_counting_eff_0");
      summary_p->GetXaxis()->SetBinLabel(60,"sk_1Rmu_ccqe_ring_counting_eff_1");
      summary_p->GetXaxis()->SetBinLabel(61,"sk_1Rmu_ccqe_ring_counting_eff_2");
      for(int i = 0; i < conv::kNumOfFSIParams; i++)
	summary_p->GetXaxis()->SetBinLabel(i+62,Form("fsi_%02d",i));
    }
    if(!summary_m) {
      summary_m = new TH1D("summary_m","summary_m",77,0,77);
      summary_m->SetLineColor(kBlue);
      for(int i = 0; i < conv::kNumOfParamsInBANFFfit; i++)
	summary_m->GetXaxis()->SetBinLabel(i+1,Form("banff_%02d",i));
      summary_m->GetXaxis()->SetBinLabel(40,"xsec_pf");
      summary_m->GetXaxis()->SetBinLabel(41,"xsec_sf");
      summary_m->GetXaxis()->SetBinLabel(42,"xsec_wshape");
      summary_m->GetXaxis()->SetBinLabel(43,"xsec_shape_ccoth");
      summary_m->GetXaxis()->SetBinLabel(44,"xsec_pilessdelta");
      summary_m->GetXaxis()->SetBinLabel(45,"xsec_enushape");
      summary_m->GetXaxis()->SetBinLabel(46,"xsec_norm_ccqeE1");
      summary_m->GetXaxis()->SetBinLabel(47,"xsec_norm_ccqeE2");
      summary_m->GetXaxis()->SetBinLabel(48,"xsec_norm_cc1piE1");
      summary_m->GetXaxis()->SetBinLabel(49,"xsec_norm_cccoh");
      summary_m->GetXaxis()->SetBinLabel(50,"xsec_norm_nc1piC");
      summary_m->GetXaxis()->SetBinLabel(51,"xsec_norm_nccoh");
      summary_m->GetXaxis()->SetBinLabel(52,"xsec_norm_ncoth");
      summary_m->GetXaxis()->SetBinLabel(53,"xsec_norm_nuetonumu");
      summary_m->GetXaxis()->SetBinLabel(54,"sk_energy_scale");
      summary_m->GetXaxis()->SetBinLabel(55,"sk_1Rmu_ccqe_other_eff");
      summary_m->GetXaxis()->SetBinLabel(56,"sk_1Rmu_ccnonqe_eff");
      summary_m->GetXaxis()->SetBinLabel(57,"sk_1Rmu_nc_eff");
      summary_m->GetXaxis()->SetBinLabel(58,"sk_1Rmu_nuecc_eff");
      summary_m->GetXaxis()->SetBinLabel(59,"sk_1Rmu_ccqe_ring_counting_eff_0");
      summary_m->GetXaxis()->SetBinLabel(60,"sk_1Rmu_ccqe_ring_counting_eff_1");
      summary_m->GetXaxis()->SetBinLabel(61,"sk_1Rmu_ccqe_ring_counting_eff_2");
      for(int i = 0; i < conv::kNumOfFSIParams; i++)
	summary_m->GetXaxis()->SetBinLabel(i+62,Form("fsi_%02d",i));
    }

    summary_p->Fill(modetitle.c_str(),frac_nsk_p);
    summary_m->Fill(modetitle.c_str(),frac_nsk_m);
    summary_p->Write();
    summary_m->Write();
  }//mode>=100

  f.Close();
}

