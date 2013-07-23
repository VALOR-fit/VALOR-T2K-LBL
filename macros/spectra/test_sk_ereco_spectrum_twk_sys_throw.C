//
// ------------------------------------------------------------------------------------------------
// Calculate SuperK Ereco spectrum error envelope due to systematic parameter variations,
// randomised using Gaussians or Cholesky decomposition (for FSI/BANFF/CCQERingCounting)
//
// There are many run modes:
//
//   0 tweak all parameters
//   1 tweak all BANFF
//   2 tweak all unconstrained xsec
//   3 tweak all sk detector
//   4 tweak all fsi
//   7 tweak BANFF flux only
//   8 tweak BANFF xsec only
//   9 tweak SKdet + FSI
//  10 tweak BANFF + unconstrained xsec
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
// 214 tweak xsec_norm_antinutonu
// 215 tweak xsec_bindE
//
// 300 tweak sk_energy_scale
// 301 tweak sk_1Rmu_ccqe_other_eff
// 302 tweak sk_1Rmu_ccnonqe_eff
// 303 tweak sk_1Rmu_nc_eff
// 304 tweak sk_1Rmu_nuecc_eff
// 305 tweak sk_1Rmu_ccqe_ring_counting
// 306 tweak sk_1Rmu_ccqe_ring_counting0
// 307 tweak sk_1Rmu_ccqe_ring_counting1
// 308 tweak sk_1Rmu_ccqe_ring_counting2
//
// 400 tweak fsi00
// 401 tweak fsi01
// ...
// 400+kNumOfFSIParams tweak last fsi parameter
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void test_sk_ereco_spectrum_twk_sys_throw
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
  bool         save_toy_expts     = false,                          /* save the ereco spectra? */
  const char * config_file        = "analysis_configuration.xml",   /* */
  const int    nexp               = 20000,                          /* number of experiments to make */
  bool         no_spectra         = true                            /* if true, won't try to get the spectra (for PRODUCTIONMODE) */
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
      modetitle = "unconstrained_xsec"; break;
    case 3:
      modetitle = "sk_detector"; break;
    case 4:
      char fsibins[40] = "";
      sprintf(fsibins, "fsi_%d", conv::kNumOfFSIParams);
      modetitle = fsibins; break;
    case 7:
      modetitle = "BANFF_flux"; break;
    case 8:
      modetitle = "BANFF_xsec"; break;
    case 9:
      modetitle = "SKDet_FSI"; break;
    case 10:
      modetitle = "BANFF_unconstrained_xsec"; break;

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
    case 209:
      modetitle = "xsec_norm_cccoh"; break;
    case 210:
      modetitle = "xsec_norm_nc1piC"; break;
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
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff"; break;
    case 306:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_0"; break;
    case 307:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_1"; break;
    case 308:
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff_2"; break;


    default:
      if(mode >= 100 && mode < 100+conv::kNumOfParamsInBANFFfit) {
	modetitle = Form("banff_%02d",mode-100); break;
      }
      else if(mode >= 400 && mode < 400+conv::kNumOfFSIParams) {
	modetitle = Form("fsi_%03d",mode-400); break;
      }
      cout << "Unknown mode " << mode << endl; exit(1); break;
  }//switch mode

  TString config_str(config_file);
  if(config_str.Contains("prefit"))
    modetitle += "prefit";

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  //Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
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
  
  // Output histograms
  TH1D * hErecoNominal = 0;
  TH1D * hEreco[nexp] = {0};
  TH1D * hNsk = new TH1D(Form("%s_nsk",modetitle.c_str()),Form("%s_nsk",modetitle.c_str()),10000,0,1000);
  int nbanff = conv::kNumOfParamsInBANFFfit;
  //int nfsi = conv::kNumOfFSIParams;
  int nfsi = conv::kNumOfSKDetFSIParams;
  const int nsys = nbanff + nfsi + 8 + 13;
  TH1D * hSysThrown[nsys];
  for(int isys = 0; isys < nsys; isys++) {
    string sysstr = "";
    if(isys < nbanff)
      sysstr = Form("banff_%02i",isys);
    else if(isys < nbanff + nfsi)
      sysstr = Form("fsi_%03i",isys - nbanff);
    else {
      int isys2 = isys - nbanff - nfsi;
      switch(isys2) {
        case(0) :
	  sysstr = "sk_energy_scale"; break;
        case(1) :
	  sysstr = "sk_ccqe_rc_eff_0"; break;
        case(2) :
	  sysstr = "sk_ccqe_rc_eff_1"; break;
        case(3) :
	  sysstr = "sk_ccqe_rc_eff_2"; break;
        case(4) :
	  sysstr = "sk_1Rmu_ccqe_other_eff"; break;
        case(5) :
	  sysstr = "sk_1Rmu_ccnonqe_eff"; break;
        case(6) :
	  sysstr = "sk_1Rmu_nc_eff"; break;
        case(7) :
	  sysstr = "sk_1Rmu_nuecc_eff"; break;
        case(8) :
	  sysstr = "xsec_pf"; break;
        case(9) :
	  sysstr = "xsec_sf"; break;
        case(10) :
	  sysstr = "xsec_wshape"; break;
        case(11) :
	  sysstr = "xsec_shape_ccoth"; break;
        case(12) :
	  sysstr = "xsec_pilessdelta"; break;
        case(14) :
	  sysstr = "xsec_norm_cccoh"; break;
        case(15) :
	  sysstr = "xsec_norm_nc1piC"; break;
        case(17) :
	  sysstr = "xsec_norm_ncoth"; break;
        case(18) :
	  sysstr = "xsec_norm_nuetonumu"; break;
        case(19) :
	  sysstr = "xsec_norm_antinutonu"; break;
        case(20) :
	  sysstr = "xsec_bindE"; break;
        default:
	  sysstr = ""; break;
      }//switch(isys2)
    }//isys>=nbanff+nfsi
    hSysThrown[isys] = new TH1D(Form("%s_tweak_throws_%s",modetitle.c_str(),sysstr.c_str()),Form("%s_tweak_throws_%s",modetitle.c_str(),sysstr.c_str()),1000,-5,5);
    hSysThrown[isys]->SetDirectory(0);
  }

  // Open output files
  string rootmode = "update";
  if(!update_rootfile)
    rootmode = "recreate";
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_sys_throw_%s.root",oscparams.c_str()),rootmode.c_str());

  if(!update_rootfile)
    ofstream outfile(Form("sk_ereco_spectrum-twk_sys_throw_%s.txt",oscparams.c_str()),ios_base::trunc);
  else
    ofstream outfile(Form("sk_ereco_spectrum-twk_sys_throw_%s.txt",oscparams.c_str()),ios_base::app);

  //remove the [] so can use it in SaveAs
  while(oscparams.find_first_of("[]") != string::npos) {
    oscparams.erase(oscparams.find_first_of("[]"), 1);
  }

  // Get the nominal
  skexp.Generate(ph);
  if(!no_spectra) {
    hErecoNominal = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
    hErecoNominal->SetDirectory(0);
    hErecoNominal->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hErecoNominal->GetYaxis()->SetTitle("Events per bin");
    hErecoNominal->SetLineWidth(2);
    if(!update_rootfile)
      hErecoNominal->Write("nominal");
  }

  TStopwatch w;
  w.Start();

  // Tweak!
  int iskip = 0;
  for(int iexp = 0; iexp < nexp; iexp++) {

    cout << endl << "Experiment " << iexp << " for systematics: " << modetitle << endl;

    if(mode < 100) {

      // All
      if(mode == 0) {
	skexp.TweakAllSystematics();
      }

      // BANFF
      if(mode == 1) {
	skexp.TweakAllBANFFSystematics();
      }

      // Uncorrelated XSec
      if(mode == 2) {
	skexp.TweakAllUncorrelatedXSecSystematics();
      }

      // SK detector
      if(mode == 3) {
	skexp.TweakAllSKDetectorSystematics();
      }

      // FSI
      if(mode == 4) {
	skexp.TweakAllFSISystematics();
	//skexp.GetFSI()->PrintSystematicVariations();
      }

      //BANFF (flux only)
      if(mode == 7) {
	TRandom3 rand(0);
	for(int ib = 0; ib < conv::kNumOfFluxParamsInBANFFfit; ib++) {
	  skexp.TweakSystematic(Form("banff_%02d",ib),rand.Gaus(0.,1.));
	}
      }

      //BANFF (xsec only)
      if(mode == 8) {
	TRandom3 rand(0);
	skexp.TweakSystematic("banff_maqe",rand.Gaus(0.,1.));
	skexp.TweakSystematic("banff_mares",rand.Gaus(0.,1.));
	skexp.TweakSystematic("banff_ccqeE0",rand.Gaus(0.,1.));
	skexp.TweakSystematic("banff_ccqeE1",rand.Gaus(0.,1.));
	skexp.TweakSystematic("banff_ccqeE2",rand.Gaus(0.,1.));
	skexp.TweakSystematic("banff_cc1piE0",rand.Gaus(0.,1.));
	skexp.TweakSystematic("banff_cc1piE1",rand.Gaus(0.,1.));
	skexp.TweakSystematic("banff_mecnorm",rand.Gaus(0.,1.));
	skexp.TweakSystematic("banff_pdd",rand.Gaus(0.,1.));
      }

      //SK det + FSI
      if(mode == 9) {
	//skexp.TweakAllSKDetectorSystematics();
	//skexp.TweakAllFSISystematics();
	skexp.TweakAllSKDetectorFSISystematics();
      }

     //unconstrained xsec + BANFF
      if(mode == 10) {
	skexp.TweakAllUncorrelatedXSecSystematics();
	skexp.TweakAllBANFFSystematics();
      }

    }//mode<100

    else {
      
      TRandom3 rand(0);
      double tweak = rand.Gaus(0.,1.);

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
	  skexp.GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(); break;
        case 306:
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_0",tweak); break;
        case 307:
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_1",tweak); break;
        case 308:
	  skexp.TweakSystematic("sk_1Rmu_ccqe_ring_counting_eff_2",tweak); break;


        default:
	  if(mode >= 100 && mode < 100+conv::kNumOfParamsInBANFFfit) {
	    TVectorD v(conv::kNumOfParamsInBANFFfit);
	    v = skexp.GetBANFF()->BestFitValues();
	    int ib = mode - 100;
	    v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
	    skexp.GetBANFF()->Tweak(true,&v,1);
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
	  cout << "Unknown mode " << mode << endl; exit(1); break;
      }//switch(mode)
    }//mode>=100

    // generate predictions
    bool applied_negevents_condition = skexp.Generate(ph);
    if(applied_negevents_condition) {
      cout << "skipping" << endl;
      iskip++;
      continue;
    }

    //store spectrum
    if(!no_spectra) {
      hEreco[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEreco[iexp]->SetDirectory(0);
      hEreco[iexp]->SetTitle(Form("%s_exp%d",modetitle.c_str(),iexp));
      if(save_toy_expts)
	hEreco[iexp]->Write(Form("%s_exp%d",modetitle.c_str(),iexp));		    
    }

    //store Nsk
    double nskint = 0;
    for(int ir = 1; ir <= conv::kNErecoBins1Rmu; ir++)
      nskint += skexp.CurrPredictedErecoSpectrumArray(s)[ir-1];
    hNsk->Fill(nskint);

    //store thrown parameter values
    for(int isys = 0; isys < nbanff; isys++) {
      hSysThrown[isys]->Fill(skexp.GetBANFF()->CurrentValue(conv::kNuisBANFF00 + isys));
    }
    //for(int isys = 0; isys < nfsi; isys++) {
      //hSysThrown[isys+nbanff]->Fill(skexp.GetFSI()->CurrentValue(util::PhysParamIdFromString(Form("fsi_%03d",isys))));
    for(int isys = 0; isys < nfsi-1; isys++) {
      hSysThrown[isys+nbanff]->Fill(skexp.GetSuperkDetectorFSI()->CurrentValue(util::PhysParamIdFromString(Form("skdetfsi_%03d",isys))));
    }
    hSysThrown[nbanff+nfsi ]->Fill(skexp.GetSuperkDetectorFSI()->CurrentValue(util::PhysParamIdFromString("sk_energy_scale")));
    /*
    hSysThrown[nbanff+nfsi+0 ]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_energy_scale")));
    hSysThrown[nbanff+nfsi+1 ]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccqe_ring_counting_eff_0")));
    hSysThrown[nbanff+nfsi+2 ]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccqe_ring_counting_eff_1")));
    hSysThrown[nbanff+nfsi+3 ]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccqe_ring_counting_eff_2")));
    hSysThrown[nbanff+nfsi+4 ]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccqe_other_eff")));
    hSysThrown[nbanff+nfsi+5 ]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccnonqe_eff")));
    hSysThrown[nbanff+nfsi+6 ]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_nc_eff")));
    hSysThrown[nbanff+nfsi+7 ]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_nuecc_eff")));
    */
    hSysThrown[nbanff+nfsi+8 ]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_pf")));
    hSysThrown[nbanff+nfsi+9 ]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_sf")));
    hSysThrown[nbanff+nfsi+10]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_wshape")));
    hSysThrown[nbanff+nfsi+11]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_shape_ccoth")));
    hSysThrown[nbanff+nfsi+12]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_pilessdelta")));
    hSysThrown[nbanff+nfsi+13]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_cccoh")));
    hSysThrown[nbanff+nfsi+14]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_nc1piC")));
    hSysThrown[nbanff+nfsi+15]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_ncoth")));
    hSysThrown[nbanff+nfsi+16]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_nuetonumu")));
    hSysThrown[nbanff+nfsi+17]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_antinutonu")));
    hSysThrown[nbanff+nfsi+18]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_bindE")));
  }//end loop over experiments
  w.Stop();
  cout << iskip << " experiments skipped due to min events" << endl;

  //save some histograms
  for(int isys = 0; isys < nsys; isys++) {
    hSysThrown[isys]->Write();
  }
  hNsk->Write();

  // Find min/max and calculate fractional deviation from nominal
  if(!no_spectra) {
    hErecoMax = new TH1D(*hErecoNominal);
    hErecoMin = new TH1D(*hErecoNominal);
    hErecoMax->SetDirectory(0);
    hErecoMin->SetDirectory(0);
    hErecoMaxEnv = new TH1D(*hErecoNominal);
    hErecoMinEnv = new TH1D(*hErecoNominal);
    hErecoMaxEnv->SetDirectory(0);
    hErecoMinEnv->SetDirectory(0);
    for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
      double nom   = hErecoNominal->GetBinContent(ibin);
      double sumsq = 0;
      for(int iexp = 0; iexp < nexp; iexp++) {
	if(hEreco[iexp]) {
	  double N = hEreco[iexp]->GetBinContent(ibin);
	  sumsq += TMath::Power(N-nom,2);
	}
      }//iexp
      double sig = TMath::Sqrt(sumsq/nexp);
      double min = nom-sig;
      double max = nom+sig;
      if(nom>0){
	hErecoMax -> SetBinContent (ibin,  (max-nom)/nom);
	hErecoMin -> SetBinContent (ibin,  (min-nom)/nom);
	hErecoMaxEnv -> SetBinContent (ibin,  max);
	hErecoMinEnv -> SetBinContent (ibin,  min);
      }
    }//ibin
    hErecoMax->Write(Form("%s_pluserr",modetitle.c_str()));
    hErecoMin->Write(Form("%s_minuserr",modetitle.c_str()));
    hErecoMaxEnv->Write(Form("%s_envpluserr",modetitle.c_str()));
    hErecoMinEnv->Write(Form("%s_envminuserr",modetitle.c_str()));
  }

  //Draw & save
  TCanvas * c = new TCanvas();
  hNsk->Draw();
  c->SaveAs(Form("%s_nsk_%s.eps",modetitle.c_str(),oscparams.c_str()));
  if(!no_spectra) {
    hErecoMax->Draw();
    c->SaveAs(Form("%s_pluserr_%s.eps",modetitle.c_str(),oscparams.c_str()));
    hErecoMin->Draw();
    c->SaveAs(Form("%s_minuserr_%s.eps",modetitle.c_str(),oscparams.c_str()));
    hErecoMaxEnv->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
    hErecoMaxEnv->SetLineStyle(2);
    hErecoMinEnv->SetLineStyle(2);
    hErecoMaxEnv->SetLineColor(kBlue+3);
    hErecoMinEnv->SetLineColor(kBlue+3);
    hErecoMaxEnv->Draw();
    hErecoMinEnv->Draw("same");
    hErecoNominal->Draw("same");
    c->SaveAs(Form("%s_envelope_%s.eps",modetitle.c_str(),oscparams.c_str()));
  }

  TCanvas * csys = new TCanvas();
  csys->Divide(ceil(sqrt(nsys)),ceil(sqrt(nsys)));
  for(int isys=0; isys < nsys; isys++) {
    csys->cd(isys+1);
    hSysThrown[isys]->Draw();
    cout << isys << "\t" << hSysThrown[isys]->GetTitle() << endl;
  }
  csys->SaveAs(Form("%s_throws_%s.eps",modetitle.c_str(),oscparams.c_str()));
  delete c;
  delete csys;

  f.Close();

  // format info to write
  stringstream streamer;
  if(!update_rootfile) {
    streamer << "\\begin{table}"   << endl;
    streamer << "\\begin{tabular}{| c | c c |}" << endl;
    streamer << "\\hline" << endl;
    streamer << "Nsk RMS/Mean &  "
	     << "BANFF pre "
	     << "BANFF post "
	     << endl;
    streamer << "\\hline" << endl;
  }
  string banffpre = "-";
  if(mode == 0 || mode == 1)
    banffpre = "?";
  streamer << modetitle << " & & " << hNsk->GetRMS() / hNsk->GetMean() << " \\\\ \\hline" << endl;
  cout << streamer.str() << endl;

  // get size of file
  long size;
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);
  // write to outfile
  outfile.write (streamer.str().c_str(),size);
  outfile.close();

  cout << "Tweaks performed in " 
       << w.CpuTime() << " sec." << endl;
}
