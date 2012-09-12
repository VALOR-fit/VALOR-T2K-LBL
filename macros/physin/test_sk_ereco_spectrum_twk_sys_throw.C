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
// 305 tweak sk_1Rmu_ccqe_ring_counting
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
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  const int nexp = 20000; //number of experiments to make

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

  //Check if the mode is valid
  string modetitle = "";
  switch(mode) {
    case 0:
      modetitle = "all"; break;
    case 1:
      modetitle = "BANFF"; break;
    case 2:
      modetitle = "unconstrained-xsec"; break;
    case 3:
      modetitle = "sk-detector"; break;
    case 4:
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
      modetitle = "sk_1Rmu_ccqe_ring_counting_eff"; break;

    default:
      cout << "Unknown mode " << mode << endl; exit(1); break;
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
  
  // Output histograms
  TH1D * hErecoNominal = 0;
  TH1D * hEreco[nexp];
  TH1D * hNsk = new TH1D(Form("%s_nsk",modetitle.c_str()),Form("%s_nsk",modetitle.c_str()),300,0,300);
  int nbanff = conv::kNumOfParamsInBANFFfit;
  int nfsi = conv::kNumOfFSIParams;
  const int nsys = nbanff + nfsi + 8 + 11;
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
	  sysstr = "sk_1Rmu_ccqe_ring_counting_eff_0"; break;
        case(2) :
	  sysstr = "sk_1Rmu_ccqe_ring_counting_eff_1"; break;
        case(3) :
	  sysstr = "sk_1Rmu_ccqe_ring_counting_eff_2"; break;
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
        case(13) :
	  sysstr = "xsec_enushape"; break;
        case(14) :
	  sysstr = "xsec_norm_cccoh"; break;
        case(15) :
	  sysstr = "xsec_norm_nc1piC"; break;
        case(16) :
	  sysstr = "xsec_norm_nccoh"; break;
        case(17) :
	  sysstr = "xsec_norm_ncoth"; break;
        case(18) :
	  sysstr = "xsec_norm_nuetonumu"; break;
        default:
	  sysstr = ""; break;
      }//switch(isys2)
    }//isys>=nbanff+nfsi
    hSysThrown[isys] = new TH1D(Form("%s_tweak_throws_%s",modetitle.c_str(),sysstr.c_str()),Form("%s_tweak_throws_%s",modetitle.c_str(),sysstr.c_str()),100,-5,5);
    hSysThrown[isys]->SetDirectory(0);
  }

  // Open output file
  string rootmode = "update";
  if(!update_rootfile)
    rootmode = "recreate";
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_sys_throw_%s.root",oscparams.c_str()),rootmode.c_str());

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
  for(int iexp = 0; iexp < nexp; iexp++) {

    cout << endl << "Experiment " << iexp << " for systematics: " << modetitle << endl;

    if(mode < 100) {

      // BANFF
      if(mode == 0 || mode == 1) {
	skexp.TweakAllBANFFSystematics();
      }

      // Uncorrelated XSec
      if(mode == 0 || mode == 2) {
	skexp.TweakAllUncorrelatedXSecSystematics();
      }

      // SK detector
      if(mode == 0 || mode == 3) {
	skexp.TweakAllSKDetectorSystematics();
      }

      // FSI
      if(mode == 0 || mode == 4) {
	skexp.TweakAllFSISystematics();
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
	  skexp.GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(); break;

        default:
	  cout << "Unknown mode " << mode << endl; exit(1); break;
      }//switch(mode)
    }//mode>=100

    // generate predictions
    skexp.Generate(ph);

    //store spectrum
    hEreco[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
    hEreco[iexp]->SetDirectory(0);
    hEreco[iexp]->SetTitle(Form("%s_exp%d",modetitle.c_str(),iexp));
    hEreco[iexp]->Write(Form("%s_exp%d",modetitle.c_str(),iexp));		    
	
    //store Nsk
    hNsk->Fill(hEreco[iexp]->Integral());

    //store thrown parameter values
    nbflux = conv::kNumOfFluxParamsInBANFFfit;
    for(int isys = 0; isys < nbanff; isys++) {
      hSysThrown[isys]->Fill(skexp.GetBANFF()->CurrentValue(conv::kNuisBANFF00 + isys));
    }
    for(int isys = 0; isys < nfsi; isys++) {
      hSysThrown[isys+nbanff]->Fill(skexp.GetFSI()->CurrentValue(util::PhysParamIdFromString(Form("fsi_%03d",isys))));
    }
    hSysThrown[nbanff+nfsi+0]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_energy_scale")));
    hSysThrown[nbanff+nfsi+1]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccqe_ring_counting_eff_0")));
    hSysThrown[nbanff+nfsi+2]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccqe_ring_counting_eff_1")));
    hSysThrown[nbanff+nfsi+3]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccqe_ring_counting_eff_2")));
    hSysThrown[nbanff+nfsi+4]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccqe_other_eff")));
    hSysThrown[nbanff+nfsi+5]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_ccnonqe_eff")));
    hSysThrown[nbanff+nfsi+6]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_nc_eff")));
    hSysThrown[nbanff+nfsi+7]->Fill(skexp.GetSuperkDetector()->CurrentValue(util::PhysParamIdFromString("sk_1Rmu_nuecc_eff")));
    hSysThrown[nbanff+nfsi+8]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_pf")));
    hSysThrown[nbanff+nfsi+9]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_sf")));
    hSysThrown[nbanff+nfsi+10]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_wshape")));
    hSysThrown[nbanff+nfsi+11]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_shape_ccoth")));
    hSysThrown[nbanff+nfsi+12]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_pilessdelta")));
    hSysThrown[nbanff+nfsi+13]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_enushape")));
    hSysThrown[nbanff+nfsi+14]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_cccoh")));
    hSysThrown[nbanff+nfsi+15]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_nc1piC")));
    hSysThrown[nbanff+nfsi+16]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_nccoh")));
    hSysThrown[nbanff+nfsi+17]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_ncoth")));
    hSysThrown[nbanff+nfsi+18]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("xsec_norm_nuetonumu")));
  }//end loop over experiments

  //save some histograms
  for(int isys = 0; isys < nsys; isys++) {
    hSysThrown[isys]->Write();
  }
  hNsk->Write();

  // Find min/max and calculate fractional deviation from nominal
#line 0
  hErecoMax = new TH1D(*hErecoNominal);
  hErecoMin = new TH1D(*hErecoNominal);
  hErecoMax->SetDirectory(0);
  hErecoMin->SetDirectory(0);
  for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
    double nom   = hErecoNominal->GetBinContent(ibin);
    double sumsq = 0;
    for(int iexp = 0; iexp < nexp; iexp++) {
       double N = hEreco[iexp]->GetBinContent(ibin);
       sumsq += TMath::Power(N-nom,2);
    }//iexp
    double sig = TMath::Sqrt(sumsq/nexp);
    double min = nom-sig;
    double max = nom+sig;
    if(nom>0){
      hErecoMax -> SetBinContent (ibin,  (max-nom)/nom);
      hErecoMin -> SetBinContent (ibin,  (min-nom)/nom);
    }
  }//ibin
  hErecoMax->GetYaxis()->SetRangeUser(-2,2);
  hErecoMin->GetYaxis()->SetRangeUser(-2,2);
  hErecoMax->Write(Form("%s_pluserr",modetitle.c_str()));
  hErecoMin->Write(Form("%s_minuserr",modetitle.c_str()));

  f.Close();
}
