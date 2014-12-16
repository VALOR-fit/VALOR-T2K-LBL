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
//   7 tweak BANFF flux only
//   8 tweak BANFF xsec only
//   9 tweak SKdet + FSI
//  10 tweak BANFF + unconstrained xsec
//  11 oscillation parameters
//
// 999 stat only
//
// 100 tweak banff00
// 101 tweak banff01
// ...
// 100+kNumOfParamsInBANFFfit tweak last banff parameter
//
// 200 tweak pfermi
// 201 tweak specfunc
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
// 900 tweak skdetfsi00
// 901 tweak skdetfsi01
// ...
// 900+kNumOfSKDetFSIParams tweak last skdetfsi parameter
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

#include<typeinfo.h>

TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}

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
  const char * mixparamchoice     = "sinsqtheta",
  const char * hierarchy          = "normal"
)
{
  using namespace std;
  using namespace valor;
  
  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

  //Check if the mode is valid
  string modetitle = "";
  switch(mode) {
    case 0:
      modetitle = "all"; break;
    case 1:
      modetitle = "BANFF"; break;
    case 2:
      modetitle = "unconstrained_xsec"; break;
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
    case 11:
      modetitle = "Osc_param"; break;

    case 200:
      modetitle = "pfermi"; break;
    case 201:
      modetitle = "specfunc"; break;
    case 202:
      modetitle = "wshape"; break;
    case 203:
      modetitle = "shapeCoth"; break;
    case 204:
      modetitle = "pilessdelta"; break;
    case 209:
      modetitle = "norm_cccoh"; break;
    case 210:
      modetitle = "norm_nc1piC"; break;
    case 212:
      modetitle = "norm_ncoth"; break;
    case 213:
      modetitle = "norm_nuetonumu"; break;
    case 214:
      modetitle = "norm_antinutonu"; break;
    case 215:
      modetitle = "bindE"; break;
      
    case 700:
      modetitle = "mixing13"; break;
    case 701:
      modetitle = "mixing12"; break;
    case 702:
      modetitle = "dm12sq"; break;
    case 703:
      modetitle = "deltacp"; break;

    case 999:
      modetitle = "statonly"; break;

    default:
      if(mode >= 100 && mode < 100+conv::kNumOfParamsInBANFFfit) {
        modetitle = Form("banff%02d",mode-100); break;
      }
      else if(mode >= 900 && mode < 900+conv::kNumOfSKDetFSIParams) {
        modetitle = Form("skdetfsi%03d",mode-900); break;
      }
      cout << "Unknown mode " << mode << endl; exit(1); break;
  }//switch mode

  TString config_str(config_file);
  if(config_str.Contains("prefit"))
    modetitle += "prefit";

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  T2KOpt::Instance()->ReadXmlPath(config);
    
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
  }
  else if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }

  //return if tweaking osc params & there are no oscillations
  if(mode >= 700 && mode <= 799
     && input_mixing23 < conv::kASmallNum 
     && input_mixing13 < conv::kASmallNum
     && input_mixing12 < conv::kASmallNum
     && input_deltacp  < conv::kASmallNum
     && input_dm12sq   < conv::kASmallNum
     && input_dm23sq   < conv::kASmallNum)
    return;

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Physics hypothesis: \n" << ph << endl;
  cout << *Opt::Instance() << endl;

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);
  //
  SuperkObservation skobs;
  
  for(int is = 0; is < conv::kns; is++) 
  {
    conv::SuperkSample_t s = conv::kSKSamples[is];
    if(s == conv::kSKsUndefined) continue;
    
    string sample_name = util::AsStringCompact(s);
    const char * sname = sample_name.c_str();
  
    // Output histograms
    TH1D * hErecoNominal = 0;
    TH1D * hEreco[nexp] = {0};
    TH1D * hNsk = new TH1D(Form("%s_nsk",modetitle.c_str()),Form("%s_%s_nsk",sname,modetitle.c_str()),10000,0,1000);
    int nbanff = conv::kNumOfParamsInBANFFfit;
    
    int nfsi = conv::kNumOfSKDetFSIParams;
    const int nsys = nbanff + nfsi + 8 + 13;
    TH1D * hSysThrown[nsys];
    for(int isys = 0; isys < nsys; isys++) {
      string sysstr = "";
      if(isys < nbanff)
        sysstr = Form("banff%02i",isys);
      else if(isys < nbanff + nfsi)
        sysstr = Form("fsi_%03i",isys - nbanff);
      else {
        int isys2 = isys - nbanff - nfsi;
        switch(isys2) {
          case(0) :
            sysstr = "sk_e_scale"; break;
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
            sysstr = "pfermi"; break;
          case(9) :
            sysstr = "specfunc"; break;
          case(10) :
            sysstr = "wshape"; break;
          case(11) :
            sysstr = "shapeCCoth"; break;
          case(12) :
            sysstr = "pilessdelta"; break;
          case(14) :
            sysstr = "norm_cccoh"; break;
          case(15) :
            sysstr = "norm_nc1piC"; break;
          case(17) :
            sysstr = "norm_ncoth"; break;
          case(18) :
            sysstr = "norm_nuetonumu"; break;
          case(19) :
            sysstr = "norm_antinutonu"; break;
          case(20) :
            sysstr = "bindE"; break;
          default:
            sysstr = ""; break;
        }//switch(isys2)
      }//isys>=nbanff+nfsi
      string hname(Form("%s_%s_tweak_throws_%s",sname,modetitle.c_str(),sysstr.c_str()));
      hSysThrown[isys] = new TH1D(hname.c_str(),hname.c_str(),1000,-5,5);
      hSysThrown[isys]->SetDirectory(0);
    }

    // Open output files
    string rootmode = "update";
    if(!update_rootfile)
      rootmode = "recreate";
    string oscparams = util::AsString(ph.Hypothesis(),
       input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
    TFile f(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s.root",sname,oscparams.c_str()),rootmode.c_str());

    if(!update_rootfile)
      ofstream outfile(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s.txt",sname,oscparams.c_str()),ios_base::trunc);
    else
      ofstream outfile(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s.txt",sname,oscparams.c_str()),ios_base::app);

    //remove the [] so can use it in SaveAs
    while(oscparams.find_first_of("[]") != string::npos) {
      oscparams.erase(oscparams.find_first_of("[]"), 1);
    }

    // Get the nominal
    skexp.Generate(ph);
    hErecoNominal = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
    for(int ir = 1; ir <= hErecoNominal->GetNbinsX(); ir++) {
      hErecoNominal->SetBinContent(ir, skexp.CurrPredictedErecoSpectrumArray(s)[ir-1]);
    }
    hErecoNominal->SetDirectory(0);
    hErecoNominal->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hErecoNominal->GetYaxis()->SetTitle("Events per bin");
    hErecoNominal->SetLineWidth(2);
    if(!update_rootfile)
      hErecoNominal->Write("nominal");
      
    TH1D * hErecoNominal_per50MeV = spectrum_to_per_50mev(hErecoNominal);
    hErecoNominal_per50MeV->Write("nominal_per_50MeV");

    TStopwatch w;
    w.Start();

    // Tweak!
    int iskip = 0;
    for(int iexp = 0; iexp < nexp; iexp++) {

      if((iexp % 1000) == 0)
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
        
        //BANFF (flux only)
        if(mode == 7) {
          TRandom3 rand(0);
          for(int ib = 0; ib < conv::kNumOfFluxParamsInBANFFfit; ib++) {
            skexp.TweakSystematic(Form("banff%02d",ib),rand.Gaus(0.,1.));
          }
        }

        //BANFF (xsec only)
        if(mode == 8) {
          TRandom3 rand(0);
          skexp.TweakSystematic("banffmaqe",rand.Gaus(0.,1.));
          skexp.TweakSystematic("banffmares",rand.Gaus(0.,1.));
          skexp.TweakSystematic("banffccqeE0",rand.Gaus(0.,1.));
          skexp.TweakSystematic("banffccqeE1",rand.Gaus(0.,1.));
          skexp.TweakSystematic("banffccqeE2",rand.Gaus(0.,1.));
          skexp.TweakSystematic("banffcc1piE0",rand.Gaus(0.,1.));
          skexp.TweakSystematic("banffcc1piE1",rand.Gaus(0.,1.));
          skexp.TweakSystematic("banffmecnorm",rand.Gaus(0.,1.));
          skexp.TweakSystematic("banffpdd",rand.Gaus(0.,1.));
        }

        //SK det + FSI
        if(mode == 9) {
          skexp.TweakAllSKDetectorFSISystematics();
        }

       //unconstrained xsec + BANFF
        if(mode == 10) {
          skexp.TweakAllUncorrelatedXSecSystematics();
          skexp.TweakAllBANFFSystematics();
        }

        //non-23 oscillation parameters
        if(mode == 11) {
          ph.SetCurrValue(conv::kPhypMixing12,  util::RandomizeMixing12());
          ph.SetCurrValue(conv::kPhypMixing13,  util::RandomizeMixing13());
          ph.SetCurrValue(conv::kPhypCPV    ,   util::RandomizeCPV());
          ph.SetCurrValue(conv::kPhypDm12sq,    util::RandomizeDm12sq());
        }

      }//mode<100

      else if(mode < 999) {//everything else but the stat
        
        TRandom3 rand(0);
        double tweak = rand.Gaus(0.,1.);

        switch(mode) {
          case 200:
            skexp.TweakSystematic("pfermi",tweak); break;
          case 201:
            skexp.TweakSystematic("specfunc",tweak); break;
          case 202:
            skexp.TweakSystematic("wshape",tweak); break;
          case 203:
            skexp.TweakSystematic("shapeCCoth",tweak); break;
          case 204:
            skexp.TweakSystematic("pilessdelta",tweak); break;
          case 209:
            skexp.TweakSystematic("norm_cccoh",tweak); break;
          case 210:
            skexp.TweakSystematic("norm_nc1piC",tweak); break;
          case 212:
            skexp.TweakSystematic("norm_ncoth",tweak); break;
          case 213:
            skexp.TweakSystematic("norm_nuetonumu",tweak); break;
          case 214:
            skexp.TweakSystematic("norm_antinutonu",tweak); break;
          case 215:
            skexp.TweakSystematic("bindE",tweak); break;

          case 300:
            skexp.TweakSystematic("sk_e_scale",tweak); break;
      
          case 700:
            ph.SetCurrValue(conv::kPhypMixing13,  util::RandomizeMixing13()); break;
          case 701:
            ph.SetCurrValue(conv::kPhypMixing12,  util::RandomizeMixing12()); break;
          case 702:
            ph.SetCurrValue(conv::kPhypDm12sq,    util::RandomizeDm12sq()); break;
          case 703:
            ph.SetCurrValue(conv::kPhypCPV    ,   util::RandomizeCPV()); break;

          default:
            if(mode >= 100 && mode < 100+conv::kNumOfParamsInBANFFfit) {
              TVectorD v(conv::kNumOfParamsInBANFFfit);
              v = skexp.GetBANFF()->BestFitValues();
              int ib = mode - 100;
              v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
              skexp.GetBANFF()->Tweak(true,&v,1);
              break;
            }
            else if(mode >= 900 && mode < 900+conv::kNumOfSKDetFSIParams) {
              TVectorD v(conv::kNumOfSKDetFSIParams);
              for (int ipf = 0; ipf < conv::kNumOfSKDetFSIParams; ipf++) {
                v[ipf] = 1;
              }
              int ib = mode - 900;
              v[ib] = 1 + (tweak * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()[ib][ib]));
              v.Print();
              skexp.GetSuperkDetectorFSI()->Tweak(true,&v);
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

      if(mode == 999) {
        bool fluctuate = true;
        skobs.CreateMockData(skexp,ph,fluctuate);
      }

      //store spectrum
      hEreco[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      for(int ir = 1; ir <= hEreco[iexp]->GetNbinsX(); ir++) {
        if(mode != 999)
          hEreco[iexp]->SetBinContent(ir, skexp.CurrPredictedErecoSpectrumArray(s)[ir-1]);
        else
          hEreco[iexp]->SetBinContent(ir, skobs.ErecoSpectrumArray(s)[ir-1]);
      }

      hEreco[iexp]->SetDirectory(0);
      hEreco[iexp]->SetTitle(Form("%s_exp%d",modetitle.c_str(),iexp));
      if(save_toy_expts)
        hEreco[iexp]->Write(Form("%s_exp%d",modetitle.c_str(),iexp));                    

      //store Nsk
      double nskint = 0;
      for(int ir = 1; ir <= util::NumErecoBins(s); ir++)
        if(mode != 999)
          nskint += skexp.CurrPredictedErecoSpectrumArray(s)[ir-1];
        else
          nskint += skobs.ErecoSpectrumArray(s)[ir-1];
      hNsk->Fill(nskint);

      //store thrown parameter values
      for(int isys = 0; isys < nbanff; isys++) {
        hSysThrown[isys]->Fill(skexp.GetBANFF()->CurrentValue(conv::kNuisBANFF00 + isys));
      }
      //for(int isys = 0; isys < nfsi; isys++) {
        //hSysThrown[isys+nbanff]->Fill(skexp.GetFSI()->CurrentValue(util::PhysParamIdFromString(Form("fsi_%03d",isys))));
      for(int isys = 0; isys < nfsi-1; isys++) {
        hSysThrown[isys+nbanff]->Fill(skexp.GetSuperkDetectorFSI()->CurrentValue(util::PhysParamIdFromString(Form("skdetfsi%03d",isys))));
      }
      hSysThrown[nbanff+nfsi ]->Fill(skexp.GetSuperkDetectorFSI()->CurrentValue(util::PhysParamIdFromString("sk_e_scale")));
      
      hSysThrown[nbanff+nfsi+8 ]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("pfermi")));
      hSysThrown[nbanff+nfsi+9 ]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("specfunc")));
      hSysThrown[nbanff+nfsi+10]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("wshape")));
      hSysThrown[nbanff+nfsi+11]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("shapeCCoth")));
      hSysThrown[nbanff+nfsi+12]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("pilessdelta")));
      hSysThrown[nbanff+nfsi+13]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("norm_cccoh")));
      hSysThrown[nbanff+nfsi+14]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("norm_nc1piC")));
      hSysThrown[nbanff+nfsi+15]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("norm_ncoth")));
      hSysThrown[nbanff+nfsi+16]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("norm_nuetonumu")));
      hSysThrown[nbanff+nfsi+17]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("norm_antinutonu")));
      hSysThrown[nbanff+nfsi+18]->Fill(skexp.GetNuXSec()->CurrentValue(util::PhysParamIdFromString("bindE")));
    }//end loop over experiments
    w.Stop();
    cout << iskip << " experiments skipped due to min events" << endl;

    //save some histograms
    for(int isys = 0; isys < nsys; isys++) {
      hSysThrown[isys]->Write();
    }
    hNsk->Write();

    // Find min/max and calculate fractional deviation from nominal
    hErecoMax = new TH1D(*hErecoNominal);
    hErecoMin = new TH1D(*hErecoNominal);
    hErecoMax->SetDirectory(0);
    hErecoMin->SetDirectory(0);
    hErecoMax->GetYaxis()->SetTitle("Ratio of events at 1 #sigma to nominal");
    hErecoMin->GetYaxis()->SetTitle("Ratio of events at 1 #sigma to nominal");
    hErecoMaxEnv = new TH1D(*hErecoNominal);
    hErecoMinEnv = new TH1D(*hErecoNominal);
    hErecoMaxEnv->SetDirectory(0);
    hErecoMinEnv->SetDirectory(0);
    for(int ibin=1; ibin <= util::NumErecoBins(s); ibin++) {
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

    //Draw & save
    TCanvas * c = new TCanvas();
    hNsk->Draw();
    c->SaveAs(Form("%s_%s_nsk_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
    hErecoMax->Draw();
    c->SaveAs(Form("%s_%s_pluserr_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
    hErecoMin->Draw();
    c->SaveAs(Form("%s_%s_minuserr_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
    hErecoMaxEnv->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
    hErecoMaxEnv->SetLineStyle(2);
    hErecoMinEnv->SetLineStyle(2);
    hErecoMaxEnv->SetLineColor(kBlue+3);
    hErecoMinEnv->SetLineColor(kBlue+3);
    hErecoMaxEnv->Draw();
    hErecoMinEnv->Draw("same");
    hErecoNominal->Draw("same");
    c->SaveAs(Form("%s_%s_envelope_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
    
    TH1D * hErecoMaxEnv_per50MeV = spectrum_to_per_50mev(hErecoMaxEnv);
    TH1D * hErecoMinEnv_per50MeV = spectrum_to_per_50mev(hErecoMinEnv);
    
    TCanvas * c = new TCanvas();
    hErecoMaxEnv_per50MeV->Draw();
    hErecoMinEnv_per50MeV->Draw("same");
    hErecoNominal_per50MeV->Draw("same");
    c->SaveAs(Form("%s_%s_envelope_per_50MeV_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
    

    TCanvas * csys = new TCanvas();
    csys->Divide(ceil(sqrt(nsys)),ceil(sqrt(nsys)));
    for(int isys=0; isys < nsys; isys++) {
      csys->cd(isys+1);
      hSysThrown[isys]->Draw();
      cout << isys << "\t" << hSysThrown[isys]->GetTitle() << endl;
    }
    csys->SaveAs(Form("%s_%s_throws_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
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
  } //samples
  cout << "Tweaks performed in " 
       << w.CpuTime() << " sec." << endl;
         
}
