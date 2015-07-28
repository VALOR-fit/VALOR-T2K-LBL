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
//   5 tweak only BANFF unstrained parameters
//   6 tweak all BANFF constrained parameters
//   7 tweak BANFF flux only
//   8 tweak BANFF xsec only
//   9 tweak SKdet + FSI
//  10 tweak BANFF + unconstrained xsec
//  11 oscillation parameters
//
// 200 tweak pfermi
// 201 tweak specfunc
// 202 tweak xsec_wshape
// 203 tweak xsec_shape_ccoth
// 204 tweak xsec_pilessdelta
// 205 tweak xsec_norm_cccoh
// 206 tweak xsec_norm_nc1piC
// 207 tweak xsec_norm_ncoth
// 208 tweak xsec_norm_nuetonumu
// 209 tweak xsec_norm_antinutonu
// 210 tweak xsec_bindE
// 211 tweak xsec_rpa (0 or 1)
//
// 400 tweak banff00
// 401 tweak banff01
// 400+kNumOfParamsInBANFFfit tweak last banff parameter
//
// 500 tweak banff00_rhc
// 501 tweak banff01_rhc
// 500+kNumOfRHCParamsInBANFFfit tweak last banff parameter
//
// 800 tweak skdetfsi00
// 801 tweak skdetfsi01
// 800+kNumOfSKDetFSIParams tweak last skdetfsi parameter
//
// 900 tweak skdetfsi00_rhc
// 901 tweak skdetfsi01_rhc
// 900+kNumOfRHCSKDetFSIParams tweak last skdetfsi parameter
//
// 999 stat only
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

//#include<typeinfo.h>



#ifndef __CINT__
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TAxis.h"
#include "TXMLEngine.h"

#include "Conventions/Opt.h"
#include "Conventions.h"
#include "Utilities.h"
#include "Utilities/Rnd.h"
#include "Utilities/Env.h"
#include "FitCore/PhysHypothesis.h"
#include "SuperkExpectation.h"
#include "SuperkObservation.h"
#include "FitCore/FitInputs.h"
#include "FitCore/FitterMinuit.h"
#include "FitCore/UtlFit.h"
#include "BANFF.h"
#include "SuperkDetectorFSI.h"
#include "NuXSec.h"
#include "UtlFitParamListFactory.h"
#include "T2KPhysicalParameterLimits.h"
#endif


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
  int          mode                     = 0,                                         /* options: see above */
  bool         update_rootfile          = false,                                        /* "update" ROOT output file? If false, will "recreate" */
  const char * sk_sample                = "RHC1Rmu",                                   /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
  const char * hypothesis               = "3f_osc",                                    /* options: "2f_numu_disap", "3f_osc", "3f_oscbar" */
  double       input_mixing12           = -99,                                         /* */
  double       input_mixing13           = -99,                                         /* */
  double       input_mixing23           = -99,                                         /* */
  double       input_deltacp            = -99,                                         /* */
  double       input_dm12sq             = -99,                                         /* */
  double       input_dm23sq             = -99,                                         /* */
  bool         save_toy_expts           = false,                                       /* save the ereco spectra? */
  bool         save_throws              = true,                                        /* save the thrown parameter values*/
  const char * config_file              = "analysis_configuration_run1to6_test.xml",   /* */
  const int    nexp                     = 50,                                          /* number of experiments to make */
  const char * mixparamchoice           = "sinsqtheta",
  const char * hierarchy                = "normal",
  bool         reject_unphysical_throws = true    ,
  bool         use_new_definition_of_dn_over_n = true
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  bool RHCFlag = false;
  TString str_sk_sample(sk_sample);
  if(str_sk_sample.Contains("RHC"))
    RHCFlag = true;

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  string sample_name = util::AsStringCompact(s);
  const char * sname = sample_name.c_str();
  
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
      std::cerr<<"FSI only no longer implemented"<<std::endl;
      exit(-1);
      //~ char fsibins[40] = "";
      //~ sprintf(fsibins, "fsi_%d", conv::kNumOfFSIParams);
      //~ modetitle = fsibins; break;
    case 5:
      modetitle = "BANFF_unconstrained"; break;
    case 6:
      modetitle = "BANFF_constrained"; break;
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
    case 205:
      modetitle = "norm_cccoh"; break;
    case 206:
      modetitle = "norm_nc1piC"; break;
    case 207:
      modetitle = "norm_ncoth";
      break;
    case 208:
      modetitle = "norm_nuetonumu"; break;
    case 209:
      modetitle = "norm_antinutonu"; break;
    case 210:
      modetitle = "bindE"; break;

    case 211:
      modetitle = "rpa"; break;
      
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

      // banff fhc
      if(mode >= 400 && mode < 400+conv::kNumOfParamsInBANFFfit) {
        modetitle = Form("banff%02d",mode-400); break;
      }
      // banff rhc
      if(mode >= 500 && mode < 500+conv::kNumOfRHCParamsInBANFFfit) {
        modetitle = Form("banff%02d_rhc",mode-500); break;
      }
      // skdetfsi fhc
      else if(mode >= 800 && mode < 800+conv::kNumOfSKDetFSIParams) {
        modetitle = Form("skdetfsi%03d",mode-800); break;
      }
      // skdetfsi rhc
      else if(mode >= 900 && mode < 900+conv::kNumOfRHCSKDetFSIParams) {
        modetitle = Form("skdetfsi%03d_rhc",mode-900); break;
      }
      cout << "Unknown mode " << mode << endl; exit(1); break;
  }//switch mode
  
  int syst_nsigma = util::GetSystNSigma();

  TString config_str(config_file);
  bool isprefit = config_str.Contains("prefit");
  if(isprefit)
    modetitle += "prefit";

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  T2KOpt::Instance()->ReadXmlPath(config);
    
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
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
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
  
      
  string oscparams = util::AsString(ph.Hypothesis(),
                       input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  
  
  // Output histograms
  TH1D * hErecoNominal = 0;
  VALORArray * obsEreco[nexp];
  for (int i = 0 ; i < nexp; i++) {obsEreco[i]=0;}
  TH1D * hNsk = new TH1D(Form("%s_%s_nsk",sname,modetitle.c_str()),Form("%s_%s_nsk",sname,modetitle.c_str()),10000,0,1000);
  
  
  int nbanff;
  if(RHCFlag) nbanff = conv::kNumOfRHCParamsInBANFFfit;
  else        nbanff = conv::kNumOfParamsInBANFFfit;
  
  int nfsi;
  if(RHCFlag) nfsi = conv::kNumOfRHCSKDetFSIParams;
  else        nfsi = conv::kNumOfRHCSKDetFSIParams;

  int nxsecuncorr;
  if(RHCFlag) nxsecuncorr = conv::kNumOfRHCXSecSyst;
  else        nxsecuncorr = conv::kNumOfXSecSyst;

  const int sys = nbanff + nfsi + nxsecuncorr;

  FitParamList * fparm = util::NamedParameterList("nuosc3f_12fix+nuisance"); 
  const int np = fparm->GetNumberParameters();
  TH1D * hSysThrown[np];
  
  if (save_throws) {
    for(int ip = 0; ip < np; ip++) {
      conv::PhysParamId_t p = fparm->WhichParam(ip);
      //if(!fparm->IsNuisance(p)) continue; 
      
      string sysstr = util::AsStringCompact(p);
      //cout << sysstr << endl;
      string hname(Form("%s_%s_tweak_throws_%s",sname,modetitle.c_str(),sysstr.c_str()));
      hSysThrown[ip] = new TH1D(hname.c_str(),hname.c_str(),1000,syst_nsigma,syst_nsigma);
      hSysThrown[ip]->SetDirectory(0);
    }
  }

  
  // Open output files
  string rootmode = "update";
  if(!update_rootfile)
    rootmode = "recreate";

  TFile f(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s.root",sname,oscparams.c_str()),rootmode.c_str());
  
  ofstream * outfilep = 0;
  
  if(!update_rootfile){
    outfilep = new ofstream(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s.txt",sname,oscparams.c_str()),ios_base::trunc);
  }
  else {
    outfilep = new ofstream(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s.txt",sname,oscparams.c_str()),ios_base::app);
  }
  
  ofstream & outfile = *outfilep;
  
  //remove the [] so can use it in SaveAs
  while(oscparams.find_first_of("[]") != string::npos) {
    oscparams.erase(oscparams.find_first_of("[]"), 1);
  }

  // Get the nominal
  skexp.Generate(ph);
  hErecoNominal = skexp.CurrPredictedErecoSpectrum(s);
  hErecoNominal->SetDirectory(0);
  hErecoNominal->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoNominal->GetYaxis()->SetTitle("Events per bin");
  hErecoNominal->SetLineWidth(2);
  
  if (modetitle.find("prefit")==string::npos) {
    hErecoNominal->Write(Form("%s_nominal",sname));
  }
  else {
    hErecoNominal->Write(Form("%s_nominalprefit",sname));
  }
  
  double nominal_total_nsk = hErecoNominal->Integral();
  
  TH1D * hErecoNominal_per50MeV = spectrum_to_per_50mev(hErecoNominal);
  hErecoNominal_per50MeV->Write(Form("%s_nominal_per_50MeV",sname));

  TStopwatch w;
  w.Start();
  TStopwatch w2;
  w2.Start();
  
  bool done_print = true;
  
  int n_ereco_bins = util::NumErecoBins(s);
  
  double * sumsq_per_bin = (double*)calloc(n_ereco_bins,sizeof(double));
  VALORArray * nsk_per_toy = new VALORArray(nexp);
  double * nsk_per_toy_arr = nsk_per_toy.GetArray();
  
  TRandom3 rand(0);
  
  PhysParamList pparm;
  
  // Tweak!
  int iskip = 0;
  for(int iexp = 0; iexp < nexp; iexp++) {

    if((iexp % 100) == 0 && (!done_print)) {
      done_print = true;
      cout << endl << "Experiment " << iexp << " for systematics: " << modetitle << endl;
      cout << "Last 100 experiments took "<<w2.RealTime()<<" seconds"<<endl;
      w2.Start();
    }
    
    if(mode < 100) {
      
      // All
      if(mode == 0) {
        skexp.TweakAllSystematics();
      }

      // BANFF
      if(mode == 1) {
        skexp.TweakAllSampleBANFFSystematics();
      }

      // Uncorrelated XSec
      if(mode == 2) {
        skexp.TweakAllUncorrelatedXSecSystematics();
        skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
      }
      
      //BANFF (unconstrained only)
      if(mode == 5) {
        skexp.GetBANFF()->TweakSampleConstrainedOnly(false,0,true); // cholesky w/ full unconstrained matrix
      }
      //BANFF (constrained only)
      if(mode == 6) {
        skexp.GetBANFF()->TweakSampleConstrainedOnly(false,0, false); // cholesky w/ full cov matrix
      }
      
      //BANFF (flux only)
      if(mode == 7) {
        skexp.GetBANFF()->TweakSampleFluxOnly(false,0); // cholesky w/ full cov matrix
      }

      //BANFF (xsec only)
      if(mode == 8) {
        skexp.GetBANFF()->TweakSampleXSecOnly(false,0); // cholesky w/ full cov matrix
      }

      //SK det + FSI
      if(mode == 9) {
        skexp.TweakAllSampleSKDetectorFSISystematics();
      }

     //unconstrained xsec + BANFF
      if(mode == 10) {
        skexp.TweakAllUncorrelatedXSecSystematics();
        skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
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
      
      double tweak = rand.Gaus(0.,1.);

      switch(mode) {
        case 200:
          if( !util::IsXSec(util::PhysParamIdFromString("pfermi")) ) break; 
          skexp.TweakSystematic("pfermi",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 201:
          if( !util::IsXSec(util::PhysParamIdFromString("specfunc")) ) break; 
          skexp.TweakSystematic("specfunc",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 202:
          if( !util::IsXSec(util::PhysParamIdFromString("wshape")) ) break; 
          skexp.TweakSystematic("wshape",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 203:
          if( !util::IsXSec(util::PhysParamIdFromString("shapeCCoth")) ) break; 
          skexp.TweakSystematic("shapeCCoth",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 204:
          if( !util::IsXSec(util::PhysParamIdFromString("pilessdelta")) ) break; 
          skexp.TweakSystematic("pilessdelta",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 205:
          if( !util::IsXSec(util::PhysParamIdFromString("norm_cccoh")) ) break; 
          skexp.TweakSystematic("norm_cccoh",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 206:
          if( !util::IsXSec(util::PhysParamIdFromString("norm_nc1piC")) ) break; 
          skexp.TweakSystematic("norm_nc1piC",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 207:
          if( !util::IsXSec(util::PhysParamIdFromString("norm_ncoth")) ) break; 
          skexp.TweakSystematic("norm_ncoth",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 208:
          if( !util::IsXSec(util::PhysParamIdFromString("norm_nuetonumu")) ) break; 
          skexp.TweakSystematic("norm_nuetonumu",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 209:
          if( !util::IsXSec(util::PhysParamIdFromString("norm_antinutonu")) ) break; 
          skexp.TweakSystematic("norm_antinutonu",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;
        case 210:
          if( !util::IsXSec(util::PhysParamIdFromString("bindE")) ) break; 
          skexp.TweakSystematic("bindE",tweak); 
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;

        case 211:
          if( !util::IsXSec(util::PhysParamIdFromString("rpa")) ) break; 
          tweak = rand.Uniform(0.,1.);
          skexp.TweakSystematic("rpa",tweak); // Automatically set at 0 or 1 in the Tweak
          skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
          break;

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
          
          // banff fhc
          if(mode >= 400 && mode < 400+conv::kNumOfParamsInBANFFfit) {
            TVectorD v(conv::kNumOfParamsInBANFFfit);
            v = skexp.GetBANFF()->BestFitValues();
            int ib = mode - 400;
            v[ib] = skexp.GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->CovMatrix()[ib][ib]));
            skexp.GetBANFF()->TweakSample(true,&v);
            break;
          }
          // banff rhc
          else if(mode >= 500 && mode < 500+conv::kNumOfRHCParamsInBANFFfit) {
            TVectorD v(conv::kNumOfRHCParamsInBANFFfit);
            v = skexp.GetBANFF()->RHCBestFitValues();
            int ib = mode - 500;
            v[ib] = skexp.GetBANFF()->RHCBestFitValues()[ib] + (tweak * sqrt(skexp.GetBANFF()->RHCCovMatrix()[ib][ib]));
            skexp.GetBANFF()->TweakSample(true,&v);
            break;
          }
          // skdetfsi fhc
          else if(mode >= 800 && mode < 800+conv::kNumOfSKDetFSIParams) {
            TVectorD v(conv::kNumOfSKDetFSIParams);
            for (int ipf = 0; ipf < conv::kNumOfSKDetFSIParams; ipf++) {
              v[ipf] = 1;
            }
            int ib = mode - 800;
            v[ib] = 1 + (tweak * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()[ib][ib]));
            skexp.GetSuperkDetectorFSI()->Tweak(true,&v);
            break;
          }
          // skdetfsi rhc
          else if(mode >= 900 && mode < 900+conv::kNumOfRHCSKDetFSIParams) {
            TVectorD v(conv::kNumOfRHCSKDetFSIParams);
            for (int ipf = 0; ipf < conv::kNumOfRHCSKDetFSIParams; ipf++) {
              v[ipf] = 1;
            }
            int ib = mode - 900;
            v[ib] = 1 + (tweak * sqrt(skexp.GetSuperkDetectorFSI()->RHCCovMatrix()[ib][ib]));
            skexp.GetSuperkDetectorFSI()->TweakRHC(true,&v);
            break;
          }

          cout << "Unknown mode " << mode << endl; exit(1); break;
      }//switch(mode)
    }//mode>=100
    
    skexp.CopyTweaksTo(pparm);
    
    if (reject_unphysical_throws) {
      bool within_limits = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(pparm);
      if (!within_limits) {
        //cout<<"Systematic thrown outside physical limits, regenerating"<<std::endl;
        iskip++;
        iexp--;
        continue;
      }
      else {done_print = false;}
    }
    else {done_print = false;}
    
    
    // generate predictions
    bool applied_negevents_condition = skexp.Generate(ph);
    if(applied_negevents_condition) {
      cout << "skipping" << endl;
      iskip++;
      iexp--;
      continue;
    }

    if(mode == 999) {
      bool fluctuate = true;
      skobs.CreateMockData(skexp,ph,fluctuate);
    }
    
    VALORArray * spectrum = 0;

    //store spectrum
    if (save_toy_expts) {
      if(mode == 999 ) {
        spectrum = new VALORArray(*skobs.ErecoSpectrumVArray(s));
      }
      else {
        spectrum = new VALORArray(*skexp.CurrPredictedErecoVArray(s));
      }
      
      spectrum->SetTitle(Form("%s_%s_exp%d",sname,modetitle.c_str(),iexp));
      spectrum->Write   (Form("%s_%s_exp%d",sname,modetitle.c_str(),iexp));
    }
    else {
      if(mode == 999 ) {
        spectrum = skobs.ErecoSpectrumVArray(s);
      }
      else {
        spectrum = skexp.CurrPredictedErecoVArray(s);
      }
    }

    //store Nsk
    double nskint = spectrum->Integral();
    hNsk->Fill(nskint);
    nsk_per_toy_arr[iexp] = nskint;
    
    // Store tweaks.
    for (int ibin = 0 ; ibin < n_ereco_bins ; ibin++) {
      double nom   = hErecoNominal->GetBinContent(ibin+1);
      double N = spectrum->GetArray()[ibin];
      sumsq_per_bin[ibin] += (N-nom)*(N-nom);
    }
    
    if (save_toy_expts) { delete spectrum; }
    
    //store thrown parameter values
    
    if (save_throws) {
      for(int ip = 0; ip < np; ip++) {
        conv::PhysParamId_t p = fparm->WhichParam(ip);
        //if(!fparm->IsNuisance(p)) continue; 
        
        string sysstr = util::AsStringCompact(p);
        
        if(util::IsBANFF(p)){
          if(RHCFlag)
            hSysThrown[ip]->Fill(skexp.GetBANFF()->RHCCurrentValue(p));
          else 
            hSysThrown[ip]->Fill(skexp.GetBANFF()->CurrentValue(p));
        }
        else if(util::IsXSec(p)){
          hSysThrown[ip]->Fill(skexp.GetNuXSec()->CurrentValue(p));
        }
        if(util::IsSKDetFSI(p)){
          if(RHCFlag)
            hSysThrown[ip]->Fill(skexp.GetSuperkDetectorFSI()->RHCCurrentValue(p));
          else 
            hSysThrown[ip]->Fill(skexp.GetSuperkDetectorFSI()->CurrentValue(p));
        }
      }
    }// save throws
    
  }//end loop over experiments
  
  w.Stop();
  cout << iskip << " experiments skipped due to limits" << endl;

  //save some histograms
  if (save_throws) {
    for(int ip = 0; ip < np; ip++) {
      conv::PhysParamId_t p = fparm->WhichParam(ip);
      //if(!fparm->IsNuisance(p)) continue; 
      hSysThrown[ip]->Write();
    }
  }
  hNsk->Write();


  // Find min/max and calculate fractional deviation from nominal
  TH1D * hErecoMax = new TH1D(*hErecoNominal);
  TH1D * hErecoMin = new TH1D(*hErecoNominal);
  hErecoMax->SetDirectory(0);
  hErecoMin->SetDirectory(0);
  hErecoMax->GetYaxis()->SetTitle("Ratio of events at 1 #sigma to nominal");
  hErecoMin->GetYaxis()->SetTitle("Ratio of events at 1 #sigma to nominal");
  TH1D * hErecoMaxEnv = new TH1D(*hErecoNominal);
  TH1D * hErecoMinEnv = new TH1D(*hErecoNominal);
  hErecoMaxEnv->SetDirectory(0);
  hErecoMinEnv->SetDirectory(0);
  for(int ibin=1; ibin <= util::NumErecoBins(s); ibin++) {
    double nom   = hErecoNominal->GetBinContent(ibin);
    double sumsq = 0;
    double sig = TMath::Sqrt(sumsq_per_bin[ibin-1]/nexp);
    double min = nom-sig;
    double max = nom+sig;
    if(nom>0){
      hErecoMax -> SetBinContent (ibin,  (max-nom)/nom);
      hErecoMin -> SetBinContent (ibin,  (min-nom)/nom);
      hErecoMaxEnv -> SetBinContent (ibin,  max);
      hErecoMinEnv -> SetBinContent (ibin,  min);
    }
  }//ibin
  hErecoMax   ->Write(Form("%s_%s_pluserr" ,sname,modetitle.c_str()));
  hErecoMin   ->Write(Form("%s_%s_minuserr",sname,modetitle.c_str()));
  hErecoMaxEnv->Write(Form("%s_%s_envpluserr" ,sname,modetitle.c_str()));
  hErecoMinEnv->Write(Form("%s_%s_envminuserr",sname,modetitle.c_str()));

  //Draw & save
  TCanvas * c = new TCanvas();
  hNsk->Draw();
  c->SaveAs(Form("%s_%s_nsk_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
  hErecoMax->Draw();
  c->SaveAs(Form("%s_%s_pluserr_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
  hErecoMin->Draw();
  c->SaveAs(Form("%s_%s_minuserr_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
  hErecoMaxEnv->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
  hErecoMaxEnv->SetLineStyle(3);
  hErecoMinEnv->SetLineStyle(3);
  hErecoMaxEnv->SetLineColor(kRed);
  hErecoMinEnv->SetLineColor(kRed);
  hErecoMaxEnv->SetFillColor(kRed);
  hErecoMinEnv->SetFillColor(10); //white, to unfill the lower area
  hErecoMaxEnv->SetFillStyle(3244); // hatched
  hErecoMinEnv->SetFillStyle(1001); // solid (white, to unfill the lower area)
  hErecoMaxEnv->Draw();
  hErecoMinEnv->Draw("same");
  hErecoNominal->Draw("same");
  c->SaveAs(Form("%s_%s_envelope_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
  
  TH1D * hErecoMaxEnv_per50MeV = spectrum_to_per_50mev(hErecoMaxEnv);
  TH1D * hErecoMinEnv_per50MeV = spectrum_to_per_50mev(hErecoMinEnv);
  
  TCanvas * c2 = new TCanvas();
  hErecoMaxEnv_per50MeV->Draw();
  hErecoMinEnv_per50MeV->Draw("same");
  hErecoNominal_per50MeV->Draw("same");
  c2->SaveAs(Form("%s_%s_envelope_per_50MeV_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
  
  
  c2->cd();
  double max_limit = 1.1*TMath::Max(hErecoMax->GetMaximum(),hErecoMin->GetMaximum());
  double min_limit = 1.1*TMath::Min(hErecoMax->GetMinimum(),hErecoMin->GetMinimum());
  max_limit = TMath::Max(fabs(max_limit),fabs(min_limit));
  min_limit = -max_limit;
  hErecoMax->SetMaximum(max_limit);
  hErecoMax->SetMinimum(min_limit);
  hErecoMax->Draw();
  hErecoMin->Draw("same");
  c2->SaveAs(Form("%s_%s_ratio_%s.eps",sname,modetitle.c_str(),oscparams.c_str()));
  
  if (save_throws) {
    TCanvas * csys = new TCanvas();
    csys->Divide(ceil(sqrt(np)),ceil(sqrt(np)));
    for(int ip=0; ip < np; ip++) {
      //if(!fparm->IsNuisance(p)) continue; 
      csys->cd(ip+1);
      hSysThrown[ip]->Draw();
      cout << ip << "\t" << hSysThrown[ip]->GetTitle() << endl;
    }
    delete csys;
  }
 
  delete c2;

  f.Close();

  // format info to write
  stringstream streamer;
  if(!update_rootfile) {
    streamer << "\\begin{table}"   << endl;
    streamer << "\\begin{tabular}{| c | c c |}" << endl;
    streamer << "\\hline" << endl;
    streamer << "Nsk RMS/Mean & "
             << "BANFF pre & "
             << "BANFF post "
             << endl;
    streamer << "\\hline" << endl;
  }
  string banffpre = "-";
  if(mode == 0 || mode == 1)
    banffpre = "?";
    
  double dn_over_n = hNsk->GetRMS() / hNsk->GetMean();
  
  if (use_new_definition_of_dn_over_n) {
    double sum_nsk = 0;
    for (int iexp = 0 ; iexp < nexp; iexp++) {
      sum_nsk += nsk_per_toy_arr[iexp];
    }
    double mean_nsk = sum_nsk / nexp;
    
    double sumsq = 0;
    for (int iexp = 0 ; iexp < nexp; iexp++) {
      sumsq += pow(nsk_per_toy_arr[iexp]-mean_nsk,2);
    }
    
    double expected_deviation = sqrt(sumsq/nexp);
    cerr<<"Expected deviation: "<<expected_deviation<<endl;
    cerr<<"Nominal total SK: "<<nominal_total_nsk<<endl;
    
    dn_over_n = expected_deviation/nominal_total_nsk;
    
  }
      
  
  if (isprefit) {
    streamer << modetitle << " & & " << dn_over_n << " \\\\ \\hline" << endl;
  }
  else {
    streamer << modetitle << " & & " << dn_over_n << " \\\\ \\hline" << endl;
  }
    
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


#ifndef __CINT__
int main()
{
  valor::Env::Instance()->AutoSet("t2xk");
  TString conf_dir  = valor::Env::Instance()->TopLevelConfDir();
  TString verbosity = conf_dir + "/verbosity.xml";
  valor::util::Verbosity(verbosity.Data());
  test_sk_ereco_spectrum_twk_sys_throw(0,false,"RHC1Rmu","3f_osc",0,0,0,0,-99,-99,false,false,"analysis_configuration_2.5E20.xml",100000,"sinsqtheta","normal");
  return 0;
}
#endif


