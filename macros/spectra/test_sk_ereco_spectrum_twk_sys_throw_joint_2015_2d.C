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

void test_sk_ereco_spectrum_twk_sys_throw_joint_2015_2d
(
 int          mode                     = 0,                                         /* options: see above */
 bool         update_rootfile          = true,                                        /* "update" ROOT output file? If false, will "recreate" */
 const char * sk_sample                = "RHC1Re",                                   /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
 const char * hypothesis               = "3f_osc",                                    /* options: "2f_numu_disap", "3f_osc", "3f_oscbar" */
 
 double       input_mixing12       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp        = -1.601, //-99,    /* deltacp, */
 double       input_dm12sq         = 7.53E-5,    /* Dm^2_{12},  */
 double       input_dm23sq         = 2.509E-3, //-99,    /* |Dm^2_{23}| */
 double       input_betabar        = -99,
 
 bool         save_toy_expts           = false,                                       /* save the ereco spectra? */
 bool         save_throws              = false,                                        /* save the thrown parameter values*/

 //const char * config_file      = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml",
 //const char * config_file      = "analysis_configuration_run1to4_joint_Nu6.57E20POT.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_Nu6.57E20POT.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_FullStat.xml",
 const char * config_file      = "analysis_configuration_test_theta_1bin.xml",

 const int    nexp                     = 200,                                          /* number of experiments to make */
 const char * dmsqscheme               = "auto",
 const char * mixparamchoice           = "sinsqtheta",
 const char * hierarchy                = "normal",
 const char * mixparamoctant           = "two",       
 bool         reject_unphysical_throws = true    ,
 bool         use_new_definition_of_dn_over_n = true,
 
 const bool   doMargMH             = false,      // marginalize Mass Hierarchy
 const double PriorProbNH          = 0.5, //0.684, // Mach3  // Prior for MH --> give the prob. of NH as input
 TString file_priors = "" // Set to "" if no user priors are used (VALOR) 
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  //

  TH1D *hMEC = new TH1D("hMEC","hMEC",200,-5,5);

  // Take prior distributions
 
  TH2D *hToy_th13_dcp_nh;
  TH2D *hToy_th13_dcp_ih;
  TH2D *hToy_23_nh;
  TH2D *hToy_23_ih;
  TH1D *hToy_MH = new TH1D("hToy_MH","0-->NH, 1-->IH",2,0,2);

  TString name_prior_th13_dcp_nh = "hPrior_th13_dcp_nh";
  TString name_prior_th13_dcp_ih = "hPrior_th13_dcp_ih";
  TString name_prior_23_nh = "hPrior_L23_nh";
  TString name_prior_23_ih = "hPrior_L23_ih";
 
  TH2D *hPrior_th13_dcp_nh = 0;
  TH2D *hPrior_23_nh = 0;
  TH2D *hPrior_th13_dcp_ih = 0;
  TH2D *hPrior_23_ih = 0;
 
  bool UsePrior_th13_dcp = false;
  bool UsePrior_23 = false;
  
  if(file_priors!=""){

    UsePrior_th13_dcp = true;
    UsePrior_23 = true;
    
    cout << endl;
    cout << "Take priors from: " << file_priors << endl;
    TFile *fPrior = new TFile(file_priors.Data(),"READ");
    if(!fPrior->IsOpen() || !fPrior){
      cerr << endl;
      cerr << "Error: the file " << file_priors << " is not open!" << endl;
      cerr << endl;
      exit(1);
    }

    hPrior_th13_dcp_nh = (TH2D*) fPrior->Get(name_prior_th13_dcp_nh.Data());        
    if(hPrior_th13_dcp_nh->IsZombie() || !hPrior_th13_dcp_nh){
      cerr << endl;
      cerr << "Error: the histogram theta13 Vs deltaCP normal hierarchy is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
    hPrior_th13_dcp_ih = (TH2D*) fPrior->Get(name_prior_th13_dcp_ih.Data());       
    if(hPrior_th13_dcp_ih->IsZombie() || !hPrior_th13_dcp_ih){
      cerr << endl;
      cerr << "Error: the histogram theta13 Vs deltaCP inverted hierarchy is not taken!" << endl;
      cerr << endl;
      exit(1);
    }

    hPrior_23_nh = (TH2D*) fPrior->Get(name_prior_23_nh.Data()); 
    if(hPrior_23_nh->IsZombie() || !hPrior_23_nh){
      cerr << endl;
      cerr << "Error: the histogram par23_hist_normal is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
    hPrior_23_ih = (TH2D*) fPrior->Get(name_prior_23_ih.Data()); 
    if(hPrior_23_ih->IsZombie() || !hPrior_23_ih){
      cerr << endl;
      cerr << "Error: the histogram par23_hist_inverted is not taken!" << endl;
      cerr << endl;
      exit(1);
    }
    

    hToy_th13_dcp_nh = (TH2D*) hPrior_th13_dcp_nh->Clone("hToy_th13_dcp_nh");
    hToy_th13_dcp_nh->Reset();
    hToy_th13_dcp_ih = (TH2D*) hPrior_th13_dcp_ih->Clone("hToy_th13_dcp_ih");
    hToy_th13_dcp_ih->Reset();
    
    hToy_23_nh = (TH2D*) hPrior_23_nh->Clone("hToy_23_nh");
    hToy_23_nh->Reset();
    hToy_23_ih = (TH2D*) hPrior_23_ih->Clone("hToy_23_ih");
    hToy_23_ih->Reset();
    
    // TCanvas *c0 = new TCanvas("c0","c0");
    // hPrior_th13_dcp_nh->Draw("colz");
    // TCanvas *c1 = new TCanvas("c1","c1");
    // hPrior_th13_dcp_ih->Draw("colz");
    // TCanvas *c2 = new TCanvas("c2","c2");
    // hPrior_23_nh->Draw("colz");
    // TCanvas *c3 = new TCanvas("c3","c3");
    // hPrior_23_ih->Draw("colz");
    // return;

  } // if(file_priors)


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
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  Opt::Instance()->SetDefaultParameterValues();

  // Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  // Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  // if(Opt::Instance()->Undefined()) {
  //   cerr << "Global options not defined correctly. Exiting..." << endl;
  //   cerr << *Opt::Instance() << endl;
  //   exit(1);
  // }
  string sample_name = util::AsStringCompact(s);
  const char * sname = sample_name.c_str();

  //

  TH1D* h_nominal;
  TH2D* h_2d = new TH2D("h_2d","Envelope",util::NumErecoBins(s),util::ErecoBinEdges(s),2000,0,3);
  
  TH1D** h_bin = new TH1D*[util::NumErecoBins(s)];
  TH1D** h_bin_coarse = new TH1D*[util::NumErecoBins(s)];
  for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){    
    char name[200];
    sprintf(name,"h_bin_%d",ibin);
    h_bin[ibin] = new TH1D(name,name,20000,0,1);
    //h_bin[ibin] = new TH1D(name,name,2000,0,1);
    sprintf(name,"h_bin_coarse_%d",ibin);
    h_bin_coarse[ibin] = new TH1D(name,name,200,0,1);
  }

  TH1D* hErecoAverage = new TH1D("hErecoAverage","hErecoAverage",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1D* h_up = new TH1D("h_up","h_up",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1D* h_down = new TH1D("h_down","h_down",util::NumErecoBins(s),util::ErecoBinEdges(s));

  TH1D* hErecoMedian = new TH1D("hErecoMedian","hErecoMedian",util::NumErecoBins(s),util::ErecoBinEdges(s)); 
  TH1D* h_up_Median = new TH1D("h_up_Median","h_up_Median",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1D* h_down_Median = new TH1D("h_down_Median","h_down_Median",util::NumErecoBins(s),util::ErecoBinEdges(s));

  //

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
    

  // Check whether some osc. params need to be randomized
  
  bool randomize_mixing12   = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13   = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_mixing23   = (TMath::Abs(input_mixing23 + 999) < kASmallNum);
  bool randomize_deltacp    = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq     = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);
  bool randomize_dm23sq     = (TMath::Abs(input_dm23sq   + 999) < kASmallNum);

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

  if(input_betabar < -1 ) input_betabar = conv::kGlobFitBetaBar;

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
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FNormBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypBetaBar,   input_betabar);
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
  
  // string oscparams = util::AsString(ph.Hypothesis(),
  //                      input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  string oscparams = util::AsString(ph.Hypothesis(),
			  input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq,input_betabar);
  
 
  // Output histograms
  TH1D * hErecoNominal = 0;
  VALORArray * obsEreco[nexp];
  for (int i = 0 ; i < nexp; i++) {obsEreco[i]=0;}
  TH1D * hNsk = new TH1D(Form("%s_%s_nsk",sname,modetitle.c_str()),Form("%s_%s_nsk",sname,modetitle.c_str()),100,0,10);
  
  
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

  FitParamList * fparm = util::NamedParameterList("nuosc3f+nuisance"); 
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
  hErecoNominal = (TH1D*) skexp.CurrPredictedErecoVsQ2Spectrum(s)->ProjectionY()->Clone("hErecoNominal");
  //hErecoNominal = skexp.CurrPredictedErecoSpectrum(s);
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
  int n_q2_bins    = util::NumQ2Bins(s);
  
  double * sumsq_per_bin = (double*)calloc(n_ereco_bins,sizeof(double));
  VALORArray * nsk_per_toy = new VALORArray(nexp);
  double * nsk_per_toy_arr = nsk_per_toy->GetArray();
  
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
      if(mode == 11 || mode == 0) {
	
       	double input_mixing12_used = conv::kBadResult;  
	double input_mixing13_used = conv::kBadResult;    
	double input_mixing23_used = conv::kBadResult;   
	double input_deltacp_used  = conv::kBadResult; 
	double input_dm12sq_used   = conv::kBadResult;   
	double input_dm23sq_used   = conv::kBadResult;   
	double input_beta_bar_used = conv::kBadResult;
	
	// Extract Mass Hierarchy
	
	conv::NuMassHierarchy_t hier;
	conv::NuMassHierarchy_t hierBar;
	
	//int iMH= util::RandomizeMH(PriorProbNH); // 0 --> NH, 1 --> IH
	
	int iMH = conv::kBadResult;  // Index of Mass Hierarchy: 0 --> NH, 1 --> IH
	
	// Don't marginalize Mass Hierarchy. Only the assumed one
	if(doMargMH == false){
	  hier    = Opt::Instance()->gMassHierarchyAssumed;
	  hierBar = Opt::Instance()->gMassHierarchyAssumed;
	  //  - iMH=0 --> kNmhNormal
	  //  - iMH=1 --> kNmhInverted	
	  if     ( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal   ) iMH = 0;
	  else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted ) iMH = 1; 
	  
	  // Mass Hierarchy has to be defined
	  else {
	    cerr << endl;
	    cerr << "Mass Hierarchy is not defined" << endl;
	    cerr << endl;
	    exit(1);
	  }
	}
	
	// Extract randomly the Mass Hierarchy 
	else	  
	  iMH = util::RandomizeMH(PriorProbNH);   // iMH = 0 --> NH, iMH = 1 --> IH

	hToy_MH->Fill(iMH+conv::kASmallNum);

	if( iMH==0 ){ // Normal Hierarchy
	  
	  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
	  hier    = kNmhNormal;
	  hierBar = kNmhNormal;
	  skexp.GetSuperkMC()->SetMassHierarchy(hier);
	  skexp.GetSuperkMC()->SetMassHierarchyBar(hierBar);
	  
	  // Randomize specified osc params for current toy MC expt
	  if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
	  if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }     
	  
	  if ( (randomize_dm23sq || randomize_mixing23) && UsePrior_23){
	    double dm=0;
	    double mix23=0;
	    hPrior_23_nh->GetRandom2(mix23,dm);      
	    if(randomize_mixing23) input_mixing23 = mix23;
	    if(randomize_dm23sq)   input_dm23sq = dm;
	  }
	  
	  if ( (randomize_mixing13 || randomize_deltacp) && UsePrior_th13_dcp){
	    double mix13=0;
	    double dcp=0;
	    hPrior_th13_dcp_nh->GetRandom2(mix13,dcp);
	    if(randomize_mixing13) input_mixing13 = mix13;
	    if(randomize_deltacp)  input_deltacp = dcp;
	  }
	}
	
	else if(iMH==1){ // Inverted Hierarchy
	  
	  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
	  hier    = kNmhInverted;
	  hierBar = kNmhInverted;
	  skexp.GetSuperkMC()->SetMassHierarchy(hier);
	  skexp.GetSuperkMC()->SetMassHierarchyBar(hierBar);

	  // Randomize specified osc params for current toy MC expt 
	  if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
	  if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }     
	  
	  if ( (randomize_dm23sq || randomize_mixing23) && UsePrior_23){
	    double dm=0;
	    double mix23=0;
	    hPrior_23_ih->GetRandom2(mix23,dm);      
	    if(randomize_mixing23) input_mixing23 = mix23;
	    if(randomize_dm23sq)   input_dm23sq = dm;
	  }
	  
	  if ( (randomize_mixing13 || randomize_deltacp) && UsePrior_th13_dcp){
	    double mix13=0;
	    double dcp=0;
	    hPrior_th13_dcp_ih->GetRandom2(mix13,dcp);
	    if(randomize_mixing13) input_mixing13 = mix13;
	    if(randomize_deltacp)  input_deltacp = dcp;
	  }
	}
	else{
	  cerr << endl;
	  cerr << "Mass Hierarchy is not defined" << endl;
	  cerr << "iMH = " << iMH << endl;
	  cerr << endl;
	  exit(1);
	}	
	  	
       	if(ph.Hypothesis() == kHypOsc3F) {  
	  ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
	  ph.SetCurrValue(conv::kPhypCPV,       input_deltacp);
	  ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
	  ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
	}
	else if(ph.Hypothesis() == kHypOsc3FBar) {
	  ph.SetCurrValue(conv::kPhypMixing13,     input_mixing13);
	  ph.SetCurrValue(conv::kPhypCPV,          input_deltacp);
	  ph.SetCurrValue(conv::kPhypMixing23,     input_mixing23);
	  ph.SetCurrValue(conv::kPhypDm23sq,       input_dm23sq);
	  ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
	  ph.SetCurrValue(conv::kPhypCPVBar,       input_deltacp);
	  ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
	  ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
	}
	else if(ph.Hypothesis() == kHypOsc3FNormBar) {
	  ph.SetCurrValue(conv::kPhypMixing13,     input_mixing13);
	  ph.SetCurrValue(conv::kPhypCPV,          input_deltacp);
	  ph.SetCurrValue(conv::kPhypMixing23,     input_mixing23);
	  ph.SetCurrValue(conv::kPhypDm23sq,       input_dm23sq);
	  ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
	  ph.SetCurrValue(conv::kPhypCPVBar,       input_deltacp);
	  ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
	  ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
	  ph.SetCurrValue(conv::kPhypBetaBar,      input_betabar);
	}
       
	//cout << "- Physics hypothesis set: " << ph << endl;

	double mix12 = ph.CurrValue(kPhypMixing12);
	double mix13 = ph.CurrValue(kPhypMixing13);
	double mix23 = ph.CurrValue(kPhypMixing23);
	double dm12 = ph.CurrValue(kPhypDm12sq);
	double dm23 = ph.CurrValue(kPhypDm23sq);
	double dcp = ph.CurrValue(kPhypCPV);
	
	// if(iMH==0){
	//   hToy_th13_dcp_nh->Fill(mix13,dcp);
	//   hToy_23_nh->Fill(mix23,dm23);
	// }
	// else if(iMH==1){
	//   hToy_th13_dcp_ih->Fill(mix13,dcp);
	//   hToy_23_ih->Fill(mix23,dm23);
	// }
  
	//cout << "mix13 = " << mix13 << " - dcp = " << dcp << endl;

      } // mode 11 (osc param)

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
    
    
    if(mode==425) hMEC->Fill(v[25]);

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
        //spectrum = new VALORArray(*skobs.ErecoSpectrumVArray(s));
        spectrum = new VALORArray(*skobs.TwoDSpectrum(s));
      }
      else {
        spectrum = new VALORArray(*skexp.CurrPredictedErecoVsQ2VArray(s));
        //spectrum = new VALORArray(*skexp.CurrPredictedErecoVArray(s));
      }
      
      spectrum->SetTitle(Form("%s_%s_exp%d",sname,modetitle.c_str(),iexp));
      spectrum->Write   (Form("%s_%s_exp%d",sname,modetitle.c_str(),iexp));
    }
    else {
      if(mode == 999 ) {
        //spectrum = skobs.ErecoSpectrumVArray(s);
        spectrum = skobs->TwoDSpectrum(s);
      }
      else {
        //spectrum = skexp.CurrPredictedErecoVArray(s);
        spectrum = skexp.CurrPredictedErecoVsQ2VArray(s);
      }
    }
  
    //store Nsk
    double nskint = spectrum->Integral();
 
    cout << "nskint = " << nskint << endl;

    hNsk->Fill(nskint);
    nsk_per_toy_arr[iexp] = nskint;
    
    //double ntot = 0.;

    // Store tweaks.
    for (int ibin = 0 ; ibin < n_ereco_bins ; ibin++) {
      
      double N=0.;
      for (int ibinQ2 = 0 ; ibinQ2 < n_q2_bins ; ibinQ2++) {
	N += spectrum->GetEntry(ibin,ibinQ2); 
      }
      //double N = spectrum->GetArray()[ibin]; // TO CHANGE

      double nom   = hErecoNominal->GetBinContent(ibin+1);
      
      //cout << "nom = " << nom << " - N = " << N << endl;

      sumsq_per_bin[ibin] += (N-nom)*(N-nom);

      h_2d->Fill(h_2d->GetXaxis()->GetBinCenter(ibin+1),N);

      h_bin[ibin]->Fill(N);
      h_bin_coarse[ibin]->Fill(N);

      //ntot += N;
    }
     
    //cout << "ntot = " << ntot << endl;
    //cout << endl;
  
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
  
  cout << endl;
  cout << "Nominal = " << hErecoNominal->Integral() << endl;
  cout << "Average = " << hNsk->GetMean() << endl;
  cout << endl;
  
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
  
  // Assume symmetric errors wrt the average (not nominal)

  for(int i=0;i<util::NumErecoBins(s);i++){
    
    // wrt average
    double average = h_bin[i]->GetMean();
    double rms = h_bin[i]->GetRMS();
    hErecoAverage->SetBinContent(i+1,average);
    h_up->SetBinContent(i+1,average+rms);
    h_down->SetBinContent(i+1,average-rms);  

    // wrt median
    int medianbin = GetMedianStatBin(h_bin[i],0.5); // median
    double medianstat = h_bin[i]->GetBinLowEdge(medianbin);
    int median_p1sigma_bin = GetMedianStatBin(h_bin[i],0.5+0.34); // median +1sigma
    double medianstat_p1sigma = h_bin[i]->GetBinLowEdge(median_p1sigma_bin); 
    int medianstat_m1sigma_int = GetMedianStatBin(h_bin[i],0.5-0.34); // median -1sigma
    double medianstat_m1sigma = h_bin[i]->GetBinCenter(medianstat_m1sigma_int); 
    hErecoMedian->SetBinContent(i+1,medianstat);
    h_up_Median->SetBinContent(i+1,medianstat_p1sigma);
    h_down_Median->SetBinContent(i+1,medianstat_m1sigma);
  }
  
  hErecoAverage->Write(Form("%s_%s_average",sname,modetitle.c_str()));
  h_up   ->Write(Form("%s_%s_up_average" ,sname,modetitle.c_str()));
  h_down ->Write(Form("%s_%s_down_average",sname,modetitle.c_str()));

  hErecoMedian->Write(Form("%s_%s_median",sname,modetitle.c_str()));
  h_up_Median   ->Write(Form("%s_%s_up_median" ,sname,modetitle.c_str()));
  h_down_Median ->Write(Form("%s_%s_down_median",sname,modetitle.c_str()));
  
  TCanvas * cAverage = new TCanvas("cAverage","Wrt average");
  h_up->GetXaxis()->SetTitle("E_{reco} (GeV)");
  h_up->GetYaxis()->SetTitle("Events per bin");
  h_down->GetXaxis()->SetTitle("E_{reco} (GeV)");
  h_down->GetYaxis()->SetTitle("Events per bin");
  h_up->SetLineStyle(3);
  h_down->SetLineStyle(3);
  h_up->SetLineColor(kRed);
  h_down->SetLineColor(kRed);
  h_up->SetFillColor(kRed);
  h_down->SetFillColor(10); //white, to unfill the lower area
  h_up->SetFillStyle(3244); // hatched
  h_down->SetFillStyle(1001); // solid (white, to unfill the lower area)
  h_up->DrawCopy("");
  h_down->DrawCopy("same");
  hErecoNominal->SetLineColor(kBlack);
  hErecoNominal->DrawCopy("same");
  hErecoAverage->SetLineColor(kRed);
  hErecoAverage->SetLineWidth(2);
  hErecoAverage->DrawCopy("same");
  hErecoMedian->SetLineColor(kBlue);
  hErecoMedian->SetLineWidth(2);
  hErecoMedian->DrawCopy("same");
  TLegend *lAverage = new TLegend(0.7,0.6,1.0,0.9);
  lAverage->SetFillColor(0);
  lAverage->SetLineColor(0);
  lAverage->SetTextFont(132);
  lAverage->SetTextSize(0.06);
  lAverage->AddEntry(hErecoNominal,"Nominal","l");
  lAverage->AddEntry(hErecoAverage,"Average","l");
  lAverage->AddEntry(hErecoMedian,"Median","l");
  lAverage->AddEntry(h_up,"Error wrt average","f");
  lAverage->Draw();
  
  TCanvas * cMedian = new TCanvas("cMedian","Wrt median");
  h_up_Median->GetXaxis()->SetTitle("E_{reco} (GeV)");
  h_up_Median->GetYaxis()->SetTitle("Events per bin");
  h_down_Median->GetXaxis()->SetTitle("E_{reco} (GeV)");
  h_down_Median->GetYaxis()->SetTitle("Events per bin");
  h_up_Median->SetLineStyle(3);
  h_down_Median->SetLineStyle(3);
  h_up_Median->SetLineColor(kRed);
  h_down_Median->SetLineColor(kRed);
  h_up_Median->SetFillColor(kRed);
  h_down_Median->SetFillColor(10); //white, to unfill the lower area
  h_up_Median->SetFillStyle(3244); // hatched
  h_down_Median->SetFillStyle(1001); // solid (white, to unfill the lower area)
  h_up_Median->DrawCopy("");
  h_down_Median->DrawCopy("same");
  hErecoNominal->SetLineColor(kBlack);
  hErecoNominal->DrawCopy("same");
  hErecoAverage->SetLineColor(kRed);
  hErecoAverage->SetLineWidth(2);
  hErecoAverage->DrawCopy("same");
  hErecoMedian->SetLineColor(kBlue);
  hErecoMedian->SetLineWidth(2);
  hErecoMedian->DrawCopy("same");
  TLegend *lMedian = new TLegend(0.7,0.6,1.0,0.9);
  lMedian->SetFillColor(0);
  lMedian->SetLineColor(0);
  lMedian->SetTextFont(132);
  lMedian->SetTextSize(0.06);
  lMedian->AddEntry(hErecoNominal,"Nominal","l");
  lMedian->AddEntry(hErecoAverage,"Average","l");
  lMedian->AddEntry(hErecoMedian,"Median","l");
  lMedian->AddEntry(h_up_Median,"Error wrt median","f");
  lMedian->Draw();
  
  TCanvas * c3 = new TCanvas();
  hNsk->Draw();
  hNsk->GetXaxis()->SetTitle("# of SK events");
  TLine *lNom = new TLine(hErecoNominal->Integral(),0,hErecoNominal->Integral(),hNsk->GetMaximum());
  lNom->SetLineColor(kRed);
  lNom->Draw();  

  TCanvas **c4 = new TCanvas*[util::NumErecoBins(s)];
  for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){

    //if(ibin+1!=10) continue;
    
    char name[200];
    sprintf(name,"ErecoBin_%d",ibin);
    c4[ibin] = new TCanvas(name,name);
    h_bin_coarse[ibin]->Draw();
    TLine *lNom_bin = new TLine(hErecoNominal->GetBinContent(ibin+1),0,hErecoNominal->GetBinContent(ibin+1),h_bin_coarse[ibin]->GetMaximum());
    lNom_bin->SetLineColor(kRed);
    lNom_bin->Draw();

    TLine *lNom_bin_up = new TLine(h_up->GetBinContent(ibin+1),0,h_up->GetBinContent(ibin+1),h_bin_coarse[ibin]->GetMaximum());
    lNom_bin_up->SetLineColor(kBlue);
    lNom_bin_up->Draw();

    TLine *lNom_bin_down = new TLine(h_down->GetBinContent(ibin+1),0,h_down->GetBinContent(ibin+1),h_bin_coarse[ibin]->GetMaximum());
    lNom_bin_down->SetLineColor(kBlue);
    lNom_bin_down->Draw();

    TLine *lMedian_bin = new TLine(hErecoMedian->GetBinContent(ibin+1),0,hErecoMedian->GetBinContent(ibin+1),h_bin_coarse[ibin]->GetMaximum());
    lMedian_bin->SetLineColor(kViolet);
    lMedian_bin->Draw();

    TLine *lNom_bin_up_Median = new TLine(h_up_Median->GetBinContent(ibin+1),0,h_up_Median->GetBinContent(ibin+1),h_bin_coarse[ibin]->GetMaximum());
    lNom_bin_up_Median->SetLineColor(kGreen);
    lNom_bin_up_Median->Draw();

    TLine *lNom_bin_down_Median = new TLine(h_down_Median->GetBinContent(ibin+1),0,h_down_Median->GetBinContent(ibin+1),h_bin_coarse[ibin]->GetMaximum());
    lNom_bin_down_Median->SetLineColor(kGreen);
    lNom_bin_down_Median->Draw();
    
    //   cout << "bin " << ibin << " : " 
    // 	 << endl 
    // 	 << " - nominal = " << hErecoNominal->GetBinContent(ibin+1)
    // 	 << endl
    // 	 // << " - mean = " 
    // 	 // << hErecoAverage->GetBinContent(ibin+1)
    // 	 // << " +" << h_up->GetBinContent(ibin+1) - hErecoAverage->GetBinContent(ibin+1)
    // 	 // << " -" << hErecoAverage->GetBinContent(ibin+1) - h_down->GetBinContent(ibin+1)
    // 	 // << endl
    //     	 // << " - median = " 
    // 	 // << hErecoMedian->GetBinContent(ibin+1) 
    // 	 // << " +" << h_up_Median->GetBinContent(ibin+1) - hErecoMedian->GetBinContent(ibin+1) 
    // 	 // << " -" << hErecoMedian->GetBinContent(ibin+1) - h_down_Median->GetBinContent(ibin+1)
    // 	 // << endl;
    // 	 << " - mean = " 
    // 	 << hErecoAverage->GetBinContent(ibin+1) << " - "
    // 	 << h_up->GetBinContent(ibin+1) << " - "
    // 	 << h_down->GetBinContent(ibin+1) 
    // 	 << endl
    //     	 << " - median = " 
    // 	 << hErecoMedian->GetBinContent(ibin+1) << " - "
    // 	 << h_up_Median->GetBinContent(ibin+1) << " - "
    // 	 << h_down_Median->GetBinContent(ibin+1)
    // 	 << endl;
  }
  cout << endl;
 
  TCanvas *c5 = new TCanvas("c5","c5");
  hMEC->Draw();

  //

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
  
  cout << "RMS / Mean = " << dn_over_n << endl;

  if (use_new_definition_of_dn_over_n) {
    double sum_nsk = 0;
    for (int iexp = 0 ; iexp < nexp; iexp++) {
      sum_nsk += nsk_per_toy_arr[iexp];
    }
    double mean_nsk = sum_nsk / nexp;
    
    double sumsq = 0;
    for (int iexp = 0 ; iexp < nexp; iexp++) {
      //cerr<<nsk_per_toy_arr[iexp]<<"\t"<<mean_nsk<<endl;
      sumsq += pow(nsk_per_toy_arr[iexp]-mean_nsk,2);
    }
    
    double expected_deviation = sqrt(sumsq/nexp);
    cerr<<"Expected deviation: "<<expected_deviation<<endl;
    cerr<<"Nominal total SK: "<<nominal_total_nsk<<endl;
    
    dn_over_n = expected_deviation/nominal_total_nsk;
    
    cout << "Use the following:" << endl;
    cout << "RMS / Nominal = " << dn_over_n << endl;
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





  // Draw histograms
  
  // TCanvas *cToy_th13_dcp_nh = new TCanvas("cToy_th13_dcp_nh","cToy_th13_dcp_nh");
  // hToy_th13_dcp_nh->Draw("colz");
  
  // TCanvas *cToy_th13_dcp_ih = new TCanvas("cToy_th13_dcp_ih","cToy_th13_dcp_ih");
  // hToy_th13_dcp_ih->Draw("colz");
  
  // TCanvas *cToy_23_nh = new TCanvas("cToy_23_nh","cToy_23_nh");
  // hToy_23_nh->Draw("colz");
  
  // TCanvas *cToy_23_ih = new TCanvas("cToy_23_ih","cToy_23_ih");
  // hToy_23_ih->Draw("colz");
  
  TCanvas *cToy_MH = new TCanvas("cToy_MH","cToy_MH");
  hToy_MH->Draw();
  
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


// double average(TH1D* h1){
//   double avg(0);
//   for(int i=1; i<=h1->GetNbinsX();i++){
//     avg += h1->GetBinCenter(i)*h1->GetBinContent(i);
//   }
//   avg = avg/h1->GetEntries();  
//   return avg;
// }

// double rms(TH1D* h1,double avg){
//   double rms(0);  
//   for(int i=1; i<=h1->GetNbinsX();i++){
//     rms += ((h1->GetBinCenter(i) - avg) **2) * (h1->GetBinContent(i));
//   }
//   rms = rms/h1->GetEntries();
//   rms = rms**0.5;
//   return rms;
// }

TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}

double GetMedianStatBin(TH1D *h,double frac){
  
  int binmedian = -999;
  double totevents=h->Integral();

  double currevt=0;
  
  for(int ibin=0;ibin<h->GetNbinsX();ibin++){
    
    currevt += h->GetBinContent(ibin+1);
    
    if( (currevt/totevents) >= frac ){
      binmedian = ibin+1;
      break;
    }
  }

  //double median = h->GetBinCenter(binmedian);
  
  return binmedian;
}
