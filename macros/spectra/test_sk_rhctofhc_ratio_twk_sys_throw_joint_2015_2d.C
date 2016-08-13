//
// ------------------------------------------------------------------------------------------------
// Calculate SuperK rhc/fhc ratio of # of events distribution due to systematic parameter variations,
// randomised using Gaussians or Cholesky decomposition (for FSI/BANFF/CCQERingCounting)
//
// There are many run modes:
//
//   0 tweak all parameters
//   1 tweak all BANFF
//   2 tweak all unconstrained xsec
//   5 tweak only BANFF unconstrained parameters
//   6 tweak all BANFF constrained parameters
//   7 tweak BANFF flux only
//   8 tweak BANFF xsec only
//   9 tweak SKdet + FSI
//  10 tweak BANFF + unconstrained xsec
//  11 oscillation parameters
//  12 sigma nue/numu (nu and nubar) parameters
//  13 tweak all parameters except oscillations
//  14 tweak FSI (use different config file)
//  15 tweak SKDet (use different config file)
//  16 tweak Escale
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
// banff fhc only
// 400 tweak banff00
// 401 tweak banff01
// 425 tweak banff25 (MEC)
// 441 tweak banff41 (MECbar)
// 400+kNumOfParamsInBANFFfit tweak last banff parameter
//
// banff rhc only
// 500 tweak banff00_rhc
// 501 tweak banff01_rhc
// 425 tweak banff25 (MEC)
// 441 tweak banff41 (MECbar)
// 500+kNumOfRHCParamsInBANFFfit tweak last banff parameter
//
// banff fhc+rhc w/ correlations
// 600 tweak banff00_fhcprhc
// 601 tweak banff01_fhcprhc
// 650 tweak banff50_fhcprhc (MEC)
// 662 tweak banff62_fhcprhc (NCother)
// 664 tweak banff64_fhcprhc (NC1gamma)
// 666 tweak banff66_fhcprhc (MECbar)
// 600+kNumOfFHCpRHCParamsInBANFFfit tweak last banff parameter
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

#include "Utilities/RunOpt.h"
#include "Utilities/OscParams.h"

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
#include "UtlForMacrosT2XK.h"

#endif


TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}

void test_sk_rhctofhc_ratio_twk_sys_throw_joint_2015_2d
(
 int          mode                     = 0, //650, //11,                                         /* options: see above */
 bool         update_rootfile          = true,                                        /* "update" ROOT output file? If false, will "recreate" */
  
 const char * hypothesis               = "3f_osc",                                    /* options: "2f_numu_disap", "3f_osc", "3f_oscbar" */

 // double       input_mixing12       = 0.304,    
 // double       input_mixing13       = 0.0217, //0.0243   
 // double       input_mixing23       = 0.528,    
 // double       input_deltacp        = -1.601,   
 // double       input_dm12sq         = 7.53E-5,   
 // double       input_dm23sq         = 2.509E-3,     
 // double       input_betabar        = 1.,

 // Asimov A
 double       input_mixing12       = -999,    
 double       input_mixing13       = -999,   
 double       input_mixing23       = 0.528,    
 double       input_deltacp        = -1.601,  
 double       input_dm12sq         = -999,   
 double       input_dm23sq         = 2.509E-3,     
 double       input_betabar        = 1.,

 // // Asimov B
 // double       input_mixing12       = -999,    
 // double       input_mixing13       = -999,   
 // double       input_mixing23       = 0.45,    
 // double       input_deltacp        = 0.,  
 // double       input_dm12sq         = -999,   
 // double       input_dm23sq         = 2.509E-3,     
 // double       input_betabar        = 1.,


 //
 // NEW CONFIGURATION FILES (4/6/16 BANFF RELEASE W/ FIT BIAS AND ADDITIONAL UNCERTAINTY ON RUN 1-7 FLUX)
 //
 const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17_FSIOnly.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17_SKDetOnly.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_Run17.xml", // Run 1-6 flux tuning

 //
 // NEW CONFIGURATION FILES (4/3/16 BANFF RELEASE)
 //
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias.xml",
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_ExpSummer2016.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_FSIOnly.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_SKDetOnly.xml",
 
 const int    nexp                     = 10000,                                          /* number of experiments to make */
 const char * dmsqscheme               = "auto",
 const char * mixparamchoice           = "sinsqtheta",
 const char * hierarchy                = "normal",
 const char * mixparamoctant           = "two",       
 bool         reject_unphysical_throws = true    ,
 bool         use_new_definition_of_dn_over_n = false,
 
 const bool   doMargMH             = false,      // marginalize Mass Hierarchy
 const double PriorProbNH          = 0.5, //0.684, // Mach3  // Prior for MH --> give the prob. of NH as input
 TString file_priors = "" // Set to "" if no user priors are used (VALOR) 
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
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
  case 12:
    modetitle = "NueNumu"; break;
  case 13:
    modetitle = "skdetfsi_banff"; break;
  case 14:
    modetitle = "FSIonly"; break;
  case 15:
    modetitle = "SKDetonly"; break;
  case 16:
    modetitle = "Escale_SK"; break;
    
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
    // banff fhc+rhc
    if(mode >= 600 && mode < 600+conv::kNumOfFHCpRHCParamsInBANFFfit) {
      modetitle = Form("banff%02d_fhcprhc",mode-600); break;
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
  
  //
  
  TString config_str(config_file);
  bool isprefit = config_str.Contains("prefit");
  if(isprefit)
    modetitle += "prefit";
  
  //
  

  string sample_name_RHC1Re = util::AsStringCompact(kSKsRHCFCFV1RingELike);
  const char * sname_RHC1Re = sample_name_RHC1Re.c_str();

  string sample_name_RHC1Rmu = util::AsStringCompact(kSKsRHCFCFV1RingMuLike);
  const char * sname_RHC1Rmu = sample_name_RHC1Rmu.c_str();

  string sample_name_1Re = util::AsStringCompact(kSKsFCFV1RingELike);
  const char * sname_1Re = sample_name_1Re.c_str();

  string sample_name_1Rmu = util::AsStringCompact(kSKsFCFV1RingMuLike);
  const char * sname_1Rmu = sample_name_1Rmu.c_str();

  TH1D * hNsk_RHC1Re = new TH1D(Form("%s_nsk",sname_RHC1Re),Form("%s_nsk",sname_RHC1Re),100,0,10);
  TH1D * hNsk_1Re = new TH1D(Form("%s_nsk",sname_1Re),Form("%s_nsk",sname_1Re),100,0,100);
  TH1D * hNsk_RHC1Rmu = new TH1D(Form("%s_nsk",sname_RHC1Rmu),Form("%s_nsk",sname_RHC1Rmu),100,0,100);
  TH1D * hNsk_1Rmu = new TH1D(Form("%s_nsk",sname_1Rmu),Form("%s_nsk",sname_1Rmu),100,0,300);      

  // TH1D * hNsk_RHC1Re = new TH1D(Form("%s_%s_nsk",sname_RHC1Re,modetitle.c_str()),Form("%s_%s_nsk",sname_RHC1Re,modetitle.c_str()),100,0,10);
  // TH1D * hNsk_1Re = new TH1D(Form("%s_%s_nsk",sname_1Re,modetitle.c_str()),Form("%s_%s_nsk",sname_1Re,modetitle.c_str()),100,0,100);
  // TH1D * hNsk_RHC1Rmu = new TH1D(Form("%s_%s_nsk",sname_RHC1Rmu,modetitle.c_str()),Form("%s_%s_nsk",sname_RHC1Rmu,modetitle.c_str()),100,0,100);
  // TH1D * hNsk_1Rmu = new TH1D(Form("%s_%s_nsk",sname_1Rmu,modetitle.c_str()),Form("%s_%s_nsk",sname_1Rmu,modetitle.c_str()),100,0,300);      
 

  TH1D * hRatio_FHCtoRHC_elike  = new TH1D(Form("%s_ratio","Re"),Form("%s_nsk","Re"),200,0,20);
  TH1D * hRatio_FHCtoRHC_mulike = new TH1D(Form("%s_ratio","Rmu"),Form("%s_nsk","Rmu"),200,0,20);
      
  cout << endl;
  cout << "Systematic mode (" << mode << ") --> " << modetitle << endl;
  cout << endl;
  
  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  string sconfig(config);  
  T2KOpt::Instance()->ReadXmlPath(sconfig);    
  
  // Check whether some osc. params need to be randomized
  
  bool randomize_mixing12   = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13   = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_mixing23   = (TMath::Abs(input_mixing23 + 999) < kASmallNum);
  bool randomize_deltacp    = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq     = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);
  bool randomize_dm23sq     = (TMath::Abs(input_dm23sq   + 999) < kASmallNum);
  
  // Set Physics Hypothesis
  PhysHypothesis ph(hypothesis);
  
  // Set oscillation parameter values and priors
  
  // util::OscParams::Instance()->SetValues(input_mixing12, input_mixing13, input_mixing23, input_deltacp, input_dm12sq, input_dm23sq,
  // 					 input_betabar,
  // 					 input_mixing12, input_mixing13, input_mixing23, input_deltacp, input_dm12sq, input_dm23sq,
  // 					 PriorProbNH, doMargMH
  // 					 );
  
  // Use same values for Nu and NuBar oscillation parameters 
  util::OscParams::Instance()->SetValues(input_mixing12, input_mixing13, input_mixing23, input_deltacp, input_dm12sq, input_dm23sq,
					 input_betabar,
					 input_mixing12, input_mixing13, input_mixing23, input_deltacp, input_dm12sq, input_dm23sq,
					 PriorProbNH, doMargMH,
					 input_mixing12, input_mixing13, input_mixing23, input_deltacp, input_dm12sq, input_dm23sq,
					 input_mixing12, input_mixing13, input_mixing23, input_deltacp, input_dm12sq, input_dm23sq
					 );

  bool MargMixing12        = true;
  string PriorMixing12       = "gaus";   
  double PriorParMixing12[2];
  PriorParMixing12[0] = conv::kPrior_Sinsq2th12_Used_Mean; 
  PriorParMixing12[1] = conv::kPrior_Sinsq2th12_Used_RMS;  
  string MixModeMixing12     = "sinsq2th";
  
  bool MargMixing13        = true;
  string PriorMixing13       = "gaus";   
  double PriorParMixing13[2];
  PriorParMixing13[0] = conv::kPrior_Sinsq2th13_Used_Mean; 
  PriorParMixing13[1] = conv::kPrior_Sinsq2th13_Used_RMS; 
  string MixModeMixing13     = "sinsq2th";
  
  bool MargMixing23        = true;
  string PriorMixing23       = "flat";   
  double PriorParMixing23[2];
  PriorParMixing23[0] = conv::kPrior_Sinsq2th23_Flat_Low; 
  PriorParMixing23[1] = conv::kPrior_Sinsq2th23_Flat_Up;  
  string MixModeMixing23     = "sinsqth";
  
  bool MargDm23sq          = true;
  string PriorDm23sq         = "flat";   
  double PriorParDm23sq[2];
  PriorParDm23sq[0]   = conv::kPrior_Dm23sq_Flat_Low; 
  PriorParDm23sq[1]   = conv::kPrior_Dm23sq_Flat_Up;  
  
  bool MargDm12sq          = true;
  string PriorDm12sq         = "gaus";   
  double PriorParDm12sq[2];
  PriorParDm12sq[0]   = conv::kPrior_Dm12sq_Used_Mean; 
  PriorParDm12sq[1]   = conv::kPrior_Dm12sq_Used_RMS;  
  
  bool MargDeltaCP         = true;
  string PriorDeltaCP        = "flat";   
  double PriorParDeltaCP[2];
  PriorParDeltaCP[0]  = conv::kPrior_deltaCP_Flat_Low; 
  PriorParDeltaCP[1]  = conv::kPrior_deltaCP_Flat_Up;  
  
  if(randomize_mixing23) 
    util::OscParams::Instance()->SetPrior(conv::kPhypMixing23, PriorMixing23, PriorParMixing23[0], PriorParMixing23[1], MixModeMixing23);      
  if(randomize_dm23sq)   
    util::OscParams::Instance()->SetPrior(conv::kPhypDm23sq,   PriorDm23sq,   PriorParDm23sq[0],   PriorParDm23sq[1],   ""              );  
  if(randomize_mixing13) 
    util::OscParams::Instance()->SetPrior(conv::kPhypMixing13, PriorMixing13, PriorParMixing13[0], PriorParMixing13[1], MixModeMixing13); // Convert to sinsqth
  if(randomize_deltacp)  
    util::OscParams::Instance()->SetPrior(conv::kPhypCPV,      PriorDeltaCP,  PriorParDeltaCP[0],  PriorParDeltaCP[1],  ""              ); 
  if(randomize_mixing12) 
    util::OscParams::Instance()->SetPrior(conv::kPhypMixing12, PriorMixing12, PriorParMixing12[0], PriorParMixing12[1], MixModeMixing12); // Convert to sinsqth
  if(randomize_dm12sq)   
    util::OscParams::Instance()->SetPrior(conv::kPhypDm12sq,   PriorDm12sq,   PriorParDm12sq[0],   PriorParDm12sq[1],   ""              ); 
  
  // Update physics hypothesis and fitparamlist with chosen value of oscillation parameters
  util::OscParams::Instance()->UpdatePhysHypothesis(ph);
  
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
  SuperkExpectation *skexp = new SuperkExpectation;
  skexp->Read(sconfig);
  
  SuperkObservation skobs;
  
  string oscparams = util::AsString(ph.Hypothesis(),input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq,input_betabar,-999,-999);

  FitParamList * fparm = util::NamedParameterList("nuosc3f+nuisance"); 
  const int np = fparm->GetNumberParameters();
  
  // Open output files
  string rootmode = "update";
  if(!update_rootfile)
    rootmode = "recreate";

  ofstream * outfilep = 0;
  if(!update_rootfile){
    outfilep = new ofstream(Form("sk_rhctofhc_ratio-twk_sys_throw_%s_%s.txt",oscparams.c_str(),modetitle.c_str()),ios_base::trunc);
  }
  else {
    outfilep = new ofstream(Form("sk_rhctofhc_ratio-twk_sys_throw_%s_%s.txt",oscparams.c_str(),modetitle.c_str()),ios_base::app);
  }  
  ofstream & outfile = *outfilep;
  
  //remove the [] so can use it in SaveAs
  while(oscparams.find_first_of("[]") != string::npos) {
    oscparams.erase(oscparams.find_first_of("[]"), 1);
  }
  
  TStopwatch w;
  w.Start();
  TStopwatch w2;
  w2.Start();
  
  bool done_print = true;
    
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
	skexp->TweakAllSystematics();
      }
      
      // BANFF
      if(mode == 1) {
	skexp->TweakAllSampleBANFFSystematics();
      }
      
      // Uncorrelated XSec
      if(mode == 2) {
	skexp->TweakAllUncorrelatedXSecSystematics();
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
      }
      
      //BANFF (unconstrained only)
      if(mode == 5) {
	skexp->GetBANFF()->TweakSampleConstrainedOnly(false,0,true); // cholesky w/ full unconstrained matrix
      }
      //BANFF (constrained only)
      if(mode == 6) {
	skexp->GetBANFF()->TweakSampleConstrainedOnly(false,0, false); // cholesky w/ full cov matrix
      }
      
      //BANFF (flux only)
      if(mode == 7) {
	skexp->GetBANFF()->TweakSampleFluxOnly(false,0); // cholesky w/ full cov matrix
      }
      
      //BANFF (xsec only)
      if(mode == 8) {
	skexp->GetBANFF()->TweakSampleXSecOnly(false,0); // cholesky w/ full cov matrix
      }
      
      //SK det + FSI, FSI only or SKDet only
      if(mode == 9 || mode == 14 || mode == 15) {
	skexp->TweakAllSampleSKDetectorFSISystematics();
      }
      
      //unconstrained xsec + BANFF
      if(mode == 10) {
	skexp->TweakAllUncorrelatedXSecSystematics();
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	skexp->TweakAllBANFFSystematics();
      }
      
      //sigma nue/numu (nu and nubar) xsec (they are correlated)
      if(mode == 12) {
	skexp->GetBANFF()->TweakSampleXSecNumuNueOnly(false,0); // cholesky w/ full cov matrix
      }

      // All except oscillations
      if(mode == 13) {
	skexp->TweakAllSystematics();
      }

      // Energy Scale
      if(mode == 16){
	// Same as done in SuperkExpectation::TweakAllFHCpRHCSKDetectorFSISystematics()
	skexp->GetSuperkDetectorFSI()->TweakEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
	skexp->GetSuperkDetectorFSI()->TweakRHCEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));

	// cout << "FHC Escale = " << skexp->GetSuperkDetectorFSI()->FHCpRHCCurrentValue(kNuisSuperkEnergyScale) << endl;
	// cout << "FHC Escale = " << skexp->GetSuperkDetectorFSI()->CurrentValue(kNuisSuperkEnergyScale) << endl;
	// cout << "RHC Escale = " << skexp->GetSuperkDetectorFSI()->FHCpRHCCurrentValue(kRHCNuisSuperkEnergyScale) << endl;
	// cout << "RHC Escale = " << skexp->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSuperkEnergyScale) << endl;
      }

      //non-23 oscillation parameters
      if(mode == 11 || mode == 0) {
	
	// Randomize specified osc params for current toy MC expt
	util::OscParams::Instance()->Randomize();
	
	//util::OscParams::Instance()->Print(cout);
	
	// Set the mass hierarchy
	//valor::util::SetMassHierarchy(skexp);
	util::SetMassHierarchy(skexp);
	
	// Update physics hypothesis and fitparamlist with chosen value of oscillation parameters
	bool print = false;
	util::OscParams::Instance()->UpdatePhysHypothesis(ph,print);
	bool set_init = false;
	util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);
	
	// Fill toy histograms of osc params
	util::OscParams::Instance()->FillOscToyHisto(ph);
	
      } // mode 11 (osc param)
      
    }//mode<100
    
    else if(mode < 999) {//everything else but the stat
      
      double tweak = rand.Gaus(0.,1.);
      
      switch(mode) {
      case 200:
	if( !util::IsXSec(util::PhysParamIdFromString("pfermi")) ) break; 
	skexp->TweakSystematic("pfermi",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 201:
	if( !util::IsXSec(util::PhysParamIdFromString("specfunc")) ) break; 
	skexp->TweakSystematic("specfunc",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 202:
	if( !util::IsXSec(util::PhysParamIdFromString("wshape")) ) break; 
	skexp->TweakSystematic("wshape",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 203:
	if( !util::IsXSec(util::PhysParamIdFromString("shapeCCoth")) ) break; 
	skexp->TweakSystematic("shapeCCoth",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 204:
	if( !util::IsXSec(util::PhysParamIdFromString("pilessdelta")) ) break; 
	skexp->TweakSystematic("pilessdelta",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 205:
	if( !util::IsXSec(util::PhysParamIdFromString("norm_cccoh")) ) break; 
	skexp->TweakSystematic("norm_cccoh",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 206:
	if( !util::IsXSec(util::PhysParamIdFromString("norm_nc1piC")) ) break; 
	skexp->TweakSystematic("norm_nc1piC",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 207:
	if( !util::IsXSec(util::PhysParamIdFromString("norm_ncoth")) ) break; 
	skexp->TweakSystematic("norm_ncoth",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 208:
	if( !util::IsXSec(util::PhysParamIdFromString("norm_nuetonumu")) ) break; 
	skexp->TweakSystematic("norm_nuetonumu",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 209:
	if( !util::IsXSec(util::PhysParamIdFromString("norm_antinutonu")) ) break; 
	skexp->TweakSystematic("norm_antinutonu",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
      case 210:
	if( !util::IsXSec(util::PhysParamIdFromString("bindE")) ) break; 
	skexp->TweakSystematic("bindE",tweak); 
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
	
      case 211:
	if( !util::IsXSec(util::PhysParamIdFromString("rpa")) ) break; 
	tweak = rand.Uniform(0.,1.);
	skexp->TweakSystematic("rpa",tweak); // Automatically set at 0 or 1 in the Tweak
	skexp->GetNuXSec()->SyncRHCToFHCSystematics(); // need to synchronize FHC and RHC syst
	break;
	
      case 700:
	ph.SetCurrValue(conv::kPhypMixing13,  util::RandomizeMixing13()); break;
      case 701:
	ph.SetCurrValue(conv::kPhypMixing12,  util::RandomizeMixing12()); break;
      case 702:
	ph.SetCurrValue(conv::kPhypDm12sq,    util::RandomizeDm12sq()); break;
      case 703:
	ph.SetCurrValue(conv::kPhypCPV    ,   util::RandomizeCPV()); break;
	
      default:
	
	// banff fhc only
	if(mode >= 400 && mode < 400+conv::kNumOfParamsInBANFFfit) {
	  TVectorD v(conv::kNumOfParamsInBANFFfit);
	  v = skexp->GetBANFF()->BestFitValues();
	  int ib = mode - 400;
	  v[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	  skexp->GetBANFF()->TweakSample(true,&v);
	  break;
	}
	// banff rhc
	else if(mode >= 500 && mode < 500+conv::kNumOfRHCParamsInBANFFfit) {
	  TVectorD v(conv::kNumOfRHCParamsInBANFFfit);
	  v = skexp->GetBANFF()->RHCBestFitValues();
	  int ib = mode - 500;
	  v[ib] = skexp->GetBANFF()->RHCBestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->RHCCovMatrix()[ib][ib]));
	  skexp->GetBANFF()->TweakSample(true,&v);
	  break;
	}
	// banff fhc+rhc
	else if(mode >= 600 && mode < 600+conv::kNumOfFHCpRHCParamsInBANFFfit) {
	  TVectorD v(conv::kNumOfFHCpRHCParamsInBANFFfit);
	  v = skexp->GetBANFF()->FHCpRHCBestFitValues();
	  int ib = mode - 600;
	  v[ib] = skexp->GetBANFF()->FHCpRHCBestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->FHCpRHCCovMatrix()[ib][ib]));
	  //cout << "Single Param Variance = " << skexp->GetBANFF()->FHCpRHCCovMatrix()[ib][ib] << endl;
	  skexp->GetBANFF()->TweakSample(true,&v);
	  break;
	}
	
	// skdetfsi fhc
	else if(mode >= 800 && mode < 800+conv::kNumOfSKDetFSIParams) {
	  TVectorD v(conv::kNumOfSKDetFSIParams);
	  for (int ipf = 0; ipf < conv::kNumOfSKDetFSIParams; ipf++) {
	    v[ipf] = 1;
	  }
	  int ib = mode - 800;
	  v[ib] = 1 + (tweak * sqrt(skexp->GetSuperkDetectorFSI()->CovMatrix()[ib][ib]));
	  skexp->GetSuperkDetectorFSI()->Tweak(true,&v);
	  break;
	}
	// skdetfsi rhc
	else if(mode >= 900 && mode < 900+conv::kNumOfRHCSKDetFSIParams) {
	  TVectorD v(conv::kNumOfRHCSKDetFSIParams);
	  for (int ipf = 0; ipf < conv::kNumOfRHCSKDetFSIParams; ipf++) {
	    v[ipf] = 1;
	  }
	  int ib = mode - 900;
	  v[ib] = 1 + (tweak * sqrt(skexp->GetSuperkDetectorFSI()->RHCCovMatrix()[ib][ib]));
	  skexp->GetSuperkDetectorFSI()->TweakRHC(true,&v);
	  break;
	}
	
	cout << "Unknown mode " << mode << endl; exit(1); break;
      }//switch(mode)
    }//mode>=100
    
    skexp->CopyTweaksTo(pparm);
    
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
    bool applied_negevents_condition = skexp->Generate(ph);
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
    
    
    VALORArray * spectrum_RHC1Re = 0;
    VALORArray * spectrum_RHC1Rmu = 0;
    VALORArray * spectrum_1Re = 0;
    VALORArray * spectrum_1Rmu = 0;
    
    //store spectrum
    if(mode == 999 ) {
      //spectrum_RHC1Rmu = skobs.ErecoSpectrumVArray(kSKsRHCFCFV1RingMuLike);
      //spectrum_RHC1Re  = skobs.ErecoSpectrumVArray(kSKsRHCFCFV1RingELike);
      //spectrum_1Rmu    = skobs.ErecoSpectrumVArray(kSKsFCFV1RingMuLike);
      //spectrum_1Re     = skobs.ErecoSpectrumVArray(kSKsFCFV1RingELike);
      spectrum_RHC1Rmu = skobs->TwoDSpectrum(kSKsRHCFCFV1RingMuLike);
      spectrum_RHC1Re  = skobs->TwoDSpectrum(kSKsRHCFCFV1RingELike);
      spectrum_1Rmu    = skobs->TwoDSpectrum(kSKsFCFV1RingMuLike);
      spectrum_1Re     = skobs->TwoDSpectrum(kSKsFCFV1RingELike);
    }
    else {
      //spectrum_RHC1Rmu = skexp->CurrPredictedErecoVArray(kSKsRHCFCFV1RingMuLike);
      //spectrum_RHC1Re  = skexp->CurrPredictedErecoVArray(kSKsRHCFCFV1RingELike);
      //spectrum_1Rmu    = skexp->CurrPredictedErecoVArray(kSKsFCFV1RingMuLike);
      //spectrum_1Re     = skexp->CurrPredictedErecoVArray(kSKsFCFV1RingELike);
      spectrum_RHC1Rmu   = skexp->CurrPredictedErecoVsQ2VArray(kSKsRHCFCFV1RingMuLike);
      spectrum_RHC1Re    = skexp->CurrPredictedErecoVsQ2VArray(kSKsRHCFCFV1RingELike);
      spectrum_1Rmu      = skexp->CurrPredictedErecoVsQ2VArray(kSKsFCFV1RingMuLike);
      spectrum_1Re       = skexp->CurrPredictedErecoVsQ2VArray(kSKsFCFV1RingELike);
    }
    
    //store Nsk
    double nskint_RHC1Rmu = spectrum_RHC1Rmu->Integral();
    double nskint_RHC1Re  = spectrum_RHC1Re ->Integral();
    double nskint_1Rmu    = spectrum_1Rmu   ->Integral();
    double nskint_1Re     = spectrum_1Re    ->Integral();
    
    hNsk_RHC1Re ->Fill(nskint_RHC1Re);
    hNsk_RHC1Rmu->Fill(nskint_RHC1Rmu);
    hNsk_1Re    ->Fill(nskint_1Re);
    hNsk_1Rmu   ->Fill(nskint_1Rmu);
    
    double ratio_mulike = nskint_1Rmu / nskint_RHC1Rmu;
    double ratio_elike = nskint_1Re / nskint_RHC1Re;
    hRatio_FHCtoRHC_elike ->Fill(ratio_elike);
    hRatio_FHCtoRHC_mulike->Fill(ratio_mulike);
    
  }//end loop over experiments
  
  w.Stop();
  cout << iskip << " experiments skipped due to limits" << endl;
  
  cout << endl;
  cout << "Systematic group: " << modetitle.c_str() << endl;
  cout << endl;
  cout << "- 1Rmu: " << endl;
  cout << "   Mean = " << hNsk_1Rmu->GetMean() << endl;
  cout << "   RMS = " << hNsk_1Rmu->GetRMS() << endl;
  cout << "   RMS/Mean = " << hNsk_1Rmu->GetRMS() / hNsk_1Rmu->GetMean() << endl;
  cout << "- 1Re: " << endl;
  cout << "   Mean = " << hNsk_1Re->GetMean() << endl;
  cout << "   RMS = " << hNsk_1Re->GetRMS() << endl;
  cout << "   RMS/Mean = " << hNsk_1Re->GetRMS() / hNsk_1Re->GetMean() << endl;
  cout << "- RHC1Rmu: " << endl;
  cout << "   Average = " << hNsk_RHC1Rmu->GetMean() << endl;
  cout << "   RMS = " << hNsk_RHC1Rmu->GetRMS() << endl;
  cout << "   RMS/Mean = " << hNsk_RHC1Rmu->GetRMS() / hNsk_RHC1Rmu->GetMean() << endl;
  cout << "- RHC1Re: " << endl;
  cout << "   Average = " << hNsk_RHC1Re->GetMean() << endl;
  cout << "   RMS = " << hNsk_RHC1Re->GetRMS() << endl;
  cout << "   RMS/Mean = " << hNsk_RHC1Re->GetRMS() / hNsk_RHC1Re->GetMean() << endl;
  cout << endl;
  cout << "e-like  FHC / RHC Mean     = " << hRatio_FHCtoRHC_elike->GetMean() << endl;
  cout << "e-like  FHC / RHC RMS      = " << hRatio_FHCtoRHC_elike->GetRMS() << endl;
  cout << "e-like  FHC / RHC RMS/Mean = " << hRatio_FHCtoRHC_elike->GetRMS() / hRatio_FHCtoRHC_elike->GetMean() << endl;
  cout << "mu-like FHC / RHC Mean     = " << hRatio_FHCtoRHC_mulike->GetMean() << endl;
  cout << "mu-like FHC / RHC RMS      = " << hRatio_FHCtoRHC_mulike->GetRMS() << endl;
  cout << "mu-like FHC / RHC RMS/Mean = " << hRatio_FHCtoRHC_mulike->GetRMS() / hRatio_FHCtoRHC_mulike->GetMean() << endl;
  cout << endl;

  outfile << "mulike Mean = " << hRatio_FHCtoRHC_mulike->GetMean() << endl;
  outfile << "mulike RMS = " << hRatio_FHCtoRHC_mulike->GetRMS()  << endl;
  outfile << "mulike RMS/Mean = " << hRatio_FHCtoRHC_mulike->GetRMS() / hRatio_FHCtoRHC_mulike->GetMean()  << endl;
  outfile << "elike Mean = " << hRatio_FHCtoRHC_elike->GetMean() << endl;
  outfile << "elike RMS = " << hRatio_FHCtoRHC_elike->GetRMS()  << endl;
  outfile << "elike RMS/Mean = " << hRatio_FHCtoRHC_elike->GetRMS() / hRatio_FHCtoRHC_elike->GetMean()  << endl;
  outfile.close();

  TCanvas *cNsk_1Rmu = new TCanvas("cNsk_1Rmu","cNsk_1Rmu");
  hNsk_1Rmu->Draw();

  TCanvas *cNsk_1Re = new TCanvas("cNsk_1Re","cNsk_1Re");
  hNsk_1Re->Draw();

  TCanvas *cNsk_RHC1Rmu = new TCanvas("cNsk_RHC1Rmu","cNsk_RHC1Rmu");
  hNsk_RHC1Rmu->Draw();

  TCanvas *cNsk_RHC1Re = new TCanvas("cNsk_RHC1Re","cNsk_RHC1Re");
  hNsk_RHC1Re->Draw();

  TCanvas *cRatio_elike = new TCanvas("cRatio_elike","cRatio_elike");
  hRatio_FHCtoRHC_elike->Draw();
  hRatio_FHCtoRHC_elike->GetXaxis()->SetTitle("Ratio # of events FHC / RHC");
  hRatio_FHCtoRHC_elike->GetYaxis()->SetTitle("# of events");
  TLegend *lRatio_elike = new TLegend(0.7,0.6,0.9,0.9);
  lRatio_elike->SetFillColor(0);
  lRatio_elike->SetLineColor(0);
  lRatio_elike->SetTextFont(132);
  lRatio_elike->SetTextSize(0.06);
  lRatio_elike->AddEntry(hRatio_FHCtoRHC_elike,modetitle.c_str(),"l");
  lRatio_elike->Draw();

  TCanvas *cRatio_mulike = new TCanvas("cRatio_mulike","cRatio_mulike");
  hRatio_FHCtoRHC_mulike->Draw();
  hRatio_FHCtoRHC_mulike->GetXaxis()->SetTitle("Ratio # of events FHC / RHC");
  hRatio_FHCtoRHC_mulike->GetYaxis()->SetTitle("# of events");
  TLegend *lRatio_mulike = new TLegend(0.7,0.6,0.9,0.9);
  lRatio_mulike->SetFillColor(0);
  lRatio_mulike->SetLineColor(0);
  lRatio_mulike->SetTextFont(132);
  lRatio_mulike->SetTextSize(0.06);
  lRatio_mulike->AddEntry(hRatio_FHCtoRHC_mulike,modetitle.c_str(),"l");
  lRatio_mulike->Draw();

  // Store histograms
  TFile *fout = new TFile(Form("sk_rhctofhc_ratio-twk_sys_throw_%s_%s.root",oscparams.c_str(),modetitle.c_str()),"RECREATE");
  fout->cd();
  hNsk_RHC1Re->Write();
  hNsk_RHC1Rmu->Write();
  hNsk_1Re->Write();
  hNsk_1Rmu->Write();
  hRatio_FHCtoRHC_mulike->Write();
  hRatio_FHCtoRHC_elike->Write();
  fout->Close();

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
