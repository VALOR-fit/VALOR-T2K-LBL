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
// 800+ number of SK params tweak last skdetfsi parameter
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

void test_sk_ereco_spectrum_twk_sys_throw_joint_2015_2d
(
 int          mode                     = 13, //13,                                          /* options: see above */
 bool         update_rootfile          = true,                                        /* "update" ROOT output file? If false, will "recreate" */
 const char * sk_sample                = "1Re",                                   /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
 const char * hypothesis               = "3f_osc",                                    /* options: "2f_numu_disap", "3f_osc", "3f_oscbar" */
 
 // Asimov data set A
 double       input_mixing12       = 0.304,//-99,    
 double       input_mixing13       = 0.0217, //-99,    
 double       input_mixing23       = 0.528,    
 double       input_deltacp        = -1.601,   
 double       input_dm12sq         = 7.53E-5,   
 double       input_dm23sq         = 2.509E-3,     
 double       input_betabar        = 1,

 // // Asimov data set B
 // double       input_mixing12       = 0.304,    
 // double       input_mixing13       = 0.0217,  
 // double       input_mixing23       = 0.45,    
 // double       input_deltacp        = 0.,   
 // double       input_dm12sq         = 7.53E-5,   
 // double       input_dm23sq         = 2.509E-3,     
 // double       input_betabar        = 1,
 
 bool         save_toy_expts           = false,                                       /* save the ereco spectra? */
 bool         save_throws              = false,                                        /* save the thrown parameter values*/
 
 //
 // NEW CONFIGURATION FILES (4/6/16 BANFF RELEASE W/ FIT BIAS AND ADDITIONAL UNCERTAINTY ON RUN 1-7 FLUX)
 //
 const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17.xml", // Run 5-7 flux tuning for RHC
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_Run17.xml", // Run 1-6 flux tuning

 //
 // NEW CONFIGURATION FILES (14/3/16 BANFF RELEASE W/ FIT BIAS)
 //
 //const char * config_file      = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FSIOnly.xml",
 //const char * config_file      = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_SKDetOnly.xml",

 const int    nexp                     = 10000,                                          /* number of experiments to make */
 const char * dmsqscheme               = "auto",
 const char * mixparamchoice           = "sinsqtheta",
 const char * hierarchy                = "normal",
 const char * mixparamoctant           = "two",       
 bool         reject_unphysical_throws = true    ,
 bool         use_new_definition_of_dn_over_n = false,
 
 const bool   doMargMH             = false,      // marginalize Mass Hierarchy
 const double PriorProbNH          = 0.5, //0.684, // Mach3  // Prior for MH --> give the prob. of NH as input
 TString file_priors = "", // Set to "" if no user priors are used (VALOR)
 const char*  tag                      = "nom",

 //TString       randomised_osc          = "ss23-ss13-dcp-ss12-dmsq23-dmsq12"
 TString       randomised_osc          = ""
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

  TString str_sk_sample(sk_sample);
  bool RHCFlag = str_sk_sample.Contains("RHC");

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
      else if(mode >= 800 && mode < 900) {
        modetitle = Form("skdetfsi%03d",mode-800); break;
      }
      cout << "Unknown mode " << mode << endl; exit(1); break;
  }//switch mode
  

  //
  
  TString config_str(config_file);
  bool isprefit = config_str.Contains("prefit");
  if(isprefit) {
    modetitle += "prefit";
  }

  //

  TH1D * hErecoNominal = 0;
  TH1D * hNsk; 
  TH1D* h_nominal;
  TH1D** h_bin = new TH1D*[util::NumErecoBins(s)];
  TH1D** h_bin_coarse = new TH1D*[util::NumErecoBins(s)];

  // RHC1Re
  if(s==conv::kSKsRHCFCFV1RingELike){
    hNsk = new TH1D(Form("%s_%s_nsk",sname,modetitle.c_str()),Form("%s_%s_nsk",sname,modetitle.c_str()),100,0,10);
  }
  // 1Re
  else if(s==conv::kSKsFCFV1RingELike){
    hNsk = new TH1D(Form("%s_%s_nsk",sname,modetitle.c_str()),Form("%s_%s_nsk",sname,modetitle.c_str()),100,0,100);
  }
  // RHC1Rmu
  else if(s==conv::kSKsRHCFCFV1RingMuLike){
    hNsk = new TH1D(Form("%s_%s_nsk",sname,modetitle.c_str()),Form("%s_%s_nsk",sname,modetitle.c_str()),100,0,100);
  }
  // 1Rmu
  else if(s==conv::kSKsFCFV1RingMuLike){
    hNsk = new TH1D(Form("%s_%s_nsk",sname,modetitle.c_str()),Form("%s_%s_nsk",sname,modetitle.c_str()),100,0,300);      
  } 

  for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){    
    char name[200];
    sprintf(name,"h_bin_%d",ibin);

    // RHC1Re
    if(s==conv::kSKsRHCFCFV1RingELike){
      h_bin[ibin] = new TH1D(name,name,20000,0,3);
      sprintf(name,"h_bin_coarse_%d",ibin);
      h_bin_coarse[ibin] = new TH1D(name,name,200,0,3);
    }
    // 1Re
    else if(s==conv::kSKsFCFV1RingELike){
      h_bin[ibin] = new TH1D(name,name,20000,0,10);
      sprintf(name,"h_bin_coarse_%d",ibin);
      h_bin_coarse[ibin] = new TH1D(name,name,200,0,10);
    }
    // RHC1Rmu
    else if(s==conv::kSKsRHCFCFV1RingMuLike){
      h_bin[ibin] = new TH1D(name,name,20000,0,20);
      sprintf(name,"h_bin_coarse_%d",ibin);
      h_bin_coarse[ibin] = new TH1D(name,name,200,0,20);
    }
    // 1Rmu
    else if(s==conv::kSKsFCFV1RingMuLike){
      h_bin[ibin] = new TH1D(name,name,20000,0,50);
      sprintf(name,"h_bin_coarse_%d",ibin);
      h_bin_coarse[ibin] = new TH1D(name,name,200,0,50);
    } 
  }

  TH1D* hErecoAverage = new TH1D("hErecoAverage","hErecoAverage",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1D* h_up = new TH1D("h_up","h_up",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1D* h_down = new TH1D("h_down","h_down",util::NumErecoBins(s),util::ErecoBinEdges(s));

  TH1D* hErecoMedian = new TH1D("hErecoMedian","hErecoMedian",util::NumErecoBins(s),util::ErecoBinEdges(s)); 
  TH1D* h_up_Median = new TH1D("h_up_Median","h_up_Median",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1D* h_down_Median = new TH1D("h_down_Median","h_down_Median",util::NumErecoBins(s),util::ErecoBinEdges(s));

  //

  int syst_nsigma = util::GetSystNSigma();

  cout << endl;
  cout << "Systematic mode (" << mode << ") --> " << modetitle << endl;
  cout << endl;

  // Analysis configuration
  TString config  = TString::Format("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  string sconfig(config.Data());  
  T2KOpt::Instance()->ReadXmlPath(sconfig);    

  // Check whether some osc. params need to be randomized

  bool randomize_mixing12   = randomised_osc.Contains("ss12");
  bool randomize_mixing13   = randomised_osc.Contains("ss13");
  bool randomize_mixing23   = randomised_osc.Contains("ss23");
  bool randomize_deltacp    = randomised_osc.Contains("dcp");
  bool randomize_dm12sq     = randomised_osc.Contains("dmsq12");
  bool randomize_dm23sq     = randomised_osc.Contains("dmsq23");
  
  // Set Physics Hypothesis
  PhysHypothesis ph(hypothesis);

  // Set oscillation parameter values and priors

  // util::OscParams::Instance()->SetValues(input_mixing12, input_mixing13, input_mixing23, input_deltacp, input_dm12sq, input_dm23sq,
  //  input_betabar,
  //  input_mixing12, input_mixing13, input_mixing23, input_deltacp, input_dm12sq, input_dm23sq,
  //  PriorProbNH, doMargMH
  //  );

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
  
  if(randomize_mixing23){ 
    util::OscParams::Instance()->SetPrior(conv::kPhypMixing23,    PriorMixing23, PriorParMixing23[0], PriorParMixing23[1], MixModeMixing23);      
    util::OscParams::Instance()->SetRandomize(conv::kPhypMixing23,1);
  }
  if(randomize_dm23sq){   
    util::OscParams::Instance()->SetPrior(conv::kPhypDm23sq,    PriorDm23sq,   PriorParDm23sq[0],   PriorParDm23sq[1],   ""              );  
    util::OscParams::Instance()->SetRandomize(conv::kPhypDm23sq,1);
  }
  if(randomize_mixing13){ 
    util::OscParams::Instance()->SetPrior(conv::kPhypMixing13,    PriorMixing13, PriorParMixing13[0], PriorParMixing13[1], MixModeMixing13); // Convert to sinsqth
    util::OscParams::Instance()->SetRandomize(conv::kPhypMixing13,1);
  }
  if(randomize_deltacp){  
    util::OscParams::Instance()->SetPrior(conv::kPhypCPV,      PriorDeltaCP,  PriorParDeltaCP[0],  PriorParDeltaCP[1],  ""              ); 
    util::OscParams::Instance()->SetRandomize(conv::kPhypCPV,1);
  }
  if(randomize_mixing12){ 
    util::OscParams::Instance()->SetPrior(conv::kPhypMixing12, PriorMixing12, PriorParMixing12[0], PriorParMixing12[1], MixModeMixing12); // Convert to sinsqth
    util::OscParams::Instance()->SetRandomize(conv::kPhypMixing12,1);
  }
  if(randomize_dm12sq){   
    util::OscParams::Instance()->SetPrior(conv::kPhypDm12sq,   PriorDm12sq,   PriorParDm12sq[0],   PriorParDm12sq[1],   ""              ); 
    util::OscParams::Instance()->SetRandomize(conv::kPhypDm12sq,1);
  }

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
  //SuperkExpectation *skexp;
  SuperkExpectation *skexp = new SuperkExpectation;
  skexp->Read(sconfig);
  
  const SuperkDetectorFSI * skdet = skexp->GetSuperkDetectorFSI();
  int n_total_skdetfsi_parameters = kNumOfSKDetFSIParams;
  
  SuperkObservation skobs;
  
  //string oscparams = util::AsString(ph.Hypothesis(),input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq,input_betabar);
  string oscparams = util::AsString(ph.Hypothesis(),
    input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq,input_betabar,-999,-999);

  cout << "input_deltacp = " << input_deltacp << endl;
  cout << oscparams << endl;
  


  // Output histograms
  //TH1D * hErecoNominal = 0;
  VALORArray * obsEreco[nexp];
  for (int i = 0 ; i < nexp; i++) {obsEreco[i]=0;}
  //TH1D * hNsk = new TH1D(Form("%s_%s_nsk",sname,modetitle.c_str()),Form("%s_%s_nsk",sname,modetitle.c_str()),100,0,10);
  
  int nbanff;
  if(RHCFlag) nbanff = conv::kNumOfRHCParamsInBANFFfit;
  else        nbanff = conv::kNumOfParamsInBANFFfit;

  int nxsecuncorr;
  if(RHCFlag) nxsecuncorr = conv::kNumOfRHCXSecSyst;
  else        nxsecuncorr = conv::kNumOfXSecSyst;

  const int sys = nbanff + n_total_skdetfsi_parameters + nxsecuncorr;

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

  TFile *fosctoyhisto = new TFile(Form("sk_osc_toy_histo_%s-twk_sys_throw_%s_%s.root",sname,oscparams.c_str(),tag),rootmode.c_str(),"recreate");
  TFile f(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s_%s.root",sname,oscparams.c_str()),rootmode.c_str(),tag);
  
  ofstream * outfilep = 0;
  
  if(!update_rootfile){
    outfilep = new ofstream(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s_%s.txt",sname,oscparams.c_str(),tag),ios_base::trunc);
  }
  else {
    outfilep = new ofstream(Form("sk_ereco_spectrum_%s-twk_sys_throw_%s_%s.txt",sname,oscparams.c_str(),tag),ios_base::app);
  }
  
  ofstream & outfile = *outfilep;
  
  //remove the [] so can use it in SaveAs
  while(oscparams.find_first_of("[]") != string::npos) {
    oscparams.erase(oscparams.find_first_of("[]"), 1);
  }

  // Get the nominal
  skexp->Generate(ph);
  hErecoNominal = (TH1D*) skexp->CurrPredictedErecoVsQ2Spectrum(s)->ProjectionY()->Clone("hErecoNominal");
  //hErecoNominal = skexp->CurrPredictedErecoSpectrum(s);
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
      
      // All expect oscillations
      if(mode == 13) {
        skexp->TweakAllSystematics();
      }
      
      // Energy Scale
      if(mode == 16){
        skexp->GetSuperkDetectorFSI()->TweakEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
      }
      
      //oscillations parameters
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

        case 300:
          skexp->TweakSystematic("sk_e_scale",tweak); break;
    
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
          else if(mode >= 800 && mode < 900) {
            TVectorD v(n_total_skdetfsi_parameters);
            for (int ipf = 0; ipf < n_total_skdetfsi_parameters; ipf++) {
              v[ipf] = 1;
            }
            int ib = mode - 800;
            PhysParamId_t p = skdet->WhichParam(ib);
            v[ib] = 1 + (tweak * skdet->OneSigmaError(p));
            skdet->Tweak(true,&v);
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
    
    
    if(mode==425) hMEC->Fill(v[25]);

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
    
    VALORArray * spectrum = 0;

    //store spectrum
    if (save_toy_expts) {
      if(mode == 999 ) {
        //spectrum = new VALORArray(*skobs.ErecoSpectrumVArray(s));
        spectrum = new VALORArray(*skobs.TwoDSpectrum(s));
      }
      else {
        spectrum = new VALORArray(*skexp->CurrPredictedErecoVsQ2VArray(s));
        //spectrum = new VALORArray(*skexp->CurrPredictedErecoVArray(s));
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
        //spectrum = skexp->CurrPredictedErecoVArray(s);
        spectrum = skexp->CurrPredictedErecoVsQ2VArray(s);
      }
    }
  
    //store Nsk
    double nskint = spectrum->Integral();
 
    //cout << "nskint = " << nskint << endl;

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

      //h_2d->Fill(h_2d->GetXaxis()->GetBinCenter(ibin+1),N);

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
            hSysThrown[ip]->Fill(skexp->GetBANFF()->RHCCurrentValue(p));
          else 
            hSysThrown[ip]->Fill(skexp->GetBANFF()->CurrentValue(p));
        }
        else if(util::IsXSec(p)){
          hSysThrown[ip]->Fill(skexp->GetNuXSec()->CurrentValue(p));
        }
        if(util::IsSKDetFSI(p)){
          hSysThrown[ip]->Fill(skexp->GetSuperkDetectorFSI()->CurrentValue(p));
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


  // Store extracted osc param toy histos
  fosctoyhisto->cd();
  util::OscParams::Instance()->Write(fosctoyhisto); // store the toy osc param histograms
  
  fosctoyhisto->Close();
  f->Close();

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
