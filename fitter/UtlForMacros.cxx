//____________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//____________________________________________________________________________


#include "UtlForMacros.h"

//#include "Utilities.h"
//#include "Enumerations.h"

#include "BANFF.h"
#include "SuperkDetectorFSI.h"

#include "Conventions/Opt.h"
#include "Conventions/PhysParamId.h"
#include "Conventions/PhysHypothesisId.h"
#include "Conventions/Constants.h"
#include "Utilities/UtlXML.h"

#include "TSystem.h"

#include <iostream>
#include <string>

using std::cout;
using std::cerr;
using std::endl;

using namespace valor;
using namespace valor::conv;

//____________________________________________________________________________
const char * valor::util::FormFilePath(const char * file, const char * dir)
{
  return  (new std::string(Form("%s/%s/%s", gSystem->Getenv("T2KOSCANA"), dir, file)))->c_str();
}
//____________________________________________________________________________
string valor::util::SysThrowIdToString(unsigned int tweak_systematics)
{
  string tweak_sys_str;
  switch(tweak_systematics)
    {
    case(0):
      tweak_sys_str = "None"; break;
    case(1):
      tweak_sys_str = "All"; break;
    case(2):
      tweak_sys_str = "Special low energy numu flux increase"; break;
    case(3):
      tweak_sys_str = "Tweak all but SK and FSI"; break;
    case(4):
      tweak_sys_str = "NC1pi norm = 0.5"; break;
    case(5):
      tweak_sys_str = "High energy >1.5 GeV flux increase"; break;
    case(6):
      tweak_sys_str = "BANFF"; break;
    case(7):
      tweak_sys_str = "XSec"; break;
    case(8):
      tweak_sys_str = "SK+FSI"; break;
    default:
      cerr << "SysThrowIdToString() tweak_systematics = " << tweak_systematics
	   << " not implemented" << endl;
      exit(-1);
    }
  return tweak_sys_str;
}
//____________________________________________________________________________
void valor::util::PerformSysThrow(unsigned int tweak_systematics, SuperkExpectation * skexp, bool print)
{
  if(tweak_systematics == 0) {
    return;
  }
  else if(tweak_systematics == 1) {
    // tweak all sources of systematic err (gaussian err)
    cout << "- Randomizing all systematic params in pdf..." << endl;
    skexp->TweakAllSystematics();
    // skexp->TweakSystematic("sk_energy_scale",0.0,true);
  }
  else if(tweak_systematics == 2) {
    // special tweak, increasing numu flux parameters <700MeV by 5%
    cout << "- Increasing numu flux parameters <700MeV by 5% in pdf..." << endl;
    TVectorD twkv(skexp->GetBANFF()->BestFitValues());
    twkv[0] += 0.05;
    twkv[1] += 0.05;
    twkv[2] += 0.05;
    twkv[3] += 0.05;
    skexp->GetBANFF()->Tweak(true,&twkv);
  }
  else if(tweak_systematics == 3) {
    cout << "- Randomizing all but SK systematic params in pdf..." << endl;
    skexp->TweakAllBANFFSystematics();
    skexp->TweakAllUncorrelatedXSecSystematics();
  }
  else if(tweak_systematics == 4) {
    cout << "- Setting xsec_norm_nc1piC to -0.5" << endl;
    skexp->TweakSystematic("xsec_norm_nc1piC",-0.5);
  }
  else if(tweak_systematics == 5) {
    cout << "- Increasing flux parameters >1500MeV by 15% in pdf..." << endl;
    TVectorD twkv(skexp->GetBANFF()->BestFitValues());
    twkv[6]  += 0.15;
    twkv[7]  += 0.15;
    twkv[8]  += 0.15;
    twkv[9]  += 0.15;
    twkv[10] += 0.15;
    twkv[14] += 0.15;
    twkv[15] += 0.15;
    skexp->GetBANFF()->Tweak(true,&twkv);
  }
  else if(tweak_systematics == 6) {
    cout << "- Randomizing all BANFF systematic params in pdf..." << endl;
    skexp->TweakAllBANFFSystematics();
  }
  else if(tweak_systematics == 7) {
    cout << "- Randomizing all uncorrelated cross section systematic params in pdf..." << endl;
    skexp->TweakAllUncorrelatedXSecSystematics();
  }
  else if(tweak_systematics == 8) {
    cout << "- Randomizing all SK+FSI systematic params in pdf..." << endl;
    skexp->TweakAllSKDetectorFSISystematics();
  }
  else {
    cerr << "PerformSysThrow() tweak_systematics = " << tweak_systematics
	 << " not implemented" << endl;
    exit(-1);
  }
  if(print)
    skexp->PrintSystematicVariations();
}
//____________________________________________________________________________
void valor::util::SetGlobalOptions(const char * dmsqscheme, const char * hierarchy,
				   const char * mixparamchoice, const char * mixparamoctant,
				   const char * cpvchoice)
{
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  Opt::Instance()->gOscDeltaCP           = util::CPVChoiceFromString(cpvchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  cout << *Opt::Instance() << endl;
}
//____________________________________________________________________________
void valor::util::SetOscParamOptions(double& input_mixing12_nh, double& input_mixing13_nh, double& input_mixing23_nh,
				     double& input_deltacp_nh,  double& input_dm12sq_nh,   double& input_dm23sq_nh,
				     bool& randomize_mixing12_nh, bool& randomize_mixing13_nh, bool& randomize_mixing23_nh,
				     bool& randomize_deltacp_nh,  bool& randomize_dm12sq_nh,   bool& randomize_dm23sq_nh,
				     double& input_beta_bar, bool& randomize_beta_bar,
				     double& input_mixing12_ih, double& input_mixing13_ih, double& input_mixing23_ih,
				     double& input_deltacp_ih,  double& input_dm12sq_ih,   double& input_dm23sq_ih,
				     bool& randomize_mixing12_ih, bool& randomize_mixing13_ih, bool& randomize_mixing23_ih,
				     bool& randomize_deltacp_ih,  bool& randomize_dm12sq_ih,   bool& randomize_dm23sq_ih
				     )
{
  // Check whether some osc. params need to be randomized (i.e. if == -999)
  //NH
  randomize_mixing12_nh = (TMath::Abs(input_mixing12_nh + 999) < kASmallNum);
  randomize_mixing13_nh = (TMath::Abs(input_mixing13_nh + 999) < kASmallNum);
  randomize_mixing23_nh = (TMath::Abs(input_mixing23_nh + 999) < kASmallNum);
  randomize_deltacp_nh  = (TMath::Abs(input_deltacp_nh  + 999) < kASmallNum);
  randomize_dm12sq_nh   = (TMath::Abs(input_dm12sq_nh   + 999) < kASmallNum);
  randomize_dm23sq_nh   = (TMath::Abs(input_dm23sq_nh   + 999) < kASmallNum);
  //IH
  randomize_mixing12_ih = (TMath::Abs(input_mixing12_ih + 999) < kASmallNum);
  randomize_mixing13_ih = (TMath::Abs(input_mixing13_ih + 999) < kASmallNum);
  randomize_mixing23_ih = (TMath::Abs(input_mixing23_ih + 999) < kASmallNum);
  randomize_deltacp_ih  = (TMath::Abs(input_deltacp_ih  + 999) < kASmallNum);
  randomize_dm12sq_ih   = (TMath::Abs(input_dm12sq_ih   + 999) < kASmallNum);
  randomize_dm23sq_ih   = (TMath::Abs(input_dm23sq_ih   + 999) < kASmallNum);
  //Beta bar
  randomize_beta_bar    = (TMath::Abs(input_beta_bar    + 999) < kASmallNum);

  //Set to defaults (only if enough negative (standard value is -99)
  //NH
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13_nh < -1 && input_mixing13_nh > -900) input_mixing13_nh = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_nh < -1 && input_mixing23_nh > -900) input_mixing23_nh = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_nh < -1 && input_mixing12_nh > -900) input_mixing12_nh = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13_nh < -1 && input_mixing13_nh > -900) input_mixing13_nh = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_nh < -1 && input_mixing23_nh > -900) input_mixing23_nh = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_nh < -1 && input_mixing12_nh > -900) input_mixing12_nh = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq_nh < -1 && input_dm23sq_nh > -900)  input_dm23sq_nh = conv::kGlobFitDm23sq;
  if(input_dm12sq_nh < -1 && input_dm12sq_nh > -900)  input_dm12sq_nh = conv::kGlobFitDm12sq;
  if(input_deltacp_nh < -90 && input_deltacp_nh > -900) input_deltacp_nh = conv::kGlobFitDeltaCP;
  //IH
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13_ih < -1 && input_mixing13_ih > -900) input_mixing13_ih = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_ih < -1 && input_mixing23_ih > -900) input_mixing23_ih = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_ih < -1 && input_mixing12_ih > -900) input_mixing12_ih = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13_ih < -1 && input_mixing13_ih > -900) input_mixing13_ih = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_ih < -1 && input_mixing23_ih > -900) input_mixing23_ih = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_ih < -1 && input_mixing12_ih > -900) input_mixing12_ih = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq_ih < -1 && input_dm23sq_ih > -900)  input_dm23sq_ih = conv::kGlobFitDm23sq;
  if(input_dm12sq_ih < -1 && input_dm12sq_ih > -900)  input_dm12sq_ih = conv::kGlobFitDm12sq;
  if(input_deltacp_ih < -90 && input_deltacp_ih > -900) input_deltacp_ih = conv::kGlobFitDeltaCP;
  //Beta bar
  if(input_beta_bar < -90 && input_beta_bar > -900) input_beta_bar = conv::kGlobFitBetaBar;

  cout << "   - True oscillation parameters: " << endl;
  cout << "     - Normal ordering: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = "; ((randomize_mixing12_nh) ? cout << "* randomized *" : cout << input_mixing12_nh) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = "; ((randomize_mixing13_nh) ? cout << "* randomized *" : cout << input_mixing13_nh) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = "; ((randomize_mixing23_nh) ? cout << "* randomized *" : cout << input_mixing23_nh) << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = "; ((randomize_deltacp_nh)  ? cout << "* randomized *" : cout << input_deltacp_nh)  << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = "; ((randomize_dm12sq_nh)   ? cout << "* randomized *" : cout << input_dm12sq_nh )  << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = "; ((randomize_dm23sq_nh)   ? cout << "* randomized *" : cout << input_dm23sq_nh )  << " eV^2/c^4" << endl;
  cout << "     - Inverted ordering: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = "; ((randomize_mixing12_ih) ? cout << "* randomized *" : cout << input_mixing12_ih) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = "; ((randomize_mixing13_ih) ? cout << "* randomized *" : cout << input_mixing13_ih) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = "; ((randomize_mixing23_ih) ? cout << "* randomized *" : cout << input_mixing23_ih) << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = "; ((randomize_deltacp_ih)  ? cout << "* randomized *" : cout << input_deltacp_ih)  << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = "; ((randomize_dm12sq_ih)   ? cout << "* randomized *" : cout << input_dm12sq_ih )  << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = "; ((randomize_dm23sq_ih)   ? cout << "* randomized *" : cout << input_dm23sq_ih )  << " eV^2/c^4" << endl;
  cout << "     - Other: " << endl;
  cout << "       > " << util::AsString(kPhypBetaBar)  << " = "; ((randomize_beta_bar)    ? cout << "* randomized *" : cout << input_beta_bar) << endl;
}
//____________________________________________________________________________
void valor::util::RandomizeOscParams(double& input_mixing12_nh, double& input_mixing13_nh, double& input_mixing23_nh,
				     double& input_deltacp_nh,  double& input_dm12sq_nh,   double& input_dm23sq_nh,
				     bool randomize_mixing12_nh, bool randomize_mixing13_nh, bool randomize_mixing23_nh,
				     bool randomize_deltacp_nh,  bool randomize_dm12sq_nh,   bool randomize_dm23sq_nh,
				     double& input_beta_bar, bool randomize_beta_bar,
				     double& input_mixing12_ih, double& input_mixing13_ih, double& input_mixing23_ih,
				     double& input_deltacp_ih,  double& input_dm12sq_ih,   double& input_dm23sq_ih,
				     bool randomize_mixing12_ih, bool randomize_mixing13_ih, bool randomize_mixing23_ih,
				     bool randomize_deltacp_ih,  bool randomize_dm12sq_ih,   bool randomize_dm23sq_ih
				     )
{
  if (randomize_mixing12_nh) { input_mixing12_nh = util::RandomizeMixing12(); }
  if (randomize_mixing13_nh) { input_mixing13_nh = util::RandomizeMixing13(); }
  if (randomize_mixing23_nh) { input_mixing23_nh = util::RandomizeMixing23(); }
  if (randomize_deltacp_nh ) { input_deltacp_nh  = util::RandomizeCPV();      }
  if (randomize_dm12sq_nh  ) { input_dm12sq_nh   = util::RandomizeDm12sq();   }
  if (randomize_dm23sq_nh  ) { input_dm23sq_nh   = util::RandomizeDm23sq();   }

  if (randomize_beta_bar   ) { input_beta_bar    = util::RandomizeBetaBar();  }

  if (randomize_mixing12_ih) { input_mixing12_ih = util::RandomizeMixing12(); }
  if (randomize_mixing13_ih) { input_mixing13_ih = util::RandomizeMixing13(); }
  if (randomize_mixing23_ih) { input_mixing23_ih = util::RandomizeMixing23(); }
  if (randomize_deltacp_ih ) { input_deltacp_ih  = util::RandomizeCPV();      }
  if (randomize_dm12sq_ih  ) { input_dm12sq_ih   = util::RandomizeDm12sq();   }
  if (randomize_dm23sq_ih  ) { input_dm23sq_ih   = util::RandomizeDm23sq();   }
}
//____________________________________________________________________________
void valor::util::SetOscParams(PhysHypothesis& ph,
			       double input_mixing12_used, double input_mixing13_used, double input_mixing23_used,
			       double input_deltacp_used,  double input_dm12sq_used,   double input_dm23sq_used,
			       double input_beta_bar_used,
			       bool print)
{
  if(ph.Hypothesis() == kHypNuMuDisapp2F) {
       ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
  }
  else if(ph.Hypothesis() == kHypOsc3F) {  
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_deltacp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
  }
  else if(ph.Hypothesis() == kHypOsc3FBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_deltacp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypMixing12Bar,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13Bar,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23Bar,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPVBar    ,   input_deltacp_used);
    ph.SetCurrValue(kPhypDm12sqBar,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sqBar,    input_dm23sq_used);
  }
  else if(ph.Hypothesis() == kHypOsc3FNormBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_deltacp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypBetaBar,    input_beta_bar_used);
  }
  else if(ph.Hypothesis() == kHypOsc3FBarNormBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_deltacp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypMixing12Bar,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13Bar,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23Bar,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPVBar    ,   input_deltacp_used);
    ph.SetCurrValue(kPhypDm12sqBar,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sqBar,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypBetaBar,    input_beta_bar_used);
  }
  else {
    cout << "SetOscParams() " << util::AsString(ph.Hypothesis()) << " not implemented" << endl;
    exit(1);
  }
  if(print)
    cout << "- Physics hypothesis set: " << ph << endl;
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
