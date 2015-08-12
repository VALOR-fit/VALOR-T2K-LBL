//____________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//____________________________________________________________________________


#include "UtlForMacros.h"

#include "SuperkExpectation.h"
#include "SuperkMC.h"
#include "BANFF.h"
#include "SuperkDetectorFSI.h"

#include "Conventions/Opt.h"
#include "Conventions/PhysParamId.h"
#include "Conventions/PhysHypothesisId.h"
#include "Conventions/Constants.h"
#include "FitCore/PhysHypothesis.h"
#include "Utilities/UtlXML.h"

#include "TSystem.h"

#include <iostream>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

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
  cout << *Opt::Instance();
}
//____________________________________________________________________________
util::OscParams valor::util::SetOscParamOptions(double input_mixing12_nh, double input_mixing13_nh, double input_mixing23_nh,
						double input_dcp_nh,      double input_dm12sq_nh,   double input_dm23sq_nh,
						double input_beta_bar,
						double input_mixing12_ih, double input_mixing13_ih, double input_mixing23_ih,
						double input_dcp_ih,      double input_dm12sq_ih,   double input_dm23sq_ih)
{
  util::OscParams inputs;
  
  // Check whether some osc. params need to be randomized (i.e. if == -999)
  //NH
  inputs.randomize_mixing12_nh = (TMath::Abs(input_mixing12_nh + 999) < kASmallNum);
  inputs.randomize_mixing13_nh = (TMath::Abs(input_mixing13_nh + 999) < kASmallNum);
  inputs.randomize_mixing23_nh = (TMath::Abs(input_mixing23_nh + 999) < kASmallNum);
  inputs.randomize_dcp_nh      = (TMath::Abs(input_dcp_nh      + 999) < kASmallNum);
  inputs.randomize_dm12sq_nh   = (TMath::Abs(input_dm12sq_nh   + 999) < kASmallNum);
  inputs.randomize_dm23sq_nh   = (TMath::Abs(input_dm23sq_nh   + 999) < kASmallNum);
  //IH
  inputs.randomize_mixing12_ih = (TMath::Abs(input_mixing12_ih + 999) < kASmallNum);
  inputs.randomize_mixing13_ih = (TMath::Abs(input_mixing13_ih + 999) < kASmallNum);
  inputs.randomize_mixing23_ih = (TMath::Abs(input_mixing23_ih + 999) < kASmallNum);
  inputs.randomize_dcp_ih      = (TMath::Abs(input_dcp_ih      + 999) < kASmallNum);
  inputs.randomize_dm12sq_ih   = (TMath::Abs(input_dm12sq_ih   + 999) < kASmallNum);
  inputs.randomize_dm23sq_ih   = (TMath::Abs(input_dm23sq_ih   + 999) < kASmallNum);
  //Beta bar
  inputs.randomize_beta_bar    = (TMath::Abs(input_beta_bar    + 999) < kASmallNum);

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
  if(input_dm23sq_nh < -1  && input_dm23sq_nh > -900) input_dm23sq_nh = conv::kGlobFitDm23sq;
  if(input_dm12sq_nh < -1  && input_dm12sq_nh > -900) input_dm12sq_nh = conv::kGlobFitDm12sq;
  if(input_dcp_nh    < -90 && input_dcp_nh    > -900) input_dcp_nh    = conv::kGlobFitDeltaCP;
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
  if(input_dm23sq_ih < -1  && input_dm23sq_ih > -900) input_dm23sq_ih = conv::kGlobFitDm23sq;
  if(input_dm12sq_ih < -1  && input_dm12sq_ih > -900) input_dm12sq_ih = conv::kGlobFitDm12sq;
  if(input_dcp_ih    < -90 && input_dcp_ih    > -900) input_dcp_ih = conv::kGlobFitDeltaCP;
  //Beta bar
  if(input_beta_bar < -90 && input_beta_bar > -900) input_beta_bar = conv::kGlobFitBetaBar;

  //Set the input values in the OscParams
  //NH
  inputs.mixing12_nh = input_mixing12_nh;
  inputs.mixing13_nh = input_mixing13_nh;
  inputs.mixing23_nh = input_mixing23_nh;
  inputs.dcp_nh      = input_dcp_nh;
  inputs.dm12sq_nh   = input_dm12sq_nh;
  inputs.dm23sq_nh   = input_dm23sq_nh;
  //IH
  inputs.mixing12_ih = input_mixing12_ih;
  inputs.mixing13_ih = input_mixing13_ih;
  inputs.mixing23_ih = input_mixing23_ih;
  inputs.dcp_ih      = input_dcp_ih;
  inputs.dm12sq_ih   = input_dm12sq_ih;
  inputs.dm23sq_ih   = input_dm23sq_ih;
  //Other
  inputs.beta_bar    = input_beta_bar;
  
  cout << "   - True oscillation parameters: " << endl;
  cout << "     - Normal ordering: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = "; ((inputs.randomize_mixing12_nh) ? cout << "* randomized *" : cout << input_mixing12_nh) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = "; ((inputs.randomize_mixing13_nh) ? cout << "* randomized *" : cout << input_mixing13_nh) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = "; ((inputs.randomize_mixing23_nh) ? cout << "* randomized *" : cout << input_mixing23_nh) << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = "; ((inputs.randomize_dcp_nh)      ? cout << "* randomized *" : cout << input_dcp_nh) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = "; ((inputs.randomize_dm12sq_nh)   ? cout << "* randomized *" : cout << input_dm12sq_nh )  << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = "; ((inputs.randomize_dm23sq_nh)   ? cout << "* randomized *" : cout << input_dm23sq_nh )  << " eV^2/c^4" << endl;
  cout << "     - Inverted ordering: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = "; ((inputs.randomize_mixing12_ih) ? cout << "* randomized *" : cout << input_mixing12_ih) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = "; ((inputs.randomize_mixing13_ih) ? cout << "* randomized *" : cout << input_mixing13_ih) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = "; ((inputs.randomize_mixing23_ih) ? cout << "* randomized *" : cout << input_mixing23_ih) << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = "; ((inputs.randomize_dcp_ih)      ? cout << "* randomized *" : cout << input_dcp_ih) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = "; ((inputs.randomize_dm12sq_ih)   ? cout << "* randomized *" : cout << input_dm12sq_ih )  << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = "; ((inputs.randomize_dm23sq_ih)   ? cout << "* randomized *" : cout << input_dm23sq_ih )  << " eV^2/c^4" << endl;
  cout << "     - Other: " << endl;
  cout << "       > " << util::AsString(kPhypBetaBar)  << " = "; ((inputs.randomize_beta_bar)    ? cout << "* randomized *" : cout << input_beta_bar) << endl;

  return inputs;
}
//____________________________________________________________________________
void valor::util::RandomizeOscParams(OscParams & inputs)
{
  if (inputs.randomize_mixing12_nh) { inputs.mixing12_nh = util::RandomizeMixing12(); }
  if (inputs.randomize_mixing13_nh) { inputs.mixing13_nh = util::RandomizeMixing13(); }
  if (inputs.randomize_mixing23_nh) { inputs.mixing23_nh = util::RandomizeMixing23(); }
  if (inputs.randomize_dcp_nh     ) { inputs.dcp_nh      = util::RandomizeCPV();      }
  if (inputs.randomize_dm12sq_nh  ) { inputs.dm12sq_nh   = util::RandomizeDm12sq();   }
  if (inputs.randomize_dm23sq_nh  ) { inputs.dm23sq_nh   = util::RandomizeDm23sq();   }

  if (inputs.randomize_beta_bar   ) { inputs.beta_bar    = util::RandomizeBetaBar();  }

  if (inputs.randomize_mixing12_ih) { inputs.mixing12_ih = util::RandomizeMixing12(); }
  if (inputs.randomize_mixing13_ih) { inputs.mixing13_ih = util::RandomizeMixing13(); }
  if (inputs.randomize_mixing23_ih) { inputs.mixing23_ih = util::RandomizeMixing23(); }
  if (inputs.randomize_dcp_ih     ) { inputs.dcp_ih      = util::RandomizeCPV();      }
  if (inputs.randomize_dm12sq_ih  ) { inputs.dm12sq_ih   = util::RandomizeDm12sq();   }
  if (inputs.randomize_dm23sq_ih  ) { inputs.dm23sq_ih   = util::RandomizeDm23sq();   }
}
//____________________________________________________________________________
void valor::util::SetOscParams(PhysHypothesis& ph, OscParams inputs, int iMH, bool print)
{
  double input_mixing12_used, input_mixing13_used, input_mixing23_used;
  double input_dcp_used,      input_dm12sq_used,   input_dm23sq_used;
  double input_beta_bar_used;

  //Get the values to use (based on MH)
  if(iMH == 0) {
    //NH
    input_mixing12_used = inputs.mixing12_nh;
    input_mixing13_used = inputs.mixing13_nh;
    input_mixing23_used = inputs.mixing23_nh;
    input_dcp_used      = inputs.dcp_nh;
    input_dm12sq_used   = inputs.dm12sq_nh;
    input_dm23sq_used   = inputs.dm23sq_nh;
  }
  else if(iMH == 1) {
    //IH
    input_mixing12_used = inputs.mixing12_ih;
    input_mixing13_used = inputs.mixing13_ih;
    input_mixing23_used = inputs.mixing23_ih;
    input_dcp_used      = inputs.dcp_ih;
    input_dm12sq_used   = inputs.dm12sq_ih;
    input_dm23sq_used   = inputs.dm23sq_ih;
  }
  else {
    cerr << "SetOscParams has unknown value of iMH: " << iMH << endl;
    exit(1);
  }
  //Other
  input_beta_bar_used = inputs.beta_bar;
  
  //Update the PhysHypothesis
  if(ph.Hypothesis() == kHypNuMuDisapp2F) {
       ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
       ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
  }
  else if(ph.Hypothesis() == kHypOsc3F) {  
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
  }
  else if(ph.Hypothesis() == kHypOsc3FBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypMixing12Bar,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13Bar,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23Bar,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPVBar    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sqBar,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sqBar,    input_dm23sq_used);
  }
  else if(ph.Hypothesis() == kHypOsc3FNormBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypBetaBar,    input_beta_bar_used);
  }
  else if(ph.Hypothesis() == kHypOsc3FBarNormBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypMixing12Bar,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13Bar,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23Bar,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPVBar    ,   input_dcp_used);
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
int valor::util::SetMassHierarchy(SuperkExpectation * skexp, bool doMargMH, double PriorProbNH, TH1D * hToy_MH)
{
  conv::NuMassHierarchy_t hier;
     
  //  - iMH=0 --> kNmhNormal
  //  - iMH=1 --> kNmhInverted	
  int iMH = conv::kBadResult;
	
  // Don't marginalize Mass Hierarchy. Use the one already set
  if(!doMargMH) {
    hier = Opt::Instance()->gMassHierarchyAssumed;
    if(hier == kNmhNormal)
      iMH = 0;
    else if(hier == kNmhInverted)
      iMH = 1; 
    else {
      cerr << "Mass Hierarchy is not undefined" << endl;
      exit(1);
    }
  }//!doMargMH     
  // Extract randomly the Mass Hierarchy 
  else {
    iMH = util::RandomizeMH(PriorProbNH);
  }

  //save the choice in a histogram
  if(hToy_MH)
    hToy_MH->Fill(iMH+conv::kASmallNum);

  //convert to NuMassHierarchy_t
  if(iMH == 0)
    hier = kNmhNormal;
  else if(iMH == 1)
    hier = kNmhInverted;
  else {
    cerr << "util::SetMassHierarchy() Unknown value of iMH " << iMH << endl;
    exit(1);
  }

  //set the hierarchy
  Opt::Instance()->gMassHierarchyAssumed = hier;
  skexp->GetSuperkMC()->SetMassHierarchy(hier);
  skexp->GetSuperkMC()->SetMassHierarchyBar(hier);

  return iMH;
}
//____________________________________________________________________________
void GetPriorHists(const char * file_priors,
		   TH2D ** hPrior_th13_dcp, TH2D ** hPrior_23,
		   TH2D ** hToy_th13_dcp, TH2D ** hToy_23)
{
  TString name_prior_th13_dcp = "hPrior_th13_dcp_";
  TString name_prior_23 = "hPrior_L23";
  TString name_mh[2];
  name_mh[0] = "nh";
  name_mh[1] = "ih";
    
  cout << "Taking priors from: " << file_priors << endl;
  TFile *fPrior = new TFile(file_priors,"READ");
  if(!fPrior->IsOpen() || !fPrior){
    cerr << "Error: the file " << file_priors << " is not open!" << endl;
      exit(1);
  }

  for(int im = 0; im < 2; im++) {
    //theta13 vs deltacp
    string name = Form("%s%s", name_prior_th13_dcp.Data(), name_mh[im].Data());
    fPrior->GetObject(name.c_str(), hPrior_th13_dcp[im]);
    if(hPrior_th13_dcp[im]->IsZombie() || !hPrior_th13_dcp[im]){
      cerr << "Error: the histogram " << name << " hasn't been found" << endl;
      exit(1);
    }
    //theta23 vs dm32sq
    name = Form("%s%s", name_prior_23.Data(), name_mh[im].Data());
    fPrior->GetObject(name.c_str(), hPrior_23[im]);
    if(hPrior_23[im]->IsZombie() || !hPrior_23[im]){
      cerr << "Error: the histogram " << name << " hasn't been found" << endl;
      exit(1);
    }

    //setup the toys histograms too
    hToy_th13_dcp[im] = new TH2D(*hPrior_th13_dcp[im]);
    hToy_th13_dcp[im]->Reset();    
    hToy_23[im] = new TH2D(*hPrior_23[im]);
    hToy_23[im]->Reset();
  }//im
  
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
