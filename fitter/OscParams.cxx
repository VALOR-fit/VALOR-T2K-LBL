//____________________________________________________________________________
/*
  \copyright (c) 2010-2014, The VALOR group.

  >> ChangeLog:

*/
//____________________________________________________________________________

#include <cstdlib>

#include "OscParams.h"

#include "UtlForMacros.h"
#include "T2KPhysParamId.h"
#include "Conventions/Opt.h"
#include "Conventions/PhysHypothesisId.h"
#include "Conventions/Constants.h"
#include "Utilities/UtlNumStat.h"
#include "Utilities/Rnd.h"
#include <TMath.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using namespace valor;
using namespace valor::conv;
  
//____________________________________________________________________________
namespace valor {
  namespace util {
    ostream & operator << (ostream & stream, const OscParams & oscparams)
    {
      oscparams.Print(stream);
      return stream;
    }
  }
}
//____________________________________________________________________________
valor::util::OscParams * valor::util::OscParams::fInstance = 0;
//____________________________________________________________________________
valor::util::OscParams::OscParams()
{
  fInstance = 0;

  this->Init();
}
//____________________________________________________________________________
valor::util::OscParams::~OscParams()
{
  fInstance = 0;

  this->Cleanup();
}
//____________________________________________________________________________
valor::util::OscParams * valor::util::OscParams::Instance()
{
  if(fInstance == 0) {
    static valor::util::OscParams::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new OscParams;
  }
  return fInstance;
}
//____________________________________________________________________________
void valor::util::OscParams::Init()
{
  mixmode_mixing12 = "";
  mixmode_mixing13 = "";
  mixmode_mixing23 = "";

  priormode_mixing12 = "";
  priormode_mixing13 = "";
  priormode_mixing23 = "";
  priormode_dcp      = "";
  priormode_dm12sq   = "";
  priormode_dm23sq   = "";

  for(int i=0;i<2;i++){
    params_mixing12[i] = conv::kBadResult;
    params_mixing13[i] = conv::kBadResult;
    params_mixing23[i] = conv::kBadResult;
    params_dcp[i]      = conv::kBadResult;
    params_dm12sq[i]   = conv::kBadResult;
    params_dm23sq[i]   = conv::kBadResult;
  }

  hPrior_th13_dcp_nh = 0;
  hPrior_th13_dcp_ih = 0;
  hPrior_23_nh = 0;
  hPrior_23_ih = 0;
  hToy_th13_dcp_nh = 0;
  hToy_th13_dcp_ih = 0;
  hToy_23_nh = 0;
  hToy_23_ih = 0;
  hToy_MH = 0;
}
//____________________________________________________________________________
void valor::util::OscParams::Cleanup()
{
  mixmode_mixing12 = "";
  mixmode_mixing13 = "";
  mixmode_mixing23 = "";

  priormode_mixing12 = "";
  priormode_mixing13 = "";
  priormode_mixing23 = "";
  priormode_dcp      = "";
  priormode_dm12sq   = "";
  priormode_dm23sq   = "";

  for(int i=0;i<2;i++){
    params_mixing12[i] = conv::kBadResult;
    params_mixing13[i] = conv::kBadResult;
    params_mixing23[i] = conv::kBadResult;
    params_dcp[i]      = conv::kBadResult;
    params_dm12sq[i]   = conv::kBadResult;
    params_dm23sq[i]   = conv::kBadResult;
  }
  
  if(hPrior_th13_dcp_nh) delete hPrior_th13_dcp_nh;
  if(hPrior_th13_dcp_ih) delete hPrior_th13_dcp_ih;
  if(hPrior_23_nh) delete hPrior_23_nh;
  if(hPrior_23_ih) delete hPrior_23_ih;
  if(hToy_th13_dcp_nh) delete hToy_th13_dcp_nh;
  if(hToy_th13_dcp_ih) delete hToy_th13_dcp_ih;
  if(hToy_23_nh) delete hToy_23_nh;
  if(hToy_23_ih) delete hToy_23_ih;
  if(hToy_MH) delete hToy_MH;
}
//____________________________________________________________________________

void valor::util::OscParams::SetPrior(conv::PhysParamId_t physpar,string mode,double par1,double par2,string mixmode)
{
  if(physpar==conv::kPhypMixing13){
    priormode_mixing13 = mode;
    mixmode_mixing13   = mixmode;
    params_mixing13[0] = par1;
    params_mixing13[1] = par2;
  }
  else if(physpar==conv::kPhypMixing23){
    priormode_mixing23 = mode;
    mixmode_mixing23   = mixmode;
    params_mixing23[0] = par1;
    params_mixing23[1] = par2;
  }
  else if(physpar==conv::kPhypMixing12){
    priormode_mixing12 = mode;
    mixmode_mixing12   = mixmode;
    params_mixing12[0] = par1;
    params_mixing12[1] = par2;  
  }
  else if(physpar==conv::kPhypCPV){
    priormode_dcp = mode;
    params_dcp[0] = par1;
    params_dcp[1] = par2;
  }
  else if(physpar==conv::kPhypDm23sq){
    priormode_dm23sq = mode;
    params_dm23sq[0] = par1;
    params_dm23sq[1] = par2;
  }
  else if(physpar==conv::kPhypDm12sq){
    priormode_dm12sq = mode;
    params_dm12sq[0] = par1;
    params_dm12sq[1] = par2;
  }
  else{ //kPhypUndefined
    cerr << endl;
    cerr << "The physics " << util::AsString(physpar) << " parameter is not defined" << endl;
    cerr << endl;
    exit(1);    
  }

  cout << endl;
  cout << "Prior for " << util::AsString(physpar) << " : " << endl;
  if(mode == "gaus"){
    cout << " - gaus : mean = " << par1 << " - sigma = " << par2 << endl;
  }
  else if(mode == "flat"){
    cout << " - flat : low = " << par1 << " - up = " << par2 << endl;
  }
  else if(mode == "user"){
    cout << " - user " << endl;
  }
  else{
    cerr << endl;
    cerr << "Prior mode " << mode << "is not defined for " << util::AsString(physpar) << endl;
    cerr << endl;
    exit(1);
  }

  if(IsMixingAngle(physpar)){
    if( Opt::Instance()->gOscMixingParamChoice == kNuMixSinsqTheta){       
      if(mixmode == "sinsq2th"){
	cout << " - Convert extracted value from sinsq2th to sinsqth" << endl;
      }
    }
    else if( Opt::Instance()->gOscMixingParamChoice == kNuMixSinsq2Theta){       
      if(mixmode == "sinsqth"){
	cout << " - Convert extracted value from sinsqth to sinsq2th" << endl;
      }
    }
    if(mixmode != "sinsq2th" && mixmode != "sinsqth"){
      cerr << endl;
      cerr << "The mixing mode " << mixmode << " is not used" << endl;
      cerr << endl;
      exit(1);
    }
  }

}
//____________________________________________________________________________
void valor::util::OscParams::SetValues(double input_mixing12_nh, double input_mixing13_nh, double input_mixing23_nh,
				       double input_dcp_nh,      double input_dm12sq_nh,   double input_dm23sq_nh,
				       double input_beta_bar,
				       double input_mixing12_ih, double input_mixing13_ih, double input_mixing23_ih,
				       double input_dcp_ih,      double input_dm12sq_ih,   double input_dm23sq_ih,
				       double input_mh_nh_prior, bool marginalize_MH
				       )
{  
  // Check whether some osc. params need to be randomized (i.e. if == -999)
  //NH
  randomize_mixing12_nh = (TMath::Abs(input_mixing12_nh + 999) < kASmallNum);
  randomize_mixing13_nh = (TMath::Abs(input_mixing13_nh + 999) < kASmallNum);
  randomize_mixing23_nh = (TMath::Abs(input_mixing23_nh + 999) < kASmallNum);
  randomize_dcp_nh      = (TMath::Abs(input_dcp_nh      + 999) < kASmallNum);
  randomize_dm12sq_nh   = (TMath::Abs(input_dm12sq_nh   + 999) < kASmallNum);
  randomize_dm23sq_nh   = (TMath::Abs(input_dm23sq_nh   + 999) < kASmallNum);
  //IH
  randomize_mixing12_ih = (TMath::Abs(input_mixing12_ih + 999) < kASmallNum);
  randomize_mixing13_ih = (TMath::Abs(input_mixing13_ih + 999) < kASmallNum);
  randomize_mixing23_ih = (TMath::Abs(input_mixing23_ih + 999) < kASmallNum);
  randomize_dcp_ih      = (TMath::Abs(input_dcp_ih      + 999) < kASmallNum);
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

  //MH
  if(input_mh_nh_prior < -1 && input_mh_nh_prior > -900) input_mh_nh_prior = 0.5;
  cerr << "Default value for input_mh_nh_prior is hard coded. This should be moved into the relevant Conventions file" << endl;
  
  //Set the input values in the OscParams
  //NH
  mixing12_nh = input_mixing12_nh;
  mixing13_nh = input_mixing13_nh;
  mixing23_nh = input_mixing23_nh;
  dcp_nh      = input_dcp_nh;
  dm12sq_nh   = input_dm12sq_nh;
  dm23sq_nh   = input_dm23sq_nh;
  //IH
  mixing12_ih = input_mixing12_ih;
  mixing13_ih = input_mixing13_ih;
  mixing23_ih = input_mixing23_ih;
  dcp_ih      = input_dcp_ih;
  dm12sq_ih   = input_dm12sq_ih;
  dm23sq_ih   = input_dm23sq_ih;
  //Other
  beta_bar       = input_beta_bar;
  marginalize_mh = marginalize_MH;
  mh_nh_prior    = input_mh_nh_prior;
  MH             = Opt::Instance()->gMassHierarchyAssumed;
  iMH            = (MH == kNmhNormal) ? 0 : 1; 
}
//____________________________________________________________________________
void valor::util::OscParams::Randomize()
{
    
  // Mass Hierarchy
  
  if(marginalize_mh){
    
    const double min = 0.;
    const double max = 1.;
    double Pnh = mh_nh_prior;
    double Pih = 1 - Pnh;
    
    if(Pnh < 0. || Pnh > 1.){
      cerr << endl;
      cerr << "Randomize" << endl;
      cerr << "The input probability for NH is not possible" << endl;
      cerr << "P(NH) = " << Pnh << endl;
      cerr << endl;
      exit(1);
    }
    if(Pih < 0. || Pih > 1.){
      cerr << endl;
      cerr << "Randomize" << endl;
      cerr << "The input probability for IH is not possible" << endl;
      cerr << "P(IH) = " << Pih << endl;
      cerr << endl;
      exit(1);
    }
    
    double rndm = Rnd::Instance()->RndToyMCosc().Uniform(min,max);
    
    iMH = conv::kBadResult;
    if(rndm <= Pnh) iMH = 0; // Normal Hierarchy
    else            iMH = 1; // Inverted Hierarchy
    
    if( iMH==0 ){
      Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
      // hier    = kNmhNormal;
      // hierBar = kNmhNormal;
      // skexp->GetSuperkMC()->SetMassHierarchy(hier);
      // skexp->GetSuperkMC()->SetMassHierarchyBar(hierBar);
      // fparm->Copy(*fparm_nh); // NEW MARG MH
    }
    else if( iMH==1 ){
      Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
      // hier    = kNmhInverted;
      // hierBar = kNmhInverted;
      // skexp->GetSuperkMC()->SetMassHierarchy(hier);
      // skexp->GetSuperkMC()->SetMassHierarchyBar(hierBar);
      // fparm->Copy(*fparm_ih); // NEW MARG MH
    }
    
  } // Mass Hierarchy
  
  double par1 = conv::kBadResult;
  double par2 = conv::kBadResult;
  double valrndm = conv::kBadResult;
  string mixmode = "";
  string priormode = "";
  // Normal Hierarchy
  
  if(iMH==0){
    
    // mixing23
    par1 = params_mixing23[0]; 
    par2 = params_mixing23[1];
    priormode = priormode_mixing23;
    mixmode = mixmode_mixing23;
    mixing23_nh = this->ExtractParam(conv::kPhypMixing23,priormode,par1,par2,mixmode);
        
    // mixing13 
    par1 = params_mixing13[0]; 
    par2 = params_mixing13[1];
    priormode = priormode_mixing13;
    mixmode = mixmode_mixing13;
    mixing13_nh = this->ExtractParam(conv::kPhypMixing13,priormode,par1,par2,mixmode);
    
    // deltaCP
    par1 = -TMath::Pi(); 
    par2 = +TMath::Pi();
    priormode = priormode_dcp;
    mixmode = "";
    dcp_nh = this->ExtractParam(conv::kPhypCPV,priormode,par1,par2,mixmode);
    
    // dmsq23
    par1 = 2E-3; 
    par2 = 3E-3;
    priormode = priormode_dm23sq;
    mixmode = "";
    dm23sq_nh = this->ExtractParam(conv::kPhypDm23sq,priormode,par1,par2,mixmode);

    // dmsq12
    par1 = 7.53E-5; 
    par2 = 0.18E-5;
    priormode = priormode_dm12sq;
    mixmode = "";
    dm12sq_nh = this->ExtractParam(conv::kPhypDm12sq,priormode,par1,par2,mixmode);

    // mixing12
    par1 = params_mixing12[0]; 
    par2 = params_mixing12[1];
    priormode = priormode_mixing13;
    mixmode = mixmode_mixing12;
    mixing12_nh = this->ExtractParam(conv::kPhypMixing12,priormode,par1,par2,mixmode);
  }
  
  // Inverted Hierarchy
  
  else if(iMH==1){

    // mixing23
    par1 = params_mixing23[0]; 
    par2 = params_mixing23[1];
    priormode = priormode_mixing23;
    mixmode = mixmode_mixing23;
    mixing23_ih = this->ExtractParam(conv::kPhypMixing23,priormode,par1,par2,mixmode);
        
    // mixing13 
    par1 = params_mixing13[0]; 
    par2 = params_mixing13[1];
    priormode = priormode_mixing13;
    mixmode = mixmode_mixing13;
    mixing13_ih = this->ExtractParam(conv::kPhypMixing13,priormode,par1,par2,mixmode);
    
    // deltaCP
    par1 = -TMath::Pi(); 
    par2 = +TMath::Pi();
    priormode = priormode_dcp;
    mixmode = "";
    dcp_ih = this->ExtractParam(conv::kPhypCPV,priormode,par1,par2,mixmode);
    
    // dmsq23
    par1 = 2E-3; 
    par2 = 3E-3;
    priormode = priormode_dm23sq;
    mixmode = "";
    dm23sq_ih = this->ExtractParam(conv::kPhypDm23sq,priormode,par1,par2,mixmode);

    // dmsq12
    par1 = 7.53E-5; 
    par2 = 0.18E-5;
    priormode = priormode_dm12sq;
    mixmode = "";
    dm12sq_ih = this->ExtractParam(conv::kPhypDm12sq,priormode,par1,par2,mixmode);

    // mixing12
    par1 = params_mixing12[0]; 
    par2 = params_mixing12[1];
    priormode = priormode_mixing13;
    mixmode = mixmode_mixing12;
    mixing12_ih = this->ExtractParam(conv::kPhypMixing12,priormode,par1,par2,mixmode);
  }
  
  else{
    cerr << endl;
    cerr << "Mass Hierarchy is not defined" << endl;
    cerr << endl;
    exit(1);
  }
      
}
//____________________________________________________________________________
double valor::util::OscParams::ExtractParam(conv::PhysParamId_t physpar,string mode,double par1,double par2,string mixmode){
  
  double valrndm = conv::kBadResult;

  if(mode=="gaus")
    valrndm = Rnd::Instance()->RndToyMCosc().Gaus(par1,par2);
  else if(mode=="flat")
    valrndm = Rnd::Instance()->RndToyMCosc().Uniform(par1,par2);
  
  if(IsMixingAngle(physpar)){    
    if(Opt::Instance()->gOscMixingParamChoice == kNuMixSinsqTheta){
      if(mixmode == "sinsqth")       return valrndm;
      else if(mixmode == "sinsq2th") return util::FromSinsq2thToSinsqth(valrndm);
    }
    else if(Opt::Instance()->gOscMixingParamChoice == kNuMixSinsq2Theta){
      if(mixmode == "sinsqth")       return util::FromSinsqthToSinsq2th(valrndm);
      else if(mixmode == "sinsq2th") return valrndm;
    }
  }
  return valrndm;
}
//____________________________________________________________________________
void valor::util::OscParams::UpdatePhysHypothesis(PhysHypothesis& ph, bool print)
{
  double input_mixing12_used, input_mixing13_used, input_mixing23_used;
  double input_dcp_used,      input_dm12sq_used,   input_dm23sq_used;
  double input_beta_bar_used;

  //Get the values to use (based on MH)
  if(iMH == 0) {
    //NH
    input_mixing12_used = mixing12_nh;
    input_mixing13_used = mixing13_nh;
    input_mixing23_used = mixing23_nh;
    input_dcp_used      = dcp_nh;
    input_dm12sq_used   = dm12sq_nh;
    input_dm23sq_used   = dm23sq_nh;
  }
  else if(iMH == 1) {
    //IH
    input_mixing12_used = mixing12_ih;
    input_mixing13_used = mixing13_ih;
    input_mixing23_used = mixing23_ih;
    input_dcp_used      = dcp_ih;
    input_dm12sq_used   = dm12sq_ih;
    input_dm23sq_used   = dm23sq_ih;
  }
  else {
    cerr << "SetOscParams has unknown value of iMH: " << iMH << endl;
    exit(1);
  }
  //Other 
  input_beta_bar_used = beta_bar;  // It's never randomized
  
  //Update the PhysHypothesis
  if(ph.Hypothesis() == kHypNuMuDisapp2F) {
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
  }//kHypNuMuDisapp2F
  else if(ph.Hypothesis() == kHypOsc3F) {  
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
  }//kHypOsc3F
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
  }//kHypOsc3FBar
  else if(ph.Hypothesis() == kHypOsc3FNormBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypBetaBar,    input_beta_bar_used);
  }//kHypOsc3FNormBar
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
  }//kHypOsc3FBarNormBar
  else {
    cout << "SetOscParams() " << util::AsString(ph.Hypothesis()) << " not implemented" << endl;
    exit(1);
  }
  
  if(print)
    cout << "- Physics hypothesis set: " << ph << endl;

}
//____________________________________________________________________________
void valor::util::OscParams::UpdateFitParamList(FitParamList& fparm, bool set_init, bool print)
{
  double input_mixing12_used, input_mixing13_used, input_mixing23_used;
  double input_dcp_used,      input_dm12sq_used,   input_dm23sq_used;
  double input_beta_bar_used;

  //Get the values to use (based on MH)
  if(iMH == 0) {
    //NH
    input_mixing12_used = mixing12_nh;
    input_mixing13_used = mixing13_nh;
    input_mixing23_used = mixing23_nh;
    input_dcp_used      = dcp_nh;
    input_dm12sq_used   = dm12sq_nh;
    input_dm23sq_used   = dm23sq_nh;
  }
  else if(iMH == 1) {
    //IH
    input_mixing12_used = mixing12_ih;
    input_mixing13_used = mixing13_ih;
    input_mixing23_used = mixing23_ih;
    input_dcp_used      = dcp_ih;
    input_dm12sq_used   = dm12sq_ih;
    input_dm23sq_used   = dm23sq_ih;
  }
  else {
    cerr << "SetOscParams has unknown value of iMH: " << iMH << endl;
    exit(1);
  }
  //Other
  input_beta_bar_used = beta_bar;
  
  //Update the FitParamList
  if(fparm.PhysicsHypothesis() == kHypNuMuDisapp2F) {
    if(set_init) {
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
    }
  }//kHypNuMuDisapp2F
  else if(fparm.PhysicsHypothesis() == kHypOsc3F) {  
    if(set_init) {
      fparm.SetInitValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
    }
  }//kHypOsc3F
  else if(fparm.PhysicsHypothesis() == kHypOsc3FBar) {
    if(set_init) {
      fparm.SetInitValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetInitValue(kPhypMixing12Bar,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13Bar,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23Bar,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPVBar    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sqBar,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sqBar,    input_dm23sq_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetCurrentValue(kPhypMixing12Bar,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13Bar,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23Bar,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPVBar    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sqBar,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sqBar,    input_dm23sq_used);
    }
  }//kHypOsc3FBar
  else if(fparm.PhysicsHypothesis() == kHypOsc3FNormBar) {
    if(set_init) {
      fparm.SetInitValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetInitValue(kPhypBetaBar,    input_beta_bar_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetCurrentValue(kPhypBetaBar,    input_beta_bar_used);
    }
  }//kHypOsc3FNormBar
  else if(fparm.PhysicsHypothesis() == kHypOsc3FBarNormBar) {
    if(set_init) {
      fparm.SetInitValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetInitValue(kPhypMixing12Bar,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13Bar,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23Bar,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPVBar    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sqBar,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sqBar,    input_dm23sq_used);
       
      fparm.SetInitValue(kPhypBetaBar,    input_beta_bar_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetCurrentValue(kPhypMixing12Bar,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13Bar,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23Bar,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPVBar    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sqBar,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sqBar,    input_dm23sq_used);
       
      fparm.SetCurrentValue(kPhypBetaBar,    input_beta_bar_used);
    }
  }//kHypOsc3FBarNormBar
  else {
    cout << "SetOscParams() " << util::AsString(fparm.PhysicsHypothesis()) << " not implemented" << endl;
    exit(1);
  }
  
  if(print)
    cout << "- Oscillation parameters in physics param list updated: " << fparm << endl;

  //save the values used
  if(iMH == 0) {
    if(hToy_th13_dcp_nh)
      hToy_th13_dcp_nh->Fill(input_mixing13_used, input_dcp_used);
    if(hToy_23_nh)
      hToy_23_nh->Fill(input_mixing23_used, input_dm23sq_used);
  }//NH
  else if(iMH == 1) {
    if(hToy_th13_dcp_ih)
      hToy_th13_dcp_ih->Fill(input_mixing13_used, input_dcp_used);
    if(hToy_23_ih)
      hToy_23_ih->Fill(input_mixing23_used, input_dm23sq_used);
  }//IH

}
//____________________________________________________________________________
void valor::util::OscParams::Print(ostream & stream) const
{
  stream << "   - True oscillation parameters: " << endl
	 << "     - Normal ordering: " << endl
	 << "       > " << util::AsString(kPhypMixing12) << " = "; ((randomize_mixing12_nh) ? cout << "* randomized *" : cout << mixing12_nh) << endl
	 << "       > " << util::AsString(kPhypMixing13) << " = "; ((randomize_mixing13_nh) ? cout << "* randomized *" : cout << mixing13_nh) << endl
	 << "       > " << util::AsString(kPhypMixing23) << " = "; ((randomize_mixing23_nh) ? cout << "* randomized *" : cout << mixing23_nh) << endl
	 << "       > " << util::AsString(kPhypCPV    )  << " = "; ((randomize_dcp_nh)      ? cout << "* randomized *" : cout << dcp_nh) << endl
	 << "       > " << util::AsString(kPhypDm12sq)   << " = "; ((randomize_dm12sq_nh)   ? cout << "* randomized *" : cout << dm12sq_nh )  << " eV^2/c^4" << endl
	 << "       > " << util::AsString(kPhypDm23sq)   << " = "; ((randomize_dm23sq_nh)   ? cout << "* randomized *" : cout << dm23sq_nh )  << " eV^2/c^4" << endl
	 << "     - Inverted ordering: " << endl
	 << "       > " << util::AsString(kPhypMixing12) << " = "; ((randomize_mixing12_ih) ? cout << "* randomized *" : cout << mixing12_ih) << endl
	 << "       > " << util::AsString(kPhypMixing13) << " = "; ((randomize_mixing13_ih) ? cout << "* randomized *" : cout << mixing13_ih) << endl
	 << "       > " << util::AsString(kPhypMixing23) << " = "; ((randomize_mixing23_ih) ? cout << "* randomized *" : cout << mixing23_ih) << endl
	 << "       > " << util::AsString(kPhypCPV    )  << " = "; ((randomize_dcp_ih)      ? cout << "* randomized *" : cout << dcp_ih) << endl
	 << "       > " << util::AsString(kPhypDm12sq)   << " = "; ((randomize_dm12sq_ih)   ? cout << "* randomized *" : cout << dm12sq_ih )  << " eV^2/c^4" << endl
	 << "       > " << util::AsString(kPhypDm23sq)   << " = "; ((randomize_dm23sq_ih)   ? cout << "* randomized *" : cout << dm23sq_ih )  << " eV^2/c^4" << endl
	 << "     - Other: " << endl
	 << "       > " << util::AsString(kPhypBetaBar)  << " = "; ((randomize_beta_bar)    ? cout << "* randomized *" : cout << beta_bar) << endl
	 << "       > Mass hierarchy = "; ((marginalize_mh) ? cout << "Margninalized (NH prior = " << mh_nh_prior << ")" : cout << util::AsString(MH)) << endl;
  if(!priors_file.empty()) {
    if((randomize_mixing13_nh && randomize_dcp_nh)
       || (randomize_mixing13_ih && randomize_dcp_ih)
       || (randomize_mixing23_nh && randomize_dm23sq_nh)
       || (randomize_mixing23_ih && randomize_dm23sq_ih)
       ) {
      stream << "       > Using priors from file: " << priors_file << " for parameters:" << endl;
      if(randomize_mixing13_nh && randomize_dcp_nh)
	stream << "          - " << util::AsString(kPhypMixing13) << " vs " << util::AsString(kPhypCPV) << " (NH)" << endl;
      if(randomize_mixing13_ih && randomize_dcp_ih)
	stream << "          - " << util::AsString(kPhypMixing13) << " vs " << util::AsString(kPhypCPV) << " (IH)" << endl;
      if(randomize_mixing23_nh && randomize_dm23sq_nh)
	stream << "          - " << util::AsString(kPhypMixing23) << " vs " << util::AsString(kPhypDm23sq) << " (NH)" << endl;
      if(randomize_mixing23_ih && randomize_dm23sq_ih)
	stream << "          - " << util::AsString(kPhypMixing23) << " vs " << util::AsString(kPhypDm23sq) << " (IH)" << endl;
    }
  }
  stream << endl;
  
}
//____________________________________________________________________________
void valor::util::OscParams::Write(TFile * f)
{
  f->cd();

  //priors
  if(hPrior_th13_dcp_nh) {
    hPrior_th13_dcp_nh->SetDirectory(0);
    hPrior_th13_dcp_nh->Write();
  }
  if(hPrior_th13_dcp_ih) {
    hPrior_th13_dcp_ih->SetDirectory(0);
    hPrior_th13_dcp_ih->Write();
  }
  if(hPrior_23_nh) {
    hPrior_23_nh->SetDirectory(0);
    hPrior_23_nh->Write();
  }
  if(hPrior_23_ih) {
    hPrior_23_ih->SetDirectory(0);
    hPrior_23_ih->Write();
  }
  //toys used
  if(hToy_th13_dcp_nh) {
    hToy_th13_dcp_nh->SetDirectory(0);
    hToy_th13_dcp_nh->Write();
  }
  if(hToy_th13_dcp_ih) {
    hToy_th13_dcp_ih->SetDirectory(0);
    hToy_th13_dcp_ih->Write();
  }
  if(hToy_23_nh) {
    hToy_23_nh->SetDirectory(0);
    hToy_23_nh->Write();
  }
  if(hToy_23_ih) {
    hToy_23_ih->SetDirectory(0);
    hToy_23_ih->Write();
  }
  if(hToy_MH) {
    hToy_MH->SetDirectory(0);
    hToy_MH->Write();
  }
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
