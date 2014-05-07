//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include <string>
#include <cstdlib>

#include "T2KOpt.h"
#include "Utilities/Rnd.h"
#include "Utilities.h"

using std::endl;
using std::string;

using namespace valor;
using namespace valor::conv;

//________________________________________________________________________________________
namespace valor {
  ostream & operator<<(ostream & stream, const T2KOpt & opt)
  {
    opt.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
T2KOpt * T2KOpt::fInstance = 0;
//____________________________________________________________________________
T2KOpt::T2KOpt()
{
  this->SetDefaults();
  fInstance = 0;
}
//________________________________________________________________________________________
T2KOpt::T2KOpt(const T2KOpt &) 
{

}
//________________________________________________________________________________________
T2KOpt::~T2KOpt()
{

}
//____________________________________________________________________________
T2KOpt * T2KOpt::Instance()
{
  if(fInstance == 0) {
    static T2KOpt::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new T2KOpt;
  }
  return fInstance;
}
//____________________________________________________________________________
void T2KOpt::SetDefaults(void) 
{

  bool fit_1Rmu = false;
  bool fit_1Re  = false;
  bool fit_RHC_1Rmu = false;
  bool fit_RHC_1Re  = false;

  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsFCFV1RingMuLike) fit_1Rmu = true;
    if(s == kSKsFCFV1RingELike ) fit_1Re  = true;
    if(s == kSKsRHCFCFV1RingMuLike) fit_RHC_1Rmu = true;
    if(s == kSKsRHCFCFV1RingELike ) fit_RHC_1Re  = true;
  }

#if defined __FHC_RUNNING__ && defined __RHC_RUNNING__
#ifdef __NUMUONLY__
  fit_1Re  = false;
  fit_RHC_1Re  = false;
#endif
#ifdef __NUEONLY__
  fit_1Rmu = false;
  fit_RHC_1Rmu = false;
#endif
#else
#ifdef __FHC_RUNNING__
#ifdef __NUMUONLY__
  fit_1Re  = false;
  fit_RHC_1Re  = false;
  fit_RHC_1Rmu = false;
#endif
#ifdef __NUEONLY__
  fit_1Rmu = false;
  fit_RHC_1Rmu = false;
  fit_RHC_1Re  = false;
#endif
#else
#ifdef __RHC_RUNNING__
#ifdef __NUMUONLY__
  fit_1Re  = false;
  fit_RHC_1Re  = false;
  fit_1Rmu = false;
#endif
#ifdef __NUEONLY__
  fit_1Rmu = false;
  fit_RHC_1Rmu = false;
  fit_1Re  = false;
#endif
#endif
#endif
#endif

  if(fit_1Rmu && fit_1Re && fit_RHC_1Rmu && fit_RHC_1Re)
    gOscAnaTypeChoice = kAnaFHCAndRHCCombined;
  else if(fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && !fit_RHC_1Re) 
    gOscAnaTypeChoice = kAnaFHCAndRHCNuMuDisappOnly;
  else if(!fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && fit_RHC_1Re)
    gOscAnaTypeChoice = kAnaFHCAndRHCNuEAppOnly;
  else if(fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
    gOscAnaTypeChoice = kAnaCombined;
  else if(fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
    gOscAnaTypeChoice = kAnaNuMuDisappOnly;
  else if(!fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
    gOscAnaTypeChoice = kAnaNuEAppOnly;
  else if(!fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && fit_RHC_1Re)
    gOscAnaTypeChoice = kAnaRHCCombined;
  else if(!fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && !fit_RHC_1Re)
    gOscAnaTypeChoice = kAnaRHCNuMuDisappOnly;
  else if(!fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && fit_RHC_1Re)
    gOscAnaTypeChoice = kAnaRHCNuEAppOnly;
  else
    gOscAnaTypeChoice = kAnaUndefined;
    
  this->SetDefaultParameterValues();
}
//____________________________________________________________________________
void T2KOpt::SetDefaultParameterValues()
{
  // Set default global parameter values
  // set empty values
  for (int ip = 0 ; ip < knppmax ; ip++ ) {
    fGlobalBestFit [ip] = kBadResult;
    fGlobalFitError[ip] = kBadResult;
    //fInitialValue  [ip] = 0.;
  }
  
  NuOscDmsqScheme_t odmsq = Opt::Instance()->gOscDmsqScheme;
  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  
  double mean;
  
  //dm12
  fGlobalBestFit  [kPhypDm12sq] = kGlobFitDm12sq;
  fGlobalFitError [kPhypDm12sq] = kGlobFitErrDm12sq;
  //fInitialValue   [kPhypDm12sq] = kInitFitDm12sq;
 
  //dm23
  switch (odmsq) {
    case(kNuDmsqFogliLisi):
      mean  =  kGlobFitDmsqFL;
      break;
    case(kNuDmsq12andAbs23):
      mean  = kGlobFitDm23sq; 
      break;
    case(kNuDmsq12andAbs31):
      if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal)
        mean  = kGlobFitDm23sq + kGlobFitDm12sq; 
      else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted)
        mean  = kGlobFitDm23sq - kGlobFitDm12sq; 
      break;
    default:
      mean  = kGlobFitDm23sq;  
  }
  fGlobalBestFit  [kPhypDm23sq] = mean;
  fGlobalFitError [kPhypDm23sq] = kGlobFitErrDm23sq;
  //fInitialValue   [kPhypDm23sq] = fGlobalBestFit[kPhypDm23sq];
  
  // mixing angles
  switch(omp) {
    case(kNuMixSinsq2Theta):
    {
      fGlobalBestFit  [kPhypMixing12] = kGlobFitSinsq2Theta12;
      fGlobalFitError [kPhypMixing12] = kGlobFitErrSinsq2Theta12;
      //fInitialValue   [kPhypMixing12] = fGlobalBestFit[kPhypMixing12];
      
      fGlobalBestFit  [kPhypMixing23] = kGlobFitSinsq2Theta23;
      fGlobalFitError [kPhypMixing23] = kGlobFitErrSinsq2Theta23;
      //fInitialValue   [kPhypMixing23] = fGlobalBestFit[kPhypMixing23];
      
      fGlobalBestFit  [kPhypMixing13] = kGlobFitSinsq2Theta13;
      fGlobalFitError [kPhypMixing13] = kGlobFitErrSinsq2Theta13;
      //fInitialValue   [kPhypMixing13] = fGlobalBestFit[kPhypMixing13];
      break;
    }
    case(kNuMixSinsqTheta):
    {
      fGlobalBestFit  [kPhypMixing12] = kGlobFitSinsqTheta12;
      fGlobalFitError [kPhypMixing12] = kGlobFitErrSinsqTheta12;
      //fInitialValue   [kPhypMixing12] = fGlobalBestFit[kPhypMixing12];
      
      fGlobalBestFit  [kPhypMixing23] = kGlobFitSinsqTheta23;
      fGlobalFitError [kPhypMixing23] = kGlobFitErrSinsqTheta23;
      //fInitialValue   [kPhypMixing23] = fGlobalBestFit[kPhypMixing23];
      
      fGlobalBestFit  [kPhypMixing13] = kGlobFitSinsqTheta13;
      fGlobalFitError [kPhypMixing13] = kGlobFitErrSinsqTheta13;
      //fInitialValue   [kPhypMixing13] = fGlobalBestFit[kPhypMixing13];
      break;
    }
    default: exit(-1);
  }
}
//____________________________________________________________________________
bool T2KOpt::Undefined(void)
{

  return false;
}
//____________________________________________________________________________
void T2KOpt::Print(ostream & stream) const
{
  stream << " - Oscillation analysis choice: " 
         << util::AsString(gOscAnaTypeChoice) 
         << endl;
}
//________________________________________________________________________________________
//________________________________________________________________________________________
//________________________________________________________________________________________
double T2KOpt::GetGlobalBestFit(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fGlobalBestFit[p];
}
//________________________________________________________________________________________
double T2KOpt::GetGlobalFitError(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fGlobalFitError[p];
}
//________________________________________________________________________________________
double T2KOpt::GetRandomGaus(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  
  double mean  = fGlobalBestFit[p]; 
  double sigma = fGlobalFitError[p]*mean;
  
  assert(mean!=kBadResult);
  double randval = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
  return randval;
}
//________________________________________________________________________________________
