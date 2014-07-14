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
  std::ostream & operator<<(std::ostream & stream, const T2KOpt & opt)
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
  bool fit_MultiRmu = false;
  bool fit_MultiRe  = false;

  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsFCFV1RingMuLike) fit_1Rmu = true;
    if(s == kSKsFCFV1RingELike ) fit_1Re  = true;
    if(s == kSKsRHCFCFV1RingMuLike) fit_RHC_1Rmu = true;
    if(s == kSKsRHCFCFV1RingELike ) fit_RHC_1Re  = true;
    if(s == kSKsFCFVMultiRingMuLike) fit_MultiRmu = true;
    if(s == kSKsFCFVMultiRingELike ) fit_MultiRe  = true;
  }

  // MultiRing FHC
#if defined __MULTIRING_FHC_RUNNING__
#ifdef __NUMUONLY__
  fit_MultiRe  = false;
#endif
#ifdef __NUEONLY__
  fit_MultiRmu = false;
#endif
  // No MultiRing FHC
#else
  fit_MultiRmu = false;
  fit_MultiRe = false;
#endif

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


  // w/ MultiRing sample
  if(fit_MultiRmu || fit_MultiRe){
    
    // MultiRing-FHC + 1Ring-FHC&RHC
    if(fit_1Rmu && fit_1Re && fit_RHC_1Rmu && fit_RHC_1Re & fit_MultiRmu && fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingFHCAnd1RingFHCAndRHCCombined;
    else if(fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && !fit_RHC_1Re & fit_MultiRmu && !fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingFHCAnd1RingFHCAndRHCNuMuDisappOnly;
    else if(!fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && fit_RHC_1Re & !fit_MultiRmu && fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingFHCAnd1RingFHCAndRHCNuEAppOnly;
    
    // MultiRing-FHC + 1Ring-FHC
    else if(fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re & fit_MultiRmu && fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingAnd1RingCombined;
    else if(fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re & fit_MultiRmu && !fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingAnd1RingNuMuDisappOnly;
    else if(!fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re & !fit_MultiRmu && fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingAnd1RingNuEAppOnly;
    
    // MultiRing-FHC 
    else if(!fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re & fit_MultiRmu && fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingCombined;
    else if(!fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re & fit_MultiRmu && !fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingNuMuDisappOnly;
    else if(!fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re & !fit_MultiRmu && fit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingNuEAppOnly;
  }
  // w/o MultiRing sample
  else {
    
    // 1Ring-FHC&RHC 
    if(fit_1Rmu && fit_1Re && fit_RHC_1Rmu && fit_RHC_1Re)
      gOscAnaTypeChoice = kAnaFHCAndRHCCombined;
    else if(fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && !fit_RHC_1Re) 
      gOscAnaTypeChoice = kAnaFHCAndRHCNuMuDisappOnly;
    else if(!fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && fit_RHC_1Re)
      gOscAnaTypeChoice = kAnaFHCAndRHCNuEAppOnly;
        
    // 1Ring-FHC 
    else if(fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
      gOscAnaTypeChoice = kAnaCombined;
    else if(fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
      gOscAnaTypeChoice = kAnaNuMuDisappOnly;
    else if(!fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
      gOscAnaTypeChoice = kAnaNuEAppOnly;
    
    // 1Ring-RHC 
    else if(!fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && fit_RHC_1Re)
      gOscAnaTypeChoice = kAnaRHCCombined;
    else if(!fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && !fit_RHC_1Re)
      gOscAnaTypeChoice = kAnaRHCNuMuDisappOnly;
    else if(!fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && fit_RHC_1Re)
      gOscAnaTypeChoice = kAnaRHCNuEAppOnly;
    
    // Nothing
    else
      gOscAnaTypeChoice = kAnaUndefined; 
  }
  

  // OLD

// #if defined __FHC_RUNNING__ && defined __RHC_RUNNING__
// #ifdef __NUMUONLY__
//   fit_1Re  = false;
//   fit_RHC_1Re  = false;
// #endif
// #ifdef __NUEONLY__
//   fit_1Rmu = false;
//   fit_RHC_1Rmu = false;
// #endif
// #else
// #ifdef __FHC_RUNNING__
// #ifdef __NUMUONLY__
//   fit_1Re  = false;
//   fit_RHC_1Re  = false;
//   fit_RHC_1Rmu = false;
// #endif
// #ifdef __NUEONLY__
//   fit_1Rmu = false;
//   fit_RHC_1Rmu = false;
//   fit_RHC_1Re  = false;
// #endif
// #else
// #ifdef __RHC_RUNNING__
// #ifdef __NUMUONLY__
//   fit_1Re  = false;
//   fit_RHC_1Re  = false;
//   fit_1Rmu = false;
// #endif
// #ifdef __NUEONLY__
//   fit_1Rmu = false;
//   fit_RHC_1Rmu = false;
//   fit_1Re  = false;
// #endif
// #endif
// #endif
// #endif

  // if(fit_1Rmu && fit_1Re && fit_RHC_1Rmu && fit_RHC_1Re)
  //   gOscAnaTypeChoice = kAnaFHCAndRHCCombined;
  // else if(fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && !fit_RHC_1Re) 
  //   gOscAnaTypeChoice = kAnaFHCAndRHCNuMuDisappOnly;
  // else if(!fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && fit_RHC_1Re)
  //   gOscAnaTypeChoice = kAnaFHCAndRHCNuEAppOnly;
  // else if(fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
  //   gOscAnaTypeChoice = kAnaCombined;
  // else if(fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
  //   gOscAnaTypeChoice = kAnaNuMuDisappOnly;
  // else if(!fit_1Rmu && fit_1Re && !fit_RHC_1Rmu && !fit_RHC_1Re)
  //   gOscAnaTypeChoice = kAnaNuEAppOnly;
  // else if(!fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && fit_RHC_1Re)
  //   gOscAnaTypeChoice = kAnaRHCCombined;
  // else if(!fit_1Rmu && !fit_1Re && fit_RHC_1Rmu && !fit_RHC_1Re)
  //   gOscAnaTypeChoice = kAnaRHCNuMuDisappOnly;
  // else if(!fit_1Rmu && !fit_1Re && !fit_RHC_1Rmu && fit_RHC_1Re)
  //   gOscAnaTypeChoice = kAnaRHCNuEAppOnly;
  // else
  //   gOscAnaTypeChoice = kAnaUndefined;
    
}
//____________________________________________________________________________
bool T2KOpt::Undefined(void)
{

  return false;
}
//____________________________________________________________________________
void T2KOpt::Print(std::ostream & stream) const
{
  stream << " - Oscillation analysis choice: " 
         << util::AsString(gOscAnaTypeChoice) 
         << endl;
}
//________________________________________________________________________________________
