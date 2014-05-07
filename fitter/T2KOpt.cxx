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
