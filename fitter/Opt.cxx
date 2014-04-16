//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include <string>
#include <cstdlib>

#include "Opt.h"
#include "Utilities.h"

using std::endl;
using std::string;

using namespace t2kana;
using namespace t2kana::conv;

ClassImp(Opt::Persistable);

//________________________________________________________________________________________
namespace t2kana {
  ostream & operator<<(ostream & stream, const Opt & opt)
  {
    opt.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
Opt * Opt::fInstance = 0;
//____________________________________________________________________________
Opt::Opt()
{
  this->SetDefaults();
  fInstance = 0;
}
//________________________________________________________________________________________
Opt::Opt(const Opt &) 
{

}
//________________________________________________________________________________________
Opt::~Opt()
{

}
//____________________________________________________________________________
Opt * Opt::Instance()
{
  if(fInstance == 0) {
    static Opt::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new Opt;
  }
  return fInstance;
}
//____________________________________________________________________________
void Opt::SetDefaults(void) 
{
  gMassHierarchyAssumed = kNmhNormal;
  gOscMixingParamChoice = kNuMixSinsq2Theta;
  gOscDmsqScheme        = kNuDmsq12andAbs23;
  gOscDeltaCP           = kNuDeltaCPPhase;
  gOscOctantChoice      = kNuOctantOne;

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
bool Opt::Undefined(void)
{
  if(gMassHierarchyAssumed == kNmhUndefined)
    return true;
  if(gOscMixingParamChoice == kNuMixUndefined)
    return true;
  if(gOscDmsqScheme        == kNuDmsqUndefined)
    return true;
  if(gOscDeltaCP           == kNuDeltaCPUndefined)
    return true;
  if(gOscOctantChoice      == kNuOctantUndefined)
    return true;

  return false;
}
//____________________________________________________________________________
void Opt::Print(ostream & stream) const
{
  stream << "[Global analysis options]" 
         << endl;
  stream << " - Hierarchy: " 
         << util::AsString(gMassHierarchyAssumed) 
         << endl;
  stream << " - Mixing param choice: " 
         << util::AsString(gOscMixingParamChoice) 
         << endl;
  stream << " - 23 octant choice: "
	 << util::AsString(gOscOctantChoice)
	 << endl;
  stream << " - Independent squared mass splittings: "
         << util::AsString(gOscDmsqScheme)
         << endl;
  stream << " - Choice for Delta CP phase: "
         << util::AsString(gOscDeltaCP)
         << endl;
  stream << " - Oscillation analysis choice: " 
         << util::AsString(gOscAnaTypeChoice) 
         << endl;
}
//________________________________________________________________________________________
void Opt::WriteTo (Opt::Persistable & optp)
{
  optp.mh    = gMassHierarchyAssumed; 
  optp.mixp  = gOscMixingParamChoice;
  optp.oct   = gOscOctantChoice; 
  optp.dmsqs = gOscDmsqScheme;
  optp.dcp   = gOscDeltaCP;
}
//________________________________________________________________________________________
void Opt::ReadFrom (Opt::Persistable & optp)
{
  gMassHierarchyAssumed = optp.mh;
  gOscMixingParamChoice = optp.mixp;
  gOscOctantChoice      = optp.oct;
  gOscDmsqScheme        = optp.dmsqs;
  gOscDeltaCP           = optp.dcp;
}
//________________________________________________________________________________________

