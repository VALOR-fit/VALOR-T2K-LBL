//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @ 18/2/13 - LE                                                                                                                               
 Added back the selection of dmsq schema,from branch 2010_v3: choice between usual Dm12sq + Dm23sq and Fogli&Lisi schema

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

//____________________________________________________________________________
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
//____________________________________________________________________________
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
#ifdef __FHC_RUNNING__
    if(s == kSKsFCFV1RingMuLike) fit_1Rmu = true;
    if(s == kSKsFCFV1RingELike ) fit_1Re  = true;
#endif
#ifdef __RHC_RUNNING__
    if(s == kSKsRHCFCFV1RingMuLike) fit_RHC_1Rmu = true;
    if(s == kSKsRHCFCFV1RingELike ) fit_RHC_1Re  = true;
#endif
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
  fit_1Re  = false;
#endif
#endif
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
//____________________________________________________________________________

