//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include <string>
#include <cstdlib>

#include "TSystem.h"

#include "T2KOpt.h"
#include "valor/Utilities/Rnd.h"
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
int T2KOpt::ReadXmlPath(std::string config_file)
{
  FileStat_t fs;
  int status = gSystem->GetPathInfo(config_file.c_str(),fs);
  if (status) {
    std::cerr<<"T2KOpt/ReadXmlPath > Warning, cannot access file: "<<config_file<<endl;
    return -1;
  }
  fConfigFileName = config_file;
  return 0;
}
//____________________________________________________________________________
void T2KOpt::SetDefaults(void) 
{
  fConfigFileName = "T2KOpt - No config file set";
  
  fFit_1Rmu     = false;
  fFit_1Re      = false;
  fFit_RHC_1Rmu = false;
  fFit_RHC_1Re  = false;
  fFit_MultiRmu = false;
  fFit_MultiRe  = false;

  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsFCFV1RingMuLike    ) fFit_1Rmu = true;
    if(s == kSKsFCFV1RingELike     ) fFit_1Re  = true;
    if(s == kSKsRHCFCFV1RingMuLike ) fFit_RHC_1Rmu = true;
    if(s == kSKsRHCFCFV1RingELike  ) fFit_RHC_1Re  = true;
    if(s == kSKsFCFVMultiRingMuLike) fFit_MultiRmu = true;
    if(s == kSKsFCFVMultiRingELike ) fFit_MultiRe  = true;
  }

  // MultiRing FHC
#if defined __MULTIRING_FHC_RUNNING__
#ifdef __NUMUONLY__
  fFit_MultiRe  = false;
#endif
#ifdef __NUEONLY__
  fFit_MultiRmu = false;
#endif
  // No MultiRing FHC
#else
  fFit_MultiRmu = false;
  fFit_MultiRe = false;
#endif

#if defined __FHC_RUNNING__ && defined __RHC_RUNNING__
#ifdef __NUMUONLY__
  fFit_1Re  = false;
  fFit_RHC_1Re  = false;
#endif
#ifdef __NUEONLY__
  fFit_1Rmu = false;
  fFit_RHC_1Rmu = false;
#endif
#else
#ifdef __FHC_RUNNING__
#ifdef __NUMUONLY__
  fFit_1Re  = false;
  fFit_RHC_1Re  = false;
  fFit_RHC_1Rmu = false;
#endif
#ifdef __NUEONLY__
  fFit_1Rmu = false;
  fFit_RHC_1Rmu = false;
  fFit_RHC_1Re  = false;
#endif
#else
#ifdef __RHC_RUNNING__
#ifdef __NUMUONLY__
  fFit_1Re  = false;
  fFit_RHC_1Re  = false;
  fFit_1Rmu = false;
#endif
#ifdef __NUEONLY__
  fFit_1Rmu = false;
  fFit_RHC_1Rmu = false;
  fFit_1Re  = false;
#endif
#endif
#endif
#endif


  // w/ MultiRing sample
  if(fFit_MultiRmu || fFit_MultiRe){
    
    // MultiRing-FHC + 1Ring-FHC&RHC
    if(fFit_1Rmu && fFit_1Re && fFit_RHC_1Rmu && fFit_RHC_1Re & fFit_MultiRmu && fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingFHCAnd1RingFHCAndRHCCombined;
    else if(fFit_1Rmu && !fFit_1Re && fFit_RHC_1Rmu && !fFit_RHC_1Re & fFit_MultiRmu && !fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingFHCAnd1RingFHCAndRHCNuMuDisappOnly;
    else if(!fFit_1Rmu && fFit_1Re && !fFit_RHC_1Rmu && fFit_RHC_1Re & !fFit_MultiRmu && fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingFHCAnd1RingFHCAndRHCNuEAppOnly;
    
    // MultiRing-FHC + 1Ring-FHC
    else if(fFit_1Rmu && fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re & fFit_MultiRmu && fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingAnd1RingCombined;
    else if(fFit_1Rmu && !fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re & fFit_MultiRmu && !fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingAnd1RingNuMuDisappOnly;
    else if(!fFit_1Rmu && fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re & !fFit_MultiRmu && fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingAnd1RingNuEAppOnly;
    
    // MultiRing-FHC 
    else if(!fFit_1Rmu && !fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re & fFit_MultiRmu && fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingCombined;
    else if(!fFit_1Rmu && !fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re & fFit_MultiRmu && !fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingNuMuDisappOnly;
    else if(!fFit_1Rmu && !fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re & !fFit_MultiRmu && fFit_MultiRe)
      gOscAnaTypeChoice = kAnaMultiRingNuEAppOnly;
  }
  // w/o MultiRing sample
  else {
    
    // 1Ring-FHC&RHC 
    if(fFit_1Rmu && fFit_1Re && fFit_RHC_1Rmu && fFit_RHC_1Re)
      gOscAnaTypeChoice = kAnaFHCAndRHCCombined;
    else if(fFit_1Rmu && !fFit_1Re && fFit_RHC_1Rmu && !fFit_RHC_1Re) 
      gOscAnaTypeChoice = kAnaFHCAndRHCNuMuDisappOnly;
    else if(!fFit_1Rmu && fFit_1Re && !fFit_RHC_1Rmu && fFit_RHC_1Re)
      gOscAnaTypeChoice = kAnaFHCAndRHCNuEAppOnly;
        
    // 1Ring-FHC 
    else if(fFit_1Rmu && fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re)
      gOscAnaTypeChoice = kAnaCombined;
    else if(fFit_1Rmu && !fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re)
      gOscAnaTypeChoice = kAnaNuMuDisappOnly;
    else if(!fFit_1Rmu && fFit_1Re && !fFit_RHC_1Rmu && !fFit_RHC_1Re)
      gOscAnaTypeChoice = kAnaNuEAppOnly;
    
    // 1Ring-RHC 
    else if(!fFit_1Rmu && !fFit_1Re && fFit_RHC_1Rmu && fFit_RHC_1Re)
      gOscAnaTypeChoice = kAnaRHCCombined;
    else if(!fFit_1Rmu && !fFit_1Re && fFit_RHC_1Rmu && !fFit_RHC_1Re)
      gOscAnaTypeChoice = kAnaRHCNuMuDisappOnly;
    else if(!fFit_1Rmu && !fFit_1Re && !fFit_RHC_1Rmu && fFit_RHC_1Re)
      gOscAnaTypeChoice = kAnaRHCNuEAppOnly;
    
    // Nothing
    else
      gOscAnaTypeChoice = kAnaUndefined; 
  }
    
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
