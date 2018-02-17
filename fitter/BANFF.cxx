//________________________________________________________________________________________
/*

 \copyright (c) 2010-2016, VaLOR group.

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__  
//#define __show_init_mesg_debug__   
//#define __show_banffsystconv_check_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include "valor/Conventions/Opt.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Registry/Registry.h"
#include "valor/Utilities/UtlSystem.h"

#include "Build.h"
#include "BANFF.h" 
#include "BANFFSystConv.h" 

using namespace valor;
using namespace valor::conv;
using namespace valor::util;

using namespace std;

ClassImp(BANFF)

//________________________________________________________________________________________
BANFF::BANFF() :
CorrelatedSystematics()
{
  // Should rewrite the Conv file so this #ifdef isn't necessary and we only use the FHC branch.
#if defined(__FHC_RUNNING__)&&defined(__RHC_RUNNING__)
  this->Init(kNumOfFHCpRHCParamsInBANFFfit, kListOfFHCpRHCBANFFParams, "syst/banff");
#elif defined(__FHC_RUNNING__)
  this->Init(kNumOfParamsInBANFFfit,kListOfBANFFParams,"syst/banff");
#elif defined(__RHC_RUNNING__)
  this->Init(kNumOfRHCParamsInBANFFfit,kListOfRHCBANFFParams, "syst/banff");
#else
  LOG("BANFF",pFATAL)<<"BANFF can't work out what samples are being used.";
  util::Die(-1);return;
#endif

  this->SetName("BANFF");

  Registry & reg = Opt::Instance()->MutableGlobalRegistry();
  reg.Set("banff", this);
  
  // Somewhat hacky handling of the fact that some of our code still thinks there are RHC specific BANFF-XSEC params.
  for (int irhcp = kRHCNuisBANFFMAQE ; irhcp <= kRHCNuisBANFFEb; irhcp++) {
    int ifhcp = irhcp - (kRHCNuisBANFFMAQE-kNuisBANFFMAQE);
    fParam2Idx[irhcp] = fParam2Idx[ifhcp];
  }
  (*fDefaultValues)[this->Param2Idx(conv::kNuisBANFFShapeCCoth)] = 0;
  
  // (*fDefaultValues)[this->Param2Idx(conv::kNuisBANFFBeRPA_A)]   = 0.59;
  // (*fDefaultValues)[this->Param2Idx(conv::kNuisBANFFBeRPA_B)]   = 1.05;
  // (*fDefaultValues)[this->Param2Idx(conv::kNuisBANFFBeRPA_D)]   = 1.13;
  // (*fDefaultValues)[this->Param2Idx(conv::kNuisBANFFBeRPA_E)]   = 0.88;
  // (*fDefaultValues)[this->Param2Idx(conv::kNuisBANFFBeRPA_U)]   = 1.20;
  
  
  LOG("BANFF",pINFO) << "Defined final BANFF matrix rows: ";
  for (unsigned int i = 0 ; i < fNSystParameters; i++ ) {
    LOG("BANFF",pINFO) << i <<" (default: "<<(*fDefaultValues)[i]<< ")" << "\t" <<util::AsString(fSystParameters[i]);
  }

}
//________________________________________________________________________________________
void BANFF::TweakSampleConstrainedOnly  (bool force_inputs, TVectorD * input_vec, bool unconstrained)
{
  
  this->Tweak(force_inputs, input_vec);
  
  for(unsigned int ip=0;ip<fNSystParameters;ip++){
    PhysParamId_t p = fSystParameters[ip];
    if (unconstrained!=util::IsBANFFUnconstrained(p)) {(*fCurrentValues)[ip]=(*fInitialValues)[ip];};
  }
  return;
}
//________________________________________________________________________________________
void BANFF::TweakSampleFluxOnly  (bool force_inputs, TVectorD * input_vec)
{
  
  this->Tweak(force_inputs, input_vec);

  for(unsigned int ip=0;ip<fNSystParameters;ip++){
    PhysParamId_t p = fSystParameters[ip];
    if (!util::IsFlux(p)) {
      (*fCurrentValues)[ip]=(*fInitialValues)[ip];
    };
  }
  return;
}
//________________________________________________________________________________________
void BANFF::TweakSampleXSecOnly  (bool force_inputs, TVectorD * input_vec)
{
  
  this->Tweak(force_inputs, input_vec);
  
  for(unsigned int ip=0;ip<fNSystParameters;ip++){
    PhysParamId_t p = fSystParameters[ip];
    if (!util::IsXSec(p)) {
      (*fCurrentValues)[ip]=(*fInitialValues)[ip];
    };
  }
  return;
}
//________________________________________________________________________________________
void BANFF::TweakSampleXSecNumuNueOnly  (bool force_inputs, TVectorD * input_vec)
{
  
  this->Tweak(force_inputs, input_vec);
  
  for(unsigned int ip=0;ip<fNSystParameters;ip++){
    PhysParamId_t p = fSystParameters[ip];
    if (!( p == kNuisBANFFNormNuEBarToNuMuBar   ||
          p == kRHCNuisBANFFNormNuEBarToNuMuBar ||
          p == kNuisBANFFNormNuEToNuMu   ||
          p == kRHCNuisBANFFNormNuEToNuMu
        ))
    {
      (*fCurrentValues)[ip]=(*fInitialValues)[ip];
    };
  }
  return;
}

//________________________________________________________________________________________
void BANFF::TweakSample2p2h  (bool force_inputs, TVectorD * input_vec)
{
  
  this->Tweak(force_inputs, input_vec);
  
  for(unsigned int ip=0;ip<fNSystParameters;ip++){
    PhysParamId_t p = fSystParameters[ip];
    if (!( p == kNuisBANFFNorm2p2h || p == kRHCNuisBANFF2p2h)) {(*fCurrentValues)[ip]=(*fInitialValues)[ip];};
  }
  return;
}
//________________________________________________________________________________________
void BANFF::TweakSample2p2hBar  (bool force_inputs, TVectorD * input_vec)
{
  
  this->Tweak(force_inputs, input_vec);
  
  for(unsigned int ip=0;ip<fNSystParameters;ip++){
    PhysParamId_t p = fSystParameters[ip];
    if (!( p == kNuisBANFFNorm2p2hBar || p == kRHCNuisBANFF2p2hBar)) {(*fCurrentValues)[ip]=(*fInitialValues)[ip];};
  }
  return;
}
//________________________________________________________________________________________
void BANFF::TweakSampleNCGamma  (bool force_inputs, TVectorD * input_vec)
{
  
  this->Tweak(force_inputs, input_vec);
  
  for(unsigned int ip=0;ip<fNSystParameters;ip++){
    PhysParamId_t p = fSystParameters[ip];
    if (!( p == kNuisBANFFNormNC1Gamma || p == kRHCNuisBANFFNormNC1Gamma)) {(*fCurrentValues)[ip]=(*fInitialValues)[ip];};
  }
  return;
}
//________________________________________________________________________________________
void BANFF::TweakSampleNCoth  (bool force_inputs, TVectorD * input_vec)
{
  
  this->Tweak(force_inputs, input_vec);
  
  for(unsigned int ip=0;ip<fNSystParameters;ip++){
    PhysParamId_t p = fSystParameters[ip];
    if (!( p == kNuisBANFFNormNCoth || p == kRHCNuisBANFFNormNCoth)) {(*fCurrentValues)[ip]=(*fInitialValues)[ip];};
  }
  return;
}
//________________________________________________________________________________________
BANFF::~BANFF()
{
  
}
//________________________________________________________________________________________
