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

#include "Conventions/Opt.h"
#include "Messenger/Messenger.h"
#include "Registry/Registry.h"
#include "Utilities/UtlSystem.h"

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
  
  LOG("BANFF",pINFO) << "Defined final BANFF matrix rows: ";
  for (unsigned int i = 0 ; i < fNSystParameters; i++ ) {
    LOG("BANFF",pINFO) << "\t"<<util::AsString(fSystParameters[i])<<" (default: "<<(*fDefaultValues)[i]<<")";
  }

  Registry & reg = Opt::Instance()->MutableGlobalRegistry();
  reg.Set("banff", this);
  
  // Somewhat hacky handling of the fact that some of our code still thinks there are RHC specific BANFF-XSEC params.
  for (int irhcp = kRHCNuisBANFFMAQE ; irhcp <= kRHCNuisBANFFMECBar; irhcp++) {
    int ifhcp = irhcp - (kRHCNuisBANFFMAQE-kNuisBANFFMAQE);
    fParam2Idx[irhcp] = fParam2Idx[ifhcp];
  }
  
  (*fDefaultValues)[this->Param2Idx(conv::kNuisBANFFShapeCCoth)] = 0;
}
//________________________________________________________________________________________
BANFF::~BANFF()
{
  
}
//________________________________________________________________________________________
