//________________________________________________________________________________________
/*
 
 Copyright (c) 2010-2014, The VALOR group. All rights reserved.
 
 ChangeLog:
 
 */

#include "Enumerations.h"
#include "SuperkModes.h"
#include "Utilities.h"
#include "Binning.h"

using namespace valor;
using namespace valor::conv;

//___________________________________________________________________________
int valor::util::SuperkMode2FluxNuPdg(valor::conv::SuperkMode_t m)
{
    if ( util::IsNuMu      (m) ) return kPdgNuMu;
    if ( util::IsOscNuMu   (m) ) return kPdgNuE; /// nue -> numu
    if ( util::IsNuMuBar   (m) ) return kPdgNuMuBar;
    if ( util::IsOscNuMuBar(m) ) return kPdgNuEBar; /// nue -> numu
    if ( util::IsNuE       (m) ) return kPdgNuE;
    if ( util::IsNuEBar    (m) ) return kPdgNuEBar;
    if ( util::IsOscNuEBar (m) ) return kPdgNuMuBar; /// nue -> numu
    if ( util::IsOscNuE    (m) ) return kPdgNuMu; /// numu -> nue
    
    // NOTE: Beam NuEBars not handled
    
    return 0;
}
//____________________________________________________________________________
int valor::util::SuperkMode2SKMCPdg(valor::conv::SuperkMode_t m)
{
    if ( util::IsNuMu      (m) ) return kPdgNuMu;
    if ( util::IsOscNuMu   (m) ) return kPdgNuMu;
    if ( util::IsNuMuBar   (m) ) return kPdgNuMuBar;
    if ( util::IsOscNuMuBar(m) ) return kPdgNuMuBar;
    if ( util::IsNuE       (m) ) return kPdgNuE;
    if ( util::IsNuEBar    (m) ) return kPdgNuEBar;
    if ( util::IsOscNuEBar (m) ) return kPdgNuEBar;
    if ( util::IsOscNuE    (m) ) return kPdgNuE;
    
    // NOTE: Beam NuEBars not handled
    
    return 0;
}
//____________________________________________________________________________
bool valor::util::IsNuMu(valor::conv::SuperkMode_t  m)
{
    
    if(m >= kSKmNuMuCCQE && m <= kSKmNuMuNC1Gamma)
        return true;
    
    return false;
}
//____________________________________________________________________________
bool valor::util::IsOscNuMu(
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
  valor::conv::SuperkMode_t  m)
#else
  valor::conv::SuperkMode_t  /*m*/)
#endif
{
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    if(m >= kSKmOscNuMuCCQE && m <= kSKmOscNuMuCCoth)
        return true;
#endif
    return false;
}
//____________________________________________________________________________
bool valor::util::IsOscNuMuBar(
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
  valor::conv::SuperkMode_t  m)
#else
  valor::conv::SuperkMode_t  /*m*/)
#endif
{
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    if(m >= kSKmOscNuMuBarCCQE && m <= kSKmOscNuMuBarCCoth)
        return true;
#endif
    return false;
}
//____________________________________________________________________________
bool valor::util::IsOscNuEBar(valor::conv::SuperkMode_t  m)
{
    if(m >= kSKmOscNuEBarCCQE && m <= kSKmOscNuEBarCCoth)
        return true;
    
    return false;
}
//____________________________________________________________________________
bool valor::util::IsNuMuBar(valor::conv::SuperkMode_t m)
{
    if(m >= kSKmNuMuBarCCQE && m <= kSKmNuMuBarNC1Gamma)
        return true;
    
    return false;
}
//____________________________________________________________________________
bool valor::util::IsNuE(valor::conv::SuperkMode_t m)
{
    if(m >= kSKmNuECCQE && m <= kSKmNuENC1Gamma)
        return true;
    
    return false;
}
//____________________________________________________________________________
bool valor::util::IsNuEBar(valor::conv::SuperkMode_t m)
{
    if(m >= kSKmNuEBarCCQE && m <= kSKmNuEBarNC1Gamma)
        return true;
    
    return false;
}
//____________________________________________________________________________
bool valor::util::IsOscNuE(valor::conv::SuperkMode_t m)
{
    if(m >= kSKmOscNuECCQE && m <= kSKmOscNuECCoth)
        return true;
    return false;
}
//____________________________________________________________________________
bool valor::util::IsCC(valor::conv::SuperkMode_t m)
{
    return !(IsNC(m));
}
//____________________________________________________________________________
bool valor::util::IsCCQE(valor::conv::SuperkMode_t m){
    switch(m) {
        case ( kSKmNuMuCCQE       ) :
        case ( kSKmNuMuBarCCQE    ) :
        case ( kSKmOscNuMuCCQE    ) :
        case ( kSKmOscNuMuBarCCQE ) :
        case ( kSKmNuECCQE        ) :
        case ( kSKmNuEBarCCQE     ) :
        case ( kSKmOscNuECCQE     ) :
        case ( kSKmOscNuEBarCCQE  ) :
            return true;
            break;
        default:
            return false;
            break;
    }
    return false;
}
//____________________________________________________________________________
bool valor::util::IsCCnQE(valor::conv::SuperkMode_t m){
    return ( IsCC(m) && !IsCCQE(m) && !IsMEC(m) );
}
//____________________________________________________________________________
bool valor::util::IsMEC(valor::conv::SuperkMode_t m)
{
    switch(m) {
        case ( kSKmNuMuCCMEC        ) :
        case ( kSKmNuMuBarCCMEC     ) :
        case ( kSKmNuECCMEC         ) :
        case ( kSKmNuEBarCCMEC      ) :
        case ( kSKmOscNuMuCCMEC        ) :
        case ( kSKmOscNuMuBarCCMEC     ) :
        case ( kSKmOscNuECCMEC         ) :
        case ( kSKmOscNuEBarCCMEC      ) :
            return true;
            break;
        default:
            return false;
            break;
    }
    return false;
}
//____________________________________________________________________________
bool valor::util::IsNC(valor::conv::SuperkMode_t m)
{
    switch(m) {
        case ( kSKmNuMuNC1piC      ) :
        case ( kSKmNuMuNC1pi0      ) :
        case ( kSKmNuMuNCcoh       ) :
        case ( kSKmNuMuNCoth       ) :
        case ( kSKmNuMuBarNC1piC   ) :
        case ( kSKmNuMuBarNC1pi0   ) :
        case ( kSKmNuMuBarNCcoh    ) :
        case ( kSKmNuMuBarNCoth    ) :
        case ( kSKmNuENC1piC       ) :
        case ( kSKmNuENC1pi0       ) :
        case ( kSKmNuENCcoh        ) :
        case ( kSKmNuENCoth        ) :
        case ( kSKmNuEBarNC1piC    ) :
        case ( kSKmNuEBarNC1pi0    ) :
        case ( kSKmNuEBarNCcoh     ) :
        case ( kSKmNuEBarNCoth     ) :
      // NC 1 gamma
        case ( kSKmNuMuNC1Gamma    ) :
        case ( kSKmNuMuBarNC1Gamma ) :
        case ( kSKmNuENC1Gamma     ) :          
        case ( kSKmNuEBarNC1Gamma  ) :
             
        return true;
        break;
     default:
        return false;
        break;
    }
    return false;
}
//____________________________________________________________________________
int valor::util::AppliesToCCQE(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  == modes.end())
    return 0;
  //effects all CCQE
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  != modes.end())
    return 9;
  //signue+signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  != modes.end())
    return 10;
  //non-signue+non-signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCQE)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCQE)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCQE)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCQE) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCQE)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCQE)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCQE)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCQE)  == modes.end())
    return 11;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToCCMEC(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  == modes.end())
    return 0;
  //effects all CCMEC
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  != modes.end())
    return 9;
  //signue+signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  != modes.end())
    return 10;
  //non-signue+non-signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCMEC)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCMEC)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCMEC)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCMEC) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCMEC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCMEC)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCMEC)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCMEC)  == modes.end())
    return 11;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToCCcoh(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  == modes.end())
    return 0;
  //effects all CCcoh
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  != modes.end())
    return 9;
  //signue+signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  != modes.end())
    return 10;
  //non-signue+non-signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCcoh)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCcoh)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCcoh)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCcoh) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCcoh)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCcoh)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCcoh)  == modes.end())
    return 11;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToCC1pi(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  == modes.end())
    return 0;
  //effects all CC1pi
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  != modes.end())
    return 9;
  //signue+signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  != modes.end())
    return 10;
  //non-signue+non-signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCC1pi)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCC1pi)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCC1pi)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCC1pi) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECC1pi)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECC1pi)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCC1pi)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCC1pi)  == modes.end())
    return 11;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToCCoth(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  == modes.end())
    return 0;
  //effects all CCoth
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       != modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    != modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) != modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  != modes.end())
    return 9;
  //signue+signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  != modes.end())
    return 10;
  //non-signue+non-signuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuCCoth)       == modes.end() &&
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuCCoth)    == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuMuBarCCoth)    == modes.end() &&
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     find(modes.begin(), modes.end(), kSKmOscNuMuBarCCoth) == modes.end() &&
#endif
     find(modes.begin(), modes.end(), kSKmNuECCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuECCoth)     == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarCCoth)     != modes.end() &&
     find(modes.begin(), modes.end(), kSKmOscNuEBarCCoth)  == modes.end())
    return 11;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToNC1pi0(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     == modes.end())
    return 0;
  //effects all NC1pi0
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1pi0)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1pi0)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1pi0)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1pi0)     != modes.end())
    return 9;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToNC1piC(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     == modes.end())
    return 0;
  //effects all NC1piC
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1piC)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1piC)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1piC)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1piC)     != modes.end())
    return 9;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToNCcoh(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     == modes.end())
    return 0;
  //effects all NCcoh
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCcoh)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCcoh)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCcoh)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCcoh)     != modes.end())
    return 9;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToNC1Gamma(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     == modes.end())
    return 0;
  //effects all NC1Gamma
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNC1Gamma)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNC1Gamma)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENC1Gamma)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNC1Gamma)     != modes.end())
    return 9;
  return 99;
}
//____________________________________________________________________________
int valor::util::AppliesToNCoth(vector<valor::conv::SuperkMode_t> modes)
{
  if(!modes.size())
    return 0;
  //no effect
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     == modes.end())
    return 0;
  //effects all NCoth
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     != modes.end())
    return 1;
  //numu+numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     == modes.end())
    return 2;
  //nue+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     != modes.end())
    return 3;
  //numu+nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     == modes.end())
    return 4;
  //numubar+nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     != modes.end())
    return 5;
  //numu
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     == modes.end())
    return 6;
  //numubar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     == modes.end())
    return 7;
  //nue
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        != modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     == modes.end())
    return 8;
  //nuebar
  if(find(modes.begin(), modes.end(), kSKmNuMuNCoth)       == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuMuBarNCoth)    == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuENCoth)        == modes.end() &&
     find(modes.begin(), modes.end(), kSKmNuEBarNCoth)     != modes.end())
    return 9;
  return 99;
}
//____________________________________________________________________________
void valor::util::PrintAppliesType(int mode, const char * type, bool verbose)
{
  switch(mode)
    {
    case 0:
      if(!verbose)
	return;
      std::cout << "Doesn't effect " << type << std::endl;
      break;
    case 1:
      std::cout << "Effects all of " << type << std::endl;
      break;
    case 2:
      std::cout << "Effects numu+numubar for " << type << std::endl;
      break;
    case 3:
      std::cout << "Effects nue+nuebar for " << type << std::endl;
      break;
    case 4:
      std::cout << "Effects numu+nue for " << type << std::endl;
      break;
    case 5:
      std::cout << "Effects numubar+nuebar for " << type << std::endl;
      break;
    case 6:
      std::cout << "Effects numu for " << type << std::endl;
      break;
    case 7:
      std::cout << "Effects numubar for " << type << std::endl;
      break;
    case 8:
      std::cout << "Effects nue for " << type << std::endl;
      break;
    case 9:
      std::cout << "Effects nuebar for " << type << std::endl;
      break;
    case 10:
      std::cout << "Effects signue+signuebar for " << type << std::endl;
      break;
    case 11:
      std::cout << "Effects non-signue+non-signuebar for " << type << std::endl;
      break;
    default:
      std::cout << "Odd effects for " << type << std::endl;
      break;
    }
}
//____________________________________________________________________________



//____________________________________________________________________________
string valor::util::AsStringCompact_KDI (valor::conv::SuperkMode_t m)
{
    // Compact string used for naming output histograms stored ROOT files
    // using the Japanese analysis naming convention
    //
    // needs updating with MEC modes
    switch(m) {
        case ( kSKmNuMuCCQE        ) : return "numu_ccqe";       break;
        case ( kSKmNuMuCC1pi       ) : return "numu_ccpi";       break;
        case ( kSKmNuMuCCcoh       ) : return "numu_cccoh";      break;
        case ( kSKmNuMuCCoth       ) : return "numu_ccoth";      break;
        case ( kSKmNuMuNC1piC      ) : return "numu_ncpi";       break;
        case ( kSKmNuMuNC1pi0      ) : return "numu_nc1pi0";     break;
        case ( kSKmNuMuNCcoh       ) : return "numu_nccoh";      break;
        case ( kSKmNuMuNCoth       ) : return "numu_ncoth";      break;
        case ( kSKmOscNuMuCCQE     ) : return "do_not_know_kdi_name_numu_ccqe";       break;
        case ( kSKmOscNuMuCC1pi    ) : return "do_not_know_kdi_name_numu_ccpi";       break;
        case ( kSKmOscNuMuCCcoh    ) : return "do_not_know_kdi_name_numu_cccoh";      break;
        case ( kSKmOscNuMuCCoth    ) : return "do_not_know_kdi_name_numu_ccoth";      break;
        case ( kSKmNuMuBarCCQE     ) : return "numb_ccqe";       break;
        case ( kSKmNuMuBarCC1pi    ) : return "numb_ccpi";       break;
        case ( kSKmNuMuBarCCcoh    ) : return "numb_cccoh";      break;
        case ( kSKmNuMuBarCCoth    ) : return "numb_ccoth";      break;
        case ( kSKmNuMuBarNC1piC   ) : return "numb_ncpi";       break;
        case ( kSKmNuMuBarNC1pi0   ) : return "numb_nc1pi0";     break;
        case ( kSKmNuMuBarNCcoh    ) : return "numb_nccoh";      break;
        case ( kSKmNuMuBarNCoth    ) : return "numb_ncoth";      break;
        case ( kSKmNuECCQE         ) : return "nue_ccqe";        break;
        case ( kSKmNuECC1pi        ) : return "nue_ccpi";        break;
        case ( kSKmNuECCcoh        ) : return "nue_cccoh";       break;
        case ( kSKmNuECCoth        ) : return "nue_ccoth";       break;
        case ( kSKmNuENC1piC       ) : return "nue_ncpi";        break;
        case ( kSKmNuENC1pi0       ) : return "nue_nc1pi0";      break;
        case ( kSKmNuENCcoh        ) : return "nue_nccoh";       break;
        case ( kSKmNuENCoth        ) : return "nue_ncoth";       break;
        case ( kSKmNuEBarCCQE      ) : return "nueb_ccqe";       break;
        case ( kSKmNuEBarCC1pi     ) : return "nueb_ccpi";       break;
        case ( kSKmNuEBarCCcoh     ) : return "nueb_cccoh";      break;
        case ( kSKmNuEBarCCoth     ) : return "nueb_ccoth";      break;
        case ( kSKmNuEBarNC1piC    ) : return "nueb_nc1piC";     break;
        case ( kSKmNuEBarNC1pi0    ) : return "nueb_nc1pi0";     break;
        case ( kSKmNuEBarNCcoh     ) : return "nueb_nccoh";      break;
        case ( kSKmNuEBarNCoth     ) : return "nueb_ncoth";      break;
        case ( kSKmOscNuECCQE      ) : return "signue_ccqe";     break;
        case ( kSKmOscNuECC1pi     ) : return "signue_ccpi";     break;
        case ( kSKmOscNuECCcoh     ) : return "signue_cccoh";    break;
        case ( kSKmOscNuECCoth     ) : return "signue_ccoth";    break;
            //NC 1 gamma
        case ( kSKmNuMuNC1Gamma    ) : return "numu_nc1gamma";   break;
        case ( kSKmNuMuBarNC1Gamma ) : return "numubar_nc1gamma";break;
        case ( kSKmNuENC1Gamma     ) : return "nue_nc1gamma";    break;   
        case ( kSKmNuEBarNC1Gamma  ) : return "nuebar_nc1gamma"; break;
        default                      : return "";
    }
    return "";
}
//____________________________________________________________________________
string valor::util::AsString(valor::conv::SuperkMode_t m)
{
    switch(m) {
        case ( kSKmNuMuCCQE        ) : return "numu CC QE";                        break;
        case ( kSKmNuMuCC1pi       ) : return "numu CC 1pi";                       break;
        case ( kSKmNuMuCCcoh       ) : return "numu CC coh";                       break;
        case ( kSKmNuMuCCMEC       ) : return "numu CC MEC";                       break;
        case ( kSKmNuMuCCoth       ) : return "numu CC oth";                       break;
        case ( kSKmNuMuNC1piC      ) : return "numu NC 1piC";                      break;
        case ( kSKmNuMuNC1pi0      ) : return "numu/nutau NC 1pi0";                break;
        case ( kSKmNuMuNCcoh       ) : return "numu/nutau NC coh";                 break;
        case ( kSKmNuMuNCoth       ) : return "numu/nutau NC oth";                 break;
        case ( kSKmOscNuMuCCQE     ) : return "oscillated numu CC QE";             break;
        case ( kSKmOscNuMuCC1pi    ) : return "oscillated numu CC 1pi";            break;
        case ( kSKmOscNuMuCCcoh    ) : return "oscillated numu CC coh";            break;
        case ( kSKmOscNuMuCCMEC    ) : return "oscillated numu CC MEC";            break;
        case ( kSKmOscNuMuCCoth    ) : return "oscillated numu CC oth";            break;
        case ( kSKmNuMuBarCCQE     ) : return "numubar CC QE";                     break;
        case ( kSKmNuMuBarCC1pi    ) : return "numubar CC 1pi";                    break;
        case ( kSKmNuMuBarCCcoh    ) : return "numubar CC coh";                    break;
        case ( kSKmNuMuBarCCMEC    ) : return "numubar CC MEC";                    break;
        case ( kSKmNuMuBarCCoth    ) : return "numubar CC oth";                    break;
        case ( kSKmNuMuBarNC1piC   ) : return "numubar/nutaubar NC 1piC";          break;
        case ( kSKmNuMuBarNC1pi0   ) : return "numubar/nutaubar NC 1pi0";          break;
        case ( kSKmNuMuBarNCcoh    ) : return "numubar/nutaubar NC coh";           break;
        case ( kSKmNuMuBarNCoth    ) : return "numubar/nutaubar NC oth";           break;
        case ( kSKmOscNuMuBarCCQE  ) : return "oscillated numubar CC QE";          break;
        case ( kSKmOscNuMuBarCC1pi ) : return "oscillated numubar CC 1pi";         break;
        case ( kSKmOscNuMuBarCCcoh ) : return "oscillated numubar CC coh";         break;
        case ( kSKmOscNuMuBarCCMEC ) : return "oscillated numubar CC MEC";         break;
        case ( kSKmOscNuMuBarCCoth ) : return "oscillated numubar CC oth";         break;
        case ( kSKmNuECCQE         ) : return "nue CC QE";                         break;
        case ( kSKmNuECC1pi        ) : return "nue CC 1pi";                        break;
        case ( kSKmNuECCcoh        ) : return "nue CC coh";                        break;
        case ( kSKmNuECCMEC        ) : return "nue CC MEC";                        break;
        case ( kSKmNuECCoth        ) : return "nue CC oth";                        break;
        case ( kSKmNuENC1piC       ) : return "nue NC 1piC";                       break;
        case ( kSKmNuENC1pi0       ) : return "nue NC 1pi0";                       break;
        case ( kSKmNuENCcoh        ) : return "nue NC coh";                        break;
        case ( kSKmNuENCoth        ) : return "nue NC oth";                        break;
        case ( kSKmOscNuECCQE      ) : return "oscillated nue CC QE";              break;
        case ( kSKmOscNuECC1pi     ) : return "oscillated nue CC 1pi";             break;
        case ( kSKmOscNuECCcoh     ) : return "oscillated nue CC coh";             break;
        case ( kSKmOscNuECCMEC     ) : return "oscillated nue CC MEC";             break;
        case ( kSKmOscNuECCoth     ) : return "oscillated nue CC oth";             break;
        case ( kSKmNuEBarCCQE      ) : return "nuebar CC QE";                      break;
        case ( kSKmNuEBarCC1pi     ) : return "nuebar CC 1pi";                     break;
        case ( kSKmNuEBarCCcoh     ) : return "nuebar CC coh";                     break;
        case ( kSKmNuEBarCCMEC     ) : return "nuebar CC MEC";                     break;
        case ( kSKmNuEBarCCoth     ) : return "nuebar CC oth";                     break;
        case ( kSKmNuEBarNC1piC    ) : return "nuebar NC 1piC";                    break;
        case ( kSKmNuEBarNC1pi0    ) : return "nuebar NC 1pi0";                    break;
        case ( kSKmNuEBarNCcoh     ) : return "nuebar NC coh";                     break;
        case ( kSKmNuEBarNCoth     ) : return "nuebar NC oth";                     break;
        case ( kSKmOscNuEBarCCQE   ) : return "oscillated nuebar CC QE";           break;
        case ( kSKmOscNuEBarCC1pi  ) : return "oscillated nuebar CC 1pi";          break;
        case ( kSKmOscNuEBarCCcoh  ) : return "oscillated nuebar CC coh";          break;
        case ( kSKmOscNuEBarCCMEC  ) : return "oscillated nuebar CC MEC";          break;
        case ( kSKmOscNuEBarCCoth  ) : return "oscillated nuebar CC oth ";         break;
            //NC 1 gamma
        case ( kSKmNuMuNC1Gamma      ) : return "numu NC 1gamma";                    break;
        case ( kSKmNuMuBarNC1Gamma   ) : return "numubar NC 1gamma";                 break;          
        case ( kSKmNuENC1Gamma       ) : return "nue NC 1gamma";                     break;
        case ( kSKmNuEBarNC1Gamma    ) : return "nuebar NC 1gamma";                  break;

        default                      : return "--";
    }
    return "";
}
//____________________________________________________________________________
string valor::util::AsStringRooTex(valor::conv::SuperkMode_t m)
{
    string s = AsStringTex(m);
    replace( s.begin(), s.end(), '\\', '#' ); // Make slashes hashes for RooTex
    s.erase( std::remove( s.begin(), s.end(), '$'), s.end() ) ; // Remove dollar signs.
    return s;
}
//________________________________________________________________________________________
string valor::util::AsStringTex(valor::conv::SuperkMode_t m)
{
    switch(m) {
            
        case ( kSKmNuMuCCQE        ) : return "$\\nu_\\mu$ CCQE";                                               break;
        case ( kSKmNuMuCC1pi       ) : return "$\\nu_\\mu$ CC 1$\\pi$";                                         break;
        case ( kSKmNuMuCCcoh       ) : return "$\\nu_\\mu$ CC coherent";                                        break;
        case ( kSKmNuMuCCMEC       ) : return "$\\nu_\\mu$ CC MEC";                                             break;
        case ( kSKmNuMuCCoth       ) : return "$\\nu_\\mu$ CC other";                                           break;
        case ( kSKmNuMuNC1piC      ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC 1$\\pi^{+/-}$";                      break;
        case ( kSKmNuMuNC1pi0      ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC 1$\\pi^0$";                          break;
        case ( kSKmNuMuNCcoh       ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC coherent";                           break;
        case ( kSKmNuMuNCoth       ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC other";                              break;
        case ( kSKmOscNuMuCCQE     ) : return "Osc. $\\nu_\\mu$ CCQE";                                          break;
        case ( kSKmOscNuMuCC1pi    ) : return "Osc. $\\nu_\\mu$ CC 1$\\pi$";                                    break;
        case ( kSKmOscNuMuCCcoh    ) : return "Osc. $\\nu_\\mu$ CC coherent";                                   break;
        case ( kSKmOscNuMuCCMEC    ) : return "Osc. $\\nu_\\mu$ CC MEC";                                        break;
        case ( kSKmOscNuMuCCoth    ) : return "Osc. $\\nu_\\mu$ CC other";                                      break;
        case ( kSKmNuMuBarCCQE     ) : return "$\\bar{\\nu}_\\mu$ CCQE";                                        break;
        case ( kSKmNuMuBarCC1pi    ) : return "$\\bar{\\nu}_\\mu$ CC 1$\\pi$";                                  break;
        case ( kSKmNuMuBarCCcoh    ) : return "$\\bar{\\nu}_\\mu$ CC coherent";                                 break;
        case ( kSKmNuMuBarCCMEC    ) : return "$\\bar{\\nu}_\\mu$ CC MEC";                                      break;
        case ( kSKmNuMuBarCCoth    ) : return "$\\bar{\\nu}_\\mu$ CC other";                                    break;
        case ( kSKmNuMuBarNC1piC   ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC 1$\\pi^{+/-}$";        break;
        case ( kSKmNuMuBarNC1pi0   ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC 1$\\pi^0$";            break;
        case ( kSKmNuMuBarNCcoh    ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC coherent";             break;
        case ( kSKmNuMuBarNCoth    ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC other";                break;
        case ( kSKmOscNuMuBarCCQE  ) : return "Osc. $\\bar{\\nu}_\\mu$ CCQE";                                   break;
        case ( kSKmOscNuMuBarCC1pi ) : return "Osc. $\\bar{\\nu}_\\mu$ CC 1$\\pi$";                             break;
        case ( kSKmOscNuMuBarCCcoh ) : return "Osc. $\\bar{\\nu}_\\mu$ CC coherent";                            break;
        case ( kSKmOscNuMuBarCCMEC ) : return "Osc. $\\bar{\\nu}_\\mu$ CC MEC";                                 break;
        case ( kSKmOscNuMuBarCCoth ) : return "Osc. $\\bar{\\nu}_\\mu$ CC other";                               break;
        case ( kSKmNuECCQE         ) : return "$\\nu_e$ CCQE";                                                  break;
        case ( kSKmNuECC1pi        ) : return "$\\nu_e$ CC 1$\\pi$";                                            break;
        case ( kSKmNuECCcoh        ) : return "$\\nu_e$ CC coherent";                                           break;
        case ( kSKmNuECCMEC        ) : return "$\\nu_e$ CC MEC";                                                break;
        case ( kSKmNuECCoth        ) : return "$\\nu_e$ CC other";                                              break;
        case ( kSKmNuENC1piC       ) : return "$\\nu_e$ NC 1$\\pi^{+/-}$";                                      break;
        case ( kSKmNuENC1pi0       ) : return "$\\nu_e$ NC 1$\\pi^0$";                                          break;
        case ( kSKmNuENCcoh        ) : return "$\\nu_e$ NC coherent";                                           break;
        case ( kSKmNuENCoth        ) : return "$\\nu_e$ NC other";                                              break;
        case ( kSKmOscNuECCQE      ) : return "Osc. $\\nu_e$ CCQE";                                             break;
        case ( kSKmOscNuECC1pi     ) : return "Osc. $\\nu_e$ CC 1$\\pi$";                                       break;
        case ( kSKmOscNuECCcoh     ) : return "Osc. $\\nu_e$ CC coherent";                                      break;
        case ( kSKmOscNuECCMEC     ) : return "Osc. $\\nu_e$ CC MEC";                                           break;
        case ( kSKmOscNuECCoth     ) : return "Osc. $\\nu_e$ CC other";                                         break;
        case ( kSKmNuEBarCCQE      ) : return "$\\bar{\\nu}_e$ CCQE";                                           break;
        case ( kSKmNuEBarCC1pi     ) : return "$\\bar{\\nu}_e$ CC 1$\\pi$";                                     break;
        case ( kSKmNuEBarCCcoh     ) : return "$\\bar{\\nu}_e$ CC coherent";                                    break;
        case ( kSKmNuEBarCCMEC     ) : return "$\\bar{\\nu}_e$ CC MEC";                                         break;
        case ( kSKmNuEBarCCoth     ) : return "$\\bar{\\nu}_e$ CC other";                                       break;
        case ( kSKmNuEBarNC1piC    ) : return "$\\bar{\\nu}_e$ NC 1$\\pi^{+/-}$";                               break;
        case ( kSKmNuEBarNC1pi0    ) : return "$\\bar{\\nu}_e$ NC 1$\\pi^0$";                                   break;
        case ( kSKmNuEBarNCcoh     ) : return "$\\bar{\\nu}_e$ NC coherent";                                    break;
        case ( kSKmNuEBarNCoth     ) : return "$\\bar{\\nu}_e$ NC other";                                       break;
        case ( kSKmOscNuEBarCCQE   ) : return "Osc. $\\bar{\\nu}_e$ CCQE";                                      break;
        case ( kSKmOscNuEBarCC1pi  ) : return "Osc. $\\bar{\\nu}_e$ CC 1$\\pi$";                                break;
        case ( kSKmOscNuEBarCCcoh  ) : return "Osc. $\\bar{\\nu}_e$ CC coherent";                               break;
        case ( kSKmOscNuEBarCCMEC  ) : return "Osc. $\\bar{\\nu}_e$ CC MEC";                                    break;
        case ( kSKmOscNuEBarCCoth  ) : return "Osc. $\\bar{\\nu}_e$ CC other";                                  break;
       // NC 1 gamma
        case ( kSKmNuMuNC1Gamma    ) : return "$\\nu_{\\mu}$ NC 1$\\gamma$";                                    break;
        case ( kSKmNuMuBarNC1Gamma ) : return "$\\bar{\\nu}_\\mu$ NC 1$\\gamma$";                               break;
        case ( kSKmNuENC1Gamma     ) : return "$\\nu_e$ NC 1$\\gamma$";                                         break;          
        case ( kSKmNuEBarNC1Gamma  ) : return "$\\bar{\\nu}_e$ NC 1$\\gamma$";                                  break;

        default                      : return "--";
    }
    return "";
}
//____________________________________________________________________________
string valor::util::AsStringCompact (valor::conv::SuperkMode_t m)
{
    // Compact string used for naming output histograms stored ROOT files
    //
    switch(m) {
        case ( kSKmNuMuCCQE        ) : return "numuccqe";          break;
        case ( kSKmNuMuCC1pi       ) : return "numucc1pi";         break;
        case ( kSKmNuMuCCcoh       ) : return "numucccoh";         break;
        case ( kSKmNuMuCCMEC       ) : return "numuccmec";         break;
        case ( kSKmNuMuCCoth       ) : return "numuccoth";         break;
        case ( kSKmNuMuNC1piC      ) : return "numunc1piC";        break;
        case ( kSKmNuMuNC1pi0      ) : return "numunc1pi0";        break;
        case ( kSKmNuMuNCcoh       ) : return "numunccoh";         break;
        case ( kSKmNuMuNCoth       ) : return "numuncoth";         break;
        case ( kSKmOscNuMuCCQE     ) : return "oscnumuccqe";       break;
        case ( kSKmOscNuMuCC1pi    ) : return "oscnumucc1pi";      break;
        case ( kSKmOscNuMuCCcoh    ) : return "oscnumucccoh";      break;
        case ( kSKmOscNuMuCCMEC    ) : return "oscnumuccmec";      break;
        case ( kSKmOscNuMuCCoth    ) : return "oscnumuccoth";      break;
        case ( kSKmNuMuBarCCQE     ) : return "numubarccqe";       break;
        case ( kSKmNuMuBarCC1pi    ) : return "numubarcc1pi";      break;
        case ( kSKmNuMuBarCCcoh    ) : return "numubarcccoh";      break;
        case ( kSKmNuMuBarCCMEC    ) : return "numubarccmec";      break;
        case ( kSKmNuMuBarCCoth    ) : return "numubarccoth";      break;
        case ( kSKmNuMuBarNC1piC   ) : return "numubarnc1piC";     break;
        case ( kSKmNuMuBarNC1pi0   ) : return "numubarnc1pi0";     break;
        case ( kSKmNuMuBarNCcoh    ) : return "numubarnccoh";      break;
        case ( kSKmNuMuBarNCoth    ) : return "numubarncoth";      break;
        case ( kSKmOscNuMuBarCCQE  ) : return "oscnumubarccqe";    break;
        case ( kSKmOscNuMuBarCC1pi ) : return "oscnumubarcc1pi";   break;
        case ( kSKmOscNuMuBarCCcoh ) : return "oscnumubarcccoh";   break;
        case ( kSKmOscNuMuBarCCMEC ) : return "oscnumubarccmec";   break;
        case ( kSKmOscNuMuBarCCoth ) : return "oscnumubarccoth";   break;
        case ( kSKmNuECCQE         ) : return "nueccqe";           break;
        case ( kSKmNuECC1pi        ) : return "nuecc1pi";          break;
        case ( kSKmNuECCcoh        ) : return "nuecccoh";          break;
        case ( kSKmNuECCMEC        ) : return "nueccmec";          break;
        case ( kSKmNuECCoth        ) : return "nueccoth";          break;
        case ( kSKmNuENC1piC       ) : return "nuenc1piC";         break;
        case ( kSKmNuENC1pi0       ) : return "nuenc1pi0";         break;
        case ( kSKmNuENCcoh        ) : return "nuenccoh";          break;
        case ( kSKmNuENCoth        ) : return "nuencoth";          break;
        case ( kSKmOscNuECCQE      ) : return "oscnueccqe";        break;
        case ( kSKmOscNuECC1pi     ) : return "oscnuecc1pi";       break;
        case ( kSKmOscNuECCcoh     ) : return "oscnuecccoh";       break;
        case ( kSKmOscNuECCMEC     ) : return "oscnueccmec";       break;
        case ( kSKmOscNuECCoth     ) : return "oscnueccoth";       break;
        case ( kSKmNuEBarCCQE      ) : return "nuebarccqe";        break;
        case ( kSKmNuEBarCC1pi     ) : return "nuebarcc1pi";       break;
        case ( kSKmNuEBarCCcoh     ) : return "nuebarcccoh";       break;
        case ( kSKmNuEBarCCMEC     ) : return "nuebarccmec";       break;
        case ( kSKmNuEBarCCoth     ) : return "nuebarccoth";       break;
        case ( kSKmNuEBarNC1piC    ) : return "nuebarnc1piC";      break;
        case ( kSKmNuEBarNC1pi0    ) : return "nuebarnc1pi0";      break;
        case ( kSKmNuEBarNCcoh     ) : return "nuebarnccoh";       break;
        case ( kSKmNuEBarNCoth     ) : return "nuebarncoth";       break;
        case ( kSKmOscNuEBarCCQE   ) : return "oscnuebarccqe";     break;
        case ( kSKmOscNuEBarCC1pi  ) : return "oscnuebarcc1pi";    break;
        case ( kSKmOscNuEBarCCcoh  ) : return "oscnuebarcccoh";    break;
        case ( kSKmOscNuEBarCCMEC  ) : return "oscnuebarccmec";    break;
        case ( kSKmOscNuEBarCCoth  ) : return "oscnuebarccoth";    break;
    // NC 1 gamma
        case ( kSKmNuMuNC1Gamma    ) : return "numunc1gamma";      break;  
        case ( kSKmNuMuBarNC1Gamma ) : return "numubarnc1gamma";   break;
        case ( kSKmNuENC1Gamma     ) : return "nuenc1gamma";       break;
        case ( kSKmNuEBarNC1Gamma  ) : return "nuebarnc1gamma";    break;

        default                      : return "";
    }
    return "";
}
//____________________________________________________________________________
valor::conv::SuperkMode_t valor::util::SuperkModeFromString(string sh)
{
    string lsh = util::ToLowerCase(sh);
    for (int i = 0 ; i < kNSuperkModes ; i ++ ) {
        SuperkMode_t m = (SuperkMode_t) i;
        string full    = util::ToLowerCase(AsString(m));
        string compact = util::ToLowerCase(AsStringCompact(m));
        if (m != kSKmUndefined && (lsh == full || lsh == compact) ) return m;
    }
    std::cerr<<"SuperkModeFromString >> Did not get mode for string: "<<sh<<std::endl;
    return kSKmUndefined;
}
//____________________________________________________________________________
string valor::util::AsStringCompact_splinehack (valor::conv::SuperkMode_t m)
{
    // TEMPORARY HACK using the old naming scheme for osc nue
    //
    switch(m) {
        case ( kSKmOscNuECCQE      ) : return "signueccqe";        break;
        case ( kSKmOscNuECC1pi     ) : return "signuecc1pi";       break;
        case ( kSKmOscNuECCcoh     ) : return "signuecccoh";       break;
        case ( kSKmOscNuECCoth     ) : return "signueccoth";       break;
            
        default                      : return AsStringCompact(m);
    }
    return AsStringCompact(m);
}
//___________________________________________________________________________
string valor::util::ModeAsString       (conv::SuperkMode_t m) { return AsString       (m); }
string valor::util::ModeAsStringCompact(conv::SuperkMode_t m) { return AsStringCompact(m); }
string valor::util::ModeAsStringTex    (conv::SuperkMode_t m) { return AsStringTex    (m); }
string valor::util::ModeAsStringRooTex (conv::SuperkMode_t m) { return AsStringRooTex (m); }
//___________________________________________________________________________
valor::conv::SuperkMode_t valor::util::GetMode(int is, int im){
    return valor::conv::kSKReactionModes[is][im];
}
//___________________________________________________________________________
