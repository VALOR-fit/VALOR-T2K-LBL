//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/

#include "Enumerations.h"
#include "SuperkModes.h"
#include "Utilities.h"
#include "Binning.h"

using namespace t2kana;
using namespace t2kana::conv;

//___________________________________________________________________________
int t2kana::util::SuperkMode2FluxNuPdg(t2kana::conv::SuperkMode_t m)
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
int t2kana::util::SuperkMode2SKMCPdg(t2kana::conv::SuperkMode_t m)
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
bool t2kana::util::IsNuMu(t2kana::conv::SuperkMode_t  m)
{

  if(m >= kSKmNuMuCCQE && m <= kSKmNuMuNCoth_e)
    return true;
  /*
  switch(m) {
     case ( kSKmNuMuCCQE        ) :
     case ( kSKmMEC             ) : 
     case ( kSKmNuMuCC1pi       ) :
     case ( kSKmNuMuCCcoh       ) :
     case ( kSKmNuMuCCoth    ) :
       //     case ( kSKmNuMuCCoth     ) :
     case ( kSKmNuMuNC1piC      ) :
     case ( kSKmNuMuNC1pi0      ) :
     case ( kSKmNuMuNCcoh       ) :
     case ( kSKmNuMuNCoth_e     ) :
     case ( kSKmNuMuNCoth_mu    ) :
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsOscNuMu(t2kana::conv::SuperkMode_t  m)
{
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
  if(m >= kSKmOscNuMuCCQE && m <= kSKmOscNuMuCCoth)
    return true;
#endif
  /*
  switch(m) {
     case ( kSKmOscNuMuCCQE     ) : 
     case ( kSKmOscNuMuCC1pi    ) :
     case ( kSKmOscNuMuCCcoh    ) :
       //     case ( kSKmOscNuMuCCoth ) :
     case ( kSKmOscNuMuCCoth  ) :
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsOscNuMuBar(t2kana::conv::SuperkMode_t  m)
{
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
  if(m >= kSKmOscNuMuBarCCQE && m <= kSKmOscNuMuBarCCoth)
    return true;
#endif
  /*
  switch(m) {
     case ( kSKmOscNuMuBarCCQE     ) : 
     case ( kSKmOscNuMuBarCC1pi    ) :
     case ( kSKmOscNuMuBarCCcoh    ) :
       //     case ( kSKmOscNuMuCCoth ) :
     case ( kSKmOscNuMuBarCCoth  ) :
       return true;
       break;
     default :                
       return false;
       break;
       }*/
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsOscNuEBar(t2kana::conv::SuperkMode_t  m)
{
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
  if(m >= kSKmOscNuEBarCCQE && m <= kSKmOscNuEBarCCoth)
    return true;
#endif
  /*switch(m) {
     case ( kSKmOscNuEBarCCQE     ) : 
     case ( kSKmOscNuEBarCC1pi    ) :
     case ( kSKmOscNuEBarCCcoh    ) :
       //     case ( kSKmOscNuMuCCoth ) :
     case ( kSKmOscNuEBarCCoth  ) :
       return true;
       break;
     default :                
       return false;
       break;
       }*/
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuMuBar(t2kana::conv::SuperkMode_t m)
{
  if(m >= kSKmNuMuBarCCQE && m <= kSKmNuMuBarNCoth_e)
    return true;
  /*
  switch(m) {
     case (kSKmNuMuBarCCQE    ):
     case (kSKmNuMuBarCC1pi   ):
     case (kSKmNuMuBarCCcoh   ):
       //     case (kSKmNuMuBarCCoth):
     case (kSKmNuMuBarCCoth ):
     case (kSKmNuMuBarNC1piC  ):
     case (kSKmNuMuBarNC1pi0  ):
     case (kSKmNuMuBarNCcoh   ):
     case (kSKmNuMuBarNCoth_e ):
     case (kSKmNuMuBarNCoth_mu):
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuE(t2kana::conv::SuperkMode_t m)
{
  if(m >= kSKmNuECCQE && m <= kSKmNuENCoth_e)
    return true;
  /*
  switch(m) {
     case (kSKmNuECCQE    ):
     case (kSKmNuECC1pi   ):
     case (kSKmNuECCcoh   ):
       //     case (kSKmNuECCoth):
     case (kSKmNuECCoth ):
     case (kSKmNuENC1piC  ):
     case (kSKmNuENC1pi0  ):
     case (kSKmNuENCcoh   ):
     case (kSKmNuENCoth_e ):
     case (kSKmNuENCoth_mu):
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuEBar(t2kana::conv::SuperkMode_t m)
{
  if(m >= kSKmNuEBarCCQE && m <= kSKmNuEBarNCoth_e)
    return true;
  /*
  switch(m) {
     case (kSKmNuEBarCCQE    ):
     case (kSKmNuEBarCC1pi   ):
     case (kSKmNuEBarCCcoh   ):
       //     case (kSKmNuEBarCCoth):
     case (kSKmNuEBarCCoth ):
     case (kSKmNuEBarNC1piC  ):
     case (kSKmNuEBarNC1pi0  ):
     case (kSKmNuEBarNCcoh   ):
     case (kSKmNuEBarNCoth_e ):
     case (kSKmNuEBarNCoth_mu):
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsOscNuE(t2kana::conv::SuperkMode_t m)
{
  if(m >= kSKmOscNuECCQE && m <= kSKmOscNuECCoth)
    return true;
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsCC(t2kana::conv::SuperkMode_t m)
{
  return !(IsNC(m));
}
//____________________________________________________________________________
bool t2kana::util::IsNC(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuNC1piC      ) : 
     case ( kSKmNuMuNC1pi0      ) : 
     case ( kSKmNuMuNCcoh       ) : 
     case ( kSKmNuMuNCoth_e     ) : 
     case ( kSKmNuMuNCoth_mu    ) : 
     case ( kSKmNuMuBarNC1piC   ) : 
     case ( kSKmNuMuBarNC1pi0   ) : 
     case ( kSKmNuMuBarNCcoh    ) : 
     case ( kSKmNuMuBarNCoth_e  ) : 
     case ( kSKmNuMuBarNCoth_mu ) : 
     case ( kSKmNuENC1piC       ) : 
     case ( kSKmNuENC1pi0       ) : 
     case ( kSKmNuENCcoh        ) : 
     case ( kSKmNuENCoth_e      ) : 
     case ( kSKmNuENCoth_mu     ) : 
     case ( kSKmNuEBarNC1piC    ) : 
     case ( kSKmNuEBarNC1pi0    ) : 
     case ( kSKmNuEBarNCcoh     ) : 
     case ( kSKmNuEBarNCoth_e   ) : 
     case ( kSKmNuEBarNCoth_mu  ) : 
        return true;
        break;
     default: 
        return false;
        break;
  }
  return false;
}
//____________________________________________________________________________



//____________________________________________________________________________
string t2kana::util::AsStringCompact_KDI (t2kana::conv::SuperkMode_t m)
{
// Compact string used for naming output histograms stored ROOT files
// using the Japanese analysis naming convention
//
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numu_ccqe";       break;
     case ( kSKmMEC             ) : return "mec";             break;
     case ( kSKmNuMuCC1pi       ) : return "numu_ccpi";       break;
     case ( kSKmNuMuCCcoh       ) : return "numu_cccoh";      break;
     case ( kSKmNuMuCCoth       ) : return "numu_ccoth";      break;
     case ( kSKmNuMuNC1piC      ) : return "numu_ncpi";       break;
     case ( kSKmNuMuNC1pi0      ) : return "numu_nc1pi0";     break;
     case ( kSKmNuMuNCcoh       ) : return "numu_nccoh";      break;
     case ( kSKmNuMuNCoth_e     ) : return "numu_ncoth";      break;
     case ( kSKmNuMuNCoth_mu    ) : return "numu_ncoth";      break;
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
     case ( kSKmNuMuBarNCoth_e  ) : return "numb_ncoth";      break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numb_ncoth";      break;
     case ( kSKmNuECCQE         ) : return "nue_ccqe";        break;
     case ( kSKmNuECC1pi        ) : return "nue_ccpi";        break;
     case ( kSKmNuECCcoh        ) : return "nue_cccoh";       break;
     case ( kSKmNuECCoth        ) : return "nue_ccoth";       break;
     case ( kSKmNuENC1piC       ) : return "nue_ncpi";        break;
     case ( kSKmNuENC1pi0       ) : return "nue_nc1pi0";      break;
     case ( kSKmNuENCcoh        ) : return "nue_nccoh";       break;
     case ( kSKmNuENCoth_e      ) : return "nue_ncoth";       break;
     case ( kSKmNuENCoth_mu     ) : return "nue_ncoth";       break;
     case ( kSKmNuEBarCCQE      ) : return "nueb_ccqe";       break;
     case ( kSKmNuEBarCC1pi     ) : return "nueb_ccpi";       break;
     case ( kSKmNuEBarCCcoh     ) : return "nueb_cccoh";      break;
     case ( kSKmNuEBarCCoth     ) : return "nueb_ccoth";      break;
     case ( kSKmNuEBarNC1piC    ) : return "nueb_ncpi";       break;
     case ( kSKmNuEBarNC1pi0    ) : return "nueb_nc1pi0";     break;
     case ( kSKmNuEBarNCcoh     ) : return "nueb_nccoh";      break;
     case ( kSKmNuEBarNCoth_e   ) : return "nueb_ncoth";      break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nueb_ncoth";      break;
     case ( kSKmOscNuECCQE      ) : return "signue_ccqe";     break;
     case ( kSKmOscNuECC1pi     ) : return "signue_ccpi";     break;
     case ( kSKmOscNuECCcoh     ) : return "signue_cccoh";    break;
     case ( kSKmOscNuECCoth     ) : return "signue_ccoth";    break;
     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numu CC QE";                        break;
     case ( kSKmMEC             ) : return "numu CC MEC";                       break;
     case ( kSKmNuMuCC1pi       ) : return "numu CC 1pi";                       break;
     case ( kSKmNuMuCCcoh       ) : return "numu CC coh";                       break;
     case ( kSKmNuMuCCoth       ) : return "numu CC oth";                       break;
     case ( kSKmNuMuNC1piC      ) : return "numu NC 1piC";                      break;
     case ( kSKmNuMuNC1pi0      ) : return "numu/nutau NC 1pi0";                break;
     case ( kSKmNuMuNCcoh       ) : return "numu/nutau NC coh";                 break;
     case ( kSKmNuMuNCoth_e     ) : return "numu/nutau NC oth (nue app)";       break;
     case ( kSKmNuMuNCoth_mu    ) : return "numu/nutau NC oth (numu disapp)";   break;
     case ( kSKmOscNuMuCCQE     ) : return "oscillated numu CC QE";             break;
     case ( kSKmOscNuMuCC1pi    ) : return "oscillated numu CC 1pi";            break;  
     case ( kSKmOscNuMuCCcoh    ) : return "oscillated numu CC coh";            break;
     case ( kSKmOscNuMuCCoth    ) : return "oscillated numu CC oth";            break;
     case ( kSKmNuMuBarCCQE     ) : return "numubar CC QE";                     break;
     case ( kSKmNuMuBarCC1pi    ) : return "numubar CC 1pi";                    break;
     case ( kSKmNuMuBarCCcoh    ) : return "numubar CC coh";                    break;
     case ( kSKmNuMuBarCCoth    ) : return "numubar CC oth";                    break;
     case ( kSKmNuMuBarNC1piC   ) : return "numubar/nutaubar NC 1piC";          break;
     case ( kSKmNuMuBarNC1pi0   ) : return "numubar/nutaubar NC 1pi0";          break;
     case ( kSKmNuMuBarNCcoh    ) : return "numubar/nutaubar NC coh";           break;
     case ( kSKmNuMuBarNCoth_e  ) : return "numubar/nutaubar NC oth (nue app)"; break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numubar/nutaubar NC oth (numu disapp)"; break;
     case ( kSKmOscNuMuBarCCQE  ) : return "oscillated numubar CC QE";          break;
     case ( kSKmOscNuMuBarCC1pi ) : return "oscillated numubar CC 1pi";         break;  
     case ( kSKmOscNuMuBarCCcoh ) : return "oscillated numubar CC coh";         break;
     case ( kSKmOscNuMuBarCCoth ) : return "oscillated numubar CC oth";         break;
     case ( kSKmNuECCQE         ) : return "nue CC QE";                         break;
     case ( kSKmNuECC1pi        ) : return "nue CC 1pi";                        break;
     case ( kSKmNuECCcoh        ) : return "nue CC coh";                        break;
     case ( kSKmNuECCoth        ) : return "nue CC oth";                        break;
     case ( kSKmNuENC1piC       ) : return "nue NC 1piC";                       break;
     case ( kSKmNuENC1pi0       ) : return "nue NC 1pi0";                       break;
     case ( kSKmNuENCcoh        ) : return "nue NC coh";                        break;
     case ( kSKmNuENCoth_e      ) : return "nue NC oth (nue app)";              break;
     case ( kSKmNuENCoth_mu     ) : return "nue NC oth (numu disapp)";          break;
     case ( kSKmOscNuECCQE      ) : return "oscillated nue CC QE";              break;
     case ( kSKmOscNuECC1pi     ) : return "oscillated nue CC 1pi";             break;
     case ( kSKmOscNuECCcoh     ) : return "oscillated nue CC coh";             break;
     case ( kSKmOscNuECCoth     ) : return "oscillated nue CC oth";             break;
     case ( kSKmNuEBarCCQE      ) : return "nuebar CC QE";                      break;
     case ( kSKmNuEBarCC1pi     ) : return "nuebar CC 1pi";                     break;
     case ( kSKmNuEBarCCcoh     ) : return "nuebar CC coh";                     break;
     case ( kSKmNuEBarCCoth     ) : return "nuebar CC oth";                     break;
     case ( kSKmNuEBarNC1piC    ) : return "nuebar NC 1piC";                    break;
     case ( kSKmNuEBarNC1pi0    ) : return "nuebar NC 1pi0";                    break;
     case ( kSKmNuEBarNCcoh     ) : return "nuebar NC coh";                     break;
     case ( kSKmNuEBarNCoth_e   ) : return "nuebar NC oth (nue app)";           break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nuebar NC oth (numu disapp)";       break;
     case ( kSKmOscNuEBarCCQE   ) : return "oscillated nuebar CC QE";           break;
     case ( kSKmOscNuEBarCC1pi  ) : return "oscillated nuebar CC 1pi";          break;  
     case ( kSKmOscNuEBarCCcoh  ) : return "oscillated nuebar CC coh";          break;
     case ( kSKmOscNuEBarCCoth  ) : return "oscillated nuebar CC oth ";         break;

     default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringRooTex(t2kana::conv::SuperkMode_t m)
{
  string s = AsStringTex(m);
  replace( s.begin(), s.end(), '\\', '#' ); // Make slashes hashes for RooTex
  s.erase( std::remove( s.begin(), s.end(), '$'), s.end() ) ; // Remove dollar signs.
  return s;
}
//________________________________________________________________________________________
string t2kana::util::AsStringTex(t2kana::conv::SuperkMode_t m)
{
  switch(m) {

     case ( kSKmNuMuCCQE        ) : return "$\\nu_\\mu$ CCQE";                                               break;
     case ( kSKmMEC             ) : return "$\\nu_\\mu$ CC MEC";                                             break;
     case ( kSKmNuMuCC1pi       ) : return "$\\nu_\\mu$ CC 1$\\pi$";                                         break;
     case ( kSKmNuMuCCcoh       ) : return "$\\nu_\\mu$ CC coherent";                                        break;
     case ( kSKmNuMuCCoth       ) : return "$\\nu_\\mu$ CC other";                                           break;
     case ( kSKmNuMuNC1piC      ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC 1$\\pi^{+/-}$";                      break;
     case ( kSKmNuMuNC1pi0      ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC 1$\\pi^0$";                          break;
     case ( kSKmNuMuNCcoh       ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC coherent";                           break;
     case ( kSKmNuMuNCoth_e     ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC other ($\\nu_e$ app)";               break;
     case ( kSKmNuMuNCoth_mu    ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC other ($\\nu_\\mu$ disap)";          break;
     case ( kSKmOscNuMuCCQE     ) : return "Osc. $\\nu_\\mu$ CCQE";                                          break;
     case ( kSKmOscNuMuCC1pi    ) : return "Osc. $\\nu_\\mu$ CC 1$\\pi$";                                    break;
     case ( kSKmOscNuMuCCcoh    ) : return "Osc. $\\nu_\\mu$ CC coherent";                                   break; 
     case ( kSKmOscNuMuCCoth    ) : return "Osc. $\\nu_\\mu$ CC other";                                      break;
     case ( kSKmNuMuBarCCQE     ) : return "$\\bar{\\nu}_\\mu$ CCQE";                                        break;
     case ( kSKmNuMuBarCC1pi    ) : return "$\\bar{\\nu}_\\mu$ CC 1$\\pi$";                                  break;
     case ( kSKmNuMuBarCCcoh    ) : return "$\\bar{\\nu}_\\mu$ CC coherent";                                 break;
     case ( kSKmNuMuBarCCoth    ) : return "$\\bar{\\nu}_\\mu$ CC other";                                    break;
     case ( kSKmNuMuBarNC1piC   ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC 1$\\pi^{+/-}$";        break;
     case ( kSKmNuMuBarNC1pi0   ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC 1$\\pi^0$";            break;
     case ( kSKmNuMuBarNCcoh    ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC coherent";             break;
     case ( kSKmNuMuBarNCoth_e  ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC other ($\\nu_e$ app)"; break;
     case ( kSKmNuMuBarNCoth_mu ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC other ($\\nu_\\mu$ disap)"; break;
     case ( kSKmOscNuMuBarCCQE  ) : return "Osc. $\\bar{\\nu}_\\mu$ CCQE";                                   break;
     case ( kSKmOscNuMuBarCC1pi ) : return "Osc. $\\bar{\\nu}_\\mu$ CC 1$\\pi$";                             break;
     case ( kSKmOscNuMuBarCCcoh ) : return "Osc. $\\bar{\\nu}_\\mu$ CC coherent";                            break; 
     case ( kSKmOscNuMuBarCCoth ) : return "Osc. $\\bar{\\nu}_\\mu$ CC other";                               break;
     case ( kSKmNuECCQE         ) : return "$\\nu_e$ CCQE";                                                  break;
     case ( kSKmNuECC1pi        ) : return "$\\nu_e$ CC 1$\\pi$";                                            break;
     case ( kSKmNuECCcoh        ) : return "$\\nu_e$ CC coherent";                                           break;
     case ( kSKmNuECCoth        ) : return "$\\nu_e$ CC other";                                              break;
     case ( kSKmNuENC1piC       ) : return "$\\nu_e$ NC 1$\\pi^{+/-}$";                                      break;
     case ( kSKmNuENC1pi0       ) : return "$\\nu_e$ NC 1$\\pi^0$";                                          break;
     case ( kSKmNuENCcoh        ) : return "$\\nu_e$ NC coherent";                                           break;
     case ( kSKmNuENCoth_e      ) : return "$\\nu_e$ NC other ($\\nu_e$ app)";                               break;
     case ( kSKmNuENCoth_mu     ) : return "$\\nu_e$ NC other ($\\nu_\\mu$ disap)";                          break;
     case ( kSKmOscNuECCQE      ) : return "Osc. $\\nu_e$ CCQE";                                             break;
     case ( kSKmOscNuECC1pi     ) : return "Osc. $\\nu_e$ CC 1$\\pi$";                                       break;
     case ( kSKmOscNuECCcoh     ) : return "Osc. $\\nu_e$ CC coherent";                                      break;
     case ( kSKmOscNuECCoth     ) : return "Osc. $\\nu_e$ CC other";                                         break;
     case ( kSKmNuEBarCCQE      ) : return "$\\bar{\\nu}_e$ CCQE";                                           break;
     case ( kSKmNuEBarCC1pi     ) : return "$\\bar{\\nu}_e$ CC 1$\\pi$";                                     break;
     case ( kSKmNuEBarCCcoh     ) : return "$\\bar{\\nu}_e$ CC coherent";                                    break;
     case ( kSKmNuEBarCCoth     ) : return "$\\bar{\\nu}_e$ CC other";                                       break;
     case ( kSKmNuEBarNC1piC    ) : return "$\\bar{\\nu}_e$ NC 1$\\pi^{+/-}$";                               break;
     case ( kSKmNuEBarNC1pi0    ) : return "$\\bar{\\nu}_e$ NC 1$\\pi^0$";                                   break;
     case ( kSKmNuEBarNCcoh     ) : return "$\\bar{\\nu}_e$ NC coherent";                                    break;
     case ( kSKmNuEBarNCoth_e   ) : return "$\\bar{\\nu}_e$ NC other ($\\nu_e$ app)";                        break;
     case ( kSKmNuEBarNCoth_mu  ) : return "$\\bar{\\nu}_e$ NC other ($\\nu_\\mu$ disap)";                   break;
     case ( kSKmOscNuEBarCCQE   ) : return "Osc. $\\bar{\\nu}_e$ CCQE";                                      break;
     case ( kSKmOscNuEBarCC1pi  ) : return "Osc. $\\bar{\\nu}_e$ CC 1$\\pi$";                                break;
     case ( kSKmOscNuEBarCCcoh  ) : return "Osc. $\\bar{\\nu}_e$ CC coherent";                               break; 
     case ( kSKmOscNuEBarCCoth  ) : return "Osc. $\\bar{\\nu}_e$ CC other";                                  break;

    default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact (t2kana::conv::SuperkMode_t m)
{
// Compact string used for naming output histograms stored ROOT files
//
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numuccqe";          break;
     case ( kSKmMEC             ) : return "mec";               break;
     case ( kSKmNuMuCC1pi       ) : return "numucc1pi";         break;
     case ( kSKmNuMuCCcoh       ) : return "numucccoh";         break;
     case ( kSKmNuMuCCoth       ) : return "numuccoth";         break;
     case ( kSKmNuMuNC1piC      ) : return "numunc1piC";        break;
     case ( kSKmNuMuNC1pi0      ) : return "numunc1pi0";        break;
     case ( kSKmNuMuNCcoh       ) : return "numunccoh";         break;
     case ( kSKmNuMuNCoth_e     ) : return "numuncoth";         break;
     case ( kSKmNuMuNCoth_mu    ) : return "numuncoth";         break;
     case ( kSKmOscNuMuCCQE     ) : return "oscnumuccqe";       break;
     case ( kSKmOscNuMuCC1pi    ) : return "oscnumucc1pi";      break;
     case ( kSKmOscNuMuCCcoh    ) : return "oscnumucccoh";      break;
     case ( kSKmOscNuMuCCoth    ) : return "oscnumuccoth";      break;
     case ( kSKmNuMuBarCCQE     ) : return "numubarccqe";       break;
     case ( kSKmNuMuBarCC1pi    ) : return "numubarcc1pi";      break;
     case ( kSKmNuMuBarCCcoh    ) : return "numubarcccoh";      break;
     case ( kSKmNuMuBarCCoth    ) : return "numubarccoth";      break;
     case ( kSKmNuMuBarNC1piC   ) : return "numubarnc1piC";     break;
     case ( kSKmNuMuBarNC1pi0   ) : return "numubarnc1pi0";     break;
     case ( kSKmNuMuBarNCcoh    ) : return "numubarnccoh";      break;
     case ( kSKmNuMuBarNCoth_e  ) : return "numubarncoth";      break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numubarncoth";      break;
     case ( kSKmOscNuMuBarCCQE  ) : return "oscnumubarccqe";    break;
     case ( kSKmOscNuMuBarCC1pi ) : return "oscnumubarcc1pi";   break;
     case ( kSKmOscNuMuBarCCcoh ) : return "oscnumubarcccoh";   break;
     case ( kSKmOscNuMuBarCCoth ) : return "oscnumubarccoth";   break;
     case ( kSKmNuECCQE         ) : return "nueccqe";           break;
     case ( kSKmNuECC1pi        ) : return "nuecc1pi";          break;
     case ( kSKmNuECCcoh        ) : return "nuecccoh";          break;
     case ( kSKmNuECCoth        ) : return "nueccoth";          break;
     case ( kSKmNuENC1piC       ) : return "nuenc1piC";         break;
     case ( kSKmNuENC1pi0       ) : return "nuenc1pi0";         break;
     case ( kSKmNuENCcoh        ) : return "nuenccoh";          break;
     case ( kSKmNuENCoth_e      ) : return "nuencoth";          break;
     case ( kSKmNuENCoth_mu     ) : return "nuencoth";          break;
     case ( kSKmOscNuECCQE      ) : return "oscnueccqe";        break;
     case ( kSKmOscNuECC1pi     ) : return "oscnuecc1pi";       break;
     case ( kSKmOscNuECCcoh     ) : return "oscnuecccoh";       break;
     case ( kSKmOscNuECCoth     ) : return "oscnueccoth";       break;
     case ( kSKmNuEBarCCQE      ) : return "nuebarccqe";        break;
     case ( kSKmNuEBarCC1pi     ) : return "nuebarcc1pi";       break;
     case ( kSKmNuEBarCCcoh     ) : return "nuebarcccoh";       break;
     case ( kSKmNuEBarCCoth     ) : return "nuebarccoth";       break;
     case ( kSKmNuEBarNC1piC    ) : return "nuebarnc1piC";      break;
     case ( kSKmNuEBarNC1pi0    ) : return "nuebarnc1pi0";      break;
     case ( kSKmNuEBarNCcoh     ) : return "nuebarnccoh";       break;
     case ( kSKmNuEBarNCoth_e   ) : return "nuebarncoth";       break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nuebarncoth";       break;
     case ( kSKmOscNuEBarCCQE   ) : return "oscnuebarccqe";     break;
     case ( kSKmOscNuEBarCC1pi  ) : return "oscnuebarcc1pi";    break;
     case ( kSKmOscNuEBarCCcoh  ) : return "oscnuebarcccoh";    break;
     case ( kSKmOscNuEBarCCoth  ) : return "oscnuebarccoth";    break;

     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
t2kana::conv::SuperkMode_t t2kana::util::SuperkModeFromString(string sh)
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
//____________________DEPRECATED______________________________________________
bool t2kana::util::IsMEC(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmMEC      ) :
        return true;
        break;
     default: 
        return false;
        break;
  }
  return false;
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_splinehack (t2kana::conv::SuperkMode_t m)
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
