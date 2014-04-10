//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include "OscChannelId.h" 
#include "PDG.h"

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
string util::AsString(OscChannelId_t osch)
{
  switch(osch) {
     case ( kOscNuE2NuE           ) : return "nue -> nue";            break;
     case ( kOscNuE2NuMu          ) : return "nue -> numu";           break;
     case ( kOscNuE2NuTau         ) : return "nue -> nutau";          break;
     case ( kOscNuMu2NuE          ) : return "numu -> nue";           break;
     case ( kOscNuMu2NuMu         ) : return "numu -> numu";          break;
     case ( kOscNuMu2NuTau        ) : return "numu -> nutau";         break;
     case ( kOscNuTau2NuE         ) : return "nutau -> nue";          break;
     case ( kOscNuTau2NuMu        ) : return "nutau -> numu";         break;
     case ( kOscNuTau2NuTau       ) : return "nutau -> nutau";        break;
     case ( kOscNuEBar2NuEBar     ) : return "nuebar -> nuebar";      break;
     case ( kOscNuEBar2NuMuBar    ) : return "nuebar -> numubar";     break;
     case ( kOscNuEBar2NuTauBar   ) : return "nuebar -> nutaubar";    break;
     case ( kOscNuMuBar2NuEBar    ) : return "numubar -> nuebar";     break;
     case ( kOscNuMuBar2NuMuBar   ) : return "numubar -> numubar";    break;
     case ( kOscNuMuBar2NuTauBar  ) : return "numubar -> nutaubar";   break;
     case ( kOscNuTauBar2NuEBar   ) : return "nutaubar -> nuebar";    break;
     case ( kOscNuTauBar2NuMuBar  ) : return "nutaubar -> numubar";   break;
     case ( kOscNuTauBar2NuTauBar ) : return "nutaubar -> nutaubar";  break;
     
     case ( kOscNuE2NuS4          ) : return "nue -> nus4";
     case ( kOscNuE2NuS5          ) : return "nue -> nus5";
     case ( kOscNuMu2NuS4         ) : return "numu -> nus4";
     case ( kOscNuMu2NuS5         ) : return "numu -> nus5";
     case ( kOscNuTau2NuS4        ) : return "nutau -> nus4";
     case ( kOscNuTau2NuS5        ) : return "nutau -> nus5";
     
     case ( kOscNuS42NuE          ) : return "nus4 -> nue";
     case ( kOscNuS52NuE          ) : return "nus5 -> nue";
     case ( kOscNuS42NuMu         ) : return "nus4 -> numu";
     case ( kOscNuS52NuMu         ) : return "nus5 -> numu";
     case ( kOscNuS42NuTau        ) : return "nus4 -> nutau";
     case ( kOscNuS52NuTau        ) : return "nus5 -> nutau";
     
     case ( kOscNuEBar2NuS4Bar    ) : return "nuebar -> nus4bar";
     case ( kOscNuEBar2NuS5Bar    ) : return "nuebar -> nus5bar";
     case ( kOscNuMuBar2NuS4Bar   ) : return "numubar -> nus4bar";
     case ( kOscNuMuBar2NuS5Bar   ) : return "numubar -> nus5bar";
     case ( kOscNuTauBar2NuS4Bar  ) : return "nutaubar -> nus4bar";
     case ( kOscNuTauBar2NuS5Bar  ) : return "nutaubar -> nus5bar";
     
     case ( kOscNuS4Bar2NuEBar    ) : return "nus4bar -> nuebar";
     case ( kOscNuS5Bar2NuEBar    ) : return "nus5bar -> nuebar";
     case ( kOscNuS4Bar2NuMuBar   ) : return "nus4bar -> numubar";
     case ( kOscNuS5Bar2NuMuBar   ) : return "nus5bar -> numubar";
     case ( kOscNuS4Bar2NuTauBar  ) : return "nus4bar -> nutaubar";
     case ( kOscNuS5Bar2NuTauBar  ) : return "nus5bar -> nutaubar";
     
     default                        : return "--";
  }
  return "--";
}
//________________________________________________________________________________________
OscChannelId_t util::GetOscChannel(int nu_start, int nu_end)
{
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuE      ) return kOscNuE2NuE;
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuMu     ) return kOscNuE2NuMu;
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuTau    ) return kOscNuE2NuTau;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuE      ) return kOscNuMu2NuE;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuMu     ) return kOscNuMu2NuMu;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuTau    ) return kOscNuMu2NuTau;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuE      ) return kOscNuTau2NuE;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuMu     ) return kOscNuTau2NuMu;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuTau    ) return kOscNuTau2NuTau;

  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuEBar   ) return kOscNuEBar2NuEBar;
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuMuBar  ) return kOscNuEBar2NuMuBar;
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuTauBar ) return kOscNuEBar2NuTauBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuEBar   ) return kOscNuMuBar2NuEBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuMuBar  ) return kOscNuMuBar2NuMuBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuTauBar ) return kOscNuMuBar2NuTauBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuEBar   ) return kOscNuTauBar2NuEBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuMuBar  ) return kOscNuTauBar2NuMuBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuTauBar ) return kOscNuTauBar2NuTauBar;

  //oscillations to or from sterile neutrinos NuS4 or NuS5
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuS4     ) return kOscNuE2NuS4;
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuS5     ) return kOscNuE2NuS5;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuS4     ) return kOscNuMu2NuS4;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuS5     ) return kOscNuMu2NuS5;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuS4     ) return kOscNuTau2NuS4;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuS5     ) return kOscNuTau2NuS5;

  if ( nu_start == kPdgNuS4     && nu_end == kPdgNuE      ) return kOscNuS42NuE;
  if ( nu_start == kPdgNuS5     && nu_end == kPdgNuE      ) return kOscNuS52NuE;
  if ( nu_start == kPdgNuS4     && nu_end == kPdgNuMu     ) return kOscNuS42NuMu;
  if ( nu_start == kPdgNuS5     && nu_end == kPdgNuMu     ) return kOscNuS52NuMu;
  if ( nu_start == kPdgNuS4     && nu_end == kPdgNuTau    ) return kOscNuS42NuTau;
  if ( nu_start == kPdgNuS5     && nu_end == kPdgNuTau    ) return kOscNuS52NuTau;

  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuS4Bar  ) return kOscNuEBar2NuS4Bar;
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuS5Bar  ) return kOscNuEBar2NuS5Bar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuS4Bar  ) return kOscNuMuBar2NuS4Bar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuS5Bar  ) return kOscNuMuBar2NuS5Bar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuS4Bar  ) return kOscNuTauBar2NuS4Bar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuS5Bar  ) return kOscNuTauBar2NuS5Bar;

  if ( nu_start == kPdgNuS4Bar  && nu_end == kPdgNuEBar   ) return kOscNuS4Bar2NuEBar;
  if ( nu_start == kPdgNuS5Bar  && nu_end == kPdgNuEBar   ) return kOscNuS5Bar2NuEBar;
  if ( nu_start == kPdgNuS4Bar  && nu_end == kPdgNuMuBar  ) return kOscNuS4Bar2NuMuBar;
  if ( nu_start == kPdgNuS5Bar  && nu_end == kPdgNuMuBar  ) return kOscNuS5Bar2NuMuBar;
  if ( nu_start == kPdgNuS4Bar  && nu_end == kPdgNuTauBar ) return kOscNuS4Bar2NuTauBar;
  if ( nu_start == kPdgNuS5Bar  && nu_end == kPdgNuTauBar ) return kOscNuS5Bar2NuTauBar;

  return kOscUndefined;
}
//________________________________________________________________________________________

