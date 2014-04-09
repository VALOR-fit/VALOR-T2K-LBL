//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include "PDG.h"

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
bool util::IsNeutrino(int pdg_code)
{
  if ( pdg_code == kPdgNuE    ||
       pdg_code == kPdgNuMu   ||
       pdg_code == kPdgNuTau  ||
       pdg_code == kPdgNuS4   ||
       pdg_code == kPdgNuS5   )
  {
    return true;
  }
  return false;
}
//________________________________________________________________________________________
bool util::IsAntiNeutrino(int pdg_code)
{
  if ( pdg_code == kPdgNuEBar    ||
       pdg_code == kPdgNuMuBar   ||
       pdg_code == kPdgNuTauBar  ||
       pdg_code == kPdgNuS4Bar   ||
       pdg_code == kPdgNuS5Bar   )
  {
    return true;
  }
  return false;
}
//________________________________________________________________________________________
string util::Pdg2String(int pdg_code)
{
  switch(pdg_code) {
     case ( kPdgNuE      ) : return "nue";      break;
     case ( kPdgNuMu     ) : return "numu";     break;
     case ( kPdgNuTau    ) : return "nutau";    break;
     case ( kPdgNuEBar   ) : return "nuebar";   break;
     case ( kPdgNuMuBar  ) : return "numubar";  break;
     case ( kPdgNuTauBar ) : return "nutaubar"; break;
     default               : return "";
  }
  return "";
}
//________________________________________________________________________________________
string util::Pdg2StringTex(int pdg_code)
{
  switch(pdg_code) {
     case ( kPdgNuE      ) : return "$\\nu_e$";      break;
     case ( kPdgNuMu     ) : return "$\\nu_\\mu$";     break;
     case ( kPdgNuTau    ) : return "$\\nu_\\tau$";    break;
     case ( kPdgNuEBar   ) : return "$\\bar{\nu_e}$";   break;
     case ( kPdgNuMuBar  ) : return "$\\bar{\\nu_\\mu}$";  break;
     case ( kPdgNuTauBar ) : return "$\\bar{\\nu_\\tau}$"; break;
     default               : return "";
  }
  return "";
}
//________________________________________________________________________________________



