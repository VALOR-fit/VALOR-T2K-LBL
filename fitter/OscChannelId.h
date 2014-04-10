//____________________________________________________________________________
/*!

\class     t2kana::conv::OscChannelId_t

\brief     Enumeration of oscillation channels

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef osc_channel_id_h
#define osc_channel_id_h

#include <string>

using std::string;

namespace t2kana {

namespace conv {

typedef enum EOscChannelId 
{   
   kOscUndefined = -1, 
   kOscNuE2NuE,
   kOscNuE2NuMu,
   kOscNuE2NuTau,
   kOscNuMu2NuE,
   kOscNuMu2NuMu,
   kOscNuMu2NuTau,
   kOscNuTau2NuE,
   kOscNuTau2NuMu,
   kOscNuTau2NuTau,
   kOscNuEBar2NuEBar,
   kOscNuEBar2NuMuBar,
   kOscNuEBar2NuTauBar,
   kOscNuMuBar2NuEBar,
   kOscNuMuBar2NuMuBar,
   kOscNuMuBar2NuTauBar,
   kOscNuTauBar2NuEBar,
   kOscNuTauBar2NuMuBar,
   kOscNuTauBar2NuTauBar,
   //oscillations to or from sterile neutrinos NuS4 or NuS5
   kOscNuE2NuS4,
   kOscNuE2NuS5,
   kOscNuMu2NuS4,
   kOscNuMu2NuS5,
   kOscNuTau2NuS4,
   kOscNuTau2NuS5,
   kOscNuS42NuE,
   kOscNuS52NuE,
   kOscNuS42NuMu,
   kOscNuS52NuMu,
   kOscNuS42NuTau,
   kOscNuS52NuTau,
   kOscNuEBar2NuS4Bar,
   kOscNuEBar2NuS5Bar,
   kOscNuMuBar2NuS4Bar,
   kOscNuMuBar2NuS5Bar,
   kOscNuTauBar2NuS4Bar,
   kOscNuTauBar2NuS5Bar,
   kOscNuS4Bar2NuEBar,
   kOscNuS5Bar2NuEBar,
   kOscNuS4Bar2NuMuBar,
   kOscNuS5Bar2NuMuBar,
   kOscNuS4Bar2NuTauBar,
   kOscNuS5Bar2NuTauBar
} 
OscChannelId_t;
   
} // conv

namespace util {

  string AsString (conv::OscChannelId_t o);

  conv::OscChannelId_t GetOscChannel(int nu_start, int nu_end);

} // util

} // t2kana
  
#endif 

