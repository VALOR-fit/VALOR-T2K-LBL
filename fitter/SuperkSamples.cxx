//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/

#include "Enumerations.h"
#include "SuperkSamples.h"
#include "Utilities.h"
#include "Binning.h"

using namespace t2kana;
using namespace t2kana::conv;

//____________________________________________________________________________
t2kana::conv::SuperkSample_t t2kana::util::SKSampleFromString(string sh)
{
  if(sh == "1Rmu") return kSKsFCFV1RingMuLike;
  if(sh == "1Re" ) return kSKsFCFV1RingELike;
  if(sh == "RHC1Rmu") return kSKsRHCFCFV1RingMuLike;
  if(sh == "RHC1Re" ) return kSKsRHCFCFV1RingELike;

  return kSKsUndefined;
}
//________________________________________________________________________________________
bool t2kana::util::IsFHCSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    return false;
  return true;
}
//___________________________________________________________________________
bool t2kana::util::IsRHCSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    return true;
  return false;
}
//___________________________________________________________________________
bool t2kana::util::Is1RMuSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsFCFV1RingMuLike)
    return true;
  return false;
}
//___________________________________________________________________________
bool t2kana::util::Is1RESample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingELike || s == kSKsFCFV1RingELike)
    return true;
  return false;
}
//___________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return "FCFV 1-ring mu-like"; break;
     case ( kSKsFCFV1RingELike  ) : return "FCFV 1-ring e-like";  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return "FCFV 1-ring mu-like (RHC running)"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "FCFV 1-ring e-like (RHC running)";  break;

     default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact (t2kana::conv::SuperkSample_t s)
{
// Compact string used for naming output histograms stored ROOT files
//
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return "1Rmu";    break;
     case ( kSKsFCFV1RingELike     ) : return "1Re";     break;
     case ( kSKsRHCFCFV1RingMuLike ) : return "RHC1Rmu"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "RHC1Re";  break;

     default                         : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringTex(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return "1 $\\mu$-like ring"; break;
     case ( kSKsFCFV1RingELike  ) : return "1 e-like ring";      break;
     case ( kSKsRHCFCFV1RingMuLike ) : return "1 $\\mu$-like ring (RHC running)"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "1 e-like ring (RHC running)";      break;

     default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringRooTex(t2kana::conv::SuperkSample_t samp)
{
  string s = AsStringTex(samp);
  replace( s.begin(), s.end(), '\\', '#' ); // Make slashes hashes for RooTex
  s.erase( std::remove( s.begin(), s.end(), '$'), s.end() ) ; // Remove dollar signs.
  return s;
}
//________________________________________________________________________________________
