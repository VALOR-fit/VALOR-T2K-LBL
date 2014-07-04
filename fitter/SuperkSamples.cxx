//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/

#include "Enumerations.h"
#include "SuperkSamples.h"
#include "Utilities.h"
#include "Binning.h"

using namespace valor;
using namespace valor::conv;

//____________________________________________________________________________
valor::conv::SuperkSample_t valor::util::SKSampleFromString(string sh)
{
  if(sh == "1Rmu") return kSKsFCFV1RingMuLike;
  if(sh == "1Re" ) return kSKsFCFV1RingELike;
  if(sh == "RHC1Rmu") return kSKsRHCFCFV1RingMuLike;
  if(sh == "RHC1Re" ) return kSKsRHCFCFV1RingELike;

  return kSKsUndefined;
}
//________________________________________________________________________________________
bool valor::util::IsFHCSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    return false;
  return true;
}
//___________________________________________________________________________
bool valor::util::IsRHCSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    return true;
  return false;
}
//___________________________________________________________________________
bool valor::util::Is1RMuSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsFCFV1RingMuLike)
    return true;
  return false;
}
//___________________________________________________________________________
bool valor::util::Is1RESample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingELike || s == kSKsFCFV1RingELike)
    return true;
  return false;
}
//___________________________________________________________________________
string valor::util::GetXMLMCTreeName(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike )    : return "one_mu_ring";
     case ( kSKsFCFV1RingELike  )    : return "one_e_ring";
     case ( kSKsRHCFCFV1RingMuLike ) : return "one_mu_ring";
     case ( kSKsRHCFCFV1RingELike  ) : return "one_e_ring";

     default                      : return "--";
  }
  return "";
}
//___________________________________________________________________________
string valor::util::GetXMLDatasetName(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike )    : return "one_mu_ring";
     case ( kSKsFCFV1RingELike  )    : return "one_e_ring";
     case ( kSKsRHCFCFV1RingMuLike ) : return "one_mu_ring_rhc";
     case ( kSKsRHCFCFV1RingELike  ) : return "one_e_ring_rhc";

     default                      : return "--";
  }
  return "";
}
//___________________________________________________________________________
string valor::util::AsString(valor::conv::SuperkSample_t s)
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
string valor::util::AsStringCompact (valor::conv::SuperkSample_t s)
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
string valor::util::AsStringTex(valor::conv::SuperkSample_t s)
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
string valor::util::AsStringRooTex(valor::conv::SuperkSample_t samp)
{
  string s = AsStringTex(samp);
  replace( s.begin(), s.end(), '\\', '#' ); // Make slashes hashes for RooTex
  s.erase( std::remove( s.begin(), s.end(), '$'), s.end() ) ; // Remove dollar signs.
  return s;
}
//________________________________________________________________________________________
string valor::util::SampleAsString       (conv::SuperkSample_t s) { return AsString       (s);}
string valor::util::SampleAsStringCompact(conv::SuperkSample_t s) { return AsStringCompact(s);}
string valor::util::SampleAsStringTex    (conv::SuperkSample_t s) { return AsStringTex    (s);}
string valor::util::SampleAsStringRooTex (conv::SuperkSample_t s) { return AsStringRooTex (s);}
//________________________________________________________________________________________
