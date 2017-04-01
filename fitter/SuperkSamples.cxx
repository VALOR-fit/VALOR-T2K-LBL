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
  if(sh == "1Rmu")     return kSKsFCFV1RingMuLike;
  if(sh == "1Re" )     return kSKsFCFV1RingELike;
  if(sh == "RHC1Rmu")  return kSKsRHCFCFV1RingMuLike;
  if(sh == "RHC1Re" )  return kSKsRHCFCFV1RingELike;
  if(sh == "MultiRmu") return kSKsFCFVMultiRingMuLike;     
  if(sh == "MultiRe" ) return kSKsFCFVMultiRingELike;      

  return kSKsUndefined;
}
//________________________________________________________________________________________
bool valor::util::IsFHCSample(conv::SuperkSample_t s)
{
  if(s == kSKsFCFV1RingMuLike || s == kSKsFCFV1RingELike ||
     s == kSKsFCFVMultiRingMuLike || s == kSKsFCFVMultiRingELike)
    return true;
  return false;
}
//___________________________________________________________________________
bool valor::util::IsRHCSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    return true;
  return false;
}
//___________________________________________________________________________
bool valor::util::IsMultiRingSample(conv::SuperkSample_t s)
{
  if(s == kSKsFCFVMultiRingMuLike || s == kSKsFCFVMultiRingELike)
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
bool valor::util::IsMultiRMuSample(conv::SuperkSample_t s)
{
  if(s == kSKsFCFVMultiRingMuLike)
    return true;
  return false;
}
//___________________________________________________________________________
bool valor::util::IsMultiRESample(conv::SuperkSample_t s)
{
  if(s == kSKsFCFVMultiRingELike)
    return true;
  return false;
}
//___________________________________________________________________________
string valor::util::GetXMLMCTreeName(valor::conv::SuperkSample_t s)
{
  switch(s) {
    case ( kSKsFCFV1RingMuLike     ) : return "one_mu_ring";
    case ( kSKsFCFV1RingELike      ) : return "one_e_ring";
      
    case ( kSKsRHCFCFV1RingMuLike  ) : return "one_mu_ring_rhc";  
    case ( kSKsRHCFCFV1RingELike   ) : return "one_e_ring_rhc" ;  
    //~ case ( kSKsRHCFCFV1RingMuLike  ) : return "one_mu_ring"; // For when using old MC
    //~ case ( kSKsRHCFCFV1RingELike   ) : return "one_e_ring";
      
      /// NEW MULTIRING
    case ( kSKsFCFVMultiRingMuLike ) : return "multi_mu_ring";    // multiring not available yet
    case ( kSKsFCFVMultiRingELike  ) : return "multi_e_ring";     // multiring not available yet
      //case ( kSKsFCFVMultiRingMuLike ) : return "one_mu_ring";    
      //case ( kSKsFCFVMultiRingELike  ) : return "one_e_ring";
      ////
      
    default                          : return "--";
  }
  return "";
}
//___________________________________________________________________________
string valor::util::GetXMLDatasetName(valor::conv::SuperkSample_t s)
{
  switch(s) {
     
     case ( kSKsFCFV1RingMuLike     ) : return "one_mu_ring";
     case ( kSKsFCFV1RingELike      ) : return "one_e_ring";
     
  case ( kSKsRHCFCFV1RingMuLike  ) : return "one_mu_ring_rhc";  // rhc available
  case ( kSKsRHCFCFV1RingELike   ) : return "one_e_ring_rhc"; // rhc available
    
    /// NEW MULTIRING
  case ( kSKsFCFVMultiRingMuLike ) : return "multi_mu_ring";  // multiring not available yet
  case ( kSKsFCFVMultiRingELike  ) : return "multi_e_ring";   // multiring not available yet
    //case ( kSKsFCFVMultiRingMuLike ) : return "one_mu_ring";
    //case ( kSKsFCFVMultiRingELike  ) : return "one_e_ring";
    ///
     
     default                          : return "--";
  }
  return "";
}
//___________________________________________________________________________
string valor::util::AsString(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return "FCFV 1-ring mu-like"; break;
     case ( kSKsFCFV1RingELike      ) : return "FCFV 1-ring e-like";  break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return "FCFV 1-ring mu-like (RHC running)"; break;
     case ( kSKsRHCFCFV1RingELike   ) : return "FCFV 1-ring e-like (RHC running)";  break;
     case ( kSKsFCFVMultiRingMuLike ) : return "FCFV multi-ring mu-like"; break;          /// NEW MULTIRING
     case ( kSKsFCFVMultiRingELike  ) : return "FCFV multi-ring e-like";  break;       /// NEW MULTIRING
     default                          : return "--";
  }
  return "";
}
//____________________________________________________________________________
string valor::util::AsStringCompact (valor::conv::SuperkSample_t s)
{
// Compact string used for naming output histograms stored ROOT files
//
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return "1Rmu";    break;
     case ( kSKsFCFV1RingELike      ) : return "1Re";     break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return "RHC1Rmu"; break;
     case ( kSKsRHCFCFV1RingELike   ) : return "RHC1Re";  break;
     case ( kSKsFCFVMultiRingMuLike ) : return "MultiRmu";    break;         /// NEW MULTIRING
     case ( kSKsFCFVMultiRingELike  ) : return "MultiRe";     break;       /// NEW MULTIRING
     default                          : return "";
  }
  return "";
}
//____________________________________________________________________________
string valor::util::AsStringTex(valor::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike     ) : return "1 $\\mu$-like ring"; break;
     case ( kSKsFCFV1RingELike      ) : return "1 e-like ring";      break;
     case ( kSKsRHCFCFV1RingMuLike  ) : return "1 $\\mu$-like ring (RHC running)"; break;
     case ( kSKsRHCFCFV1RingELike   ) : return "1 e-like ring (RHC running)";      break;
     case ( kSKsFCFVMultiRingMuLike ) : return "Multi $\\mu$-like ring"; break;       /// NEW MULTIRING
     case ( kSKsFCFVMultiRingELike  ) : return "Multi e-like ring";      break;       /// NEW MULTIRING
     default                          : return "--";
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
conv::SuperkSample_t valor::util::GetSKSample (int s)
{
  if (s < conv::kNSuperkSamples) {
    return conv::kSKSamples[s];
  }
  else return conv::kSKsUndefined;
}
//________________________________________________________________________________________
