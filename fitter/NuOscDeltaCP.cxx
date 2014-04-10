//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include "NuOscDeltaCP.h"

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
string t2kana::util::AsString(NuOscDeltaCP_t d)
{
  switch(d) {
  case(kNuDeltaCPPhase) :
    return "d_{CP}";
    break;
  case(kNuDeltaCPSin) :
    return "sin(d_{CP})";
    break;
  default :
    return "Undefined";
    break;
  }
}
//________________________________________________________________________________________
NuOscDeltaCP_t t2kana::util::CPVChoiceFromString (string d)
{
  if(d == "sindelta" ) return kNuDeltaCPSin;
  if(d == "delta"    ) return kNuDeltaCPPhase;

  return kNuDeltaCPUndefined;
}
//________________________________________________________________________________________

