//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include "NuOscOctant.h"

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
string t2kana::util::AsString(NuOscOctant_t o)
{
  switch(o) {
   case (kNuOctantOne) :
     return "sin^2(theta_{23}) <= 0.5";
     break;
   case (kNuOctantTwo) :
     return "sin^2(theta_{23}) >= 0.5";
     break;
   default :
     return "Undefined";
     break;
  }
  return "";
}
//________________________________________________________________________________________
NuOscOctant_t t2kana::util::OctantFromString(string h)
{
  if(h == "one") return kNuOctantOne;
  if(h == "two") return kNuOctantTwo;

  return kNuOctantUndefined;
}
//________________________________________________________________________________________

