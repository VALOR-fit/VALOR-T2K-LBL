//____________________________________________________________________________
/*

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________

#include "NuOscMixParam.h"

using namespace t2kana;
using namespace t2kana::conv;

//____________________________________________________________________________
string t2kana::util::AsString(conv::NuOscMixParam_t p)
{
  switch(p) {
   case (kNuMixSinsq2Theta) :
     return "sin^2(2*theta_{ij})";
     break;
   case (kNuMixSinsqTheta) :
     return "sin^2(theta_{ij})";
     break;
   default :
     return "Undefined";
     break;
  }
  return "";
}
//____________________________________________________________________________
NuOscMixParam_t t2kana::util::MixingTypeFromString(string h)
{
  if(h == "sinsq2theta") return kNuMixSinsq2Theta;
  if(h == "sinsqtheta" ) return kNuMixSinsqTheta;

  return kNuMixUndefined;
}
//____________________________________________________________________________

