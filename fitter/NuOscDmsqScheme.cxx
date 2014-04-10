//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include "NuOscDmsqScheme.h"

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
string util::AsString(conv::NuOscDmsqScheme_t dmsqs)
{
  switch(dmsqs) {
  case(kNuDmsq12andAbs23) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |Dm^2_{32}|";
    break;
  case(kNuDmsq12andAbs31) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |Dm^2_{31}|";
    break;
  case(kNuDmsqFogliLisi ) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|";
    break;
  default :
    return "Undefined";
    break;
  }
}
//________________________________________________________________________________________
NuOscDmsqScheme_t util::DmsqSchemeFromString(string h)
{
  if(h == "abs23"    ) return kNuDmsq12andAbs23;
  if(h == "abs31"    ) return kNuDmsq12andAbs31;
  if(h == "fogli"    ) return kNuDmsqFogliLisi;

  return kNuDmsqUndefined;
}
//________________________________________________________________________________________

