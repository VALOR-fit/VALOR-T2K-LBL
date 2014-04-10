//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include "NuMassHierarchy.h"

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
string t2kana::util::AsString(NuMassHierarchy_t m)
{
  switch(m) {
   case (kNmhNormal) :
     return "Normal";
     break;
   case (kNmhInverted) :
     return "Inverted";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//________________________________________________________________________________________
NuMassHierarchy_t t2kana::util::HierarchyFromString(string h)
{
  if(h == "normal"   ) return kNmhNormal;
  if(h == "inverted" ) return kNmhInverted;

  return kNmhUndefined;
}
//________________________________________________________________________________________

