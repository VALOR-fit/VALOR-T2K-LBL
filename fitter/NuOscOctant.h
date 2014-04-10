//____________________________________________________________________________
/*!

\class     t2kana::conv::NuOscOctant_t

\brief     Choice of theta_{23} octant

\author    Lorena Escudero <loressa@ific.uv.es>
           IFIC Valencia
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef nuosc_octant_h
#define nuosc_octant_h

#include <string>

using std::string;

namespace t2kana {

namespace conv {

typedef enum ENuOscOctant
{   
   kNuOctantUndefined = -1, 
   kNuOctantOne,      // sin^2(theta_{23}) <= 0.5
   kNuOctantTwo       // sin^2(theta_{23}) >= 0.5
} 
NuOscOctant_t;

} // conv

namespace util {

  string AsString(conv::NuOscOctant_t p);

  conv::NuOscOctant_t OctantFromString(string h);

} // util

} // t2kana
  
#endif 

