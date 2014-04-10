//____________________________________________________________________________
/*!

\class    t2kana::conv::NuOscMixParam_t

\brief    Choices of oscillation mixing parameter

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

*/
//____________________________________________________________________________

#ifndef nuosc_mixing_parameter_h
#define nuosc_mixing_parameter_h

#include <string>

using std::string;

namespace t2kana {

namespace conv {

typedef enum ENuOscMixParam
{   
   kNuMixUndefined = -1, 
   kNuMixSinsq2Theta,      // sin^2(2*theta)
   kNuMixSinsqTheta        // sin^2(theta)
} 
NuOscMixParam_t;
   
} // conv

namespace util {

   string AsString (conv::NuOscMixParam_t o);

   conv::NuOscMixParam_t MixingTypeFromString (string h);

} // util

} // t2kana
  
#endif 

