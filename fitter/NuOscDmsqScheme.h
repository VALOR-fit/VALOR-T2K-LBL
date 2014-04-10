//____________________________________________________________________________
/*!

\class     t2kana::conv::NuOscDmsqScheme_t

\brief     Choice of squared mass splitting parameter

\author    Lorena Escudero <loressa@ific.uv.es>
           IFIC Valencia
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group
*/
//____________________________________________________________________________

#ifndef nuosc_dmsq_scheme_h
#define nuosc_dmsq_scheme_h

#include <string>

using std::string;

namespace t2kana {

namespace conv {

typedef enum ENuOscDmsqScheme
{
  kNuDmsqUndefined = -1,
  kNuDmsq12andAbs23,  // use |Dm^2_{23}| and Dm^2_{12} as the 2 independent squared mass splittings                                                                            
  kNuDmsqFogliLisi,   // scheme used by Fogli and Lisi, dm^2 = Dm^2_{12} and Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|                                                            
  kNuDmsq12andAbs31   // use |Dm^2_{31}| and Dm^2_{12} as the 2 independent squared mass splittings                                                                            
}
NuOscDmsqScheme_t;
   
} // conv

namespace util {

  string AsString(conv::NuOscDmsqScheme_t dmsqs);

  conv::NuOscDmsqScheme_t DmsqSchemeFromString (string h);

} // util

} // t2kana
  
#endif 

