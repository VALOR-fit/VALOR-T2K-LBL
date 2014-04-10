//____________________________________________________________________________
/*!

\enum      t2kana::conv::NuMassHierarchy_t

\brief     Neutrino mass hierarchy enumeration

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef nu_mass_hierarchy_h
#define nu_mass_hierarchy_h

#include <string>

using std::string;

namespace t2kana {

namespace conv {

typedef enum ENuMassHierarchy 
{   
   kNmhUndefined = -1, 
   kNmhNormal,
   kNmhInverted
} 
NuMassHierarchy_t;
   
} // conv

namespace util {

  string AsString (conv::NuMassHierarchy_t  m);

  conv::NuMassHierarchy_t HierarchyFromString (string h);

} // util

} // t2kana
  
#endif 

