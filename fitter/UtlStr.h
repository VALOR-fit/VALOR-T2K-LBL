//____________________________________________________________________________
/*!

\class     t2kana::utils

\brief     Utilities for string manipulation

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_str_utils_h_
#define _t2kana_str_utils_h_

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace t2kana {
namespace util {

  string         TrimSpaces (string input);	
  vector<string> Tokens     (string input, string delimiter);
  
} // util
} // t2kana
  
#endif


