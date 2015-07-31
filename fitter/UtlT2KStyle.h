//____________________________________________________________________________
/*!

\class     valor::utils

\brief     Utilities for formatting and style

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _valor_t2kstyle_utils_h_
#define _valor_t2kstyle_utils_h_

#include <iostream>
#include <string>
#include <vector>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TString.h>

#include "Enumerations.h"

class TGraph;
class TH1;
class TH2;

using std::string;
using std::vector;
using std::cout;

namespace valor {
namespace util {
  
  // T2K specific
 
   void SetDefaultSKSpectrumStyle (conv::SuperkSample_t s, conv::SuperkMode_t m, TH1D * histogram);
   void SetDefaultSKSpectrumStyle (conv::SuperkSample_t s, conv::SuperkMode_t m, TH2D * histogram); // NEW 2D TEMPLATES
   // When making an ereco hist, set directory/axis titles/line colour
   void FormatErecoHist(TH1 * hist);
   void FormatErecoVsQ2Hist(TH2 * hist); // NEW 2D TEMPLATES

} // util
} // valor
  
#endif


