//____________________________________________________________________________
/*!

\class     t2kana::utils

\brief     Utilities for formatting and style

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_style_utils_h_
#define _t2kana_style_utils_h_

#include <iostream>
#include <string>
#include <vector>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TString.h>

#include "Enumerations.h"

class TGraph;
class TH1;

using std::string;
using std::vector;
using std::cout;

namespace t2kana {
namespace util {

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr,
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst,
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);

  void PrintBanner(void);
  void PrintBanner(string filename, UInt_t wait_msec);

  //~ void PrintCompact    (const TVectorD & v, int w=15, int prec=9, Option_t * opt="", ostream & stream = cout);
  //~ void PrintCompact    (const TVectorD & current, const TVectorD & nominal, int w, int prec, Option_t * option, ostream & stream = cout);
  //~ void PrintCompact    (const TMatrixD & m, int w=15, int prec=9, ostream & stream = cout);
  //~ void PrintCompactTex (const TMatrixD & m, int w=20, int prec=9, ostream & stream = cout, TString opt = "");

  // T2K specific
 
   void SetDefaultSKSpectrumStyle (conv::SuperkSample_t s, conv::SuperkMode_t m, TH1D * histogram);
   // When making an ereco hist, set directory/axis titles/line colour
   void FormatErecoHist(TH1 * hist);

} // util
} // t2kana
  
#endif


