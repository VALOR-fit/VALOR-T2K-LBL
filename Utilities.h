//____________________________________________________________________________
/*!

\class    t2kana::yr1numu::utils

\brief    Simple util

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_year1_numu_analysis_util_h_
#define _t2kana_year1_numu_analysis_util_h_

#include <string>

class TGraph;
class TH1;

#include "year1_numu_osc_analysis/Conventions.h"

using std::string;

namespace t2kana {
namespace yr1numu {
namespace util {

  double DeltaChisquare (t2kana::yr1numu::conv::ConfLevel_t cl, int nfitparams);
  double ChisqGauss     (double ndata, double nmc);
  double ChisqPoisson   (double ndata, double nmc);
  double NuMuSurvProb2D (double sinsq2theta, double dmsq, double Et);

  string AsString(t2kana::yr1numu::conv::PhysParamId_t p);

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);

} // util
} // yr1numu
} // t2kana
  
#endif


