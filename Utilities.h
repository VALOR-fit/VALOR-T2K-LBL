//____________________________________________________________________________
/*!

\class    t2kana::year1_numuosc::utils

\brief    Simple utilities

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_year1_numu_analysis_utilities_h_
#define _t2kana_year1_numu_analysis_utilities_h_

#include <string>

class TGraph;
class TH1;

#include "year1_numu_osc_analysis/Conventions.h"

using std::string;

namespace t2kana {
namespace year1_numuosc {
namespace utilities {

  double DeltaChisquare (t2kana::year1_numuosc::conventions::ConfLevel_t cl, int nfitparams);
  double ChisqGauss     (double ndata, double nmc);
  double ChisqPoisson   (double ndata, double nmc);
  double NuMuSurvProb2D (double sinsq2theta, double dmsq, double Et);

  string AsString(t2kana::year1_numuosc::conventions::PhysParamId_t p);

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);

} // utilities
} // year1_numuosc
} // t2kana
  
#endif


