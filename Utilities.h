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
#include <vector>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TXMLEngine.h>

class TGraph;
class TH1;
class TH1D;

#include "year1_numu_osc_analysis/Conventions.h"

using std::string;
using std::vector;

namespace t2kana {
namespace yr1numu {
namespace util {

  double DeltaChisquare   (t2kana::yr1numu::conv::ConfLevel_t cl, int nfitparams);
  double ChisqGauss       (double ndata, double nmc);
  double ChisqPoisson     (double ndata, double nmc);
  double NuMuSurvProb2D   (double sinsq2theta, double dmsq, double Et);
  double NuMuSurvProb2D   (double sinsq2theta, double dmsq, double Et, t2kana::yr1numu::conv::SuperkMode_t m);

  TVectorD RunCholeskyThingy(TMatrixD& cov_matrix, TVectorD& parm);

  string AsString(t2kana::yr1numu::conv::PhysParamId_t p);

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);

  vector<string> Tokens     (string input, string delimiter);
  string         TrimSpaces (string input);

  TH1D *           XML2TH1D(TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  XMLNodePointer_t FindNode(TXMLEngine * eng, XMLNodePointer_t topnode, string path);

} // util
} // yr1numu
} // t2kana
  
#endif


