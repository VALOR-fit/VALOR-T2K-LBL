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
#include <TH1D.h>

class TGraph;
class TH1;
class TH1F;
class TH1D;
class TH2D;
class TAxis;

#include "year1_numu_osc_analysis/Conventions.h"

using std::string;
using std::vector;

namespace t2kana {
namespace yr1numu {
namespace util {

  class TH1DErrEnv;

  double DeltaChisquare   (t2kana::yr1numu::conv::ConfLevel_t cl, int nfitparams);
  double ChisqGauss       (double ndata, double nmc);
  double ChisqPoisson     (double ndata, double nmc);
  double NuMuSurvProb2D   (double sinsq2theta, double dmsq, double Et);
  double NuMuSurvProb2D   (double sinsq2theta, double dmsq, double Et, t2kana::yr1numu::conv::SuperkMode_t m);

  // Cholesky decomposition

  TMatrixD CholeskyDecomposition      (TMatrixD& cov_matrix); ///< Cholesky decomposition. Return lower triangular matrix.
  TVectorD CholeskyGenerateCorrelated (TMatrixD& cholesky_triangular, TVectorD& mean_parm); ///< Generate a vector of correlated params

  string AsString(t2kana::yr1numu::conv::PhysParamId_t  p);
  string AsString(t2kana::yr1numu::conv::SuperkSample_t s);
  string AsString(t2kana::yr1numu::conv::SuperkMode_t   m);

  // string utilities

  vector<string> Tokens     (string input, string delimiter);
  string         TrimSpaces (string input);

  // XML utilities

  XMLNodePointer_t FindNode      (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TH1F *           XML2TH1F      (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TH1D *           XML2TH1D      (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TH2D *           XML2TH2D      (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TAxis *          XML2BinInfo   (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  double           XML2Dbl       (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TH1DErrEnv *     XML2TH1DErrEnv(TXMLEngine * eng, XMLNodePointer_t topnode, string path);

  // ROOT style/formatting

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);

  // etc

  int Sign(double tweaking_dial);

  // utility classes

  class TH1DErrEnv
  {
   public:
     TH1DErrEnv() : fLow(0), fHigh(0) { }
     TH1DErrEnv(TH1D * low, TH1D * high) : fLow(low), fHigh(high) { }
    ~TH1DErrEnv() { 
       if(fLow ) delete fLow;
       if(fHigh) delete fHigh;
     }
     TH1D * fLow;
     TH1D * fHigh;
  };

} // util
} // yr1numu
} // t2kana
  
#endif


