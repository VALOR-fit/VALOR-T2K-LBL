//____________________________________________________________________________
/*!

\class    t2kana::utils

\brief    Utilities

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_utils_h_
#define _t2kana_utils_h_

#include <string>
#include <vector>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TXMLEngine.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TTree.h>

class TGraph;
class TH1;
class TH1F;
class TH1D;
class TH2D;
class TAxis;

#include "Conventions.h"

using std::string;
using std::vector;

namespace t2kana {
namespace util {

  double DeltaChisquare   (t2kana::conv::ConfLevel_t cl, int nfitparams);
  double ChisqGauss       (double ndata, double nmc);
  double ChisqPoisson     (double ndata, double nmc);
  double NuMuSurvProbOsc2D(double sinsq2theta, double dmsq, double Et);
  double NuMuSurvProbOsc2D(double sinsq2theta, double dmsq, double Et, t2kana::conv::SuperkMode_t m);
  double NuMuSurvProbDecay(double sinsq2theta, double alpha, double Et);
  double NuMuSurvProbDecay(double sinsq2theta, double alpha, double Et, t2kana::conv::SuperkMode_t m);
  double NuMuSurvProbDecoh(double sinsq2theta, double musq, double Et);
  double NuMuSurvProbDecoh(double sinsq2theta, double musq, double Et, t2kana::conv::SuperkMode_t m);

  // Cholesky decomposition

  TMatrixD CholeskyDecomposition                     (TMatrixD& cov_matrix); ///< Cholesky decomposition. Return lower triangular matrix.
  TVectorD CholeskyGenerateCorrelatedParams          (TMatrixD& cholesky_triangular, TVectorD& mean_parm); ///< Generate a vector of correlated params
  TVectorD CholeskyGenerateCorrelatedParamVariations (TMatrixD& cholesky_triangular); ///< Generate a vector of correlated param variations
  TVectorD CholeskyCalculateCorrelatedParamVariations(TMatrixD& cholesky_triangular, TVectorD& g_uncorrelated);

  // Functions randomizing oscillation params (for use when params are marginalized)
  double RandomizeSinsq2Theta12 (void); ///< pull random from gaussian with mean and sigma taken from global fit
  double RandomizeSinsq2Theta13 (void); ///< uniform from 0 up to CHOOZ limit 
  double RandomizeDeltaCP       (void); ///< uniform from 0 to 2*pi
  double RandomizeDm12sq        (void); ///< pull random from gaussian with mean and sigma taken from global fit

  //

  string AsString (string hypothesis, double sinsq2theta12, double sinsq2theta13, double sinsq2theta23, double deltacp, double dm12sq, double dm23sq, Option_t * option="all");
  string AsString (t2kana::conv::NuMassHierarchy_t m);
  string AsString (t2kana::conv::PhysHypothesis_t  h);
  string AsString (t2kana::conv::PhysParamId_t     p);
  string AsString (t2kana::conv::SuperkSample_t    s);
  string AsString (t2kana::conv::SuperkMode_t      m);
  string AsString (t2kana::conv::OscChannelId_t    o);
  string AsStringTex (t2kana::conv::SuperkSample_t    s);
  string AsStringTex (t2kana::conv::SuperkMode_t      m);

  t2kana::conv::PhysHypothesis_t FromString(string h);
  t2kana::conv::SuperkSample_t   SKSampleFromString(string h);

  t2kana::conv::OscChannelId_t GetOscChannel(int nu_start, int nu_end);

  bool IsNuMu         (t2kana::conv::SuperkMode_t  m);
  bool IsNuMuBar      (t2kana::conv::SuperkMode_t  m);
  bool IsBeamNuE      (t2kana::conv::SuperkMode_t  m);
  bool IsSignalNuE    (t2kana::conv::SuperkMode_t  m);
  bool IsCC           (t2kana::conv::SuperkMode_t  m);
  bool IsNC           (t2kana::conv::SuperkMode_t  m);
  bool IsSystematic   (t2kana::conv::PhysParamId_t p);
  bool IsNeutrino     (int pdg_code);
  bool IsAntiNeutrino (int pdg_code);

  // string utilities

  vector<string> Tokens     (string input, string delimiter);
  string         TrimSpaces (string input);

  // XML utilities

  XMLNodePointer_t FindNode        (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  bool             NodeExists      (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TH1F *           XML2TH1F        (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TH1D *           XML2TH1D        (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TH2D *           XML2TH2D        (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TTree *          XML2TTree       (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TAxis *          XML2BinInfo     (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  double           XML2Dbl         (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  vector<double>   XML2DblArray    (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  string           XML2String      (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TMatrixDSym *    XML2TMatrixDSym (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  TMatrixD *       XML2TMatrixD    (TXMLEngine * eng, XMLNodePointer_t topnode, string path);
  
  // reaction modes etc

  bool UseEventToFillInTemplate(t2kana::conv::SuperkMode_t m, int nu_pdg_code, int neut_reaction_mode);

  bool IsCCQE     (int neut_reaction_mode); 
  bool IsCC1pi    (int neut_reaction_mode);
  bool IsCCcoh    (int neut_reaction_mode); 
  bool IsCCoth_e  (int neut_reaction_mode); 
  bool IsCCoth_mu (int neut_reaction_mode);
  bool IsNC1pi0   (int neut_reaction_mode);
  bool IsNCcoh    (int neut_reaction_mode); 
  bool IsNCoth_e  (int neut_reaction_mode);
  bool IsNCall    (int neut_reaction_mode); 

  // ROOT style/formatting

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);


  // handle variable binning in fluc covariance matrices

  int    FluxFix2VarBinNumber             (int fixBinNo);
  double FluxFractionsFixedBinRelToVarBin (int fixBinNo, TH1D* fixBinFluxHisto, TH1D* varBinFluxHisto);

  // etc

  int Sign (double tweaking_dial);

} // util
} // t2kana
  
#endif


