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

  //
  // Statistics
  //

  double DeltaChisquare   (t2kana::conv::ConfLevel_t cl, int nfitparams);
  double ChisqGauss       (double ndata, double nmc);
  double ChisqPoisson     (double ndata, double nmc);

  // Cholesky decomposition

  TMatrixD CholeskyDecomposition                     (TMatrixD& cov);                 ///< Cholesky decomposition. Return lower triangular matrix.
  TVectorD CholeskyGenerateCorrelatedParams          (TMatrixD& Lch, TVectorD& mean); ///< Generate a vector of correlated params
  TVectorD CholeskyGenerateCorrelatedParamVariations (TMatrixD& Lch);                 ///< Generate a vector of correlated param variations
  TVectorD CholeskyCalculateCorrelatedParamVariations(TMatrixD& Lch, TVectorD& g);    ///<

  // Functions randomizing oscillation params (for use when params are marginalized)
  double RandomizeMixing23 (void); ///< pull from gaussian with mean & sigma from global fit (mixing param is either sinsq2theta or sinsqtheta depending on config)
  double RandomizeMixing12 (void); ///< pull from gaussian with mean & sigma from global fit ( -//- )
  double RandomizeMixing13 (void); ///< uniform from 0 up to ~ small_factor x CHOOZ limit    ( -//- )
  double RandomizeDeltaCP  (void); ///< uniform from -pi to pi
  double RandomizeDm12sq   (void); ///< pull from gaussian with mean & sigma from global fit
  double RandomizeDm23sq   (void); ///< pull from gaussian with mean & sigma from global fit

  //

  int FindSampleSlot(t2kana::conv::SuperkSample_t s);
  int FindModeSlot  (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);

  t2kana::conv::FluxNormErrIdx_t FindFluxNormErrIdx (t2kana::conv::SuperkSample_t s);
  t2kana::conv::FluxNormErrIdx_t FindFluxNormErrIdx (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);

  string AsString    (t2kana::conv::PhysHypothesis_t hyp, double mix12, double mix13, double mix23, double deltacp, double dm12sq, double dm23sq);
  string AsString    (t2kana::conv::NuMassHierarchy_t m);
  string AsString    (t2kana::conv::NuOscMixParam_t   p);
  string AsString    (t2kana::conv::PhysHypothesis_t  h);
  string AsString    (t2kana::conv::PhysParamId_t     p);
  string AsString    (t2kana::conv::SuperkSample_t    s);
  string AsString    (t2kana::conv::SuperkMode_t      m);
  string AsString    (t2kana::conv::OscChannelId_t    o);
  string AsString    (t2kana::conv::FluxNormErrIdx_t  i);
  string AsStringTex (t2kana::conv::SuperkSample_t    s);
  string AsStringTex (t2kana::conv::SuperkMode_t      m);

  t2kana::conv::PhysHypothesis_t   PhysHypFromString(string h);
  t2kana::conv::SuperkSample_t     SKSampleFromString(string h);
  t2kana::conv::NuMassHierarchy_t  HierarchyFromString(string h);
  t2kana::conv::NuOscMixParam_t    MixingTypeFromString(string h);

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

  const int      NumErecoBins  (t2kana::conv::SuperkSample_t s);
  const double * ErecoBinEdges (t2kana::conv::SuperkSample_t s);

  //
  // String utilities
  //

  vector<string> Tokens     (string input, string delimiter);
  string         TrimSpaces (string input);

  //
  // XML utilities
  //

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
  
  //
  // Reaction mode utilities etc
  //

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

  //
  // Simple physics hypotheses (effective 2-flavour oscillation, decay and decoherence).
  // Code for 3-flavour oscillations incl. matter effects is in NuOsc3Flavour.cxx
  //

  double NuMuSurvProbOsc2Flavour (double mixing23, double dm23sq, double Et);
  double NuMuSurvProbOsc2Flavour (double mixing23, double dm23sq, double Et, t2kana::conv::SuperkMode_t m);
  double NuMuSurvProbDecay       (double sinsq2thetadecay, double alpha, double Et);
  double NuMuSurvProbDecay       (double sinsq2thetadecay, double alpha, double Et, t2kana::conv::SuperkMode_t m);
  double NuMuSurvProbDecoh       (double sinsq2thetadecoh, double musq, double Et);
  double NuMuSurvProbDecoh       (double sinsq2thetadecoh, double musq, double Et, t2kana::conv::SuperkMode_t m);

  //
  // ROOT style & formatting
  //

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);

  //
  // Handle variable binning in fluc covariance matrices
  //

  int    FluxFix2VarBinNumber             (int fixBinNo);
  double FluxFractionsFixedBinRelToVarBin (int fixBinNo, TH1D* fixBinFluxHisto, TH1D* varBinFluxHisto);

  //
  // Etc
  //

  int Sign (double tweaking_dial);

} // util
} // t2kana
  
#endif


