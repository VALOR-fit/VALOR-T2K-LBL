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

#include <TXMLEngine.h>
#include <TFile.h>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TH1D.h>
#include <TSpline.h>

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

  double DeltaChisquare   (conv::ConfLevel_t cl, int nfitparams);
  double ChisqGauss       (double ndata, double nmc);
  double ChisqPoisson     (double ndata, double nmc);

  //
  // Cholesky decomposition
  //

  TMatrixD CholeskyDecomposition                     (TMatrixD& cov);                 ///< Cholesky decomposition. Return lower triangular matrix.
  TVectorD CholeskyGenerateCorrelatedParams          (TMatrixD& Lch, TVectorD& mean); ///< Generate a vector of correlated params
  TVectorD CholeskyGenerateCorrelatedParamVariations (TMatrixD& Lch);                 ///< Generate a vector of correlated param variations
  TVectorD CholeskyCalculateCorrelatedParamVariations(TMatrixD& Lch, TVectorD& g);    ///<

  //
  // Functions randomizing oscillation params (for use when params are marginalized)
  //

  double RandomizeMixing23 (void); ///< pull from gaussian with mean & sigma from global fit (mixing param is either sinsq2theta or sinsqtheta depending on config)
  double RandomizeMixing12 (void); ///< pull from gaussian with mean & sigma from global fit ( -//- )
  double RandomizeMixing13 (void); ///< uniform from 0 up to ~ small_factor x CHOOZ limit    ( -//- )
  double RandomizeDeltaCP  (void); ///< uniform from -pi to pi
  double RandomizeDm12sq   (void); ///< pull from gaussian with mean & sigma from global fit
  double RandomizeDm23sq   (void); ///< pull from gaussian with mean & sigma from global fit

  //
  // Utility functions for T2KAna enumerations
  //

  void SetDefaultSKSpectrumStyle (conv::SuperkSample_t s, conv::SuperkMode_t m, TH1D * histogram);

  string AsString        (conv::PhysHypothesis_t hyp, double mix12, double mix13, double mix23, double deltacp, double dm12sq, double dm23sq);
  string AsString        (conv::OscAnaType_t       a);
  string AsString        (conv::NuMassHierarchy_t  m);
  string AsString        (conv::NuOscMixParam_t    p);
  string AsString        (conv::PhysHypothesis_t   h);
  string AsString        (conv::PhysParamId_t      p);
  string AsString        (conv::SuperkSample_t     s);
  string AsString        (conv::SuperkMode_t       m);
  string AsString        (conv::SKMCTemplateType_t t);
  string AsString        (conv::FluxType_t         f);
  string AsString        (conv::OscChannelId_t     o);
  string AsStringTex     (conv::SuperkSample_t     s);
  string AsStringTex     (conv::SuperkMode_t       m);
  string AsStringTex     (conv::PhysParamId_t      p);
  string AsStringRooTex  (conv::SuperkSample_t s, conv::SuperkMode_t m, string append_begin="", string append_end="");
  string AsStringRooTex  (conv::SuperkSample_t     s);
  string AsStringRooTex  (conv::SuperkMode_t       m);
  string AsStringCompact (conv::SuperkSample_t s, conv::SuperkMode_t m, string append_begin="", string append_end="");
  string AsStringCompact (conv::SuperkSample_t     s);
  string AsStringCompact (conv::SuperkMode_t       m);
  string AsStringCompact_Japan (conv::SuperkMode_t m);
  string AsStringCompact (conv::PhysParamId_t      p);
  string AsStringCompact_Japan (conv::PhysParamId_t p);
  string Pdg2String      (int pdg_code);

  conv::PhysHypothesis_t   PhysHypFromString         (string h);
  conv::SuperkSample_t     SKSampleFromString        (string h);
  conv::OscAnaType_t       OscAnaFromString          (string h);
  conv::NuMassHierarchy_t  HierarchyFromString       (string h);
  conv::NuOscMixParam_t    MixingTypeFromString      (string h);
  conv::PhysParamId_t      PhysParamIdFromString     (string h);
  conv::SKMCTemplateType_t SKMCTemplateTypeFromString(string h);

  conv::OscChannelId_t GetOscChannel(int nu_start, int nu_end);

  int SuperkMode2FluxNuPdg(conv::SuperkMode_t m);
  int SuperkMode2SKMCPdg  (conv::SuperkMode_t m);

  bool IsNuMu         (conv::SuperkMode_t  m);
  bool IsNuMuBar      (conv::SuperkMode_t  m);
  bool IsNuE          (conv::SuperkMode_t  m);
  bool IsNuEBar       (conv::SuperkMode_t  m);
  bool IsSignalNuE    (conv::SuperkMode_t  m);
  bool IsCC           (conv::SuperkMode_t  m);
  bool IsNC           (conv::SuperkMode_t  m);
  bool IsSystematic   (conv::PhysParamId_t p);
  bool IsNeutrino     (int pdg_code);
  bool IsAntiNeutrino (int pdg_code);

  const int      NumErecoBins  (conv::SuperkSample_t s);
  const double * ErecoBinEdges (conv::SuperkSample_t s);
  const int      NumVariableXsEtrueBins(conv::SuperkSample_t s, conv::SuperkMode_t m);
  const double * VariableXsEtrueBinEdges(conv::SuperkSample_t s, conv::SuperkMode_t m);
  const int      NumVariableMCEtrueBins(conv::SuperkSample_t s, conv::SuperkMode_t m);
  const double * VariableMCEtrueBinEdges(conv::SuperkSample_t s, conv::SuperkMode_t m);

  //
  // String utilities
  //

  vector<string> Tokens     (string input, string delimiter);
  string         TrimSpaces (string input);

  //
  // XML utilities
  //

  XMLNodePointer_t FindNode        (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  bool             NodeExists      (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);

  double           XML2Dbl         (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  vector<double>   XML2DblArray    (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  string           XML2String      (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);

  string           XML2RootFileName(TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  string           XML2RootObjName (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);

  TFile *          XML2TFile       (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TH1F *           XML2TH1F        (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TH1D *           XML2TH1D        (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TH2D *           XML2TH2D        (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TTree *          XML2TTree       (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TVectorD *       XML2TVectorD    (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TMatrixDSym *    XML2TMatrixDSym (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TMatrixD *       XML2TMatrixD    (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  TSpline3 *       XML2TSpline3    (TXMLEngine * xmleng, XMLNodePointer_t topnode, string path);
  
  //
  // Reaction mode utilities etc
  //

  bool UseEventToFillInTemplate(conv::SuperkMode_t m, int nu_pdg_code, int neut_reaction_mode);

  bool IsCCQE     (int neut_reaction_mode); 
  bool IsCC1pi    (int neut_reaction_mode);
  bool IsCCcoh    (int neut_reaction_mode); 
  bool IsCCoth    (int neut_reaction_mode); // this is a catch all. if new categories are added, please alter this method too 
  bool IsCCoth_e  (int neut_reaction_mode); 
  bool IsCCoth_mu (int neut_reaction_mode);
  bool IsNC1pi0   (int neut_reaction_mode);
  bool IsNC1piC   (int neut_reaction_mode);
  bool IsNCcoh    (int neut_reaction_mode); 
  bool IsNCoth    (int neut_reaction_mode); // this is a catch all. if new categories are added, please alter this method too 
  bool IsNCoth_e  (int neut_reaction_mode);
  bool IsNCoth_mu (int neut_reaction_mode);

  //
  // Simple physics hypotheses (effective 2-flavour oscillation, decay and decoherence).
  // Code for 3-flavour oscillations incl. matter effects is in NuOsc3Flavour.cxx
  //

  double NuMuSurvProbOsc2Flavour (double mixing23, double dm23sq, double Et);
  double NuMuSurvProbOsc2Flavour (double mixing23, double dm23sq, double Et, conv::SuperkMode_t m);
  double NuMuSurvProbDecay       (double sinsq2thetadecay, double alpha, double Et);
  double NuMuSurvProbDecay       (double sinsq2thetadecay, double alpha, double Et, conv::SuperkMode_t m);
  double NuMuSurvProbDecoh       (double sinsq2thetadecoh, double musq, double Et);
  double NuMuSurvProbDecoh       (double sinsq2thetadecoh, double musq, double Et, conv::SuperkMode_t m);

  //
  // ROOT style & formatting
  //

  void SetDefaultStyle(bool black_n_white=false);

  void Format(TGraph* gr, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
  void Format(TH1* hst, 
        int lcol, int lsty, int lwid, int mcol, int msty, double msiz);

  //
  // Etc
  //

  int Sign(double x);
  
  void PrintCompact (const TVectorD & v, int w=15, int prec=9, Option_t * opt="");
  void PrintCompact (const TVectorD & current, const TVectorD & nominal, int w, int prec, Option_t * option);
  void PrintCompact (const TMatrixD & m, int w=15, int prec=9);

  //
  //
  //

  class ArrIdx {
  public:
    static ArrIdx * Instance();
    int SuperkSampleIdx (conv::SuperkSample_t s)                       { return fSuperkSampleIdx[s];    }
    int SuperkModeIdx   (conv::SuperkSample_t s, conv::SuperkMode_t m) { return fSuperkModeIdx  [s][m]; }
    int XSecSystIdx     (conv::PhysParamId_t p)                        { return fXSecSystIdx    [p];    }
    int BANFFSystIdx    (conv::PhysParamId_t p)                        { return fBANFFSystIdx   [p];    }
    int FSISystIdx      (conv::PhysParamId_t p)                        { return fFSISystIdx     [p];    }
  private:
    ArrIdx();
    ArrIdx(const ArrIdx &) {}
   ~ArrIdx() {}
    void Init();
    int fSuperkSampleIdx [2*conv::kns];                  ///< maps sample to kSKSamples array index
    int fSuperkModeIdx   [2*conv::kns][2*conv::knm];     ///< maps sample/mode to kSKReactionModes array index
    int fXSecSystIdx     [conv::knppmax];                ///< maps physics param to kListOfXSecSyst array index
    int fBANFFSystIdx    [conv::knppmax];                ///< maps physics param to BANFF covariance matrix
    int fFSISystIdx    [conv::knppmax];                  ///< maps physics param to FSI covariance matrix
    static ArrIdx * fInstance;
    struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (ArrIdx::fInstance !=0) {
            delete ArrIdx::fInstance;
            ArrIdx::fInstance = 0;
         }
      }
    };
    friend struct Cleaner;
  };

} // util
} // t2kana
  
#endif


