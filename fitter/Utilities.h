//____________________________________________________________________________
/*!

\class    t2kana::utils

\brief    Utilities

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_utils_h_
#define _t2kana_utils_h_

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <complex>
#include <cassert>

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
#include "Enumerations.h"
#include "ConfLevel.h"
#include "PDG.h"
#include "UtlNumStat.h"
#include "UtlStr.h"
#include "UtlStyle.h"
#include "UtlXML.h"


using std::string;
using std::vector;

namespace t2kana {
namespace util {

  //
  // Utility functions for T2KAna enumerations
  //
  
  string AsString           (conv::OscAnaType_t       a);
  string AsString           (conv::SKMCTemplateType_t t);
  string AsString           (conv::FluxType_t         f);
/*string AsStringPresentTex (conv::PhysParamId_t      p);*/
  string AsStringRooTex     (conv::SuperkSample_t s, conv::SuperkMode_t m, string append_begin="", string append_end="");
  string AsStringCompact    (conv::SuperkSample_t s, conv::SuperkMode_t m, string append_begin="", string append_end="");
  string AsStringCompact_KDI(conv::SuperkMode_t m);

  //temporary methods to get around spline naming
  string AsStringCompact_splinehack(conv::SuperkSample_t s);
  string AsStringCompact_splinehack(conv::SuperkMode_t   m);
  
  conv::OscAnaType_t       OscAnaFromString          (string h);
  conv::SKMCTemplateType_t SKMCTemplateTypeFromString(string h);

  conv::OscChannelId_t GetOscChannel(int nu_start, int nu_end);
  
  //
  // Reaction mode utilities etc
  //

  bool UseEventToFillInTemplate(conv::SuperkMode_t m, int nu_pdg_code, int neut_reaction_mode);

  bool IsCCQE     (int neut_reaction_mode); 
  bool IsCC1pi    (int neut_reaction_mode);
  bool IsCCcoh    (int neut_reaction_mode); 
  bool IsCCoth    (int neut_reaction_mode); // this is a catch all. if new categories are added, please alter this method too 
  //  bool IsCCoth_e  (int neut_reaction_mode); 
  //bool IsCCoth_mu (int neut_reaction_mode);
  bool IsNC1pi0   (int neut_reaction_mode);
  bool IsNC1piC   (int neut_reaction_mode);
  bool IsNCcoh    (int neut_reaction_mode); 
  bool IsNCoth    (int neut_reaction_mode); // this is a catch all. if new categories are added, please alter this method too 
  bool IsNCoth_e  (int neut_reaction_mode);
  bool IsNCoth_mu (int neut_reaction_mode);

  //
  // Etc
  //
  
  void PrintCompact (const TVectorD & v, int w=15, int prec=9, Option_t * opt="", ostream & stream = std::cout);
  void PrintCompact (const TVectorD & current, const TVectorD & nominal, int w, int prec, Option_t * option, ostream & stream = std::cout);
  void PrintCompact (const TMatrixD & m, int w=15, int prec=9, ostream & stream = std::cout);
  void PrintCompactTex (const TMatrixD & m, int w=20, int prec=9, ostream & stream = std::cout, TString opt = "");

  //
  //
  //

  class ArrIdx {
  public:
    static ArrIdx * Instance();
    int SuperkSampleIdx (conv::SuperkSample_t s)                       { return fSuperkSampleIdx[s];    }
    int SuperkModeIdx   (conv::SuperkSample_t s, conv::SuperkMode_t m) { return fSuperkModeIdx  [s][m]; }
    int XSecSystIdx     (conv::PhysParamId_t p)                        { return fXSecSystIdx    [p];    }
    int RHCXSecSystIdx  (conv::PhysParamId_t p)                        { return fRHCXSecSystIdx    [p]; }
    int BANFFSystIdx    (conv::PhysParamId_t p)                        { return fBANFFSystIdx   [p];    }
    int RHCBANFFSystIdx (conv::PhysParamId_t p)                        { return fRHCBANFFSystIdx   [p]; }
    int FSISystIdx      (conv::PhysParamId_t p)                        { return fFSISystIdx     [p];    }
    int SKDetFSISystIdx (conv::PhysParamId_t p)                        { return fSKDetFSISystIdx[p];    }
    int RHCSKDetFSISystIdx (conv::PhysParamId_t p)                     { return fRHCSKDetFSISystIdx[p]; }
    int SKDet1ReSystIdx (conv::PhysParamId_t p)                        { return fSKDet1ReSystIdx[p];    }
    int PhysParamIdx    (conv::PhysParamId_t p)                        { return fPhysParamIdx   [p];    }
  private:
    ArrIdx();
    ArrIdx(const ArrIdx &) {}
   ~ArrIdx() {}
    void Init();
    int fSuperkSampleIdx [2*conv::kns];                  ///< maps sample to kSKSamples array index
    int fSuperkModeIdx   [2*conv::kns][2*conv::knm];     ///< maps sample/mode to kSKReactionModes array index
    int fXSecSystIdx     [conv::knppmax];                ///< maps physics param to kListOfXSecSyst array index
    int fRHCXSecSystIdx  [conv::knppmax];                ///< maps physics param to kListOfRHCXSecSyst array index
    int fBANFFSystIdx    [conv::knppmax];                ///< maps physics param to BANFF covariance matrix
    int fRHCBANFFSystIdx [conv::knppmax];                ///< maps physics param to BANFF covariance matrix
    int fFSISystIdx      [conv::knppmax];                ///< maps physics param to FSI covariance matrix
    int fSKDetFSISystIdx [conv::knppmax];                ///< maps physics param to FSI covariance matrix
    int fRHCSKDetFSISystIdx [conv::knppmax];             ///< maps physics param to FSI covariance matrix
    int fSKDet1ReSystIdx [conv::knppmax];                ///< maps physics param to 1Re eff covariance matrix
    int fPhysParamIdx    [conv::knppmax];                ///< maps physics param
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


