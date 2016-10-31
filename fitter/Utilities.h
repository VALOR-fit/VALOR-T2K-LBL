//____________________________________________________________________________
/*!
 
 \class    valor::utils
 
 \brief    Utilities
 
 \author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
 STFC, Rutherford Appleton Laboratory
 
 \created  Oct 20, 2010
 
 \copyright (c) 2010-2013, VaLOR group.
 
 */
//____________________________________________________________________________

#ifndef _valor_t2k_utils_h_
#define _valor_t2k_utils_h_

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
#include "T2KOpt.h"

#include "UtlT2KStyle.h"
#include "MakeFakeData.h"

#include "valor/Conventions/Opt.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/UtlPrint.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/UtlXML.h"
#include "valor/VArray/VALORArrayD.h"


using std::string;
using std::vector;

namespace valor {
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
        
        // Apply the energy scale systematic (used in both SKMC and SKToyArray)
        void ApplyEnergyScaleSystematic(conv::SuperkSample_t s, double * events_per_bin, double escale_tweak);
        
        void ApplyEnergyScaleSystematic2D(conv::SuperkSample_t s, double *events_per_bin, double escale_tweak); // NEW 2D TEMPLATES
        
        
        //
        // Reaction mode utilities etc
        //
        
        bool UseEventToFillInTemplate(conv::SuperkMode_t m, int nu_pdg_code, int neut_reaction_mode);
        
        bool IsCCQE     (int neut_reaction_mode);
        bool IsCC1pi    (int neut_reaction_mode);
        bool IsMEC      (int neut_reaction_mode);
        bool IsCCcoh    (int neut_reaction_mode);
        bool IsCCoth    (int neut_reaction_mode); // this is a catch all. if new categories are added, please alter this method too
        //  bool IsCCoth_e  (int neut_reaction_mode);
        //bool IsCCoth_mu (int neut_reaction_mode);
        bool IsNC1pi0   (int neut_reaction_mode);
        bool IsNC1piC   (int neut_reaction_mode);
        bool IsNCcoh    (int neut_reaction_mode);
        bool IsNCoth    (int neut_reaction_mode);
        bool IsNC1Gamma (int neut_reaction_mode); // NC 1 GAMMA
        
        //
        // Etc
        //
        
        // These internally decide what params to use. Deprecated - use the generic ones in fmwk UtlPrint.
        void PrintCompact (const TVectorD & v, int w=15, int prec=9, Option_t * opt="", std::ostream & stream = std::cout);
        void PrintCompact (const TVectorD & current, const TVectorD & nominal, int w, int prec, Option_t * option, std::ostream & stream = std::cout);
        
        // For python, hack
        conv::PhysParamId_t GetBANFFParam(int FHCpRHC_param_no);
        
        // Access compiler macro from interpreter
        int GetSystNSigma();
        
        //
        //
        //
        
        class ArrIdx {
        public:
            static ArrIdx * Instance();
            int SuperkSampleIdx     (conv::SuperkSample_t s)                       { return fSuperkSampleIdx[s];    }
            int SuperkModeIdx       (conv::SuperkSample_t s, conv::SuperkMode_t m) { return fSuperkModeIdx  [s][m]; }
            int XSecSystIdx         (conv::PhysParamId_t p)                        { return fXSecSystIdx        [p]; }
            int RHCXSecSystIdx      (conv::PhysParamId_t p)                        { return fRHCXSecSystIdx     [p]; }
            int BANFFSystIdx        (conv::PhysParamId_t p)                        { return fBANFFSystIdx       [p]; }
            int BANFFXSecSystIdx    (conv::PhysParamId_t p)                        { return fBANFFXSecSystIdx   [p]; }
            int RHCBANFFSystIdx     (conv::PhysParamId_t p)                        { return fRHCBANFFSystIdx    [p]; }
            int RHCBANFFXSecSystIdx (conv::PhysParamId_t p)                        { return fRHCBANFFXSecSystIdx[p]; }
            int FSISystIdx          (conv::PhysParamId_t p)                        { return fFSISystIdx         [p]; }
            int SKDetFSISystIdx     (conv::PhysParamId_t p)                        { return fSKDetFSISystIdx    [p]; }
            int RHCSKDetFSISystIdx  (conv::PhysParamId_t p)                        { return fRHCSKDetFSISystIdx [p]; }
            
            // NEW FHC+RHC JOINT
            int FHCpRHCSKDetFSISystIdx  (conv::PhysParamId_t p)                    { return fFHCpRHCSKDetFSISystIdx [p];    }
            int FHCpRHCBANFFSystIdx     (conv::PhysParamId_t p)                    { return fFHCpRHCBANFFSystIdx[p];        }
            //
            
            /// NEW MULTIRING 2
            int MultiRingSKDetFSISystIdx  (conv::PhysParamId_t p)               { return fMultiRingSKDetFSISystIdx [p]; }
            int FHCpRHCpMRSKDetFSISystIdx (conv::PhysParamId_t p)               { return fFHCpRHCpMRSKDetFSISystIdx[p]; }
            ///
            int SKDet1ReSystIdx (conv::PhysParamId_t p)                        { return fSKDet1ReSystIdx[p];    }
            int PhysParamIdx    (conv::PhysParamId_t p)                        { return fPhysParamIdx   [p];    }
        private:
            ArrIdx();
            ArrIdx(const ArrIdx &) {}
            ~ArrIdx() {}
            void Init();
            int fSuperkSampleIdx [2*conv::kns];                  ///< maps sample to kSKSamples array index
            int fSuperkModeIdx   [2*conv::kns][2*conv::knm];     ///< maps sample/mode to kSKReactionModes array index
            int fXSecSystIdx         [conv::knppmax];            ///< maps physics param to kListOfXSecSyst array index
            int fRHCXSecSystIdx      [conv::knppmax];            ///< maps physics param to kListOfRHCXSecSyst array index
            int fBANFFSystIdx        [conv::knppmax];            ///< maps physics param to BANFF covariance matrix
            int fRHCBANFFSystIdx     [conv::knppmax];            ///< maps physics param to BANFF covariance matrix
            int fBANFFXSecSystIdx    [conv::knppmax];            ///< maps physics param to BANFF covariance matrix, xsec only
            int fRHCBANFFXSecSystIdx [conv::knppmax];            ///< maps physics param to BANFF covariance matrix, xsec only
            int fFSISystIdx          [conv::knppmax];            ///< maps physics param to FSI covariance matrix
            int fSKDetFSISystIdx     [conv::knppmax];            ///< maps physics param to FSI covariance matrix
            int fRHCSKDetFSISystIdx  [conv::knppmax];            ///< maps physics param to FSI covariance matrix
            
            // NEW FHC+RHC JOINT
            int fFHCpRHCSKDetFSISystIdx  [conv::knppmax];        ///< maps physics param to FSI covariance matrix
            int fFHCpRHCBANFFSystIdx     [conv::knppmax];        ///< maps physics param to BANFF covariance matrix
            //
            
            /// NEW MULTIRING 2
            int fMultiRingSKDetFSISystIdx [conv::knppmax];       ///< maps physics param to FSI covariance matrix
            int fFHCpRHCpMRSKDetFSISystIdx[conv::knppmax];       ///< maps physics param to FSI covariance matrix
            ///
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
} // valor

#endif


