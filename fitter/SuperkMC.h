//____________________________________________________________________________
/*!

\class    t2kana::SuperkMC

\brief    SuperK Monte Carlo templates

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 15, 2011

*/
//____________________________________________________________________________

#ifndef _t2kana_sk_mc_templates_h_
#define _t2kana_sk_mc_templates_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TMatrixD.h>

#include "Conventions.h"
#include "PhysHypothesis.h"
#include "NuOsc3Flavour.h"

class TVirtualPad;
class TTree;
class TH1D;
class TH2D;

using t2kana::conv::kns;
using t2kana::conv::knm;

namespace t2kana {

class NuFlux;
class NuXSec;
class FSI;
class SuperkDetector;

class SuperkMC: public TObject
{
public:
  SuperkMC();
 ~SuperkMC();

  // Apply systematic variations + oscillations (,decay,decoherence) to MC templates
  void UpdateTemplates(void);

  // Integrate current MC templates and get Ereco spectra (needed for fitting the data)
  void  UpdateErecoSpectra  (void);
  TH1D* CurrErecoSpectrum   (conv::SuperkSample_t s) const;
  TH1D* CurrErecoSpectrum   (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
      
  // MC event samples used to fill-in MC templates and POT exposure for each
  TTree * MCTree (conv::SuperkSample_t s, conv::SuperkMode_t m);
  double  MCPOT  (conv::SuperkSample_t s, conv::SuperkMode_t m);

  // Access MC templates
  TH2D * Template (conv::SuperkSample_t s, conv::SuperkMode_t m, conv::SKMCTemplateType_t type) const;

  // POT of analysed data sample (MC is scaled to that POT value)
  double AnalyzedSamplePOT(void) const;

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
      
  // Physics hypothesis 
  PhysHypothesis & Hypothesis(void);

  // Hooks for obtaining necessary flux, cross-section and detector information
  void TellMeAboutTheFlux       (NuFlux *         flux ) { fNuFlux = flux; }
  void TellMeAboutCrossSections (NuXSec *         xsec ) { fNuXSec = xsec; }
  void TellMeAboutFSI           (FSI *            fsi  ) { fFSI    = fsi;  }
  void TellMeAboutTheDetector   (SuperkDetector * skdet) { fSKDet  = skdet; }

  // Draw inputs:
  // ...
  void Draw(TVirtualPad * c, Option_t * opt = "") const;

private:

  // Init and clean-up
  void Init    (void);
  void CleanUp (void);
  void Book    (void);

  // Fill MC templates from the input event trees
  void FillTemplate (conv::SuperkSample_t s, conv::SuperkMode_t m);

  // Apply any flux tuning based on external (NA61,...) data to nominal MC
  void ApplyFluxExtTuningToTemplates (void); 

  // Apply BANFF flux + cross section tuning (based on ND280 data) on top of NA61-based flux tuned MC
  void ApplyBANFFTuningToTemplates (void); 

  // Apply systematic variations to MC templates
  void ApplySystematicsToTemplates                   (void); 
  void ApplySystematicsToTemplates_Multiplicative    (void); 
  void ApplySystematicsToTemplates_NonMultiplicative (void); 

  // Apply physics hypothesis (oscillations,...) to MC templates
  void ApplyPhysHypothesisToTemplates (void); 

  // Ignore nue_signal templates for all but the 3-flavour oscillation hypothesis
  bool SkipTemplate(conv::SuperkMode_t m);

  //
  // Private data members
  //

  double fPOT;  ///< POT for actual data-set (MC is scaled to that exposure)

  // input MC templates (nominal, nominal + systematics, nominal + systematics + physics hypothesis) and Ereco spectra
  TH2D *     fMCTemplNom[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning
  TH2D *     fMCTemplTun[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning + nd280/tuning
  TH2D *     fMCTemplSys[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning + nd280/tuning + systematics 
  TH2D *     fMCTemplFin[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
  TH1D *     fMCErecoSpectrum[kns][knm];   ///< Ereco MC spectrum (sample,mode)     : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
  TH1D *     fMCErecoSpectrumTot[kns];     ///< Ereco MC spectrum (sample)          : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis

  // input event trees and normalisation factors
  TTree*     fMCTreeNuMu1Rmu;              ///< input MC event tree /  sample: numu                         / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuMuBar1Rmu;           ///< input MC event tree /  sample: numubar                      / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuE1Rmu;               ///< input MC event tree /  sample: nue (intrinsic)              / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuEBar1Rmu;            ///< input MC event tree /  sample: nuebar (intrinsic)           / cuts: 1-ring mu-like 
  TTree*     fMCTreeSigNuE1Rmu;            ///< input MC event tree /  sample: nue (signal, 100% numu->nue) / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuMu1Re;               ///< input MC event tree /  sample: numu                         / cuts: 1-ring e-like 
  TTree*     fMCTreeNuMuBar1Re;            ///< input MC event tree /  sample: numubar                      / cuts: 1-ring e-like 
  TTree*     fMCTreeNuE1Re;                ///< input MC event tree /  sample: nue (intrinsic)              / cuts: 1-ring e-like 
  TTree*     fMCTreeNuEBar1Re;             ///< input MC event tree /  sample: nuebar (intrinsic)           / cuts: 1-ring e-like 
  TTree*     fMCTreeSigNuE1Re;             ///< input MC event tree /  sample: nue (signal, 100% numu->nue) / cuts: 1-ring e-like 
  double     fMCPOTNuMu;                   ///< input MC event sample POT : numu
  double     fMCPOTNuMuBar;                ///< input MC event sample POT : numubar
  double     fMCPOTNuE;                    ///< input MC event sample POT : nue 
  double     fMCPOTNuEBar;                 ///< input MC event sample POT : nuebar
  double     fMCPOTSigNuE;                 ///< input MC event sample POT : nue (signal)
  double     fBrEreco;                     ///< branch var for reading input MC mini-DST
  double     fBrEtrue;                     ///< branch var for reading input MC mini-DST
  int        fBrPdg;                       ///< branch var for reading input MC mini-DST
  int        fBrGenCode;                   ///< branch var for reading input MC mini-DST
      
  // utility objects for computing cross-section and flux weights
  NuFlux *         fNuFlux;                ///< Neutrino flux inputs
  NuXSec *         fNuXSec;                ///< Neutrino cross-section inputs
  FSI *            fFSI;                   ///< FSI inputs
  SuperkDetector * fSKDet;                 ///< SuperK detector efficiency error etc

  // utility objects for defining physics hypothesis (oscillations, decay, ...) & computing probabilities
  PhysHypothesis fPhysHypothesis;          ///< Physics hypothesis 
  NuOsc3Flavour  fNuOsc3Flavour;           ///< 3-flavour neutrino osc. probability calculations

ClassDef(SuperkMC, 5);
};

} //t2kana

#endif
