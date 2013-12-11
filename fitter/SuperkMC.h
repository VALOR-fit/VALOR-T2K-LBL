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
#include "Utilities.h"

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
class SuperkDetectorFSI;

class SuperkMC: public TObject
{
public:
  SuperkMC();
 ~SuperkMC();

  // Apply systematic variations + oscillations (,decay,decoherence) to MC templates
  void UpdateTemplates(void);

  // Integrate current MC templates and get Ereco spectra (needed for fitting the data)
  bool  UpdateErecoSpectra  (void);
  TH1D* CurrErecoSpectrum   (conv::SuperkSample_t s) const;
  TH1D* CurrErecoSpectrum   (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  double* CurrErecoSpectrumArray   (conv::SuperkSample_t s);
  double* CurrErecoSpectrumArray   (conv::SuperkSample_t s, conv::SuperkMode_t m);
      
  // MC event samples used to fill-in MC templates and POT exposure for each
  TTree * MCTree (conv::SuperkSample_t s, conv::SuperkMode_t m);
  double  MCPOT  (conv::SuperkSample_t s, conv::SuperkMode_t m);

  // Access MC templates
  TH2D * Template (conv::SuperkSample_t s, conv::SuperkMode_t m, conv::SKMCTemplateType_t type) const;
  TH2D * VarTemplate (conv::SuperkSample_t s, conv::SuperkMode_t m) const;

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
  void TellMeAboutTheDetectorFSI(SuperkDetectorFSI * skdetfsi) { fSKDetFSI  = skdetfsi; }

  // Calculate error in each 2D template bin due to finite MC statistics & fluctuate to study effect of MC stat err
  void CalculateMCStatError(void); // called at init, using the original templates
  void ApplyMCStatFluctationToTemplates(bool shape_only=true); // fluctuate contents of MC templates using MC stat error

  // Draw inputs:
  // ...
  void Draw(TVirtualPad * c, Option_t * opt = "") const;

  // Hack for x_mec until we have a parameter we can tweak
  void SetXMEC(double xmec)              { fXMEC = xmec; }

  // Allow osc prob > 1 in 2 flavour oscillation?
  void SetUnphysicalOscProb(bool unphys) { fUnphysicalOscProb = unphys; }

private:

  // Init and clean-up
  void Init    (void);
  void CleanUp (void);
  void Book    (void);

  // Fill MC templates from the input SK event tree or from other sources
  void FillTemplate (conv::SuperkSample_t s, conv::SuperkMode_t m);
  void FillSpecialTemplates (void);

  // Apply any flux tuning based on external (NA61,...) data to nominal MC
  void ApplyFluxExtTuningToTemplates (void); 

  // Apply BANFF flux + cross section tuning (based on ND280 data) on top of NA61-based flux tuned MC
  void ApplyBANFFTuningToTemplates (void); 

  // Apply systematic variations to MC templates
  void ApplySystematicsToTemplates                   (void); 
  void ApplySystematicsToTemplates_Multiplicative    (void); 
  void ApplyMECRenormalisationToTemplates            (void); 
  void ApplySystematicsToTemplates_NonMultiplicative (void); 
  bool ApplyEnergyScaleSystematicToTotalSpectrum     (void);

  // Apply physics hypothesis (oscillations,...) to MC templates
  void ApplyPhysHypothesisToTemplates (void); 

  // Ignore nue_signal templates for all but the 3-flavour oscillation hypothesis
  bool SkipTemplate(conv::SuperkMode_t m);

  // Debug method, prints out some number of events, and also saves templates & spectra
  void NevDebug(conv::SKMCTemplateType_t t, string name);
  void MCStatErrDebug(void);

  // Looks to see if there are any events for the given Etrue in the sample
  bool IsThereAnEvent(conv::SuperkSample_t s, double Etrue) const;

  // Integrate current MC templates and get Etrue spectra (needed for applying physics hypothesis)
  void FillEtrueSpectra(void);

  //
  // Private data members
  //

  double fPOT;  ///< POT for actual data-set (MC is scaled to that exposure)

  // input MC templates (nominal, nominal + systematics, nominal + systematics + physics hypothesis) and Ereco spectra
  TH2D *     fMCTemplOri[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal (doesn't even have pot weighting applied)
  TH2D *     fMCTemplNom[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning
#ifndef PRODUCTIONMODE
  TH2D *     fMCTemplTun[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning + nd280/tuning
  TH2D *     fMCTemplSys[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning + nd280/tuning + systematics 
  TH2D *     fMCTemplFin[kns][knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
#endif
  TH2D *     fMCTemplMECMartini;           ///< Martini MEC template
  TH2D *     fMCTemplMECNieves;            ///< Nieves MEC template
  double     fMCTemplArrayMECMartini[conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu]; ///< Martini MEC template
  double     fMCTemplArrayMECNieves [conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu]; ///< Nieves MEC template
  TH2D *     fMCVariableTempl[kns][knm];
  int        fMCVariableTemplArray1Rmu[kns][knm][conv::kMCNEtrueBins][conv::kXsNErecoBins1Rmu];
  int        fMCVariableTemplArray1Re [kns][knm][conv::kMCNEtrueBins][conv::kXsNErecoBins1Re];
  int        fMCTemplArrayOri1Rmu[knm][conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu];
  int        fMCTemplArrayOri1Re [knm][conv::kNEtrueBins1Re][conv::kNErecoBins1Re];
  double     fMCTemplArrayNom1Rmu[knm][conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu];
  double     fMCTemplArrayNom1Re [knm][conv::kNEtrueBins1Re][conv::kNErecoBins1Re];
  double     fMCTemplArrayTun1Rmu[knm][conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu];
  double     fMCTemplArrayTun1Re [knm][conv::kNEtrueBins1Re][conv::kNErecoBins1Re];
  //double     fMCTemplArraySys[kns][knm][conv::kNEtrueBins][conv::kNErecoBins1Rmu];
  double     fMCTemplArrayFin1Rmu[knm][conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu];
  double     fMCTemplArrayFin1Re [knm][conv::kNEtrueBins1Re][conv::kNErecoBins1Re];
  double     fMCStatFracErr1Rmu[knm][conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu];
  double     fMCStatFracErr1Re [knm][conv::kNEtrueBins1Re][conv::kNErecoBins1Re];

  TH1D *     fMCErecoSpectrum[kns][knm];   ///< Ereco MC spectrum (sample,mode)     : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
  TH1D *     fMCErecoSpectrumTot[kns];     ///< Ereco MC spectrum (sample)          : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
  TH1D *     fMCEtrueSpectrumTot[kns];     ///< Etrue MC spectrum (sample)          : nominal
  double     fMCEtrueSpectrumTotArray1Rmu[conv::kNEtrueBins1Rmu];
  double     fMCEtrueSpectrumTotArray1Re[conv::kNEtrueBins1Re];
  double     fMCErecoSpectrumTotArray1Rmu[conv::kNErecoBins1Rmu];
  double     fMCErecoSpectrumTotArray1Re[conv::kNErecoBins1Re];
  double     fMCErecoSpectrumArray1Rmu[conv::knm][conv::kNErecoBins1Rmu];
  double     fMCErecoSpectrumArray1Re [conv::knm][conv::kNErecoBins1Re];

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
  
  TH1D*      fMCPOTHistNuMu;               ///< input MC event sample POT histogram : numu
  TH1D*      fMCPOTHistNuMuBar;            ///< input MC event sample POT histogram : numubar
  TH1D*      fMCPOTHistNuE;                ///< input MC event sample POT histogram : nue 
  TH1D*      fMCPOTHistNuEBar;             ///< input MC event sample POT histogram : nuebar
  TH1D*      fMCPOTHistSigNuE;             ///< input MC event sample POT histogram : nue (signal)
  double     fMCPOTNuMu;                   ///< MC event sample POT : numu
  double     fMCPOTNuMuBar;                ///< MC event sample POT : numubar
  double     fMCPOTNuE;                    ///< MC event sample POT : nue 
  double     fMCPOTNuEBar;                 ///< MC event sample POT : nuebar
  double     fMCPOTSigNuE;                 ///< MC event sample POT : nue (signal)
  float      fBrEreco;                     ///< branch var for reading input MC mini-DST
  float      fBrEtrue;                     ///< branch var for reading input MC mini-DST
  int        fBrGenCode;                   ///< branch var for reading input MC mini-DST
  double     fBrWght;                      ///< branch var for reading input MC mini-DST
  float      fBrWght1Re;                      ///< branch var for reading input MC mini-DST
  int        fBrEvFlag;                      ///< branch var for reading input MC mini-DST

  // utility objects for computing cross-section and flux weights
  NuFlux *         fNuFlux;                ///< Neutrino flux inputs
  NuXSec *         fNuXSec;                ///< Neutrino cross-section inputs
  FSI *            fFSI;                   ///< FSI inputs
  SuperkDetector * fSKDet;                 ///< SuperK detector efficiency error etc
  SuperkDetectorFSI * fSKDetFSI;           ///< SuperK detector efficiency error etc

  // utility objects for defining physics hypothesis (oscillations, decay, ...) & computing probabilities
  PhysHypothesis fPhysHypothesis;          ///< Physics hypothesis 
  NuOsc3Flavour  fNuOsc3Flavour;           ///< 3-flavour neutrino osc. probability calculations

  double fXMEC; //hack x_mec until we have a parameter we can use

  bool fUnphysicalOscProb; // Allow unphysical region (probs not in range [0,1]) in 2 flavour case
  
  util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject

ClassDef(SuperkMC, 5);
};

} //t2kana

#endif
