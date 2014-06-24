//____________________________________________________________________________
/*!

\class    valor::SuperkMC

\brief    SuperK Monte Carlo templates

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 15, 2011

*/
//____________________________________________________________________________

#ifndef _valor_sk_mc_templates_h_
#define _valor_sk_mc_templates_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TMatrixD.h>
#include <vector>

#include "Conventions.h"
#include "FitCore/PhysHypothesis.h"
#include "Utilities/NuOsc3Flavour.h"
#include "Utilities/NuOscSterile.h"
#include "Utilities.h"
#include "Utilities/VALORArray.h"

class TVirtualPad;
class TTree;
class TH1D;
class TH2D;

namespace valor {

class NuFlux;
class NuXSec;
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
  TH1D* CurrErecoSpectrum   (conv::SuperkSample_t s) const; // Returns new object, mut be deleted
  TH1D* CurrErecoSpectrum   (conv::SuperkSample_t s, conv::SuperkMode_t m) const; // Returns new object, mut be deleted
  double*       CurrErecoSpectrumArray   (conv::SuperkSample_t s);
  const double* CurrErecoSpectrumArray   (conv::SuperkSample_t s) const;
  double*       CurrErecoSpectrumArray   (conv::SuperkSample_t s, conv::SuperkMode_t m);
  const double* CurrErecoSpectrumArray   (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
      
  // MC event samples used to fill-in MC templates and POT exposure for each
  TTree * MCTree (conv::SuperkSample_t s, conv::SuperkMode_t m);
  double  MCPOT  (conv::SuperkSample_t s, conv::SuperkMode_t m);
  
  // Access MC templates
  TH2D * Template (conv::SuperkSample_t s, conv::SuperkMode_t m, conv::SKMCTemplateType_t type) const;
  TH2D * VarTemplate (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  
  TH2D Get2DHist(conv::SuperkSample_t s, conv::SKMCTemplateType_t type) const;
  TH2D Get2DHist(conv::SuperkSample_t s, conv::SuperkMode_t m , conv::SKMCTemplateType_t type) const;

  // POT of analysed data sample (MC is scaled to that POT value)
  double AnalyzedSamplePOT(conv::SuperkSample_t s) const;

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
      
  // Physics hypothesis 
  PhysHypothesis & Hypothesis(void);

  // Hooks for obtaining necessary flux, cross-section and detector information
  void TellMeAboutTheFlux       (NuFlux *         flux )       { fNuFlux = flux; }
  void TellMeAboutCrossSections (NuXSec *         xsec )       { fNuXSec = xsec; }
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

  // Allow changing the mass hierarchy
  void SetMassHierarchy(conv::NuMassHierarchy_t hier) { fNuOsc3Flavour.SetMassHierarchy(hier); }
    
  void NevDebug(conv::SKMCTemplateType_t t, string name);

private:

  // Init and clean-up
  void Init    (void);
  void CleanUp (void);
  void Book    (void);

  // Fill MC templates from the input SK event tree or from other sources
  void FillTemplate (conv::SuperkSample_t s, conv::SuperkMode_t m);
  void FillSpecialTemplates (void);
  
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
 // void NevDebug(conv::SKMCTemplateType_t t, string name);
  void MCStatErrDebug(void);

  // Looks to see if there are any events for the given Etrue in the sample
  bool IsThereAnEvent(conv::SuperkSample_t s, double Etrue) const;

  // Integrate current MC templates and get Etrue spectra (needed for applying physics hypothesis)
  void FillEtrueSpectra(void);
  
  //
  // Private data members
  //
 
  double fPOT;     ///< POT for actual data-set (MC is scaled to that exposure)
  double fPOTRHC;  ///< POT for actual data-set in RHC sample (MC is scaled to that exposure)

  // input MC templates (nominal, nominal + systematics, nominal + systematics + physics hypothesis) and Ereco spectra
  TH2D *     fMCTemplOri[conv::kns][conv::knm];        ///< Ereco-Etrue MC template(sample,mode): nominal (doesn't even have pot weighting applied)
  TH2D *     fMCTemplNom[conv::kns][conv::knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning
  
  TH2D *     fMCTemplMECMartini;           ///< Martini MEC template
  TH2D *     fMCTemplRHCMECMartini;           ///< Martini MEC template    // ANTINU

  TH2D *     fMCTemplMECNieves;            ///< Nieves MEC template
  TH2D *     fMCTemplRHCMECNieves;            ///< Nieves MEC template     // ANTINU

  double     fMCTemplArrayMECMartini[conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu]; ///< Martini MEC template
  double     fMCTemplArrayRHCMECMartini[conv::kRHCNEtrueBins1Rmu][conv::kRHCNErecoBins1Rmu]; ///< Martini MEC template     // ANTINU

  double     fMCTemplArrayMECNieves [conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu]; ///< Nieves MEC template
  double     fMCTemplArrayRHCMECNieves [conv::kRHCNEtrueBins1Rmu][conv::kRHCNErecoBins1Rmu]; ///< Nieves MEC template     // ANTINU

  TH2D *     fMCVariableTempl[conv::kns][conv::knm];
  
  int        fMCVariableTemplArray1Rmu[conv::knm][conv::kMCNEtrueBins][conv::kXsNErecoBins1Rmu];
  int        fMCVariableTemplArray1Re [conv::knm][conv::kMCNEtrueBins][conv::kXsNErecoBins1Re];
  int        fMCVariableTemplArrayRHC1Rmu[conv::knm][conv::kMCNEtrueBins][conv::kRHCXsNErecoBins1Rmu];     // ANTINU
  int        fMCVariableTemplArrayRHC1Re [conv::knm][conv::kMCNEtrueBins][conv::kRHCXsNErecoBins1Re];      // ANTINU

  int        fMCTemplArrayOri1Rmu[conv::knm][conv::kNEtrueBins1Rmu][conv::kNErecoBins1Rmu];
  int        fMCTemplArrayOri1Re [conv::knm][conv::kNEtrueBins1Re][conv::kNErecoBins1Re];
  int        fMCTemplArrayOriRHC1Rmu[conv::knm][conv::kRHCNEtrueBins1Rmu][conv::kRHCNErecoBins1Rmu];    // ANTINU
  int        fMCTemplArrayOriRHC1Re [conv::knm][conv::kRHCNEtrueBins1Re][conv::kRHCNErecoBins1Re];     // ANTINU
  
  std::vector<VALORArray> fMCTemplArrayNom; //sample -> (mode, etrue, observable)
  std::vector<VALORArray> fMCTemplArrayTun; //sample -> (mode, etrue, observable)
  std::vector<VALORArray> fMCTemplArrayFin; //sample -> (mode, etrue, observable)
  
  std::vector<VALORArray> fMCStatFracErr;           // sample -> (mode, etrue, observable)
  std::vector<VALORArray> fMCEtrueSpectrumTotArray; // sample -> etrue
  std::vector<VALORArray> fMCErecoSpectrumTotArray; // sample -> observable
  std::vector<VALORArray> fMCErecoSpectrumArray   ; // sample -> (mode, observable)
  
  TH1D *     fMCErecoSpectrum[conv::kns][conv::knm];   ///< Ereco MC spectrum (sample,mode)     : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
  TH1D *     fMCErecoSpectrumTot[conv::kns];     ///< Ereco MC spectrum (sample)          : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
  TH1D *     fMCEtrueSpectrumTot[conv::kns];     ///< Etrue MC spectrum (sample)          : nominal
  
  // input event trees and normalisation factors
  TTree*     fMCTreeNuMu1Rmu;              ///< input MC event tree /  sample: numu                         / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuMuBar1Rmu;           ///< input MC event tree /  sample: numubar                      / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuE1Rmu;               ///< input MC event tree /  sample: nue (intrinsic)              / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuEBar1Rmu;            ///< input MC event tree /  sample: nuebar (intrinsic)           / cuts: 1-ring mu-like 
  TTree*     fMCTreeOscNuE1Rmu;            ///< input MC event tree /  sample: nue (signal, 100% numu->nue) / cuts: 1-ring mu-like 
  
  TTree*     fMCTreeNuMu1Re;               ///< input MC event tree /  sample: numu                         / cuts: 1-ring e-like 
  TTree*     fMCTreeNuMuBar1Re;            ///< input MC event tree /  sample: numubar                      / cuts: 1-ring e-like 
  TTree*     fMCTreeNuE1Re;                ///< input MC event tree /  sample: nue (intrinsic)              / cuts: 1-ring e-like 
  TTree*     fMCTreeNuEBar1Re;             ///< input MC event tree /  sample: nuebar (intrinsic)           / cuts: 1-ring e-like 
  TTree*     fMCTreeOscNuE1Re;             ///< input MC event tree /  sample: nue (signal, 100% numu->nue) / cuts: 1-ring e-like 
 
  TTree*     fMCTreeNuMuRHC1Rmu;              ///< input MC event tree /  sample: numu                         / cuts: 1-ring mu-like    // ANTINU
  TTree*     fMCTreeNuMuBarRHC1Rmu;           ///< input MC event tree /  sample: numubar                      / cuts: 1-ring mu-like    // ANTINU
  TTree*     fMCTreeNuERHC1Rmu;               ///< input MC event tree /  sample: nue (intrinsic)              / cuts: 1-ring mu-like    // ANTINU
  TTree*     fMCTreeNuEBarRHC1Rmu;            ///< input MC event tree /  sample: nuebar (intrinsic)           / cuts: 1-ring mu-like    // ANTINU
  TTree*     fMCTreeOscNuERHC1Rmu;            ///< input MC event tree /  sample: nue (signal, 100% numu->nue) / cuts: 1-ring mu-like    // ANTINU
  
  TTree*     fMCTreeNuMuRHC1Re;               ///< input MC event tree /  sample: numu                         / cuts: 1-ring e-like    // ANTINU
  TTree*     fMCTreeNuMuBarRHC1Re;            ///< input MC event tree /  sample: numubar                      / cuts: 1-ring e-like    // ANTINU
  TTree*     fMCTreeNuERHC1Re;                ///< input MC event tree /  sample: nue (intrinsic)              / cuts: 1-ring e-like    // ANTINU
  TTree*     fMCTreeNuEBarRHC1Re;             ///< input MC event tree /  sample: nuebar (intrinsic)           / cuts: 1-ring e-like    // ANTINU
  TTree*     fMCTreeOscNuERHC1Re;             ///< input MC event tree /  sample: nue (signal, 100% numu->nue) / cuts: 1-ring e-like    // ANTINU
 
  TH1D*      fMCPOTHistNuMu;               ///< input MC event sample POT histogram : numu
  TH1D*      fMCPOTHistNuMuBar;            ///< input MC event sample POT histogram : numubar
  TH1D*      fMCPOTHistNuE;                ///< input MC event sample POT histogram : nue 
  TH1D*      fMCPOTHistNuEBar;             ///< input MC event sample POT histogram : nuebar
  TH1D*      fMCPOTHistOscNuE;             ///< input MC event sample POT histogram : nue (signal)
  double     fMCPOTNuMu;                   ///< MC event sample POT : numu
  double     fMCPOTNuMuBar;                ///< MC event sample POT : numubar
  double     fMCPOTNuE;                    ///< MC event sample POT : nue 
  double     fMCPOTNuEBar;                 ///< MC event sample POT : nuebar
  double     fMCPOTOscNuE;                 ///< MC event sample POT : nue (signal)

  // ANTINU
  TH1D*      fMCPOTHistRHCNuMu;            ///< input MC event sample POT histogram : numu
  TH1D*      fMCPOTHistRHCNuMuBar;         ///< input MC event sample POT histogram : numubar
  TH1D*      fMCPOTHistRHCNuE;             ///< input MC event sample POT histogram : nue 
  TH1D*      fMCPOTHistRHCNuEBar;          ///< input MC event sample POT histogram : nue
  TH1D*      fMCPOTHistRHCOscNuE;          ///< input MC event sample POT histogram : nuebar (signal)
  double     fMCPOTRHCNuMu;                ///< MC event sample POT : numu
  double     fMCPOTRHCNuMuBar;             ///< MC event sample POT : numubar
  double     fMCPOTRHCNuE;                 ///< MC event sample POT : nue 
  double     fMCPOTRHCNuEBar;              ///< MC event sample POT : nue
  double     fMCPOTRHCOscNuE;              ///< MC event sample POT : nuebar (signal)
  ////

  float      fBrEreco;                     ///< branch var for reading input MC mini-DST
  float      fBrEtrue;                     ///< branch var for reading input MC mini-DST
  int        fBrGenCode;                   ///< branch var for reading input MC mini-DST
  double     fBrWght;                      ///< branch var for reading input MC mini-DST
  float      fBrWght1Re;                      ///< branch var for reading input MC mini-DST
  int        fBrEvFlag;                      ///< branch var for reading input MC mini-DST
  
  static const int kMaxVarBins = 1000;
  
  // Arrays for remembering which big bins have events, and which little bins for those big bigs have events.
  int fNNonEmptyVariableBins[conv::kns][conv::knm];               //!
  int fNNonEmptyFixedBins   [conv::kns][conv::knm][kMaxVarBins];  //!
  
  int fNonEmptyVariableBins [conv::kns][conv::knm][kMaxVarBins];  //! These contain both the reco and true bin numbers, masked with conv::kBigIndex;
  int * fNonEmptyFixedBins  [conv::kns][conv::knm][kMaxVarBins];  //! These contain both the reco and true bin numbers, masked with conv::kBigIndex;

  // utility objects for computing cross-section and flux weights
  NuFlux *         fNuFlux;                ///< Neutrino flux inputs
  NuXSec *         fNuXSec;                ///< Neutrino cross-section inputs
  SuperkDetectorFSI * fSKDetFSI;           ///< SuperK detector efficiency error etc

  // utility objects for defining physics hypothesis (oscillations, decay, ...) & computing probabilities
  PhysHypothesis fPhysHypothesis;          ///< Physics hypothesis 
  NuOsc3Flavour  fNuOsc3Flavour;           ///< 3-flavour neutrino osc. probability calculations

  double fXMEC; //hack x_mec until we have a parameter we can use

  bool fUnphysicalOscProb; // Allow unphysical region (probs not in range [0,1]) in 2 flavour case
  
  util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject

ClassDef(SuperkMC, 7);
};

} //valor

#endif
