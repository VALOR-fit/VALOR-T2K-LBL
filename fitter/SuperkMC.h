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

class SuperkMC: public TObject
{
public:
  SuperkMC();
 ~SuperkMC();

  // Apply tuning, systematic variations + oscillations (,decay,decoherence) to nominal MC templates
  void UpdateTemplates(void);

  // Integrate current MC templates and get Ereco spectra (needed for fitting the data)
  void  UpdateErecoSpectra  (void);
  TH1D* CurrErecoSpectrum   (t2kana::conv::SuperkSample_t s) const;
  TH1D* CurrErecoSpectrum   (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;

  // Integrate current MC templates and get Etrue spectra (needed for flux normalisation error calculation)
  void  UpdateEtrueSpectra  (void);
  TH1D* Curr1RmuNuMuEtrue   (void) const { return f1RmNuMuEt;    } ///< current 1Rmu (allowing for applied weights) Etrue distribution for numu
  TH1D* Curr1RmuNuMuBarEtrue(void) const { return f1RmNuMuBarEt; } ///<      -//-  for numubar
  TH1D* Curr1RmuNuEEtrue    (void) const { return f1RmNuEEt;     } ///<      -//-  for nue
  TH1D* Curr1ReNuMuEtrue    (void) const { return f1ReNuMuEt;    } ///< current 1Re (allowing for applied weights) Etrue distribution for numu
  TH1D* Curr1ReNuMuBarEtrue (void) const { return f1ReNuMuBarEt; } ///<      -//-  for numubar
  TH1D* Curr1ReNuEEtrue     (void) const { return f1ReNuEEt;     } ///<      -//-  for nue
  TH1D* Curr1ReNuESigEtrue  (void) const { return f1ReNuESigEt;  } ///<      -//-  for nue signal
      
  // MC event samples used to fill-in MC templates and POT exposure for each
  TTree * MCTree (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);
  double  MCPOT  (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);

  // Actual MC templates
  TH2D * Template (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, t2kana::conv::SKMCTemplateType_t type) const;

  // POT of analysed data sample (MC is scaled to that POT value)
  double AnalyzedSamplePOT(void) const;

  // Systematic variations
  void TweakEnergyScale                    (double twkdial);
  void TweakEfficiency1RmuCCQEOther        (double twkdial); ///< "other": all but ring-counting
  void TweakEfficiency1RmuCCQERingCounting (bool force_input_twkdials=false, double twkdial1=0, double twkdial2=0, double twkdial3=0, bool correlated=true); 
  void TweakEfficiency1RmuCCnonQE          (double twkdial);
  void TweakEfficiency1RmuNC               (double twkdial);
  void TweakEfficiency1RmuNuECC            (double twkdial);
  void TweakEfficiency1ReSig               (double twkdial);
  void TweakEfficiency1ReBkg               (double twkdial);

  // Calculate weight due to efficiency tweaks
  double  EfficiencyWeight (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, double Ereco) const;

  // Info about the energy scale error and current energy scale tweaking factor
  double  EnergyScaleRelErr   (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const; ///< energy scale err
  double  EnergyScale         (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const; ///< energy scale

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
      
  // Physics hypothesis 
  PhysHypothesis & Hypothesis(void);

  // Hooks for obtaining necessary flux and cross-section information
  void TellMeAboutTheFlux       (NuFlux * flux) { fNuFlux     = flux; }
  void TellMeAboutCrossSections (NuXSec * xsec) { fGenInputs  = xsec; }

  // Draw inputs:
  // ...
  void Draw(TVirtualPad * c, Option_t * opt = "") const;

  // Print      
  void PrintSystematicVariations (void) const;

private:

  void Init         (void);
  void CleanUp      (void);
  void Book         (void);
  void FillTemplate (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);

  ///< Apply systematic variations to MC templates
  void ApplySystematicsToTemplates (void); 

  // Apply physics hypothesis (oscillations,...) to MC templates
  void ApplyPhysHypothesisToTemplates (void); 

  // Calculate the number of events in given (Ereco, Etrue) bin (applying all multiplicative systematics; i.e. all but energy scale)
  double NumOfEvents (double Ereco, double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;

  // Ignore nue_signal templates for all but the 3-flavour oscillation hypothesis
  bool SkipTemplate(t2kana::conv::SuperkMode_t m);
      
  //
  // private data members
  //

  double     fPOT;                         ///< POT for actual data-set (MC is scaled to that exposure)

  // systematic uncertainties: 1 sigma fractional errors & covariance matrices
  // > 1Rmu
  double     fEnergyScale1RmuRelErr;       ///< 1Rmu energy scale relative error
  double     fEff1RmuCCQEOthRelErr;        ///< 1Rmu CCQE other (all but ring-counting) efficiency relative error 
  TMatrixD * fEff1RmuCCQERCCovarRelErr;    ///< 1Rmu CCQE ring-counting efficiency relative error covariance matrix (3x3)
  double     fEff1RmuCCQERCCovarE1RelErr;  ///<  -  ~0.4 GeV (above RC cov matrix is 3x3, 3 Ereco bins: (0, fEff1RmuCCQERCCovarE1RelErr, fEff1RmuCCQERCCovarE2RelErr, oo))
  double     fEff1RmuCCQERCCovarE2RelErr;  ///<  -  ~1.1 GeV (above RC cov matrix is 3x3, 3 Ereco bins: (0, fEff1RmuCCQERCCovarE1RelErr, fEff1RmuCCQERCCovarE2RelErr, oo))
  TMatrixD * fEff1RmuCCQERCCovarChoLRelErr;///<  -  Cholesky decomposition (trangular matrix L) of fEff1RmuCCQERCCovarRelErr covariance matrix
  double     fEff1RmuCCnonQERelErr;        ///< 1Rmu CCnonQE efficiency relative error
  double     fEff1RmuNCRelErr;             ///< 1Rmu NC efficiency relative error
  double     fEff1RmuNuECCRelErr;          ///< 1Rmu nue CC efficiency relative error
  // > 1Re
  double     fEnergyScale1ReSigRelErr;     ///< 1Re signal energy scale relative error
  double     fEnergyScale1ReBkgRelErr;     ///< 1Re background energy scale relative error
  double     fEff1ReSigRelErr;             ///< 1Re signal efficiency relative error
  double     fEff1ReBkgRelErr;             ///< 1Re background efficiency relative error

  // systematic uncertainties: tweaking dials
  // > 1Rmu + 1Re
  double     fEnergyScaleTwkDial;          ///< 1Rmu/1Re Energy scale tweaking dial
  // > 1Rmu 
  double     fEff1RmuCCQEOthTwkDial;       ///< 1Rmu CCQE other (all-but-ring-counting) efficiency tweaking dial
  double     fEff1RmuCCQERC0TwkDial;       ///< 1Rmu CCQE ring-counting efficiency tweaking dial 0 (Ereco < fEff1RmuCCQERCCovarE1RelErr)
  double     fEff1RmuCCQERC1TwkDial;       ///< 1Rmu CCQE ring-counting efficiency tweaking dial 1 (fEff1RmuCCQERCCovarE1RelErr <= Ereco <= fEff1RmuCCQERCCovarE2RelErr)
  double     fEff1RmuCCQERC2TwkDial;       ///< 1Rmu CCQE ring-counting efficiency tweaking dial 2 (Ereco > fEff1RmuCCQERCCovarE2RelErr)
  double     fEff1RmuCCnonQETwkDial;       ///< 1Rmu CCnonQE efficiency tweaking dial
  double     fEff1RmuNCTwkDial;            ///< 1Rmu NC efficiency tweaking dial
  double     fEff1RmuNuECCTwkDial;         ///< 1Rmu nue CC efficiency tweaking dial
  // > 1Re
  double     fEff1ReSigTwkDial;            ///< 1Re signal effiency tweaking dial
  double     fEff1ReBkgTwkDial;            ///< 1Re background effiency tweaking dial

  // input MC templates (nominal, nominal + systematics, nominal + systematics + physics hypothesis) and Ereco spectra
  TH2D *     fMCTemplNom[kns][knm];        ///< MC template(sample,mode): Ereco-Etrue (nominal)
  TH2D *     fMCTemplSys[kns][knm];        ///< MC template(sample,mode): Ereco-Etrue + systematic variations only
  TH2D *     fMCTemplFin[kns][knm];        ///< MC template(sample,mode): Ereco-Etrue + systematic variations + physics hypothesis (final)
  TH1D *     fMCErecoSpectrum[kns][knm];   ///< MC Ereco spectrum(sample,mode)        + systematic variations + physics hypothesis (final)
  TH1D *     fMCErecoSpectrumTot[kns];     ///< MC total Ereco spectrum (sample)      + systematic variations + physics hypothesis (final)

  // input event trees and normalisation factors
  TTree*     fMCTreeNuMu1Rmu;              ///< input MC event tree /  sample: numu                        / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuMuBar1Rmu;           ///< input MC event tree /  sample: numubar                     / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuE1Rmu;               ///< input MC event tree /  sample: nue (intrinsic)             / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuESig1Rmu;            ///< input MC event tree /  sample: nue (signal 100% numu->nue) / cuts: 1-ring mu-like 
  TTree*     fMCTreeNuMu1Re;               ///< input MC event tree /  sample: numu                        / cuts: 1-ring e-like 
  TTree*     fMCTreeNuMuBar1Re;            ///< input MC event tree /  sample: numubar                     / cuts: 1-ring e-like 
  TTree*     fMCTreeNuE1Re;                ///< input MC event tree /  sample: nue (intrinsic)             / cuts: 1-ring e-like 
  TTree*     fMCTreeNuESig1Re;             ///< input MC event tree /  sample: nue (signal 100% numu->nue) / cuts: 1-ring e-like 
  double     fMCPOTNuMu;                   ///< input MC event sample POT : numu
  double     fMCPOTNuMuBar;                ///< input MC event sample POT : numubar
  double     fMCPOTNuE;                    ///< input MC event sample POT : nue (bkg)
  double     fMCPOTNuESig;                 ///< input MC event sample POT : nue (sig)
  double     fBrEreco;                     ///< branch var for reading input event sample tree
  double     fBrEtrue;                     ///< branch var for reading input event sample tree
  int        fBrPdg;                       ///< branch var for reading input event sample tree
  int        fBrGenCode;                   ///< branch var for reading input event sample tree

  // histograms for Nsk/Nnd280 (beam norm) error calc
  TH1D *     f1RmNuMuEt;                   ///< Etrue distribution for 1-Ring mu-like numu    
  TH1D *     f1RmNuMuBarEt;                ///< Etrue distribution for 1-Ring mu-like numubar 
  TH1D *     f1RmNuEEt;                    ///< Etrue distribution for 1-Ring mu-like nue     
  TH1D *     f1ReNuMuEt;                   ///< Etrue distribution for 1-Ring mu-like numu
  TH1D *     f1ReNuMuBarEt;                ///< Etrue distribution for 1-Ring mu-like numubar
  TH1D *     f1ReNuEEt;                    ///< Etrue distribution for 1-Ring mu-like nue (bkg)
  TH1D *     f1ReNuESigEt;                 ///< Etrue distribution for 1-Ring mu-like nue (sig)
      
  NuFlux *   fNuFlux;                      ///< Neutrino flux analysis inputs
  NuXSec *   fGenInputs;                   ///< Neutrino cross section analysis inputs

  PhysHypothesis fPhysHypothesis;          ///< Hypothesis for observed deficit: disappearance, decay, decoherence
  NuOsc3Flavour  fNuOsc3Flavour;           ///< 3-flavour neutrino osc. probability calculations

ClassDef(SuperkMC, 4);
};

} //t2kana

#endif
