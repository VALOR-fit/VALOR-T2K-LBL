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
  void    UpdateErecoSpectra  (void);
  TH1D*   CurrErecoSpectrum   (t2kana::conv::SuperkSample_t s) const                               { return fMCErecoSpectrumTot[s]; }
  TH1D*   CurrErecoSpectrum   (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const { return fMCErecoSpectrum[s][m]; }

  // Integrate current MC templates and get Etrue spectra (needed for flux normalisation error calculation)
  void    UpdateEtrueSpectra  (void);
  TH1D*   Curr1RmuNuMuEtrue   (void) const { return f1RmNuMuEt;    } ///< current (allowing for applied weights) Etrue distribution for numu
  TH1D*   Curr1RmuNuMuBarEtrue(void) const { return f1RmNuMuBarEt; } ///<      -//-  for numubar
  TH1D*   Curr1RmuNuEEtrue    (void) const { return f1RmNuEEt;     } ///<      -//-  for nue
      
  // MC event samples used to fill-in MC templates and POT exposure for each
  TTree * MCTree (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);
  double  MCPOT  (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);

  // POT of analysed sample (MC scaled to that POT-figure)
  double AnalyzedSamplePOT(void) const;

  // Systematic variations
  void TweakEnergyScale                 (double twkdial);
  void TweakOtherEfficiencyCCQE         (double twkdial); ///< all but ring-counting
  void TweakRingCountingEfficiencyCCQE  (bool force_input_twkdials=false, double twkdial1=0, double twkdial2=0, double twkdial3=0, bool correlated=true); 
  void TweakEfficiencyCCnonQE           (double twkdial);
  void TweakEfficiencyNC                (double twkdial);
  void TweakEfficiencyNuECC             (double twkdial);

  // Calculate weight due to efficiency tweaks
  double  EfficiencyWeight (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, double Ereco) const;

  // Info about the energy scale error and current energy scale tweaking factor
  double  EnergyScaleRelErr   (void) const; ///< energy scale err
  double  EnergyScale         (void) const; ///< energy scale

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
      
  // Physics hypothesis and calculation of survival probability for each hypothesis
  PhysHypothesis & Hypothesis(void);
  //double SurvivalProbability(double Etrue, t2kana::conv::SuperkMode_t m) const;

  // Hooks for obtaining necessary flux and cross-section information
  void TellMeAboutTheFlux       (NuFlux * flux) { fNuFlux     = flux; }
  void TellMeAboutCrossSections (NuXSec * xsec) { fGenInputs  = xsec; }

  // Draw inputs:
  // ...
  void Draw(TVirtualPad * c, Option_t * opt = "") const;
      
private:

  void Init           (void);
  void CleanUp        (void);
  void CreateTemplate (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);

  void ApplySystematicsToTemplates    (void); ///< apply systematic variations to MC templates
  void ApplyPhysHypothesisToTemplates (void); ///< apply physics hypothesis (oscillations,...) to MC templates

  // Loop over templates and calculate the number of events in given (Ereco, Etrue) bin
  double  NumOfEvents (double Ereco, double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;
      
  //
  // private data members
  //

  double     fPOT;                         ///< POT for actual data-set (MC is scaled to that exposure)

  double     fEnergyScaleRelErr;           ///< Energy scale relative error
  double     fOthEffCCQERelErr;            ///< CCQE other (all but ring-counting) efficiency relative error 
  TMatrixD * fEffCCQERelErrCovarRC;        ///< CCQE ring-counting efficiency relative error covariance matrix (3x3)
  double     fEreco1RCCCQE;                ///< ~0.4 GeV (RC cov matrix is 3x3, 3 Ereco bins: (0, fEreco1RCCCQE, fEreco2RCCCQE, oo))
  double     fEreco2RCCCQE;                ///< ~1.1 GeV (RC cov matrix is 3x3, 3 Ereco bins: (0, fEreco1RCCCQE, fEreco2RCCCQE, oo))
  double     fEffCCnonQERelErr;            ///< CCnonQE efficiency relative error
  double     fEffNCRelErr;                 ///< NC efficiency relative error
  double     fEffNuECCRelErr;              ///< nue CC efficiency relative error

  double     fEnergyScaleTwkDial;          ///< Energy scale tweaking dial
  double     fOthEffCCQETwkDial;           ///< CCQE other (all-but-ring-counting) efficiency tweaking dial
  double     fRC0EffCCQETwkDial;           ///< CCQE ring-counting efficiency tweaking dial 0 (Ereco < fEreco1RCCCQE)
  double     fRC1EffCCQETwkDial;           ///< CCQE ring-counting efficiency tweaking dial 1 (fEreco1RCCCQE <= Ereco <= fEreco2RCCCQE)
  double     fRC2EffCCQETwkDial;           ///< CCQE ring-counting efficiency tweaking dial 2 (Ereco > fEreco2RCCCQE)
  double     fEffCCnonQETwkDial;           ///< CCnonQE efficiency tweaking dial
  double     fEffNCTwkDial;                ///< NC efficiency tweaking dial
  double     fEffNuECCTwkDial;             ///< nue CC efficiency tweaking dial
  TMatrixD * fEffCCQERelErrCovarCholeskyL; ///< Cholesky decomposition (trangular matrix L) of fEffCCQERelErrCovarRC covariance matrix

  TH2D *     fMCTemplNom[kns][knm];        ///< MC template(sample,mode): Ereco-Etrue (nominal)
  TH2D *     fMCTemplSys[kns][knm];        ///< MC template(sample,mode): Ereco-Etrue (+ systematic variation)
  TH2D *     fMCTemplSysHyp[kns][knm];     ///< MC template(sample,mode): Ereco-Etrue (+ systematic variation + physics hypothesis)
  TH1D *     fMCErecoSpectrum[kns][knm];   ///< MC Ereco spectrum(sample,mode)        (+ systematic variation + physics hypothesis)
  TH1D *     fMCErecoSpectrumTot[kns];     ///< MC total Ereco spectrum (sample)      (+ systematic variation + physics hypothesis)

  TTree*     fMC1Rmu_numu;                 ///< input numu MC event sample
  double     fMCPOT_numu;                  ///< numu MC sample POT
  TTree*     fMC1Rmu_numubar;              ///< input numubar MC event sample
  double     fMCPOT_numubar;               ///< numubar MC sample POT
  TTree*     fMC1Rmu_nue;                  ///< input nue MC event sample
  double     fMCPOT_nue;                   ///< nue MC sample POT
  double     fBrEreco;                     ///< branch var for reading input event sample tree
  double     fBrEtrue;                     ///< branch var for reading input event sample tree
  int        fBrPdg;                       ///< branch var for reading input event sample tree
  int        fBrGenCode;                   ///< branch var for reading input event sample tree

  TH1D *     f1RmNuMuEt;                   ///< Etrue distribution for 1-Ring mu-like numu    - needed for beam norm err calculation
  TH1D *     f1RmNuMuBarEt;                ///< Etrue distribution for 1-Ring mu-like numubar - needed for beam norm err calculation
  TH1D *     f1RmNuEEt;                    ///< Etrue distribution for 1-Ring mu-like nue     - needed for beam norm err calculation
      
  NuFlux *   fNuFlux;                      ///< Neutrino flux analysis inputs
  NuXSec *   fGenInputs;                   ///< Neutrino cross section analysis inputs

  PhysHypothesis fPhysHypothesis;          ///< Hypothesis for observed deficit: disappearance, decay, decoherence
  NuOsc3Flavour  fNuOsc3Flavour;           ///< 3-flavour neutrino osc. probability calculations

ClassDef(SuperkMC, 3);
};

} //t2kana

#endif
