//____________________________________________________________________________
/*!

\class    t2kana::SuperkMC

\brief    SuperK Monte Carlo templates

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>

\created  March 15, 2011

*/
//____________________________________________________________________________

#ifndef _t2kana_sk_mc_templates_h_
#define _t2kana_sk_mc_templates_h_

#include <TObject.h>
#include <TXMLEngine.h>

#include "Conventions.h"

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

  // Get information needed from SK for predicting the SK spectrum

  void    UseNeutrinoOscillationModel(bool on_off) { fHaveNuOsc   = on_off; }
  void    UseNeutrinoDecayModel      (bool on_off) { fHaveNuDecay = on_off; }
  void    UseNeutrinoDecoherenceModel(bool on_off) { fHaveNuDecoh = on_off; }

  void    UseOscillationParams(double sinsq2theta23, double dm23sq) { fCurrSinsq2Theta23 = sinsq2theta23; fCurrDm23sq     = dm23sq; }
  void    UseDecayParams      (double sinsq2theta23, double alpha ) { fCurrSinsq2Theta23 = sinsq2theta23; fCurrAlphaDecay = alpha;  }
  void    UseDecoherenceParams(double sinsq2theta23, double musq  ) { fCurrSinsq2Theta23 = sinsq2theta23; fCurrMusqDecoh  = musq;   }
      
  double  AnalyzedSamplePOT   (void) const;
  TTree * MCTree              (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);
  double  MCPOT               (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);
  double  NumOfEvents_Ereco   (double Ereco, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, bool ignore_energy_scale_twk=false);
  double  NumOfEvents_Etrue   (double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);

  void    UpdateErecoSpectra  (void);
  TH1D*   CurrErecoSpectrum   (void) const { return fMCErecoSpectrumGrandTot; }
  TH1D*   CurrErecoSpectrum   (t2kana::conv::SuperkSample_t s) const { return fMCErecoSpectrumTot[s]; }
  TH1D*   CurrErecoSpectrum   (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const { return fMCErecoSpectrum[s][m]; }

  void    UpdateEtrueSpectra  (void);
  TH1D*   Curr1RmuNuMuEtrue   (void) const { return f1RmNuMuEt;    } ///< current (allowing for applied weights) Etrue distribution for numu
  TH1D*   Curr1RmuNuMuBarEtrue(void) const { return f1RmNuMuBarEt; } ///<      -//-  for numubar
  TH1D*   Curr1RmuNuEEtrue    (void) const { return f1RmNuEEt;     } ///<      -//-  for nue
      
  double  EfficiencyRelErr    (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const; ///< rel err on efficiency
  double  EfficiencyTwkDial   (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const; ///< curr efficiency twk dial
  double  EfficiencyWeight    (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const; ///< efficiency weight
  double  EnergyScaleRelErr   (void) const; ///< energy scale err
  double  EnergyScale         (void) const; ///< energy scale

  // Draw inputs:
  // ...
  void Draw(TVirtualPad * c, Option_t * opt = "") const;

  // Systematic variations
  void TweakEnergyScale       (double twkdial) { fEnergyScaleTwkDial = twkdial; }
  void TweakEfficiencyCCQE    (double twkdial) { fEffCCQETwkDial     = twkdial; }
  void TweakEfficiencyCCnonQE (double twkdial) { fEffCCnonQETwkDial  = twkdial; }
  void TweakEfficiencyNC      (double twkdial) { fEffNCTwkDial       = twkdial; }
  void TweakEfficiencyNuECC   (double twkdial) { fEffNuECCTwkDial    = twkdial; }
      
  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
      
  //
  void TellMeAboutTheFlux       (NuFlux *      bi) { fNuFlux = bi; }
  void TellMeAboutCrossSections (NuXSec * gi) { fGenInputs  = gi; }
      

private:

  void Init           (void);
  void CleanUp        (void);
  void CreateTemplate (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);
      
  double NumOfEventsReadFromTree(double Ereco, double ErecoBinSz,
          t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);
      
  double  fCurrSinsq2Theta23;
  double  fCurrDm23sq;
  double  fCurrAlphaDecay;
  double  fCurrMusqDecoh;

  double  fPOT;                       ///< Number of POT
  double  fEnergyScaleRelErr;         ///< Energy scale relative error
  double  fEffCCQERelErr;             ///< CCQE efficiency relative error
  double  fEffCCnonQERelErr;          ///< CCnonQE efficiency relative error
  double  fEffNCRelErr;               ///< NC efficiency relative error
  double  fEffNuECCRelErr;            ///< nue CC efficiency relative error

  double  fEnergyScaleTwkDial;        ///< Energy scale tweaking dial
  double  fEffCCQETwkDial;            ///< CCQE efficiency tweaking dial
  double  fEffCCnonQETwkDial;         ///< CCnonQE efficiency tweaking dial
  double  fEffNCTwkDial;              ///< NC efficiency tweaking dial
  double  fEffNuECCTwkDial;           ///< nue CC efficiency tweaking dial

  TH2D *  fMCTemplates[kns][knm];     ///< Ereco-Etrue
  TH1D *  fMCEtrueSpectrum[kns][knm]; ///< Etrue
  TH1D *  fMCErecoSpectrum[kns][knm]; ///< Ereco
  TH1D *  fMCErecoSpectrumTot[kns];   ///< Ereco
  TH1D *  fMCErecoSpectrumGrandTot;   ///< Ereco
  TTree*  fMC1Rmu_numu;               ///<
  double  fMCPOT_numu;                ///<
  TTree*  fMC1Rmu_numubar;            ///<
  double  fMCPOT_numubar;             ///<
  TTree*  fMC1Rmu_nue;                ///<
  double  fMCPOT_nue;                 ///<
  double  fBrEreco;                   ///<
  double  fBrEtrue;                   ///<
  int     fBrPdg;                     ///<
  int     fBrGenCode;                 ///<

  TH1D *  f1RmNuMuEt;                 ///< fEtrue distribution for 1-Ring mu-like numu    - needed for beam norm err calculation
  TH1D *  f1RmNuMuBarEt;              ///< fEtrue distribution for 1-Ring mu-like numubar - needed for beam norm err calculation
  TH1D *  f1RmNuEEt;                  ///< fEtrue distribution for 1-Ring mu-like nue     - needed for beam norm err calculation
      
  NuFlux * fNuFlux;                   ///<
  NuXSec * fGenInputs;                ///<
      
  bool fHaveNuOsc;
  bool fHaveNuDecay;
  bool fHaveNuDecoh;

ClassDef(SuperkMC, 2);
};

} //t2kana

#endif
