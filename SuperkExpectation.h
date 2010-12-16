//____________________________________________________________________________
/*!

\class    t2kana::yr1numu::SuperkExpectation

\brief    Code for predicting the SuperK reconstructed energy spectrum.

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk> 
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_year1_numu_analysis_sk_expectation_h_
#define _t2kana_year1_numu_analysis_sk_expectation_h_

#include <string>
#include <vector>

#include <TObject.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpline.h>
#include <TVirtualPad.h>
#include <TMatrixTSym.h>
#include <TVectorD.h>
#include <TXMLEngine.h>

#include "year1_numu_osc_analysis/PhysParamList.h"
#include "year1_numu_osc_analysis/Conventions.h"
#include "year1_numu_osc_analysis/Utilities.h"

using std::string;
using std::vector;
using t2kana::yr1numu::conv::kns;
using t2kana::yr1numu::conv::knm;
using t2kana::yr1numu::conv::knmnd;

class HadErrCalc;

namespace t2kana {
namespace yr1numu {

class SuperkExpectation: public TObject
{
public:
   // forward declaration of nested utility classes
   class ND280Inputs;      ///< ND280 inputs 
   class SuperkInputs;     ///< SuperK inputs
   class BeamInputs;       ///< Beam inputs
   class GeneratorInputs;  ///< Neutrino generator (GENIE/NEUT/...) inputs

   // ctors & dtor
   SuperkExpectation();
   SuperkExpectation(const vector<t2kana::yr1numu::conv::SuperkSample_t> & samples);
  ~SuperkExpectation();

   // Read inputs from input XML file
   bool Read(string config_file);

   // Generate SK energy spectrum expectation for the current inputs
   void Generate(const PhysParamList & fp);

   // Access calculated SK energy spectrum 
   TH1D* CurrPredictedErecoSpectrum (void) const;
   TH1D* CurrPredictedErecoSpectrum (t2kana::yr1numu::conv::SuperkSample_t s) const;
   TH1D* CurrPredictedErecoSpectrum (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const;

   // Access inputs used for building-up the expected spectrum at SK
   ND280Inputs *     GetND280Inputs     (void) { return fND280Inputs;     }
   SuperkInputs *    GetSuperkInputs    (void) { return fSuperkInputs;    }
   BeamInputs *      GetBeamInputs      (void) { return fBeamInputs;      }
   GeneratorInputs * GetGeneratorInputs (void) { return fGeneratorInputs; }

   // Reset all systematic variations
   void ResetSystematicVariations (void);
   void TweakSystematic           (string systematic, double twkdial=0.);
   void TweakAllSystematics       (bool do_tweak=true);

   //.........................................................................
   // nested class ND280Inputs
   //.........................................................................
   class ND280Inputs 
   {
    public:
      ND280Inputs();
     ~ND280Inputs();

      // Get nd280 oscillation analysis inputs
      double CorrectionFactor (void) const;
      double NObsCCLike       (void) const; ///< observed number of CC-like events
      double NExpNCBkg        (bool nominal=false) const; ///< expected NC contamination of CC-like sample (by default report tweaked value)
      double NExpCC           (bool nominal=false) const; ///< expected number of CC events  (by default report tweaked value)
      double NumOfEvents      (t2kana::yr1numu::conv::ND280Mode_t mode, bool nominal=false) const; 
      double NumOfEvents      (t2kana::yr1numu::conv::ND280Mode_t mode, int bin, bool nominal=false) const;
      double NumOfEvents      (t2kana::yr1numu::conv::ND280Mode_t mode, double Etrue, bool nominal=false) const;
      double EfficiencyRelErr (double Et, t2kana::yr1numu::conv::ND280Mode_t mode, int sign) const; ///< relative uncertainty on efficiency

      void   UpdateEtrueSpectra (void);
      TH1D*  CurrNuMuEtrue      (void) { return fNuMuEt;    } ///< current (allowing for applied weights) Etrue distribution for numu
      TH1D*  CurrNuMuBarEtrue   (void) { return fNuMuBarEt; } ///<      -//-  for numubar
      TH1D*  CurrNuEEtrue       (void) { return fNuEEt;     } ///<      -//-  for nue

      // Draw inputs:
      // "Saclay/EtrueSpectrum/NuMu/All"
      // "Saclay/EtrueSpectrum/NuMu/CC"
      // "Saclay/EtrueSpectrum/NuMu/NC"
      // "Saclay/EtrueSpectrum/NuMuBar/CC"
      // "Saclay/EtrueSpectrum/NuMuBar/NC"
      // "Saclay/EtrueSpectrum/NuMuBar/CC"
      // "Saclay/EtrueSpectrum/NuMuBar/NC"
      // "Saclay/Efficiency/NuMu/All"
      // "Saclay/Efficiency/NuMu/CC"
      // "Saclay/Efficiency/NuMu/NC"
      // "Saclay/Efficiency/NuMuBar/CC"
      // "Saclay/Efficiency/NuMuBar/NC"
      // "Saclay/Efficiency/NuMuBar/CC"
      // "Saclay/Efficiency/NuMuBar/NC"
      void Draw(TVirtualPad * c, Option_t * opt = "") const;

      // Systematic variations
      void TweakEfficiency (t2kana::yr1numu::conv::ND280Mode_t mode, double twkdial);
      void TweakEfficiency (double twkdial);

      //
      void TellMeAboutCrossSections(GeneratorInputs * geninp) { fGenInputs = geninp; }

      // Read inputs from XML file
      bool Read (TXMLEngine* xml, XMLNodePointer_t root_node); 

    private:

      void   Init    (void);
      void   CleanUp (void);

      double  fPOTData;                ///< number of POT (data)
      double  fNumObsCCLike;           ///< Observed number of CC-like events 
      double  fNumExpNCBkgDeflt;       ///<
      double  fNumNExpCCDeflt;         ///<
      double  fPOTMC;                  ///< number of POT (MC)
      TH1F *  fNumSelEt[knmnd];        ///< number of events passing cuts as function of Etrue for various modes
      TH1F *  fSelEffEt[knmnd];        ///< selection efficiency as function of Etrue for various modes
      double  fSelEffTwkDial[knmnd];   ///< selection efficiency tweaking dials for various modes
      double  fSelEffRelErr;           ///< relative error on selection effciencies

      TH1D *  fNuMuEt;                 ///< fEtrue distribution for numu    - needed for beam norm err calculation
      TH1D *  fNuMuBarEt;              ///< fEtrue distribution for numubar - needed for beam norm err calculation
      TH1D *  fNuEEt;                  ///< fEtrue distribution for nue     - needed for beam norm err calculation

      double  fNomCorr;                ///< nominal correction factor (without efficiency, cross section tweaks)

      GeneratorInputs * fGenInputs;
   };

   //.........................................................................
   // nested class SuperkInputs
   //.........................................................................
   class SuperkInputs 
   {
    public:
      SuperkInputs();
     ~SuperkInputs();

      // Get information needed from SK for predicting the SK spectrum

      void    UseOscillationParams(double sinsq2theta23, double dm23sq) { fCurrSinsq2Theta23 = sinsq2theta23; fCurrDm23sq = dm23sq; }

      double  AnalyzedSamplePOT   (void) const;   
      TTree * MCTree              (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);
      double  MCPOT               (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);
      double  NumOfEvents_Ereco   (double Ereco, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m, bool ignore_energy_scale_twk=false);
      double  NumOfEvents_Etrue   (double Etrue, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);

      void    UpdateErecoSpectra  (void);
      TH1D*   CurrErecoSpectrum   (void) const { return fMCErecoSpectrumGrandTot; }
      TH1D*   CurrErecoSpectrum   (t2kana::yr1numu::conv::SuperkSample_t s) const { return fMCErecoSpectrumTot[s]; }
      TH1D*   CurrErecoSpectrum   (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const { return fMCErecoSpectrum[s][m]; }

      void    UpdateEtrueSpectra  (void);
      TH1D*   Curr1RmuNuMuEtrue   (void) const { return f1RmNuMuEt;    } ///< current (allowing for applied weights) Etrue distribution for numu
      TH1D*   Curr1RmuNuMuBarEtrue(void) const { return f1RmNuMuBarEt; } ///<      -//-  for numubar
      TH1D*   Curr1RmuNuEEtrue    (void) const { return f1RmNuEEt;     } ///<      -//-  for nue

      double  EfficiencyRelErr    (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const; ///< rel err on efficiency
      double  EnergyScaleRelErr   (void) const; ///< energy scale err
      double  EnergyScale         (void) const; ///< energy scale 

      // Draw inputs:
      // ...
      void Draw(TVirtualPad * c, Option_t * opt = "") const;

      // Systematic variations
      void TweakEnergyScale (double twkdial);
      void TweakEfficiency  (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m, double twkdial); 

      // Read inputs from XML file
      bool Read (TXMLEngine* xml, XMLNodePointer_t root_node); 

      //
      void TellMeAboutTheFlux       (BeamInputs *      bi) { fBeamInputs = bi; }
      void TellMeAboutCrossSections (GeneratorInputs * gi) { fGenInputs  = gi; }

    private:

      void Init           (void);
      void CleanUp        (void);
      void CreateTemplate (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);

      double NumOfEventsReadFromTree(double Ereco, double ErecoBinSz, 
          t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);

      double  fCurrSinsq2Theta23;
      double  fCurrDm23sq;

      double  fPOT;                       ///< Number of POT
      double  fEnergyScaleRelErr;         ///< Energy scale relative error
      double  fEnergyScaleTwkDial;        ///< Energy scale tweaking dial
      double  fEffRelErr[kns][knm];       ///< Efficiency relative error
      double  fEffTwkDial[kns][knm];      ///< Efficiency tweaking dial

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

      BeamInputs *      fBeamInputs;      ///<
      GeneratorInputs * fGenInputs;       ///<
   };

   //.........................................................................
   // nested class BeamInputs
   //.........................................................................
   class BeamInputs 
   {
    public:
      BeamInputs();
     ~BeamInputs();

      // Get beam oscillation analysis inputs
      double              SuperkNuMuFlux10aTo10dWght    (double Etrue) const;
      double              SuperkNuMuBarFlux10aTo10dWght (double Etrue) const;
      double              SuperkNuEFlux10aTo10dWght     (double Etrue) const;
      double              SuperkNuMuFlux                (double Etrue, bool wght_to_10d=true) const;
      double              SuperkNuMuFluxPOT             (void) const; 
      const TMatrixD &    SuperkNuMuFluxCov             (void) const; ///< Superk flux covariance matrix
      double              SuperkNuMuFluxShapeWeight     (double Etrue);
      double              SuperkFluxNormErr             (TH1D * hnumu1Rmu, TH1D * hnumubar1Rmu, TH1D * hnue1Rm, 
                                                         TH1D * hnumuND, TH1D * hnumubarND, TH1D * hnueND) const; ///< normalization err
      double              SuperkFluxNorm                (TH1D * hnumu1Rmu, TH1D * hnumubar1Rmu, TH1D * hnue1Rm, 
                                                         TH1D * hnumuND, TH1D * hnumubarND, TH1D * hnueND) const; ///< normalization 
      bool                UsingKuboSansCode             (void) { return (fHEC!=0); }

      // Drawing option:
      // - "SKFlux/10a/NuMu"              : nominal SuperK numu flux (10a)
      // - "SKFlux/10d/NuMu"              : weighted SuperK numu flux (10d)
      // - "SKFlux/Both10a10d/NuMu"       : nominal (10a) and weighted SuperK numu flux (10d)
      // - "SKFlux/Weights/10aTo10d/NuMu" : 10a->10d numu flux weights
      // - "SKFlux/CovarianceMatrix/NuMu" : SK numu flux covariance matrix
      void Draw(TVirtualPad * c, Option_t * opt = "SKFlux/10d") const;

      // Systematic variations
      void TweakNorm  (double twkdial) { fSuperkNuMuFluxNormTwkDial = twkdial; }
      void TweakShape (bool enable_tweaking=true);

      // Read inputs from XML file
      bool Read (TXMLEngine* xml, XMLNodePointer_t root_node); 

    private:

      void Init    (void);
      void CleanUp (void);

      TH1D *        fSuperkNuMuFlux;               ///< SuperK numu flux estimation (50 MeV binning)
      TH1D *        fSuperkNuMuFluxVarBin;         ///< SuperK numu flux estimation (beam group var binning scheme)
      double        fSuperkNuMuFluxPOT;            ///< integrated POT for numu flux simulation
      TMatrixD *    fSuperkNuMuFluxCov;            ///< SuperK numu flux estimation
      TH1D *        fSuperkNuMuWght10aTo10d;       ///< SuperK numu    flux weights (10a -> 10d)
      TH1D *        fSuperkNuMuBarWght10aTo10d;    ///< SuperK numubar flux weights (10a -> 10d)
      TH1D *        fSuperkNuEWght10aTo10d;        ///< SuperK nue     flux weights (10a -> 10d)
      double        fSuperkNuMuFluxNormTwkDial;    ///<
      bool          fSuperkNuMuFluxHaveCholeskyL;
      TMatrixD *    fSuperkNuMuFluxCovCholeskyL;
      TH1D *        fSuperkNuMuFluxTweaked;        ///< SuperK numu flux estimation (tweaked / 50 MeV binning)
      bool          fEnableFluxShapeTweaking;
      HadErrCalc *  fHEC;                          ///< Kubo-san / BEAM WG norm err calc code
      double        fSuperkNormRelErr_simple;      ///< a fixed norm rel error (to use when not using Kubo-san's code)
   };

   //.........................................................................
   // nested class GeneratorInputs
   //.........................................................................
   class GeneratorInputs 
   {
    public:
      GeneratorInputs();
     ~GeneratorInputs();
  
      // Get neutrino generator oscillation analysis inputs
      double ND280XSecRatioWght  (double Etrue, t2kana::yr1numu::conv::ND280Mode_t  m) const; 
      double SuperkXSecRatioWght (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double SuperkCCQEXSecWght  (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double SuperkNuclEffectWght(double Ereco, t2kana::yr1numu::conv::SuperkMode_t m) const;

      // Draw inputs:
      //
      void Draw(TVirtualPad * c, Option_t * opt = "") const;

      // Systematic variations
      void TweakCCQEXSec             (double twkdial) { fCCQEXSecTwkDial             = twkdial; } ///< tweak CCQE (far only)
      void TweakCC1piToCCQEXSecRatio (double twkdial) { fCC1piToCCQEXSecRatioTwkDial = twkdial; } ///< tweak CC1pi/CCQE (both near & far)
      void TweakCCOthToCCQEXSecRatio (double twkdial) { fCCOthToCCQEXSecRatioTwkDial = twkdial; } ///< tweak CCoth/CCQE (both near & far)
      void TweakNCToCCQEXSecRatio    (double twkdial) { fNCToCCQEXSecRatioTwkDial    = twkdial; } ///< tweak NC/CCQE (both near & far)
      void TweakNuclEffect           (double twkdial) { fNuclEffectTwkDial           = twkdial; } ///< tweak nuclear effect (far only)

      bool TweakAffectsND280 (void) const;
      bool TweakAffectsSuperk(void) const;

      // Read inputs from XML file
      bool Read (TXMLEngine* xml, XMLNodePointer_t root_node); 

    private:

      void Init    (void);
      void CleanUp (void);

      double CC1piToCCQEXSecRatioWght(double Etrue) const;
      double CCOthToCCQEXSecRatioWght(double Etrue) const;
      double NCToCCQEXSecRatioWght   (double Etrue) const;
      double XSecRatioWght           (double Etrue, double twkdial, TH1D * hRelErr) const;

      double fCCQEXSecTwkDial;               ///<
      double fCC1piToCCQEXSecRatioTwkDial;   ///<
      double fCCOthToCCQEXSecRatioTwkDial;   ///<
      double fNCToCCQEXSecRatioTwkDial;      ///<
      double fNuclEffectTwkDial;             ///<
      TH1D * fCCQEXSecRelErr;                ///<
      TH1D * fCC1piToCCQEXSecRatioRelErr;    ///<
      TH1D * fCCOthToCCQEXSecRatioRelErr;    ///<
      TH1D * fNCToCCQEXSecRatioRelErr;       ///<
      TH1D * fSuperkNuclearEffectRelErr;     ///<
   };

private:

   void Init    (void);
   void CleanUp (void);

   bool              fIsConfigured;          ///<
   ND280Inputs *     fND280Inputs;           ///<
   SuperkInputs *    fSuperkInputs;          ///<
   BeamInputs *      fBeamInputs;            ///<
   GeneratorInputs * fGeneratorInputs;       ///<

   TH1D *  fPredictErecoSpectrum[kns][knm];  ///< Ereco
   TH1D *  fPredictErecoSpectrumTot[kns];    ///< Ereco
   TH1D *  fPredictErecoSpectrumGrandTot;    ///< Ereco

   vector<t2kana::yr1numu::conv::SuperkSample_t> fSamples;

ClassDef(SuperkExpectation, 1)
};

} // yr1numu
} // t2kana

#endif // _t2kana_year1_numu_analysis_sk_expectation_h_

