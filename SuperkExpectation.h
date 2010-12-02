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
   double Generate(double Ereco, const PhysParamList & fp);
   double Generate(double Ereco, const PhysParamList & fp, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);

   TH1D * AsTH1D(const PhysParamList & fp);
   TH1D * AsTH1D(const PhysParamList & fp, t2kana::yr1numu::conv::SuperkSample_t s);
   TH1D * AsTH1D(const PhysParamList & fp, t2kana::yr1numu::conv::SuperkMode_t m);
   TH1D * AsTH1D(const PhysParamList & fp, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);

   // Access inputs used for building-up the expected spectrum at SK
   ND280Inputs *     GetND280Inputs     (void) { return fND280Inputs;     }
   SuperkInputs *    GetSuperkInputs    (void) { return fSuperkInputs;    }
   BeamInputs *      GetBeamInputs      (void) { return fBeamInputs;      }
   GeneratorInputs * GetGeneratorInputs (void) { return fGeneratorInputs; }

   // Reset all systematic variations
   void ResetSystematicVariations(void);

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
      double NExpNCBkg        (void) const; ///< expected NC contamination of CC-like sample
      double NExpCC           (void) const; ///< expected number of CC events
      double NumOfEvents      (t2kana::yr1numu::conv::ND280Mode_t mode) const;
      double EfficiencyRelErr (double Et, t2kana::yr1numu::conv::ND280Mode_t mode, int sign) const; ///< relative uncertainty on efficiency

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
      double  AnalyzedSamplePOT (void) const;   
      TTree * MCTree            (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);
      double  MCPOT             (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);
      double  NumOfEvents       (double Ereco, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m, double sinsq2theta23, double dmsq23);
      double  EfficiencyRelErr  (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const; ///< rel err on efficiency
      double  EnergyScaleRelErr (void) const; ///< energy scale err
      double  EnergyScale       (void) const; ///< energy scale 

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

      void Init          (void);
      void CleanUp       (void);
      void FillHistogram (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);

      double NumOfEventsReadFromTree(double Ereco, double ErecoBinSz, 
          t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m, double sinsq2theta23, double dmsq23);

      double  fPOT;                     ///< Number of POT
      double  fEnergyScaleRelErr;       ///< Energy scale relative error
      double  fEnergyScaleTwkDial;      ///< Energy scale tweaking dial
      double  fEffRelErr[kns][knm];     ///< Efficiency relative error
      double  fEffTwkDial[kns][knm];    ///< Efficiency tweaking dial

      TH2D *  fMCSpectrum[kns][knm];    ///<
      TTree*  fMC1Rmu_numu;             ///<
      double  fMCPOT_numu;              ///<
      TTree*  fMC1Rmu_numubar;          ///<
      double  fMCPOT_numubar;           ///<
      TTree*  fMC1Rmu_nue;              ///<
      double  fMCPOT_nue;               ///<
      double  fBrEreco;                 ///<
      double  fBrEtrue;                 ///<
      int     fBrPdg;                   ///<
      int     fBrGenCode;               ///<

      BeamInputs * fBeamInputs;         ///<
      GeneratorInputs * fGenInputs;     ///<
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
      double              SuperkNuMuFluxNormErr         (double sinsq2theta23, double dmsq23) const; ///< normalization err
      double              SuperkNuMuFluxNorm            (double sinsq2theta23, double dmsq23) const; ///< normalization 
      double              SuperkNuMuFluxShapeWeight     (double Etrue);

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

      // Draw inputs:
      //
      void Draw(TVirtualPad * c, Option_t * opt = "") const;

      // Systematic variations
      void TweakCC1piToCCQEXSecRatio (double twkdial) { fCC1piToCCQEXSecRatioTwkDial = twkdial; } ///< tweak CC1pi/CCQE (both near & far)
      void TweakCCOthToCCQEXSecRatio (double twkdial) { fCCOthToCCQEXSecRatioTwkDial = twkdial; } ///< tweak CCoth/CCQE (both near & far)
      void TweakNCToCCQEXSecRatio    (double twkdial) { fNCToCCQEXSecRatioTwkDial    = twkdial; } ///< tweak NC/CCQE (both near & far)

      // Read inputs from XML file
      bool Read (TXMLEngine* xml, XMLNodePointer_t root_node); 

    private:

      void Init    (void);
      void CleanUp (void);

      double CC1piToCCQEXSecRatioWght(double Etrue) const;
      double CCOthToCCQEXSecRatioWght(double Etrue) const;
      double NCToCCQEXSecRatioWght   (double Etrue) const;
      double XSecRatioWght           (double Etrue, double twkdial, TH1D * hRelErr) const;

//      double SigmaModeSK (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m, bool neglect_twk=false) const; ///< cross section for input mode
//      double SigmaTotSK  (double Etrue) const;
//      double FractModeSK (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const;
//
//      TH1D * fSigmaTotSK;        ///< numu+H20 total cross section as a function of energy
//      TH1D * fFractModeSK[knm];  ///< fraction of event modes in numu+H20 as a function of energy

      double fCC1piToCCQEXSecRatioTwkDial;
      double fCCOthToCCQEXSecRatioTwkDial;
      double fNCToCCQEXSecRatioTwkDial;
      TH1D * fCCQEXSecRelErr;
      TH1D * fCC1piToCCQEXSecRatioRelErr;
      TH1D * fCCOthToCCQEXSecRatioRelErr;
      TH1D * fNCToCCQEXSecRatioRelErr;
   };

private:

   void Init    (void);
   void CleanUp (void);

   bool              fIsConfigured;
   ND280Inputs *     fND280Inputs;
   SuperkInputs *    fSuperkInputs;
   BeamInputs *      fBeamInputs;
   GeneratorInputs * fGeneratorInputs;

   vector<t2kana::yr1numu::conv::SuperkSample_t> fSamples;

ClassDef(SuperkExpectation, 1)
};

} // yr1numu
} // t2kana

#endif // _t2kana_year1_numu_analysis_sk_expectation_h_

