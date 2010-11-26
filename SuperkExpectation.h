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
    friend class SuperkExpectation;

    public:
      ND280Inputs();
     ~ND280Inputs();

      // Get nd280 oscillation analysis inputs
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
      void TweakEfficiency   (t2kana::yr1numu::conv::ND280Mode_t mode, double twkdial);
      void TweakCrossSection (t2kana::yr1numu::conv::ND280Mode_t mode, double twkdial);

    private:

      void   Init    (void);
      void   CleanUp (void);
      bool   Read    (TXMLEngine* xml, XMLNodePointer_t root_node); ///< read inputs from XML file

      double  fPOTData;              ///< number of POT (data)
      double  fNumObsCCLike;         ///< Observed number of CC-like events 
      double  fPOTMC;                ///< number of POT (MC)
      TH1F *  fNumSelEt[knmnd];      ///< number of events passing cuts as function of Etrue for various modes
      TH1F *  fSelEffEt[knmnd];      ///< selection efficiency as function of Etrue for various modes
      double  fSelEffTwkDial[knmnd]; ///< selection efficiency tweaking dials for various modes
      double  fSelEffRelErr;         ///< relative error on selection effciencies
   };

   //.........................................................................
   // nested class SuperkInputs
   //.........................................................................
   class SuperkInputs 
   {
    friend class SuperkExpectation;

    public:
      SuperkInputs();
     ~SuperkInputs();

      // Get information needed from SK for predicting the SK spectrum
      double       EffectiveTargetMassNumber(void) const;
      double       FiducialMass             (void) const;
      double       AnalyzedSamplePOT        (void) const;
      double       RecoEfficiency           (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double       SelectionEfficiency      (double Ereco, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double       EnergyTransferMatrix     (double Etrue, double Ereco, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double       FiducialMassRelErr       (void) const; ///< relative err on fiducial mass 
      double       RecoEfficiencyRelErr     (void) const; ///< relative err on "reconstruction efficiency" (reduction, OD cut, visible energy cut)
      double       EnergyScaleRelErr        (void) const; ///< energy scale err
      double       EnergyScale              (void) const; ///< energy scale 
      const TAxis& ErecoBinning             (void) const { return *fErBinning; }
      const TAxis& EtrueBinning             (void) const { return *fEtBinning; }

      // Draw inputs:
      // "ReconstructionEfficiency"            :
      // "SelectionEfficiency/OneMuRingTight"  :
      // "EnergyTransferMatrix/OneMuRingTight" :
      void Draw(TVirtualPad * c, Option_t * opt = "") const;

      // Systematic variations
      void TweakFiducialMass   (double twkdial) { fFidMassTwkDial     = twkdial; }
      void TweakRecoEfficiency (double twkdial) { fRecoEffTwkDial     = twkdial; }
      void TweakEnergyScale    (double twkdial) { fEnergyScaleTwkDial = twkdial; }
      void TweakEfficiency   (void) { } ///TODO

    private:

      void Init               (void);
      void CleanUp            (void);
      bool Read               (TXMLEngine* xml, XMLNodePointer_t root_node); ///< read inputs from XML file
      void NormTransferMatrix (t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m);

      double  fAeff;                 ///< Effective mass number
      double  fFidMass;              ///< Fiducial mass
      double  fFidMassRelErr;        ///< Fiducial mass relative error
      double  fFidMassTwkDial;       ///< Fiducial mass tweaking dial
      double  fPOT;                  ///< Number of POT
      TAxis * fErBinning;            ///< E_{reco} binning (both fixed and variable bin widths supported)
      TAxis * fEtBinning;            ///< E_{true} binning (both fixed and variable bin widths supported)
      TH1D  * fRecoEff[knm];         ///< Reconstruction efficiency = f(E_{true}, mode)
      double  fRecoEffRelErr;        ///< Reconstruction efficiency relative error
      double  fRecoEffTwkDial;       ///< Reconstruction efficiency tweaking dial
      TH1D  * fSelEff[kns][knm];     ///< Selection efficiency = f(E_{reco}, sample, mode)
      TH2D  * fETransMtrx[kns][knm]; ///< Energy transfer matrix = f(E_{true}, E_{reco}, sample, mode)
      double  fEnergyScaleRelErr;    ///< Energy scale relative error
      double  fEnergyScaleTwkDial;   ///< Energy scale tweaking dial
   };

   //.........................................................................
   // nested class BeamInputs
   //.........................................................................
   class BeamInputs 
   {
    friend class SuperkExpectation;

    public:
      BeamInputs();
     ~BeamInputs();

      // Get beam oscillation analysis inputs
      double              SuperkNuMuFlux        (double Etrue, bool wght_to_10d=true) const;
      double              SuperkNuMuFluxPOT     (void) const; 
      const TMatrixDSym & SuperkNuMuFluxCov     (void) const; ///< Superk flux covariance matrix
      double              SuperkNuMuFluxNormErr (double sinsq2theta23, double dmsq23) const; ///< normalization err
      double              SuperkNuMuFluxNorm    (double sinsq2theta23, double dmsq23) const; ///< normalization 

      // Drawing option:
      // - "SKFlux/10a/NuMu"              : nominal SuperK numu flux (10a)
      // - "SKFlux/10d/NuMu"              : weighted SuperK numu flux (10d)
      // - "SKFlux/Both10a10d/NuMu"       : nominal (10a) and weighted SuperK numu flux (10d)
      // - "SKFlux/Weights/10aTo10d/NuMu" : 10a->10d numu flux weights
      // - "SKFlux/CovarianceMatrix/NuMu" : SK numu flux covariance matrix
      void Draw(TVirtualPad * c, Option_t * opt = "SKFlux/10d") const;

      //
      void TweakNorm  (double twkdial) { fSuperkNuMuFluxNormTwkDial = twkdial; }
      void TweakShape (bool reset=false);

    private:

      void Init    (void);
      void CleanUp (void);
      bool Read    (TXMLEngine* xml, XMLNodePointer_t root_node); ///< read inputs from XML file

      TH1D *        fSuperkNuMuFlux;            ///< SuperK numu flux estimation
      double        fSuperkNuMuFluxPOT;         ///< integrated POT for numu flux simulation
      TMatrixDSym * fSuperkNuMuFluxCov;         ///< SuperK numu flux estimation
      TH1D *        fSuperkNuMuWght10aTo10d;    ///< SuperK flux weights (10a -> 10d)0
      double        fSuperkNuMuFluxNormTwkDial; ///<
      bool          fSuperkNuMuFluxHaveCholeskyL;
      TMatrixD      fSuperkNuMuFluxCovCholeskyL;
      TH1D *        fSuperkNuMuFluxTwks;
   };

   //.........................................................................
   // nested class GeneratorInputs
   //.........................................................................
   class GeneratorInputs 
   {
    friend class SuperkExpectation;

    public:
      GeneratorInputs();
     ~GeneratorInputs();
  
      // Get neutrino generator oscillation analysis inputs
      double SigmaModeSK      (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const; ///< cross section for input mode
      double SigmaModeSKRelErr(double Etrue, t2kana::yr1numu::conv::SuperkMode_t m, int sign) const; /// relative err on cross section for input mode

      // Draw inputs:
      // - "SKSigTotal" : draw the numu+H20 total cross section
      // - "SKSigMode"  : draw cross section = f(Etrue) for various event modes in numu+H20
      // - "SKSigAll"   : numu+H20 total cross section + cross section for various event modes
      // - "SKFracMode" : draw fraction = f(Etrue) for various event modes in numu+H20
      void Draw(TVirtualPad * c, Option_t * opt = "") const;

      // Systematic variations
      void TweakSigmaModeSK(t2kana::yr1numu::conv::SuperkMode_t m, double twk_dial);

    private:

      void Init    (void);
      void CleanUp (void);
      bool Read    (TXMLEngine* xml, XMLNodePointer_t root_node); ///< read inputs from XML file

      double SigmaTotSK  (double Etrue) const;
      double FractModeSK (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const;

      TH1D *             fSigmaTotSK;              ///< numu+H20 total cross section as a function of energy
      TH1D *             fFractModeSK[knm];        ///< fraction of event modes in numu+H20 as a function of energy
      util::TH1DErrEnv * fSigmaModeSKRelErr[knm];  /// 
      double             fSigmaModeSKTwkDial[knm]; ///<
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

