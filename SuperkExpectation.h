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

#include <TObject.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpline.h>
#include <TVirtualPad.h>
#include <TXMLEngine.h>

#include "year1_numu_osc_analysis/PhysParamList.h"
#include "year1_numu_osc_analysis/SuperkObservation.h"
#include "year1_numu_osc_analysis/Conventions.h"

using std::string;
using t2kana::yr1numu::conv::kns;
using t2kana::yr1numu::conv::knm;

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
  ~SuperkExpectation();

   // Generate SK energy spectrum expectation for the current inputs
   double Generate(double Ereco, const PhysParamList & fp);

   // Read inputs from input XML file
   bool Read(string config_file);

   // Access inputs used for building-up the expected spectrum at SK
   ND280Inputs *     GetND280Inputs     (void) { return fND280Inputs;     }
   SuperkInputs *    GetSuperkInputs    (void) { return fSuperkInputs;    }
   BeamInputs *      GetBeamInputs      (void) { return fBeamInputs;      }
   GeneratorInputs * GetGeneratorInputs (void) { return fGeneratorInputs; }

   //.........................................................................
   // nested class ND280Inputs
   //.........................................................................
   class ND280Inputs 
   {
    friend class SuperkExpectation;

    public:
      ND280Inputs();
     ~ND280Inputs();

      // Get information needed from nd280 for predicting the SK spectrum
      double       EffectiveTargetMassNumber(void) const;
      double       FiducialMass             (void) const;
      double       AnalyzedSamplePOT        (void) const;
      double       NumOfObsCCLikeEvents     (void) const; ///< observed number of CC-like events
      double       NumOfEstNCBkgEvents      (void) const; ///< estimated NC background
      double       Efficiency               (double Etrue) const;
      const TAxis& EtrueBinning             (void) const { return *fEtBinning; }

    private:

      void Init    (void);
      void CleanUp (void);
      bool Read    (TXMLEngine* xml, XMLNodePointer_t root_node); ///< read inputs from XML file

      double  fAeff;            ///< target mass number
      double  fFidMass;         ///< fiducial mass
      double  fPOT;             ///< number of POT
      TAxis * fEtBinning;       ///< E_{true} binning (both fixed and variable bin widths supported)
      double  fNumObsCCLike;    ///< Observed number of CC-like events 
      double  fNumEstNCBkg;     ///< Estimated number of NC-bkg events
      TH1D  * fEfficiency;      ///< Reconstruction efficiency = f(E_{reco})
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
      double       ReconstructionEfficiency (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double       SelectionEfficiency      (double Ereco, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double       EnergyResolutionMatrix   (double Etrue, double Ereco, t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const;
      const TAxis& ErecoBinning             (void) const { return *fErBinning; }
      const TAxis& EtrueBinning             (void) const { return *fEtBinning; }

      // Draw inputs:
      void Draw(TVirtualPad * c, Option_t * opt = "") const;

    private:

      void Init    (void);
      void CleanUp (void);
      bool Read    (TXMLEngine* xml, XMLNodePointer_t root_node); ///< read inputs from XML file

      double  fAeff;                 ///< effective mass number
      double  fFidMass;              ///< fiducial mass
      double  fPOT;                  ///< number of POT
      TAxis * fErBinning;            ///< E_{reco} binning (both fixed and variable bin widths supported)
      TAxis * fEtBinning;            ///< E_{true} binning (both fixed and variable bin widths supported)
      TH1D  * fRecoEff[knm];         ///< Reconstruction efficiency = f(E_{true}, mode)
      TH1D  * fSelEff[kns][knm];     ///< Selection efficiency = f(E_{reco}, sample, mode)
      TH2D  * fResolMtrx[kns][knm];  ///< Energy resolution matrix = f(E_{true}, E_{reco}, sample, mode)
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

      // Get information needed from beam for predicting the SK spectrum
      double       FNRatio    (double Etrue) const;
      double       SuperkFlux (double Etrue) const;
      double       ND280Flux  (double Etrue) const;
      const TAxis& EtBinning  (void) const { return *fEtBinning; }

      // Draw inputs:
      // Set option to "SKFlux"       to draw the SK flux
      // Set option to "ND280Flux"    to draw the ND280 flux
      // Set option to "OffAxisRatio" to draw the SK/ND280
      void Draw(TVirtualPad * c, Option_t * opt = "OffAxisRatio") const;

    private:

      void Init    (void);
      void CleanUp (void);
      bool Read    (TXMLEngine* xml, XMLNodePointer_t root_node); ///< read inputs from XML file

      TAxis * fEtBinning;   // E_{true} binning (both fixed and variable bin widths supported)
      TH1D *  fSuperkFlux;  // SuperK flux estimation
      TH1D *  fND280Flux;   // ND280 flux estimation
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
  
      // get information needed from neutrino generator for predicting the SK spectrum
      double SigmaTotSK     (double Etrue) const;
      double FractModeSK    (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double SigmaModeSK    (double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const;
      double SigmaTotCCND280(double Etrue) const;

      // Draw inputs:
      // - "SKSigTotal" : draw the numu+H20 total cross section
      // - "SKSigMode"  : draw cross section = f(Etrue) for various event modes in numu+H20
      // - "SKSigAll"   : numu+H20 total cross section + cross section for various event modes
      // - "SKFracMode" : draw fraction = f(Etrue) for various event modes in numu+H20
      void Draw(TVirtualPad * c, Option_t * opt = "") const;

    private:

      void Init    (void);
      void CleanUp (void);
      bool Read    (TXMLEngine* xml, XMLNodePointer_t root_node); ///< read inputs from XML file

      TH1D * fFractModeSK[knm];
      TH1D * fSigmaTotSK;
      TH1D * fSigmaTotCCND280;
   };

private:

   void Init    (void);
   void CleanUp (void);

   PhysParamList *   fPhysParamList;
   ND280Inputs *     fND280Inputs;
   SuperkInputs *    fSuperkInputs;
   BeamInputs *      fBeamInputs;
   GeneratorInputs * fGeneratorInputs;

ClassDef(SuperkExpectation, 1)
};

} // yr1numu
} // t2kana

#endif // _t2kana_year1_numu_analysis_sk_expectation_h_

