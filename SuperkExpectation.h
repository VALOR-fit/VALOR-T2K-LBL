//____________________________________________________________________________
/*!

\class    t2kana::year1_numuosc::SuperkExpectation

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

#include "year1_numu_osc_analysis/PhysParamList.h"
#include "year1_numu_osc_analysis/SuperkObservation.h"
#include "year1_numu_osc_analysis/Conventions.h"

class TAxis;
class TH1D;
class TH2D;
class TSpline3;

using std::string;
using t2kana::year1_numuosc::conventions::SuperkSample_t;
using t2kana::year1_numuosc::conventions::SuperkMode_t;
using t2kana::year1_numuosc::conventions::kns;
using t2kana::year1_numuosc::conventions::knm;

namespace t2kana {
namespace year1_numuosc {

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

   // generate SK energy spectrum expectation for the current inputs
   bool Generate(double Er);

   // update SK energy spectrum expectation 
   void Update(PhysParamList & fp);

   // judge the agreement between the input data and the current prediction
   double Chisquare(SuperkObservation & skobs, Option_t * opt = "poisson");

   // access inputs used for building-up the expected spectrum at SK
   ND280Inputs *     GetND280Inputs     (void) { return fND280Inputs;     }
   SuperkInputs *    GetSuperkInputs    (void) { return fSuperkInputs;    }
   BeamInputs *      GetBeamInputs      (void) { return fBeamInputs;      }
   GeneratorInputs * GetGeneratorInputs (void) { return fGeneratorInputs; }

   //.........................................................................
   // nested class ND280Inputs
   //.........................................................................
   class ND280Inputs 
   {
    public:
      ND280Inputs();
     ~ND280Inputs();

      // get information needed from nd280 for predicting the SK spectrum
      double NumOfObsCCLikeEvents     (double Er); ///< observed number of CC-like events
      double NumOfEstNCBkgEvents      (double Er); ///< estimated NC background
      double FiducialMass             (void);
      double AnalyzedSamplePOT        (void);
      double ReconstructionEfficiency (double Er);
      double SelectionEfficiency      (double Et);
      double EnergyResolutionMatrix   (double Et, double Er);
      TAxis& ErBinning                (void) { return *fErBinning; }
      TAxis& EtBinning                (void) { return *fEtBinning; }

      // read inputs
      bool Read(string config_file);

    private:

      void Init    (void);
      void CleanUp (void);

      double  fPOT;             ///< number of POT
      double  fFidMass;         ///< fiducial mass
      TAxis * fErBinning;       ///< E_{reco} binning (both fixed and variable bin widths supported)
      TAxis * fEtBinning;       ///< E_{true} binning (both fixed and variable bin widths supported)
      TH1D  * fObsCCLikeEr;     ///< Observed CC-like inclusive dN/dE_{reco} 
      TH1D  * fEstNCBkgEr;      ///< Estimated NC-bkg dN/dE_{reco} 
      TH1D  * fRecoEffCC;       ///< Reconstruction efficiency = f(E_{true})
      TH1D  * fSelEffCC;        ///< Reconstruction efficiency = f(E_{reco})
      TH2D  * fResolMtrxCC;     ///< Energy resolution matrix = f(E_{true}, E_{reco})
   };

   //.........................................................................
   // nested class SuperkInputs
   //.........................................................................
   class SuperkInputs 
   {
    public:
      SuperkInputs();
     ~SuperkInputs();

      // get information needed from SK for predicting the SK spectrum
      double FiducialMass             (void);
      double AnalyzedSamplePOT        (void);
      double ReconstructionEfficiency (double Er, SuperkSample_t s, SuperkMode_t m);
      double SelectionEfficiency      (double Et, SuperkMode_t m);
      double EnergyResolutionMatrix   (double Et, double Er, SuperkSample_t s, SuperkMode_t m);
      TAxis& ErBinning                (void) { return *fErBinning; }
      TAxis& EtBinning                (void) { return *fEtBinning; }

      // read inputs
      bool Read(string config_file);

    private:

      void Init    (void);
      void CleanUp (void);

      double  fPOT;                  ///< number of POT
      double  fFidMass;              ///< fiducial mass
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
    public:
      BeamInputs();
     ~BeamInputs();

      // get information needed from beam for predicting the SK spectrum
      double FNRatio    (double Et) { return this->SuperkFlux(Et)/this->ND280Flux(Et); }
      double SuperkFlux (double Et);
      double ND280Flux  (double Et);
      TAxis& EtBinning  (void)      { return *fEtBinning; }

      // read inputs
      bool Read(string config_file);

    private:

      void Init    (void);
      void CleanUp (void);

      TAxis * fEtBinning;   // E_{true} binning (both fixed and variable bin widths supported)
      TH1D *  fSuperkFlux;  // SuperK flux estimation
      TH1D *  fND280Flux;   // ND280 flux estimation
   };

   //.........................................................................
   // nested class GeneratorInputs
   //.........................................................................
   class GeneratorInputs 
   {
    public:
      GeneratorInputs();
     ~GeneratorInputs();
  
      // get information needed from neutrino generator for predicting the SK spectrum
      double ModeFraction   (double Et, SuperkMode_t m);
      double SigmaTotWater  (double Et);
      double SigmaTotCCND280(double Et);

      // read inputs
      bool Read(string config_file);

    private:

      void Init    (void);
      void CleanUp (void);

      TSpline3 * fModeFraction[knm];
      TSpline3 * fSigmaTotWater;
      TSpline3 * fSigmaTotCCND280;
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

} // year1_numuosc
} // t2kana

#endif // _t2kana_year1_numu_analysis_sk_expectation_h_

