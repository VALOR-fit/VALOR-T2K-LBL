//____________________________________________________________________________
/*!

\class    valor::SuperkExpectation

\brief    Code for predicting the SuperK reconstructed energy spectrum.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________

#ifndef _valor_sk_expectation_h_
#define _valor_sk_expectation_h_

#include <string>

#include <TObject.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpline.h>
#include <TVirtualPad.h>
#include <TMatrixTSym.h>
#include <TVectorD.h>
#include <TXMLEngine.h>

#include "Enumerations.h"
#include "Conventions.h"
#include "Utilities.h"

#include "FitCore/FitModelABC.h"
#include "FitCore/FitParamList.h"
#include "FitCore/PhysParamList.h"
#include "Utilities/VALORArray.h"
#include "FitCore/PhysHypothesis.h"

using std::string;

namespace valor {

class SuperkDetectorFSI; ///< SuperK detector systematics + FSI
class SuperkMC;          ///< SuperK MC templates
class NuFlux;            ///< Neutrino flux info
class NuXSec;            ///< Neutrino generator / cross-section info
class BANFF;             ///< BANFF near-detector flux and cross-section fit info

class SuperkExpectation: public FitModelABC
{
public:
   // ctors & dtor
   SuperkExpectation();
  ~SuperkExpectation();

   // *** FitModelABC Interface ***
   string Name(void);
   
   // Read inputs from input XML file
   bool Read(string config_file, Option_t * opt = 0);

   // Methods to generate SuperK spectra for a physics hypothesis (oscillations,decay,decoherence) 
   // and tweak systematics (tweak the p.d.fs used for generating the SuperK spectra).
   bool Generate             (const PhysHypothesis & ph, Option_t * opt = 0);
   void RandomizeSystematics (void); ///< random variation of all systematics
   void TweakSystematic      (string syst_name, double twkdial=0., Option_t * opt = 0);
   void TweakSystematic      (conv::PhysParamId_t syst, double twkdial=0., Option_t * opt = 0);
   // Additional related methods to reset the systematics and enquire their values.
   void ResetSystematics (void);
   void PrintSystematics (std::ostream & stream = std::cout) const;
   
   PhysParamList Tweaks(void);
   
   bool Generate (const FitParamList & fp, Option_t * opt = 0);
   
   // Write all predictions to the input container or to the input file
  bool Write (TObjArray* collection, const char * tag = "model", Option_t * write_opt = 0);
  bool Write (const char * filename, const char * tag = "model", Option_t * write_opt = 0, Option_t * file_opt = 0);
   
  // Get penalty term for current systematics
  double GetFitStatisticPenalty(Option_t * opt = 0) const;
   
   // ***T2K Legacy interface ***
   
   void ResetSystematicVariations (void) {ResetSystematics();}    ///< reset all systematics
   
   void TweakSystematicLegacy    (string systematic, double twkdial=0., bool best_banff=true); ///< set systematic variation (from top level)
   //                     best_banff refers to whether to use bestfit or param=1 when fitting a single banff parameter

   void TweakAllSystematics                       (void) {RandomizeSystematics();}

   void TweakAllBANFFSystematics                  (void);        ///< random variation of all BANFF systematics
   void TweakAllRHCBANFFSystematics               (void);        ///< random variation of all BANFF systematics
   void TweakAllSampleBANFFSystematics            (void);        ///< random variation of all BANFF systematics

   void TweakAllUncorrelatedXSecSystematics       (void);        ///< random variation of all Uncorrelated XSec systematics

   void TweakAllSKDetectorFSISystematics          (void);        ///< random variation of all SK Detector systematics
   void TweakAllRHCSKDetectorFSISystematics       (void);        ///< random variation of all SK Detector systematics   // NEW FLUX
   void TweakAllMultiRingSKDetectorFSISystematics (void);        ///< random variation of all SK Detector systematics   /// NEW MULTIRING 2
   void TweakAllSampleSKDetectorFSISystematics    (void);        ///< random variation of all SK Detector systematics   /// NEW MULTIRING 2
   

   // Interface method for fitters. The fit parameter list can contain info for both the 
   // physics hypothesis and the systematic variations.
   bool Generate (const FitParamList & fp, bool all_banff_fit, bool all_1Re_eff);

   // Access calculated SK reconstructed energy spectrum for each sample/mode
   const VALORArray * CurrPredictedErecoVArray     (valor::conv::SuperkSample_t s) const;
   const VALORArray * CurrPredictedErecoVArrayMode (valor::conv::SuperkSample_t s) const;
   VALORArray       * CurrPredictedErecoVArray     (valor::conv::SuperkSample_t s);
   VALORArray       * CurrPredictedErecoVArrayMode (valor::conv::SuperkSample_t s);
   
   // As histograms
   TH1D*   CurrPredictedErecoSpectrum      (valor::conv::SuperkSample_t s) const; // Returns new object, mut be deleted
   TH1D*   CurrPredictedErecoSpectrum      (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const; // Returns new object, mut be deleted
   
   // Double array getters
   double* CurrPredictedErecoSpectrumArray     (valor::conv::SuperkSample_t s);
   double* CurrPredictedErecoSpectrumArrayMode (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m);
   
   const double* CurrPredictedErecoSpectrumArray     (valor::conv::SuperkSample_t s) const;
   const double* CurrPredictedErecoSpectrumArrayMode (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const;
   double  CurrPredictedErecoSpectrumArray           (valor::conv::SuperkSample_t s, int i) const;
   
   // Access inputs used for building-up the expected spectrum at SK
   SuperkDetectorFSI * GetSuperkDetectorFSI (void) { return fSuperkDetectorFSI; }
   SuperkMC          * GetSuperkMC          (void) { return fSuperkMC;          }
   NuFlux            * GetNuFlux            (void) { return fNuFlux;            }
   NuXSec            * GetNuXSec            (void) { return fNuXSec;            }
   BANFF             * GetBANFF             (void) { return fBANFF;             }

   // Print
   void PrintSystematicVariations(std::ostream & stream = std::cout) const;
   
private:

   // Initialise and clean-up
   void Init    (void);
   void Book    (void);
   void CleanUp (void);

   // Use the nuisance param values in the input FitParamList to tweak the pdfs.
   void SetNuisanceParams (const FitParamList & fp, bool all_banff_fit, bool all_1Re_eff);

   bool             fIsConfigured;                ///< Is this object properly configured?
   SuperkDetectorFSI * fSuperkDetectorFSI;        ///< SuperK detector info + FSI
   SuperkMC          * fSuperkMC;                 ///< SuperK MC templates
   NuFlux            * fNuFlux;                   ///< Neutrino flux simulation inputs
   NuXSec            * fNuXSec;                   ///< Neutrino generator / cross section inputs
   BANFF             * fBANFF;                    ///< BANFF fit info
   
   VALORArray * fPredictErecoSpectrumArrayTot [conv::kns];
   VALORArray * fPredictErecoSpectrumArray    [conv::kns];
   
ClassDef(SuperkExpectation, 6)
};

} // valor

#endif 

