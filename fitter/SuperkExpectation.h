 //____________________________________________________________________________
/*!

\class    t2kana::SuperkExpectation

\brief    Code for predicting the SuperK reconstructed energy spectrum.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_sk_expectation_h_
#define _t2kana_sk_expectation_h_

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

#include "FitModelABC.h"
#include "FitParamList.h"
#include "PhysParamList.h"
#include "PhysHypothesis.h"
#include "Conventions.h"
#include "Utilities.h"
#include "VALORArray.h"

using std::string;
using t2kana::conv::kns;
using t2kana::conv::knm;

namespace t2kana {

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
   string Name(void) const {return "SuperkExpectation";};
   
   // Read inputs from input XML file
   bool Read(string config_file, Option_t * opt = 0);

   // Methods to generate SuperK spectra for a physics hypothesis (oscillations,decay,decoherence) 
   // and tweak systematics (tweak the p.d.fs used for generating the SuperK spectra).
   bool Generate             (const PhysHypothesis & ph);
   void RandomizeSystematics (void); ///< random variation of all systematics
   void TweakSystematic      (string syst_name, double twkdial=0., Option_t * opt = 0);
   void TweakSystematic      (conv::PhysParamId_t syst, double twkdial=0., Option_t * opt = 0);
   // Additional related methods to reset the systematics and enquire their values.
   void ResetSystematics (void);
   void PrintSystematics (ostream & stream = std::cout) const;
   
   PhysParamList Tweaks(void);
   
   bool Generate (const FitParamList & fp);
   
   // Write all predictions to the input container or to the input file
  bool Write (TObjArray* collection, const char * tag = "model", Option_t * write_opt = 0) const;
  bool Write (const char * filename, const char * tag = "model", Option_t * write_opt = 0, Option_t * file_opt = 0) const;
   
   // Get penalty term for current systematics
   double GetFitStatisticPenalty(Option_t * opt = 0) const;
   
   // ***T2K Legacy interface ***
   
   void ResetSystematicVariations (void) {ResetSystematics();}    ///< reset all systematics
   
   void TweakSystematicLegacy    (string systematic, double twkdial=0., bool best_banff=true); ///< set systematic variation (from top level)
   //                     best_banff refers to whether to use bestfit or param=1 when fitting a single banff parameter
   void TweakAllSystematics                  (void) {RandomizeSystematics();}
   void TweakAllBANFFSystematics             (void);           ///< random variation of all BANFF systematics
   void TweakAllUncorrelatedXSecSystematics  (void);           ///< random variation of all Uncorrelated XSec systematics
   void TweakAllSKDetectorFSISystematics     (void);           ///< random variation of all SK Detector systematics

   // Interface method for fitters. The fit parameter list can contain info for both the 
   // physics hypothesis and the systematic variations.
   bool Generate (const FitParamList & fp, bool all_banff_fit, bool all_1Re_eff);

   // Access calculated SK reconstructed energy spectrum for each sample/mode
   const VALORArray * CurrPredictedErecoVArray     (t2kana::conv::SuperkSample_t s) const;
   const VALORArray * CurrPredictedErecoVArrayMode (t2kana::conv::SuperkSample_t s) const;
   VALORArray       * CurrPredictedErecoVArray     (t2kana::conv::SuperkSample_t s);
   VALORArray       * CurrPredictedErecoVArrayMode (t2kana::conv::SuperkSample_t s);
   
   // As histograms
   TH1D*   CurrPredictedErecoSpectrum      (t2kana::conv::SuperkSample_t s) const; // Returns new object, mut be deleted
   TH1D*   CurrPredictedErecoSpectrum      (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const; // Returns new object, mut be deleted
   
   // Double array getters
   double* CurrPredictedErecoSpectrumArray     (t2kana::conv::SuperkSample_t s);
   double* CurrPredictedErecoSpectrumArrayMode (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);
   
   const double* CurrPredictedErecoSpectrumArray     (t2kana::conv::SuperkSample_t s) const;
   const double* CurrPredictedErecoSpectrumArrayMode (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;
   double  CurrPredictedErecoSpectrumArray           (t2kana::conv::SuperkSample_t s, int i) const;
   
   // Access inputs used for building-up the expected spectrum at SK
   SuperkDetectorFSI * GetSuperkDetectorFSI (void) { return fSuperkDetectorFSI; }
   SuperkMC          * GetSuperkMC          (void) { return fSuperkMC;          }
   NuFlux            * GetNuFlux            (void) { return fNuFlux;            }
   NuXSec            * GetNuXSec            (void) { return fNuXSec;            }
   BANFF             * GetBANFF             (void) { return fBANFF;             }

   // Print
   void PrintSystematicVariations(ostream & stream = std::cout) const;
   
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
   
   VALORArray * fPredictErecoSpectrumArrayTot1Rmu   ;
   VALORArray * fPredictErecoSpectrumArrayTot1Re    ;
   VALORArray * fPredictErecoSpectrumArray1Rmu      ;
   VALORArray * fPredictErecoSpectrumArray1Re       ;
   VALORArray * fPredictRHCErecoSpectrumArrayTot1Rmu;
   VALORArray * fPredictRHCErecoSpectrumArrayTot1Re ;
   VALORArray * fPredictRHCErecoSpectrumArray1Rmu   ;
   VALORArray * fPredictRHCErecoSpectrumArray1Re    ;
   
ClassDef(SuperkExpectation, 6)
};

} // t2kana

#endif 

