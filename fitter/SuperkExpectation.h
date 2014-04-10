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

#include "FitParamList.h"
#include "PhysParamList.h"
#include "PhysHypothesis.h"
#include "Conventions.h"
#include "Utilities.h"

using std::string;
using t2kana::conv::kns;
using t2kana::conv::knm;

namespace t2kana {

class SuperkDetectorFSI; ///< SuperK detector systematics + FSI
class SuperkMC;          ///< SuperK MC templates
class NuFlux;            ///< Neutrino flux info
class NuXSec;            ///< Neutrino generator / cross-section info
class BANFF;             ///< BANFF near-detector flux and cross-section fit info

class SuperkExpectation: public TObject
{
public:
   // ctors & dtor
   SuperkExpectation();
  ~SuperkExpectation();

   // Methods to generate SuperK spectra for a physics hypothesis (oscillations,decay,decoherence) 
   // and tweak systematics (tweak the p.d.fs used for generating the SuperK spectra).
   bool Generate                  (const PhysHypothesis & ph);
   void ResetSystematicVariations (void);                                 ///< reset all systematics
   void TweakSystematic           (string systematic, double twkdial=0., bool best_banff=true); ///< set systematic variation (from top level)
   //                     best_banff refers to whether to use bestfit or param=1 when fitting a single banff parameter
   void TweakAllSystematics                  (void);           ///< random variation of all systematics
   void TweakAllBANFFSystematics             (void);           ///< random variation of all BANFF systematics
   void TweakAllUncorrelatedXSecSystematics  (void);           ///< random variation of all Uncorrelated XSec systematics
   void TweakAllSKDetectorFSISystematics     (void);           ///< random variation of all SK Detector systematics

   PhysParamList Tweaks(void);

   // Interface method for fitters. The fit parameter list can contain info for both the 
   // physics hypothesis and the systematic variations.
   bool Generate (const FitParamList & fp, bool all_banff_fit, bool all_1Re_eff);

   // Access calculated SK reconstructed energy spectrum for each sample/mode
   TH1D*   CurrPredictedErecoSpectrum      (t2kana::conv::SuperkSample_t s) const; // Returns new object, mut be deleted
   TH1D*   CurrPredictedErecoSpectrum      (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const; // Returns new object, mut be deleted
   
   double* CurrPredictedErecoSpectrumArray     (t2kana::conv::SuperkSample_t s);
   double* CurrPredictedErecoSpectrumArrayMode (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m);
   
   const double* CurrPredictedErecoSpectrumArray     (t2kana::conv::SuperkSample_t s) const;
   const double* CurrPredictedErecoSpectrumArrayMode (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;
   double  CurrPredictedErecoSpectrumArray           (t2kana::conv::SuperkSample_t s, int i) const;

   // Read inputs from input XML file
   bool Read(string config_file);

   // Access inputs used for building-up the expected spectrum at SK
   SuperkDetectorFSI * GetSuperkDetectorFSI (void) { return fSuperkDetectorFSI; }
   SuperkMC *       GetSuperkMC       (void) { return fSuperkMC;       }
   NuFlux *         GetNuFlux         (void) { return fNuFlux;         }
   NuXSec *         GetNuXSec         (void) { return fNuXSec;         }
   BANFF *          GetBANFF          (void) { return fBANFF;          }

   // Print
   void PrintSystematicVariations(ostream & stream = std::cout) const;

   // Write all predictions to the input container or to the input file
   bool Write(TObjArray * data_collection, const char * tag = "pdf");
   bool Write(const char * filename, const char * tag = "pdf", const char * option_file="update");

private:

   // Initialise and clean-up
   void Init    (void);
   void Book    (void);
   void CleanUp (void);

   // Use the nuisance param values in the input FitParamList to tweak the pdfs.
   void SetNuisanceParams (const FitParamList & fp, bool all_banff_fit, bool all_1Re_eff);

   bool             fIsConfigured;                ///< Is this object properly configured?
   SuperkDetectorFSI * fSuperkDetectorFSI;              ///< SuperK detector info + FSI
   SuperkMC *       fSuperkMC;                    ///< SuperK MC templates
   NuFlux *         fNuFlux;                      ///< Neutrino flux simulation inputs
   NuXSec *         fNuXSec;                      ///< Neutrino generator / cross section inputs
   BANFF *          fBANFF;                       ///< BANFF fit info
   
   double  fPredictErecoSpectrumArrayTot1Rmu   [conv::kNErecoBins1Rmu];
   double  fPredictErecoSpectrumArrayTot1Re    [conv::kNErecoBins1Re];
   double  fPredictErecoSpectrumArray1Rmu      [conv::knm][conv::kNErecoBins1Rmu];
   double  fPredictErecoSpectrumArray1Re       [conv::knm][conv::kNErecoBins1Re];
   double  fPredictRHCErecoSpectrumArrayTot1Rmu[conv::kRHCNErecoBins1Rmu];
   double  fPredictRHCErecoSpectrumArrayTot1Re [conv::kRHCNErecoBins1Re];
   double  fPredictRHCErecoSpectrumArray1Rmu   [conv::knm][conv::kRHCNErecoBins1Rmu];
   double  fPredictRHCErecoSpectrumArray1Re    [conv::knm][conv::kRHCNErecoBins1Re];
   
ClassDef(SuperkExpectation, 5)
};

} // t2kana

#endif 

