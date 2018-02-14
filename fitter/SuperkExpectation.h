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

#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/VArray/VALORArray.h"
#include "valor/FitCore/PhysHypothesis.h"

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
  string Name(void) const;
  
  // Read inputs from input XML file and root file
  bool Read(string config_file, Option_t * opt = 0);
  bool Read(string config_file, int data_set_id, bool is_xml = true);


  // Methods to generate SuperK spectra for a physics hypothesis (oscillations,decay,decoherence)
  // and tweak systematics (tweak the p.d.fs used for generating the SuperK spectra).
  bool Generate             (const PhysHypothesis & ph, Option_t * opt = 0);
  void RandomizeSystematics (Option_t * opt = 0); ///< random variation of all systematics
  void TweakSystematic      (string syst_name, double twkdial=0., Option_t * opt = 0);
  void TweakSystematic      (conv::PhysParamId_t syst, double twkdial=0., Option_t * opt = 0);
  // Additional related methods to reset the systematics and enquire their values.
  void ResetSystematics (Option_t * opt = 0);
  void PrintSystematics (std::ostream & stream = std::cout, Option_t * opt = 0) const;
  
  PhysParamList Tweaks(Option_t * opt = 0) const;
  void CopyTweaksTo(PhysParamList & pl, Option_t * opt = 0) const;
  
  bool Generate (const FitParamList & fp, Option_t * opt = 0);
   
   // Write all predictions to the input container or to the input file
  bool Write (TObjArray* collection, const char * tag = "model", Option_t * write_opt = 0);
  bool Write (const char * filename, const char * tag = "model", Option_t * write_opt = 0, Option_t * file_opt = 0);
   
  // Get penalty term for current systematics
  double GetFitStatisticPenalty(Option_t * opt = 0) const;
   
  // ***T2K Legacy interface ***
  
  void ResetSystematicVariations (void) {ResetSystematics();}    ///< reset all systematics
  
  void TweakAllSystematics                       (void) {RandomizeSystematics();}
  void TweakAllBANFFSystematics                  (void);        ///< random variation of all BANFF systematics
  void TweakAllUncorrelatedXSecSystematics       (void){};        ///< random variation of all Uncorrelated XSec systematics
  void TweakAllSKDetectorFSISystematics          (void);        ///< random variation of all SK Detector systematics
  
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
  TH1D*   CurrPredictedQ2Spectrum      (valor::conv::SuperkSample_t s) const; // Returns new object, mut be deleted
  TH1D*   CurrPredictedQ2Spectrum      (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const; // Returns new object, mut be deleted

  // Double array getters
  double* CurrPredictedErecoSpectrumArray     (valor::conv::SuperkSample_t s);
  double* CurrPredictedErecoSpectrumArrayMode (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m);
  const double* CurrPredictedErecoSpectrumArray     (valor::conv::SuperkSample_t s) const;
  const double* CurrPredictedErecoSpectrumArrayMode (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const;
  double  CurrPredictedErecoSpectrumArray           (valor::conv::SuperkSample_t s, int i) const;


  // NEW 2D TEMPLATES

  // Access calculated SK reconstructed energy Vs Q2 spectrum for each sample/mode
  const VALORArray * CurrPredictedErecoVsQ2VArray     (valor::conv::SuperkSample_t s) const;
  const VALORArray * CurrPredictedErecoVsQ2VArrayMode (valor::conv::SuperkSample_t s) const;
  VALORArray       * CurrPredictedErecoVsQ2VArray     (valor::conv::SuperkSample_t s);
  VALORArray       * CurrPredictedErecoVsQ2VArrayMode (valor::conv::SuperkSample_t s);
  // As histograms
  TH2D*   CurrPredictedErecoVsQ2Spectrum      (valor::conv::SuperkSample_t s) const; // Returns new object, mut be deleted
  TH2D*   CurrPredictedErecoVsQ2Spectrum      (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const; // Returns new object, mut be deleted
  // Double array getters
  double* CurrPredictedErecoVsQ2SpectrumArray     (valor::conv::SuperkSample_t s);
  double* CurrPredictedErecoVsQ2SpectrumArrayMode (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m);
  const double* CurrPredictedErecoVsQ2SpectrumArray     (valor::conv::SuperkSample_t s) const;
  const double* CurrPredictedErecoVsQ2SpectrumArrayMode (valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const;
  double  CurrPredictedErecoVsQ2SpectrumArray           (valor::conv::SuperkSample_t s, int ir, int iq) const;
  //   







  // Access inputs used for building-up the expected spectrum at SK
  SuperkDetectorFSI * GetSuperkDetectorFSI (void) { return fSuperkDetectorFSI; }
  SuperkMC          * GetSuperkMC          (void) { return fSuperkMC;          }
  NuFlux            * GetNuFlux            (void) { return fNuFlux;            }
  NuXSec            * GetNuXSec            (void) { return fNuXSec;            }
  BANFF             * GetBANFF             (void) { return fBANFF;             }
  // Const of above
  const SuperkDetectorFSI * GetSuperkDetectorFSI (void) const { return fSuperkDetectorFSI; }
  const SuperkMC          * GetSuperkMC          (void) const { return fSuperkMC;          }
  const NuFlux            * GetNuFlux            (void) const { return fNuFlux;            }
  const NuXSec            * GetNuXSec            (void) const { return fNuXSec;            }
  const BANFF             * GetBANFF             (void) const { return fBANFF;             }

  // Print
  void PrintSystematicVariations(std::ostream & stream = std::cout) const;
  
private:
  
  // Initialise and clean-up
  void Init    (void);
  void Book    (void);
  void CleanUp (void);

  void CheckPhysParamId(void);

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

  // NEW 2D TEMPLATES
  VALORArray * fPredictErecoVsQ2SpectrumArrayTot [conv::kns];
  VALORArray * fPredictErecoVsQ2SpectrumArray    [conv::kns];
  //
  
  
  // DEPRECATED INTERFACE, DO NOT USE
public:
  void TweakAllRHCSKDetectorFSISystematics       (void);
  void TweakAllFHCpRHCSKDetectorFSISystematics   (void);
  void TweakAllMultiRingSKDetectorFSISystematics (void);
  void TweakAllSampleSKDetectorFSISystematics    (void);
  void TweakAllRHCBANFFSystematics               (void);
  void TweakAllFHCpRHCBANFFSystematics           (void);
  void TweakAllMultiRingBANFFSystematics         (void);
  void TweakAllSampleBANFFSystematics            (void);
  
ClassDef(SuperkExpectation, 6)
};

} // valor

#endif 

