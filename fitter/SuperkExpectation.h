//____________________________________________________________________________
/*!

\class    t2kana::SuperkExpectation

\brief    Code for predicting the SuperK reconstructed energy spectrum.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

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
#include "PhysHypothesis.h"
#include "Conventions.h"
#include "Utilities.h"

using std::string;
using t2kana::conv::kns;
using t2kana::conv::knm;
using t2kana::conv::knmnd;

namespace t2kana {

class SuperkMC;          ///< SuperK MC templates
class NuFlux;            ///< Neutrino flux info
class NuXSec;            ///< Neutrino generator / cross section info
class BANFF;             ///< BANFF fit info

class SuperkExpectation: public TObject
{
public:
   // ctors & dtor
   SuperkExpectation();
  ~SuperkExpectation();

   // Read inputs from input XML file
   bool Read(string config_file);

   // Methods to generate SuperK spectra for a physics hypothesis (oscillations,decay,decoherence) 
   // and tweak systematics (tweak the p.d.fs used for generating the SuperK spectra).
   void Generate                  (const PhysHypothesis & ph);
   void ResetSystematicVariations (void);                                 ///< reset all systematics
   void TweakSystematic           (string systematic, double twkdial=0.); ///< set systematic variation (from top level)
   void TweakAllSystematics       (bool do_tweak=true);                   ///< random variation of all systematics

   // Interface method for fitters. The fit parameter list can contain info for both the 
   // physics hypothesis and the systematic variations.
   void Generate (const FitParamList & fp);

   // Access calculated SK reconstructed energy spectrum for each sample/mode
   TH1D* CurrPredictedErecoSpectrum (t2kana::conv::SuperkSample_t s) const;
   TH1D* CurrPredictedErecoSpectrum (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;

   // Access inputs used for building-up the expected spectrum at SK
   SuperkMC *         GetSuperkMC          (void) { return fSuperkMC;          }
   NuFlux *           GetNuFlux            (void) { return fNuFlux;            }
   NuXSec *           GetNuXSec            (void) { return fNuXSec;            }
   BANFF *            GetBANFF             (void) { return fBANFF;             }

   // Print
   void PrintSystematicVariations(void) const;

private:

   // Initialise and clean-up
   void Init    (void);
   void Book    (void);
   void CleanUp (void);

   // Use the nuisance param values in the input FitParamList to tweak the pdfs.
   void SetNuisanceParams (const FitParamList & fp);

   bool               fIsConfigured;         ///< Is this object properly configured?
   SuperkMC *         fSuperkMC;             ///< SuperK MC templates
   NuFlux *           fNuFlux;               ///< Neutrino flux simulation inputs
   NuXSec *           fNuXSec;               ///< Neutrino generator / cross section inputs
   BANFF *            fBANFF;                ///< BANFF fit info

   TH1D *  fPredictErecoSpectrum[kns][knm];  ///< Ereco spectrum = f(sample,mode)
   TH1D *  fPredictErecoSpectrumTot[kns];    ///< Ereco spectrum = f(sample), sum over modes

ClassDef(SuperkExpectation, 5)
};

} // t2kana

#endif 

