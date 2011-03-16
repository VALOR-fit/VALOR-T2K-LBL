//____________________________________________________________________________
/*!

\class    t2kana::SuperkExpectation

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

#ifndef _t2kana_sk_expectation_h_
#define _t2kana_sk_expectation_h_

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

#include "PhysParamList.h"
#include "Conventions.h"
#include "Utilities.h"

using std::string;
using std::vector;
using t2kana::conv::kns;
using t2kana::conv::knm;
using t2kana::conv::knmnd;

namespace t2kana {

class ND280Measurement;  ///< ND280 measurement
class SuperkMC;          ///< SuperK MC templates
class NuFlux;            ///< Neutrino flux info
class NuXSec;            ///< Neutrino generator / cross section info

class SuperkExpectation: public TObject
{
public:
   // ctors & dtor
   SuperkExpectation();
   SuperkExpectation(const vector<t2kana::conv::SuperkSample_t> & samples);
  ~SuperkExpectation();

   // Read inputs from input XML file
   bool Read(string config_file);

   // Generate SK energy spectrum expectation for the current inputs
   void Generate(const PhysParamList & fp);

   // Access calculated SK energy spectrum 
   TH1D* CurrPredictedErecoSpectrum (void) const;
   TH1D* CurrPredictedErecoSpectrum (t2kana::conv::SuperkSample_t s) const;
   TH1D* CurrPredictedErecoSpectrum (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;

   // Access inputs used for building-up the expected spectrum at SK
   ND280Measurement * GetND280Measurement  (void) { return fND280Measurement;  }
   SuperkMC *         GetSuperkMC          (void) { return fSuperkMC;          }
   NuFlux *           GetNuFlux            (void) { return fNuFlux;            }
   NuXSec *           GetNuXSec            (void) { return fNuXSec;            }

   // Reset all systematic variations
   void ResetSystematicVariations (void);
   void TweakSystematic           (string systematic, double twkdial=0.);
   void TweakAllSystematics       (bool do_tweak=true);

private:

   void Init    (void);
   void CleanUp (void);

   bool               fIsConfigured;      ///< properly configured?
   ND280Measurement * fND280Measurement;  ///< ND280 inputs
   SuperkMC *         fSuperkMC;          ///< SuperK MC templates
   NuFlux *           fNuFlux;            ///< Neutrino flux simulation inputs
   NuXSec *           fNuXSec;            ///< Neutrino generator / cross section inputs

   TH1D *  fPredictErecoSpectrum[kns][knm];  ///< Ereco
   TH1D *  fPredictErecoSpectrumTot[kns];    ///< Ereco
   TH1D *  fPredictErecoSpectrumGrandTot;    ///< Ereco

   vector<t2kana::conv::SuperkSample_t> fSamples;

ClassDef(SuperkExpectation, 2)
};

} // t2kana

#endif 

