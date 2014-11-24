//____________________________________________________________________________
/*!

\class    valor::SuperkToyArray

\brief    Code for using many toy predictions of the Superk Ereco spectrum to calculate a likelihood.

\author   Steve Dennis <stephen.dennis \at warwick.ac.uk>
          University of Warwick, Rutherford Appleton Laboratory

\created  November 2014

\copyright (c) 2010-2015, VaLOR group.

*/
//____________________________________________________________________________

#ifndef _valor_sk_toyarray_h_
#define _valor_sk_toyarray_h_

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

#include "AnalysisOptions.h"
#include "Enumerations.h"
#include "Conventions.h"
#include "Utilities.h"
#include "SuperkExpectation.h"

#include "FitCore/FitModelABC.h"
#include "FitCore/FitParamList.h"
#include "FitCore/PhysParamList.h"
#include "Utilities/VALORArray.h"
#include "FitCore/PhysHypothesis.h"

using std::string;

namespace valor {

class SuperkToyArray: public FitModelABC
{
public:
  // ctors & dtor
  SuperkToyArray(SuperkExpectation * master_skexp);
  ~SuperkToyArray();

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
  
  // Get full mode/etrue/ereco array.
  const VALORArray * CurrPredictedVArrayMode (int itoy, valor::conv::SuperkSample_t s) const;
  VALORArray *       CurrPredictedVArrayMode (int itoy, valor::conv::SuperkSample_t s);

  // Access calculated SK reconstructed energy spectrum for each sample/mode
  const VALORArray * CurrPredictedErecoVArray     (int itoy, valor::conv::SuperkSample_t s) const;
  const VALORArray * CurrPredictedErecoVArrayMode (int itoy, valor::conv::SuperkSample_t s) const;
  VALORArray       * CurrPredictedErecoVArray     (int itoy, valor::conv::SuperkSample_t s);
  VALORArray       * CurrPredictedErecoVArrayMode (int itoy, valor::conv::SuperkSample_t s);
  
  // As histograms
  TH1D*   CurrPredictedErecoSpectrum      (int itoy, valor::conv::SuperkSample_t s) const; // Returns new object, mut be deleted
  TH1D*   CurrPredictedErecoSpectrum      (int itoy, valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const; // Returns new object, mut be deleted
  
  // Double array getters
  double* CurrPredictedErecoSpectrumArray     (int itoy, valor::conv::SuperkSample_t s);
  double* CurrPredictedErecoSpectrumArrayMode (int itoy, valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m);
  
  const double* CurrPredictedErecoSpectrumArray     (int itoy, valor::conv::SuperkSample_t s) const;
  const double* CurrPredictedErecoSpectrumArrayMode (int itoy, valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const;

  // Print
  void PrintSystematicVariations(std::ostream & stream = std::cout) const;
  
  int GetNToys(void) const;
  
  SuperkExpectation * GetMasterExpectation(void);
  const SuperkExpectation * GetMasterExpectation(void) const;
  
  void UpdateErecoSpectra(void);
  
private:

  // Initialise and clean-up
  void Init    (SuperkExpectation * master_skexp);
  void Book    (void);
  void CleanUp (void);
  
  SuperkExpectation * fMasterExpectation;
  
  VALORArray * fPredictErecoSpectrumArrayTot   [NUMOFTOYSTOMARGINALIZE][conv::kns];
  VALORArray * fPredictErecoSpectrumArrayMode  [NUMOFTOYSTOMARGINALIZE][conv::kns];
  VALORArray * fPredictErecoSpectrumArrayUnosc [NUMOFTOYSTOMARGINALIZE][conv::kns];
  
  VALORArray * fMasterSpectrumUnoscillated     [conv::kns];
  VALORArray * fMasterErecoSpectrumUnoscillated[conv::kns];
  
  PhysParamList * fToyPhysParams[NUMOFTOYSTOMARGINALIZE];
  
ClassDef(SuperkToyArray, 1)
};

} // valor

#endif 

