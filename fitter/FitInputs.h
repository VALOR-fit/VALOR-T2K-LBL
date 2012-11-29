//____________________________________________________________________________
/*!

\class    t2kana::FitInputs

\brief    Bundles fitter inputs into a single TObject-derived class so that 
          they can be passed to the minuit fitting function via 
          TMinuit::Set/GetObjectFit();

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_fit_inputs_h_
#define _t2kana_fit_inputs_h_

#include <Minuit2/FCNBase.h>

#include "FitParamList.h"
#include "Conventions.h"

namespace t2kana {

class SuperkObservation;
class SuperkExpectation;

class FitInputs: public ROOT::Minuit2::FCNBase, public TObject
{
public :
  FitInputs(SuperkObservation * o, SuperkExpectation * e, FitParamList * p, string opt = "poisson", double errdef = 1, bool goodnessoffit = false);
 ~FitInputs();

  SuperkObservation & Observation(void);
  SuperkExpectation & Expectation(void);
  FitParamList &      FitParams  (void);

  // set the chisq definition 5 bins vs 73 bins
  void SetChisqDefinition (bool goodnessoffit) { fGoodnessOfFit = goodnessoffit; }

  // judge the agreement between the input data and the current prediction
  double Chisquare (Option_t * opt = "poisson", conv::SuperkSample_t s = conv::kSKsUndefined);

  // get chisquare for the best fit prediction
  double BestFitChisquare (Option_t * opt = "poisson");

  // Set the chisquare type (e.g. poisson) to be used in the fit
  void SetChisquareFitOption(Option_t * opt) { fOpt = opt; }

  // Minuit2: Minimum chi2 fit (delta chi2 = 1 for 1 sigma error)
  double Up          (void) const    { return fErrDef;   }
  double ErrorDef    (void) const    { return Up();   }
  void   SetErrorDef (double errdef) { fErrDef = errdef; }

  // Minuit2: Returns Chi square
  double operator() (const std::vector<double>& x) const;
      
private:

  double ChisquarePenaltyTerm(void);

  SuperkObservation * fObs;
  SuperkExpectation * fExp;
  FitParamList *      fFitParams;
  string              fOpt;
  double              fErrDef;
  bool                fAllBANFFFit;
  bool                fGoodnessOfFit;

ClassDef(FitInputs, 4)

};

//
// MINUIT fit functions
//

void FitFunc (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);

} // t2kana

#endif

