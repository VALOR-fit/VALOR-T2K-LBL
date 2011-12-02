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

#include <TObject.h>

#include "FitParamList.h"
#include "Conventions.h"

namespace t2kana {

class SuperkObservation;
class SuperkExpectation;

class FitInputs: public TObject 
{
public :
  FitInputs(SuperkObservation * o, SuperkExpectation * e, FitParamList * p);
 ~FitInputs();

  SuperkObservation & Observation(void);
  SuperkExpectation & Expectation(void);
  FitParamList &      FitParams  (void);

  // judge the agreement between the input data and the current prediction
  double Chisquare (Option_t * opt = "poisson", conv::SuperkSample_t s = conv::kSKsUndefined);

  // get chisquare for the best fit prediction
  double BestFitChisquare (Option_t * opt = "poisson");
      
private:

  SuperkObservation * fObs;
  SuperkExpectation * fExp;
  FitParamList *      fFitParams;

ClassDef(FitInputs, 4)

};

//
// MINUIT fit functions
//

void FitFunc (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);

} // t2kana

#endif

