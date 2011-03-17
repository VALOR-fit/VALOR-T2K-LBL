//____________________________________________________________________________
/*!

\class    t2kana::FitInputs

\brief    Bundles fitter inputs into a single TObject-derived class so that 
          they can be passed to the minuit fitting function via 
          TMinuit::Set/GetObjectFit();

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_fit_inputs_h_
#define _t2kana_fit_inputs_h_

#include <TObject.h>

#include "PhysParamList.h"

namespace t2kana {

class SuperkObservation;
class SuperkExpectation;

class FitInputs: public TObject 
{
public :
  FitInputs(SuperkObservation * o, SuperkExpectation * e, PhysParamList * p);
 ~FitInputs();

  SuperkObservation & Observation(void);
  SuperkExpectation & Expectation(void);
  PhysParamList &     FitParams  (void);

  // judge the agreement between the input data and the current prediction
  double Chisquare(Option_t * opt = "poisson");
      
private:

  SuperkObservation * fObs;
  SuperkExpectation * fExp;
  PhysParamList *     fFitParams;

ClassDef(FitInputs, 2)

};

//
// MINUIT fit functions
//

void FitFunc (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);

} // t2kana

#endif

