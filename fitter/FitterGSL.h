//____________________________________________________________________________
/*!

\class    t2kana::FitterGSL

\brief    T2KAna fitter based on GSL implementation of the following algoritms:

          ....................................................................
           Algorithm name                                 | GSL fitter option
          ....................................................................
	   Simulated annealing                            | SimAnneal
	   Levemberg-Marquardt                            | LevMarq
           Fletcher-Reeves conjugate gradient             | ConjugateFR
           Polak-Ribiere conjugate gradient               | ConjugatePR
           Vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) | VectorBFGS
           Steepest Descent                               | SteepestDescent
          ....................................................................

\reference
          http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 06, 2013

*/
//____________________________________________________________________________

#ifndef _t2kana_fitter_gsl_h_
#define _t2kana_fitter_gsl_h_

#include "Build.h"

#include <TString.h>
#include <TObject.h>
#include <TMatrixD.h>

#ifdef __GSL_ENABLED__
#include <Math/Minimizer.h>
#include <Math/IFunction.h>
#endif

#include "Conventions.h"
#include "FitterABC.h"
#include "FitInputs.h"

namespace t2kana {

class FitterGSL: public FitterABC
{
public:
   FitterGSL();
   FitterGSL(FitInputs * inputs, Option_t * option);
   virtual ~FitterGSL();

   virtual bool       DoTheFit      (void);
   virtual int        FitStatusCode (void)        { return fStatusCode;                 }
   virtual double     FitTime       (void)        { return fFitTime;                    }
   virtual int        NDOF          (void)        { return fNFitBins - fNFreeFitParams; }
   virtual int        NBinsInFit    (void)        { return fNFitBins;                   }
   virtual int        NFitParams    (void)        { return fNFreeFitParams;             }

//rm
   virtual int        FitMigradStatusCode (void)  { return -1;   }
   virtual int        FitHesseStatusCode  (void)  { return -1;   }

protected:

   void Init    (void);
   void Init(FitInputs * fit_inputs, Option_t * option);

   FitInputs*  fFitInputs;        ///< fit inputs
   TString     fOption;           ///< fit option
   int         fStatusCode;       ///< fit status code
   double      fFitTime;          ///< fit time
   int         fNFitBins;         ///< number of bins in fit
   int         fNFreeFitParams;   ///< number of free params in fit

#ifdef __GSL_ENABLED__
   ROOT::Math::IMultiGenFunction * fFunc;
   ROOT::Math::FitMethodFunction * fFitMethodFunc;
   ROOT::Math::Minimizer *         fMinimizer;
#endif
};


} // t2kana

#endif

