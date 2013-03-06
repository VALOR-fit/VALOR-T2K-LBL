//____________________________________________________________________________
/*!

\class    t2kana::FitterLevembergMarquardt

\brief    T2KAna fitter based on GSL Levemberg-Marquardt

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 06, 2013

*/
//____________________________________________________________________________

#ifndef _t2kana_fitter_levemberg_marquardt_h_
#define _t2kana_fitter_levemberg_marquardt_h_

#include <TString.h>
#include <TObject.h>
#include <TMatrixD.h>
#include <Math/GSLNLSMinimizer.h>
#include <Math/IFunction.h>

#include "Conventions.h"
#include "FitterABC.h"
#include "FitInputs.h"

namespace t2kana {

class FitterLevembergMarquardt: public FitterABC
{
public:
   FitterLevembergMarquardt();
   FitterLevembergMarquardt(FitInputs * inputs, Option_t * option);
   virtual ~FitterLevembergMarquardt();

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

   ROOT::Math::IBaseFunctionMultiDim * fFunc;
   ROOT::Math::GSLNLSMinimizer * fMinimizer;
};


} // t2kana

#endif

