//____________________________________________________________________________
/*!

\class    t2kana::FitterKernel

\brief    Fitter kernel

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_fitter_kernel_h_
#define _t2kana_fitter_kernel_h_

#include <TMinuit.h>
#include <TString.h>

#include "Conventions.h"
#include "FitterABC.h"

namespace t2kana {

class FitInputs;

class FitterKernel: public FitterABC
{
public:
   virtual ~FitterKernel();

   virtual bool DoTheFit      (void);
   virtual int  FitStatusCode (void) { return fStatusCode;     }
   virtual int  NBinsInFit    (void) { return fNFitBins;       }
   virtual int  NFitParams    (void) { return fNFreeFitParams; }
   virtual void Write         (Option_t *, const char *) { }

protected:
   FitterKernel();
   FitterKernel(FitInputs * inputs, Option_t * option = "minuit");

   void  Init    (void);
   void  Init    (FitInputs * fit_inputs, Option_t * option = "minuit");
   bool  RunOnce (void);

   FitInputs*  fFitInputs;        ///< fit inputs
   TString     fOption;           ///< fit option
   int         fStatusCode;       ///< fit status code
   int         fNFitBins;         ///< number of bins in fit
   int         fNFreeFitParams;   ///< number of free params in fit
   int         fNdof;             ///< number of degrees of freedom
};

} // t2kana

#endif

