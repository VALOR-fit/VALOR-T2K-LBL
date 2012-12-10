//____________________________________________________________________________
/*!

\class    t2kana::FitterKernelMinuit2

\brief    Fitter kernel for Minuit2

\author   Thomas Dealtry <thomas.dealtry \at physics.ox.ac.uk>
          STFC, Rutherford Appleton Laboratory
	  Subdepartment of Particle Physics, University of Oxford

\created  October 21, 2012

*/
//____________________________________________________________________________

#ifndef _t2kana_fitter_kernel_minuit2_h_
#define _t2kana_fitter_kernel_minuit2_h_

#include <TString.h>
#include <TObject.h>
#include <TMatrixD.h>
#include <Minuit2/MnUserCovariance.h>
#include <Minuit2/FunctionMinimum.h>

#include "Conventions.h"
#include "FitterABC.h"
#include "FitInputs.h"

namespace t2kana {

class FitterKernelMinuit2: public FitterABC
{

public:
   virtual ~FitterKernelMinuit2();

   virtual bool       DoTheFit      (void);
   virtual int        FitStatusCode (void)        { return fStatusCode;                 }
   virtual double     FitTime       (void)        { return fFitTime;                    }
   virtual int        NDOF          (void)        { return fNFitBins - fNFreeFitParams; }
   virtual int        NBinsInFit    (void)        { return fNFitBins;                   }
   virtual int        NFitParams    (void)        { return fNFreeFitParams;             }
/*
   virtual TObjArray& GetResults    (void) const  { return *fResults;                   }
   virtual void       Analyze       (Option_t*)   { }
   virtual void       Save          (const char*) { }
*/
protected:
   FitterKernelMinuit2 ();
   FitterKernelMinuit2 (FitInputs * inputs, Option_t * option = "minuit2");

   void  Init    (void);
   void  Init    (FitInputs * fit_inputs, Option_t * option = "minuit2");
   bool  RunOnce (void);
   void  ConstructCovMatrix(void);
   void  GetStatusCode(ROOT::Minuit2::FunctionMinimum bestfit);

   FitInputs*  fFitInputs;        ///< fit inputs
   TString     fOption;           ///< fit option
   int         fStatusCode;       ///< fit status code
   double      fFitTime;          ///< fit time
   int         fNFitBins;         ///< number of bins in fit
   int         fNFreeFitParams;   ///< number of free params in fit
/*
   TObjArray*  fResults;          ///<
*/
   TMatrixD*   fCovMtrx;          ///< covariance matrix of all free parameters in the fit

   ROOT::Minuit2::MnUserParameters fUserParams;   ///< store the user parameters
};


} // t2kana

#endif

