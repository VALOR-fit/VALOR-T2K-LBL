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

   virtual bool       DoTheFit      (void);
   virtual int        FitStatusCode (void)        { return fStatusCode;                 }
   virtual int        FitMigradStatusCode (void)  { return fMigradStatusCode;           }
   virtual int        FitHesseStatusCode  (void)  { return fHesseStatusCode;            }
   virtual double     FitTime       (void)        { return fFitTime;                    }
   virtual int        NDOF          (void)        { return fNFitBins - fNFreeFitParams; }
   virtual int        NBinsInFit    (void)        { return fNFitBins;                   }
   virtual int        NFitParams    (void)        { return fNFreeFitParams;             }
   virtual TObjArray& GetResults    (void) const  { return *fResults;                   }
   virtual void       Analyze       (Option_t*)   { }
   virtual void       Save          (const char*) { }
   void               SetHesseErrors(bool evaluate_hesse) { fGetHesseErrors = evaluate_hesse; }

protected:
   FitterKernel();
   FitterKernel(FitInputs * inputs, Option_t * option = "minuit");

   void  Init    (void);
   void  Init    (FitInputs * fit_inputs, Option_t * option = "minuit");
   bool  RunOnce (bool use_simplex);

   FitInputs*  fFitInputs;        ///< fit inputs
   TString     fOption;           ///< fit option
   int         fStatusCode;       ///< fit status code
   int         fMigradStatusCode; ///< return value of the MIGRAD command
   int         fHesseStatusCode;  ///< return value of the HESSE command
   double      fFitTime;          ///< time taken to do the fit
   int         fNFitBins;         ///< number of bins in fit
   int         fNFreeFitParams;   ///< number of free params in fit
   TObjArray*  fResults;          ///<
   bool        fGetHesseErrors;   ///< get hesse errors?
};

} // t2kana

#endif

