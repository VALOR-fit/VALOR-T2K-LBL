//____________________________________________________________________________
/*!

\class    t2kana::FitterMinuit

\brief    T2KAna fitter based on Minuit

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_fitter_minuit_h_
#define _t2kana_fitter_minuit_h_

#include <TMinuit.h>
#include <TString.h>

//#include "Conventions.h"
#include "Formatting.h"
#include "FitterABC.h"

namespace t2kana {

class FitInputs;

class FitterMinuit: public FitterABC
{
public:
   FitterMinuit();
   FitterMinuit(FitInputs * inputs, Option_t * option = "minuit"); // option: "minuit","minuit+hesse_err","minuit+minos_err" (migrad_err by default)
   virtual ~FitterMinuit();

   virtual bool       DoTheFit      (void);
   virtual int        FitStatusCode (Option_t * option = ""); // option: "migrad_min", "migrad_err", "hesse_err", "minos_err"
   virtual double     FitTime       (void)              { return fFitTime;                    }
   virtual int        FitNTries     (void)              { return fFitTries;                   }
   virtual int        NDOF          (void)              { return fNFitBins - fNFreeFitParams; }
   virtual int        NBinsInFit    (void)              { return fNFitBins;                   }
   virtual int        NFitParams    (void)              { return fNFreeFitParams;             }
   virtual void       SetTol        (double tol)        { fTol  = tol ;                       }
   virtual void       SetStr        (int str   )        { fStr  = str ;                       }
   
protected:

   void  Init    (void);
   void  Init    (FitInputs * fit_inputs, Option_t * option = "minuit");
   bool  RunOnce (bool use_simplex);

   FitInputs*  fFitInputs;           ///< fit inputs
   TString     fOption;              ///< fit option
   int         fMigradMinStatusCode; ///< MIGRAD minimization status code
   int         fMigradErrStatusCode; ///< MIGRAD err status code
   int         fHesseErrStatusCode;  ///< HESSE err status code
   int         fMinosErrStatusCode;  ///< MINOS err status code
   double      fFitTime;             ///< time taken to do the fit
   int         fFitTries;            ///< number of tries taken to do the fit
   int         fNFitBins;            ///< number of bins in fit
   int         fNFreeFitParams;      ///< number of free params in fit
   bool        fCalcHesseErrors;     ///< calculate HESSE errors?
   bool        fCalcMinosErrors;     ///< calculate MINOS errors?
   bool        fCalcEScaleErrors;    ///< calculate MINOS errors for Energy scale paramter only
 //bool        fSimplex;             ///< 
   double      fTol;                 ///< set tolerance
   int         fStr;
};

} // t2kana

#endif

