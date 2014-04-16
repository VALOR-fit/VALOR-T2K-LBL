//____________________________________________________________________________
/*!

\class     t2kana::FitterMinuit

\brief     Minuit-based VALOR fitter.

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
           Thomas Dealtry <Thomas.Dealtry \at physics.ox.ac.uk>
           Universty of Oxford
           Lorena Escudero <loressa \at ific.uv.es>
           IFIC Valencia
           Nick Grant <n.grant \at lancaster.ac.uk>
           University of Lancaster

\created   October 20, 2010

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef _t2kana_fitter_minuit_h_
#define _t2kana_fitter_minuit_h_

#include <TMinuit.h>
#include <TString.h>

#include "Conventions.h"
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
   double      fTol;                 ///< set tolerance
   int         fStr;
};

} // t2kana

#endif

