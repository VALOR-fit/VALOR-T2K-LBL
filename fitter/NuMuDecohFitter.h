//____________________________________________________________________________
/*!

\class    t2kana::NuMuDecohFitter

\brief    Muon neutrino decoherence fitter

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_decoh_fitter_h_
#define _t2kana_decoh_fitter_h_

#include <TMinuit.h>
#include <TString.h>

class TVirtualPad;
class TH1D;
class TH2D;

#include "Conventions.h"
#include "FitterABC.h"

namespace t2kana {

class SuperkObservation;
class SuperkExpectation;
class FitInputs;

class NuMuDecohFitter : public FitterABC
{
public:
   NuMuDecohFitter();
   NuMuDecohFitter(FitInputs * inputs, Option_t * option = "minuit");
  ~NuMuDecohFitter();

   bool      DoTheFit          (void);
   int       FitStatusCode     (void) { return fStatusCode; }
   void      Draw              (TVirtualPad * c, Option_t * option = "results", const char * outtag = 0);
   void      ActualDecohParams (double sinsq2theta23, double alpha);
   TMinuit * Minuit            (void) { return gMinuit; }
   TH2D *    GetChisqSurface   (void) { return fChisqSurface; }

private:

   void  Init                    (void);
   void  SetInputs               (FitInputs * inputs);
   void  PrepareFitter           (Option_t * option = "minuit");
   bool  DoTheFitCustom          (void);
   bool  DoTheFitMinuit          (void);
   void  BuildChisqSurface       (void);
   TH1D* MarginalizeChisqSurface (Option_t * option);

   FitInputs*  fFitInputs;        ///<
   bool        fKnowActualParams; ///<
   double      fSinsq2Theta23;    ///<
   double      fMusq;             ///<
   int         fNdof;             ///<
   int         fNFreeFitParams;   ///<
   int         fNFitBins;         ///<
   TH2D *      fChisqSurface;     ///<
   TString     fOption;           ///<
   int         fStatusCode;       ///<

};

} // t2kana

#endif

