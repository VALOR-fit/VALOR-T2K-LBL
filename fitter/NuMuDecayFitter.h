//____________________________________________________________________________
/*!

\class    t2kana::NuMuDecayFitter

\brief    Muon neutrino decay fitter

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_numu_decay_fitter_h_
#define _t2kana_numu_decay_fitter_h_

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

class NuMuDecayFitter: public FitterABC
{
public:
   NuMuDecayFitter();
   NuMuDecayFitter(FitInputs * inputs, Option_t * option = "minuit");
  ~NuMuDecayFitter();

   bool      DoTheFit          (void);
   void      Draw              (TVirtualPad * c, Option_t * option = "results", const char * outtag = 0);
   void      ActualDecayParams (double sinsq2theta23, double alpha);
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


   FitInputs*           fFitInputs;
   bool                 fKnowActualParams; ///<
   double               fSinsq2Theta23;    ///<
   double               fAlpha;            ///<
   int                  fNdof;             ///<
   int                  fNFreeFitParams;   ///<
   int                  fNFitBins;         ///<
   TH2D *               fChisqSurface;
   TString              fOption;
};

} // t2kana

#endif

