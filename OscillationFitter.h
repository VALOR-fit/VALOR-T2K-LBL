//____________________________________________________________________________
/*!

\class    t2kana::OscillationFitter

\brief    Muon neutrino disapperance fitter

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_disapp_fitter_h_
#define _t2kana_disapp_fitter_h_

#include <TMinuit.h>
#include <TString.h>

class TVirtualPad;

#include "Conventions.h"
#include "PhysParamList.h"

namespace t2kana {

class SuperkObservation;
class SuperkExpectation;
class FitInputs;

class OscillationFitter
{
public:
   OscillationFitter();
   OscillationFitter(FitInputs * inputs, Option_t * option = "minuit");
  ~OscillationFitter();

   bool      DoTheFit        (void);
   void      Draw            (TVirtualPad * c, Option_t * option = "results", const char * outtag = 0);
   void      ActualOscParams (double sinsq2theta23, double dmsq23);
   TMinuit * Minuit          (void) { return gMinuit; }
   TH2D *    GetChisqSurface (void) { return fChisqSurface; }

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
   double               fDm23sq;           ///<
   int                  fNdof;             ///<
   int                  fNFreeFitParams;   ///<
   int                  fNFitBins;         ///<
   TH2D *               fChisqSurface;
   TString              fOption;
};

} // t2kana

#endif

