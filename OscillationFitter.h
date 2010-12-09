//____________________________________________________________________________
/*!

\class    t2kana::yr1numu::OscillationFitter

\brief    MINUIT-based numu oscillation param fitter

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_year1_numu_analysis_fitter_minuit_h_
#define _t2kana_year1_numu_analysis_fitter_minuit_h_

#include <TMinuit.h>
#include <TString.h>

class TVirtualPad;

#include "year1_numu_osc_analysis/Conventions.h"
#include "year1_numu_osc_analysis/PhysParamList.h"

namespace t2kana {
namespace yr1numu {

class SuperkObservation;
class SuperkExpectation;

class OscillationFitter
{
public:
   class FitInputs;

   OscillationFitter();
   OscillationFitter(OscillationFitter::FitInputs * inputs, Option_t * option = "minuit");
  ~OscillationFitter();

   bool      DoTheFit        (void);
   void      Draw            (TVirtualPad * c, Option_t * option = "results", const char * outtag = 0);
   void      ActualOscParams (double sinsq2theta23, double dmsq23);
   TMinuit * Minuit          (void) { return gMinuit; }
   TH2D *    GetChisqSurface (void) { return fChisqSurface; }

   //
   // nested utility classes
   //

   // Bundle fitter inputs into a single TObject-derived class so that they can be 
   // passed to the minuit fitting function via TMinuit::Set/GetObjectFit();
   class FitInputs: public TObject 
   {
    public :
      FitInputs(SuperkObservation * skobs, SuperkExpectation * skexp, PhysParamList * fpar) :
       fObs(skobs), fExp(skexp), fFitParams(fpar) { }
     ~FitInputs()  { }

      SuperkObservation & Observation(void) { return *fObs;       }
      SuperkExpectation & Expectation(void) { return *fExp;       }
      PhysParamList &     FitParams  (void) { return *fFitParams; }

      // judge the agreement between the input data and the current prediction
      double Chisquare(Option_t * opt = "poisson");
      
    private:
      SuperkObservation * fObs;
      SuperkExpectation * fExp;
      PhysParamList *     fFitParams;

    ClassDef(FitInputs,1)
   };

private:

   void  Init                    (void);
   void  SetInputs               (OscillationFitter::FitInputs * inputs);
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
   int                  fNonZeroBins;      ///<
   TH2D *               fChisqSurface;
   TString              fOption;
};

//
// MINUIT fit functions
//
void FitFunc1 (
       Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);

} // yr1numu
} // t2kana

#endif

