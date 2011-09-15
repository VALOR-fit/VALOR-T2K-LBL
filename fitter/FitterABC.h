//____________________________________________________________________________
/*!

\class    t2kana::FitterABC

\brief    Fitter ABC

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 05, 2011

*/
//____________________________________________________________________________

#ifndef _t2kana_fitter_abc_h_
#define _t2kana_fitter_abc_h_

#include <TString.h>

class TVirtualPad;

namespace t2kana {

class FitterABC
{
public:
   virtual ~FitterABC();

   virtual bool DoTheFit      (void) = 0;
   virtual int  FitStatusCode (void) = 0;
   virtual int  NDOF          (void) { return this->NBinsInFit() - this->NFitParams(); }
   virtual int  NBinsInFit    (void) = 0;
   virtual int  NFitParams    (void) = 0;
   virtual void Draw          (TVirtualPad * c, Option_t * option, const char * outtag) = 0;

protected:
   FitterABC();

};

} // t2kana

#endif

