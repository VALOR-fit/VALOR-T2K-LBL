//____________________________________________________________________________
/*!

\class    t2kana::Superk1ReNFitter

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_superk_1ring_elike_n_fitter_h_
#define _t2kana_superk_1ring_elike_n_fitter_h_

#include "FitterKernel.h"

class TH2D;

namespace t2kana {

class Superk1ReNFitter: public FitterKernel
{
public:
   Superk1ReNFitter();
   Superk1ReNFitter(FitInputs * inputs, Option_t * option = "minuit");
  ~Superk1ReNFitter();

   void Analyze (Option_t * option);
   void Save    (const char * filename);

private:

   TH2D * BuildChisqSurfaceNuEApp2ParamFit(void);

};

} // t2kana

#endif

