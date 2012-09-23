//____________________________________________________________________________
/*!

\class    t2kana::Superk1RmuErecoFitter

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_superk_1ring_mulike_ereco_fitter_h_
#define _t2kana_superk_1ring_mulike_ereco_fitter_h_

#include "FitterKernel.h"

class TH2D;

namespace t2kana {

class Superk1RmuErecoFitter: public FitterKernel
{
public:
   Superk1RmuErecoFitter();
   Superk1RmuErecoFitter(FitInputs * inputs, Option_t * option = "minuit");
  ~Superk1RmuErecoFitter();

   void Analyze (Option_t * option);
   void Save    (const char * filename);
   void Save    (const char * filename, const char * option_file);

private:

  TH2D * BuildChisqSurfaceNuMuDisap2ParamFit(void);
  TH2D * BuildChisqSurfaceNuMuDisapWithSystParamFit(void);
  TH2D * BuildChisqSurfaceNuMuDecay2ParamFit(void);
  TH2D * BuildChisqSurfaceNuMuDecoh2ParamFit(void);

};

} // t2kana

#endif

