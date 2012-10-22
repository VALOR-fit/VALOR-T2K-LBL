//____________________________________________________________________________
/*!

\class    t2kana::Superk1RmuErecoFitterMinuit2

\brief    

\author   Thomas Dealtry <thomas.dealtry \at physics.ox.ac.uk>
          STFC, Rutherford Appleton Laboratory
	  Subdepartment of Particle Physics, University of Oxford

\created  October 21, 2012

*/
//____________________________________________________________________________

#ifndef _t2kana_superk_1ring_mulike_ereco_fitter_minuit2_h_
#define _t2kana_superk_1ring_mulike_ereco_fitter_minuit2_h_

#include "FitterKernelMinuit2.h"

class TH2D;

namespace t2kana {

class Superk1RmuErecoFitterMinuit2: public FitterKernelMinuit2
{
public:
   Superk1RmuErecoFitterMinuit2();
   Superk1RmuErecoFitterMinuit2(FitInputs * inputs, Option_t * option = "minuit2");
  ~Superk1RmuErecoFitterMinuit2();

   void Analyze (Option_t * option);
   void Save    (const char * filename);
   void Save    (const char * filename, const char * option_file);

private:

  TH2D * BuildChisqSurfaceNuMuDisap2ParamFit(void);
  TH2D * BuildChisqSurfaceNuMuDisapWithSystParamFit(void);
  TH2D * BuildChisqSurfaceNuMuDecay2ParamFit(void);
  TH2D * BuildChisqSurfaceNuMuDecoh2ParamFit(void);

  TH2D * fDChisq_surface;

};

} // t2kana

#endif

