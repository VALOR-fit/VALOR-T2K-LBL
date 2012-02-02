//____________________________________________________________________________
/*!

\class    t2kana::SuperkCombinedFitter

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_superk_combined_fitter_h_
#define _t2kana_superk_combined_fitter_h_

#include "FitterKernel.h"
#include "Superk1RmuErecoFitter.h"
#include "Superk1ReNFitter.h"

namespace t2kana {

  class SuperkCombinedFitter: public FitterKernel
{
public:
   SuperkCombinedFitter();
   SuperkCombinedFitter(FitInputs * inputs, Option_t * option = "minuit");
  ~SuperkCombinedFitter();

   void Analyze (Option_t * option);
   void Save    (const char * filename);
   void Save    (const char * filename, const char * opt_file);

private:

   FitterABC * fitter1Rmu; 
   FitterABC * fitter1Re; 

};

} // t2kana

#endif

