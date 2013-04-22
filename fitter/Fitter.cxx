#include <TString.h>

#include "FitterABC.h"
#include "FitterGSL.h"
#include "FitterMinuit.h"
#include "FitterMInuit2.h"

#include "Fitter.h"

using namespace t2kana;

FitterABC * GetFitter (FitInputs * inputs, Option_t * option)
{
  if(!inputs) return 0;
  if(!option) return 0;

  TString opt(option);

  if ( opt.Contains("minuit2") ) 
  {
     return new FitterMinuit2(inputs,option);
  }

  if ( opt.Contains("minuit") ) 
  {
     return new FitterMinuit(inputs,option);
  }

  if ( opt.Contains("SimAnneal")   ||
       opt.Contains("LevMarq")     ||
       opt.Contains("ConjugateFR") || 
       opt.Contains("ConjugatePR") ||
       opt.Contains("VectorBFGS")  ||
       opt.Contains("SteepestDescent")
     )
  {
     return new FitterGSL(inputs,option);     
  }

  return 0;
}


