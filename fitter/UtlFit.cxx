#include <iostream>

#include <TString.h>

#include "FitterABC.h"
#include "FitterGSL.h"
#include "FitterMinuit.h"
#include "FitterMinuit2.h"

#include "UtlFit.h"

using namespace t2kana;

using std::cerr;
using std::cout;
using std::endl;

FitterABC * t2kana::util::GetFitter (FitInputs * inputs, Option_t * option)
{
  if(!inputs) return 0;
  if(!option) return 0;

  TString opt(option);

  if ( opt.Contains("minuit2") ) 
  {
     cout << "Fitter >> Instantianting MINUIT2 fitter with option: " << option << endl; 
     return new FitterMinuit2(inputs,option);
  }

  if ( opt.Contains("minuit") ) 
  {
     cout << "Fitter >> Instantianting MINUIT fitter with option: " << option << endl; 
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
     cout << "Fitter >> Instantianting GSL fitter with option: " << option << endl; 
     return new FitterGSL(inputs,option);     
  }

  cerr << "Fitter >> No fitter was instantiated. Input fitter option was: " << option << endl; 
  return 0;
}


