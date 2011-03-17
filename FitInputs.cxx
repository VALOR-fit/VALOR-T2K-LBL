//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Created from the old NuMuDisappFitter::FitInputs nested class. Also replaces similar
   nested classes which were located in the NuMuDecayFitter and NuMuDecohFitter classes.
*/
//________________________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TMinuit.h>
#include <TString.h>

#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "FitInputs.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;

ClassImp(FitInputs)

//________________________________________________________________________________________
FitInputs::FitInputs(
  SuperkObservation * skobs, SuperkExpectation * skexp, PhysParamList * fpar) :
fObs(skobs), 
fExp(skexp), 
fFitParams(fpar) 
{ 

}
//________________________________________________________________________________________
FitInputs::~FitInputs()  
{ 

}
//________________________________________________________________________________________
SuperkObservation & FitInputs::Observation(void) 
{ 
  return *fObs;       
}
//________________________________________________________________________________________
SuperkExpectation & FitInputs::Expectation(void) 
{ 
  return *fExp;       
}
//________________________________________________________________________________________
PhysParamList & FitInputs::FitParams(void) 
{ 
  return *fFitParams; 
}
//________________________________________________________________________________________
double FitInputs::Chisquare(Option_t * option)
{
  double chisq = 0;
  TString opt = option;

  this->Expectation().Generate(this->FitParams());
  TH1D * obs_ereco_spectrum = this->Observation().ErecoSpectrum();
  TH1D * exp_ereco_spectrum = this->Expectation().CurrPredictedErecoSpectrum();

  int nbins = exp_ereco_spectrum->GetNbinsX();

#ifdef __debug_level_3__
   cout << "FitInputs::Chisquare()> nbins = " << nbins << endl;
#endif

  for(int ibin = 1; ibin <= nbins; ibin++) {

     double nobs = obs_ereco_spectrum->GetBinContent(ibin);
     double nexp = exp_ereco_spectrum->GetBinContent(ibin);

     double delta = 0.;
     if      (opt.Contains("poisson")) { delta = util::ChisqPoisson(nobs,nexp); }
     else if (opt.Contains("gauss")  ) { delta = util::ChisqGauss  (nobs,nexp); }

     chisq += delta;

#ifdef __debug_level_3__
     cout << "FitInputs::Chisquare()> chisq (bin = "
          << ibin << " / Ereco = " << Ereco << " GeV) = " << delta << endl;
#endif

  }//ibin

#ifdef __debug_level_2__
  cout << "FitInputs::Chisquare()> chisq = " << chisq << endl;
#endif

  return chisq;
}
//____________________________________________________________________________ 
// MINUIT fit functions 
// 
void t2kana::FitFunc (
  Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/) 
{
  FitInputs * inputs = (FitInputs *) gMinuit->GetObjectFit();

  // get fit params
  PhysParamList & fpar = inputs->FitParams(); 

  // synchronize PhysParamList with the paramerer array used by MINUIT
  fpar.Sync(par);

  // judge the agreement between the data and the current prediction
  f = inputs->Chisquare("poisson");
}
//____________________________________________________________________________

      

