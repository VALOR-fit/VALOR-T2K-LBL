//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Created from the old NuMuDisappFitter::FitInputs nested class. Also replaces similar
   nested classes which were located in the NuMuDecayFitter and NuMuDecohFitter classes.
 @ 28/03/2011 - CA
   In Chisquare(), make it clear that we are fitting only the 1-R mu-like sample at this
   point.
 @ 05/07/2011 - CA
   Adapt code since PhysParamList was renamed to FitParamList. Added penalty term in
   chi^2 since fitting for systematics in now enabled.
*/
//________________________________________________________________________________________

#define __show_bin_contributions_to_chisq__
#define __show_chisq__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TMinuit.h>
#include <TString.h>

#include "Conventions.h"
#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "FitInputs.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;

ClassImp(FitInputs)

//________________________________________________________________________________________
FitInputs::FitInputs(
  SuperkObservation * skobs, SuperkExpectation * skexp, FitParamList * fpar) :
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
FitParamList & FitInputs::FitParams(void) 
{ 
  return *fFitParams; 
}
//________________________________________________________________________________________
double FitInputs::Chisquare(Option_t * option)
{
  double chisq = 0;
  TString opt = option;

  this->Expectation().Generate(this->FitParams());
  TH1D * obs_ereco_spectrum = this->Observation().ErecoSpectrum(kSKsFCFV1RingMuLike);
  TH1D * exp_ereco_spectrum = this->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike);

  int nbins = exp_ereco_spectrum->GetNbinsX();

#ifdef __show_bin_contributions_to_chisq__
   cout << "FitInputs::Chisquare()> nbins = " << nbins << endl;
#endif

  for(int ibin = 1; ibin <= nbins; ibin++) {

     double nobs = obs_ereco_spectrum->GetBinContent(ibin);
     double nexp = exp_ereco_spectrum->GetBinContent(ibin);

     double delta_chisq = 0.;
     if      (opt.Contains("poisson")) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
     else if (opt.Contains("gauss")  ) { delta_chisq = util::ChisqGauss  (nobs,nexp); }

     chisq += delta_chisq;

#ifdef __show_bin_contributions_to_chisq__
     double xlo = obs_ereco_spectrum->GetXaxis()->GetBinLowEdge(ibin);  
     double xup = obs_ereco_spectrum->GetXaxis()->GetBinUpEdge(ibin);
     cout << "FitInputs::Chisquare()> delta_chisq (bin = " << ibin 
          << "; [" << xlo << ", " << xup << "] GeV) = " << delta_chisq << endl;
#endif

  }//ibin

  // add penalty term for parameters declared as nuisance (and having values
  // in units of sigma)
  double chisq_penalty = this->FitParams().ChisquarePenaltyTerm();
  chisq += chisq_penalty;

#ifdef __show_chisq__
  cout << "FitInputs::Chisquare()> chisq = " << chisq 
       << " (contribution from penalty term = " << chisq_penalty << ")" 
       << endl;
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
  FitParamList & fpar = inputs->FitParams(); 

  // synchronize FitParamList with the paramerer array used by MINUIT
  fpar.Sync(par);

  // judge the agreement between the data and the current prediction
  f = inputs->Chisquare("poisson");
}
//____________________________________________________________________________

      

