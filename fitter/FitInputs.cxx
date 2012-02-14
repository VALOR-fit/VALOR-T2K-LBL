//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis) 
 see code in the corresponding SVN branch.
 
 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

//#define __show_bin_contributions_to_chisq__
//#define __show_chisq__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TMinuit.h>
#include <TString.h>

#include "Conventions.h"
#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "FitInputs.h"
#include "Opt.h" 

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
double FitInputs::Chisquare(Option_t * option, conv::SuperkSample_t s0)
{
  TString opt = option;

  //set the oscillation analysis type for chisquare computation
  OscAnaType_t oscana = Opt::Instance()->gOscAnaTypeChoice;
  
  double chisq = 0;

  // loop over samples included in the analysis
  for(int is = 0; is < kns; is++) 
    {
      SuperkSample_t s = kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;

      // If input is s0 = kSKsUndefined then go ahead and calculate a chisq which
      // is the sum of chisq for all samples included in the analysis.
      // If a s0 is specified, then calculate the chisq only for this sample
      if(s0 != conv::kSKsUndefined && s != s0) continue;
      

      //
      // 1Rmu
      //
      if ((s == conv::kSKsFCFV1RingMuLike) && (oscana != conv::kAnaNuEAppOnly))
      //check if we only want the nue appearance analysis, then continue 
	{
	  this->Expectation().Generate(this->FitParams());
	  TH1D * obs_ereco_spectrum = this->Observation().ErecoSpectrum(s);
	  TH1D * exp_ereco_spectrum = this->Expectation().CurrPredictedErecoSpectrum(s);
	  int nbins = exp_ereco_spectrum->GetNbinsX();
	  
#ifdef __show_bin_contributions_to_chisq__
	  cout << "FitInputs::Chisquare()> nbins = " << nbins << endl;
#endif
	  //loop over Ereco bins
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
	}
      
      //
      // 1Re
      //
      else 
	if ((s == conv::kSKsFCFV1RingELike) && (oscana != conv::kAnaNuMuDisappOnly) )
	  {
	    this->Expectation().Generate(this->FitParams());
	    // Note: Counting experiment for the time-being
	    double nobs = this->Observation().ErecoSpectrum(s)->Integral();
	    double nexp = this->Expectation().CurrPredictedErecoSpectrum(s)->Integral();
	    
	    double delta_chisq = 0.;
	    if      (opt.Contains("poisson")) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	    else if (opt.Contains("gauss")  ) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
	    
	    chisq += delta_chisq;

#ifdef __show_chisq__
	    cout << "FitInputs::Chisquare()> chisq = " << chisq 
		 << endl;
#endif
	  }
  
    }//loop over samples
  
  //
  // Add penalty term for parameters declared as nuisance 
  //
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
double FitInputs::BestFitChisquare(Option_t * opt)
{
  fFitParams->CopyBestFitToCurrValues();
  return this->Chisquare(opt);
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
      
      TString smpl = fpar.GetTheSample();

      // judge the agreement between the data and the current prediction
      if(smpl.Contains("numu")) {
            f = inputs->Chisquare("poisson",conv::kSKsFCFV1RingMuLike); 
      }else if(smpl.Contains("nue")){
            f = inputs->Chisquare("poisson",conv::kSKsFCFV1RingELike); 
      }else{
            f = inputs->Chisquare("poisson"); 
      }

}
//____________________________________________________________________________

      

