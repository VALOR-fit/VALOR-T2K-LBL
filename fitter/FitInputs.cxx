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
//#define __show_contributions_to_penalty_term__
//#define __show_chisq__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TMinuit.h>
#include <TString.h>

#include "Conventions.h"
#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "SuperkDetector.h"
#include "BANFF.h"
#include "NuXSec.h"
#include "FSI.h"
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
// If input is s0 = kSKsUndefined then go ahead and calculate a chisq which is the sum of 
// chisq for all SK samples included in the analysis.
// If a s0 is specified, then calculate the chisq only for this sample.
// Add penalty terms as appropriate given which nuisance parameters were included in the
// fit and which oscillation parameters are marginalized.
// The input option can be "poisson" or "gauss".
//

  double chisq = 0; // total chisq to be returned

  TString opt = option;

  // Get the oscillation analysis type for chisquare computation
  OscAnaType_t oscana = Opt::Instance()->gOscAnaTypeChoice;
  
  // loop over samples included in the analysis
  for(int is = 0; is < kns; is++) 
  {
      SuperkSample_t s = kSKSamples[is];
      // current sample is undefined: skip
      if (s == conv::kSKsUndefined) continue;
      // an input sample was specified and the current one is not it: skip
      if(s0 != conv::kSKsUndefined && s != s0) continue;      

      // Add chisq for 1Rmu sample
      // If we only want the nue appearance analysis only, then skip
      if ((s == conv::kSKsFCFV1RingMuLike) && (oscana != conv::kAnaNuEAppOnly))
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
      
      // Add chisq for 1Re sample
      // If we only want the numu disappearance analysis only, then skip
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
  // Add penalty term for nuisance parameters
  //
  double chisq_penalty = this->ChisquarePenaltyTerm();
  chisq += chisq_penalty;      

#ifdef __show_chisq__
  cout << "FitInputs::Chisquare()> chisq = " << chisq 
       << " (contribution from penalty term = " << chisq_penalty << ")" 
       << endl;
#endif

  return chisq; 
}
//____________________________________________________________________________ 
double FitInputs::ChisquarePenaltyTerm(void)
{
// Add penalty term for nuisance parameters

  double chisq_penalty = 0;

  //
  // Add penalty term for BANFF params 
  // First, make sure that all BANFF params were included in the fit param
  // list and were marked as nuisance parameters.
  //
  bool incl_banff_penalty = true;
  const int np = conv::kNumOfParamsInBANFFfit;
  for(int ip = 0; ip < np; ip++) {
     PhysParamId_t p = PhysParamId_t ( (int) kNuisBANFF00 + ip );
     incl_banff_penalty = incl_banff_penalty && 
       (this->FitParams().Exists(p) && this->FitParams().IsNuisance(p));
     if(!incl_banff_penalty) break;
  }
  if(incl_banff_penalty) {
    double dchisq_penalty = this->Expectation().GetBANFF()->ChisqPenalty();
    chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
      cout << "FitInputs::ChisquarePenaltyTerm >> Current BANFF param variation"
           << " => \delta\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }

  //
  // Add penalty term for SuperkDetector params
  // First check all the parms considered by SuperkDetector are included in the fit
  // param list and set as nuisance
  // ...

  bool incl_skdet_penalty = true;
  int np_skdet = 8; // 7 eff + energy scale
  for(int ip = 0; ip < np_skdet; ip++) {
    PhysParamId_t p = PhysParamId_t ( (int) kNuisSuperkEnergyScale  + ip );
    incl_skdet_penalty = incl_skdet_penalty &&
      (this->FitParams().Exists(p) && this->FitParams().IsNuisance(p));
    if(!incl_skdet_penalty) break;
  }
  if(incl_skdet_penalty) {
    double dchisq_penalty = this->Expectation().GetSuperkDetector()->ChisqPenalty();
    chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
    cout << "FitInputs::ChisquarePenaltyTerm >> Current variation of params considered by SuperkDetector"
         << " => \delta\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }


  //
  // Add penalty term for NuXSec params (other than the ones incl. in BANFF)
  //
  // ...

  //
  // Add penalty term for FSI params 
  //
  // ...

  //
  // Add penalty term for marginalised oscillation parameters.
  // Penalize deviations from global analysis best-fit value.
  // ...

  if( this->FitParams().Exists    (kPhypMixing23) && 
      this->FitParams().IsNuisance(kPhypMixing23) )
  {
    double mix = this->FitParams().CurrValue(kPhypMixing23);
    double nom = 0;
    double err = 0;
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    switch(omp) {
     case(kNuMixSinsq2Theta):
          nom = kGlobFitSinsq2Theta23;
          err = kGlobFitSinsq2Theta23 * kGlobFitErrSinsq2Theta23;
          break;
     case(kNuMixSinsqTheta):
          nom = kGlobFitSinsqTheta23;
          err = kGlobFitSinsqTheta23 * kGlobFitErrSinsqTheta23;
          break;
     default:
          exit(1);
    }
    double twk = (mix - nom) / err;
    double dchisq_penalty = twk*twk;
    chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
    cout << "FitInputs::ChisquarePenaltyTerm >> mixing23 = " << mix
         << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
         << " => \delta\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }//mix23

  if( this->FitParams().Exists    (kPhypMixing13) && 
      this->FitParams().IsNuisance(kPhypMixing13) )
  {
     double mix = this->FitParams().CurrValue(kPhypMixing13);
     double nom = 0;
     double err = 0;
     NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
     switch(omp) {
      case(kNuMixSinsq2Theta):
            nom = kGlobFitSinsq2Theta13;
            err = kGlobFitSinsq2Theta13 * kGlobFitErrSinsq2Theta13;
            break;
      case(kNuMixSinsqTheta):
            nom = kGlobFitSinsqTheta13;
            err = kGlobFitSinsqTheta13 * kGlobFitErrSinsqTheta13;
            break;
      default:
            exit(1);
     }
     double twk = (mix - nom) / err;
     double dchisq_penalty = twk*twk;
     chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
     cout << "FitInputs::ChisquarePenaltyTerm >> mixing13 = " << mix
          << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
          << " => \delta\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }//mix13


  if( this->FitParams().Exists    (kPhypMixing12) && 
      this->FitParams().IsNuisance(kPhypMixing12) )
  {
     double mix = this->FitParams().CurrValue(kPhypMixing12);
     double nom = 0;
     double err = 0;
     NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
     switch(omp) {
      case(kNuMixSinsq2Theta):
            nom = kGlobFitSinsq2Theta12;
            err = kGlobFitSinsq2Theta12 * kGlobFitErrSinsq2Theta12;
            break;
      case(kNuMixSinsqTheta):
            nom = kGlobFitSinsqTheta12;
            err = kGlobFitSinsqTheta12 * kGlobFitErrSinsqTheta12;
            break;
      default:
            exit(1);
     }
     double twk = (mix - nom) / err;
     double dchisq_penalty = twk*twk;
     chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
     cout << "FitInputs::ChisquarePenaltyTerm >> mixing12 = " << mix
          << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
          << " => \delta\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }//mix12

  if( this->FitParams().Exists    (kPhypDm12sq) && 
      this->FitParams().IsNuisance(kPhypDm12sq) )
  {
     double dmsq = this->FitParams().CurrValue(kPhypDm12sq);
     double nom  = kGlobFitDm12sq;
     double err  = kGlobFitDm12sq * kGlobFitErrDm12sq;
     double twk  = (dmsq - nom) / err;
     double dchisq_penalty = twk*twk;
     chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
     cout << "FitInputs::ChisquarePenaltyTerm >> dm12sq = " << dmsq
          << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
          << " => \delta\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }//dm12sq

  if( this->FitParams().Exists    (kPhypDm23sq) && 
      this->FitParams().IsNuisance(kPhypDm23sq) )
  {
      double dmsq = this->FitParams().CurrValue(kPhypDm23sq);
      double nom  = kGlobFitDm23sq;
      double err  = kGlobFitDm23sq * kGlobFitErrDm23sq;
      double twk  = (dmsq - nom) / err;
      double dchisq_penalty = twk*twk;
      chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
      cout << "FitInputs::ChisquarePenaltyTerm >> |dm23sq| = " << dmsq
           << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
           << " => \delta\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }//dm23sq


  return chisq_penalty;
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

     // judge the agreement between the data and the current prediction
     f = inputs->Chisquare("poisson"); 
      
//      TString smpl = fpar.GetTheSample();
//      // judge the agreement between the data and the current prediction
//      if(smpl.Contains("numu")) {
//            f = inputs->Chisquare("poisson",conv::kSKsFCFV1RingMuLike); 
//      }else if(smpl.Contains("nue")){
//            f = inputs->Chisquare("poisson",conv::kSKsFCFV1RingELike); 
//      }else{
//            f = inputs->Chisquare("poisson"); 
//      }

}
//____________________________________________________________________________

      

