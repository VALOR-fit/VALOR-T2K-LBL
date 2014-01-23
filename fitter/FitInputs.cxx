//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis) 
 see code in the corresponding SVN branch.
 
 Changes since 2011_v0.9.0:

 @ 18/2/13 LE
 Added back Fogli and Lisi schema in the penalty term.

 30-12-2013 Nick
 Added samples for reverse horn current running

*/
//________________________________________________________________________________________

//#define __show_bin_contributions_to_chisq__
//#define __show_contributions_to_penalty_term__
//#define __show_chisq__
#define __include_all_systematics_in_banff_penalty_term__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TMinuit.h>
#include <TString.h>
#include <TMath.h>
#include "Math/FitMethodFunction.h"
#include "HFitInterface.h"

#include "Conventions.h"
#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "SuperkDetector.h"
#include "SuperkDetectorFSI.h"
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
	  SuperkObservation * skobs, SuperkExpectation * skexp, FitParamList * fpar, string opt, double errdef, bool goodnessoffit) :
fObs(skobs), 
fExp(skexp), 
fFitParams(fpar),
fOpt(opt),
fErrDef(errdef),
fGoodnessOfFit(goodnessoffit)
{ 
  fNPoints = kNErecoBins[0]+kNErecoBins[1];
  fNDim = (*fFitParams).GetNumberParameters();
  if(fFitParams->GetSystMode() == 2)
    {
      fAllBANFFFit = false;
      fAll1ReEff = false;
    }
  else{
    fAllBANFFFit = true;
    fAll1ReEff = true;
  }
  TString option = opt.c_str();
  if(option.Contains("poisson"))
    fMode = 1;
  else if(option.Contains("gauss"))
    fMode = 2;
  else
    fMode = -1;
}
//________________________________________________________________________________________
FitInputs::~FitInputs()  
{ 

}
//________________________________________________________________________________________
void FitInputs::UpdateParamList(FitParamList * fpar)
{
  fFitParams = fpar;
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
double FitInputs::Chisquare(Option_t * option, conv::SuperkSample_t s0, int point) const
{
  TString opt = option;
  if(opt.Contains("poisson"))
    return Chisquare(1,s0,point);
  else if(opt.Contains("gauss"))
    return Chisquare(2,s0,point);
  return Chisquare(0,s0,point);
}
//________________________________________________________________________________________
double FitInputs::Chisquare(int mode, conv::SuperkSample_t s0, int point) const
// If input is s0 = kSKsUndefined then go ahead and calculate a chisq which is the sum of 
// chisq for all SK samples included in the analysis.
// If a s0 is specified, then calculate the chisq only for this sample.
// Add penalty terms as appropriate given which nuisance parameters were included in the
// fit and which oscillation parameters are marginalized.
// The input option can be 1 ("poisson") or 2 ("gauss").
//
{
  double chisq = 0; // total chisq to be returned

  //TString opt = option;

  // Get the oscillation analysis type for chisquare computation
  OscAnaType_t oscana = Opt::Instance()->gOscAnaTypeChoice;
  bool applied_negevents_condition = false;

  // loop over samples included in the analysis
  for(int is = 0; is < kns; is++) 
  {
      SuperkSample_t s = kSKSamples[is];
      // current sample is undefined: skip
      if (s == conv::kSKsUndefined) continue;
      // an input sample was specified and the current one is not it: skip
      if(s0 != conv::kSKsUndefined && s != s0) continue;      

      // Add chisq for each sample
      if ((s == conv::kSKsFCFV1RingMuLike) && (oscana == conv::kAnaNuEAppOnly || oscana == conv:: kAnaRHCCombined
	   || oscana == conv::kAnaRHCNuMuDisappOnly || oscana == conv::kAnaRHCNuEAppOnly || oscana == conv::kAnaFHCAndRHCNuEAppOnly))
	continue;
      if ((s == conv::kSKsRHCFCFV1RingMuLike) && (oscana == conv:: kAnaCombined || oscana == conv::kAnaNuMuDisappOnly 
           || oscana == conv::kAnaNuEAppOnly || oscana == conv::kAnaRHCNuEAppOnly || oscana == conv::kAnaFHCAndRHCNuEAppOnly))
	continue;
      if ((s == conv::kSKsFCFV1RingELike) && (oscana == conv::kAnaNuMuDisappOnly || oscana == conv:: kAnaRHCCombined 
           || oscana == conv::kAnaRHCNuMuDisappOnly || oscana == conv::kAnaRHCNuEAppOnly || oscana == conv::kAnaFHCAndRHCNuMuDisappOnly))
	continue;
      if ((s == conv::kSKsRHCFCFV1RingELike) && (oscana == conv:: kAnaCombined || oscana == conv::kAnaNuMuDisappOnly 
           || oscana == conv::kAnaNuEAppOnly || oscana == conv::kAnaRHCNuMuDisappOnly || oscana == conv::kAnaFHCAndRHCNuMuDisappOnly))
	continue;

      int n_ereco_bins = kNErecoBins[is];

      applied_negevents_condition = fExp->Generate(*fFitParams, fAllBANFFFit, fAll1ReEff) || applied_negevents_condition;
      if(!fGoodnessOfFit) {
#ifdef __show_bin_contributions_to_chisq__
	cout << "FitInputs::Chisquare()> nbins = " << nbins << endl;
#endif
	//loop over Ereco bins
	for(int ibin = 1; ibin <= n_ereco_bins; ibin++) {
	  if((point!=-1) && (point != ibin)) continue;
	  double nobs = fObs->ErecoSpectrumArray             (s)[ibin-1];
	  double nexp = fExp->CurrPredictedErecoSpectrumArray(s)[ibin-1];
	  
	  double delta_chisq = 0.;
	  //if      (opt.Contains("poisson")) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	  //else if (opt.Contains("gauss")  ) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
	  if      (mode == 1) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	  else if (mode == 2) { delta_chisq = util::ChisqGauss  (nobs,nexp); }

	  if((point!=-1) && (point == ibin)) return delta_chisq;
	  chisq += delta_chisq;
	  
#ifdef __show_bin_contributions_to_chisq__
	  double xlo = util::ErecoBinEdges(s,ibin-1);
	  double xup = util::ErecoBinEdges(s,ibin);
	    cout << "FitInputs::Chisquare()> delta_chisq (bin = " << ibin 
		 << "; [" << xlo << ", " << xup << "] GeV) = " << delta_chisq << endl;
	    cout << "\tNobs: " << nobs << "\tNexp:" << nexp << endl;
#endif
	}//ibin
	
      }//!goodnessoffit
      else {
	const int nbins = 5;
	double binedges[nbins+1] = {0.0,0.4,0.7,1.0,2.0,30.0}; 
	
#ifdef __show_bin_contributions_to_chisq__
	cout << "FitInputs::Chisquare()> nbins = " << nbins << endl;
#endif
	
	  //loop over Ereco bins
	  for(int ibin = 1; ibin <= nbins; ibin++) {

	    int ErecoMinBin = 1 + TMath::BinarySearch(n_ereco_bins,util::ErecoBinEdges(s),binedges[ibin-1]);
	    int ErecoMaxBin = 1 + TMath::BinarySearch(n_ereco_bins,util::ErecoBinEdges(s),binedges[ibin]);
	    double nobs = 0, nexp = 0;
	    for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	      nobs += fObs->ErecoSpectrumArray             (s)[ir-1];
	      nexp += fExp->CurrPredictedErecoSpectrumArray(s)[ir-1];
	    }
	    double delta_chisq = 0.;
	    //if      (opt.Contains("poisson")) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	    //else if (opt.Contains("gauss")  ) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
	    if      (mode == 1) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	    else if (mode == 2) { delta_chisq = util::ChisqGauss  (nobs,nexp); }

	    chisq += delta_chisq;
	      
#ifdef __show_bin_contributions_to_chisq__
	    double xlo = util::ErecoBinEdges(s,ibin-1);
	    double xup = util::ErecoBinEdges(s,ibin);
	    cout << "FitInputs::Chisquare()> delta_chisq (bin = " << ibin 
		 << "; [" << xlo << ", " << xup << "] GeV) = " << delta_chisq << endl;
	    cout << "\tNobs: " << nobs << "\tNexp:" << nexp << endl;
#endif
	  }//ibin
      }//goodnessoffit
      //      }//1Rmu
      /*
      // Add chisq for 1Re sample
      // If we only want the numu disappearance analysis, then skip
      else 
      if ((s == conv::kSKsFCFV1RingELike) && (oscana != conv::kAnaNuMuDisappOnly) )
      {
	    applied_negevents_condition = fExp->Generate(*fFitParams, fAllBANFFFit) || applied_negevents_condition;
	    // Note: Counting experiment for the time-being
	    double nobs = fObs->ErecoSpectrum(s)->Integral();
	    double nexp = fExp->CurrPredictedErecoSpectrum(s)->Integral();
	    
	    double delta_chisq = 0.;
	    //if      (opt.Contains("poisson")) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	    //else if (opt.Contains("gauss")  ) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
	    if      (mode == 1) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	    else if (mode == 2) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
	    
	    chisq += delta_chisq;

#ifdef __show_chisq__
	    cout << "FitInputs::Chisquare()> chisq = " << chisq 
		 << endl;
#endif
      }
      */  

  }//loop over samples
  
  //
  // Add penalty term for nuisance parameters
  //
  double chisq_penalty = ChisquarePenaltyTerm();
  chisq += chisq_penalty;      
  
#ifdef __show_chisq__
  cout << "FitInputs::Chisquare()> chisq = " << chisq 
       << " (contribution from penalty term = " << chisq_penalty << ")" 
       << endl;
#endif
  
  //if(applied_negevents_condition) chisq += 1E8;
  
  return chisq; 
}
//____________________________________________________________________________ 
double FitInputs::ChisquarePenaltyTerm(void) const
{
// Add penalty term for nuisance parameters

  double chisq_penalty = 0;

  // FitParamList's without nuisance parameters don't have penalty terms
  int systmode = fFitParams->GetSystMode();
  if(systmode == 0)
    return chisq_penalty;

  // Check whether its a full or partial set of systematics
  bool all_syst = true;
  if(systmode == 2)
    all_syst = false;
  //check whether we only want the oscillation parameter tweaks
  bool osc_only = false;
  if(systmode == 3 || systmode == 4)
    osc_only = true;
  bool hk_hack = false;
  if(systmode == 5 || systmode == 6)
    hk_hack = true;
  
  //
  // Add penalty term for marginalised oscillation parameters.
  // Penalize deviations from global analysis best-fit value.
  // ...
  if(systmode <= 3 || systmode == 6) {
#ifdef HK
    if( fFitParams->Exists    (kPhypMixing23) && 
	fFitParams->IsNuisance(kPhypMixing23) )
      {
	double mix = fFitParams->CurrValue(kPhypMixing23);
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
	     << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      }//mix23
#endif//HK

    if( fFitParams->Exists    (kPhypMixing13) && 
	fFitParams->IsNuisance(kPhypMixing13) )
      {
	double mix = fFitParams->CurrValue(kPhypMixing13);
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
	     << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      }//mix13


    if( fFitParams->Exists    (kPhypMixing12) && 
	fFitParams->IsNuisance(kPhypMixing12) )
      {
	double mix = fFitParams->CurrValue(kPhypMixing12);
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
	     << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      }//mix12

    if( fFitParams->Exists    (kPhypDm12sq) && 
	fFitParams->IsNuisance(kPhypDm12sq) )
      {
	double dmsq = fFitParams->CurrValue(kPhypDm12sq);
	double nom  = kGlobFitDm12sq;
	double err  = kGlobFitDm12sq * kGlobFitErrDm12sq;
	double twk  = (dmsq - nom) / err;
	double dchisq_penalty = twk*twk;
	chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
	cout << "FitInputs::ChisquarePenaltyTerm >> dm12sq = " << dmsq
	     << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
	     << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      }//dm12sq
  

#ifdef HK
    if( fFitParams->Exists    (kPhypDm23sq) && 
	fFitParams->IsNuisance(kPhypDm23sq) )
      {
	double dmsq = fFitParams->CurrValue(kPhypDm23sq);
	double nom = 0.;
	NuOscDmsqScheme_t odmsq = Opt::Instance()->gOscDmsqScheme;
	switch (odmsq) {
	case(kNuDmsqFogliLisi):
	  nom  =  kGlobFitDmsqFL;
	  break;
	case(kNuDmsq12andAbs23):
	  nom  = kGlobFitDm23sq;
	  break;
	case(kNuDmsq12andAbs31):
	  if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal)
	    nom  = kGlobFitDm23sq + kGlobFitDm12sq; 
	  else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted)
	    nom  = kGlobFitDm23sq - kGlobFitDm12sq; 
	  break;
	default:
	  nom  = kGlobFitDm23sq;
	}
	double err  = nom * kGlobFitErrDm23sq;
	double twk  = (dmsq - nom) / err;
	double dchisq_penalty = twk*twk;
	chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
	cout << "FitInputs::ChisquarePenaltyTerm >> |dm23sq| = " << dmsq
	     << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
	     << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      }//dm23sq
#endif//HK
  }//CURRENT osc params
  else {//FUTURE osc params
    if( fFitParams->Exists    (kPhypMixing13) && 
	fFitParams->IsNuisance(kPhypMixing13) )
      {
	double mix = fFitParams->CurrValue(kPhypMixing13);
	double nom = 0;
	double err = 0;
	NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
	switch(omp) {
	case(kNuMixSinsq2Theta):
	  nom = kGlobFitSinsq2Theta13;
	  err = kFutureAbsErrSinsq2Theta13;
	  break;
	case(kNuMixSinsqTheta):
	  nom = kGlobFitSinsqTheta13;
	  err = kFutureAbsErrSinsqTheta13;
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
	     << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      }//mix13

    if( fFitParams->Exists    (kPhypMixing23) && 
	fFitParams->IsNuisance(kPhypMixing23) )
      {
	double mix = fFitParams->CurrValue(kPhypMixing23);
	double nom = 0;
	double err = 0;
	NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
	switch(omp) {
	case(kNuMixSinsq2Theta):
	  nom = kGlobFitSinsq2Theta23;
	  err = kFutureAbsErrSinsq2Theta23;
	  break;
	case(kNuMixSinsqTheta):
	  nom = kGlobFitSinsqTheta23;
	  err = kFutureAbsErrSinsqTheta23;
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
	     << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      }//mix23

    if( fFitParams->Exists    (kPhypDm23sq) && 
	fFitParams->IsNuisance(kPhypDm23sq) )
      {
	double dmsq = fFitParams->CurrValue(kPhypDm23sq);
	double nom = 0.;
	NuOscDmsqScheme_t odmsq = Opt::Instance()->gOscDmsqScheme;
	switch (odmsq) {
	case(kNuDmsqFogliLisi):
	  nom  =  kGlobFitDmsqFL;
	  break;
	case(kNuDmsq12andAbs23):
	  nom  = kGlobFitDm23sq;
	  break;
	case(kNuDmsq12andAbs31):
	  if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal)
	    nom  = kGlobFitDm23sq + kGlobFitDm12sq; 
	  else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted)
	    nom  = kGlobFitDm23sq - kGlobFitDm12sq; 
	  break;
	default:
	  nom  = kGlobFitDm23sq;
	}
	double err  = kFutureAbsErrDm23sq;
	double twk  = (dmsq - nom) / err;
	double dchisq_penalty = twk*twk;
	chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
	cout << "FitInputs::ChisquarePenaltyTerm >> |dm23sq| = " << dmsq
	     << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
	     << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      }//dm23sq
  }

  //got all the osc param penalty. check if we need to exit
  if(osc_only)
    return chisq_penalty;

  if(hk_hack) {
    double dchisq_penalty_xsec = fExp->GetNuXSec()->ChisqPenalty(all_syst,true);
    chisq_penalty += dchisq_penalty_xsec;
#ifdef __show_contributions_to_penalty_term__
    cout << "FitInputs::ChisquarePenaltyTerm >> Current variation of params considered by NuXSec"
	 << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_xsec << endl;
#endif
    return chisq_penalty;
  }

  //
  // Add penalty term for BANFF params 
  //
#ifdef __include_all_systematics_in_banff_penalty_term__
  double dchisq_penalty_banff = fExp->GetBANFF()->ChisqPenalty();
#else
  double dchisq_penalty_banff = fExp->GetBANFF()->ChisqPenalty(all_syst);
#endif
  chisq_penalty += dchisq_penalty_banff;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::ChisquarePenaltyTerm >> Current BANFF param variation"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_banff << endl;
#endif

  //
  // Add penalty term for NuXSec params (other than the ones incl. in BANFF)
  //
  double dchisq_penalty_xsec = fExp->GetNuXSec()->ChisqPenalty(all_syst);
  chisq_penalty += dchisq_penalty_xsec;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::ChisquarePenaltyTerm >> Current variation of params considered by NuXSec"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_xsec << endl;
#endif

  //
  // Add penalty term for SuperkDetector params
  //
#ifndef SKFSITOGETHER
  double dchisq_penalty_skdet = fExp->GetSuperkDetector()->ChisqPenalty(all_syst);
  chisq_penalty += dchisq_penalty_skdet;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::ChisquarePenaltyTerm >> Current variation of params considered by SuperkDetector"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_skdet << endl;
#endif
#endif

  //
  // Add penalty term for FSI params 
  //
#ifndef SKFSITOGETHER
  double dchisq_penalty_fsi = fExp->GetFSI()->ChisqPenalty();
  chisq_penalty += dchisq_penalty_fsi;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::ChisquarePenaltyTerm >> Current FSI param variation"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_fsi << endl;
#endif
#endif
    
  //
  // Add penalty term for SuperkDetector + FSI params
  //
#ifdef SKFSITOGETHER
  double dchisq_penalty_skdetfsi = fExp->GetSuperkDetectorFSI()->ChisqPenalty();
  chisq_penalty += dchisq_penalty_skdetfsi;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::ChisquarePenaltyTerm >> Current SK+FSI param variation"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_skdetfsi << endl;
#endif
#endif

  return chisq_penalty;
}
//____________________________________________________________________________ 
double FitInputs::BestFitChisquare(Option_t * opt)
{
  fFitParams->CopyBestFitToCurrValues();
  return Chisquare(opt);
}
//____________________________________________________________________________ 
// MINUIT fit functions 
// 
//____________________________________________________________________________ 
double FitInputs::operator() (const std::vector<double>& x) const
{
// To be used with Minuit2

 fFitParams->Sync(x);

 //return  this->Chisquare(fOpt.c_str());
 return  this->Chisquare(fMode);
}
//____________________________________________________________________________ 
void t2kana::FitFunc (
  Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/) 
{
// To be used with a TVirtualFitter

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
FitFunctor::FitFunctor(FitInputs * inputs) 
{ 
  fFitInputs = inputs; 
}
FitFunctor::~FitFunctor() 
{

}
double FitFunctor::operator() (const double * x) const
{
  return this->DoEval(x);
}
double FitFunctor::DoEval(const double * x) const 
{
  fFitInputs->FitParams().Sync(x);
  return fFitInputs->Chisquare("poisson");
}
ROOT::Math::IBaseFunctionMultiDim * FitFunctor::Clone(void) const 
{
  return new FitFunctor(fFitInputs);
}
unsigned int FitFunctor::NDim(void) const 
{ 
  return fFitInputs->FitParams().GetNumberParameters();
}
//____________________________________________________________________________
FitMethodFunc::~FitMethodFunc() 
{

}
double FitMethodFunc::operator() (const double * x) const
{
  return this->DoEval(x);
}
double FitMethodFunc::DoEval(const double * x) const 
{
  fFitInputs->FitParams().Sync(x);
  return fFitInputs->Chisquare("poisson");
}
ROOT::Math::BasicFitMethodFunction<ROOT::Math::IBaseFunctionMultiDim>::Type_t FitMethodFunc::Type(void) const
{
  return ROOT::Math::FitMethodFunction::kLeastSquare;
}


double  FitMethodFunc::DataElement(const double *par, unsigned int ipoint, double * /* g */) const 
{
  fFitInputs->FitParams().Sync(par);
  double chisqbin = 0;
  //Todo: make this work when both __FHC_RUNNING__ and __RHC_RUNNING__ are defined.
#ifdef __FHC_RUNNING__
  if (ipoint < (unsigned int) conv::kNErecoBins1Rmu) chisqbin = fFitInputs->Chisquare("poisson",conv::kSKsFCFV1RingMuLike,ipoint); 
  else chisqbin = fFitInputs->Chisquare("poisson",conv::kSKsFCFV1RingELike,ipoint-conv::kNErecoBins1Rmu); 
#endif
#ifdef __RHC_RUNNING__
  if (ipoint < (unsigned int) conv::kRHCNErecoBins1Rmu) chisqbin = fFitInputs->Chisquare("poisson",conv::kSKsRHCFCFV1RingMuLike,ipoint);
  else chisqbin = fFitInputs->Chisquare("poisson",conv::kSKsRHCFCFV1RingELike,ipoint-conv::kRHCNErecoBins1Rmu);
#endif
  //  if (g!=0) ;
  //  cout << "DataElement: chisqbin = " << chisqbin << " for point = " << ipoint << endl;
  return chisqbin;
}

