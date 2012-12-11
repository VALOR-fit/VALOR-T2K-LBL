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
#define __include_all_systematics_in_banff_penalty_term__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TMinuit.h>
#include <TString.h>
#include <TMath.h>

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
	  SuperkObservation * skobs, SuperkExpectation * skexp, FitParamList * fpar, string opt, double errdef, bool goodnessoffit) :
fObs(skobs), 
fExp(skexp), 
fFitParams(fpar),
fOpt(opt),
fErrDef(errdef),
fGoodnessOfFit(goodnessoffit)
{ 
  if(fFitParams->GetSystMode() == 2)
    fAllBANFFFit = false;
  else
    fAllBANFFFit = true;
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
double FitInputs::Chisquare(Option_t * option, conv::SuperkSample_t s0) const
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
  bool applied_negevents_condition = false;

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
      if ((s == conv::kSKsFCFV1RingMuLike) && (oscana != conv::kAnaNuEAppOnly)) {

	applied_negevents_condition = fExp->Generate(*fFitParams, fAllBANFFFit) || applied_negevents_condition;
	if(!fGoodnessOfFit) {
	  int nbins = kNErecoBins1Rmu;
#ifdef __show_bin_contributions_to_chisq__
	  cout << "FitInputs::Chisquare()> nbins = " << nbins << endl;
#endif
	  //loop over Ereco bins
	  for(int ibin = 1; ibin <= nbins; ibin++) {
	    double nobs = fObs->ErecoSpectrumArray             (s)[ibin-1];
	    double nexp = fExp->CurrPredictedErecoSpectrumArray(s)[ibin-1];
	    
	    double delta_chisq = 0.;
	    if      (opt.Contains("poisson")) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	    else if (opt.Contains("gauss")  ) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
	    chisq += delta_chisq;

#ifdef __show_bin_contributions_to_chisq__
	    double xlo = kErecoBinEdges1Rmu[ibin-1];
	    double xup = kErecoBinEdges1Rmu[ibin];
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
	    int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins1Rmu,kErecoBinEdges1Rmu,binedges[ibin-1]);
	    int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins1Rmu,kErecoBinEdges1Rmu,binedges[ibin]);
	    double nobs = 0, nexp = 0;
	    for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	      nobs += fObs->ErecoSpectrumArray             (s)[ir-1];
	      nexp += fExp->CurrPredictedErecoSpectrumArray(s)[ibin-1];
	    }
	    double delta_chisq = 0.;
	    if      (opt.Contains("poisson")) { delta_chisq = util::ChisqPoisson(nobs,nexp); }
	    else if (opt.Contains("gauss")  ) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
	    chisq += delta_chisq;
	      
#ifdef __show_bin_contributions_to_chisq__
	    double xlo = kErecoBinEdges1Rmu[ibin-1];
	    double xup = kErecoBinEdges1Rmu[ibin];
	    cout << "FitInputs::Chisquare()> delta_chisq (bin = " << ibin 
		 << "; [" << xlo << ", " << xup << "] GeV) = " << delta_chisq << endl;
	    cout << "\tNobs: " << nobs << "\tNexp:" << nexp << endl;
#endif
	  }//ibin
	}//goodnessoffit
      }//1Rmu
      
      // Add chisq for 1Re sample
      // If we only want the numu disappearance analysis only, then skip
      else 
      if ((s == conv::kSKsFCFV1RingELike) && (oscana != conv::kAnaNuMuDisappOnly) )
      {
	    applied_negevents_condition = fExp->Generate(*fFitParams, fAllBANFFFit) || applied_negevents_condition;
	    // Note: Counting experiment for the time-being
	    double nobs = fObs->ErecoSpectrum(s)->Integral();
	    double nexp = fExp->CurrPredictedErecoSpectrum(s)->Integral();
	    
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

  //
  // Add penalty term for BANFF params 
  // First, make sure that all BANFF params were included in the fit param
  // list and were marked as nuisance parameters.
  //
  bool incl_banff_penalty = true;
  /* this check is unnecessary, as long as the fitparamlist is correctly implemented
  if(all_syst) {
    const int np_banff = conv::kNumOfParamsInBANFFfit;
    for(int ip = 0; ip < np_banff; ip++) {
      PhysParamId_t p = PhysParamId_t ( (int) kNuisBANFF00 + ip );
      incl_banff_penalty = incl_banff_penalty && 
	(fFitParams->Exists(p) && fFitParams->IsNuisance(p));
      if(!incl_banff_penalty) break;
    }
  }
  else {
    PhysParamId_t p = kNuisBANFFMAQE;
    incl_banff_penalty = incl_banff_penalty && 
      (fFitParams->Exists(p) && fFitParams->IsNuisance(p));
    p = kNuisBANFFCCQEE0;
    incl_banff_penalty = incl_banff_penalty && 
      (fFitParams->Exists(p) && fFitParams->IsNuisance(p));
  }
  */
  if(incl_banff_penalty) {
#ifdef __include_all_systematics_in_banff_penalty_term__
    double dchisq_penalty = fExp->GetBANFF()->ChisqPenalty();
#else
    double dchisq_penalty = fExp->GetBANFF()->ChisqPenalty(all_syst);
#endif
    chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
      cout << "FitInputs::ChisquarePenaltyTerm >> Current BANFF param variation"
           << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }

  //
  // Add penalty term for SuperkDetector params
  // First check all the parms considered by SuperkDetector are included in the fit
  // param list and set as nuisance
  //
  bool incl_skdet_penalty = true;
  /* this check is unnecessary, as long as the fitparamlist is correctly implemented
  if(all_syst) {
    int np_skdet = 8; // 7 eff + energy scale
    for(int ip = 0; ip < np_skdet; ip++) {
      PhysParamId_t p = PhysParamId_t ( (int) kNuisSuperkEnergyScale  + ip );
      incl_skdet_penalty = incl_skdet_penalty &&
	(fFitParams->Exists(p) && fFitParams->IsNuisance(p));
      if(!incl_skdet_penalty) break;
    }
  }
  else {
    PhysParamId_t p = kNuisSuperkEff1RmuCCnonQE;
    incl_skdet_penalty = incl_skdet_penalty &&
      (fFitParams->Exists(p) && fFitParams->IsNuisance(p));
    p = kNuisSuperkEff1RmuNC;
    incl_skdet_penalty = incl_skdet_penalty &&
      (fFitParams->Exists(p) && fFitParams->IsNuisance(p));
    p = kNuisSuperkEnergyScale;
    incl_skdet_penalty = incl_skdet_penalty &&
      (fFitParams->Exists(p) && fFitParams->IsNuisance(p));
  }
  */
  if(incl_skdet_penalty) {
    double dchisq_penalty = fExp->GetSuperkDetector()->ChisqPenalty(all_syst);
    chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
    cout << "FitInputs::ChisquarePenaltyTerm >> Current variation of params considered by SuperkDetector"
         << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }

  //
  // Add penalty term for NuXSec params (other than the ones incl. in BANFF)
  //
  bool incl_xsec_penalty = true;
  /* this check is unnecessary, as long as the fitparamlist is correctly implemented
  if(all_syst) {
    int np_xsec = 11;
    for(int ip = 0; ip < np_xsec; ip++) {
      PhysParamId_t p = PhysParamId_t ( (int) kNuisXSecFermiMom  + ip );
      incl_xsec_penalty = incl_xsec_penalty &&
	(fFitParams->Exists(p) && fFitParams->IsNuisance(p));
      if(!incl_xsec_penalty) break;
    }
  }
  else {
    PhysParamId_t p = kNuisXSecPiLessDelta;
    incl_xsec_penalty = incl_xsec_penalty &&
      (fFitParams->Exists(p) && fFitParams->IsNuisance(p));
    p = kNuisXSecSpecFunc;
    incl_xsec_penalty = incl_xsec_penalty &&
      (fFitParams->Exists(p) && fFitParams->IsNuisance(p));
    p = kNuisXSecBindE;
    incl_xsec_penalty = incl_xsec_penalty &&
      (fFitParams->Exists(p) && fFitParams->IsNuisance(p));
  }
  */
  if(incl_xsec_penalty) {
    double dchisq_penalty = fExp->GetNuXSec()->ChisqPenalty(all_syst);
    chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
    cout << "FitInputs::ChisquarePenaltyTerm >> Current variation of params considered by NuXSec"
         << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }


  //
  // Add penalty term for FSI params 
  //
  if(all_syst) {
    bool incl_fsi_penalty = true;
    /* this check is unnecessary, as long as the fitparamlist is correctly implemented
    const int np_fsi = conv::kNumOfFSIParams;
    for(int ip = 0; ip < np_fsi; ip++) {
      PhysParamId_t p = PhysParamId_t ( (int) kNuisFSI000 + ip );
      incl_fsi_penalty = incl_fsi_penalty && 
	(fFitParams->Exists(p) && fFitParams->IsNuisance(p));
      if(!incl_fsi_penalty) break;
    }
    */
    if(incl_fsi_penalty) {
      double dchisq_penalty = fExp->GetFSI()->ChisqPenalty();
      chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
      cout << "FitInputs::ChisquarePenaltyTerm >> Current FSI param variation"
           << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
    }
  }


  //
  // Add penalty term for marginalised oscillation parameters.
  // Penalize deviations from global analysis best-fit value.
  // ...
  /* don't actually need these for the current analysis
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

  if( fFitParams->Exists    (kPhypDm23sq) && 
      fFitParams->IsNuisance(kPhypDm23sq) )
  {
      double dmsq = fFitParams->CurrValue(kPhypDm23sq);
      double nom  = kGlobFitDm23sq;
      double err  = kGlobFitDm23sq * kGlobFitErrDm23sq;
      double twk  = (dmsq - nom) / err;
      double dchisq_penalty = twk*twk;
      chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
      cout << "FitInputs::ChisquarePenaltyTerm >> |dm23sq| = " << dmsq
           << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
           << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty << endl;
#endif
  }//dm23sq
  */

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

  return this->Chisquare(fOpt.c_str());
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
