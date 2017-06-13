
//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 30-12-2013  Nick
 Added samples for reverse horn current running

 21-7-2014 Davide
 Added samples for multiring forward horn current

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
// #define __show_init_mesg_debug__
#define __show_skdetfsiconv_check_mesg_info__
// #define __show_skdetfsiconv_check_mesg_debug__
// #define __show_skdetfsi_cholesky_mesg_info__
// #define __show_skdetfsi_cholesky_mesg_debug__
// #define __show_skdetfsi_tweak_mesg_debug__
// #define __show_chisqpenalty_mesg_debug__
// #define __show_skdetfsi_applysyst_mesg_debug__
// #define __use_1Rsyst_for_MultiR_sample__


#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>

#include <TMath.h>

// FMWK
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Registry/Registry.h"
#include "valor/Utilities/UtlPrint.h"
#include "valor/Utilities/UtlSystem.h"

#include "SuperkDetectorFSI.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

using namespace valor;
using namespace valor::conv;
using namespace valor::util;

ClassImp(SuperkDetectorFSI)

//________________________________________________________________________________________
SuperkDetectorFSI::SuperkDetectorFSI()
{
  
  arr_index = ArrIdx::Instance();
  
  // Decide what matrix definition we want.
  unsigned int numOfActiveSKDetFSIParams = 0;
  
  for (int is = 0 ; is < conv::kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    switch (s) {
      case (kSKsFCFV1RingMuLike     ): numOfActiveSKDetFSIParams += kNumOfSKDetFSIParams1Rmu;     break;
      case (kSKsFCFV1RingELike      ): numOfActiveSKDetFSIParams += kNumOfSKDetFSIParams1Re;      break;
      case (kSKsRHCFCFV1RingMuLike  ): numOfActiveSKDetFSIParams += kNumOfRHCSKDetFSIParams1Rmu;  break;
      case (kSKsRHCFCFV1RingELike   ): numOfActiveSKDetFSIParams += kNumOfRHCSKDetFSIParams1Re;   break;
      case (kSKsFCFVMultiRingMuLike ): numOfActiveSKDetFSIParams += kNumOfSKDetFSIParamsMultiRmu; break;
      case (kSKsFCFVMultiRingELike  ): numOfActiveSKDetFSIParams += kNumOfSKDetFSIParamsMultiRe;  break;
      case (kSKsUndefined           ): break;
    }
  }

  LOG("SuperkDetectorFSI",pINFO) << "Number of Active SKDetFSI parameters: " << numOfActiveSKDetFSIParams;
  
  int * number_of_relevant_systematics =  new int[conv::kns]();
  
  for (int is = 0 ; is < conv::kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    switch (s) {
      case (kSKsFCFV1RingMuLike     ): number_of_relevant_systematics[is] += kNumOfSKDetFSIParams1Rmu;     break;
      case (kSKsFCFV1RingELike      ): number_of_relevant_systematics[is] += kNumOfSKDetFSIParams1Re;      break;
      case (kSKsRHCFCFV1RingMuLike  ): number_of_relevant_systematics[is] += kNumOfRHCSKDetFSIParams1Rmu;  break;
      case (kSKsRHCFCFV1RingELike   ): number_of_relevant_systematics[is] += kNumOfRHCSKDetFSIParams1Re;   break;
      case (kSKsFCFVMultiRingMuLike ): number_of_relevant_systematics[is] += kNumOfSKDetFSIParamsMultiRmu; break;
      case (kSKsFCFVMultiRingELike  ): number_of_relevant_systematics[is] += kNumOfSKDetFSIParamsMultiRe;  break;
      case (kSKsUndefined           ): break;
    }
  }

  for (int is = 0 ; is < conv::kns; is++) {
    LOG("SuperkDetectorFSI",pINFO) << "Number of Relevant SKDetFSI systematics: " << number_of_relevant_systematics[is];
  }
  
  numOfActiveSKDetFSIParams = std::accumulate(number_of_relevant_systematics,number_of_relevant_systematics+conv::kns,0);
  numOfActiveSKDetFSIParams++; // Energy scale
  
  LOG("SuperkDetectorFSI",pINFO) << "Number of Active SKDetFSI Parameters: " << numOfActiveSKDetFSIParams;
  
  PhysParamId_t activeSKDetFSIParams[numOfActiveSKDetFSIParams];
  
  // This logic needs careful checking. It should probably be replaced entirely by getting rid of the separate ListOf... arrays.
  // The intention here is that we iterate through our samples and append the relevant systematics in turn to our variable
  // This is complicated by the fact that our FSI arrays are defined in pairs of samples, not one individual sample
  int ip = 0; // Will be used to iterate our position in the numOfActiveSKDetFSIParams variables
  int ip_by_type[3] = {}; // individually count our three types of systematic (FHC, RHC, MR)
  for (int is = 0 ; is < conv::kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    for (int ipar_for_this_sample = 0 ; ipar_for_this_sample < number_of_relevant_systematics[is] ; ipar_for_this_sample++) {
      switch (s) {
        case (kSKsFCFV1RingMuLike     ): 
        case (kSKsFCFV1RingELike      ): 
          activeSKDetFSIParams[ip++] = kListOfSKDetFSIParams         [ip_by_type[0]++];  break;
        case (kSKsRHCFCFV1RingMuLike  ): 
        case (kSKsRHCFCFV1RingELike   ):
          activeSKDetFSIParams[ip++] = kListOfRHCSKDetFSIParams      [ip_by_type[1]++];  break;
        case (kSKsFCFVMultiRingMuLike ): 
        case (kSKsFCFVMultiRingELike  ):
          activeSKDetFSIParams[ip++] = kListOfMultiRingSKDetFSIParams[ip_by_type[2]++];  break;
        case (kSKsUndefined           ):
          util::Die(-1);
          return;
      }
    }
  }
  
  // Add energy scale
  activeSKDetFSIParams[numOfActiveSKDetFSIParams-1] = kNuisSuperkEnergyScale;

  // Check that SKFSISystConv.h is set up corretly, and we have interpreted it right.
  CheckListOfSKDetFSISyst(numOfActiveSKDetFSIParams, activeSKDetFSIParams,
        "kListOfSKDetFSIParams", true);

  LOG("SuperkDetectorFSI",pINFO) << "Checking ERange: ";
  
  // Need a minus one here since these checks don't apply to Escale
  // Need to get the num of params that applies to the sample and the sample array and name since they are separate

  // CheckERange(numOfActiveSKDetFSIParams-1, activeSKDetFSIParams,
  // 	      (const double *)kSKDetFSIParamsErecoCut, "kSKDetFSIParamsErecoCut");

  if (kNumOfSKDetFSIParams){
  CheckERange(kNumOfSKDetFSIParams, kListOfSKDetFSIParams,
  	      (const double *)kSKDetFSIParamsErecoCut, "kSKDetFSIParamsErecoCut");
  }
  if (kNumOfRHCSKDetFSIParams){
  CheckERange(kNumOfRHCSKDetFSIParams, kListOfRHCSKDetFSIParams,
  	      (const double *)kRHCSKDetFSIParamsErecoCut, "kRHCSKDetFSIParamsErecoCut");
  }
  if (kNumOfSKDetFSIParamsMultiRing){
  CheckERange(kNumOfSKDetFSIParamsMultiRing, kListOfMultiRingSKDetFSIParams,
  	      (const double *)kMultiRingSKDetFSIParamsErecoCut, "kMultiRingSKDetFSIParamsErecoCut");
  }
  
  CheckSKMode(numOfActiveSKDetFSIParams-1, activeSKDetFSIParams,
              (const bool *)kSKDetFSIParamToSKMode, "activeSKDetFSIParams", false);

  LOG("SuperkDetectorFSI",pINFO) << "Active parameters: ";
  for (unsigned int i = 0; i < numOfActiveSKDetFSIParams; i++){
    cout << activeSKDetFSIParams[i] << " " << util::AsString(activeSKDetFSIParams[i]) << "\n";
  }
  
  this->Init(numOfActiveSKDetFSIParams, activeSKDetFSIParams, "syst/skdetfsi");
  this->SetName("SuperkDetectorFSI");
  
  (*fDefaultValues)[this->Param2Idx(conv::kNuisSuperkEnergyScale)] = 0;//0; p-theta set to one, event discrepancies during validation 
  
  LOG("SuperkDetectorFSI",pINFO) << "Defined final SKDetFSI matrix rows: ";
  for (unsigned int i = 0 ; i < fNSystParameters; i++ ) {
    LOG("SuperkDetectorFSI",pINFO) << "\t" << fSystParameters[i] << " " <<util::AsString(fSystParameters[i])<<" (default: "<<(*fDefaultValues)[i]<<")";
  }
  
  Registry & reg = Opt::Instance()->MutableGlobalRegistry();
  reg.Set("skdetfsi", this);
}
//________________________________________________________________________________________
SuperkDetectorFSI::~SuperkDetectorFSI()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkDetectorFSI::CheckListOfSKDetFSISyst(int nparams, const conv::PhysParamId_t * params, const char * vname, bool allow_gaps)
{
  int id = 0;
  PhysParamId_t curr, prev;

  // check that SKDETFSI parameters are monotonically increasing
  // Nsamples-1 gaps are now allowed due to FHC/RHC/MR each having individual numu/nue options.
  int nsamples=0;
  for (int is = 0 ; is < conv::kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    if (s!=conv::kSKsUndefined) {
      nsamples++;
    }
  }
  int allowed_gaps = 2*nsamples;//nsamples-1;
  for(id = 1; id < nparams; id++) {
    prev = params[id-1];
    curr = params[id];
#ifdef __show_skdetfsiconv_check_mesg_debug__
    if(id == 1) {
      LOG("SuperkDetectorFSI", pINFO) << vname << "[" << id-1 << "] = " << prev << " " << util::AsString(prev);
    }
    LOG("SuperkDetectorFSI", pINFO) << vname << "[" << id << "] = " << curr << " " << util::AsString(curr);
#endif

    if((int)curr != ((int)prev + 1)) {
      //~ //this is to allow jumps between sample parameters
      if(allow_gaps && allowed_gaps>0) {
        allowed_gaps--;
	LOG("SuperkDetectorFSI",pDEBUG) << allowed_gaps;
        continue;
      }
      LOG("SuperkDetectorFSI", pERROR) << vname << "[" << id << "] = " << curr << " " << util::AsString(curr) << "\n"
			      << " comes before " << vname << "[" << id-1 << "] = " << prev << " " << util::AsString(prev);
      util::Die(-1); return;
    }
  }//id
}
//________________________________________________________________________________________
void SuperkDetectorFSI::CheckERange(int nparams, const PhysParamId_t * params, const double * erange, const char * vname)
{
  PhysParamId_t curr, prev;
  //check that XSec parameters affect the correct energy range
  for(int ix = 1; ix < nparams; ix++) {
    prev = params[ix-1];
    curr = params[ix];
    double this_lo = erange[ix*2];
    double this_hi = erange[ix*2 + 1];
    double prev_lo = erange[(ix-1)*2];
    double prev_hi = erange[(ix-1)*2 + 1];
#ifdef __show_skdetfsiconv_check_mesg_debug__
    if(ix == 1){
      LOG("SuperkDetectorFSI",pINFO) << vname << "[" << ix-1 << "] = " << prev << " " << prev_lo << "," << prev_hi << " " << util::AsString(prev) ;
    }
    LOG("SuperkDetectorFSI",pINFO)   << vname << "[" << ix << "] = " << curr << " " << this_lo << "," << this_hi << " " << util::AsString(curr) ;
#endif
    if(TMath::Abs(this_lo - prev_hi) > 1E-6) {
      if(TMath::Abs(prev_hi - 30) < 1E-6 && TMath::Abs(this_lo) < 1E-6) {} //it's ok! (numu modes valid 0-30)
      else if(TMath::Abs(prev_hi - 1.25) < 1E-6 && TMath::Abs(this_lo) < 1E-6) {} //it's ok! (nue modes valid 0-1.25)
      else if(TMath::Abs(prev_hi - 1.25 ) < 1E-6 && TMath::Abs(this_lo - 0.3) < 1E-6) {} // allow gaps in energy bins for the CC nue 1pi sample
      else {
	LOG("SuperkDetectorFSI",pERROR) <<  "Wrong energy range " << "\n"
			       << vname << "[" << ix-1 << "] = " << prev << " " << prev_lo << "," << prev_hi << " " << util::AsString(prev) << "\n"
			       << vname << "[" << ix << "] = " << curr << " " << this_lo << "," << this_hi << " " << util::AsString(curr) << "\n"
			       << "Should be previous high and this low should be identical (or 30 and 0, or 1.25 and 0)";
	util::Die(-1);
      }
    }
    if(TMath::Abs(this_lo + 1) < 1E-6 || TMath::Abs(this_hi + 1) < 1E-6) {
      LOG("SuperkDetectorFSI",pERROR) << "Wrong energy range " << curr << "\n"
			     << vname << "[" << ix << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << "\n"
			     << "XSec parameters shouldn't have -1 here";	
    }
  }//ix
}
//________________________________________________________________________________________
void SuperkDetectorFSI::CheckSKMode(int nparams, const conv::PhysParamId_t * params, const bool * /*modes*/, const char * vname, bool fatal_no_effect)
{
  bool found_good_sample = false;
  bool failure = false;
  //check that SK Det FSI parameters effect the correct energy range
  for(int ix = 0; ix < nparams; ix++) {
    PhysParamId_t curr = params[ix];
#if defined(__show_skdetfsiconv_check_mesg_debug__) || defined(__show_skdetfsiconv_check_mesg_info__)
    LOG("SuperkDetectorFSI",pINFO) << vname << " Parameter " << ix << " " << util::AsString(curr);
#endif
    std::vector<int> n_modes_applied_to(10);
    for(int is = 0; is <conv::kns; is++) {
      n_modes_applied_to[is] = 0;
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
      found_good_sample = true;
      
#ifdef __show_skdetfsiconv_check_mesg_debug__
      LOG("SuperkDetectorFSI",pINFO) << "\t" << util::AsString(s) << "\n";
#endif
      for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        bool applies_to_mode = this->ApplySystToSKMode(curr,s,m);
        if(applies_to_mode) {
          n_modes_applied_to[is]++;
#ifdef __show_skdetfsiconv_check_mesg_debug__
         LOG("SuperkDetectorFSI",pINFO) << util::AsString(m) << "\t\t";
#endif
        }
      }//im
#ifdef __show_skdetfsiconv_check_mesg_debug__
      cout << endl;
#endif
    }//is
    if(!found_good_sample) {
#ifdef __show_skdetfsiconv_check_mesg_info__
      cout << "SuperkDetectorFSI/CheckSKMode >> Relevant samples not active. Not performing " << vname << " check" << endl;
#endif
      return;
    }
    
    SuperkSample_t prev_sample_applied_to = conv::kSKsUndefined;
    for(int is = 0; is <conv::kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
      
      if ((n_modes_applied_to[is]) > 0) {
        if (prev_sample_applied_to!=conv::kSKsUndefined) {
          LOG("SuperkDetectorFSI",pERROR) << "Systematic "<<util::AsString(curr)<<" applies to both sample "<<util::AsString(prev_sample_applied_to)<<" and "<<util::AsString(s);
          failure = true;
        }
        prev_sample_applied_to = s;
      }
    } //is
  }//ix
  
  if (failure && !fatal_no_effect) {
    util::Die(1);return;
  }
}
//________________________________________________________________________________________
double SuperkDetectorFSI::SystWght(
  conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double /*Etrue*/) const
{
#ifdef __use_1Rsyst_for_MultiR_sample__
  // ~#{ Hack for testing set multiring splines to single ring splines
  if (util::IsMultiRMuSample(s)) s = kSKsFCFV1RingMuLike;
  if (util::IsMultiReSample(s))  s = kSKsFCFV1RingELike;
  // }#~
#endif


  //util::PrintMatrix(CovMatrix(), 7, 5);

  // find SuperkDetector+FSI param which corresponds to the current Ereco param
  int iparam = -1;
  int np = fNSystParameters;
  conv::PhysParamId_t p;
  
  for(int ip = 0; ip < np; ip++) {
    p = fSystParameters[ip];
    
    // Is the current systematic parameter applied to current sample/mode end energy bin?
    bool apply = this->ApplySyst(p,s,m,Ereco);
    if(apply) {iparam = ip; break;}
    //~# if(!apply) continue;
    
    //iparam = ip;
  } //ip
  
  if (iparam==-1) {
    LOG("SuperkDetectorFSI", pERROR) << "Failed to get SK Systematic applied for "<<util::AsString(s)<<" with "<<util::AsString(m) <<" at energy "<<Ereco<<" GeV";
  }
  //LOG("SuperkDetectorFSI",pINFO) << "Current values for SKDetFSI";
  //~# fCurrentValues->Print();
  //~# fUntunedValues->Print();
  //~# fInitialValues->Print();

  double skdetfsi_param_value = (*fCurrentValues)[iparam];

  // param tweaked
  bool tweaked = fabs(skdetfsi_param_value - 1.) > conv::kASmallNum;
  if(!tweaked) {
    return 1.;
  }
  double skdetfsi_wght = skdetfsi_param_value;
/*  if(p >= kNuisSKDetFSI006 && p <= kNuisSKDetFSI017)
  {  // Fix an SK weight for testing
     skdetfsi_wght = 1.0;
  }*/
  
  return skdetfsi_wght;
}
//________________________________________________________________________________________
bool SuperkDetectorFSI::ApplySyst(
  conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const
{
  bool apply2mode = this->ApplySystToSKMode(p,s,m);
  if (not apply2mode) return false;
  bool apply2bin  = this->ApplySystToEnergyBin(p,Ereco);
  
  return apply2bin;
}
//________________________________________________________________________________________
bool SuperkDetectorFSI::ApplySystToSKMode(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const
{
  if (p==conv::kNuisSuperkEnergyScale) {return true;}
  int idx = -999;
  bool apply = false;  

  int is = arr_index->SuperkSampleIdx(s);
  if (is == -1) return false;
  int im = arr_index->SuperkModeIdx(s,m);
  if (im == -1) return false;
 
  int idxRing = s % 2;  // 0 --> RingMu, 1--> RingE
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike) {
    idx = arr_index->RHCSKDetFSISystIdx(p);
    if (idx == -1) return false;
    apply = conv::kSKDetFSIParamToSKMode[idxRing][im][idx];             // FULLY CORRELATED RHC FOR NOW;
  }
  else if(s == kSKsFCFVMultiRingMuLike || s == kSKsFCFVMultiRingELike) {
    idx = arr_index->MultiRingSKDetFSISystIdx(p);
    if (idx == -1) return false;
    apply = conv::kMultiRingSKDetFSIParamToSKMode[idxRing][im][idx]; 
  }
  else {
    idx = arr_index->SKDetFSISystIdx(p);
    if (idx == -1) return false;
    apply = conv::kSKDetFSIParamToSKMode[idxRing][im][idx]; 
  }
  ///


#ifdef __show_skdetfsi_applysyst_mesg_debug__
  // string applystr;
  // applystr = apply ? "Applies" : "DoesNotApply";
  // LOG("",pINFO) << applystr << " " << util::AsString(p) << " " << util::AsString(s) << " " << util::AsString(m);
#endif
  return apply;
}
//________________________________________________________________________________________
bool SuperkDetectorFSI::ApplySystToEnergyBin(conv::PhysParamId_t p, double Ereco) const
{
  if (p==conv::kNuisSuperkEnergyScale) {return true;}
  
  int idx;
  
  bool isRHC = util::IsRHCSystematic(p);
  bool isFHC = util::IsFHCSystematic(p);
  bool isMRE = false;

  if(isRHC)
    idx = arr_index->RHCSKDetFSISystIdx(p);
  else if(isFHC)
    idx = arr_index->SKDetFSISystIdx(p);
  else{
    idx = arr_index->MultiRingSKDetFSISystIdx(p);
    isMRE = true;
  }
  ///

  if (idx == -1) return false;

  bool apply;
  
  if(isRHC)
    apply = Ereco >= conv::kRHCSKDetFSIParamsErecoCut[idx][0] &&
      Ereco <  conv::kRHCSKDetFSIParamsErecoCut[idx][1];
  else if(isFHC)
    apply = Ereco >= conv::kSKDetFSIParamsErecoCut[idx][0] &&
      Ereco <  conv::kSKDetFSIParamsErecoCut[idx][1];
  else
    apply = Ereco >= conv::kMultiRingSKDetFSIParamsErecoCut[idx][0] &&
      Ereco <  conv::kMultiRingSKDetFSIParamsErecoCut[idx][1];
  ///

#ifdef __show_skdetfsi_applysyst_mesg_debug__
  if (Ereco < 0.4 && isMRE){
  string applystr;
  applystr = apply ? "SysBinPassed" : "SysBinFailed";
  LOG("SuperkDetectorFSI",pINFO) << applystr << "\t" << util::AsString(p) << "\t" << Ereco;
  }
#endif
  return apply;
}
//________________________________________________________________________________________
double SuperkDetectorFSI::EnergyScale(conv::SuperkSample_t /*s*/) const
{
  PhysParamId_t p = conv::kNuisSuperkEnergyScale;
  double twk     = this->CurrentValue(p);
  double central = 1.;
  return central + twk;
}
//________________________________________________________________________________________
