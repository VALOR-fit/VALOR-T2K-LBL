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
//~ #define __show_init_mesg_debug__
#define __show_skdetfsiconv_check_mesg_info__
//~ #define __show_skdetfsiconv_check_mesg_debug__
//~ #define __show_skdetfsi_cholesky_mesg_info__
//~ #define __show_skdetfsi_cholesky_mesg_debug__
//~ #define __show_skdetfsi_tweak_mesg_debug__
//~ #define __show_chisqpenalty_mesg_debug__
//~ #define __show_skdetfsi_applysyst_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>

#include <TMath.h>

// FMWK
#include "FitCore/FitParamList.h"
#include "FitCore/PhysParamList.h"
#include "Messenger/Messenger.h"
#include "Utilities/UtlPrint.h"
#include "Utilities/UtlSystem.h"

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
  this->Init();

  // Check that SKFSISystConv.h is set up corretly, and we have interpreted it right.
  CheckListOfSKDetFSISyst(fNumOfActiveSKDetFSIParams, fActiveSKDetFSIParams,
        "kListOfSKDetFSIParams", true);

  CheckERange(fNumOfActiveSKDetFSIParams, fActiveSKDetFSIParams,
              (const double *)kSKDetFSIParamsErecoCut, "kSKDetFSIParamsErecoCut");
              
  //~ CheckSKMode(fNumOfActiveSKDetFSIParams, fActiveSKDetFSIParams,
              //~ (const bool *)kSKDetFSIParamToSKMode, "fActiveSKDetFSIParams", false);
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
  int allowed_gaps = nsamples-1;
  for(id = 1; id < nparams; id++) {
    prev = params[id-1];
    curr = params[id];
#ifdef __show_skdetfsiconv_check_mesg_debug__
    if(id == 1) {
      cout << "SuperkDetectorFSI/CheckListOfSKDetFSISyst >> " << vname << "[" << id-1 << "] = " << prev << " " << util::AsString(prev) << endl;
    }
    cout << "SuperkDetectorFSI/CheckListOfSKDetFSISyst >> " << vname << "[" << id   << "] = " << curr << " " << util::AsString(curr) << endl;
#endif
    if((int)curr != ((int)prev + 1)) {
      //~ //this is to allow jumps between sample parameters
      if(allow_gaps && allowed_gaps>0) {
        allowed_gaps--;
        continue;
      }
      cerr << "SuperkDetectorFSI/CheckLIstOfSKDetFSISyst >> " << vname << "[" << id   << "] = " << curr << " " << util::AsString(curr) << endl
           << " comes before " << vname << "[" << id-1 << "] = " << prev << " " << util::AsString(prev)
           << endl;
      util::Die(-1); return;
    }
  }//id
}
//________________________________________________________________________________________
void SuperkDetectorFSI::CheckERange(int nparams, const PhysParamId_t * params, const double * erange, const char * vname)
{
  PhysParamId_t curr, prev;
  //check that XSec parameters effect the correct energy range
  for(int ix = 1; ix < nparams; ix++) {
    prev = params[ix-1];
    curr = params[ix];
    double this_lo = erange[ix*2];
    double this_hi = erange[ix*2 + 1];
    double prev_lo = erange[(ix-1)*2];
    double prev_hi = erange[(ix-1)*2 + 1];
#ifdef __show_skdetfsiconv_check_mesg_debug__
    if(ix == 1)
      cout << "SuperkDetectorFSI/CheckERange >> " << vname << "[" << ix-1 << "] = " << prev_lo << "," << prev_hi << " " << util::AsString(prev) << endl;
    cout << "SuperkDetectorFSI/CheckERange >> " << vname << "[" << ix << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << endl;
#endif
    if(TMath::Abs(this_lo - prev_hi) > 1E-6) {
      if(TMath::Abs(prev_hi - 30) < 1E-6 && TMath::Abs(this_lo) < 1E-6) {} //it's ok! (numu modes valid 0-30)
      else if(TMath::Abs(prev_hi - 1.25) < 1E-6 && TMath::Abs(this_lo) < 1E-6) {} //it's ok! (nue modes valid 0-1.25)
      else {
        cerr << "SuperkDetectorFSI/CheckERange >> Wrong energy range" << endl
             << vname << "[" << ix-1 << "] = " << prev_lo << "," << prev_hi << " " << util::AsString(prev) << endl
             << vname << "[" << ix << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << endl
             << "Should be previous high and this low should be identical (or 30 and 0, or 1.25 and 0)" << endl;
        exit(-1);
      }
    }
    if(TMath::Abs(this_lo + 1) < 1E-6 || TMath::Abs(this_hi + 1) < 1E-6) {
      cerr << "SuperkDetectorFSI/CheckERange >> Wrong energy range" << endl
           << vname << "[" << ix << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << endl
           << "XSec parameters shouldn't have -1 here" << endl;	
    }
  }//ix
}
//________________________________________________________________________________________
void SuperkDetectorFSI::CheckSKMode(int nparams, const conv::PhysParamId_t * params, const bool * /*modes*/, const char * vname, bool fatal_no_effect)
{
  // Check this still does what we wanted.
  PhysParamId_t curr;
  
  bool found_good_sample = false;
  bool failure = false;
  //check that SK Det FSI parameters effect the correct energy range
  for(int ix = 0; ix < nparams; ix++) {
    curr = params[ix];
#if defined(__show_skdetfsiconv_check_mesg_debug__) || defined(__show_skdetfsiconv_check_mesg_info__)
    cout << "SuperkDetectorFSI/CheckSKMode >> " << vname << " Parameter " << ix << " " << util::AsString(curr) << endl;
#endif
    std::vector<int> n_modes_applied_to;
    for(int is = 0; is <conv::kns; is++) {
      n_modes_applied_to[is] = 0;
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
      found_good_sample = true;
      
#ifdef __show_skdetfsiconv_check_mesg_debug__
      cout << "\t" << util::AsString(s) << endl;
#endif
      for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        bool applies_to_mode = this->ApplySystToSKMode(curr,s,m);
        if(applies_to_mode) {
          n_modes_applied_to[is]++;
          //~ modes_set[is].push_back(m);
#ifdef __show_skdetfsiconv_check_mesg_debug__
          cout << util::AsString(m) << "\t\t";
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
      //check that they apply to all e.g. CCQE
      // this must be done by printing to screen & analyser
#ifdef __show_skdetfsiconv_check_mesg_debug__
      cout << "Effects " << util::AsString(s) << endl;
      PrintAppliesType(util::AppliesToCCQE    (modes_set[is]), "CCQE");
      PrintAppliesType(util::AppliesToCCMEC   (modes_set[is]), "CCMEC");
      PrintAppliesType(util::AppliesToCCcoh   (modes_set[is]), "CCcoh");
      PrintAppliesType(util::AppliesToCC1pi   (modes_set[is]), "CC1pi");
      PrintAppliesType(util::AppliesToCCoth   (modes_set[is]), "CCoth");
      PrintAppliesType(util::AppliesToNC1pi0  (modes_set[is]), "NC1pi0");
      PrintAppliesType(util::AppliesToNC1piC  (modes_set[is]), "NC1piC");
      PrintAppliesType(util::AppliesToNCcoh   (modes_set[is]), "NCcoh");
      PrintAppliesType(util::AppliesToNC1Gamma(modes_set[is]), "NC1Gamma");
      PrintAppliesType(util::AppliesToNCoth   (modes_set[is]), "NCoth");
#endif
    } //is
  }//ix
  
  if (failure && !fatal_no_effect) {
    util::Die(1);return;
  }
}
//________________________________________________________________________________________
void SuperkDetectorFSI::Init(void)
{
  fCurrentValues          = 0;
  fNominalValues          = 0;
  fCovMtrx                = 0;
  fCovMtrxChoL            = 0;
  fEnergyScaleRelErrMu    = 0.;
  fEnergyScaleRelErrE     = 0.;
  fEnergyScaleTwkDial     = 0.;
  
  // Decide what matrix definition we want.
  fNumOfActiveSKDetFSIParams = 0;
  
  for (int is = 0 ; is < conv::kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    switch (s) {
      case (kSKsFCFV1RingMuLike     ): fNumOfActiveSKDetFSIParams += kNumOfSKDetFSIParams1Rmu;     break;
      case (kSKsFCFV1RingELike      ): fNumOfActiveSKDetFSIParams += kNumOfSKDetFSIParams1Re;      break;
      case (kSKsRHCFCFV1RingMuLike  ): fNumOfActiveSKDetFSIParams += kNumOfRHCSKDetFSIParams1Rmu;  break;
      case (kSKsRHCFCFV1RingELike   ): fNumOfActiveSKDetFSIParams += kNumOfRHCSKDetFSIParams1Re;   break;
      case (kSKsFCFVMultiRingMuLike ): fNumOfActiveSKDetFSIParams += kNumOfSKDetFSIParamsMultiRmu; break;
      case (kSKsFCFVMultiRingELike  ): fNumOfActiveSKDetFSIParams += kNumOfSKDetFSIParamsMultiRe;  break;
      case (kSKsUndefined           ): break;
    }
  }
  
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
  fNumOfActiveSKDetFSIParams = std::accumulate(number_of_relevant_systematics,number_of_relevant_systematics+conv::kns,0);
  fActiveSKDetFSIParams = new conv::PhysParamId_t[fNumOfActiveSKDetFSIParams];
  
  // This logic needs careful checking
  // The intention here is that we iterate through our samples and append the relevant systematics in turn to our variable
  // This is complicated by the fact that our FSI arrays are defined in pairs of samples, not one indiividual sample
  int ip = 0; // Will be used to iterate our position in the fNumOfActiveSKDetFSIParams variables
  int ip_by_type[3] = {}; // individually count our three types of systematic (FHC, RHC, MR)
  for (int is = 0 ; is < conv::kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    for (int ipar_for_this_sample = 0 ; ipar_for_this_sample < number_of_relevant_systematics[is] ; ipar_for_this_sample++) {
      switch (s) {
        case (kSKsFCFV1RingMuLike     ): 
        case (kSKsFCFV1RingELike      ): 
          fActiveSKDetFSIParams[ip++] = kListOfSKDetFSIParams         [ip_by_type[0]++];  break;
        case (kSKsRHCFCFV1RingMuLike  ): 
        case (kSKsRHCFCFV1RingELike   ):
          fActiveSKDetFSIParams[ip++] = kListOfRHCSKDetFSIParams      [ip_by_type[1]++];  break;
        case (kSKsFCFVMultiRingMuLike ): 
        case (kSKsFCFVMultiRingELike  ):
          fActiveSKDetFSIParams[ip++] = kListOfMultiRingSKDetFSIParams[ip_by_type[2]++];  break;
        case (kSKsUndefined           ):
          util::Die(-1);
          return;
      }
    }
  }
  
  LOG("SuperkDetectorFSI",pINFO) << "Defined final SKDetFSI matrix rows: ";
  for (int i = 0 ; i < fNumOfActiveSKDetFSIParams; i++ ) {
    LOG("SuperkDetectorFSI",pINFO) << "\t"<<util::AsString(fActiveSKDetFSIParams[i]);
  }
  arr_index = util::ArrIdx::Instance();
}
//________________________________________________________________________________________
void SuperkDetectorFSI::CleanUp(void)
{
  if ( fCurrentValues       ) delete fCurrentValues;
  if ( fNominalValues       ) delete fNominalValues;
  if ( fCovMtrx             ) delete fCovMtrx;
  if ( fCovMtrxChoL         ) delete fCovMtrxChoL;
  if ( fActiveSKDetFSIParams) delete fActiveSKDetFSIParams;
}
//________________________________________________________________________________________
bool SuperkDetectorFSI::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "SuperkDetectorFSI >> Reading SuperkDetector+FSI inputs from input XML file." << endl;
#endif

  fEnergyScaleRelErrMu = util::XML2Dbl(xml,root_node,T2KOpt::Instance()->fPathEscale1RmuErr.c_str());
  fEnergyScaleRelErrE  = util::XML2Dbl(xml,root_node,T2KOpt::Instance()->fPathEscale1ReErr.c_str());

  const int np = fNumOfActiveSKDetFSIParams;

  fCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathSKDetFSICovMatr.c_str());    
    
  if(!fCovMtrx) {
    cerr << "SuperkDetectorFSI >> No SuperkDetector+FSI covariance matrix was found." << endl;
    exit(1);
  }

  const TMatrixD & cov = *fCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np)
  {
    cerr << "SuperkDetectorFSI >> SuperkDetector+FSI covariance matrix has unexpected size!" << endl;
    cerr << "SuperkDetectorFSI >> Expecting: " <<  np << " x " << np << " matrix" << endl;
    cerr << "SuperkDetectorFSI >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
    exit(1);
  }

  double det;
  fCovMtrxInv = new TMatrixD(*fCovMtrx);
  (*fCovMtrxInv).SetTol(1E-23);
  (*fCovMtrxInv).Invert(&det);
  //  cout << "Determinant = " << det << endl;
  const TMatrixD &   covInv  = *fCovMtrxInv;
  TMatrixD prod = cov*covInv;

#ifdef __show_init_mesg_debug__
  cout << "SuperkDetectorFSI >> " << "Error {EnergyScale;mu}    : " << fEnergyScaleRelErrMu  << endl;
  cout << "SuperkDetectorFSI >> " << "Error {EnergyScale;e}    : "  << fEnergyScaleRelErrE  << endl;
  cout << "SuperkDetectorFSI >> Input SuperkDetector+FSI covariance matrix: " << endl; 
  //cov.Print();
  util::PrintMatrix(cov,7,4);
  cout << "SuperkDetectorFSI >> Inverse SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintMatrix(*fCovMtrxInv,7,4);
  cout << "SuperkDetectorFSI >> check of product cov*(cov)^-1 : "<< endl;
  util::PrintMatrix(prod,7,4);
#endif

  // Create/init current SuperkDetector+FSI params
  fCurrentValues = new TVectorD(np);
  fNominalValues = new TVectorD(np);
  for(int ip = 0; ip < np; ip++) 
  {
     (*fCurrentValues)[ip] = 1.;
     (*fNominalValues)[ip] = 1.;
  }

#ifdef __show_init_mesg_info__
  cout << "SuperkDetectorFSI >> Done reading SuperkDetector+FSI inputs!" << endl;
#endif

  return true;  
}
//________________________________________________________________________________________                                                                                        
void SuperkDetectorFSI::TweakEnergyScale(double twkdial)
{
  fEnergyScaleTwkDial = twkdial;
}
//________________________________________________________________________________________
double SuperkDetectorFSI::SystWght(
  conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double /*Etrue*/) const
{
  // find SuperkDetector+FSI param which corresponds to the current Ereco param
  int iparam = -1;
  int np = fNumOfActiveSKDetFSIParams;
  conv::PhysParamId_t p;
  
  for(int ip = 0; ip < np; ip++) {
    p = fActiveSKDetFSIParams[ip];
    
    // Is the current systematic parameter applied to current sample/mode end energy bin?
    bool apply = this->ApplySyst(p,s,m,Ereco);
    if(!apply) continue;
    
    iparam = ip;
  } //ip
  
  if (iparam==-1) {
    LOG("SuperkDetectorFSI",pERROR) << "Failed to get SK Systematic applied for "<<util::AsString(s)<<" with "<<util::AsString(m) <<" at energy "<<Ereco<<" GeV";
  }
  
  double skdetfsi_param_value = (*fCurrentValues)[iparam];

  // param tweaked
  bool tweaked = fabs(skdetfsi_param_value - 1.) > conv::kASmallNum;
  if(!tweaked) {
    return 1.;
  }
  double skdetfsi_wght = skdetfsi_param_value;
  return skdetfsi_wght;
}
// //________________________________________________________________________________________
double SuperkDetectorFSI::CurrentValue(conv::PhysParamId_t p) const
{
  if(p == kNuisSuperkEnergyScale)
    return fEnergyScaleTwkDial;
  const TVectorD & current_values = *fCurrentValues;
  int idx = arr_index->SKDetFSISystIdx(p);
  return current_values[idx];
}
//________________________________________________________________________________________
double SuperkDetectorFSI::OneSigmaError(PhysParamId_t p) const
{
  if(p == kNuisSuperkEnergyScale)
    return 1;
  const TMatrixD & cov = *fCovMtrx;
  int idx = arr_index->SKDetFSISystIdx(p);
  if (idx == -1) {
    std::cerr<<"SuperkDetectorFSI/OneSigmaError >> Failed for "<<util::AsString(p)<<std::endl;
    return 1;
  }
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
void SuperkDetectorFSI::Tweak(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_skdetfsi_params, /* input vector of FSI params - Set to NULL to reset to nominal ones */
  Option_t * /*option*/
)
{
#ifdef __show_skdetfsi_cholesky_mesg_debug__
  cout << "SuperkDetectorFSI/Tweak >> Tweaking SuperkDetector+FSI parameters" << endl;
#endif
  
  TVectorD & current_values = *fCurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
  {
     // Input vector is NULL.
     // This resets the FSI params to their best-fit values.
     if(input_skdetfsi_params==0) {
#ifdef __show_skdetfsi_tweak_mesg_debug__
        cout << "SuperkDetectorFSI/Tweak >> Reset SuperkDetector+FSI params to nominal values." << endl;
#endif
        const int np = fNumOfActiveSKDetFSIParams;
        for(int ip = 0; ip < np; ip++) {
           current_values[ip] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_skdetfsi_tweak_mesg_debug__
        cout << "SuperkDetectorFSI/Tweak >> Forcing input SuperkDetector+FSI parameters. " << endl;
#endif
        TVectorD & vec = *input_skdetfsi_params;
        current_values = vec;
        return;
     }
  }//force_inputs?

  //
  // If this is the first time this function is called, then run the Cholesky
  // decomposition of the flux covariance matrix
  //

  if(fCovMtrxChoL == 0)
  {
#ifdef __show_skdetfsi_cholesky_mesg_info__
    cout << "SuperkDetectorFSI/Tweak >> "
         << "Running Cholesky's decomposition method on the SuperkDetector+FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_skdetfsi_cholesky_mesg_debug__
    cout << "SuperkDetectorFSI/Tweak >> SuperkDetector+FSI covariance matrix: " << endl;
    util::PrintMatrix(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_skdetfsi_cholesky_mesg_debug__
    cout << "SuperkDetectorFSI/Tweak >> Cholesky decomposition of SuperkDetector+FSI covariance matrix: " << endl;
    util::PrintMatrix(*fCovMtrxChoL,7,4);
#endif
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
  {   
#ifdef __show_skdetfsi_tweak_mesg_debug__
     cout << "SuperkDetectorFSI/Tweak >> Generating random SuperkDetector+FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD & nominal_values = *fNominalValues;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
     current_values = vec;
     return;
  }
}
//________________________________________________________________________________________
double SuperkDetectorFSI::ChisqPenalty(void) const
{
  //Vector with deviations of current values from  best fit values
  const TVectorD & bestfit_values   = *fNominalValues;
  const TVectorD & current_values   = *fCurrentValues;
  const int np = fNumOfActiveSKDetFSIParams;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "SuperkDetectorFSI/ChisqPenalty >> SuperkDetector+FSI parameters: " << endl;
  util::PrintVectorPair(current_values, bestfit_values, fActiveSKDetFSIParams, 12,6,"SKDetFSI fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in SuperkDetector+FSI parameters: " << endl;
  util::PrintVector(TMatrixDRow(dev,0), fActiveSKDetFSIParams, 12,6,"SKDetFSI fit params");
#endif

  //covariance matrix
  const TMatrixD cov    (this->CovMatrix());
  const TMatrixD covInv (this->CovMatrixInv());

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "SuperkDetectorFSI/ChisqPenalty >> SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintMatrix(cov,7,4);
  cout << "SuperkDetectorFSI/ChisqPenalty >> Inverse SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintMatrix(covInv,7,4);
#endif

  //Check the inverse:
  TMatrixD prod = cov*covInv;
#ifdef __show_chisqpenalty_mesg_debug__
  cout << "SuperkDetectorFSI/ChisqPenalty >> check of product cov*(cov)^-1 : "<< endl;
  util::PrintMatrix(prod,7,4);
#endif

  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

  //now add the contribution from energy scale
  if (fEnergyScaleTwkDial   ) chisq_penalty += fEnergyScaleTwkDial*fEnergyScaleTwkDial;

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "SuperkDetectorFSI/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;
#endif

  if(chisq_penalty < 0) 
    cout << "SuperkDetectorFSI/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;

  if(chisq_penalty) return chisq_penalty;
  
  return 0;
}
//________________________________________________________________________________________
void SuperkDetectorFSI::PrintSystematicVariations(std::ostream & stream, Option_t * /*opt*/) const
{
  const TVectorD & curr_skdetfsi_params = *fCurrentValues;
  stream << "SuperkDetectorFSI >> Energy scale twk        = "  << fEnergyScaleTwkDial      << endl;
  stream << "SuperkDetectorFSI >> SuperkDetector+FSI parameters: " << endl;
  util::PrintVector(curr_skdetfsi_params,fActiveSKDetFSIParams,15,6,"SKDetFSI params",stream);
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
  /// NEW MULTIRING 2 
  int idx;
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
  string applystr;
  applystr = apply ? "Passed" : "Failed";
  cout << applystr << "\t" << util::AsString(p) << "\t" << util::AsString(s) << "\t" << util::AsString(m) << endl;
#endif
  return apply;
}
//________________________________________________________________________________________
bool SuperkDetectorFSI::ApplySystToEnergyBin(conv::PhysParamId_t p, double Ereco) const
{
  int idx;
  
  bool isRHC = util::IsRHCSystematic(p);
  /// NEW MULTIRING 2
  bool isFHC = util::IsFHCSystematic(p);
  ///

  // OLD
  // if(isRHC)
  //   idx = arr_index->RHCSKDetFSISystIdx(p);
  // else
  //   idx = arr_index->SKDetFSISystIdx(p);

  if(isRHC)
    idx = arr_index->RHCSKDetFSISystIdx(p);
  else if(isFHC)
    idx = arr_index->SKDetFSISystIdx(p);
  /// NEW MULTIRING 2
  else
    idx = arr_index->MultiRingSKDetFSISystIdx(p);
  ///

  if (idx == -1) return false;

  bool apply;

  // OLD
  // if(isRHC)
  //   apply = Ereco >= conv::kRHCSKDetFSIParamsErecoCut[idx][0] &&
  //     Ereco <  conv::kRHCSKDetFSIParamsErecoCut[idx][1];
  // else
  //   apply = Ereco >= conv::kSKDetFSIParamsErecoCut[idx][0] &&
  //     Ereco <  conv::kSKDetFSIParamsErecoCut[idx][1];
  
  if(isRHC)
    apply = Ereco >= conv::kRHCSKDetFSIParamsErecoCut[idx][0] &&
      Ereco <  conv::kRHCSKDetFSIParamsErecoCut[idx][1];
  else if(isFHC)
    apply = Ereco >= conv::kSKDetFSIParamsErecoCut[idx][0] &&
      Ereco <  conv::kSKDetFSIParamsErecoCut[idx][1];
  /// NEW MULTIRING 2
  else
    apply = Ereco >= conv::kMultiRingSKDetFSIParamsErecoCut[idx][0] &&
      Ereco <  conv::kMultiRingSKDetFSIParamsErecoCut[idx][1];
  ///

#ifdef __show_skdetfsi_applysyst_mesg_debug__
  string applystr;
  applystr = apply ? "Passed" : "Failed";
  cout << applystr << "\t" << util::AsString(p) << "\t" << Ereco << endl;
#endif
  return apply;
}
//________________________________________________________________________________________
double SuperkDetectorFSI::EnergyScaleRelErr(SuperkSample_t s, SuperkMode_t /* m */) const
{
  return this->EnergyScaleRelErr(s);
}
//________________________________________________________________________________________
double SuperkDetectorFSI::EnergyScaleRelErr(SuperkSample_t s) const
{
  // energy scale relative error
  switch(s) {
    case(kSKsFCFV1RingMuLike):
    case(kSKsRHCFCFV1RingMuLike):
    case(kSKsFCFVMultiRingMuLike):
      return fEnergyScaleRelErrMu;
    case(kSKsFCFV1RingELike):
    case(kSKsRHCFCFV1RingELike):
    case(kSKsFCFVMultiRingELike):
      return fEnergyScaleRelErrE;
    case(kSKsUndefined):
      cerr << "EnergyScaleRelErr --> No sample is identified!" << endl;
      util::Die(-1); return -1;
  }
  return -1;
}
//________________________________________________________________________________________                                                                                       
double SuperkDetectorFSI::EnergyScale(SuperkSample_t s, SuperkMode_t m) const
{
  double twkdial = fEnergyScaleTwkDial;
  // energy scale
  double err     = this->EnergyScaleRelErr(s,m);
  double scale   = 1. + twkdial * err;
  return scale;
}
//________________________________________________________________________________________                                                                                       
double SuperkDetectorFSI::EnergyScale(SuperkSample_t s) const
{
  double twkdial = fEnergyScaleTwkDial;
  // energy scale
  double err     = this->EnergyScaleRelErr(s);
  double scale   = 1. + twkdial * err;
  return scale;
}
//________________________________________________________________________________________                                                                                        
void SuperkDetectorFSI::CopyTweaksTo(valor::PhysParamList & ppl) const
{
  for (int i = 0 ; i < fCurrentValues->GetNoElements() ; i++) {
    PhysParamId_t p = fActiveSKDetFSIParams[i];
    ppl.SetValue(p,(*fCurrentValues)[i]);
  }
  
  PhysParamId_t p = kNuisSuperkEnergyScale;
  
  ppl.SetValue(p,fEnergyScaleTwkDial);
  
  return;
}
//________________________________________________________________________________________
void SuperkDetectorFSI::CopyTweaksTo(valor::FitParamList & fpl) const
{
  for (int i = 0 ; i < fCurrentValues->GetNoElements() ; i++) {
    PhysParamId_t p = fActiveSKDetFSIParams[i];
    if (fpl.Exists(p)) {
      fpl.SetCurrentValue(p,(*fCurrentValues)[i]);
    }
    else {
      LOG("SuperkDetectorFSI",pWARN) 
        << "FPL does not contain parameter " << util::AsString(p);
    }
  }
  
  PhysParamId_t p = kNuisSuperkEnergyScale;
  if (fpl.Exists(p)) {
    fpl.SetCurrentValue(p,fEnergyScaleTwkDial);
  }
  else {
    LOG("SuperkDetectorFSI",pWARN) 
      << "FPL does not contain parameter " << util::AsString(p);
  }
    
  return;
}
//________________________________________________________________________________________
void SuperkDetectorFSI::UpdateFrom(const valor::FitParamList  & fpl)
{
  for (int i = 0 ; i < fNumOfActiveSKDetFSIParams ; i++) {
    PhysParamId_t p = fActiveSKDetFSIParams[i];
    if (fpl.Exists(p)) {
      (*fCurrentValues)[i] = fpl.CurrValue(p);
    }
    else {
      LOG("SuperkDetectorFSI",pWARN) 
        << "Unable to update flux parameter " << util::AsString(p) << " from FPL.";
    }
  }
  
  PhysParamId_t p = kNuisSuperkEnergyScale;
  fEnergyScaleTwkDial = fpl.CurrValue(p);
  //~ fIsTweaked = true;
  return;
}
//________________________________________________________________________________________
void SuperkDetectorFSI::SetCurrentValue(PhysParamId_t p, double tweak)
{
  if (p==kNuisSuperkEnergyScale) {
    fEnergyScaleTwkDial = tweak; return;
  }
  for (int i = 0 ; i < fNumOfActiveSKDetFSIParams ; i++) {
    if (p==fActiveSKDetFSIParams[i]) {
      (*fCurrentValues)[i] = tweak; return;
    }
  }
  LOG("SuperkDetectorFSI",pWARN)<<"Failed to update SKDetFSI matrix for unused parameter "<<util::AsString(p);
  return;
}
//________________________________________________________________________________________
