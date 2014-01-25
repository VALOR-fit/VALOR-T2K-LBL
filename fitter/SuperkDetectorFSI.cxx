//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 30-12-2013  Nick
 Added samples for reverse horn current running

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_init_mesg_debug__
//#define __show_skdetfsi_cholesky_mesg_info__
//#define __show_skdetfsi_cholesky_mesg_debug__
//#define __show_skdetfsi_tweak_mesg_debug__
//#define __show_chisqpenalty_mesg_debug__
//#define __show_skdetfsi_applysyst_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TMath.h>

#include "SuperkDetectorFSI.h"
#include "Utilities.h"
#include "Opt.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(SuperkDetectorFSI)

//________________________________________________________________________________________
SuperkDetectorFSI::SuperkDetectorFSI()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkDetectorFSI::~SuperkDetectorFSI()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkDetectorFSI::Init(void)
{
  fCurrentValues = 0;
  fNominalValues = 0;
  fCovMtrx       = 0;
  fCovMtrxChoL   = 0;
  fEnergyScale1RmuRelErr      = 0.;
  fEnergyScale1ReRelErr       = 0.;
  fEnergyScaleTwkDial         = 0.;
  fRHCCurrentValues = 0;
  fRHCNominalValues = 0;
  fRHCCovMtrx       = 0;
  fRHCCovMtrxChoL   = 0;
  fRHCEnergyScale1RmuRelErr      = 0.;
  fRHCEnergyScale1ReRelErr       = 0.;
  fRHCEnergyScaleTwkDial         = 0.;
  
  arr_index = util::ArrIdx::Instance();
}
//________________________________________________________________________________________
void SuperkDetectorFSI::CleanUp(void)
{
  if ( fCurrentValues ) delete fCurrentValues;
  if ( fNominalValues ) delete fNominalValues;
  if ( fCovMtrx       ) delete fCovMtrx;
  if ( fCovMtrxChoL   ) delete fCovMtrxChoL;
  if ( fRHCCurrentValues ) delete fRHCCurrentValues;
  if ( fRHCNominalValues ) delete fRHCNominalValues;
  if ( fRHCCovMtrx       ) delete fRHCCovMtrx;
  if ( fRHCCovMtrxChoL   ) delete fRHCCovMtrxChoL;
}
//________________________________________________________________________________________
bool SuperkDetectorFSI::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/SuperkDetectorFSI >> Reading SuperkDetector+FSI inputs from input XML file." << endl;
#endif

  // Get relative error on energy scale 
  fEnergyScale1RmuRelErr =
    util::XML2Dbl(xml,root_node,
		  "detector_systematics/superk/energy_scale_relative_error/one_mu_ring");
  fEnergyScale1ReRelErr =
    util::XML2Dbl(xml,root_node,
		  "detector_systematics/superk/energy_scale_relative_error/one_e_ring");
  

  const int np = conv::kNumOfSKDetFSIParams;

  // Get SuperkDetectorFSI covariance matrix
#if FSI_BINS == 1
  fCovMtrx= util::XML2TMatrixD(xml,root_node,"skdetfsi/covariance_matrix");
#elif FSI_BINS == 2
  fCovMtrx= util::XML2TMatrixD(xml,root_node,"skdetfsi/disap/covariance_matrix");
#elif FSI_BINS == 3
  fCovMtrx= util::XML2TMatrixD(xml,root_node,"skdetfsi/app/covariance_matrix");
#endif
  if(!fCovMtrx) {
    cerr << "Exp/SuperkDetectorFSI >> No SuperkDetector+FSI covariance matrix was found." << endl;
    exit(1);
  }

  const TMatrixD & cov = *fCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np)
  {
    cerr << "Exp/SuperkDetectorFSI >> SuperkDetector+FSI covariance matrix has unexpected size!" << endl;
    cerr << "Exp/SuperkDetectorFSI >> Expecting: " <<  np << " x " << np << " matrix" << endl;
    cerr << "Exp/SuperkDetectorFSI >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
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
  cout << "Exp/SuperkDetectorFSI >> " << "Error {EnergyScale;1Rmu}    : " << fEnergyScale1RmuRelErr  << endl;
  cout << "Exp/SuperkDetectorFSI >> " << "Error {EnergyScale;1Re}    : " << fEnergyScale1ReRelErr  << endl;
  cout << "Exp/SuperkDetectorFSI >> Input SuperkDetector+FSI covariance matrix: " << endl; 
  //cov.Print();
  util::PrintCompact(cov,7,4);
  cout << "Exp/SuperkDetectorFSI >> Inverse SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintCompact(*fCovMtrxInv,7,4);
  cout << "Exp/SuperkDetectorFSI >> check of product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,4);
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
  cout << "Exp/SuperkDetectorFSI >> Done reading SuperkDetector+FSI inputs!" << endl;
#endif

  return true;  
}
//________________________________________________________________________________________
bool SuperkDetectorFSI::ReadRHC(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/SuperkDetectorFSI >> Reading RHC SuperkDetector+FSI inputs from input XML file." << endl;
#endif

  // Get relative error on energy scale 
  fRHCEnergyScale1RmuRelErr =
    util::XML2Dbl(xml,root_node,
                  "detector_systematics/superk/energy_scale_relative_error/one_mu_ring");
  fRHCEnergyScale1ReRelErr =
    util::XML2Dbl(xml,root_node,
                  "detector_systematics/superk/energy_scale_relative_error/one_e_ring");
  

  const int np = conv::kNumOfRHCSKDetFSIParams;

  // Get SuperkDetectorFSI covariance matrix
#if FSI_BINS == 1
  fRHCCovMtrx= util::XML2TMatrixD(xml,root_node,"skdetfsi/covariance_matrix");
#elif FSI_BINS == 2
  fRHCCovMtrx= util::XML2TMatrixD(xml,root_node,"skdetfsi/disap/covariance_matrix");
#elif FSI_BINS == 3
  fRHCCovMtrx= util::XML2TMatrixD(xml,root_node,"skdetfsi/app/covariance_matrix");
#endif
  if(!fRHCCovMtrx) {
    cerr << "Exp/SuperkDetectorFSI >> No RHC SuperkDetector+FSI covariance matrix was found." << endl;
    exit(1);
  }

  const TMatrixD & cov = *fRHCCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np)
    {
      cerr << "Exp/SuperkDetectorFSI >> RHC SuperkDetector+FSI covariance matrix has unexpected size!" << endl;
      cerr << "Exp/SuperkDetectorFSI >> Expecting: " <<  np << " x " << np << " matrix" << endl;
      cerr << "Exp/SuperkDetectorFSI >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
      exit(1);
    }

  double det;
  fRHCCovMtrxInv = new TMatrixD(*fRHCCovMtrx);
  (*fRHCCovMtrxInv).SetTol(1E-23);
  (*fRHCCovMtrxInv).Invert(&det);
  //  cout << "Determinant = " << det << endl;
  const TMatrixD &   covInv  = *fRHCCovMtrxInv;
  TMatrixD prod = cov*covInv;


#ifdef __show_init_mesg_debug__
  cout << "Exp/SuperkDetectorFSI >> " << "Error {RHC EnergyScale;1Rmu}    : " << fRHCEnergyScale1RmuRelErr  << endl;
  cout << "Exp/SuperkDetectorFSI >> " << "Error {RHC EnergyScale;1Re}    : " << fRHCEnergyScale1ReRelErr  << endl;
  cout << "Exp/SuperkDetectorFSI >> Input RHC SuperkDetector+FSI covariance matrix: " << endl; 
  //cov.Print();
  util::PrintCompact(cov,7,4);
  cout << "Exp/SuperkDetectorFSI >> Inverse RHC SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintCompact(*fRHCCovMtrxInv,7,4);
  cout << "Exp/SuperkDetectorFSI >> check of RHC product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,4);
#endif

  // Create/init current SuperkDetector+FSI params
  fRHCCurrentValues = new TVectorD(np);
  fRHCNominalValues = new TVectorD(np);
  for(int ip = 0; ip < np; ip++) 
    {
      (*fRHCCurrentValues)[ip] = 1.;
      (*fRHCNominalValues)[ip] = 1.;
    }

#ifdef __show_init_mesg_info__
  cout << "Exp/SuperkDetectorFSI >> Done reading RHC SuperkDetector+FSI inputs!" << endl;
#endif

  return true;  
}
//________________________________________________________________________________________                                                                                        
void SuperkDetectorFSI::TweakEnergyScale(double twkdial)
{
  fEnergyScaleTwkDial = twkdial;
}
//________________________________________________________________________________________
void SuperkDetectorFSI::TweakRHCEnergyScale(double twkdial)
{
  fRHCEnergyScaleTwkDial = twkdial;
}
//________________________________________________________________________________________
double SuperkDetectorFSI::SystWght(
  conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double /*Etrue*/) const
{

  // find SuperkDetector+FSI param which corresponds to the current Ereco param
  int iparam = -1;
  int nSKDetFSIParams;
  conv::PhysParamId_t p;

  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    nSKDetFSIParams = conv::kNumOfRHCSKDetFSIParams;
  else
    nSKDetFSIParams = conv::kNumOfSKDetFSIParams;

  for(int ip = 0; ip < nSKDetFSIParams; ip++) {
    if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
      p = conv::kListOfRHCSKDetFSIParams[ip];
    else
      p = conv::kListOfSKDetFSIParams[ip];
     // Is the current systematic parameter applied to current sample/mode end energy bin?
     bool apply = this->ApplySyst(p,s,m,Ereco);
     if(!apply) continue;
     if(iparam != -1) {
       if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
	 {
	 cerr << "Exp/SuperkDetectorFSI/SystWght >> Found multiple RHC SuperkDetector+FSI bins for sample, mode, Ereco combination!" << endl
		<< util::AsString(s) << "\t" << util::AsString(m) << "\t" << Ereco << endl
#if FSI_BINS != 3
	      << util::AsString((conv::PhysParamId_t)(iparam+conv::kRHCNuisSKDetFSI000)) << "\t" <<  util::AsString(p) << endl;
#else
	      << util::AsString((conv::PhysParamId_t)(iparam+conv::kRHCNuisSKDetFSI006)) << "\t" <<  util::AsString(p) << endl;
#endif
	 }
       else
	 {
         cerr << "Exp/SuperkDetectorFSI/SystWght >> Found multiple SuperkDetector+FSI bins for sample, mode, Ereco combination!" << endl
	    << util::AsString(s) << "\t" << util::AsString(m) << "\t" << Ereco << endl
#if FSI_BINS != 3
	    << util::AsString((conv::PhysParamId_t)(iparam+conv::kNuisSKDetFSI000)) << "\t" <<  util::AsString(p) << endl;
#else
	    << util::AsString((conv::PhysParamId_t)(iparam+conv::kNuisSKDetFSI006)) << "\t" <<  util::AsString(p) << endl;
#endif
	 }
       exit(1);
     }
     if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
       iparam = arr_index->RHCSKDetFSISystIdx(p);
     else
       iparam = arr_index->SKDetFSISystIdx(p);
  }

  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    {
    int np = conv::kNumOfSKDetFSIParams;
    if(iparam < 0 || iparam >= np) {
      cerr << "Exp/SuperkDetectorFSI >> Weird SuperkDetector+FSI Ereco index: " << iparam << " for Ereco value = " << Ereco 
	 << " sample " << util::AsString(s) << " mode " << util::AsString(m) << endl;
      exit(1);
      }
    } 
  else
    {
    int np = conv::kNumOfRHCSKDetFSIParams;
    if(iparam < 0 || iparam >= np) {
      cerr << "Exp/SuperkDetectorFSI >> Weird RHC SuperkDetector+FSI Ereco index: " << iparam << " for Ereco value = " << Ereco
	   << " sample " << util::AsString(s) << " mode " << util::AsString(m) << endl;
      exit(1);
      }
    }

  // get current value of iparam^th FSI param
  double skdetfsi_param_value; 
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    skdetfsi_param_value = (*fRHCCurrentValues)[iparam];
  else
    skdetfsi_param_value = (*fCurrentValues)[iparam];

  // param tweaked
  bool tweaked = fabs(skdetfsi_param_value - 1.) > conv::kASmallNum;
  if(!tweaked) {
    return 1.;
  }

  double skdetfsi_wght = skdetfsi_param_value;
  return skdetfsi_wght;
}
//________________________________________________________________________________________
double SuperkDetectorFSI::CurrentValue(conv::PhysParamId_t p) const
{
  if(p == kNuisSuperkEnergyScale)
    return fEnergyScaleTwkDial;
  const TVectorD & current_values = *fCurrentValues;
  int idx = arr_index->SKDetFSISystIdx(p);
  return current_values[idx];
}
//________________________________________________________________________________________
double SuperkDetectorFSI::RHCCurrentValue(conv::PhysParamId_t p) const
{
  if(p == kRHCNuisSuperkEnergyScale)
    return fRHCEnergyScaleTwkDial;
  const TVectorD & current_values = *fRHCCurrentValues;
  int idx = arr_index->RHCSKDetFSISystIdx(p);
  return current_values[idx];
}
//________________________________________________________________________________________
double SuperkDetectorFSI::OneSigmaError(PhysParamId_t p) const
{
  if(p == kNuisSuperkEnergyScale)
    return 1;
  const TMatrixD & cov = *fCovMtrx;
  int idx = arr_index->SKDetFSISystIdx(p);
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double SuperkDetectorFSI::RHCOneSigmaError(PhysParamId_t p) const
{
  if(p == kRHCNuisSuperkEnergyScale)
    return 1;
  const TMatrixD & cov = *fRHCCovMtrx;
  int idx = arr_index->RHCSKDetFSISystIdx(p);
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
void SuperkDetectorFSI::SyncRHCToFHCSystematics()
{
  fRHCEnergyScaleTwkDial = fEnergyScaleTwkDial;
  this->TweakRHC(true, fCurrentValues);
}
//________________________________________________________________________________________
void SuperkDetectorFSI::Tweak(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_skdetfsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_skdetfsi_cholesky_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/Tweak >> Tweaking SuperkDetector+FSI parameters" << endl;
#endif

  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuMuDisappOnly) {
    this->TweakOnly1Rmu(force_inputs,input_skdetfsi_params);
    return;
  }

  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuEAppOnly) {
    this->TweakOnly1Re(force_inputs,input_skdetfsi_params);
    return;
  }
  
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
        cout << "Exp/SuperkDetectorFSI/Tweak >> Reset SuperkDetector+FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfSKDetFSIParams;
        for(int ip = 0; ip < np; ip++) {
           current_values[ip] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_skdetfsi_tweak_mesg_debug__
        cout << "Exp/SuperkDetectorFSI/Tweak >> Forcing input SuperkDetector+FSI parameters. " << endl;
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
    cout << "Exp/SuperkDetectorFSI/Tweak >> "
         << "Running Cholesky's decomposition method on the SuperkDetector+FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_skdetfsi_cholesky_mesg_debug__
    cout << "Exp/SuperkDetectorFSI/Tweak >> SuperkDetector+FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_skdetfsi_cholesky_mesg_debug__
    cout << "Exp/SuperkDetectorFSI/Tweak >> Cholesky decomposition of SuperkDetector+FSI covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
  {   
#ifdef __show_skdetfsi_tweak_mesg_debug__
     cout << "Exp/SuperkDetectorFSI/Tweak >> Generating random SuperkDetector+FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD & nominal_values = *fNominalValues;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
     current_values = vec;
     return;
  }
}
//________________________________________________________________________________________ 
void SuperkDetectorFSI::TweakRHC(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_skdetfsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_skdetfsi_cholesky_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/Tweak >> Tweaking RHC SuperkDetector+FSI parameters" << endl;
#endif

  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaRHCNuMuDisappOnly) {
    this->TweakOnlyRHC1Rmu(force_inputs,input_skdetfsi_params);
    return;
  }

  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaRHCNuEAppOnly) {
    this->TweakOnlyRHC1Re(force_inputs,input_skdetfsi_params);
    return;
  }
  
  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
    {
      // Input vector is NULL.
      // This resets the FSI params to their best-fit values.
      if(input_skdetfsi_params==0) {
#ifdef __show_skdetfsi_tweak_mesg_debug__
        cout << "Exp/SuperkDetectorFSI/Tweak >> Reset RHC SuperkDetector+FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfRHCSKDetFSIParams;
        for(int ip = 0; ip < np; ip++) {
	  current_values[ip] = 1.;
        }
        return;
      }
      else
	{
#ifdef __show_skdetfsi_tweak_mesg_debug__
	  cout << "Exp/SuperkDetectorFSI/Tweak >> Forcing RHC input SuperkDetector+FSI parameters. " << endl;
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

  if(fRHCCovMtrxChoL == 0)
    {
#ifdef __show_skdetfsi_cholesky_mesg_info__
      cout << "Exp/SuperkDetectorFSI/Tweak >> "
	   << "Running Cholesky's decomposition method on the RHC SuperkDetector+FSI covariance matrix." << endl;
#endif
    
      TMatrixD cov(this->RHCCovMatrix());
#ifdef __show_skdetfsi_cholesky_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/Tweak >> RHC SuperkDetector+FSI covariance matrix: " << endl;
      util::PrintCompact(cov,7,4);
#endif

      TMatrixD L = util::CholeskyDecomposition(cov);
      fRHCCovMtrxChoL = new TMatrixD(L);
#ifdef __show_skdetfsi_cholesky_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/Tweak >> Cholesky decomposition of RHC SuperkDetector+FSI covariance matrix: " << endl;
      util::PrintCompact(*fRHCCovMtrxChoL,7,4);
#endif
    }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
    {   
#ifdef __show_skdetfsi_tweak_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/Tweak >> Generating random RHC SuperkDetector+FSI params (correlated)! " << endl;
#endif
      TMatrixD & L = *fRHCCovMtrxChoL; 
      TVectorD & nominal_values = *fRHCNominalValues;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
      current_values = vec;
      return;
    }
}
//________________________________________________________________________________________
void SuperkDetectorFSI::TweakOnly1Rmu(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_skdetfsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_skdetfsi_cholesky_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Tweaking SuperkDetector+FSI parameters" << endl;
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
        cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Reset SuperkDetector+FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfSKDetFSIParams1Rmu;
        for(int ip = 0; ip < np; ip++) {
           current_values[ip] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_skdetfsi_tweak_mesg_debug__
        cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Forcing input SuperkDetector+FSI parameters. " << endl;
#endif
        TVectorD & vec = *input_skdetfsi_params;
	const int np = conv::kNumOfSKDetFSIParams1Rmu;
	for(int ip = 0; ip < np; ip++) {
	  current_values[ip] = vec[ip];
	}      
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
    cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> "
         << "Running Cholesky's decomposition method on the SuperkDetector+FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_skdetfsi_cholesky_mesg_debug__
    cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> SuperkDetector+FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_skdetfsi_cholesky_mesg_debug__
    cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Cholesky decomposition of SuperkDetector+FSI covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
  {   
#ifdef __show_skdetfsi_tweak_mesg_debug__
     cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Generating random SuperkDetector+FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD & nominal_values = *fNominalValues;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
     const int np = conv::kNumOfSKDetFSIParams1Rmu;
     for(int ip = 0; ip < np; ip++) {
       current_values[ip] = vec[ip];
     }
     return;
  }
}
//________________________________________________________________________________________
void SuperkDetectorFSI::TweakOnlyRHC1Rmu(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_skdetfsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_skdetfsi_cholesky_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Tweaking RHC SuperkDetector+FSI parameters" << endl;
#endif

  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
    {
      // Input vector is NULL.
      // This resets the FSI params to their best-fit values.
      if(input_skdetfsi_params==0) {
#ifdef __show_skdetfsi_tweak_mesg_debug__
        cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Reset RHC SuperkDetector+FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfRHCSKDetFSIParams1Rmu;
        for(int ip = 0; ip < np; ip++) {
	  current_values[ip] = 1.;
        }
        return;
      }
      else
	{
#ifdef __show_skdetfsi_tweak_mesg_debug__
	  cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Forcing RHC input SuperkDetector+FSI parameters. " << endl;
#endif
	  TVectorD & vec = *input_skdetfsi_params;
	  const int np = conv::kNumOfRHCSKDetFSIParams1Rmu;
	  for(int ip = 0; ip < np; ip++) {
	    current_values[ip] = vec[ip];
	  }      
	  return;
	}
    }//force_inputs?

  //
  // If this is the first time this function is called, then run the Cholesky
  // decomposition of the flux covariance matrix
  //

  if(fRHCCovMtrxChoL == 0)
    {
#ifdef __show_skdetfsi_cholesky_mesg_info__
      cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> "
	   << "Running Cholesky's decomposition method on the RHC SuperkDetector+FSI covariance matrix." << endl;
#endif
    
      TMatrixD cov(this->RHCCovMatrix());
#ifdef __show_skdetfsi_cholesky_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> RHC SuperkDetector+FSI covariance matrix: " << endl;
      util::PrintCompact(cov,7,4);
#endif

      TMatrixD L = util::CholeskyDecomposition(cov);
      fRHCCovMtrxChoL = new TMatrixD(L);
#ifdef __show_skdetfsi_cholesky_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Cholesky decomposition of RHC SuperkDetector+FSI covariance matrix: " << endl;
      util::PrintCompact(*fRHCCovMtrxChoL,7,4);
#endif
    }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
    {   
#ifdef __show_skdetfsi_tweak_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/TweakOnly1Rmu >> Generating random RHC SuperkDetector+FSI params (correlated)! " << endl;
#endif
      TMatrixD & L = *fRHCCovMtrxChoL; 
      TVectorD & nominal_values = *fRHCNominalValues;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
      const int np = conv::kNumOfRHCSKDetFSIParams1Rmu;
      for(int ip = 0; ip < np; ip++) {
	current_values[ip] = vec[ip];
      }
      return;
    }
}
//________________________________________________________________________________________
void SuperkDetectorFSI::TweakOnly1Re(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_skdetfsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_skdetfsi_cholesky_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Tweaking SuperkDetector+FSI parameters" << endl;
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
        cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Reset SuperkDetector+FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfSKDetFSIParams1Re;
        for(int ip = 0; ip < np; ip++) {
	  current_values[ip+conv::kNumOfSKDetFSIParams1Rmu] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_skdetfsi_tweak_mesg_debug__
        cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Forcing input SuperkDetector+FSI parameters. " << endl;
#endif
        TVectorD & vec = *input_skdetfsi_params;
	const int np = conv::kNumOfSKDetFSIParams1Re;
	for(int ip = 0; ip < np; ip++) {
	  current_values[ip+conv::kNumOfSKDetFSIParams1Rmu] = vec[ip];
	}
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
    cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> "
         << "Running Cholesky's decomposition method on the SuperkDetector+FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_skdetfsi_cholesky_mesg_debug__
    cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> SuperkDetector+FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_skdetfsi_cholesky_mesg_debug__
    cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Cholesky decomposition of SuperkDetector+FSI covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
  {   
#ifdef __show_skdetfsi_tweak_mesg_debug__
     cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Generating random SuperkDetector+FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD & nominal_values = *fNominalValues;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
     const int np = conv::kNumOfSKDetFSIParams1Re;
     for(int ip = conv::kNumOfSKDetFSIParams1Rmu; ip < np+conv::kNumOfSKDetFSIParams1Rmu; ip++) {
       current_values[ip] = vec[ip];
     }
     return;
  }
}
//________________________________________________________________________________________
void SuperkDetectorFSI::TweakOnlyRHC1Re(
 bool       force_inputs,    /* force the input params or generate random ones internally  */
 TVectorD * input_skdetfsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_skdetfsi_cholesky_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Tweaking RHC SuperkDetector+FSI parameters" << endl;
#endif

  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
    {
      // Input vector is NULL.
      // This resets the FSI params to their best-fit values.
      if(input_skdetfsi_params==0) {
#ifdef __show_skdetfsi_tweak_mesg_debug__
        cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Reset RHC SuperkDetector+FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfRHCSKDetFSIParams1Re;
        for(int ip = 0; ip < np; ip++) {
          current_values[ip+conv::kNumOfRHCSKDetFSIParams1Rmu] = 1.;
        }
        return;
      }
      else
	{
#ifdef __show_skdetfsi_tweak_mesg_debug__
	  cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Forcing RHC input SuperkDetector+FSI parameters. " << endl;
#endif
	  TVectorD & vec = *input_skdetfsi_params;
	  const int np = conv::kNumOfRHCSKDetFSIParams1Re;
	  for(int ip = 0; ip < np; ip++) {
	    current_values[ip+conv::kNumOfRHCSKDetFSIParams1Rmu] = vec[ip];
	  }
	  return;
	}
    }//force_inputs?

  //
  // If this is the first time this function is called, then run the Cholesky
  // decomposition of the flux covariance matrix
  //

  if(fRHCCovMtrxChoL == 0)
    {
#ifdef __show_skdetfsi_cholesky_mesg_info__
      cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> "
	   << "Running Cholesky's decomposition method on the RHC SuperkDetector+FSI covariance matrix." << endl;
#endif
    
      TMatrixD cov(this->RHCCovMatrix());
#ifdef __show_skdetfsi_cholesky_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> RHC SuperkDetector+FSI covariance matrix: " << endl;
      util::PrintCompact(cov,7,4);
#endif

      TMatrixD L = util::CholeskyDecomposition(cov);
      fRHCCovMtrxChoL = new TMatrixD(L);
#ifdef __show_skdetfsi_cholesky_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Cholesky decomposition of RHC SuperkDetector+FSI covariance matrix: " << endl;
      util::PrintCompact(*fRHCCovMtrxChoL,7,4);
#endif
    }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
    {   
#ifdef __show_skdetfsi_tweak_mesg_debug__
      cout << "Exp/SuperkDetectorFSI/TweakOnly1Re >> Generating random RHC SuperkDetector+FSI params (correlated)! " << endl;
#endif
      TMatrixD & L = *fRHCCovMtrxChoL; 
      TVectorD & nominal_values = *fRHCNominalValues;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
      const int np = conv::kNumOfRHCSKDetFSIParams1Re;
      for(int ip = conv::kNumOfRHCSKDetFSIParams1Rmu; ip < np+conv::kNumOfRHCSKDetFSIParams1Rmu; ip++) {
	current_values[ip] = vec[ip];
      }
      return;
    }
}
//________________________________________________________________________________________
double SuperkDetectorFSI::ChisqPenalty(void) const
//just a straight copy of BANFF::ChisqPenalty
{
  //Vector with deviations of current values from  best fit values
  const TVectorD & bestfit_values   = *fNominalValues;
  const TVectorD & current_values   = *fCurrentValues;
  const int np = conv::kNumOfSKDetFSIParams;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> SuperkDetector+FSI parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"SKDetFSI fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in SuperkDetector+FSI parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"SKDetFSI fit params");
#endif

  //covariance matrix
  const TMatrixD cov    (this->CovMatrix());
  const TMatrixD    covInv (this->CovMatrixInv());

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintCompact(cov,7,4);
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> Inverse SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintCompact(covInv,7,4);
#endif

  //Check the inverse:
  TMatrixD prod = cov*covInv;
#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> check of product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,4);
#endif

  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

  //now add the contribution from energy scale
  if (fEnergyScaleTwkDial   ) chisq_penalty += fEnergyScaleTwkDial*fEnergyScaleTwkDial;

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;
#endif

  if(chisq_penalty < 0) 
    cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;

  if(chisq_penalty) return chisq_penalty;
  
  return 0;
}
//________________________________________________________________________________________
double SuperkDetectorFSI::RHCChisqPenalty(void) const
//just a straight copy of BANFF::ChisqPenalty
{
  //Vector with deviations of current values from  best fit values
  const TVectorD & bestfit_values   = *fRHCNominalValues;
  const TVectorD & current_values   = *fRHCCurrentValues;
  const int np = conv::kNumOfRHCSKDetFSIParams;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> RHC SuperkDetector+FSI parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"RHC SKDetFSI fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in RHC SuperkDetector+FSI parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"RHC SKDetFSI fit params");
#endif

  //covariance matrix
  const TMatrixD cov    (this->RHCCovMatrix());
  const TMatrixD    covInv (this->RHCCovMatrixInv());

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> RHC SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintCompact(cov,7,4);
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> RHC inverse SuperkDetector+FSI covariance matrix: " << endl;
  util::PrintCompact(covInv,7,4);
#endif

  //Check the inverse:
  TMatrixD prod = cov*covInv;
#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> check of RHC product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,4);
#endif

  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

  //now add the contribution from energy scale
  if (fRHCEnergyScaleTwkDial   ) chisq_penalty += fRHCEnergyScaleTwkDial*fRHCEnergyScaleTwkDial;

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> RHC chisqPenalty = " << chisq_penalty << endl;
#endif

  if(chisq_penalty < 0) 
    cout << "Exp/SuperkDetectorFSI/ChisqPenalty >> RHC chisqPenalty = " << chisq_penalty << endl;

  if(chisq_penalty) return chisq_penalty;
  
  return 0;
}
//________________________________________________________________________________________
void SuperkDetectorFSI::PrintSystematicVariations(ostream & stream) const
{
  const TVectorD & curr_skdetfsi_params = *fCurrentValues;
  stream << "Exp/SuperkDetectorFSI >> Energy scale twk        = "  << fEnergyScaleTwkDial      << endl;
  stream << "Exp/SuperkDetectorFSI >> SuperkDetector+FSI parameters: " << endl;
  util::PrintCompact(curr_skdetfsi_params,15,6,"SKDetFSI params",stream);

}
//________________________________________________________________________________________
void SuperkDetectorFSI::PrintRHCSystematicVariations(ostream & stream) const
{
  const TVectorD & curr_skdetfsi_params = *fRHCCurrentValues;
  stream << "Exp/SuperkDetectorFSI >> RHC energy scale twk        = "  << fRHCEnergyScaleTwkDial      << endl;
  stream << "Exp/SuperkDetectorFSI >> RHC SuperkDetector+FSI parameters: " << endl;
  util::PrintCompact(curr_skdetfsi_params,15,6,"RHC SKDetFSI params",stream);

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
  int idx;

  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    idx = arr_index->RHCSKDetFSISystIdx(p);
  else
    idx = arr_index->SKDetFSISystIdx(p);
  if (idx == -1) return false;
  int is = arr_index->SuperkSampleIdx(s);
  if (is == -1) return false;
  int im = arr_index->SuperkModeIdx(s,m);
  if (im == -1) return false;

  bool apply = conv::kSKDetFSIApplySKDetFSIParamToSKMode[is][im][idx];
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
  TString paramName = AsString(p);
  int idx;

  if(paramName.Contains("RHC"))
    idx = arr_index->RHCSKDetFSISystIdx(p);
  else
    idx = arr_index->SKDetFSISystIdx(p);
  if (idx == -1) return false;

  bool apply;
  if(paramName.Contains("RHC"))
    apply = Ereco >= conv::kRHCSKDetFSIParamsErecoCut[idx][0] &&
      Ereco <  conv::kRHCSKDetFSIParamsErecoCut[idx][1];
  else
    apply = Ereco >= conv::kSKDetFSIParamsErecoCut[idx][0] &&
               Ereco <  conv::kSKDetFSIParamsErecoCut[idx][1];

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
  // energy scale relative error
  //Todo: add separate energy scale error for RHC running if necessary
  if(s == kSKsFCFV1RingMuLike)
    return fEnergyScale1RmuRelErr;
  else if (s == kSKsFCFV1RingELike) 
    return fEnergyScale1RmuRelErr;
  else if(s == kSKsRHCFCFV1RingMuLike)
    return fEnergyScale1RmuRelErr;
  else if (s == kSKsRHCFCFV1RingELike)
    return fEnergyScale1RmuRelErr;

      //lorena - so far never seen different energy scale for signal/background                                                                                                  
      //    if( util::IsSignalNuE(m)) {                                                                                                                                          
      /////////      return fEnergyScale1ReSigRelErr;                                                                                                                             
      // }                                                                                                                                                                      
      //else {                                                                                                                                                                   
      /////////      return fEnergyScale1ReBkgRelErr;                                                                                                                            
      //}                                                                                                                                                                         
    
  return 0;
}
//________________________________________________________________________________________                                                                                       
double SuperkDetectorFSI::EnergyScale(SuperkSample_t s, SuperkMode_t m) const
{
  double twkdial;
  // energy scale
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    twkdial = fRHCEnergyScaleTwkDial;
  else
    twkdial = fEnergyScaleTwkDial;
  double err     = this->EnergyScaleRelErr(s,m);
  double scale   = 1. + twkdial * err;
  return scale;
}
//________________________________________________________________________________________                                                                                        
double SuperkDetectorFSI::EnergyScale(SuperkSample_t s) const
{
  double twkdial;
  // energy scale
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    twkdial = fRHCEnergyScaleTwkDial;
  else       
    twkdial = fEnergyScaleTwkDial;
  double err     = this->EnergyScaleRelErr(s,kSKmNuMuCCQE);
  double scale   = 1. + twkdial * err;
  return scale;
}
//________________________________________________________________________________________                                                                                          
