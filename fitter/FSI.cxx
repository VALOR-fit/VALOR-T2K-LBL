//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_init_mesg_debug__
//#define __show_fsi_cholesky_mesg_info__
//#define __show_fsi_cholesky_mesg_debug__
//#define __show_fsi_tweak_mesg_debug__
//#define __show_chisqpenalty_mesg_debug__
//#define __show_fsi_applysyst_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TMath.h>

#include "FSI.h"
#include "Utilities.h"
#include "Opt.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

using namespace t2kana;

ClassImp(FSI)

//________________________________________________________________________________________
FSI::FSI()
{
  this->Init();
}
//________________________________________________________________________________________
FSI::~FSI()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void FSI::Init(void)
{
  fCurrentValues = 0;
  fNominalValues = 0;
  fCovMtrx       = 0;
  fCovMtrxChoL   = 0;
}
//________________________________________________________________________________________
void FSI::CleanUp(void)
{
  if ( fCurrentValues ) delete fCurrentValues;
  if ( fNominalValues ) delete fNominalValues;
  if ( fCovMtrx       ) delete fCovMtrx;
  if ( fCovMtrxChoL   ) delete fCovMtrxChoL;
}
//________________________________________________________________________________________
bool FSI::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/FSI >> Reading FSI inputs from input XML file." << endl;
#endif

  const int np = conv::kNumOfFSIParams;

  // Get FSI covariance matrix
#if FSI_BINS == 1
  fCovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/covariance_matrix");
#elif FSI_BINS == 2
  fCovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/disap/covariance_matrix");
#elif FSI_BINS == 3
  fCovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/app/covariance_matrix");
#endif
  if(!fCovMtrx) {
    cerr << "Exp/FSI >> No FSI covariance matrix was found." << endl;
    exit(1);
  }

  const TMatrixD & cov = *fCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np)
  {
    cerr << "Exp/FSI >> FSI covariance matrix has unexpected size!" << endl;
    cerr << "Exp/FSI >> Expecting: " <<  np << " x " << np << " matrix" << endl;
    cerr << "Exp/FSI >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
    exit(1);
  }

  double det;
  fCovMtrxInv = new TMatrixD(*fCovMtrx);
  (*fCovMtrxInv).Invert(&det);

#ifdef __show_init_mesg_debug__
  cout << "Exp/FSI >> Input FSI covariance matrix: " << endl; 
  util::PrintCompact(cov,7,4);
  cout << "Exp/FSI >> Inverse FSI covariance matrix: " << endl;
  util::PrintCompact(*fCovMtrxInv,7,4);
#endif

  // Create/init current FSI params
  fCurrentValues = new TVectorD(np);
  fNominalValues = new TVectorD(np);
  for(int ip = 0; ip < np; ip++) 
  {
     (*fCurrentValues)[ip] = 1.;
     (*fNominalValues)[ip] = 1.;
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/FSI >> Done reading neutrino cross-section inputs!" << endl;
#endif

  return true;  
}
//________________________________________________________________________________________
double FSI::SystWght(
  conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double /*Etrue*/) const
{
  // find FSI param which corresponds to the current Ereco param
  int iparam = -1;
  for(int ip = 0; ip < conv::kNumOfFSIParams; ip++) {
     conv::PhysParamId_t p = conv::kListOfFSIParams[ip];
     // Is the current systematic parameter applied to current sample/mode end energy bin?
     bool apply = this->ApplySyst(p,s,m,Ereco);
     if(!apply) continue;
     if(iparam != -1) {
       cerr << "Exp/FSI/SystWght >> Found multiple FSI bins for sample, mode, Ereco combination!" << endl
	    << util::AsString(s) << "\t" << util::AsString(m) << "\t" << Ereco << endl
#if FSI_BINS != 3
	    << util::AsString((conv::PhysParamId_t)(iparam+conv::kNuisFSI000)) << "\t" <<  util::AsString(p) << endl;
#else
	    << util::AsString((conv::PhysParamId_t)(iparam+conv::kNuisFSI006)) << "\t" <<  util::AsString(p) << endl;
#endif
       exit(1);
     }
     iparam = util::ArrIdx::Instance()->FSISystIdx(p);
  }
  const int np = conv::kNumOfFSIParams;
  if(iparam < 0 || iparam >= np) {
    cerr << "Exp/FSI >> Weird FSI Ereco index: " << iparam << " for Ereco value = " << Ereco 
	 << " sample " << util::AsString(s) << " mode " << util::AsString(m) << endl;
     exit(1);
  }
 
  // get current value of iparam^th FSI param
  double fsi_param_value = (*fCurrentValues)[iparam];

  // param tweaked
  bool tweaked = TMath::Abs(fsi_param_value - 1.) > conv::kASmallNum;
  if(!tweaked) {
    return 1.;
  }

  double fsi_wght = fsi_param_value;
  return fsi_wght;
}
//________________________________________________________________________________________
double FSI::CurrentValue(conv::PhysParamId_t p) const
{
  const TVectorD & current_values = *fCurrentValues;
  int idx = util::ArrIdx::Instance()->FSISystIdx(p);
  return current_values[idx];
}
//________________________________________________________________________________________
void FSI::Tweak(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_fsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_fsi_cholesky_mesg_debug__
  cout << "Exp/FSI/Tweak >> Tweaking FSI parameters" << endl;
#endif

  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuMuDisappOnly) {
    this->TweakOnly1Rmu(force_inputs,input_fsi_params);
    return;
  }

  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuEAppOnly) {
    this->TweakOnly1Re(force_inputs,input_fsi_params);
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
     if(input_fsi_params==0) {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/Tweak >> Reset FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfFSIParams;
        for(int ip = 0; ip < np; ip++) {
           current_values[ip] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/Tweak >> Forcing input FSI parameters. " << endl;
#endif
        TVectorD & vec = *input_fsi_params;
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
#ifdef __show_fsi_cholesky_mesg_info__
    cout << "Exp/FSI/Tweak >> "
         << "Running Cholesky's decomposition method on the FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/Tweak >> FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/Tweak >> Cholesky decomposition of FSI covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
  {   
#ifdef __show_fsi_tweak_mesg_debug__
     cout << "Exp/FSI/Tweak >> Generating random FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD & nominal_values = *fNominalValues;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
     current_values = vec;
     return;
  }
}
//________________________________________________________________________________________
void FSI::TweakOnly1Rmu(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_fsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_fsi_cholesky_mesg_debug__
  cout << "Exp/FSI/TweakOnly1Rmu >> Tweaking FSI parameters" << endl;
#endif

  TVectorD & current_values = *fCurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
  {
     // Input vector is NULL.
     // This resets the FSI params to their best-fit values.
     if(input_fsi_params==0) {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/TweakOnly1Rmu >> Reset FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfFSIParams1Rmu;
        for(int ip = 0; ip < np; ip++) {
           current_values[ip] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/TweakOnly1Rmu >> Forcing input FSI parameters. " << endl;
#endif
        TVectorD & vec = *input_fsi_params;
	const int np = conv::kNumOfFSIParams1Rmu;
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
#ifdef __show_fsi_cholesky_mesg_info__
    cout << "Exp/FSI/TweakOnly1Rmu >> "
         << "Running Cholesky's decomposition method on the FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/TweakOnly1Rmu >> FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/TweakOnly1Rmu >> Cholesky decomposition of FSI covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
  {   
#ifdef __show_fsi_tweak_mesg_debug__
     cout << "Exp/FSI/TweakOnly1Rmu >> Generating random FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD & nominal_values = *fNominalValues;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
     const int np = conv::kNumOfFSIParams1Rmu;
     for(int ip = 0; ip < np; ip++) {
       current_values[ip] = vec[ip];
     }
     return;
  }
}
//________________________________________________________________________________________
void FSI::TweakOnly1Re(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_fsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_fsi_cholesky_mesg_debug__
  cout << "Exp/FSI/TweakOnly1Re >> Tweaking FSI parameters" << endl;
#endif

  TVectorD & current_values = *fCurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
  {
     // Input vector is NULL.
     // This resets the FSI params to their best-fit values.
     if(input_fsi_params==0) {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/TweakOnly1Re >> Reset FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfFSIParams1Re;
        for(int ip = 0; ip < np; ip++) {
	  current_values[ip+conv::kNumOfFSIParams1Rmu] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/TweakOnly1Re >> Forcing input FSI parameters. " << endl;
#endif
        TVectorD & vec = *input_fsi_params;
	const int np = conv::kNumOfFSIParams1Re;
	for(int ip = 0; ip < np; ip++) {
	  current_values[ip+conv::kNumOfFSIParams1Rmu] = vec[ip];
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
#ifdef __show_fsi_cholesky_mesg_info__
    cout << "Exp/FSI/TweakOnly1Re >> "
         << "Running Cholesky's decomposition method on the FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/TweakOnly1Re >> FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/TweakOnly1Re >> Cholesky decomposition of FSI covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  if(!force_inputs)
  {   
#ifdef __show_fsi_tweak_mesg_debug__
     cout << "Exp/FSI/TweakOnly1Re >> Generating random FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD & nominal_values = *fNominalValues;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
     const int np = conv::kNumOfFSIParams1Re;
     for(int ip = conv::kNumOfFSIParams1Rmu; ip < np+conv::kNumOfFSIParams1Rmu; ip++) {
       current_values[ip] = vec[ip];
     }
     return;
  }
}
//________________________________________________________________________________________
double FSI::ChisqPenalty(void) const
//just a straight copy of BANFF::ChisqPenalty
{
  //Vector with deviations of current values from  best fit values
  const TVectorD & bestfit_values   = *fNominalValues;
  const TVectorD & current_values   = *fCurrentValues;
  const int np = conv::kNumOfFSIParams;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/FSI/ChisqPenalty >> FSI parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"FSI fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in FSI parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"FSI fit params");
#endif

  //covariance matrix
  const TMatrixD cov    (this->CovMatrix());
  const TMatrixD    covInv (this->CovMatrixInv());

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/FSI/ChisqPenalty >> FSI covariance matrix: " << endl;
  util::PrintCompact(cov,7,4);
  cout << "Exp/FSI/ChisqPenalty >> Inverse FSI covariance matrix: " << endl;
  util::PrintCompact(covInv,7,4);
#endif

  //Check the inverse:
  TMatrixD prod = cov*covInv;
#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/FSI/ChisqPenalty >> check of product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,4);
#endif

  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/FSI/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;
#endif

  if(chisq_penalty) return chisq_penalty;
  
  return 0;
}
//________________________________________________________________________________________
void FSI::PrintSystematicVariations(ostream & stream) const
{
  const TVectorD & curr_fsi_params = *fCurrentValues;
  stream << "Exp/FSI >> FSI parameters: " << endl;
  util::PrintCompact(curr_fsi_params,15,6,"FSI params",stream);
}
//________________________________________________________________________________________
bool FSI::ApplySyst(
  conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const
{
  bool apply2mode = this->ApplySystToSKMode(p,s,m);
  bool apply2bin  = this->ApplySystToEnergyBin(p,Ereco);
  bool apply = apply2mode && apply2bin;
  return apply;
}
//________________________________________________________________________________________
bool FSI::ApplySystToSKMode(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const
{
  int idx = util::ArrIdx::Instance()->FSISystIdx(p);
  if (idx == -1) return false;
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if (is == -1) return false;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if (im == -1) return false;

  bool apply = conv::kFSIApplyFSIParamToSKMode[is][im][idx];
#ifdef __show_fsi_applysyst_mesg_debug__
  string applystr;
  applystr = apply ? "Passed" : "Failed";
  cout << applystr << "\t" << util::AsString(p) << "\t" << util::AsString(s) << "\t" << util::AsString(m) << endl;
#endif
  return apply;
}
//________________________________________________________________________________________
bool FSI::ApplySystToEnergyBin(conv::PhysParamId_t p, double Ereco) const
{
  int idx = util::ArrIdx::Instance()->FSISystIdx(p);
  if (idx == -1) return false;

  bool apply = Ereco >= conv::kFSIParamsErecoCut[idx][0] &&
               Ereco <  conv::kFSIParamsErecoCut[idx][1];
#ifdef __show_fsi_applysyst_mesg_debug__
  string applystr;
  applystr = apply ? "Passed" : "Failed";
  cout << applystr << "\t" << util::AsString(p) << "\t" << Ereco << endl;
#endif
  return apply;
}
//________________________________________________________________________________________
