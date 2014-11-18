//________________________________________________________________________________________
/*

 \copyright (c) 2010-2013, VaLOR group.

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__  
//#define __show_init_mesg_debug__   
//#define __show_cholesky_mesg_info__
//#define __show_cholesky_mesg_debug__
//#define __show_tweak_mesg_debug__
//#define __show_update_correlated_mesg_debug__
//#define __show_chisqpenalty_mesg_debug__
//#define __show_chisqpenalty_mesg_debug2__
//#define __BANFF_BESTFIT_IS_ONE__ //debug mode to turn the best fit vector into a vector of 1's (i.e. BANFF tuned == flux tuned)

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip> 
#include <TRandom3.h>
#include <TMath.h>

#include "Enumerations.h"
#include "Utilities.h"
#include "BANFF.h" 
#include "Conventions.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace valor;
using namespace valor::conv;
using namespace valor::util;

ClassImp(BANFF)

//________________________________________________________________________________________
BANFF::BANFF() :
TObject()
{
  this->Init();
}
//________________________________________________________________________________________
BANFF::~BANFF()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void BANFF::Init(void)
{
  fBestFitValues        = 0; 
  fCurrentValues        = 0; 
  fCovMtrx              = 0; 
  fCovMtrxChoL          = 0; 
  fRHCBestFitValues     = 0;
  fRHCCurrentValues     = 0;
  fRHCCovMtrx           = 0;
  fRHCCovMtrxChoL       = 0;
 
  // NEW FLUX FHC+RHC
  fFHCpRHCBestFitValues     = 0;
  fFHCpRHCCurrentValues     = 0;
  fFHCpRHCCovMtrx           = 0;
  fFHCpRHCCovMtrxChoL       = 0;
  //

  fIsTweaked            = false;
  fIsFluxWghtUpdated    = true;
  fIsTweakedRHC         = false;
  fIsFluxWghtUpdatedRHC = true;

  // NEW FLUX FHC+RHC
  fIsTweakedFHCpRHC         = false;
  fIsFluxWghtUpdatedFHCpRHC = true;
  //
  
  for (int ip = 0 ; ip < conv::knppmax ; ip++) {
    PhysParamId_t p = (PhysParamId_t) ip;
    fParam2Idx[ip]=-1;
    for (int ibp = 0 ; ibp < kNumOfParamsInBANFFfit ; ibp++) {
      if (kListOfBANFFParams[ibp] == p) {
        fParam2Idx[ip] = ibp;
      }
    }
    for (int ibp = 0 ; ibp < kNumOfRHCParamsInBANFFfit ; ibp++) {
      if (kListOfRHCBANFFParams[ibp] == p) {
        fParam2Idx[ip] = ibp;
      }
    }
  }
  
  // NEW FLUX FHC+RHC

#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__)
  for (int ip = 0 ; ip < conv::knppmax ; ip++) {
    PhysParamId_t p = (PhysParamId_t) ip;
    fParam2Idx[ip]=-1;
    for (int ibp = 0 ; ibp < kNumOfParamsInBANFFfit ; ibp++) {
      if (kListOfBANFFParams[ibp] == p) {
        fParam2Idx[ip] = ibp;
      }
    }
    for (int ibp = 0 ; ibp < kNumOfRHCParamsInBANFFfit ; ibp++) {
      if (kListOfRHCBANFFParams[ibp] == p) {
        fParam2Idx[ip] = ibp;
      }
    }
  }
  #endif
#endif
  //
  
}
//________________________________________________________________________________________
void BANFF::CleanUp(void)
{
  if ( fBestFitValues             ) delete fBestFitValues;
  if ( fCurrentValues             ) delete fCurrentValues;
  if ( fCovMtrx                   ) delete fCovMtrx;
  if ( fCovMtrxChoL               ) delete fCovMtrxChoL;
  if ( fRHCBestFitValues          ) delete fRHCBestFitValues;
  if ( fRHCCurrentValues          ) delete fRHCCurrentValues;
  if ( fRHCCovMtrx                ) delete fRHCCovMtrx;
  if ( fRHCCovMtrxChoL            ) delete fRHCCovMtrxChoL;

  // NEW FLUX FHC+RHC
  if ( fFHCpRHCBestFitValues          ) delete fFHCpRHCBestFitValues;
  if ( fFHCpRHCCurrentValues          ) delete fFHCpRHCCurrentValues;
  if ( fFHCpRHCCovMtrx                ) delete fFHCpRHCCovMtrx;
  if ( fFHCpRHCCovMtrxChoL            ) delete fFHCpRHCCovMtrxChoL;
  //
}
//________________________________________________________________________________________
bool BANFF::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read BANFF inputs to oscillation analysis
//
#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Reading BANFF inputs from XML configuration file." << endl;
#endif

  const int np = conv::kNumOfParamsInBANFFfit;

  // Get BANFF fit central values

#ifndef __BANFF_BESTFIT_IS_ONE__

  // #if BANFF_BINS == 1
  //   fBestFitValues = util::XML2TVectorD(xml,root_node,"banff/best_fit_values");
  // #elif BANFF_BINS == 2
  //   fBestFitValues = util::XML2TVectorD(xml,root_node,"banff/disap/best_fit_values");
  // #elif BANFF_BINS == 3
  //   fBestFitValues = util::XML2TVectorD(xml,root_node,"banff/app/best_fit_values");
  // #endif
  fBestFitValues = util::XML2TVectorD(xml,root_node,T2KOpt::Instance()->fPathBanffParam.c_str());

  if(!fBestFitValues) {
    cerr << "Exp/BANFF >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
#else
  double init[np];
  for(int ip = 0; ip < np; ip++)
    init[ip] = 1;
  fBestFitValues = new TVectorD(np, init);
#endif

  const TVectorD & bestfit = *fBestFitValues;
  if(bestfit.GetNoElements() != np) 
  {
    cerr << "Exp/BANFF >> BANFF best-fit parameter vector has unexpected number of elements!" << endl;
    cerr << "Exp/BANFF >> Expecting: " <<  np << " elements" << endl;
    cerr << "Exp/BANFF >> Found: " << bestfit.GetNoElements() << " elements" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input BANFF best-fit parameters: " << endl; 
  util::PrintCompact(bestfit,12,6,"BANFF fit params");
#endif

  // init current values to best-fit values
  fCurrentValues = new TVectorD(bestfit);

  // Get BANFF fit covariance matrix

  // #if BANFF_BINS == 1
  //   fCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/covariance_matrix");
  // #elif BANFF_BINS == 2
  //   fCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/disap/covariance_matrix");
  // #elif BANFF_BINS == 3
  //   fCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/app/covariance_matrix");
  // #endif  
    fCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathBanffCovMatr.c_str());
 
  if(!fCovMtrx) {
    cerr << "Exp/BANFF >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  // Invert BANFF fit covariance matrix
  double det;
  fCovMtrxInv = new TMatrixD(*fCovMtrx);
  (*fCovMtrxInv).Invert(&det);

  const TMatrixD & cov = *fCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np) 
  {
    cerr << "Exp/BANFF >> BANFF covariance matrix has unexpected size!" << endl;
    cerr << "Exp/BANFF >> Expecting: " <<  np << " x " << np << " matrix" << endl;
    cerr << "Exp/BANFF >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
  }
#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input BANFF covariance matrix: " << endl; 
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF >> Inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(*fCovMtrxInv,7,5);
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading BANFF inputs!" << endl;
#endif
  
  fBestFitArray = fBestFitValues->GetMatrixArray();
  fCurrentArray = fCurrentValues->GetMatrixArray();
  return true;
}
//________________________________________________________________________________________
bool BANFF::ReadRHC(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  // Read BANFF inputs to oscillation analysis
  //
#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Reading RHC BANFF inputs from XML configuration file." << endl;
#endif

  const int np = conv::kNumOfRHCParamsInBANFFfit;

  // Get BANFF fit central values

#ifndef __BANFF_BESTFIT_IS_ONE__

  // #if BANFF_BINS == 1
  //   fRHCBestFitValues = util::XML2TVectorD(xml,root_node,"banff/best_fit_values");
  // #elif BANFF_BINS == 2
  //   fRHCBestFitValues = util::XML2TVectorD(xml,root_node,"banff/disap/best_fit_values");
  // #elif BANFF_BINS == 3
  //   fRHCBestFitValues = util::XML2TVectorD(xml,root_node,"banff/app/best_fit_values");
  // #endif
  fRHCBestFitValues = util::XML2TVectorD(xml,root_node,T2KOpt::Instance()->fPathRHCBanffParam.c_str());
    
  if(!fRHCBestFitValues) {
    cerr << "Exp/BANFF >> No RHC BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
#else
  double init[np];
  for(int ip = 0; ip < np; ip++)
    init[ip] = 1;
  fRHCBestFitValues = new TVectorD(np, init);
#endif

  const TVectorD & bestfit = *fRHCBestFitValues;
  if(bestfit.GetNoElements() != np) 
    {
      cerr << "Exp/BANFF >> RHC BANFF best-fit parameter vector has unexpected number of elements!" << endl;
      cerr << "Exp/BANFF >> Expecting: " <<  np << " elements" << endl;
      cerr << "Exp/BANFF >> Found: " << bestfit.GetNoElements() << " elements" << endl;
      exit(1);
    }

#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> RHC input BANFF best-fit parameters: " << endl; 
  util::PrintCompact(bestfit,12,6,"RHC BANFF fit params");
#endif

  // init current values to best-fit values
  fRHCCurrentValues = new TVectorD(bestfit);

  // Get BANFF fit covariance matrix

  //#if BANFF_BINS == 1
  // #ifdef NEW_RHC_BANFF
  //     fRHCCovMtrx = util::XML2TMatrixD(xml,root_node,"banff_rhc/covariance_matrix");
  // #else
  //     fRHCCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/covariance_matrix");
  // #endif
  // #elif BANFF_BINS == 2
  //   fRHCCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/disap/covariance_matrix");
  // #elif BANFF_BINS == 3
  //   fRHCCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/app/covariance_matrix");
  // #endif
  fRHCCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathRHCBanffCovMatr.c_str());

  if(!fRHCCovMtrx) {
    cerr << "Exp/BANFF >> No RHC BANFF covariance matrix was found." << endl;
    exit(1);
  }
  // Invert BANFF fit covariance matrix
  double det;
  fRHCCovMtrxInv = new TMatrixD(*fRHCCovMtrx);
  (*fRHCCovMtrxInv).Invert(&det);

  const TMatrixD & cov = *fRHCCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np) 
    {
      cerr << "Exp/BANFF >> RHC BANFF covariance matrix has unexpected size!" << endl;
      cerr << "Exp/BANFF >> Expecting: " <<  np << " x " << np << " matrix" << endl;
      cerr << "Exp/BANFF >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
    }
#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> RHC input BANFF covariance matrix: " << endl; 
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF >> RHC inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(*fRHCCovMtrxInv,7,5);
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading RHC BANFF inputs!" << endl;
#endif
  
  fRHCBestFitArray = fRHCBestFitValues->GetMatrixArray();
  fRHCCurrentArray = fRHCCurrentValues->GetMatrixArray();
  return true;
}










// NEW FLUX FHC+RHC

bool BANFF::ReadFHCpRHC(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read BANFF inputs to oscillation analysis
//
#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Reading BANFF inputs from XML configuration file." << endl;
#endif
  const int np = conv::kNumOfFHCpRHCParamsInBANFFfit;

  // Get BANFF fit central values

#ifndef __BANFF_BESTFIT_IS_ONE__
  fFHCpRHCBestFitValues = util::XML2TVectorD(xml,root_node,T2KOpt::Instance()->fPathBanffParam.c_str());
  
  cout << T2KOpt::Instance()->fPathBanffParam.c_str() << endl;
  
  if(!fFHCpRHCBestFitValues) {
    cerr << "Exp/BANFF >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
#else
  double init[np];
  for(int ip = 0; ip < np; ip++)
    init[ip] = 1;
  fFHCpRHCBestFitValues = new TVectorD(np, init);
#endif

  const TVectorD & bestfit = *fFHCpRHCBestFitValues;
  if(bestfit.GetNoElements() != np) 
  {
    cerr << "Exp/BANFF >> FHC+RHC BANFF best-fit parameter vector has unexpected number of elements!" << endl;
    cerr << "Exp/BANFF >> Expecting: " <<  np << " elements" << endl;
    cerr << "Exp/BANFF >> Found: " << bestfit.GetNoElements() << " elements" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input FHC+RHC BANFF best-fit parameters: " << endl; 
  util::PrintCompact(bestfit,12,6,"FHC+RHC BANFF fit params");
#endif

  // init current values to best-fit values
  fFHCpRHCCurrentValues = new TVectorD(bestfit);

  // Get BANFF fit covariance matrix
  fFHCpRHCCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathBanffCovMatr.c_str());
 
  if(!fFHCpRHCCovMtrx) {
    cerr << "Exp/BANFF >> No FHC+RHC BANFF covariance matrix was found." << endl;
    exit(1);
  }
  // Invert BANFF fit covariance matrix
  double det;
  fFHCpRHCCovMtrxInv = new TMatrixD(*fFHCpRHCCovMtrx);
  (*fFHCpRHCCovMtrxInv).Invert(&det);

  const TMatrixD & cov = *fFHCpRHCCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np) 
  {
    cerr << "Exp/BANFF >> FHC+RHC BANFF covariance matrix has unexpected size!" << endl;
    cerr << "Exp/BANFF >> Expecting: " <<  np << " x " << np << " matrix" << endl;
    cerr << "Exp/BANFF >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
  }
#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input FHC+RHC BANFF covariance matrix: " << endl; 
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF >> Inverse FHC+RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(*fFHCpRHCCovMtrxInv,7,5);
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading BANFF inputs!" << endl;
#endif
  
  fFHCpRHCBestFitArray = fFHCpRHCBestFitValues->GetMatrixArray();
  fFHCpRHCCurrentArray = fFHCpRHCCurrentValues->GetMatrixArray();
  
  // Fill FHC and RHC parameter arrays
  // for(int ipar=0;ipar<kNumOfFluxParamsInBANFFfit;ipar++)
  //   fBestFitValues(ipar) = fFHCpRHCBestFitValues(ipar);
  // for(int ipar=0;ipar<kNumOfRHCFluxParamsInBANFFfit;ipar++)
  //   fRHCBestFitValues(ipar) = fFHCpRHCBestFitValues(kNumOfFluxParamsInBANFFfit+ipar);
  // for(int ipar=0;ipar<kNumOfFHCpRHCXSecParamsInBANFFfit;ipar++){
  //   fBestFitValues(kNumOfFluxParamsInBANFFfit+ipar)    = fFHCpRHCBestFitValues(kNumOfFluxParamsInBANFFfit+kNumOfFluxParamsInBANFFfit+ipar);
  //   fRHCBestFitValues(kNumOfFluxParamsInBANFFfit+ipar) = fFHCpRHCBestFitValues(kNumOfFluxParamsInBANFFfit+kNumOfFluxParamsInBANFFfit+ipar);
  // }
  
  return true;
}

//




//________________________________________________________________________________________
double BANFF::BestFitValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx(p);
  return fBestFitArray[idx];
}
//________________________________________________________________________________________
double BANFF::BestFitValue(int idx) const
{
  return fBestFitArray[idx];
}
//________________________________________________________________________________________
double BANFF::RHCBestFitValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx(p);
  return fRHCBestFitArray[idx];
}
//________________________________________________________________________________________
double BANFF::RHCBestFitValue(int idx) const
{
  return fRHCBestFitArray[idx];
}
//________________________________________________________________________________________

// NEW FLUX FHC+RHC
double BANFF::FHCpRHCBestFitValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx(p);
  return fFHCpRHCBestFitArray[idx];
}
//________________________________________________________________________________________
double BANFF::FHCpRHCBestFitValue(int idx) const
{
  return fFHCpRHCBestFitArray[idx];
}
//


//________________________________________________________________________________________
double BANFF::CurrentValue(int idx) const
{
  return fCurrentArray[idx];
}
//________________________________________________________________________________________
double BANFF::CurrentValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx(p);
  return fCurrentArray[idx];
}
//________________________________________________________________________________________
double BANFF::RHCCurrentValue(int idx) const
{
  return fRHCCurrentArray[idx];
}
//________________________________________________________________________________________
double BANFF::RHCCurrentValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx(p);
  return fRHCCurrentArray[idx];
}
//________________________________________________________________________________________


// NEW FLUX FHC+RHC
double BANFF::FHCpRHCCurrentValue(int idx) const
{
  return fFHCpRHCCurrentArray[idx];
}
//________________________________________________________________________________________
double BANFF::FHCpRHCCurrentValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx(p);
  return fFHCpRHCCurrentArray[idx];
}
//________________________________________________________________________________________
//


double BANFF::OneSigmaError(PhysParamId_t p) const
{
  const TMatrixD & cov = *fCovMtrx;
  int idx = this->Param2Idx(p);
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::OneSigmaError(int idx) const
{
  const TMatrixD & cov = *fCovMtrx;
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::RHCOneSigmaError(PhysParamId_t p) const
{
  const TMatrixD & cov = *fRHCCovMtrx;
  int idx = this->Param2Idx(p);
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::RHCOneSigmaError(int idx) const
{
  const TMatrixD & cov = *fRHCCovMtrx;
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________


// NEW FLUX FHC+RHC
double BANFF::FHCpRHCOneSigmaError(PhysParamId_t p) const
{
  const TMatrixD & cov = *fFHCpRHCCovMtrx;
  int idx = this->Param2Idx(p);
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::FHCpRHCOneSigmaError(int idx) const
{
  const TMatrixD & cov = *fFHCpRHCCovMtrx;
  return TMath::Sqrt(cov[idx][idx]);
}
//


//________________________________________________________________________________________
void BANFF::SyncRHCToFHCSystematics()
{
  this->TweakRHC(true, fCurrentValues);
}
//________________________________________________________________________________________
void BANFF::Tweak(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //
  if(fCovMtrxChoL == 0) 
    this->CholeskyDecomposeMatrix();

  //
  // Get best-fit and current BANFF params
  //

  
  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & current_values = *fCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset BANFF params to best-fit values." << endl;
#endif
        fIsTweaked = false;
        current_values = bestfit_values;
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Forcing input BANFF parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
        assert(vec.GetNoElements() == conv::kNumOfParamsInBANFFfit);
        current_values = vec;
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/Tweak >> Generating random twk dials (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
     current_values = vec;
     return;
  }

}
//________________________________________________________________________________________
void BANFF::TweakRHC(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking RHC BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdatedRHC = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fRHCBestFitValues;
  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset RHC BANFF params to best-fit values." << endl;
#endif
        fIsTweakedRHC = false;
        current_values = bestfit_values;
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/Tweak >> Forcing RHC input BANFF parameters. " << endl;
#endif
	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfRHCParamsInBANFFfit);
	  current_values = vec;
	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/Tweak >> Generating RHC random twk dials (correlated)! " << endl;
#endif
      TMatrixD & L = *fRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      current_values = vec;
      return;
    }

}
//________________________________________________________________________________________




// NEW FLUX FHC+RHC
void BANFF::TweakFHCpRHC(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking FHC+RHC BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweakedFHCpRHC         = true;
  fIsFluxWghtUpdatedFHCpRHC = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fFHCpRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeFHCpRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset FHC+RHC BANFF params to best-fit values." << endl;
#endif
        fIsTweakedFHCpRHC = false;
        current_values = bestfit_values;
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/Tweak >> Forcing FHC+RHC input BANFF parameters. " << endl;
#endif
	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfFHCpRHCParamsInBANFFfit);
	  current_values = vec;
	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/Tweak >> Generating RHC random twk dials (correlated)! " << endl;
#endif
      TMatrixD & L = *fFHCpRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      current_values = vec;
      return;
    }
}

//________________________________________________________________________________________
void BANFF::TweakFluxOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakFluxOnly >> Tweaking BANFF flux parameters" << endl;
#endif

  // Switch flags
  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fCovMtrxChoL == 0) 
    this->CholeskyDecomposeMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & current_values = *fCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakFluxOnly >> Reset BANFF flux params to best-fit values." << endl;
#endif
        fIsTweaked = false;
        current_values = bestfit_values;
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakFluxOnly >> Forcing input BANFF flux parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
        assert(vec.GetNoElements() == conv::kNumOfFluxParamsInBANFFfit);
	for (int i=0; i< conv::kNumOfFluxParamsInBANFFfit; i++){
	  current_values[i] = vec[i];
	}
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/TweakFluxOnly >> Generating random twk dials (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
     for (int i=0; i< conv::kNumOfParamsInBANFFfit; i++){
       if(i<conv::kNumOfFluxParamsInBANFFfit) current_values[i] = vec[i];
       else current_values[i] = bestfit_values[i];
     }
     return;
  }

}
//________________________________________________________________________________________
void BANFF::TweakRHCFluxOnly(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakFluxOnly >> Tweaking RHC BANFF flux parameters" << endl;
#endif

  // Switch flags
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdatedRHC = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fRHCBestFitValues;
  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakFluxOnly >> Reset RHC BANFF flux params to best-fit values." << endl;
#endif
        fIsTweakedRHC = false;
        current_values = bestfit_values;
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/TweakFluxOnly >> Forcing RHC input BANFF flux parameters. " << endl;
#endif
	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfRHCFluxParamsInBANFFfit);
	  for (int i=0; i< conv::kNumOfRHCFluxParamsInBANFFfit; i++){
	    current_values[i] = vec[i];
	  }
	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/TweakFluxOnly >> Generating RHC random twk dials (correlated)! " << endl;
#endif
      TMatrixD & L = *fRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      for (int i=0; i< conv::kNumOfRHCParamsInBANFFfit; i++){
	if(i<conv::kNumOfRHCFluxParamsInBANFFfit) current_values[i] = vec[i];
	else current_values[i] = bestfit_values[i];
      }
      return;
    }
}
//________________________________________________________________________________________





// NEW FLUX FHC+RHC

void BANFF::TweakFHCpRHCFluxOnly(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakFluxOnly >> Tweaking FHC+RHC BANFF flux parameters" << endl;
#endif

  // Switch flags
  fIsTweakedFHCpRHC         = true;
  fIsFluxWghtUpdatedFHCpRHC = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fFHCpRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeFHCpRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakFluxOnly >> Reset FHC+RHC BANFF flux params to best-fit values." << endl;
#endif
        fIsTweakedFHCpRHC = false;
        current_values = bestfit_values;
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/TweakFluxOnly >> Forcing FHC+RHC input BANFF flux parameters. " << endl;
#endif
	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfFHCpRHCFluxParamsInBANFFfit);
	  for (int i=0; i< conv::kNumOfFHCpRHCFluxParamsInBANFFfit; i++){
	    current_values[i] = vec[i];
	  }
	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/TweakFluxOnly >> Generating FHC+RHC random twk dials (correlated)! " << endl;
#endif
      TMatrixD & L = *fFHCpRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      for (int i=0; i< conv::kNumOfFHCpRHCParamsInBANFFfit; i++){
	if(i<conv::kNumOfFHCpRHCFluxParamsInBANFFfit) current_values[i] = vec[i];
	else current_values[i] = bestfit_values[i];
      }
      return;
    }
}
//

//________________________________________________________________________________________
void BANFF::CholeskyDecomposeMatrix()
{
#ifdef __show_cholesky_mesg_info__  
    cout << "Exp/BANFF/CholeskyDecomposeMatrix >> "
         << "Running Cholesky's decomposition method on the BANFF covariance matrix." << endl;
#endif

    TMatrixD cov(this->CovMatrix());

#ifdef __show_cholesky_mesg_debug__
    cout << "Exp/BANFF/CholeskyDecomposeMatrix >> BANFF covariance matrix: " << endl; 
    util::PrintCompact(cov,7,5);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_cholesky_mesg_debug__
    cout << "Exp/BANFF/CholeskyDecomposeMatrix >> Cholesky decomposition of BANFF covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,5);
#endif
}
//________________________________________________________________________________________
void BANFF::CholeskyDecomposeRHCMatrix()
{
#ifdef __show_cholesky_mesg_info__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> "
       << "Running Cholesky's decomposition method on the RHC BANFF covariance matrix." << endl;
#endif

  TMatrixD cov(this->RHCCovMatrix());
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(cov,7,5);
#endif

  TMatrixD L = util::CholeskyDecomposition(cov);
  fRHCCovMtrxChoL = new TMatrixD(L);
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> Cholesky decomposition of RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(*fRHCCovMtrxChoL,7,5);
#endif
}
//________________________________________________________________________________________



// NEW FLUX FHC+RHC
void BANFF::CholeskyDecomposeFHCpRHCMatrix()
{
#ifdef __show_cholesky_mesg_info__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> "
       << "Running Cholesky's decomposition method on the FHC+RHC BANFF covariance matrix." << endl;
#endif

  TMatrixD cov(this->FHCpRHCCovMatrix());
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> FHC+RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(cov,7,5);
#endif

  TMatrixD L = util::CholeskyDecomposition(cov);
  fFHCpRHCCovMtrxChoL = new TMatrixD(L);
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> Cholesky decomposition of FHC+RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(*fFHCpRHCCovMtrxChoL,7,5);
#endif
}
//

//________________________________________________________________________________________
void BANFF::RestoreNominalMC(void)
{
// A special tweak method

  cout << "Exp/BANFF >> ** Applying special tweak cancelling the effect of BANFF tuning "
       << "and restoring nominal MC..." << endl; 

  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //this is doing nothing at end, because NuFlux and NuXSec
  //are still reading bestfit values from fBestFitValues
  //we need to modify BestFitValue(p) 
  TVectorD & current_values = *fCurrentValues;
  //TVectorD & bestfit_values = *fBestFitValues;

  const int np = conv::kNumOfParamsInBANFFfit;
  for(int i=0; i<np; i++) {
    current_values[i] = 1.;
    //bestfit_values[i] = 1.;
  }
}
//________________________________________________________________________________________
void BANFF::RestoreRHCNominalMC(void)
{
  // A special tweak method

  cout << "Exp/BANFF >> ** Applying special tweak cancelling the effect of RHC BANFF tuning "
       << "and restoring RHC nominal MC..." << endl; 

  fIsTweakedRHC         = true;
  fIsFluxWghtUpdatedRHC = false;

  //this is doing nothing at end, because NuFlux and NuXSec
  //are still reading bestfit values from fBestFitValues
  //we need to modify BestFitValue(p) 
  TVectorD & current_values = *fRHCCurrentValues;
  //TVectorD & bestfit_values = *fRHCBestFitValues;

  const int np = conv::kNumOfRHCParamsInBANFFfit;
  for(int i=0; i<np; i++) {
    current_values[i] = 1.;
    //bestfit_values[i] = 1.;
  }
}
//________________________________________________________________________________________



// NEW FLUX FHC+RHC
void BANFF::RestoreFHCpRHCNominalMC(void)
{
  // A special tweak method

  cout << "Exp/BANFF >> ** Applying special tweak cancelling the effect of FHC+RHC BANFF tuning "
       << "and restoring FHC+RHC nominal MC..." << endl; 

  fIsTweakedFHCpRHC         = true;
  fIsFluxWghtUpdatedFHCpRHC = false;

  //this is doing nothing at end, because NuFlux and NuXSec
  //are still reading bestfit values from fBestFitValues
  //we need to modify BestFitValue(p) 
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  //TVectorD & bestfit_values = *fFHCpRHCBestFitValues;

  const int np = conv::kNumOfFHCpRHCParamsInBANFFfit;
  for(int i=0; i<np; i++) {
    current_values[i] = 1.;
    //bestfit_values[i] = 1.;
  }
}

//________________________________________________________________________________________
double BANFF::ChisqPenalty(bool /*all_syst*/) const
{
#if !defined(__FHC_RUNNING__) && !defined(__MULTIRING_FHC_RUNNING__)
  return 0;
#endif

  //Vector of current values & best fit values
  const TVectorD & bestfit_values   = *fBestFitValues;
  const TVectorD & current_values   = *fCurrentValues;
#ifdef __show_chisqpenalty_mesg_debug2__
    cout << "Exp/BANFF/ChisqPenalty >> BANFF parameters: " << endl;
    util::PrintCompact(current_values, bestfit_values, 12,6,"BANFF fit params");
#endif

  //covariance matrix
  const TMatrixD cov    (this->CovMatrix());
  const TMatrixD covInv (this->CovMatrixInv());
  /*
  if(!all_syst) {
    if (T2KOpt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined){
      //Vector with deviations of current values from  best fit values
      const int np = 3;
      TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
      TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
      int imaqe    = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFMAQE);
      int iccqeE0  = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFCCQEE0);
      int iBANFF1  = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFF01);
      dev [0][0] = current_values[imaqe]   - bestfit_values[imaqe];
      devT[0][0] = current_values[imaqe]   - bestfit_values[imaqe];
      dev [0][1] = current_values[iccqeE0] - bestfit_values[iccqeE0];
      devT[1][0] = current_values[iccqeE0] - bestfit_values[iccqeE0];
      dev [0][2] = current_values[iBANFF1] - bestfit_values[iBANFF1];
      devT[2][0] = current_values[iBANFF1] - bestfit_values[iBANFF1];
      
#ifdef __show_chisqpenalty_mesg_debug2__
      cout << "Exp/ChisqPenalty >> deviations in BANFF parameters: " << endl;
      util::PrintCompact(TMatrixDRow(dev,0), 12,6,"BANFF selected fit params");
#endif
      
      //fill sub cov matrix
      TMatrixD covSub(np, np);
      covSub[0][0] = cov[imaqe]   [imaqe];
      covSub[1][1] = cov[iccqeE0] [iccqeE0];
      covSub[2][2] = cov[iBANFF1] [iBANFF1];
      covSub[0][1] = cov[imaqe]   [iccqeE0];
      covSub[0][2] = cov[imaqe]   [iBANFF1];
      covSub[1][0] = cov[iccqeE0] [imaqe];
      covSub[1][2] = cov[iccqeE0] [iBANFF1];
      covSub[2][0] = cov[iBANFF1] [imaqe];
      covSub[2][1] = cov[iBANFF1] [iccqeE0];
      //invert the sub matrix
      double det;
      TMatrixD covInvSub(covSub);
      covInvSub.Invert(&det);
      
#ifdef __show_chisqpenalty_mesg_debug__
      cout << "Exp/BANFF/ChisqPenalty >> BANFF covariance sub matrix : " << endl;
      util::PrintCompact(covSub,7,5);
      cout << "Exp/BANFF/ChisqPenalty >> Inverse BANFF covariance sub matrix : " << endl;
      util::PrintCompact(covInvSub,7,5);
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
      //check the inverse
      TMatrixD prod = covSub * covInvSub;
      cout << "Exp/BANFF/ChisqPenalty >> check of product covsub*(covsub)^-1 : "<< endl;
      util::PrintCompact(prod,7,4);
#endif
      
      TMatrixD chisqpenalty =  dev * covInvSub * devT;
      double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number
      return chisq_penalty;
    }
    else{
  //Vector with deviations of current values from  best fit values
      const int np = 2;
      TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
      TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
      int imaqe    = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFMAQE);
      int iccqeE0  = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFCCQEE0);
      dev [0][0] = current_values[imaqe]   - bestfit_values[imaqe];
      devT[0][0] = current_values[imaqe]   - bestfit_values[imaqe];
      dev [0][1] = current_values[iccqeE0] - bestfit_values[iccqeE0];
      devT[1][0] = current_values[iccqeE0] - bestfit_values[iccqeE0];
      
#ifdef __show_chisqpenalty_mesg_debug2__
      cout << "Exp/ChisqPenalty >> deviations in BANFF parameters: " << endl;
      util::PrintCompact(TMatrixDRow(dev,0), 12,6,"BANFF selected fit params");
#endif
      
      //fill sub cov matrix
      TMatrixD covSub(np, np);
      covSub[0][0] = cov[imaqe]   [imaqe];
      covSub[1][1] = cov[iccqeE0] [iccqeE0];
      covSub[0][1] = cov[imaqe]   [iccqeE0];
      covSub[1][0] = cov[iccqeE0] [imaqe];
      //invert the sub matrix
      double det;
      TMatrixD covInvSub(covSub);
      covInvSub.Invert(&det);
      
#ifdef __show_chisqpenalty_mesg_debug__
      cout << "Exp/BANFF/ChisqPenalty >> BANFF covariance sub matrix : " << endl;
      util::PrintCompact(covSub,7,5);
      cout << "Exp/BANFF/ChisqPenalty >> Inverse BANFF covariance sub matrix : " << endl;
      util::PrintCompact(covInvSub,7,5);
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
      //check the inverse
      TMatrixD prod = covSub * covInvSub;
      cout << "Exp/BANFF/ChisqPenalty >> check of product covsub*(covsub)^-1 : "<< endl;
      util::PrintCompact(prod,7,4);
#endif
      
      TMatrixD chisqpenalty =  dev * covInvSub * devT;
      double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number
      return chisq_penalty;
    
    }
  }//!all_syst
  else {
  */
    //Vector with deviations of current values from  best fit values
    const int np = conv::kNumOfParamsInBANFFfit;
    TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
    TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
    for(int i = 0; i < np; i++){
      dev[0][i]  = current_values[i]-bestfit_values[i];
      devT[i][0] = current_values[i]-bestfit_values[i];
    }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in BANFF parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"BANFF fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
    cout << "Exp/BANFF/ChisqPenalty >> BANFF covariance matrix: " << endl;
    util::PrintCompact(cov,7,5);
    cout << "Exp/BANFF/ChisqPenalty >> Inverse BANFF covariance matrix: " << endl;
    util::PrintCompact(covInv,7,5);
#endif

    //Check the inverse:
#ifdef __show_chisqpenalty_mesg_debug2__
    TMatrixD prod = cov*covInv;
    cout << "Exp/BANFF/ChisqPenalty >> check of product cov*(cov)^-1 : "<< endl;
    util::PrintCompact(prod,7,5);
#endif

    TMatrixD chisqpenalty =  dev * covInv * devT;
    double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

#ifdef __show_chisqpenalty_mesg_debug__
    cout << "Exp/BANFF/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;
#endif

    if(chisq_penalty) return chisq_penalty;
    //}//all_syst

  return 0;
}
//________________________________________________________________________________________ 
double BANFF::RHCChisqPenalty(bool /*all_syst*/) const
{
#ifndef __RHC_RUNNING__
  return 0;
#endif
  //Vector of current values & best fit values
  const TVectorD & bestfit_values   = *fRHCBestFitValues;
  const TVectorD & current_values   = *fRHCCurrentValues;
#ifdef __show_chisqpenalty_mesg_debug2__
  cout << "Exp/BANFF/ChisqPenalty >> RHC BANFF parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"RHC BANFF fit params");
#endif

  //covariance matrix 
  const TMatrixD cov    (this->RHCCovMatrix());
  const TMatrixD covInv (this->RHCCovMatrixInv());

  //Vector with deviations of current values from  best fit values
  const int np = conv::kNumOfRHCParamsInBANFFfit;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in RHC BANFF parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"RHC BANFF fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
  cout << "Exp/BANFF/ChisqPenalty >> RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF/ChisqPenalty >> RHC inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(covInv,7,5);
#endif

  //Check the inverse:
#ifdef __show_chisqpenalty_mesg_debug2__
  TMatrixD prod = cov*covInv;
  cout << "Exp/BANFF/ChisqPenalty >> check of RHC product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,5);
#endif

  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/BANFF/ChisqPenalty >> RHC ChisqPenalty = " << chisq_penalty << endl;
#endif

  if(chisq_penalty) return chisq_penalty;
  //}//all_syst

  return 0;
}
//________________________________________________________________________________________ 





// NEW FLUX FHC+RHC
double BANFF::FHCpRHCChisqPenalty(bool /*all_syst*/) const
{
#if !(defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__))
  return 0;
#endif
  //Vector of current values & best fit values
  const TVectorD & bestfit_values   = *fFHCpRHCBestFitValues;
  const TVectorD & current_values   = *fFHCpRHCCurrentValues;
#ifdef __show_chisqpenalty_mesg_debug2__
  cout << "Exp/BANFF/ChisqPenalty >> FHC+RHC BANFF parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"FHC+RHC BANFF fit params");
#endif

  //covariance matrix 
  const TMatrixD cov    (this->FHCpRHCCovMatrix());
  const TMatrixD covInv (this->FHCpRHCCovMatrixInv());

  //Vector with deviations of current values from  best fit values
  const int np = conv::kNumOfFHCpRHCParamsInBANFFfit;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in FHC+RHC BANFF parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"FHC+RHC BANFF fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
  cout << "Exp/BANFF/ChisqPenalty >> FHC+RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF/ChisqPenalty >> FHC+RHC inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(covInv,7,5);
#endif

  //Check the inverse:
#ifdef __show_chisqpenalty_mesg_debug2__
  TMatrixD prod = cov*covInv;
  cout << "Exp/BANFF/ChisqPenalty >> check of FHC+RHC product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,5);
#endif

  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/BANFF/ChisqPenalty >> FHC+RHC ChisqPenalty = " << chisq_penalty << endl;
#endif

  if(chisq_penalty) return chisq_penalty;
  //}//all_syst

  return 0;
}










//________________________________________________________________________________________
void BANFF::PrintSystematicVariations(std::ostream & stream) const
{
  //assert(fCurrentValues->GetNoElements() == conv::kNumOfParamsInBANFFfit);

  const TVectorD & current_values = *fCurrentValues;
  const TVectorD & bestfit_values = *fBestFitValues;

  stream << "Exp/BANFF >> BANFF parameters: " << endl; 
  util::PrintCompact(current_values, bestfit_values, 12,6,"BANFF fit params", stream);
}
//________________________________________________________________________________________
void BANFF::PrintRHCSystematicVariations(std::ostream & stream) const
{
  //assert(fCurrentValues->GetNoElements() == conv::kNumOfRHCParamsInBANFFfit);

  const TVectorD & current_values = *fRHCCurrentValues;
  const TVectorD & bestfit_values = *fRHCBestFitValues;

  stream << "Exp/BANFF >> RHC BANFF parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"RHC BANFF fit params", stream);
}
//________________________________________________________________________________________



// NEW FLUX FHC+RHC
void BANFF::PrintFHCpRHCSystematicVariations(std::ostream & stream) const
{
  //assert(fCurrentValues->GetNoElements() == conv::kNumOfRHCParamsInBANFFfit);

  const TVectorD & current_values = *fFHCpRHCCurrentValues;
  const TVectorD & bestfit_values = *fFHCpRHCBestFitValues;

  stream << "Exp/BANFF >> FHC+RHC BANFF parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"FHC+RHC BANFF fit params", stream);
}
//



//________________________________________________________________________________________
int BANFF::Param2Idx(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fParam2Idx[p];  
}
//________________________________________________________________________________________
