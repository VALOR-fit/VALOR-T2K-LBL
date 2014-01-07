//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

  @ 18/2/13 - LE                                                                                                              
  Added methods to handle efficiencies for 1Re (covariance matrix)

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_init_mesg_debug__
#define __show_rc_cholesky_mesg_info__
//#define __show_rc_cholesky_mesg_debug__ 
//#define __show_effic_weight_mesg_debug__
//#define __show_skdet_tweak_mesg_debug__
//#define __show_skdet_cholesky_mesg_debug__
//#define __show_skdet_cholesky_mesg_info__
//#define __show_skdet_applysyst_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TMath.h>

#include "SuperkDetector.h"
#include "Utilities.h"
#include "Opt.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(SuperkDetector)

//________________________________________________________________________________________
SuperkDetector::SuperkDetector()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkDetector::~SuperkDetector()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkDetector::Init(void)
{
  fEnergyScale1RmuRelErr      = 0.;
  fEnergyScale1ReRelErr       = 0.;
  fEff1RmuCCQEOthRelErr       = 0.;
  fEff1RmuCCQERCCovarRelErr   = 0;
  fEff1RmuCCQERCInvCovarRelErr= 0;
  fEff1RmuCCQERCCovarE1RelErr = 0.;
  fEff1RmuCCQERCCovarE2RelErr = 0.;      
  fEff1RmuCCnonQERelErr       = 0.;
  fEff1RmuNCRelErr            = 0.;
  fEff1RmuNuECCRelErr         = 0.;
  
  fEnergyScaleTwkDial         = 0.;
  fEff1RmuCCQEOthTwkDial      = 0.;
  fEff1RmuCCQERC0TwkDial      = 0.;
  fEff1RmuCCQERC1TwkDial      = 0.;
  fEff1RmuCCQERC2TwkDial      = 0.;
  fEff1RmuCCnonQETwkDial      = 0.;
  fEff1RmuNCTwkDial           = 0.;
  fEff1RmuNuECCTwkDial        = 0.;

  fEff1RmuCCQERCCovarChoLRelErr = 0;  

  fEff1ReCovarRelErr   = 0;
  fEff1ReInvCovarRelErr= 0;
  fEff1ReCovarChoLRelErr = 0;  

  f1ReCurrentValues = 0;
  f1ReNominalValues = 0;

}
//________________________________________________________________________________________
void SuperkDetector::CleanUp(void)
{

}
//________________________________________________________________________________________
bool SuperkDetector::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/SKDet >> Reading SuperK detector inputs from input XML file." << endl;
#endif


  // Get relative error on energy scale
  fEnergyScale1RmuRelErr = 
    util::XML2Dbl(xml,root_node,
     "detector_systematics/superk/energy_scale_relative_error/one_mu_ring");
  fEnergyScale1ReRelErr = 
    util::XML2Dbl(xml,root_node,
     "detector_systematics/superk/energy_scale_relative_error/one_e_ring");
  // Get relative error on various efficiencies
  fEff1RmuCCQEOthRelErr =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/other_error");
  fEff1RmuCCQERCCovarRelErr = 
    util::XML2TMatrixD(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_error");
  fEff1RmuCCQERCCovarE1RelErr =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_ereco_high_bin1");
  fEff1RmuCCQERCCovarE2RelErr =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_ereco_high_bin2");
  fEff1RmuCCnonQERelErr   =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccnonqe");
  fEff1RmuNCRelErr        =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/nc");
  fEff1RmuNuECCRelErr =
    util::XML2Dbl(xml,root_node, 
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/nuecc");
  fEff1ReCovarRelErr = 
    util::XML2TMatrixD(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_e_ring");

  if(!fEff1RmuCCQERCCovarRelErr) {
    cerr << "Exp/SKDet >> No input SuperK CCQE ring-counting error covariance matrix!" << endl;
    exit(1);
  }
  if (fEff1RmuCCQERCCovarRelErr){
    double det;
    fEff1RmuCCQERCInvCovarRelErr = new TMatrixD(*fEff1RmuCCQERCCovarRelErr);
    (*fEff1RmuCCQERCInvCovarRelErr).Invert(&det);
  }

  if(!fEff1ReCovarRelErr) {
    cerr << "Exp/SKDet >> No input SuperK 1Re efficiency covariance matrix!" << endl;
    exit(1);
  }
  if (fEff1ReCovarRelErr){
    double det;
    fEff1ReInvCovarRelErr = new TMatrixD(*fEff1ReCovarRelErr);
    (*fEff1ReInvCovarRelErr).Invert(&det);
  }
 
  // Create/init current 1Re eff params                                                                                                                                                                     
  const int np = conv::kNumOfSKDetParams1Re;
  f1ReCurrentValues = new TVectorD(np);
  f1ReNominalValues = new TVectorD(np);
  for(int ip = 0; ip < np; ip++)
    {
      (*f1ReCurrentValues)[ip] = 1.;
      (*f1ReNominalValues)[ip] = 1.;
    }


 
#ifdef __show_init_mesg_info__      
  cout << "Exp/SKDet >> " << "Error {EnergyScale;1Rmu}    : " << fEnergyScale1RmuRelErr  << endl;
  cout << "Exp/SKDet >> " << "Error {EnergyScale;1Re}    : " << fEnergyScale1ReRelErr  << endl;
  cout << "Exp/SKDet >> " << "Error {Efficiency;1Rmu/CCQE/Oth} : " << fEff1RmuCCQEOthRelErr << endl;
  cout << "Exp/SKDet >> " << "Error {Efficiency;1Rmu/CCQE/RC}  : " << endl;
  cout << " - The ring-counting efficiency error is given for the following 3 Ereco bins : ";
  cout << "(0, " << fEff1RmuCCQERCCovarE1RelErr << " GeV), "
       << "[" << fEff1RmuCCQERCCovarE1RelErr << " GeV, " << fEff1RmuCCQERCCovarE2RelErr << " GeV], "
       << "(" << fEff1RmuCCQERCCovarE2RelErr << " GeV, oo)" << endl;
  cout << " - The ring-counting error matrix is : " << endl;
  util::PrintCompact(*fEff1RmuCCQERCCovarRelErr,10,6);
  cout << "Exp/SKDet >> " << "Error {Efficiency;1Rmu/CCnonQE} : " << fEff1RmuCCnonQERelErr  << endl;
  cout << "Exp/SKDet >> " << "Error {Efficiency;1Rmu/NC} : "      << fEff1RmuNCRelErr       << endl;
  cout << "Exp/SKDet >> " << "Error {Efficiency;1Rmu/NuECC} : "   << fEff1RmuNuECCRelErr    << endl;
  cout << "Exp/SKDet >> " << "The efficiency error matrix for 1Re sample is: " << endl;
  util::PrintCompact(*fEff1ReCovarRelErr,12,12);

#endif
  
#ifdef __show_init_mesg_info__      
  cout << "Exp/SKDet >> Done reading SuperK detector inputs!" << endl;
#endif
  
  return true;
}
//________________________________________________________________________________________
void SuperkDetector::TweakEnergyScale(double twkdial) 
{ 
  fEnergyScaleTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkDetector::TweakEfficiency1RmuCCQEOther(double twkdial) 
{ 
  fEff1RmuCCQEOthTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkDetector::TweakEfficiency1RmuCCQERingCounting
(
   bool   force_input_twkdials, /* Force input twk dials? If not, generate using the Cholesky method */
   double twkdial0,             /* Value of 1st twk dial, needed if values are forced */
   double twkdial1,             /* Value of 2nd twk dial, needed if values are forced */
   double twkdial2,             /* Value of 3rd twk dial, needed if values are forced */
   bool   corellated            /* Are forced twk dials the correlated or uncorrelated ones */
)
{
#ifdef __show_rc_cholesky_mesg_debug__
  cout << "Exp/SKDet >> Tweaking CCQE ring-counting efficiency" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the ring-counting error matrix (covariance matrix)
  //
  if(!fEff1RmuCCQERCCovarChoLRelErr) 
  {
#ifdef __show_rc_cholesky_mesg_info__
    cout << "Exp/SKDet/TweakRingCountingEffCCQE >> "
         << "Run Cholesky's decomposition method on the CCQE RC error matrix." << endl;
#endif

    TMatrixD ccqe_rc_err_matrix(*fEff1RmuCCQERCCovarRelErr);
#ifdef __show_rc_cholesky_mesg_debug__
    cout << "Exp/SKDet/TweakRingCountingEffCCQE >> RC error matrix:" << endl;
    ccqe_rc_err_matrix.Print();
#endif

    TMatrixD L = util::CholeskyDecomposition(ccqe_rc_err_matrix);
    fEff1RmuCCQERCCovarChoLRelErr = new TMatrixD(L);
#ifdef __show_rc_cholesky_mesg_info__
    cout << "Exp/SKDet/TweakRingCountingEffCCQE >> "
         << "Cholesky decomposition of RC error matrix:" << endl;
    fEff1RmuCCQERCCovarChoLRelErr->Print();
#endif
  }


  //
  // Check whether the values of the tweaking dials are forced. 
  // If the forced tweaking dials are supposed to be the correlated ones then we don't
  // need to use the Cholesky decomposition of the  ring-counting error matrix.
  //
  if(force_input_twkdials && corellated) {
    fEff1RmuCCQERC0TwkDial = twkdial0;
    fEff1RmuCCQERC1TwkDial = twkdial1;
    fEff1RmuCCQERC2TwkDial = twkdial2;
    return;
  }

  TMatrixD & L = *fEff1RmuCCQERCCovarChoLRelErr;
  TVectorD twkv(3);

  //
  // Check whether the values of the tweaking dials are forced. 
  // If the forced tweaking dials are supposed to be the un-correlated ones then use
  // the Cholesky decomposition of the ring-counting error matrix to generate
  // correlated tweaking dials.
  //

  if(force_input_twkdials && !corellated) {
    TVectorD twkv_uncorrelated(3);
    twkv_uncorrelated[0] = twkdial0;
    twkv_uncorrelated[1] = twkdial1;
    twkv_uncorrelated[2] = twkdial2;
    twkv = util::CholeskyCalculateCorrelatedParamVariations(L,twkv_uncorrelated);
  } 

  //
  // If the tweaking dials were not forced, generate gaussian correlated tweaking dials.
  // Then use the Cholesky decomposition of the ring-counting error matrix to generate
  // correlated tweaking dials.
  //

  else {
    twkv = util::CholeskyGenerateCorrelatedParamVariations(L);
  }

  //
  // These generated tweaking dials are distributed around 0 and have the spread and 
  // correlation described by the input fEffCCQERelErrCovarRC error matrix.
  // Divide these tweaking dials with the sqrt(diagonal error matrix elements) so that
  // they have a `canonical' definition (similar with all other systematic tweaking dials 
  // used in GENIE: mean = 0, sigma = 1.
  //

  TMatrixD & rc_err_matrix = *fEff1RmuCCQERCCovarRelErr;
  fEff1RmuCCQERC0TwkDial = twkv(0) / TMath::Sqrt(rc_err_matrix(0,0)); 
  fEff1RmuCCQERC1TwkDial = twkv(1) / TMath::Sqrt(rc_err_matrix(1,1));
  fEff1RmuCCQERC2TwkDial = twkv(2) / TMath::Sqrt(rc_err_matrix(2,2));
}
//________________________________________________________________________________________
void SuperkDetector::TweakEfficiencies1Re(
     bool force_inputs, 
     TVectorD * input_params
)
{
#ifdef __show_skdet_cholesky_mesg_debug__
  cout << "Exp/SKDet/Tweak >> Tweaking Efficiencies for 1Re" << endl;
#endif

  TVectorD & current_values = *f1ReCurrentValues;

  //                                                                                                                                                                                              
  // Force input 1Re eff parameters                                                                                                                                                                        
  //                                                                                                                                                                                                    
  if(force_inputs)
    {
      // Input vector is NULL.     
      // This resets the 1Re eff params to their best-fit values.                                                                                                                                           
      if(input_params==0) {
#ifdef __show_skdet_tweak_mesg_debug__
        cout << "Exp/SKDet/Tweak >> Reset Efficiency params for 1Re to nominal values." << endl;
#endif
        const int np = conv::kNumOfSKDetParams1Re;
        for(int ip = 0; ip < np; ip++) {
	  current_values[ip] = 1.;
        }
        return;
      }
      else
	{
#ifdef __show_skdet_tweak_mesg_debug__
	  cout << "Exp/SKDet/Tweak >> Forcing input 1Re efficiency parameters. " << endl;
#endif
	  TVectorD & vec = *input_params;
	  current_values = vec;
	  return;
	}
    }//force_inputs?           

  //  
  // If this is the first time this function is called, then run the Cholesky  
  // decomposition of the flux covariance matrix 
  // 

  if(fEff1ReCovarChoLRelErr == 0)
    {
#ifdef __show_skdet_cholesky_mesg_info__
      cout << "Exp/SKDet/Tweak >> "
	   << "Running Cholesky's decomposition method on the 1Re efficiency covariance matrix." << endl;
#endif

      TMatrixD cov(*fEff1ReCovarRelErr);
#ifdef __show_skdet_cholesky_mesg_debug__
      cout << "Exp/SKDet/Tweak >> SKDet 1Re covariance matrix: " << endl;
      util::PrintCompact(cov,12,12);
#endif

      TMatrixD L = util::CholeskyDecomposition(cov);
      fEff1ReCovarChoLRelErr = new TMatrixD(L);
#ifdef __show_skdet_cholesky_mesg_debug__
      cout << "Exp/SKDet/Tweak >> Cholesky decomposition of 1Re efficiencies covariance matrix: " << endl;
      util::PrintCompact(*fEff1ReCovarChoLRelErr,12,12);
#endif
    }

  //                                     
  // If the tweaking dials were not forced, generate correlated tweaking dials 
  // using the Cholesky method 
  //              
  if(!force_inputs)
    {
#ifdef __show_skdet_tweak_mesg_debug__
      cout << "Exp/SKDet/Tweak >> Generating random 1Re efficiency params (correlated)! " << endl;
#endif
      TMatrixD & L = *fEff1ReCovarChoLRelErr;
      TVectorD & nominal_values = *f1ReNominalValues;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
      current_values = vec;
      return;
    }

}
//________________________________________________________________________________________
void SuperkDetector::TweakEfficiency1RmuCCnonQE(double twkdial) 
{ 
  fEff1RmuCCnonQETwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkDetector::TweakEfficiency1RmuNC(double twkdial) 
{ 
  fEff1RmuNCTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkDetector::TweakEfficiency1RmuNuECC(double twkdial) 
{ 
  fEff1RmuNuECCTwkDial = twkdial; 
}
//________________________________________________________________________________________
double SuperkDetector::ChisqPenalty(bool all_syst) const
{
  double chisq_penalty = 0;

  // when only fitting for most important systematics
  if(!all_syst) {
    if (fEff1RmuCCnonQETwkDial) chisq_penalty += fEff1RmuCCnonQETwkDial*fEff1RmuCCnonQETwkDial;
    if (fEff1RmuNCTwkDial     ) chisq_penalty += fEff1RmuNCTwkDial*fEff1RmuNCTwkDial;
    if (fEnergyScaleTwkDial   ) chisq_penalty += fEnergyScaleTwkDial*fEnergyScaleTwkDial;

    if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined){
      //add the contribution from the 1Re efficiency 
      TMatrixD & cov     = *fEff1ReCovarRelErr;
      const TVectorD & bestfit_values   = *f1ReNominalValues;
      const TVectorD & current_values   = *f1ReCurrentValues;
      
      const int np = 3;
      TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication 
      TMatrixD devT (np,1); //vector of one row in matrix form for multiplication 
      int isignue2    = 2;
      int inue2 = 8;
      int inc2  = 11;
      dev [0][0] = current_values[isignue2]   - bestfit_values[isignue2];
      devT[0][0] = current_values[isignue2]   - bestfit_values[isignue2];
      dev [0][1] = current_values[inue2] - bestfit_values[inue2];
      devT[1][0] = current_values[inue2] - bestfit_values[inue2];
      dev [0][2] = current_values[inc2] - bestfit_values[inc2];
      devT[2][0] = current_values[inc2] - bestfit_values[inc2];
      
      TMatrixD covSub(np, np);
      covSub[0][0] = cov[isignue2]   [isignue2];
      covSub[1][1] = cov[inue2] [inue2];
      covSub[2][2] = cov[inc2] [inc2];
      covSub[0][1] = cov[isignue2]   [inue2];
      covSub[0][2] = cov[isignue2]   [inc2];
      covSub[1][0] = cov[inue2] [isignue2];
      covSub[1][2] = cov[inue2] [inc2];
      covSub[2][0] = cov[inc2] [isignue2];
      covSub[2][1] = cov[inc2] [inue2];
      //invert the sub matrix 
      double det;
      TMatrixD covInvSub(covSub);
      covInvSub.Invert(&det);
      
      TMatrixD chisqpenalty1Re =  dev * covInvSub * devT;
      double dchisq_1Re = chisqpenalty1Re[0][0]; //The final matrix is 1x1, just a number
      if (dchisq_1Re) chisq_penalty += dchisq_1Re;
    } 
    return chisq_penalty;
  }

  //add the contributions from the simple multiplicative systematics
  if (fEff1RmuCCQEOthTwkDial) chisq_penalty += fEff1RmuCCQEOthTwkDial*fEff1RmuCCQEOthTwkDial;
  if (fEff1RmuCCnonQETwkDial) chisq_penalty += fEff1RmuCCnonQETwkDial*fEff1RmuCCnonQETwkDial;
  if (fEff1RmuNCTwkDial     ) chisq_penalty += fEff1RmuNCTwkDial*fEff1RmuNCTwkDial;
  if (fEff1RmuNuECCTwkDial  ) chisq_penalty += fEff1RmuNuECCTwkDial*fEff1RmuNuECCTwkDial;
  if (fEnergyScaleTwkDial   ) chisq_penalty += fEnergyScaleTwkDial*fEnergyScaleTwkDial;

  //add the contribution from the ring counting efficiency
  TMatrixD & rc_err_matrix     = *fEff1RmuCCQERCCovarRelErr;
  TMatrixD & rc_err_matrix_inv = *fEff1RmuCCQERCInvCovarRelErr;
  //  cout << "RC error matrix:" << endl;
  //rc_err_matrix.Print();
  //cout << "RC inverse error matrix:" << endl;
  //rc_err_matrix_inv.Print();
  TMatrixD prod = rc_err_matrix * rc_err_matrix_inv;
  //prod.Print();
  TMatrixD dev  (1,3); //vector of deviations in one column matrix form for multiplication
  TMatrixD devT (3,1); //vector of deviations in one row matrix form for multiplication
  dev [0][0] = fEff1RmuCCQERC0TwkDial*TMath::Sqrt(rc_err_matrix(0,0));
  dev [0][1] = fEff1RmuCCQERC1TwkDial*TMath::Sqrt(rc_err_matrix(1,1));
  dev [0][2] = fEff1RmuCCQERC2TwkDial*TMath::Sqrt(rc_err_matrix(2,2));
  for (int i = 0; i <3 ; i++){
    devT[i][0] = dev[0][i];
  }
  //cout << "Deviations" << endl;
  //dev.Print();
  //devT.Print();
  TMatrixD chisqpenalty =  dev * rc_err_matrix_inv * devT;
  double dchisq = chisqpenalty[0][0];
  if (dchisq) chisq_penalty += dchisq;

  //add the contribution from the 1Re efficiency
  TMatrixD & err_matrix_1Re     = *fEff1ReCovarRelErr;
  TMatrixD & err_matrix_inv_1Re = *fEff1ReInvCovarRelErr;
  TMatrixD prod_1Re = err_matrix_1Re * err_matrix_inv_1Re;
  //prod.Print();
  const TVectorD & bestfit_values   = *f1ReNominalValues;
  const TVectorD & current_values   = *f1ReCurrentValues;
  const int np = conv::kNumOfSKDetParams1Re;
  TMatrixD dev_1Re  (1,np); //vector of deviations in one column matrix form for multiplication
  TMatrixD devT_1Re (np,1); //vector of deviations in one row matrix form for multiplication

  for (int i = 0; i < np; i++){
    dev_1Re[0][i]  = current_values[i]-bestfit_values[i];
    devT_1Re[i][0] = current_values[i]-bestfit_values[i];
  }
  //cout << "Deviations" << endl;
  //dev.Print();
  //devT.Print();
  TMatrixD chisqpenalty_1Re =  dev_1Re * err_matrix_inv_1Re * devT_1Re;
  double dchisq_1Re = chisqpenalty_1Re[0][0];
  if (dchisq_1Re) chisq_penalty += dchisq_1Re;

  //final chisq penalty value
  return chisq_penalty;
}
//________________________________________________________________________________________
void SuperkDetector::PrintSystematicVariations(ostream & stream) const
{
  const TVectorD & curr_params = *f1ReCurrentValues;
  stream << "Exp/SKDet >> SuperK detector systematics: " << endl
	 << " > 1Rmu energy scale twk        = "  << fEnergyScaleTwkDial      << endl
	 << " > 1Rmu CCQE RC0 efficiency twk = "  << fEff1RmuCCQERC0TwkDial   << endl
	 << " > 1Rmu CCQE RC1 efficiency twk = "  << fEff1RmuCCQERC1TwkDial   << endl
	 << " > 1Rmu CCQE RC2 efficiency twk = "  << fEff1RmuCCQERC2TwkDial   << endl
	 << " > 1Rmu CCQE Oth efficiency twk = "  << fEff1RmuCCQEOthTwkDial   << endl
	 << " > 1Rmu CCnonQE efficiency twk  = "  << fEff1RmuCCnonQETwkDial   << endl
	 << " > 1Rmu NC efficiency twk       = "  << fEff1RmuNCTwkDial        << endl
	 << " > 1Rmu nu_e CC efficiency twk  = "  << fEff1RmuNuECCTwkDial     << endl;
  util::PrintCompact(curr_params,12,12," > 1Re efficiency covariance matrix");
}
//________________________________________________________________________________________
double SuperkDetector::EfficiencyWeight(
    SuperkSample_t s, SuperkMode_t m, double Ereco) const
{

  double weight = 0;

  //
  // 1-ring mu-like sample
  //
  if(s == kSKsFCFV1RingMuLike) {
#if FSI_BINS != 3
    switch(m) {

       // numu CCQE, numubar CCQE:
       // (1) Apply the CCQE other efficiency error. No energy dependence.
       // (2) Apply the CCQE ring-counting efficiency error. Ereco-dependence.
       //
       case (kSKmNuMuCCQE)    : 
       case (kSKmNuMuBarCCQE) : 
       case (kSKmMEC)         :
       {
          // other efficiency err
          double relerr_oth  = fEff1RmuCCQEOthRelErr;
          double twkdial_oth = fEff1RmuCCQEOthTwkDial;
          double weight_oth  = (1. + twkdial_oth*relerr_oth);
          // ring-counting efficiency err
          TMatrixD & rc_err_matrix = *fEff1RmuCCQERCCovarRelErr;
          double weight_rc0 = 1 + fEff1RmuCCQERC0TwkDial * TMath::Sqrt(rc_err_matrix(0,0));
          double weight_rc1 = 1 + fEff1RmuCCQERC1TwkDial * TMath::Sqrt(rc_err_matrix(1,1));
          double weight_rc2 = 1 + fEff1RmuCCQERC2TwkDial * TMath::Sqrt(rc_err_matrix(2,2));
          double weight_rc  = 1.;
          if      (Ereco < fEff1RmuCCQERCCovarE1RelErr) { weight_rc = weight_rc0; }
          else if (Ereco > fEff1RmuCCQERCCovarE2RelErr) { weight_rc = weight_rc2; }
          else                            { weight_rc = weight_rc1; }
          // combined weight
          weight     = weight_oth * weight_rc;
          break;
       }

       //
       // numu CC1pi, numu CCcoh, numu CCoth, numubar CC1pi, numubar CCcoh, numubar CCoth:
       // Apply the CCnonQE efficiency error. No energy dependence.
       //
       case (kSKmNuMuCC1pi)       :
       case (kSKmNuMuCCcoh)       :
       case (kSKmNuMuCCoth)    : 
       case (kSKmNuMuBarCC1pi)    : 
       case (kSKmNuMuBarCCcoh)    : 
       case (kSKmNuMuBarCCoth) : 
       {
          double relerr  = fEff1RmuCCnonQERelErr;
          double twkdial = fEff1RmuCCnonQETwkDial;
          weight  = (1. + twkdial*relerr);
          break;
       }

       //
       // numu/nutau NC1piC, numu/nutauNCoth_mu,
       // numubar/nutaubar NC1piC, numubar/nutaubarNCoth_mu,
       // nue NC1piC, nueNCoth_mu,
       // nuebar NC1piC, nuebarNCoth_mu,
       // Apply the NC efficiency error. No energy dependence.
       //
       case (kSKmNuMuNC1piC)      : 
       case (kSKmNuMuNCoth_mu)    : 
       case (kSKmNuMuBarNC1piC)   : 
       case (kSKmNuMuBarNCoth_mu) : 
       case (kSKmNuENC1piC)       :
       case (kSKmNuENCoth_mu)     :
       case (kSKmNuEBarNC1piC)    :
       case (kSKmNuEBarNCoth_mu)  :
       {
          double relerr  = fEff1RmuNCRelErr;
          double twkdial = fEff1RmuNCTwkDial;
          weight  = (1. + twkdial*relerr);
          break;
       }

       //
       // nue CCQE, nue CC1pi, nue CCcoh, nue CCoth:
       // Apply the nue CC efficiency error. No energy dependence.
       //
       case (kSKmNuECCQE)        : 
       case (kSKmNuECC1pi)       : 
       case (kSKmNuECCcoh)       : 
       case (kSKmNuECCoth)    : 
       case (kSKmNuEBarCCQE)     : 
       case (kSKmNuEBarCC1pi)    : 
       case (kSKmNuEBarCCcoh)    : 
       case (kSKmNuEBarCCoth) : 
       case (kSKmOscNuECCQE)     : 
       case (kSKmOscNuECC1pi)    : 
       case (kSKmOscNuECCcoh)    : 
       case (kSKmOscNuECCoth) : 
       {
          double relerr  = fEff1RmuNuECCRelErr;
          double twkdial = fEff1RmuNuECCTwkDial;
          weight  = (1. + twkdial*relerr);
          break;
       }

       default : 
       {
          return 0.;                 
          break;
       }
       
       return 0;

    }//mode
#endif
  }//1-ring mu-like sample
  //
  // 1-ring e-like sample
  //
  else if ( s == kSKsFCFV1RingELike ) {
#if FSI_BINS != 2
    // find which sk det param we have to apply
    int iparam = -1;
    for(int ip = 0; ip < conv::kNumOfSKDetParams1Re; ip++) {
      conv::PhysParamId_t p = conv::kListOfSKDet1ReParams[ip];
      //      cout << "Param = " << util::AsString(p) << endl;
      // Is the current systematic parameter applied to current sample/mode end energy bin? 
      bool apply = this->ApplySyst1Re(p,s,m,Ereco);
      if(!apply) continue;
      if(iparam != -1) {
	cerr << "Exp/SKDet/EfficiencyWght >> Found multiple SK 1Re eff params for sample, mode, Ereco combination!" << endl
	     << util::AsString(s) << "\t" << util::AsString(m) << "\t" << Ereco << endl
	     << util::AsString((conv::PhysParamId_t)(iparam+conv::kNuisSuperkEff1ReNueSigCC0)) << "\t" <<  util::AsString(p) << endl;
	exit(1);
      }
      iparam = util::ArrIdx::Instance()->SKDet1ReSystIdx(p);
    }
    //    cout << "iparam =" << iparam << endl;
    const int np = conv::kNumOfSKDetParams1Re;
    if(iparam < 0 || iparam >= np) {
      cerr << "Exp/SKDet/EfficiencyWght >> Weird index: " << iparam << " for Ereco value = " << Ereco
	   << " sample " << util::AsString(s) << " mode " << util::AsString(m) << endl;
      exit(1);
    }
    // get current value of iparam^th SK 1Re eff param 
    const TVectorD & current_values   = *f1ReCurrentValues;
    //    current_values.Print();
    double param_value = current_values[iparam];
    //    cout << "param_value=" << param_value << endl;

    // param tweaked           
    bool tweaked = TMath::Abs(param_value - 1.) > conv::kASmallNum;
    if(!tweaked) {
      weight  = 1.0;
    }

    double eff_1Re_wght = param_value;
    weight  = eff_1Re_wght;
    
    //    cout << "1Re eff weight = " << weight << endl;
#endif
  }//1-ring e-like sample

  return weight;

}
//________________________________________________________________________________________
double SuperkDetector::EnergyScaleRelErr(SuperkSample_t s, SuperkMode_t /*m*/) const
{
// energy scale relative error

  if(s == kSKsFCFV1RingMuLike) 
  {
     return fEnergyScale1RmuRelErr;
  }
  else 
  if (s == kSKsFCFV1RingELike) {
    return fEnergyScale1RmuRelErr;    
  }
  return 0; 
}
//________________________________________________________________________________________
double SuperkDetector::EnergyScale(SuperkSample_t s, SuperkMode_t m) const
{
// energy scale

  double twkdial = fEnergyScaleTwkDial;
  double err     = this->EnergyScaleRelErr(s,m);
  double scale   = 1. + twkdial * err;
  return scale;
}
//________________________________________________________________________________________
double SuperkDetector::EnergyScale(SuperkSample_t s) const
{
// energy scale

  double twkdial = fEnergyScaleTwkDial;
  double err     = this->EnergyScaleRelErr(s,kSKmNuMuCCQE);
  double scale   = 1. + twkdial * err;
  return scale;
}
//________________________________________________________________________________________
double SuperkDetector::CurrentValue(PhysParamId_t p) const
{
  switch(p) {
    case kNuisSuperkEnergyScale:
      return fEnergyScaleTwkDial;
      break;
#if FSI_BINS != 3
    case kNuisSuperkEff1RmuCCQEoth:
      return fEff1RmuCCQEOthTwkDial;
      break;
    case kNuisSuperkEff1RmuCCQErc0u:
      return fEff1RmuCCQERC0TwkDial;
      break;
    case kNuisSuperkEff1RmuCCQErc1u:
      return fEff1RmuCCQERC1TwkDial;
      break;
    case kNuisSuperkEff1RmuCCQErc2u:
      return fEff1RmuCCQERC2TwkDial;
      break;
    case kNuisSuperkEff1RmuCCnonQE:
      return fEff1RmuCCnonQETwkDial;
      break;
    case kNuisSuperkEff1RmuNC:
      return fEff1RmuNCTwkDial;
      break;
    case kNuisSuperkEff1RmuNueCC:
      return fEff1RmuNuECCTwkDial;
      break;
#endif
#if FSI_BINS != 2
    case kNuisSuperkEff1ReNueSigCC0:
    case kNuisSuperkEff1ReNueSigCC1:
    case kNuisSuperkEff1ReNueSigCC2:
    case kNuisSuperkEff1ReNumuCC0:
    case kNuisSuperkEff1ReNumuCC1:
    case kNuisSuperkEff1ReNumuCC2:
    case kNuisSuperkEff1ReNueCC0:
    case kNuisSuperkEff1ReNueCC1:
    case kNuisSuperkEff1ReNueCC2:
    case kNuisSuperkEff1ReNC0:
    case kNuisSuperkEff1ReNC1:
    case kNuisSuperkEff1ReNC2:
      {
      const TVectorD & current_values = *f1ReCurrentValues;
      int idx = util::ArrIdx::Instance()->SKDet1ReSystIdx(p);
      return current_values[idx];
      break;
      }
#endif
    default:
      return 0;
  }
  return 0;
}
//________________________________________________________________________________________
#if FSI_BINS != 2
bool SuperkDetector::ApplySyst1Re(
		    conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const
{
  bool apply2mode = this->ApplySystToSKMode(p,s,m);
  bool apply2bin  = this->ApplySystToEnergyBin(p,Ereco);
  bool apply = apply2mode && apply2bin;
  return apply;
}
//#endif
//________________________________________________________________________________________                                                                                                               
bool SuperkDetector::ApplySystToSKMode(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const
{
  int idx = util::ArrIdx::Instance()->SKDet1ReSystIdx(p);
  if (idx == -1) return false;
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if (is == -1) return false;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if (im == -1) return false;

  bool apply = conv::kSKDet1ReApplyParamToSKMode[is][im][idx];
#ifdef __show_skdet_applysyst_mesg_debug__
  string applystr;
  applystr = apply ? "Passed" : "Failed";
  cout << applystr << "\t" << util::AsString(p) << "\t" << util::AsString(s) << "\t" << util::AsString(m) << endl;
#endif
  return apply;
}
//________________________________________________________________________________________                                                                                                               
bool SuperkDetector::ApplySystToEnergyBin(conv::PhysParamId_t p, double Ereco) const
{
  int idx = util::ArrIdx::Instance()->SKDet1ReSystIdx(p);
  if (idx == -1) return false;

  bool apply = Ereco >= conv::kSKDet1ReParamsErecoCut[idx][0] &&
    Ereco <  conv::kSKDet1ReParamsErecoCut[idx][1];
#ifdef __show_skdet_applysyst_mesg_debug__
  string applystr;
  applystr = apply ? "Passed" : "Failed";
  cout << applystr << "\t" << util::AsString(p) << "\t" << Ereco << endl;
#endif
  return apply;
}
#endif
//________________________________________________________________________________________                                                                                                               
