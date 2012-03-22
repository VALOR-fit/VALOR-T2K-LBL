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
#define __show_init_mesg_debug__
#define __show_rc_cholesky_mesg_info__
//#define __show_rc_cholesky_mesg_debug__ 

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TMath.h>

#include "SuperkDetector.h"
#include "Utilities.h"

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
  fEff1RmuCCQEOthRelErr       = 0.;
  fEff1RmuCCQERCCovarRelErr   = 0;
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

  if(!fEff1RmuCCQERCCovarRelErr) {
    cerr << "Exp/SKDet >> No input SuperK CCQE ring-counting error covariance matrix!" << endl;
    exit(1);
  }
      
#ifdef __show_init_mesg_info__      
  cout << "Exp/SKDet >> " << "Error {EnergyScale;1Rmu}    : " << fEnergyScale1RmuRelErr  << endl;
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
double SuperkDetector::ChisqPenalty(void) const
{
  return 0;
}
//________________________________________________________________________________________
void SuperkDetector::PrintSystematicVariations(void) const
{
  cout << "Exp/SKDet >> SuperK detector systematics: " << endl
       << " > 1Rmu energy scale twk        = "  << fEnergyScaleTwkDial      << endl
       << " > 1Rmu CCQE RC0 efficiency twk = "  << fEff1RmuCCQERC0TwkDial   << endl
       << " > 1Rmu CCQE RC1 efficiency twk = "  << fEff1RmuCCQERC1TwkDial   << endl
       << " > 1Rmu CCQE RC2 efficiency twk = "  << fEff1RmuCCQERC2TwkDial   << endl
       << " > 1Rmu CCQE Oth efficiency twk = "  << fEff1RmuCCQEOthTwkDial   << endl
       << " > 1Rmu CCnonQE efficiency twk  = "  << fEff1RmuCCnonQETwkDial   << endl
       << " > 1Rmu NC efficiency twk       = "  << fEff1RmuNCTwkDial        << endl
       << " > 1Rmu nu_e CC efficiency twk  = "  << fEff1RmuNuECCTwkDial     << endl;
}
//________________________________________________________________________________________
double SuperkDetector::EfficiencyWeight(
    SuperkSample_t s, SuperkMode_t m, double Ereco) const
{
  //
  // 1-ring mu-like sample
  //
  if(s == kSKsFCFV1RingMuLike) {
    switch(m) {

       // numu CCQE, numubar CCQE:
       // (1) Apply the CCQE other efficiency error. No energy dependence.
       // (2) Apply the CCQE ring-counting efficiency error. Ereco-dependence.
       //
       case (kSKmNuMuCCQE)    : 
       case (kSKmNuMuBarCCQE) : 
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
          double weight     = weight_oth * weight_rc;
          return weight;
          break;
       }

       //
       // numu CC1pi, numu CCoth, numubar CC1pi, numubar CCoth:
       // Apply the CCnonQE efficiency error. No energy dependence.
       //
       case (kSKmNuMuCC1pi)    :
       case (kSKmNuMuCCoth_mu)    : 
       case (kSKmNuMuBarCC1pi) : 
       case (kSKmNuMuBarCCoth_mu) : 
       {
          double relerr  = fEff1RmuCCnonQERelErr;
          double twkdial = fEff1RmuCCnonQETwkDial;
          double weight  = (1. + twkdial*relerr);
          return weight;
          break;
       }

       //
       // numu/nutau NC, numubar/nutaubar NC, nue NC:
       // Apply the NC efficiency error. No energy dependence.
       //
       case (kSKmNuMuNC)        : 
       case (kSKmNuMuBarNC) : 
       case (kSKmNuENC)     :
       case (kSKmNuEBarNC)     :
       {
          double relerr  = fEff1RmuNCRelErr;
          double twkdial = fEff1RmuNCTwkDial;
          double weight  = (1. + twkdial*relerr);
          return weight;
          break;
       }

       //
       // nue CCQE, nue CC1pi, nue CCoth:
       // Apply the nue CC efficiency error. No energy dependence.
       //
       case (kSKmNuECCQE)     : 
       case (kSKmNuECC1pi)    : 
       case (kSKmNuECCoth_mu)    : 
       case (kSKmNuEBarCCQE)     : 
       case (kSKmNuEBarCC1pi)    : 
       case (kSKmNuEBarCCoth_mu)    : 
       case (kSKmSigNuECCQE)  : 
       case (kSKmSigNuECC1pi) : 
       case (kSKmSigNuECCoth) : 
       {
          double relerr  = fEff1RmuNuECCRelErr;
          double twkdial = fEff1RmuNuECCTwkDial;
          double weight  = (1. + twkdial*relerr);
          return weight;
          break;
       }

       default : 
       {
          return 0.;                 
          break;
       }
    }//mode
  }//1-ring mu-like sample
  //
  // 1-ring e-like sample
  //
  else if ( s == kSKsFCFV1RingELike ) {
    switch(m) {
       case ( kSKmNuMuCCQE        ) : 
       case ( kSKmNuMuCC1pi       ) : 
       case ( kSKmNuMuCCcoh       ) : 
       case ( kSKmNuMuCCoth_e     ) : 
       case ( kSKmNuMuNC1pi0      ) : 
       case ( kSKmNuMuNCcoh       ) : 
       case ( kSKmNuMuNCoth_e     ) : 
       case ( kSKmNuMuBarCCQE     ) : 
       case ( kSKmNuMuBarCC1pi    ) : 
       case ( kSKmNuMuBarCCcoh    ) : 
       case ( kSKmNuMuBarCCoth_e  ) : 
       case ( kSKmNuMuBarNC1pi0   ) : 
       case ( kSKmNuMuBarNCcoh    ) : 
       case ( kSKmNuMuBarNCoth_e  ) : 
       case ( kSKmNuECCQE         ) : 
       case ( kSKmNuECC1pi        ) : 
       case ( kSKmNuECCcoh        ) : 
       case ( kSKmNuECCoth_e      ) : 
       case ( kSKmNuENC1pi0       ) : 
       case ( kSKmNuENCcoh        ) : 
       case ( kSKmNuENCoth_e      ) : 
       case ( kSKmNuEBarCCQE      ) : 
       case ( kSKmNuEBarCC1pi     ) : 
       case ( kSKmNuEBarCCcoh     ) : 
       case ( kSKmNuEBarCCoth_e   ) : 
       case ( kSKmNuEBarNC1pi0    ) : 
       case ( kSKmNuEBarNCcoh     ) : 
       case ( kSKmNuEBarNCoth_e   ) : 
/////////	 {
/////////	   double relerr  = fEff1ReBkgRelErr;
/////////	   double twkdial = fEff1ReBkgTwkDial;
/////////	   double weight  = (1. + twkdial*relerr);
/////////	   return weight;
/////////	   break;
/////////	 }
       case ( kSKmSigNuECCQE      ) : 
       case ( kSKmSigNuECC1pi     ) : 
       case ( kSKmSigNuECCoth     ) : 
/////////	 {
/////////	   double relerr  = fEff1ReSigRelErr;
/////////	   double twkdial = fEff1ReSigTwkDial;
/////////	   double weight  = (1. + twkdial*relerr);
/////////	   return weight;
/////////	   break;
/////////	 }
       default:
	 return 0.;
    }//mode
  }//1-ring e-like sample

  return 0.;
}
//________________________________________________________________________________________
double SuperkDetector::EnergyScaleRelErr(SuperkSample_t s, SuperkMode_t m) const
{
// energy scale relative error

  if(s == kSKsFCFV1RingMuLike) 
  {
     return fEnergyScale1RmuRelErr;
  }
  else 
  if (s == kSKsFCFV1RingELike) {
    if( util::IsSignalNuE(m)) {
/////////      return fEnergyScale1ReSigRelErr;
    }
    else {
/////////      return fEnergyScale1ReBkgRelErr;
    }
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
