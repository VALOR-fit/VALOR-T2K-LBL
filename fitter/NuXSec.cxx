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
#define __show_wght_mesg_debug__
#define __show_fsi_cholesky_mesg_info__
#define __show_fsi_cholesky_mesg_debug__
#define __show_fsi_tweak_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TH1D.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TMath.h>

#include "NuXSec.h"
#include "BANFF.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(NuXSec)

//________________________________________________________________________________________
NuXSec::NuXSec()
{
  this->Init();
}
//________________________________________________________________________________________
NuXSec::~NuXSec()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void NuXSec::Init(void)
{
  fBANFF            = 0;
  fFSICurrentValues = 0;
  fFSICovMtrx       = 0;
  fFSICovMtrxChoL   = 0;
}
//________________________________________________________________________________________
void NuXSec::CleanUp(void)
{
  if ( fFSICurrentValues ) delete fFSICurrentValues;
  if ( fFSICovMtrx       ) delete fFSICovMtrx;
  if ( fFSICovMtrxChoL   ) delete fFSICovMtrxChoL;
}
//________________________________________________________________________________________
bool NuXSec::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Reading neutrino generator inputs from input XML file." << endl;
#endif

/*
  const int np = conv::kNumOfFSIParams;

  fCurrentValues = new TVectorD(np);
  for(int ip = 0; ip < np; ip++) {
    (*fCurrentValues)[ip] = 1.;
  }

  // Get BANFF fit covariance matrix
  fCovMtrx = util::XML2TMatrixD(xml,root_node,"xsec/fsi/covariance_matrix");
  if(!fCovMtrx) {
    cerr << "Exp/NuXSec >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Input FSI covariance matrix: " << endl;
  util::PrintCompact(*fCovMtrx,7,4);
#endif

*/

#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Done reading neutrino cross-section inputs!" << endl;
#endif

  return true;  
}
//________________________________________________________________________________________
void NuXSec::TellMeAboutTheBANFFfit(BANFF * banff)
{
  fBANFF = banff;
}
//________________________________________________________________________________________
double NuXSec::SuperkXSecFinTuningWght(
   t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, double Etrue, double Ereco) const
{
//
// This method provides weights to reweight MC to the BANFF best-fit parameters.
// It handles the 4 BANFF cross-section parameters which are correlated with the BANFF 
// flux parameters. In the current incarnation of the BANFF fit, these parameters are:
// - BANFF param # 22, M_{A} QE
// - BANFF param # 23, M_{A} RES
// - BANFF param # 24, CCQE  normalization, 0 - 1.5 GeV
// - BANFF param # 25, CCRES normalization, 0 - 2.5 GeV
//

  // A check to make sure the sample-mode combination exists
  if(s == kSKsUndefined) return 1.;
  if(m == kSKmUndefined) return 1.;
  if(util::FindModeSlot(s, m) < 0) return 1.;

  double wght_maqe    = 1.;
  double wght_mares   = 1.;
  double wght_normqe  = 1.;
  double wght_normres = 1.;
/*
  bool is_ccqe  = util::IsCCQE (m);
  bool is_ccres = util::IsCCRES(m)

  // BANFF param # 22, M_{A} QE
  if(is_ccqe) {
     double maqe = conv::kNomMAQE * fBANFF->BestFitValue(kNuisBANFF22);
     maqe_wght = this->Spline(s,m,Etrue)->Eval(maqe);
  }

  // BANFF param # 23, M_{A} RES
  if(is_ccres) {
     double mares = conv::kNomMARES * fBANFF->BestFitValue(kNuisBANFF23);
     mares_wght = this->Spline(s,m,Etrue)->Eval(mares);
  }

  // BANFF param # 24, CCQE  normalization, 0 - 1.5 GeV
  if(is_ccqe && Etrue >= 0. && Etrue < 1.5) {
     wght_normqe = fBANFF->BestFitValue(kNuisBANFF24);
  }

  // BANFF param # 25, CCRES normalization, 0 - 2.5 GeV
  if(is_ccqe && Etrue >= 0. && Etrue < 1.5) {  
    wght_normres = fBANFF->BestFitValue(kNuisBANFF25);
  }
*/
  double wght = 
       wght_maqe   *
       wght_mares  *
       wght_normqe *
       wght_normres;

#ifdef __show_wght_mesg_debug__
  if ( wght_maqe    != 1. ||
       wght_mares   != 1. || 
       wght_normqe  != 1. ||
       wght_normres != 1. )
  {
     cout << "Exp/NuXSec >> "
          << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
          << ", Etrue = " << Etrue << " GeV => Weights: "
          << "MAQE = "      << wght_maqe
          << ", MARES = "   << wght_mares
          << ", NormQE = "  << wght_normqe
          << ", NormRES = " << wght_normres
          << "  (total = "  << wght << ")"
          << endl;   
  }
#endif

  return wght;
}
//________________________________________________________________________________________
double NuXSec::SuperkXSecSystWght(
   t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, double Etrue, double Ereco) const
{
  // A check to make sure the sample-mode combination exists
  if(s == kSKsUndefined) return 1.;
  if(m == kSKmUndefined) return 1.;
  if(util::FindModeSlot(s, m) < 0) return 1.;



  return 1.;
}
//________________________________________________________________________________________
double NuXSec::FSIParamCurrentValue (conv::PhysParamId_t /*p*/) const
{
  return 0;
}
//________________________________________________________________________________________
void NuXSec::TweakFSIParams(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec     /* input vector of params - Set to 0 to reset to nominal ones */
)
{
#ifdef __show_fsi_cholesky_mesg_debug__
  cout << "Exp/NuXSec/FSI/Tweak >> Tweaking BANFF fit parameters" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky
  // decomposition of the flux covariance matrix
  //

  if(fFSICovMtrxChoL == 0)
  {
#ifdef __show_fsi_cholesky_mesg_info__
    cout << "Exp/NuXSec/FSI/Tweak >> "
         << "Running Cholesky's decomposition method on the FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->FSICovMatrix());
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/NuXSec/FSI/Tweak >> FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fFSICovMtrxChoL = new TMatrixD(L);
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/NuXSec/FSI/Tweak >> Cholesky decomposition of FSI covariance matrix: " << endl;
    util::PrintCompact(*fFSICovMtrxChoL,7,4);
#endif
  }

  TVectorD & current_values = *fFSICurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
  {
     // Input vector is NULL.
     // This resets the FSI params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/NuXSec/FSI/Tweak >> Reset FSI params to nominal values." << endl;
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
        cout << "Exp/NuXSec/FSI/Tweak >> Forcing input BANFF parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
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
#ifdef __show_fsi_tweak_mesg_debug__
     cout << "Exp/NuXSec/FSI/Tweak >> Generating random FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fFSICovMtrxChoL; 
     TVectorD vec = util::CholeskyGenerateCorrelatedParamVariations(L);
     current_values = vec;
     return;
  }

}
//________________________________________________________________________________________
void NuXSec::PrintSystematicVariations(void) const
{

}
//________________________________________________________________________________________
void NuXSec::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  //
  // ...
  //
  
  c->Update();
}
//________________________________________________________________________________________
