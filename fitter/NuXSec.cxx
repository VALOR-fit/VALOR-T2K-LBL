//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Created from the old SuperkExpectation::GeneratorInputs nested class. Removed the 
   numu10a namespace.
 @ 27/07/2011 - CA
   Take into account the newly added `nue signal' templates when systematic parameters
   are tweaked. Added `CCQEWght(double Etrue) const' to avoid code-repetition in 
   `SuperkCCQEXSecWght(double Etrue, SuperkMode_t m) const'
 @ 06/10/11 - CA
   Added `PrintSystematicVariations()'. Tweaked output messages.

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_xsec_ratio_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TH1D.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TMath.h>

#include "NuXSec.h"
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
  fCCQEXSecTwkDial             = 0;
  fCC1piToCCQEXSecRatioTwkDial = 0;
  fCCCohToCCQEXSecRatioTwkDial = 0;
  fCCOthToCCQEXSecRatioTwkDial = 0;
  fNCToCCQEXSecRatioTwkDial    = 0;
  fNC1pi0ToCCQEXSecRatioTwkDial= 0;
  fNCCohToCCQEXSecRatioTwkDial = 0;
  fNCOthToCCQEXSecRatioTwkDial = 0;
  fNuEToNuMuXSecRatioTwkDial   = 0;
  fNuclEffectTwkDial           = 0;
  fCCQEXSecRelErr              = 0;
  fCC1piToCCQEXSecRatioRelErr  = 0;
  fCCCohToCCQEXSecRatioRelErr  = 0;
  fCCOthToCCQEXSecRatioRelErr  = 0;
  fNCToCCQEXSecRatioRelErr     = 0;
  fNC1pi0ToCCQEXSecRatioRelErr = 0;
  fNCCohToCCQEXSecRatioRelErr  = 0;
  fNCOthToCCQEXSecRatioRelErr  = 0;
  fNuEToNuMuXSecRatioRelErr    = 0;
  fSuperkNuclearEffectRelErr1RmuT   = 0;
  fSuperkNuclearEffectRelErr1ReTbkg = 0;
  fSuperkNuclearEffectRelErr1ReTsig = 0;
  fCC1piInCCRES                = 0;
}
//________________________________________________________________________________________
void NuXSec::CleanUp(void)
{
  if ( fCCQEXSecRelErr             ) delete fCCQEXSecRelErr;
  if ( fCC1piToCCQEXSecRatioRelErr ) delete fCC1piToCCQEXSecRatioRelErr;
  if ( fCCCohToCCQEXSecRatioRelErr ) delete fCCCohToCCQEXSecRatioRelErr;
  if ( fCCOthToCCQEXSecRatioRelErr ) delete fCCOthToCCQEXSecRatioRelErr;
  if ( fNCToCCQEXSecRatioRelErr    ) delete fNCToCCQEXSecRatioRelErr;
  if ( fNC1pi0ToCCQEXSecRatioRelErr) delete fNC1pi0ToCCQEXSecRatioRelErr;
  if ( fNCCohToCCQEXSecRatioRelErr ) delete fNCCohToCCQEXSecRatioRelErr;
  if ( fNCOthToCCQEXSecRatioRelErr ) delete fNCOthToCCQEXSecRatioRelErr;
  if ( fNuEToNuMuXSecRatioRelErr   ) delete fNuEToNuMuXSecRatioRelErr;
  if ( fSuperkNuclearEffectRelErr1RmuT   ) delete fSuperkNuclearEffectRelErr1RmuT;
  if ( fSuperkNuclearEffectRelErr1ReTbkg ) delete fSuperkNuclearEffectRelErr1ReTbkg;
  if ( fSuperkNuclearEffectRelErr1ReTsig ) delete fSuperkNuclearEffectRelErr1ReTsig;
}
//________________________________________________________________________________________
bool NuXSec::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/XSec >> Reading neutrino generator inputs from input XML file." << endl;
#endif

  fCCQEXSecRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_relative_error/ccqe");
  fCC1piToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/cc1pi_ccqe");
  fCCCohToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/cccoh_ccqe");
  fCCOthToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/ccother_ccqe");
  fNCToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nc_ccqe");
  fNC1pi0ToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nc1pi0_ccqe");
  fNCCohToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nccoh_ccqe");
  fNCOthToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/ncother_ccqe");
  fNCOthToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nue_numu");
  fSuperkNuclearEffectRelErr1RmuT =
    util::XML2TH1D(xml,root_node,"niwg/systematics/superk_nuclear_effect_relative_error_1RmuT");
  fSuperkNuclearEffectRelErr1ReTbkg =
    util::XML2TH1D(xml,root_node,"niwg/systematics/superk_nuclear_effect_relative_error_1ReTbkg");
  fSuperkNuclearEffectRelErr1ReTsig =
    util::XML2TH1D(xml,root_node,"niwg/systematics/superk_nuclear_effect_relative_error_1ReTsig");

  fCC1piInCCRES = util::XML2Dbl(xml,root_node,"neutrino_generator/cc1pi_in_ccres");
#ifdef __show_init_mesg_info__
  cout << "Exp/XSec >> CC1pi in CCRES = " <<  fCC1piInCCRES << endl;
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/XSec >> Done reading neutrino generator inputs!" << endl;
#endif
  return true;  
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
double NuXSec::ND280XSecRatioWght(
    double Etrue, t2kana::conv::ND280Mode_t mnd) const
{
   switch(mnd) {
      case ( kNDmNuMuCCQE     ) : return 1.; break;
      case ( kNDmNuMuCCRES    ) : 
        //
        // There is a mismatch between the event categories in the nd280 sample and the 
        // event categories for which errors are given.
        // I need to propagate errors on CCOther/CCQE and CC1pi/CCQE.
        // The CCRES category contains both CC1pi and CCOther.
        // I estimate that 5% of the CCRES at the T2K energies are not CC1pi,
        // but something else (multi-pi, N+gamma, etc.)
        // Use that number ...to weight the weight... so that in every true energy bin
        // we apply the CC1pi/CCQE weight to 95% of events and the CCOth/CCQE to the
        // other 5%. Although approximate, this helps assigning the proper CCOth/CCQE
        // error on nd280
	// This doesn't apply for NEUT.
        //
	return (    fCC1piInCCRES  * this->CC1piToCCQEXSecRatioWght(Etrue) +
                 (1-fCC1piInCCRES) * this->CCOthToCCQEXSecRatioWght(Etrue)
               );
	break;

      case (kNDmNuMuCCCOH) : 
        return this->CCOthToCCQEXSecRatioWght(Etrue); 
        break;

      case (kNDmNuMuCCOth) : 
        return this->CCOthToCCQEXSecRatioWght(Etrue); 
        break;

      case (kNDmNuMuNCRES) : 
        return this->NCToCCQEXSecRatioWght(Etrue);  
        break;

      case (kNDmNuMuNCCOH) : 
        return this->NCToCCQEXSecRatioWght(Etrue);  
        break;

      case (kNDmNuMuNCOth) : 
        return this->NCToCCQEXSecRatioWght(Etrue);  
        break;

      case (kNDmNuMuBarCCQE) : 
        return 1.; 
        break;

      case (kNDmNuMuBarCCRES) : 
	return (    fCC1piInCCRES  * this->CC1piToCCQEXSecRatioWght(Etrue) +
                 (1-fCC1piInCCRES) * this->CCOthToCCQEXSecRatioWght(Etrue)
               );
	break;

      case (kNDmNuMuBarCCCOH) :
        return this->CCOthToCCQEXSecRatioWght(Etrue); 
        break;

      case (kNDmNuMuBarCCOth) : 
        return this->CCOthToCCQEXSecRatioWght(Etrue); 
        break;

      case (kNDmNuMuBarNCRES) : 
        return this->NCToCCQEXSecRatioWght(Etrue); 
        break;

      case (kNDmNuMuBarNCCOH) : 
        return this->NCToCCQEXSecRatioWght(Etrue); 
        break;

      case (kNDmNuMuBarNCOth) : 
        return this->NCToCCQEXSecRatioWght(Etrue); 
        break;

      case (kNDmNuECCQE)  : 
        return 1.; 
        break;

      case (kNDmNuECCRES) :
	return (    fCC1piInCCRES  * this->CC1piToCCQEXSecRatioWght(Etrue) +
                 (1-fCC1piInCCRES) * this->CCOthToCCQEXSecRatioWght(Etrue)
               );
        break;

      case (kNDmNuECCCOH) :
         return this->CCOthToCCQEXSecRatioWght(Etrue); 
         break;

      case (kNDmNuECCOth) : 
         return this->CCOthToCCQEXSecRatioWght(Etrue); 
         break;

      case (kNDmNuENCRES) : 
         return this->NCToCCQEXSecRatioWght(Etrue); 
         break;

      case (kNDmNuENCCOH) : 
         return this->NCToCCQEXSecRatioWght(Etrue); 
         break;

      case (kNDmNuENCOth) : 
         return this->NCToCCQEXSecRatioWght(Etrue); 
         break;

      default: 
        return 1.; 
        break;
   }
   return 1.;

}
//________________________________________________________________________________________
double NuXSec::SuperkXSecRatioWght(
   double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const
{
   // a check to make sure the sample-mode combination exists
   bool mode_in_sample = false;
   for ( int im=0; im<knm; im++ )
     if ( m == kSKReactionModes[s][im] ) {
       mode_in_sample = true;
       break;
     }
   if ( !mode_in_sample )
     return 1.;
     
   switch(m) {
      case ( kSKmNuMuCCQE        ) : return 1.; break;
      case ( kSKmNuMuCC1pi       ) : return this->CC1piToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuCCcoh       ) : return this->CCCohToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuCCoth_mu    ) :
      case ( kSKmNuMuCCoth_e     ) : return this->CCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuNC          ) : return this->NCToCCQEXSecRatioWght(Etrue);    break;
      case ( kSKmNuMuNC1pi0      ) : return this->NC1pi0ToCCQEXSecRatioWght(Etrue);break;
      case ( kSKmNuMuNCcoh       ) : return this->NCCohToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuNCoth_e     ) : return this->NCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuBarCCQE     ) : return 1.; break;
      case ( kSKmNuMuBarCC1pi    ) : return this->CC1piToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuBarCCcoh    ) : return this->CCCohToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuBarCCoth_mu ) :
      case ( kSKmNuMuBarCCoth_e  ) : return this->CCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuBarNC       ) : return this->NCToCCQEXSecRatioWght(Etrue);    break;
      case ( kSKmNuMuBarNC1pi0   ) : return this->NC1pi0ToCCQEXSecRatioWght(Etrue);break;
      case ( kSKmNuMuBarNCcoh    ) : return this->NCCohToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuBarNCoth_e  ) : return this->NCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuECCQE         ) : return this->NuEToNuMuXSecRatioWght(Etrue);   break;
      case ( kSKmNuECC1pi        ) :
	{
	  double cc1pi   = this->CC1piToCCQEXSecRatioWght(Etrue);
	  double nuenumu = this->NuEToNuMuXSecRatioWght(Etrue);
	  return TMath::Sqrt( (cc1pi*cc1pi) + (nuenumu*nuenumu) );
	  break;
	}
      case ( kSKmNuECCcoh        ) :
	{
	  double cccoh = this->CCCohToCCQEXSecRatioWght(Etrue);
	  double nuenumu  = this->NuEToNuMuXSecRatioWght(Etrue);
	  return TMath::Sqrt( (cccoh*cccoh) + (nuenumu*nuenumu) );
	  break;
	}
      case ( kSKmNuECCoth_mu     ) : 
      case ( kSKmNuECCoth_e      ) :
	{
	  double ccoth = this->CCOthToCCQEXSecRatioWght(Etrue);
	  double nuenumu  = this->NuEToNuMuXSecRatioWght(Etrue);
	  return TMath::Sqrt( (ccoth*ccoth) + (nuenumu*nuenumu) );
	  break;
	}
      case ( kSKmNuENC           ) : return this->NCToCCQEXSecRatioWght(Etrue);    break;
      case ( kSKmNuENC1pi0       ) : return this->NC1pi0ToCCQEXSecRatioWght(Etrue);break;
      case ( kSKmNuENCcoh        ) : return this->NCCohToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuENCoth_e      ) : return this->NCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmSigNuECCQE      ) : return this->NuEToNuMuXSecRatioWght(Etrue);   break;
      case ( kSKmSigNuECC1pi     ) :
	{
	  double cc1pi = this->CC1piToCCQEXSecRatioWght(Etrue);
	  double nuenumu  = this->NuEToNuMuXSecRatioWght(Etrue);
	  return TMath::Sqrt( (cc1pi*cc1pi) + (nuenumu*nuenumu) );
	  break;
	}
      case ( kSKmSigNuECCoth     ) :
	{
	  double ccoth = this->CCOthToCCQEXSecRatioWght(Etrue);
	  double nuenumu  = this->NuEToNuMuXSecRatioWght(Etrue);
	  return TMath::Sqrt( (ccoth*ccoth) + (nuenumu*nuenumu) );
	  break;
	}
      default                      : return 1.; break;
   }
   return 1.;
}
//________________________________________________________________________________________
double NuXSec::SuperkCCQEXSecWght(
    double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const
{
   // a check to make sure the sample-mode combination exists
   bool mode_in_sample = false;
   for ( int im=0; im<knm; im++ )
     if ( m == kSKReactionModes[s][im] ) {
       mode_in_sample = true;
       break;
     }
   if ( !mode_in_sample )
     return 1.;

   switch(m) {
      case ( kSKmNuMuCCQE        ) : return this->CCQEWght(Etrue); break;
      case ( kSKmNuMuCC1pi       ) : return 1.; break;
      case ( kSKmNuMuCCcoh       ) : return 1.; break;
      case ( kSKmNuMuCCoth_mu    ) : return 1.; break;
      case ( kSKmNuMuCCoth_e     ) : return 1.; break;
      case ( kSKmNuMuNC          ) : return 1.; break;
      case ( kSKmNuMuNC1pi0      ) : return 1.; break;
      case ( kSKmNuMuNCcoh       ) : return 1.; break;
      case ( kSKmNuMuNCoth_e     ) : return 1.; break;
      case ( kSKmNuMuBarCCQE     ) : return this->CCQEWght(Etrue); break;
      case ( kSKmNuMuBarCC1pi    ) : return 1.; break;
      case ( kSKmNuMuBarCCcoh    ) : return 1.; break;
      case ( kSKmNuMuBarCCoth_mu ) : return 1.; break;
      case ( kSKmNuMuBarCCoth_e  ) : return 1.; break;
      case ( kSKmNuMuBarNC       ) : return 1.; break;
      case ( kSKmNuMuBarNC1pi0   ) : return 1.; break;
      case ( kSKmNuMuBarNCcoh    ) : return 1.; break;
      case ( kSKmNuMuBarNCoth_e  ) : return 1.; break;
      case ( kSKmNuECCQE         ) : return this->CCQEWght(Etrue); break;
      case ( kSKmNuECC1pi        ) : return 1.; break;
      case ( kSKmNuECCcoh        ) : return 1.; break;
      case ( kSKmNuECCoth_mu     ) : return 1.; break;
      case ( kSKmNuECCoth_e      ) : return 1.; break;
      case ( kSKmNuENC           ) : return 1.; break;
      case ( kSKmNuENC1pi0       ) : return 1.; break;
      case ( kSKmNuENCcoh        ) : return 1.; break;
      case ( kSKmNuENCoth_e      ) : return 1.; break;
      case ( kSKmSigNuECCQE      ) : return this->CCQEWght(Etrue); break;
      case ( kSKmSigNuECC1pi     ) : return 1.; break;
      case ( kSKmSigNuECCoth     ) : return 1.; break;
      default                      : return 1.; break;
   }
   return 1.;
}
//________________________________________________________________________________________
double NuXSec::SuperkNuclEffectWght(
   double Ereco, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const
{
   double mev2gev = 1E-3; //<-- note: units in MeV, all other NIWG histos have units in GeV
   TH1D * nucleff = 0;
   //get the appropriate sk fsi histogram for the sample,mode pair
   switch(s) {
      case ( kSKsFCFV1RingMuLike ) :
	assert(fSuperkNuclearEffectRelErr1RmuT); 
	nucleff = fSuperkNuclearEffectRelErr1RmuT;
	break;
      case ( kSKsFCFV1RingELike ) :
	switch(m) {
	   case ( kSKmNuMuCCQE        ) : 
   	   case ( kSKmNuMuCC1pi       ) : 
   	   case ( kSKmNuMuCCcoh       ) : 
           case ( kSKmNuMuCCoth_mu    ) : 
           case ( kSKmNuMuCCoth_e     ) : 
           case ( kSKmNuMuNC          ) : 
           case ( kSKmNuMuNC1pi0      ) : 
           case ( kSKmNuMuNCcoh       ) : 
           case ( kSKmNuMuNCoth_e     ) : 
           case ( kSKmNuMuBarCCQE     ) : 
           case ( kSKmNuMuBarCC1pi    ) : 
           case ( kSKmNuMuBarCCcoh    ) : 
           case ( kSKmNuMuBarCCoth_mu ) : 
           case ( kSKmNuMuBarCCoth_e  ) : 
           case ( kSKmNuMuBarNC       ) : 
           case ( kSKmNuMuBarNC1pi0   ) : 
           case ( kSKmNuMuBarNCcoh    ) : 
           case ( kSKmNuMuBarNCoth_e  ) : 
           case ( kSKmNuECCQE         ) : 
           case ( kSKmNuECC1pi        ) : 
           case ( kSKmNuECCcoh        ) : 
           case ( kSKmNuECCoth_mu     ) : 
           case ( kSKmNuECCoth_e      ) : 
           case ( kSKmNuENC           ) : 
           case ( kSKmNuENC1pi0       ) : 
           case ( kSKmNuENCcoh        ) : 
           case ( kSKmNuENCoth_e      ) : 
	     assert(fSuperkNuclearEffectRelErr1ReTbkg); 
	     nucleff = fSuperkNuclearEffectRelErr1ReTbkg;
	     break;
	   case ( kSKmSigNuECCQE      ) : 
	   case ( kSKmSigNuECC1pi     ) : 
	   case ( kSKmSigNuECCoth     ) : 
	     assert(fSuperkNuclearEffectRelErr1ReTsig); 
	     nucleff = fSuperkNuclearEffectRelErr1ReTsig;
	     break;	     
   	   default                      : return 1.; break;
	}//m
     default : return 1.; break;
   }//s

   double E       = TMath::Min( Ereco, nucleff->GetXaxis()->GetXmax() * mev2gev);
   int    bin     = nucleff->FindBin(E/mev2gev);
   double relerr  = nucleff->GetBinContent(bin);
   double twk     = fNuclEffectTwkDial;
   double wght    = (1. + twk*relerr);
   return wght;
}
//________________________________________________________________________________________
double NuXSec::CCQEWght(double Etrue) const
{
   assert(fCCQEXSecRelErr);
   Etrue = TMath::Min( Etrue, fCCQEXSecRelErr->GetXaxis()->GetXmax() );
   int    bin    = fCCQEXSecRelErr->FindBin(Etrue);
   double relerr = fCCQEXSecRelErr->GetBinContent(bin);
   double wght   = (1. + fCCQEXSecTwkDial*relerr);
   return wght;
}
//________________________________________________________________________________________
double NuXSec::CC1piToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fCC1piToCCQEXSecRatioTwkDial, fCC1piToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double NuXSec::CCCohToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fCCOthToCCQEXSecRatioTwkDial, fCCOthToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double NuXSec::CCOthToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fCCOthToCCQEXSecRatioTwkDial, fCCOthToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double NuXSec::NCToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNCToCCQEXSecRatioTwkDial, fNCToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double NuXSec::NC1pi0ToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNC1pi0ToCCQEXSecRatioTwkDial, fNC1pi0ToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double NuXSec::NCCohToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNCCohToCCQEXSecRatioTwkDial, fNCCohToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double NuXSec::NCOthToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNCOthToCCQEXSecRatioTwkDial, fNCOthToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double NuXSec::NuEToNuMuXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNuEToNuMuXSecRatioTwkDial, fNuEToNuMuXSecRatioRelErr);
}
//________________________________________________________________________________________
double NuXSec::XSecRatioWght(
   double E, double twk, TH1D * hRelErr) const
{
   assert(hRelErr);

   E = TMath::Min( E, hRelErr->GetXaxis()->GetXmax() );

   int    bin   = hRelErr->FindBin(E);
   double relerr = hRelErr->GetBinContent(bin);

   double wght = (1. + twk*relerr);

#ifdef __show_xsec_ratio_mesg_debug__
   cout << "Exp/XSec/XSecRatioWght >> E = " << E 
        << ", relerr = " << relerr << ", twkdial = " << twk << ", wght = " << wght 
        << endl;
#endif

   return wght; 
}
//________________________________________________________________________________________
void NuXSec::PrintSystematicVariations(void) const
{
  cout << "Exp/XSec >> Twk CCQE xsec = "             << fCCQEXSecTwkDial             << endl;
  cout << "Exp/XSec >> Twk CC1pi/CCQE xsec ratio = " << fCC1piToCCQEXSecRatioTwkDial << endl;
  cout << "Exp/XSec >> Twk CCcoh/CCQE xsec ratio = " << fCCCohToCCQEXSecRatioTwkDial << endl;
  cout << "Exp/XSec >> Twk CCoth/CCQE xsec ratio = " << fCCOthToCCQEXSecRatioTwkDial << endl;
  cout << "Exp/XSec >> Twk NC/CCQE xsec ratio = "    << fNCToCCQEXSecRatioTwkDial    << endl;
  cout << "Exp/XSec >> Twk NC1pi0/CCQE xsec ratio = "<< fNC1pi0ToCCQEXSecRatioTwkDial<< endl;
  cout << "Exp/XSec >> Twk NCCoh/CCQE xsec ratio = " << fNCCohToCCQEXSecRatioTwkDial << endl;
  cout << "Exp/XSec >> Twk NCOth/CCQE xsec ratio = " << fNCOthToCCQEXSecRatioTwkDial << endl;
  cout << "Exp/XSec >> Twk NCOth/CCQE xsec ratio = " << fNuEToNuMuXSecRatioTwkDial   << endl;
  cout << "Exp/XSec >> Twk FSI = "                   << fNuclEffectTwkDial           << endl;
}
//________________________________________________________________________________________
bool NuXSec::TweakAffectsND280 (void) const
{
  if ( TMath::Abs(fCC1piToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fCCCohToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fCCOthToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNCToCCQEXSecRatioTwkDial    ) > kASmallNum ) return true;
  if ( TMath::Abs(fNC1pi0ToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fNCCohToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNCOthToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNuEToNuMuXSecRatioTwkDial   ) > kASmallNum ) return true;
  return false;
}
//________________________________________________________________________________________
bool NuXSec::TweakAffectsSuperk(void) const
{
  if ( TMath::Abs(fCCQEXSecTwkDial             ) > kASmallNum ) return true;
  if ( TMath::Abs(fCC1piToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fCCCohToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fCCOthToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNCToCCQEXSecRatioTwkDial    ) > kASmallNum ) return true;
  if ( TMath::Abs(fNC1pi0ToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fNCCohToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNCOthToCCQEXSecRatioTwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNuEToNuMuXSecRatioTwkDial   ) > kASmallNum ) return true;
  if ( TMath::Abs(fNuclEffectTwkDial           ) > kASmallNum ) return true;
  return false;
}
//________________________________________________________________________________________


