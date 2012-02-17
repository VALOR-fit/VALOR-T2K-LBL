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
 @ 25/11/11 - TD
   Over a few iterations have:
    added new sample (e.g. NC1pi0) xsec tweaks
    changed methods 'Superk*Wght()' to take in samples as well as modes
    changed 'SuperkNuclEffectWght()' to act with different histograms for different samples/modes
    changed names: dropped XSec from CCQEXSec, dropped XSecRatio from e.g. CC1piToCCQEXSecRatio
                   changed NuclEffect to SKFSI

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
  fCCQETwkDial        = 0;
  fCC1piToCCQETwkDial = 0;
  fCCCohToCCQETwkDial = 0;
  fCCOthToCCQETwkDial = 0;
  fNCToCCQETwkDial    = 0;
  fNC1pi0ToCCQETwkDial= 0;
  fNCCohToCCQETwkDial = 0;
  fNCOthToCCQETwkDial = 0;
  fNuEToNuMuTwkDial   = 0;
  fSKFSITwkDial       = 0;
  fCCQERelErr         = 0;
  fCC1piToCCQERelErr  = 0;
  fCCCohToCCQERelErr  = 0;
  fCCOthToCCQERelErr  = 0;
  fNCToCCQERelErr     = 0;
  fNC1pi0ToCCQERelErr = 0;
  fNCCohToCCQERelErr  = 0;
  fNCOthToCCQERelErr  = 0;
  fNuEToNuMuRelErr    = 0;
  fSKFSIRelErr1Rmu    = 0;
  fSKFSIRelErr1Rebkg  = 0;
  fSKFSIRelErr1Resig  = 0;
  fCC1piInCCRES       = 0;
}
//________________________________________________________________________________________
void NuXSec::CleanUp(void)
{
  if ( fCCQERelErr        ) delete fCCQERelErr;
  if ( fCC1piToCCQERelErr ) delete fCC1piToCCQERelErr;
  if ( fCCCohToCCQERelErr ) delete fCCCohToCCQERelErr;
  if ( fCCOthToCCQERelErr ) delete fCCOthToCCQERelErr;
  if ( fNCToCCQERelErr    ) delete fNCToCCQERelErr;
  if ( fNC1pi0ToCCQERelErr) delete fNC1pi0ToCCQERelErr;
  if ( fNCCohToCCQERelErr ) delete fNCCohToCCQERelErr;
  if ( fNCOthToCCQERelErr ) delete fNCOthToCCQERelErr;
  if ( fNuEToNuMuRelErr   ) delete fNuEToNuMuRelErr;
  if ( fSKFSIRelErr1Rmu   ) delete fSKFSIRelErr1Rmu;
  if ( fSKFSIRelErr1Rebkg ) delete fSKFSIRelErr1Rebkg;
  if ( fSKFSIRelErr1Resig ) delete fSKFSIRelErr1Resig;
}
//________________________________________________________________________________________
bool NuXSec::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/XSec >> Reading neutrino generator inputs from input XML file." << endl;
#endif

  fCCQERelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_relative_error/ccqe");
  fCC1piToCCQERelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/cc1pi_ccqe");
  fCCCohToCCQERelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/cccoh_ccqe");
  fCCOthToCCQERelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/ccother_ccqe");
  fNCToCCQERelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nc_ccqe");
  fNC1pi0ToCCQERelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nc1pi0_ccqe");
  fNCCohToCCQERelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nccoh_ccqe");
  fNCOthToCCQERelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/ncother_ccqe");
  fNuEToNuMuRelErr =
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nue_numu");
  fSKFSIRelErr1Rmu =
    util::XML2TH1D(xml,root_node,"niwg/systematics/superk_nuclear_effect_relative_error/one_mu_ring");
  fSKFSIRelErr1Rebkg =
    util::XML2TH1D(xml,root_node,"niwg/systematics/superk_nuclear_effect_relative_error/one_e_ring/bkg");
  fSKFSIRelErr1Resig =
    util::XML2TH1D(xml,root_node,"niwg/systematics/superk_nuclear_effect_relative_error/one_e_ring/sig");

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
	return (    fCC1piInCCRES  * this->CC1piToCCQEWght(Etrue) +
                 (1-fCC1piInCCRES) * this->CCOthToCCQEWght(Etrue)
               );
	break;

      case (kNDmNuMuCCCOH) : 
        return this->CCOthToCCQEWght(Etrue); 
        break;

      case (kNDmNuMuCCOth) : 
        return this->CCOthToCCQEWght(Etrue); 
        break;

      case (kNDmNuMuNCRES) : 
        return this->NCToCCQEWght(Etrue);  
        break;

      case (kNDmNuMuNCCOH) : 
        return this->NCToCCQEWght(Etrue);  
        break;

      case (kNDmNuMuNCOth) : 
        return this->NCToCCQEWght(Etrue);  
        break;

      case (kNDmNuMuBarCCQE) : 
        return 1.; 
        break;

      case (kNDmNuMuBarCCRES) : 
	return (    fCC1piInCCRES  * this->CC1piToCCQEWght(Etrue) +
                 (1-fCC1piInCCRES) * this->CCOthToCCQEWght(Etrue)
               );
	break;

      case (kNDmNuMuBarCCCOH) :
        return this->CCOthToCCQEWght(Etrue); 
        break;

      case (kNDmNuMuBarCCOth) : 
        return this->CCOthToCCQEWght(Etrue); 
        break;

      case (kNDmNuMuBarNCRES) : 
        return this->NCToCCQEWght(Etrue); 
        break;

      case (kNDmNuMuBarNCCOH) : 
        return this->NCToCCQEWght(Etrue); 
        break;

      case (kNDmNuMuBarNCOth) : 
        return this->NCToCCQEWght(Etrue); 
        break;

      case (kNDmNuECCQE)  : 
        return this->NuEToNuMuWght(Etrue);
        break;

      case (kNDmNuECCRES) :
	{
	  double ccres   = (     fCC1piInCCRES  * this->CC1piToCCQEWght(Etrue) +
			      (1-fCC1piInCCRES) * this->CCOthToCCQEWght(Etrue) );
	  double nuenumu = this->NuEToNuMuWght(Etrue);
	  return ccres *  nuenumu;
	  break;
	}

      case (kNDmNuECCCOH) :
	 return ( this->CCOthToCCQEWght(Etrue) *
	          this->NuEToNuMuWght  (Etrue) );
         break;

      case (kNDmNuECCOth) : 
	 return ( this->CCOthToCCQEWght(Etrue) *
	          this->NuEToNuMuWght  (Etrue) );
         break;

      case (kNDmNuENCRES) : 
         return this->NCToCCQEWght(Etrue); 
         break;

      case (kNDmNuENCCOH) : 
         return this->NCToCCQEWght(Etrue); 
         break;

      case (kNDmNuENCOth) : 
         return this->NCToCCQEWght(Etrue); 
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
   if(s == kSKsUndefined) return 1.;
   if(m == kSKmUndefined) return 1.;
   if(util::FindModeSlot(s, m) < 0) return 1.;

   switch(m) {
      case ( kSKmNuMuCCQE        ) : return 1.; break;
      case ( kSKmNuMuCC1pi       ) : return this->CC1piToCCQEWght(Etrue); break;
      case ( kSKmNuMuCCcoh       ) : return this->CCCohToCCQEWght(Etrue); break;
      case ( kSKmNuMuCCoth_mu    ) :
      case ( kSKmNuMuCCoth_e     ) : return this->CCOthToCCQEWght(Etrue); break;
      case ( kSKmNuMuNC          ) : return this->NCToCCQEWght(Etrue);    break;
      case ( kSKmNuMuNC1pi0      ) : return this->NC1pi0ToCCQEWght(Etrue);break;
      case ( kSKmNuMuNCcoh       ) : return this->NCCohToCCQEWght(Etrue); break;
      case ( kSKmNuMuNCoth_e     ) : return this->NCOthToCCQEWght(Etrue); break;
      case ( kSKmNuMuBarCCQE     ) : return 1.; break;
      case ( kSKmNuMuBarCC1pi    ) : return this->CC1piToCCQEWght(Etrue); break;
      case ( kSKmNuMuBarCCcoh    ) : return this->CCCohToCCQEWght(Etrue); break;
      case ( kSKmNuMuBarCCoth_mu ) :
      case ( kSKmNuMuBarCCoth_e  ) : return this->CCOthToCCQEWght(Etrue); break;
      case ( kSKmNuMuBarNC       ) : return this->NCToCCQEWght(Etrue);    break;
      case ( kSKmNuMuBarNC1pi0   ) : return this->NC1pi0ToCCQEWght(Etrue);break;
      case ( kSKmNuMuBarNCcoh    ) : return this->NCCohToCCQEWght(Etrue); break;
      case ( kSKmNuMuBarNCoth_e  ) : return this->NCOthToCCQEWght(Etrue); break;
      case ( kSKmNuECCQE         ) : return this->NuEToNuMuWght(Etrue);   break;
      case ( kSKmNuECC1pi        ) :
	{
	  double cc1pi   = this->CC1piToCCQEWght(Etrue);
	  double nuenumu = this->NuEToNuMuWght(Etrue);
	  return cc1pi * nuenumu;
	  break;
	}
      case ( kSKmNuECCcoh        ) :
	{
	  double cccoh    = this->CCCohToCCQEWght(Etrue);
	  double nuenumu  = this->NuEToNuMuWght(Etrue);
	  return cccoh * nuenumu;
	  break;
	}
      case ( kSKmNuECCoth_mu     ) : 
      case ( kSKmNuECCoth_e      ) :
	{
	  double ccoth    = this->CCOthToCCQEWght(Etrue);
	  double nuenumu  = this->NuEToNuMuWght(Etrue);
	  return ccoth * nuenumu;
	  break;
	}
      case ( kSKmNuENC           ) : return this->NCToCCQEWght(Etrue);    break;
      case ( kSKmNuENC1pi0       ) : return this->NC1pi0ToCCQEWght(Etrue);break;
      case ( kSKmNuENCcoh        ) : return this->NCCohToCCQEWght(Etrue); break;
      case ( kSKmNuENCoth_e      ) : return this->NCOthToCCQEWght(Etrue); break;
      case ( kSKmSigNuECCQE      ) : return this->NuEToNuMuWght(Etrue);   break;
      case ( kSKmSigNuECC1pi     ) :
	{
	  double cc1pi    = this->CC1piToCCQEWght(Etrue);
	  double nuenumu  = this->NuEToNuMuWght(Etrue);
	  return cc1pi * nuenumu;
	  break;
	}
      case ( kSKmSigNuECCoth     ) :
	{
	  double ccoth    = this->CCOthToCCQEWght(Etrue);
	  double nuenumu  = this->NuEToNuMuWght(Etrue);
	  return ccoth * nuenumu;
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
   if(s == kSKsUndefined) return 1.;
   if(m == kSKmUndefined) return 1.;
   if(util::FindModeSlot(s, m) < 0) return 1.;

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
	assert(fSKFSIRelErr1Rmu); 
	nucleff = fSKFSIRelErr1Rmu;
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
	     assert(fSKFSIRelErr1Rebkg); 
	     nucleff = fSKFSIRelErr1Rebkg;
	     break;
	   case ( kSKmSigNuECCQE      ) : 
	   case ( kSKmSigNuECC1pi     ) : 
	   case ( kSKmSigNuECCoth     ) : 
	     assert(fSKFSIRelErr1Resig); 
	     nucleff = fSKFSIRelErr1Resig;
	     break;	     
 	   default                      : return 1.; break;
	}//m
	break;
     default : return 1.; break;
   }//s

   double E       = TMath::Min( Ereco, nucleff->GetXaxis()->GetXmax() * mev2gev);
   int    bin     = nucleff->FindBin(E/mev2gev);
   double relerr  = nucleff->GetBinContent(bin);
   double twk     = fSKFSITwkDial;
   double wght    = (1. + twk*relerr);
   return wght;
}
//________________________________________________________________________________________
double NuXSec::CCQEWght(double Etrue) const
{
   assert(fCCQERelErr);
   Etrue = TMath::Min( Etrue, fCCQERelErr->GetXaxis()->GetXmax() );
   int    bin    = fCCQERelErr->FindBin(Etrue);
   double relerr = fCCQERelErr->GetBinContent(bin);
   double wght   = (1. + fCCQETwkDial*relerr);
   return wght;
}
//________________________________________________________________________________________
double NuXSec::CC1piToCCQEWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fCC1piToCCQETwkDial, fCC1piToCCQERelErr);
}
//________________________________________________________________________________________
double NuXSec::CCCohToCCQEWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fCCCohToCCQETwkDial, fCCCohToCCQERelErr);
}
//________________________________________________________________________________________
double NuXSec::CCOthToCCQEWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fCCOthToCCQETwkDial, fCCOthToCCQERelErr);
}
//________________________________________________________________________________________
double NuXSec::NCToCCQEWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNCToCCQETwkDial, fNCToCCQERelErr);
}
//________________________________________________________________________________________
double NuXSec::NC1pi0ToCCQEWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNC1pi0ToCCQETwkDial, fNC1pi0ToCCQERelErr);
}
//________________________________________________________________________________________
double NuXSec::NCCohToCCQEWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNCCohToCCQETwkDial, fNCCohToCCQERelErr);
}
//________________________________________________________________________________________
double NuXSec::NCOthToCCQEWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNCOthToCCQETwkDial, fNCOthToCCQERelErr);
}
//________________________________________________________________________________________
double NuXSec::NuEToNuMuWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNuEToNuMuTwkDial, fNuEToNuMuRelErr);
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
  cout << "Exp/XSec >> Twk CCQE xsec = "             << fCCQETwkDial        << endl;
  cout << "Exp/XSec >> Twk CC1pi/CCQE xsec ratio = " << fCC1piToCCQETwkDial << endl;
  cout << "Exp/XSec >> Twk CCcoh/CCQE xsec ratio = " << fCCCohToCCQETwkDial << endl;
  cout << "Exp/XSec >> Twk CCoth/CCQE xsec ratio = " << fCCOthToCCQETwkDial << endl;
  cout << "Exp/XSec >> Twk NC/CCQE xsec ratio = "    << fNCToCCQETwkDial    << endl;
  cout << "Exp/XSec >> Twk NC1pi0/CCQE xsec ratio = "<< fNC1pi0ToCCQETwkDial<< endl;
  cout << "Exp/XSec >> Twk NCCoh/CCQE xsec ratio = " << fNCCohToCCQETwkDial << endl;
  cout << "Exp/XSec >> Twk NCOth/CCQE xsec ratio = " << fNCOthToCCQETwkDial << endl;
  cout << "Exp/XSec >> Twk NCOth/CCQE xsec ratio = " << fNuEToNuMuTwkDial   << endl;
  cout << "Exp/XSec >> Twk FSI = "                   << fSKFSITwkDial       << endl;
}
//________________________________________________________________________________________
bool NuXSec::TweakAffectsND280 (void) const
{
  if ( TMath::Abs(fCC1piToCCQETwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fCCOthToCCQETwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNCToCCQETwkDial    ) > kASmallNum ) return true;
  if ( TMath::Abs(fNuEToNuMuTwkDial   ) > kASmallNum ) return true;
  return false;
}
//________________________________________________________________________________________
bool NuXSec::TweakAffectsSuperk(void) const
{
  if ( TMath::Abs(fCCQETwkDial        ) > kASmallNum ) return true;
  if ( TMath::Abs(fCC1piToCCQETwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fCCCohToCCQETwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fCCOthToCCQETwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNCToCCQETwkDial    ) > kASmallNum ) return true;
  if ( TMath::Abs(fNC1pi0ToCCQETwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fNCCohToCCQETwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNCOthToCCQETwkDial ) > kASmallNum ) return true;
  if ( TMath::Abs(fNuEToNuMuTwkDial   ) > kASmallNum ) return true;
  if ( TMath::Abs(fSKFSITwkDial       ) > kASmallNum ) return true;
  return false;
}
//________________________________________________________________________________________


