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
*/
//________________________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TH1D.h>
#include <TVirtualPad.h>
#include <TLegend.h>

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
  fCCOthToCCQEXSecRatioTwkDial = 0;
  fNCToCCQEXSecRatioTwkDial    = 0;
  fNuclEffectTwkDial           = 0;
  fCCQEXSecRelErr              = 0;
  fCC1piToCCQEXSecRatioRelErr  = 0;
  fCCOthToCCQEXSecRatioRelErr  = 0;
  fNCToCCQEXSecRatioRelErr     = 0;
  fSuperkNuclearEffectRelErr   = 0;
  fCC1piInCCRES                = 0;
}
//________________________________________________________________________________________
void NuXSec::CleanUp(void)
{
  if ( fCCQEXSecRelErr             ) delete fCCQEXSecRelErr;
  if ( fCC1piToCCQEXSecRatioRelErr ) delete fCC1piToCCQEXSecRatioRelErr;
  if ( fCCOthToCCQEXSecRatioRelErr ) delete fCCOthToCCQEXSecRatioRelErr;
  if ( fNCToCCQEXSecRatioRelErr    ) delete fNCToCCQEXSecRatioRelErr;
  if ( fSuperkNuclearEffectRelErr  ) delete fSuperkNuclearEffectRelErr;
}
//________________________________________________________________________________________
bool NuXSec::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  cout << "Exp/NuGen >> Reading neutrino generator inputs from input XML file." << endl;

  fCCQEXSecRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_relative_error/ccqe");
  fCC1piToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/cc1pi_ccqe");
  fCCOthToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/ccother_ccqe");
  fNCToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nc_ccqe");
  fSuperkNuclearEffectRelErr =
    util::XML2TH1D(xml,root_node,"niwg/systematics/superk_nuclear_effect_relative_error");

  fCC1piInCCRES = util::XML2Dbl(xml,root_node,"neutrino_generator/cc1pi_in_ccres");
  cout << "Exp/NuGen >> CC1pi in CCRES = " <<  fCC1piInCCRES << endl;

  cout << "Exp/NuGen >> Done reading neutrino generator inputs!" << endl;

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
   double Etrue, t2kana::conv::SuperkMode_t m) const
{
   switch(m) {
      case ( kSKmNuMuCCQE       ) : return 1.; break;
      case ( kSKmNuMuCC1pi      ) : return this->CC1piToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuCCoth      ) : return this->CCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNC             ) : return this->NCToCCQEXSecRatioWght(Etrue);    break;
      case ( kSKmNuMuBarCCQE    ) : return 1.; break;
      case ( kSKmNuMuBarCC1pi   ) : return this->CC1piToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuBarCCoth   ) : return this->CCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuBarNC      ) : return this->NCToCCQEXSecRatioWght(Etrue);    break;
      case ( kSKmNuECCQE        ) : return 1.; break;
      case ( kSKmNuECC1pi       ) : return this->CC1piToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuECCoth       ) : return this->CCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuENC          ) : return this->NCToCCQEXSecRatioWght(Etrue);    break;
      case ( kSKmSigNuECCQE     ) : return 1.; break;
      case ( kSKmSigNuECC1pi    ) : return this->CC1piToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmSigNuECCoth    ) : return this->CCOthToCCQEXSecRatioWght(Etrue); break;
      default                     : return 1.; break;
   }
   return 1.;
}
//________________________________________________________________________________________
double NuXSec::SuperkCCQEXSecWght(
    double Etrue, t2kana::conv::SuperkMode_t m) const
{
   switch(m) {
      case ( kSKmNuMuCCQE     ) : return this->CCQEWght(Etrue); break;
      case ( kSKmNuMuCC1pi    ) : return 1.;                    break;
      case ( kSKmNuMuCCoth    ) : return 1.;                    break;
      case ( kSKmNC           ) : return 1.;                    break;
      case ( kSKmNuMuBarCCQE  ) : return this->CCQEWght(Etrue); break;
      case ( kSKmNuMuBarCC1pi ) : return 1.;                    break;
      case ( kSKmNuMuBarCCoth ) : return 1.;                    break;
      case ( kSKmNuMuBarNC    ) : return 1.;                    break;
      case ( kSKmNuECCQE      ) : return this->CCQEWght(Etrue); break;
      case ( kSKmNuECC1pi     ) : return 1.;                    break;
      case ( kSKmNuECCoth     ) : return 1.;                    break;
      case ( kSKmNuENC        ) : return 1.;                    break;
      case ( kSKmSigNuECCQE   ) : return this->CCQEWght(Etrue); break;
      case ( kSKmSigNuECC1pi  ) : return 1.;                    break;
      case ( kSKmSigNuECCoth  ) : return 1.;                    break;
      default                   : return 1.;                    break;
   }
   return 1.;
}
//________________________________________________________________________________________
double NuXSec::SuperkNuclEffectWght(
   double Ereco, t2kana::conv::SuperkMode_t /*m*/) const
{
   assert(fSuperkNuclearEffectRelErr); //<-- note: units in MeV, all other NIWG histos have units in GeV
   double mev2gev = 1E-3;
   double E       = TMath::Min( Ereco, fSuperkNuclearEffectRelErr->GetXaxis()->GetXmax() * mev2gev);
   int    bin     = fSuperkNuclearEffectRelErr->FindBin(E/mev2gev);
   double relerr  = fSuperkNuclearEffectRelErr->GetBinContent(bin);
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
double NuXSec::XSecRatioWght(
   double E, double twk, TH1D * hRelErr) const
{
   assert(hRelErr);

   E = TMath::Min( E, hRelErr->GetXaxis()->GetXmax() );

   int    bin   = hRelErr->FindBin(E);
   double relerr = hRelErr->GetBinContent(bin);

   double wght = (1. + twk*relerr);

#ifdef __debug_level_3__
   cout << "Exp/NuGen/XSecRatioWght >> E = " << E 
        << ", relerr = " << relerr << ", twkdial = " << twk << ", wght = " << wght 
        << endl;
#endif

   return wght; 
}
//________________________________________________________________________________________
bool NuXSec::TweakAffectsND280 (void) const
{
  if ( TMath::Abs(fCC1piToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fCCOthToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fNCToCCQEXSecRatioTwkDial   ) > kASmallNum ) return true;
  return false;
}
//________________________________________________________________________________________
bool NuXSec::TweakAffectsSuperk(void) const
{
  if ( TMath::Abs(fCCQEXSecTwkDial            ) > kASmallNum ) return true;
  if ( TMath::Abs(fCC1piToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fCCOthToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fNCToCCQEXSecRatioTwkDial   ) > kASmallNum ) return true;
  if ( TMath::Abs(fNuclEffectTwkDial          ) > kASmallNum ) return true;
  return false;
}
//________________________________________________________________________________________


