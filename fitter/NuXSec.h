//____________________________________________________________________________
/*!
   
\class    t2kana::NuXSec
      
\brief    Neutrino generator / cross section info used in the oscillation fit
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________
   
#ifndef _t2kana_nu_xsec_h_
#define _t2kana_nu_xsec_h_

#include <TObject.h>
#include <TXMLEngine.h>

#include "Conventions.h"

class TH1D;
class TVirtualPad;

namespace t2kana {

class NuXSec: public TObject
{
public:
   NuXSec();
  ~NuXSec();

   // Get neutrino generator oscillation analysis inputs
   double ND280XSecRatioWght  (double Etrue, t2kana::conv::ND280Mode_t  m) const;
   double SuperkXSecRatioWght (double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;
   double SuperkCCQEXSecWght  (double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;
   double SuperkNuclEffectWght(double Ereco, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;

   // Systematic variations
   void TweakCCQEXSec             (double twkdial) { fCCQEXSecTwkDial             = twkdial; } ///< tweak CCQE (far only)
   void TweakCC1piToCCQEXSecRatio (double twkdial) { fCC1piToCCQEXSecRatioTwkDial = twkdial; } ///< tweak CC1pi/CCQE (both near & far)
   void TweakCCCohToCCQEXSecRatio (double twkdial) { fCCCohToCCQEXSecRatioTwkDial = twkdial; } ///< tweak CCcoh/CCQE (both near & far)
   void TweakCCOthToCCQEXSecRatio (double twkdial) { fCCOthToCCQEXSecRatioTwkDial = twkdial; } ///< tweak CCoth/CCQE (both near & far)
   void TweakNCToCCQEXSecRatio    (double twkdial) { fNCToCCQEXSecRatioTwkDial    = twkdial; } ///< tweak NC/CCQE (both near & far)
   void TweakNC1pi0ToCCQEXSecRatio(double twkdial) { fNC1pi0ToCCQEXSecRatioTwkDial= twkdial; } ///< tweak NC1pi0/CCQE (both near & far)
   void TweakNCCohToCCQEXSecRatio (double twkdial) { fNCCohToCCQEXSecRatioTwkDial = twkdial; } ///< tweak NCCoh/CCQE (both near & far)
   void TweakNCOthToCCQEXSecRatio (double twkdial) { fNCOthToCCQEXSecRatioTwkDial = twkdial; } ///< tweak NCOth/CCQE (both near & far)
   void TweakNuEToNuMuXSecRatio   (double twkdial) { fNuEToNuMuXSecRatioTwkDial   = twkdial; } ///< tweak sigma(nue)/sigma(numu) (both near & far) - only apply to nueCC templates (TN52, pg20)
   void TweakNuclEffect           (double twkdial) { fNuclEffectTwkDial           = twkdial; } ///< tweak nuclear effect (far only)
   bool TweakAffectsND280         (void) const;
   bool TweakAffectsSuperk        (void) const;
  
   // Read inputs from XML file
   bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

   // Draw inputs:
   //
   void Draw(TVirtualPad * c, Option_t * opt = "") const;

   // Print
   void PrintSystematicVariations (void) const;

   // Etc 
   double OnePionInRES(void) { return fCC1piInCCRES; }
   
private:

   void Init    (void);
   void CleanUp (void);

   double CCQEWght                  (double Etrue) const;
   double CC1piToCCQEXSecRatioWght  (double Etrue) const;
   double CCCohToCCQEXSecRatioWght  (double Etrue) const;
   double CCOthToCCQEXSecRatioWght  (double Etrue) const;
   double NCToCCQEXSecRatioWght     (double Etrue) const;
   double NC1pi0ToCCQEXSecRatioWght (double Etrue) const;
   double NCCohToCCQEXSecRatioWght  (double Etrue) const;
   double NCOthToCCQEXSecRatioWght  (double Etrue) const;
   double NuEToNuMuXSecRatioWght    (double Etrue) const;
   double XSecRatioWght             (double Etrue, double twkdial, TH1D * hRelErr) const;

   double fCCQEXSecTwkDial;               ///<
   double fCC1piToCCQEXSecRatioTwkDial;   ///<
   double fCCCohToCCQEXSecRatioTwkDial;   ///<
   double fCCOthToCCQEXSecRatioTwkDial;   ///<
   double fNCToCCQEXSecRatioTwkDial;      ///<
   double fNC1pi0ToCCQEXSecRatioTwkDial;  ///<
   double fNCCohToCCQEXSecRatioTwkDial;   ///<
   double fNCOthToCCQEXSecRatioTwkDial;   ///<
   double fNuEToNuMuXSecRatioTwkDial;     ///<
   double fNuclEffectTwkDial;             ///<
   TH1D * fCCQEXSecRelErr;                ///<
   TH1D * fCC1piToCCQEXSecRatioRelErr;    ///<
   TH1D * fCCCohToCCQEXSecRatioRelErr;    ///<
   TH1D * fCCOthToCCQEXSecRatioRelErr;    ///<
   TH1D * fNCToCCQEXSecRatioRelErr;       ///<
   TH1D * fNC1pi0ToCCQEXSecRatioRelErr;   ///<
   TH1D * fNCCohToCCQEXSecRatioRelErr;    ///<
   TH1D * fNCOthToCCQEXSecRatioRelErr;    ///<
   TH1D * fNuEToNuMuXSecRatioRelErr;      ///<
   TH1D * fSuperkNuclearEffectRelErr1RmuT;///<
   TH1D * fSuperkNuclearEffectRelErr1ReTbkg;///<
   TH1D * fSuperkNuclearEffectRelErr1ReTsig;///<
   double fCC1piInCCRES;                  ///< correct ND280 interaction mode / NWIG err mis-match

ClassDef(NuXSec, 3);
};

} // t2kana

#endif
