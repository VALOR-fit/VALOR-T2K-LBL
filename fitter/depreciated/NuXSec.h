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
   void TweakCCQE          (double twkdial) { fCCQETwkDial        = twkdial; } ///< tweak CCQE (far only)
   void TweakCC1piToCCQE   (double twkdial) { fCC1piToCCQETwkDial = twkdial; } ///< tweak CC1pi/CCQE (both near & far)
   void TweakCCCohToCCQE   (double twkdial) { fCCCohToCCQETwkDial = twkdial; } ///< tweak CCcoh/CCQE (both near & far)
   void TweakCCOthToCCQE   (double twkdial) { fCCOthToCCQETwkDial = twkdial; } ///< tweak CCoth/CCQE (both near & far)
   void TweakNCToCCQE      (double twkdial) { fNCToCCQETwkDial    = twkdial; } ///< tweak NC/CCQE (both near & far)
   void TweakNC1pi0ToCCQE  (double twkdial) { fNC1pi0ToCCQETwkDial= twkdial; } ///< tweak NC1pi0/CCQE (both near & far)
   void TweakNCCohToCCQE   (double twkdial) { fNCCohToCCQETwkDial = twkdial; } ///< tweak NCCoh/CCQE (both near & far)
   void TweakNCOthToCCQE   (double twkdial) { fNCOthToCCQETwkDial = twkdial; } ///< tweak NCOth/CCQE (both near & far)
   void TweakNuEToNuMu     (double twkdial) { fNuEToNuMuTwkDial   = twkdial; } ///< tweak sigma(nue)/sigma(numu) (both near & far) - only apply to nueCC templates (TN52, pg20)
   void TweakSKFSI         (double twkdial) { fSKFSITwkDial       = twkdial; } ///< tweak nuclear effect (far only)
   bool TweakAffectsND280  (void) const;
   bool TweakAffectsSuperk (void) const;
  
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

   double CCQEWght         (double Etrue) const;
   double CC1piToCCQEWght  (double Etrue) const;
   double CCCohToCCQEWght  (double Etrue) const;
   double CCOthToCCQEWght  (double Etrue) const;
   double NCToCCQEWght     (double Etrue) const;
   double NC1pi0ToCCQEWght (double Etrue) const;
   double NCCohToCCQEWght  (double Etrue) const;
   double NCOthToCCQEWght  (double Etrue) const;
   double NuEToNuMuWght    (double Etrue) const;
   double XSecRatioWght    (double Etrue, double twkdial, TH1D * hRelErr) const;

   double fCCQETwkDial;          ///<
   double fCC1piToCCQETwkDial;   ///<
   double fCCCohToCCQETwkDial;   ///<
   double fCCOthToCCQETwkDial;   ///<
   double fNCToCCQETwkDial;      ///<
   double fNC1pi0ToCCQETwkDial;  ///<
   double fNCCohToCCQETwkDial;   ///<
   double fNCOthToCCQETwkDial;   ///<
   double fNuEToNuMuTwkDial;     ///<
   double fSKFSITwkDial;         ///<
   TH1D * fCCQERelErr;           ///<
   TH1D * fCC1piToCCQERelErr;    ///<
   TH1D * fCCCohToCCQERelErr;    ///<
   TH1D * fCCOthToCCQERelErr;    ///<
   TH1D * fNCToCCQERelErr;       ///<
   TH1D * fNC1pi0ToCCQERelErr;   ///<
   TH1D * fNCCohToCCQERelErr;    ///<
   TH1D * fNCOthToCCQERelErr;    ///<
   TH1D * fNuEToNuMuRelErr;      ///<
   TH1D * fSKFSIRelErr1Rmu;      ///<
   TH1D * fSKFSIRelErr1Rebkg;    ///<
   TH1D * fSKFSIRelErr1Resig;    ///<
   double fCC1piInCCRES;         ///< correct ND280 interaction mode / NWIG err mis-match

ClassDef(NuXSec, 3);
};

} // t2kana

#endif
