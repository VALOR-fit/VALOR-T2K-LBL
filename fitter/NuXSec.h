//____________________________________________________________________________
/*!
   
\class    t2kana::NuXSec
      
\brief    Neutrino cross-section info used in the oscillation fit
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________
   
#ifndef _t2kana_nu_xsec_h_
#define _t2kana_nu_xsec_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TH2.h>
#include <TSpline.h>
#include <TFile.h>

#include "Conventions.h"
#include "Utilities.h"

using t2kana::conv::kns;
using t2kana::conv::knm;
using t2kana::conv::knwspl1Rmu;
using t2kana::conv::knwspl1Re;
using t2kana::conv::kRHCnwspl1Rmu;
using t2kana::conv::kRHCnwspl1Re;

namespace t2kana {

class BANFF;

class NuXSec: public TObject
{
public:

  NuXSec();
 ~NuXSec();

  // Reweight to best-fit values of the BANFF cross-section related params
  double FinTuningWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Calculate weight in response to systematic parameter variations of cross-section parameters
  // For the BANFF parameters, the weight is relative to the best-fit (tuned MC) parameter, not the nominal MC one
  double SystWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Tweak systematic parameters held in this object.
  // Note: BANFF cross-section systematic parameters are tweaked via the BANFF object and they are 
  // correlated with the flux systematic params.
  void TweakFermiMom      (double twkdial) { fFermiMomTwkDial      = twkdial; }
  void TweakSpecFunc      (double twkdial) { fSpecFuncTwkDial      = twkdial; }
  void TweakWshape        (double twkdial) { fWshapeTwkDial        = twkdial; }
  void TweakShapeCCoth    (double twkdial) { fShapeCCothTwkDial    = twkdial; }
  void TweakPiLessDelta   (double twkdial) { fPiLessDeltaTwkDial   = twkdial; }
  void TweakBindE         (double twkdial) { fBindETwkDial         = twkdial; }
  void TweakNormCCcoh     (double twkdial) { fNormCCcohTwkDial     = twkdial; }
  void TweakNormNC1piC    (double twkdial) { fNormNC1piCTwkDial    = twkdial; }
  void TweakNormNCcoh     (double twkdial) { fNormNCcohTwkDial     = twkdial; }
  void TweakNormNCoth     (double twkdial) { fNormNCothTwkDial     = twkdial; }
  void TweakNormNuEToNuMu (double twkdial) { fNormNuEToNuMuTwkDial = twkdial; }
  void TweakNormAntiNuToNu(double twkdial) { fNormAntiNuToNuTwkDial= twkdial; }
  
  void Tweak(conv::PhysParamId_t p , double twkdial);

  // Get current value of systematic params handled by this class
  double CurrentValue        (conv::PhysParamId_t p) const;
  // Used by SuperkMC to get hold of the BANFF best fit value for MEC norm parameter
  double MECBestFitValue     () const;

  // Get the relative error of systematic params handled by this class
  double RelativeError (conv::PhysParamId_t p) const;

  // Get weight spline for each sample / mode / parameters / Ereco bin / Etrue bin
  TSpline3 * FindWghtSpline(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Print
  void PrintSystematicVariations (ostream & stream = std::cout) const;

  // Read inputs from XML file
  bool Read        (TXMLEngine* xml, XMLNodePointer_t root_node);
  bool ReadSplines (TXMLEngine* xml, XMLNodePointer_t root_node);

  // Set BANFF fit info
  void TellMeAboutTheBANFFfit(BANFF * banff);

  //
  void BuildXsBinningMCTemplates(TH2D* input_mctempl_array[kns][knm]);

  //
  TH2D * Template(conv::SuperkSample_t s, conv::SuperkMode_t m) const;

  // Penalty term included in the chisq if NuXSec nuisance terms are included in the fit
  double ChisqPenalty(bool all_syst = true, bool hk_syst = false) const;

  //ask for splines stored - to plot them
  /*
    TSpline3 * XSecSpl(conv::SuperkSample_t s, int idx){ 
    if (s == conv::kSKsFCFV1RingMuLike){
    //if ( *fXSecWghtSpl1Rmu[idx]) 
    return *fXSecWghtSpl1Rmu[idx];
    //else return 0.;
    }
    else{
    //if(*fXSecWghtSpl1Re[idx]) 
    return *fXSecWghtSpl1Re[idx];
      //else return 0.;
      }
      }*/


private:

  void Init    (void);
  void CleanUp (void);

  // Reweight to best-fit value of input BANFF cross-section relared param. Weight is calculated for a single (input) parameter.
  // The calculated weight can not be applied directly to the MC template bins due to binning differences.
  // In FinTuningWght(), the weights for all systematic params are transformed --so that they can be applied 
  // to MC template bins-- and they are combined together.
  double FinTuningWght_1p(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Calculate weight in response to systematic parameter variations. Weight is calculated for a single (input) parameter.
  // Comments similar to the ones above apply
  double SystWght_1p(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Calculate weight spline index
  int  WghtSplineIdx        (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Apply systematic to inmput sample/mode and energy bins?
  bool ApplySyst            (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const; 
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco, double Etrue) const;

  //
  bool IsThereAnEvent       (conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  int  HowManyEvents        (conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  
  //
  // Private data members
  //

  BANFF *     fBANFF;               ///< BANFF fit info
  TSpline3 *  fXSecWghtSpl1Rmu[knwspl1Rmu]; ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  TSpline3 *  fXSecWghtSpl1Re[knwspl1Re]; ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  TSpline3 *  fRHCXSecWghtSpl1Rmu[kRHCnwspl1Rmu]; ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  TSpline3 *  fRHCXSecWghtSpl1Re[kRHCnwspl1Re]; ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  TH2D *      fXsBinningTempl[kns][knm]; ///< Ereco-Etrue templates in XSec binning scheme

  // tweaking dials
  double fFermiMomTwkDial;
  double fSpecFuncTwkDial;
  double fWshapeTwkDial;
  double fShapeCCothTwkDial;
  double fPiLessDeltaTwkDial;
  double fBindETwkDial;
  double fNormCCcohTwkDial;
  double fNormNC1piCTwkDial;
  double fNormNCcohTwkDial;
  double fNormNCothTwkDial;
  double fNormNuEToNuMuTwkDial;
  double fNormAntiNuToNuTwkDial;

  // relative errors
  double fNormCCcohRelErr;
  double fNormNC1piCRelErr;
  double fNormNCcohRelErr;
  double fNormNCothRelErr;
  double fNormNuEToNuMuRelErr;
  double fNormAntiNuToNuRelErr;
  double fSplineMAQERelErr;
  double fSpline1overMAQE2RelErr;
  double fSplineMARESRelErr;
  
  util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject

ClassDef(NuXSec, 4);
};

} // t2kana

#endif
