//____________________________________________________________________________
/*!
   
\class    valor::NuXSec
      
\brief    Neutrino cross-section info used in the oscillation fit
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________
   
#ifndef _valor_nu_xsec_h_
#define _valor_nu_xsec_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TH2.h>
#include <TH3.h>
#include <TSpline.h>
#include <TFile.h>

#include "Binning.h"
#include "Conventions.h"
#include "Enumerations.h"
#include "Utilities.h"

namespace valor {

class BANFF;

class NuXSec: public TObject
{
public:

  NuXSec();
 ~NuXSec();
 
  using TObject::Read; // Suppress clang warnings
  using TObject::Draw; // Suppress clang warnings

  // Reweight to best-fit values of the BANFF cross-section related params
  double FinTuningWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  double FinTuningWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, double Q2) const;

  // Calculate weight in response to systematic parameter variations of cross-section parameters
  // For the BANFF parameters, the weight is relative to the best-fit (tuned MC) parameter, not the nominal MC one
  double SystWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  double SystWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, double Q2) const;

  // Tweak systematic parameters held in this object.
  // Note: BANFF cross-section systematic parameters are tweaked via the BANFF object and they are 
  // correlated with the flux systematic params.

  // NEW XSEC
  void TweakRPA              (double twkdial) { this->Tweak(conv::kNuisXSecRPA, twkdial); } 

  void TweakFermiMom         (double twkdial) { fFermiMomTwkDial         = twkdial; } 
  void TweakSpecFunc         (double twkdial) { fSpecFuncTwkDial         = twkdial; } 
  void TweakWshape           (double twkdial) { fWshapeTwkDial           = twkdial; }
  void TweakShapeCCoth       (double twkdial) { fShapeCCothTwkDial       = twkdial; } 
  void TweakPiLessDelta      (double twkdial) { fPiLessDeltaTwkDial      = twkdial; }
  void TweakBindE            (double twkdial) { fBindETwkDial            = twkdial; } 
  void TweakNormCCcoh        (double twkdial) { fNormCCcohTwkDial        = twkdial; } 
  void TweakNormNC1piC       (double twkdial) { fNormNC1piCTwkDial       = twkdial; }
  void TweakNormNCcoh        (double twkdial) { fNormNCcohTwkDial        = twkdial; } 
  void TweakNormNCoth        (double twkdial) { fNormNCothTwkDial        = twkdial; } 
  void TweakNormNuEToNuMu    (double twkdial) { fNormNuEToNuMuTwkDial    = twkdial; } 
  void TweakNormAntiNuToNu   (double twkdial) { fNormAntiNuToNuTwkDial   = twkdial; }

  // NEW XSEC
  void TweakRHCRPA           (double twkdial) { this->Tweak(conv::kRHCNuisXSecRPA, twkdial); } 

  void TweakRHCFermiMom      (double twkdial) { fRHCFermiMomTwkDial      = twkdial; } 
  void TweakRHCSpecFunc      (double twkdial) { fRHCSpecFuncTwkDial      = twkdial; }
  void TweakRHCWshape        (double twkdial) { fRHCWshapeTwkDial        = twkdial; }
  void TweakRHCShapeCCoth    (double twkdial) { fRHCShapeCCothTwkDial    = twkdial; } 
  void TweakRHCPiLessDelta   (double twkdial) { fRHCPiLessDeltaTwkDial   = twkdial; }
  void TweakRHCBindE         (double twkdial) { fRHCBindETwkDial         = twkdial; } 
  void TweakRHCNormCCcoh     (double twkdial) { fRHCNormCCcohTwkDial     = twkdial; } 
  void TweakRHCNormNC1piC    (double twkdial) { fRHCNormNC1piCTwkDial    = twkdial; }
  void TweakRHCNormNCcoh     (double twkdial) { fRHCNormNCcohTwkDial     = twkdial; } 
  void TweakRHCNormNCoth     (double twkdial) { fRHCNormNCothTwkDial     = twkdial; } 
  void TweakRHCNormNuEToNuMu (double twkdial) { fRHCNormNuEToNuMuTwkDial = twkdial; } 
  void TweakRHCNormAntiNuToNu(double twkdial) { fRHCNormAntiNuToNuTwkDial= twkdial; }
  
  void Tweak(conv::PhysParamId_t p , double twkdial);

  void SyncRHCToFHCSystematics();

  // Get current value of systematic params handled by this class
  double CurrentValue        (conv::PhysParamId_t p) const;

  // Get the relative error of systematic params handled by this class
  double RelativeError (conv::PhysParamId_t p) const {assert(false); return 0;};

  // Get weight spline for each sample / mode / parameters / Ereco bin / Etrue bin
  TSpline3 * FindWghtSpline(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  TSpline3 * FindWghtSpline(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, double Q2) const;

  // Print
  void PrintSystematicVariations (std::ostream & stream = std::cout) const {};
  void PrintRHCSystematicVariations (std::ostream & stream = std::cout) const {};

  // Read inputs from XML file
  bool Read          (TXMLEngine* xml, XMLNodePointer_t root_node);
  bool ReadSplines   (TXMLEngine* xml, XMLNodePointer_t root_node);
  bool ReadSplines2D (TXMLEngine* xml, XMLNodePointer_t root_node);

  // Set BANFF fit info
  void TellMeAboutTheBANFFfit   (BANFF * banff);
  void TellMeAboutTheRHCBANFFfit(BANFF * banff);

  //
  void BuildXsBinningMCTemplates(TH2D* input_mctempl_array[valor::conv::kns][valor::conv::knm]);
  void BuildXsBinningMCTemplates(TH3D* input_mctempl_array[valor::conv::kns][valor::conv::knm]);

  //
  TH2D * Template  (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  TH3D * Template2D(conv::SuperkSample_t s, conv::SuperkMode_t m) const;

  // Penalty term included in the chisq if NuXSec nuisance terms are included in the fit
  double ChisqPenalty(bool all_syst = true, bool hk_syst = false) const {return 0;};

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
    }
  */
      
  // Apply systematic to inmput sample/mode and energy bins?
  bool ApplySyst            (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const; 
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco, double Etrue) const;
  void PrintXSecSystApplication();
 
private:

  // Check that XSecSystConv.h is setup correctly
  void CheckListOfXSecSyst(int nparams, const conv::PhysParamId_t * params, const char * vname, bool is_fatal);
  void CheckIsSplineParam (int nparams, const conv::PhysParamId_t * params, const bool * spline, const char * vname);
  void CheckIsBANFFParam  (int nparams, int nbanffxsecparams, const conv::PhysParamId_t * params, const bool * banff, const char * vname);
  void CheckERange        (int nparams, const conv::PhysParamId_t * params, const double * erange, const char * vname);
  void CheckSKMode        (int nparams, const conv::PhysParamId_t * params, const bool * modes, const int array_order_mode, const char * vname);

  void Init    (void);
  void CleanUp (void);


  // Reweight to best-fit value of input BANFF cross-section relared param. Weight is calculated for a single (input) parameter.
  // The calculated weight can not be applied directly to the MC template bins due to binning differences.
  // In FinTuningWght(), the weights for all systematic params are transformed --so that they can be applied 
  // to MC template bins-- and they are combined together.
  double FinTuningWght_1p(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  double FinTuningWght_1p(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, double Q2) const;

  // Calculate weight in response to systematic parameter variations. Weight is calculated for a single (input) parameter.
  // Comments similar to the ones above apply
  double SystWght_1p(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  double SystWght_1p(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, double Q2) const;

  // Calculate weight spline index
  int  WghtSplineIdx(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  int  WghtSplineIdx(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, double Q2) const;

  //
  bool IsThereAnEvent(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  bool IsThereAnEvent(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, double Q2) const;
  
  int  HowManyEvents(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  int  HowManyEvents(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, double Q2) const;
  
  std::string BuildSplineName(conv::SuperkSample_t s, conv::SuperkMode_t m, conv::PhysParamId_t p, int reco_bin, int true_bin) const;  
  std::string BuildSplineName(conv::SuperkSample_t s, conv::SuperkMode_t m, conv::PhysParamId_t p, int reco_bin, int true_bin, int q2bin) const;
  std::string BuildSplineName(conv::SuperkSample_t s, conv::SuperkMode_t m, conv::PhysParamId_t p) const;
  
  //
  // Private data members
  //

  BANFF *     fBANFF;               ///< BANFF fit info
  TSpline3 *  fXSecWghtSpl1Rmu   [conv::knwspl1Rmu];    ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  TSpline3 *  fXSecWghtSpl1Re    [conv::knwspl1Re];     ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  TSpline3 *  fRHCXSecWghtSpl1Rmu[conv::kRHCnwspl1Rmu]; ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  TSpline3 *  fRHCXSecWghtSpl1Re [conv::kRHCnwspl1Re];  ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin

  /// NEW MULTIRING
  TSpline3 *  fXSecWghtSplMultiRmu[conv::knwsplMultiRmu]; ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  TSpline3 *  fXSecWghtSplMultiRe [conv::knwsplMultiRe];  ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin
  
  ////

  TH2D *      fXsBinningTempl[conv::kns][conv::knm]; ///< Ereco-Etrue templates in XSec binning scheme
  TH3D *      fXsBinningTempl2D[conv::kns][conv::knm]; ///< Ereco-Etrue templates in XSec binning scheme

  // tweaking dials
  double fRPATwkDial; // NEW XSEC
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
  double fRHCRPATwkDial; // NEW XSEC
  double fRHCFermiMomTwkDial;
  double fRHCSpecFuncTwkDial;
  double fRHCWshapeTwkDial;
  double fRHCShapeCCothTwkDial;
  double fRHCPiLessDeltaTwkDial;
  double fRHCBindETwkDial;
  double fRHCNormCCcohTwkDial;
  double fRHCNormNC1piCTwkDial;
  double fRHCNormNCcohTwkDial;
  double fRHCNormNCothTwkDial;
  double fRHCNormNuEToNuMuTwkDial;
  double fRHCNormAntiNuToNuTwkDial;
  double fNormRHCToFHCTwkDial;

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
  double fRHCNormCCcohRelErr;
  double fRHCNormNC1piCRelErr;
  double fRHCNormNCcohRelErr;
  double fRHCNormNCothRelErr;
  double fRHCNormNuEToNuMuRelErr;
  double fRHCNormAntiNuToNuRelErr;
  double fRHCSplineMAQERelErr;
  double fRHCSpline1overMAQE2RelErr;
  double fRHCSplineMARESRelErr;
  double fNormRHCToFHCRelErr;
  
  TVectorD * fUntunedCentralValues;
  TVectorD * fRHCUntunedCentralValues;
  
  util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject

ClassDef(NuXSec, 6);
};

} // valor

#endif
