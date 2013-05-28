//____________________________________________________________________________
/*!

\class    t2kana::SuperkDetector

\brief    Holds info for SuperK detector systematics (E scale, efficiencies)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 15, 2011

*/
//____________________________________________________________________________

#ifndef _t2kana_sk_mc_detector_h_
#define _t2kana_sk_mc_detector_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TMatrixD.h>

#include "Conventions.h"
#include "PhysHypothesis.h"
#include "NuOsc3Flavour.h"

class TVirtualPad;
class TTree;
class TH1D;
class TH2D;

using t2kana::conv::kns;
using t2kana::conv::knm;

namespace t2kana {

class SuperkDetector: public TObject
{
public:
  SuperkDetector();
 ~SuperkDetector();

  // Tweak systematic parameters
  void TweakEnergyScale                    (double twkdial);
  void TweakEfficiency1RmuCCQEOther        (double twkdial); ///< "other": all but ring-counting
  void TweakEfficiency1RmuCCQERingCounting (bool force_input_twkdials=false, double twkdial1=0, double twkdial2=0, double twkdial3=0, bool correlated=true); 
  void TweakEfficiency1RmuCCnonQE          (double twkdial);
  void TweakEfficiency1RmuNC               (double twkdial);
  void TweakEfficiency1RmuNuECC            (double twkdial);
  // Tweak systematic parameters for 1Re - they go in a single covariance matrix
  void TweakEfficiencies1Re                (bool force_inputs = false, TVectorD * input_params = 0);

  // Calculate weight due to efficiency tweaks
  double  EfficiencyWeight (conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const;

  // Info about the energy scale error and current energy scale tweaking factor
  double  EnergyScaleRelErr   (conv::SuperkSample_t s, conv::SuperkMode_t m) const; ///< energy scale err
  double  EnergyScale         (conv::SuperkSample_t s, conv::SuperkMode_t m) const; ///< energy scale
  double  EnergyScale         (conv::SuperkSample_t s) const; ///< energy scale

  // Penalty term included in the chisq if SuperK detector nuisance terms are included in the fit
  double ChisqPenalty(bool all_syst = true) const;

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
      
  // Print systematic parameters
  void PrintSystematicVariations (ostream & stream = std::cout) const;

  // Get current value of systematic params handled by this class
  double CurrentValue (conv::PhysParamId_t p) const;

  // Access the elements of the ring counting matrix
  const TMatrixD & Eff1RmuCCQERCCovarRelErr (void) const { return *fEff1RmuCCQERCCovarRelErr; }
  const TMatrixD & Eff1RmuCCQERCInvCovarRelErr (void) const { return *fEff1RmuCCQERCInvCovarRelErr; }

  // Access the elements of the efficiency matrix for 1Re
  const TMatrixD & Eff1ReCovarRelErr (void) const { return *fEff1ReCovarRelErr; }
  const TMatrixD & Eff1ReInvCovarRelErr (void) const { return *fEff1ReInvCovarRelErr; }

private:

  // Init and clean-up
  void Init    (void);
  void CleanUp (void);

  // Apply the given SK 1Re eff parameter to the given sample, mode, Ereco?                                                                                                                              
#if FSI_BINS == 1
  bool ApplySyst1Re         (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const;
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco) const;
#endif

  //
  // Private data members
  //

  // systematic uncertainties: 1 sigma fractional errors & covariance matrices
  // > 1Rmu
  double     fEnergyScale1RmuRelErr;       ///< 1Rmu energy scale relative error
  double     fEnergyScale1ReRelErr;        ///< 1Re energy scale relative error
  double     fEff1RmuCCQEOthRelErr;        ///< 1Rmu CCQE other (all but ring-counting) efficiency relative error 
  TMatrixD * fEff1RmuCCQERCCovarRelErr;    ///< 1Rmu CCQE ring-counting efficiency relative error covariance matrix (3x3)
  TMatrixD * fEff1RmuCCQERCInvCovarRelErr; ///< 1Rmu CCQE ring-counting efficiency inverse of relative error covariance matrix (3x3)
  double     fEff1RmuCCQERCCovarE1RelErr;  ///<  -  ~0.4 GeV (above RC cov matrix is 3x3, 3 Ereco bins: (0, fEff1RmuCCQERCCovarE1RelErr, fEff1RmuCCQERCCovarE2RelErr, oo))
  double     fEff1RmuCCQERCCovarE2RelErr;  ///<  -  ~1.1 GeV (above RC cov matrix is 3x3, 3 Ereco bins: (0, fEff1RmuCCQERCCovarE1RelErr, fEff1RmuCCQERCCovarE2RelErr, oo))
  TMatrixD * fEff1RmuCCQERCCovarChoLRelErr;///<  -  Cholesky decomposition (trangular matrix L) of fEff1RmuCCQERCCovarRelErr covariance matrix
  double     fEff1RmuCCnonQERelErr;        ///< 1Rmu CCnonQE efficiency relative error
  double     fEff1RmuNCRelErr;             ///< 1Rmu NC efficiency relative error
  double     fEff1RmuNuECCRelErr;          ///< 1Rmu nue CC efficiency relative error
  // > 1Re
  TMatrixD * fEff1ReCovarRelErr;    ///< 1Re efficiencies relative error covariance matrix (12x12)
  TMatrixD * fEff1ReInvCovarRelErr; ///< 1Re efficiencies inverse of relative error covariance matrix (12x12)
  TMatrixD * fEff1ReCovarChoLRelErr;///<  -  Cholesky decomposition (trangular matrix L) of fEff1ReCovarRelErr covariance matrix

  // systematic uncertainties: tweaking dials
  // > 1Rmu + 1Re
  double     fEnergyScaleTwkDial;          ///< 1Rmu/1Re Energy scale tweaking dial
  // > 1Rmu 
  double     fEff1RmuCCQEOthTwkDial;       ///< 1Rmu CCQE other (all-but-ring-counting) efficiency tweaking dial
  double     fEff1RmuCCQERC0TwkDial;       ///< 1Rmu CCQE ring-counting efficiency tweaking dial 0 (Ereco < fEff1RmuCCQERCCovarE1RelErr)
  double     fEff1RmuCCQERC1TwkDial;       ///< 1Rmu CCQE ring-counting efficiency tweaking dial 1 (fEff1RmuCCQERCCovarE1RelErr <= Ereco <= fEff1RmuCCQERCCovarE2RelErr)
  double     fEff1RmuCCQERC2TwkDial;       ///< 1Rmu CCQE ring-counting efficiency tweaking dial 2 (Ereco > fEff1RmuCCQERCCovarE2RelErr)
  double     fEff1RmuCCnonQETwkDial;       ///< 1Rmu CCnonQE efficiency tweaking dial
  double     fEff1RmuNCTwkDial;            ///< 1Rmu NC efficiency tweaking dial
  double     fEff1RmuNuECCTwkDial;         ///< 1Rmu nue CC efficiency tweaking dial
  // > 1Re
  TVectorD * f1ReCurrentValues;    ///< Current values of 1Re efficiencies
  TVectorD * f1ReNominalValues;    ///< Nominal values of 1Re efficiencies  

ClassDef(SuperkDetector, 1);
};

} //t2kana

#endif
