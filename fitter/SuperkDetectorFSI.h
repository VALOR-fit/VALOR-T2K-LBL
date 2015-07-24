//____________________________________________________________________________
/*!
   
\class    valor::FSI
      
\brief    Final state interaction (FSI) info used in the oscillation fit
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 21, 2012

*/
//____________________________________________________________________________
   
#ifndef _valor_skdetfsi_h_
#define _valor_skdetfsi_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>

#include "Enumerations.h"
#include "Conventions.h"
#include "Utilities.h"

namespace valor {

class SuperkDetectorFSI: public TObject
{
public:

  SuperkDetectorFSI();
 ~SuperkDetectorFSI();
 
  using TObject::Read; // Suppress clang warnings
  using TObject::Draw; // Suppress clang warnings

  // Calculate weight in response to systematic parameter variations of FSI parameters
  double SystWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  
  // Info about the energy scale error and current energy scale tweaking factor              
  double  EnergyScaleRelErr   (conv::SuperkSample_t s, conv::SuperkMode_t m) const; ///< energy scale err 
  double  EnergyScaleRelErr   (conv::SuperkSample_t s) const; ///< energy scale err 
  double  EnergyScale         (conv::SuperkSample_t s, conv::SuperkMode_t m) const; ///< energy scale 
  double  EnergyScale         (conv::SuperkSample_t s) const; ///< energy scale              

  // Tweak FSI params (either force the input params or generate correlated params using the Cholesky method)
  void TweakEnergyScale                    (double twkdial);
  void TweakRHCEnergyScale                 (double twkdial);
  // NEW FHC+RHC JOINT
  void TweakFHCpRHCEnergyScale             (double twkdial);
  //
  void TweakMultiRingEnergyScale           (double twkdial); /// NEW MULTIRING 3
  void TweakSampleEnergyScale              (double twkdial); 

  void Tweak (bool force_inputs = false, TVectorD * input_params = 0);
  void TweakRHC (bool force_inputs = false, TVectorD * input_params = 0);
  // NEW FHC+RHC JOINT
  void TweakFHCpRHC (bool force_inputs = false, TVectorD * input_params = 0); 
  //
  void TweakMultiRing (bool force_inputs = false, TVectorD * input_params = 0);  /// NEW MULTIRING 2
  void TweakSample (bool force_inputs = false, TVectorD * input_params = 0); 
  
  void TweakOnly1Rmu (bool force_inputs = false, TVectorD * input_params = 0);
  void TweakOnlyRHC1Rmu (bool force_inputs = false, TVectorD * input_params = 0);
  void TweakOnlyMultiRmu (bool force_inputs = false, TVectorD * input_params = 0); /// NEW MULTIRING 2
  
  void TweakOnly1Re (bool force_inputs = false, TVectorD * input_params = 0);
  void TweakOnlyRHC1Re (bool force_inputs = false, TVectorD * input_params = 0);
  void TweakOnlyMultiRe (bool force_inputs = false, TVectorD * input_params = 0); /// NEW MULTIRING 2
 
  // NEW FHC+RHC JOINT ESCALE
  void SyncEscaleRHCToFHCSystematics(); // synchronize only Escale
  //
  void SyncRHCToFHCSystematics();
  void SyncMultiRingToFHCSystematics(); /// NEW MULTIRING 2

  // Get current value of systematic params handled by this class
  double CurrentValue (conv::PhysParamId_t p) const;
  double RHCCurrentValue (conv::PhysParamId_t p) const;
  // NEW FHC+RHC JOINT
  double FHCpRHCCurrentValue (conv::PhysParamId_t p) const;
  //
  double MultiRingCurrentValue (conv::PhysParamId_t p) const; /// NEW MULTIRING 2
  double SampleCurrentValue(conv::PhysParamId_t p) const; // Picks the right one of the above methods
  
  double OneSigmaError(conv::PhysParamId_t p) const; //square root of cov matrix
  double RHCOneSigmaError(conv::PhysParamId_t p) const; //square root of cov matrix
  // NEW FHC+RHC JOINT
  double FHCpRHCOneSigmaError(conv::PhysParamId_t p) const; //square root of cov matrix
  //
  double MultiRingOneSigmaError(conv::PhysParamId_t p) const; //square root of cov matrix   /// NEW MULTIRING 2
  double SampleOneSigmaError(conv::PhysParamId_t p) const; // Picks the right one of the above methods

  // Get input FSI covariance matrix and its Cholesky decomposition
  const TMatrixD & CovMatrix         (void) const { return *fCovMtrx;       }
  const TMatrixD & CovMatrixCholesky (void) const { return *fCovMtrxChoL;   }
  const TMatrixD & CovMatrixInv      (void) const { return *fCovMtrxInv;    }
  const TMatrixD & RHCCovMatrix         (void) const { return *fRHCCovMtrx; }
  const TMatrixD & RHCCovMatrixCholesky (void) const { return *fRHCCovMtrxChoL; }
  const TMatrixD & RHCCovMatrixInv      (void) const { return *fRHCCovMtrxInv; }
  // NEW FHC+RHC JOINT
  const TMatrixD & FHCpRHCCovMatrix         (void) const { return *fFHCpRHCCovMtrx; }
  const TMatrixD & FHCpRHCCovMatrixCholesky (void) const { return *fFHCpRHCCovMtrxChoL; }
  const TMatrixD & FHCpRHCCovMatrixInv      (void) const { return *fFHCpRHCCovMtrxInv; }
  //
  /// NEW MULTIRING 2
  const TMatrixD & MultiRingCovMatrix         (void) const { return *fMultiRingCovMtrx; }
  const TMatrixD & MultiRingCovMatrixCholesky (void) const { return *fMultiRingCovMtrxChoL; }
  const TMatrixD & MultiRingCovMatrixInv      (void) const { return *fMultiRingCovMtrxInv; }
  ///

  // Penalty term included in the chisq if FSI nuisance terms are included in the fit       
  double ChisqPenalty(void) const;
  double RHCChisqPenalty(void) const;
  // NEW FHC+RHC JOINT
  double FHCpRHCChisqPenalty(void); //const;
  //
  /// NEW MULTIRING 2
  double MultiRingChisqPenalty(void) const;
  ///

  // Print
  void PrintSystematicVariations (std::ostream & stream = std::cout) const;
  void PrintRHCSystematicVariations (std::ostream & stream = std::cout) const;
  // NEW FHC+RHC JOINT
  void PrintFHCpRHCSystematicVariations (std::ostream & stream = std::cout) const;
  //
  /// NEW MULTIRING 2
  void PrintMultiRingSystematicVariations (std::ostream & stream = std::cout) const;
  ///

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
  bool ReadRHC (TXMLEngine* xml, XMLNodePointer_t root_node);
  bool ReadFHCpRHC (TXMLEngine* xml, XMLNodePointer_t root_node);  // NEW FHC+RHC JOINT
  /// NEW MULTIRING 2
  bool ReadMultiRing (TXMLEngine* xml, XMLNodePointer_t root_node);
  ///
  
  
  // Apply the given FSI parameter to the given sample, mode, Ereco?
  bool ApplySyst            (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const;
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco) const;

  // NEW FHC+RHC JOINT
  void SplitFHCpRHCvectors(TVectorD bestfit_values_fhcprhc,TVectorD &bestfit_values_fhc,TVectorD &bestfit_values_rhc);
  void MergeFHCpRHCvectors(TVectorD &bestfit_values_fhcprhc,TVectorD bestfit_values_fhc,TVectorD bestfit_values_rhc);
  void UpdateFHCpRHCCurrentValues();
  void UpdateFHCpRHCBestFitValues();
  //

private:

  void Init      (void);
  void CleanUp   (void);
  //
  // Private data members
  //

  TVectorD * fCurrentValues;    ///< Current values of FSI params
  TVectorD * fNominalValues;    ///< Nominal values of FSI params
  TMatrixD * fCovMtrx;          ///< FSI covariance matrix
  TMatrixD * fCovMtrxInv;       ///< Inverse FSI covariance matrix
  TMatrixD * fCovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix

  TVectorD * fRHCCurrentValues;    ///< Current values of FSI params
  TVectorD * fRHCNominalValues;    ///< Nominal values of FSI params
  TMatrixD * fRHCCovMtrx;          ///< FSI covariance matrix
  TMatrixD * fRHCCovMtrxInv;       ///< Inverse FSI covariance matrix
  TMatrixD * fRHCCovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix
  

  // NEW FHC+RHC JOINT
  TVectorD * fFHCpRHCCurrentValues;    ///< Current values of FSI params
  TVectorD * fFHCpRHCNominalValues;    ///< Nominal values of FSI params
  TMatrixD * fFHCpRHCCovMtrx;          ///< FSI covariance matrix
  TMatrixD * fFHCpRHCCovMtrxInv;       ///< Inverse FSI covariance matrix
  TMatrixD * fFHCpRHCCovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix
  //

  /// NEW MULTIRING 2
  TVectorD * fMultiRingCurrentValues;    ///< Current values of FSI params
  TVectorD * fMultiRingNominalValues;    ///< Nominal values of FSI params
  TMatrixD * fMultiRingCovMtrx;          ///< FSI covariance matrix
  TMatrixD * fMultiRingCovMtrxInv;       ///< Inverse FSI covariance matrix
  TMatrixD * fMultiRingCovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix
  ///

  double     fEnergyScale1RmuRelErr;       ///< 1Rmu energy scale relative error   
  double     fEnergyScale1ReRelErr;        ///< 1Re energy scale relative error
  double     fRHCEnergyScale1RmuRelErr;    ///< RHC 1Rmu energy scale relative error   
  double     fRHCEnergyScale1ReRelErr;     ///< RHC 1Re energy scale relative error    
  /// NEW MULTIRING 3
  double     fEnergyScaleMultiRmuRelErr;       ///< MultiRmu energy scale relative error   
  double     fEnergyScaleMultiReRelErr;        ///< MultiRe energy scale relative error
  ///
  // > 1Rmu + 1Re 
  double     fEnergyScaleTwkDial;          ///< 1Rmu/1Re Energy scale tweaking dial         
  double     fRHCEnergyScaleTwkDial;          ///< 1Rmu/1Re Energy scale tweaking dial
  /// NEW MULTIRING 3
  double     fMultiRingEnergyScaleTwkDial;          ///< MultiRmu/MultiRe Energy scale tweaking dial
  ///

  util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject

ClassDef(SuperkDetectorFSI, 1);
};

} // valor

#endif
