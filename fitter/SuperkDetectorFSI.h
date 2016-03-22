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
  
class FitParamList;
class PhysParamList;

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

  void Tweak (bool force_inputs = false, TVectorD * input_params = 0, Option_t * opt = "");
  
  // Save tweaks in a VALOR fit object
  void CopyTweaksTo(valor::PhysParamList & ppl) const;
  void CopyTweaksTo(valor::FitParamList & fpl) const;
  
  // Get tweaks from a VALOR fit object
  void UpdateFrom(const valor::FitParamList  & fpl);
  
  // Or just set one
  void SetCurrentValue(conv::PhysParamId_t p, double tweak);
  
  // Get current value of systematic params handled by this class
  double CurrentValue (conv::PhysParamId_t p) const;
  double OneSigmaError(conv::PhysParamId_t p) const; //square root of cov matrix

  // Get input FSI covariance matrix and its Cholesky decomposition
  const TMatrixD & CovMatrix         (void) const { return *fCovMtrx;       }
  const TMatrixD & CovMatrixCholesky (void) const { return *fCovMtrxChoL;   }
  const TMatrixD & CovMatrixInv      (void) const { return *fCovMtrxInv;    }
  ///

  // Penalty term included in the chisq if FSI nuisance terms are included in the fit       
  double ChisqPenalty(void) const;
  ///

  // Print
  void PrintSystematicVariations    (std::ostream & stream = std::cout, Option_t * opt = "") const;
  ///

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
  
  // Apply the given FSI parameter to the given sample, mode, Ereco?
  bool ApplySyst            (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const;
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco) const;

private:

  // Check that SKFSISystConv.h is setup correctly
  void CheckListOfSKDetFSISyst(int nparams, const conv::PhysParamId_t * params, const char * vname, bool allow_gap);
  void CheckERange            (int nparams, const conv::PhysParamId_t * params, const double * erange, const char * vname);
  void CheckSKMode            (int nparams, const conv::PhysParamId_t * params, const bool * modes, const char * vname, bool fatal_no_effect);

  void Init      (void);
  void CleanUp   (void);
  //
  // Private data members
  //
  
  int fNumOfActiveSKDetFSIParams;
  conv::PhysParamId_t * fActiveSKDetFSIParams; //![fNumOfActiveSKDetFSIParams] - the parameters in use. 
  
  double     fEnergyScaleRelErrMu; ///< energy scale relative error   
  double     fEnergyScaleRelErrE;  ///< energy scale relative error   
  double     fEnergyScaleTwkDial;  ///< Energy scale tweaking dial         
  
  TVectorD * fCurrentValues;    ///< Current values of FSI params
  TVectorD * fNominalValues;    ///< Nominal values of FSI params
  TMatrixD * fCovMtrx;          ///< FSI covariance matrix
  TMatrixD * fCovMtrxInv;       ///< Inverse FSI covariance matrix
  TMatrixD * fCovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix
  ///

  util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject
  
  
  
  
  
  
  
  //____________________________________________________________________
  // DEPRECATED INTERACE, DO NOT USE
 public:

  const TMatrixD & RHCCovMatrix         (void) const { return *fCovMtrx;       }
  const TMatrixD & RHCCovMatrixCholesky (void) const { return *fCovMtrxChoL;   }
  const TMatrixD & RHCCovMatrixInv      (void) const { return *fCovMtrxInv;    }
  // NEW FHC+RHC JOINT
  const TMatrixD & FHCpRHCCovMatrix         (void) const { return *fCovMtrx;       }
  const TMatrixD & FHCpRHCCovMatrixCholesky (void) const { return *fCovMtrxChoL;   }
  const TMatrixD & FHCpRHCCovMatrixInv      (void) const { return *fCovMtrxInv;    }
  /// NEW MULTIRING 2
  const TMatrixD & MultiRingCovMatrix         (void) const { return *fCovMtrx;       }
  const TMatrixD & MultiRingCovMatrixCholesky (void) const { return *fCovMtrxChoL;   }
  const TMatrixD & MultiRingCovMatrixInv      (void) const { return *fCovMtrxInv;    }
  ///
  
  void PrintRHCSystematicVariations (std::ostream & stream = std::cout) const{stream<<"CALLING DEPRECATED FUNCTION SuperkDetectorFSI/PrintRHCSystematicVariations"<<std::endl;}
  // NEW FHC+RHC JOINT
  void PrintFHCpRHCSystematicVariations (std::ostream & stream = std::cout) const{stream<<"CALLING DEPRECATED FUNCTION SuperkDetectorFSI/PrintFHCpRHCSystematicVariations"<<std::endl;}
  //
  /// NEW MULTIRING 2
  void PrintMultiRingSystematicVariations (std::ostream & stream = std::cout) const{stream<<"CALLING DEPRECATED FUNCTION SuperkDetectorFSI/PrintMultiRingSystematicVariations"<<std::endl;}
  ///
  //
  /// NEW MULTIRING 2
  double MultiRingChisqPenalty(void) const {return 9E9;}
  
  //~ bool ReadRHC (TXMLEngine* , XMLNodePointer_t ) {return 1;}
  //~ bool ReadFHCpRHC (TXMLEngine* , XMLNodePointer_t ) {return 1;}
  //~ /// NEW MULTIRING 2
  //~ bool ReadMultiRing (TXMLEngine* , XMLNodePointer_t ) {return 1;}
  ///
  //____________________________________________________________________
  
ClassDef(SuperkDetectorFSI, 1);
};

} // valor

#endif
