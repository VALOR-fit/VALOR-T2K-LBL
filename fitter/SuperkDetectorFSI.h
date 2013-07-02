//____________________________________________________________________________
/*!
   
\class    t2kana::FSI
      
\brief    Final state interaction (FSI) info used in the oscillation fit
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 21, 2012

*/
//____________________________________________________________________________
   
#ifndef _t2kana_skdetfsi_h_
#define _t2kana_skdetfsi_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>

#include "Conventions.h"

using t2kana::conv::kns;
using t2kana::conv::knm;

namespace t2kana {

class SuperkDetectorFSI: public TObject
{
public:

  SuperkDetectorFSI();
 ~SuperkDetectorFSI();

  // Calculate weight in response to systematic parameter variations of FSI parameters
  double SystWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  
  // Info about the energy scale error and current energy scale tweaking factor              
  double  EnergyScaleRelErr   (conv::SuperkSample_t s, conv::SuperkMode_t m) const; ///< energy scale err 
  double  EnergyScale         (conv::SuperkSample_t s, conv::SuperkMode_t m) const; ///< energy scale 
  double  EnergyScale         (conv::SuperkSample_t s) const; ///< energy scale              

  // Tweak FSI params (either force the input params or generate correlated params using the Cholesky method)
  void TweakEnergyScale                    (double twkdial);
  void Tweak (bool force_inputs = false, TVectorD * input_params = 0);
  void TweakOnly1Rmu (bool force_inputs = false, TVectorD * input_params = 0);
  void TweakOnly1Re (bool force_inputs = false, TVectorD * input_params = 0);

  // Get current value of systematic params handled by this class
  double CurrentValue (conv::PhysParamId_t p) const;
  double OneSigmaError(conv::PhysParamId_t p) const; //square root of cov matrix

  // Get input FSI covariance matrix and its Cholesky decomposition
  const TMatrixD & CovMatrix         (void) const { return *fCovMtrx;       }
  const TMatrixD & CovMatrixCholesky (void) const { return *fCovMtrxChoL;   }
  const TMatrixD & CovMatrixInv      (void) const { return *fCovMtrxInv;    }

  // Penalty term included in the chisq if FSI nuisance terms are included in the fit       
  double ChisqPenalty(void) const;

  // Print
  void PrintSystematicVariations (ostream & stream = std::cout) const;

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
   
private:

  void Init      (void);
  void CleanUp   (void);
  // Apply the given FSI parameter to the given sample, mode, Ereco?
  bool ApplySyst            (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const;
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco) const;
  //
  // Private data members
  //

  TVectorD * fCurrentValues;    ///< Current values of FSI params
  TVectorD * fNominalValues;    ///< Nominal values of FSI params
  TMatrixD * fCovMtrx;          ///< FSI covariance matrix
  TMatrixD * fCovMtrxInv;       ///< Inverse FSI covariance matrix
  TMatrixD * fCovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix

  double     fEnergyScale1RmuRelErr;       ///< 1Rmu energy scale relative error   
  double     fEnergyScale1ReRelErr;        ///< 1Re energy scale relative error    
  // > 1Rmu + 1Re 
  double     fEnergyScaleTwkDial;          ///< 1Rmu/1Re Energy scale tweaking dial         

ClassDef(SuperkDetectorFSI, 1);
};

} // t2kana

#endif
