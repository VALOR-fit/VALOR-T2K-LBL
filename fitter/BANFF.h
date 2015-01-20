
//____________________________________________________________________________
/*!

\class    valor::BANFF

\brief    Class to hold analysis inputs from BANFF (Beam and ND280 Flux
          Extrapolation Task Force)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  February 09, 2012

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________

#ifndef _valor_banff_h_
#define _valor_banff_h_

#include <TObject.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMatrixTSym.h>
#include <TXMLEngine.h>

#include <iostream>

#include <Enumerations.h>

class TH1D;
class TH2D;
class TVirtualPad;

namespace valor {

class BANFF: public TObject
{
public:
  BANFF();   
 ~BANFF();
 
  using TObject::Read; // Suppress clang warnings
  using TObject::Draw; // Suppress clang warnings

  // Get BANFF fit param values (nominal + tweaked)
  double BestFitValue (conv::PhysParamId_t p) const;
  double BestFitValue (int idx) const;
  double RHCBestFitValue (conv::PhysParamId_t p) const;
  double RHCBestFitValue (int idx) const;

  // NEW FLUX FHC+RHC
  double FHCpRHCBestFitValue (conv::PhysParamId_t p) const;
  double FHCpRHCBestFitValue (int idx) const;
  //
  
  double CurrentValue (conv::PhysParamId_t p) const;
  double CurrentValue (int idx) const;
  double RHCCurrentValue (conv::PhysParamId_t p) const;
  double RHCCurrentValue (int idx) const;

  // NEW FLUX FHC+RHC
  double FHCpRHCCurrentValue (conv::PhysParamId_t p) const;
  double FHCpRHCCurrentValue (int idx) const;
  //

  double OneSigmaError(conv::PhysParamId_t p) const;
  double OneSigmaError(int idx) const;
  double RHCOneSigmaError(conv::PhysParamId_t p) const;
  double RHCOneSigmaError(int idx) const;

  // NEW FLUX FHC+RHC
  double FHCpRHCOneSigmaError(conv::PhysParamId_t p) const;
  double FHCpRHCOneSigmaError(int idx) const;
  //




  // Get BANFF best fit values, covariance matrix and cholesky decomposition of covariance matrix
  const TVectorD & BestFitValues              (void) const { return *fBestFitValues;    }
  const TMatrixD & CovMatrix                  (void) const { return *fCovMtrx;          }
  const TMatrixD & CovMatrixInv               (void) const { return *fCovMtrxInv;       }
  const TMatrixD & CovMatrixCholesky          (void) const { return *fCovMtrxChoL;      }
  const TVectorD & RHCBestFitValues           (void) const { return *fRHCBestFitValues; }
  const TMatrixD & RHCCovMatrix               (void) const { return *fRHCCovMtrx;       }
  const TMatrixD & RHCCovMatrixInv            (void) const { return *fRHCCovMtrxInv;    }

  // NEW FLUX FHC+RHC
  const TVectorD & FHCpRHCBestFitValues           (void) const { return *fFHCpRHCBestFitValues; }
  const TMatrixD & FHCpRHCCovMatrix               (void) const { return *fFHCpRHCCovMtrx;       }
  const TMatrixD & FHCpRHCCovMatrixInv            (void) const { return *fFHCpRHCCovMtrxInv;    }
  //

  // Systematic variations (wrt the best-fit BANFF params)
  void Tweak             (bool force_inputs = false, TVectorD * input_vec = 0); 
  void TweakRHC          (bool force_inputs = false, TVectorD * input_vec = 0);  
  void TweakFHCpRHC      (bool force_inputs = false, TVectorD * input_vec = 0); // NEW FLUX FHC+RHC
  void TweakSample       (bool force_inputs = false, TVectorD * input_vec = 0);   

  void TweakFluxOnly        (bool force_inputs = false, TVectorD * input_vec = 0);
  void TweakRHCFluxOnly     (bool force_inputs = false, TVectorD * input_vec = 0); 
  void TweakFHCpRHCFluxOnly (bool force_inputs = false, TVectorD * input_vec = 0); // NEW FLUX FHC+RHC
  void TweakSampleFluxOnly  (bool force_inputs = false, TVectorD * input_vec = 0);

  void TweakXSecOnly        (bool force_inputs = false, TVectorD * input_vec = 0);
  void TweakRHCXSecOnly     (bool force_inputs = false, TVectorD * input_vec = 0); 
  void TweakFHCpRHCXSecOnly (bool force_inputs = false, TVectorD * input_vec = 0); // NEW FLUX FHC+RHC
  void TweakSampleXSecOnly  (bool force_inputs = false, TVectorD * input_vec = 0);
    

  void SyncRHCToFHCSystematics();

  // A special systematic variation which cancels the effect of BANFF tuning and reverts back to the nominal MC
  void RestoreNominalMC(void);
  void RestoreRHCNominalMC(void);
  void RestoreFHCpRHCNominalMC(void);
 
  void RestoreTunedMC(void);
  void RestoreRHCTunedMC(void);
  void RestoreFHCpRHCTunedMC(void);


  // Penalty term included in the chisq if the BANFF fit parameters are included in the oscillation fit.
  double ChisqPenalty             (bool all_syst = true) const;
  double RHCChisqPenalty          (bool all_syst = true) const;
  
  // NEW FLUX FHC+RHC
  double FHCpRHCChisqPenalty          (bool all_syst = true) const;
  //

  // Get & tweak flags	
  bool IsTweaked              (void) const { return fIsTweaked;            } ///< BANFF params tweaked?
  bool IsTweakedRHC           (void) const { return fIsTweakedRHC;         } ///< BANFF params tweaked?

  // NEW FLUX FHC+RHC
  bool IsTweakedFHCpRHC           (void) const { return fIsTweakedFHCpRHC;         } ///< BANFF params tweaked?
  //

  bool IsFluxWghtUpdated      (void) const { return fIsFluxWghtUpdated;    } ///< NuFlux object updated flux in response?
  bool IsFluxWghtUpdatedRHC   (void) const { return fIsFluxWghtUpdatedRHC; } ///< NuFlux object updated flux in response?

  // NEW FLUX FHC+RHC
  bool IsFluxWghtUpdatedFHCpRHC   (void) const { return fIsFluxWghtUpdatedFHCpRHC; } ///< NuFlux object updated flux in response?
  //

  void MarkFluxWghtUpdated    (void)       { fIsFluxWghtUpdated = true;    }
  void MarkFluxWghtUpdatedRHC (void)       { fIsFluxWghtUpdatedRHC = true; }

  // NEW FLUX 
  void MarkFluxWghtUpdatedFHCpRHC (void)       { fIsFluxWghtUpdatedFHCpRHC = true; }
  //

  // Print
  void PrintSystematicVariations (std::ostream & stream = std::cout) const;
  void PrintRHCSystematicVariations (std::ostream & stream = std::cout) const;
  
  // NEW FHC+RHC
  void PrintFHCpRHCSystematicVariations (std::ostream & stream = std::cout) const;
  //

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
  bool ReadRHC (TXMLEngine* xml, XMLNodePointer_t root_node);
  
  // NEW FLUX FHC+RHC
  bool ReadFHCpRHC (TXMLEngine* xml, XMLNodePointer_t root_node);
  //

private:

  void Init      (void);
  void CleanUp   (void);
  int  Param2Idx   (conv::PhysParamId_t p) const; ///< convert BANFF PhysParamId_t to a covariance matrix or vector index
  void CholeskyDecomposeMatrix (void);    ///< Cholesky decompose the covariance matrix
  void CholeskyDecomposeRHCMatrix (void);    ///< Cholesky decompose the covariance matrix
  
  // NEW FLUX FHC+RHC
  void CholeskyDecomposeFHCpRHCMatrix (void);    ///< Cholesky decompose the covariance matrix
  //

  TVectorD *  fBestFitValues;     ///< BANFF best fit values
  TVectorD *  fCurrentValues;     ///< Current values of BANFF parameters 
  TMatrixD *  fCovMtrx;           ///< BANFF covariance matrix
  TMatrixD *  fCovMtrxInv;        ///< Inverse BANFF covariance matrix
  TMatrixD *  fCovMtrxChoL;       ///< Triangular matrix obtained from Cholesky decomposition of BANFF covariance matrix
  bool        fIsTweaked;         ///< Were BANFF params tweaked?
  bool        fIsFluxWghtUpdated; ///< Were Flux weights updated as a result of the last BANFF param variation

  TVectorD *  fRHCBestFitValues;     ///< BANFF best fit values
  TVectorD *  fRHCCurrentValues;     ///< Current values of BANFF parameters 
  TMatrixD *  fRHCCovMtrx;           ///< BANFF covariance matrix
  TMatrixD *  fRHCCovMtrxInv;        ///< Inverse BANFF covariance matrix
  TMatrixD *  fRHCCovMtrxChoL;       ///< Triangular matrix obtained from Cholesky decomposition of BANFF covariance matrix
  bool        fIsTweakedRHC;         ///< Were BANFF params tweaked?
  bool        fIsFluxWghtUpdatedRHC; ///< Were Flux weights updated as a result of the last BANFF param variation

  // NEW FLUX FHC+RHC
  TVectorD *  fFHCpRHCBestFitValues;     ///< BANFF best fit values
  TVectorD *  fFHCpRHCCurrentValues;     ///< Current values of BANFF parameters 
  TMatrixD *  fFHCpRHCCovMtrx;           ///< BANFF covariance matrix
  TMatrixD *  fFHCpRHCCovMtrxInv;        ///< Inverse BANFF covariance matrix
  TMatrixD *  fFHCpRHCCovMtrxChoL;       ///< Triangular matrix obtained from Cholesky decomposition of BANFF covariance matrix
  bool        fIsTweakedFHCpRHC;         ///< Were BANFF params tweaked?
  bool        fIsFluxWghtUpdatedFHCpRHC; ///< Were Flux weights updated as a result of the last BANFF param variation
  //

  double * fBestFitArray; //! ///< Pointer to the actual memory location of the best fit values
  double * fCurrentArray; //! ///< Pointer to the actual memory location of the current values
  
  double * fRHCBestFitArray; //! ///< Pointer to the actual memory location of the best fit values
  double * fRHCCurrentArray; //! ///< Pointer to the actual memory location of the current values

  // NEW FLUX FHC+RHC
  double * fFHCpRHCBestFitArray; //! ///< Pointer to the actual memory location of the best fit values
  double * fFHCpRHCCurrentArray; //! ///< Pointer to the actual memory location of the current values
  //

  int fParam2Idx[conv::knppmax];
  
  TVectorD * fUntunedCentralValues;    // Store the untuned values of the params (for RestoreNominalMC methods)
  TVectorD * fRHCUntunedCentralValues; // Store the untuned values of the params (for RestoreNominalMC methods)

ClassDef(BANFF, 1)
};

} //valor

#endif
