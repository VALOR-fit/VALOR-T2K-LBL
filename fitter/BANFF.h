//____________________________________________________________________________
/*!

\class    t2kana::BANFF

\brief    Class to hold analysis inputs from BANFF (Beam and ND280 Flux
          Extrapolation Task Force)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  February 09, 2012

*/
//____________________________________________________________________________

#ifndef _t2kana_banff_h_
#define _t2kana_banff_h_

#include <TObject.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMatrixTSym.h>
#include <TXMLEngine.h>

#include <iostream>

#include "Conventions.h"

class TH1D;
class TH2D;
class TVirtualPad;

namespace t2kana {

class BANFF: public TObject
{
public:
  BANFF();   
 ~BANFF();

  // Get BANFF fit param values (nominal + tweaked)
  double BestFitValue (conv::PhysParamId_t p) const;
  double CurrentValue (conv::PhysParamId_t p) const;
  double CurrentUncorrelatedValue (conv::PhysParamId_t p) const;
  double CurrentEigenDecompValue  (conv::PhysParamId_t p) const;
  double BestFitEigenDecompValue  (conv::PhysParamId_t p) const;
  double EigenDecompEigenValue    (conv::PhysParamId_t p) const;
   
  // Get BANFF best fit values, covariance matrix and cholesky decomposition of covariance matrix
  const TVectorD & BestFitValues              (void) const { return *fBestFitValues;             }
  const TVectorD & CurrentUncorrelatedValues  (void) const { return *fUncorrelatedCurrentValues; }
  const TVectorD & BestFitEigenDecompValues   (void) const { return *fEigenBestFitValues;        }
  const TVectorD & CurrentEigenDecompValues   (void) const { return *fEigenCurrentValues;        }
  const TVectorD & EigenDecompEigenValues     (void) const { return *fEigenValues;               }
  const TMatrixD & CovMatrix                  (void) const { return *fCovMtrx;                   }
  const TMatrixD & CovMatrixInv               (void) const { return *fCovMtrxInv;                }
  const TMatrixD & CovMatrixCholesky          (void) const { return *fCovMtrxChoL;               }

  // Systematic variations (wrt the best-fit BANFF params)
  void Tweak             (bool force_inputs = false, TVectorD * input_vec = 0, int mode = 1); ///> if correlated=false, still updates fCurrentValues using UpdateUncorrelated()
  void TweakCorrelated   (bool force_inputs = false, TVectorD * input_vec = 0); ///> updates fCurrentValues
  void TweakUncorrelated (bool force_inputs = false, TVectorD * input_vec = 0); ///> updates fUncorrelatedCurrentValues
  void TweakEigenDecomp  (bool force_inputs = false, TVectorD * input_vec = 0); ///> updates fEigenCurrentValues

  // A special systematic variation which cancels the effect of BANFF tuning and reverts back to the nominal MC
  void RestoreNominalMC(void);

  // Penalty term included in the chisq if the BANFF fit parameters are included in the oscillation fit.
  double ChisqPenalty             (bool all_syst = true, int mode = 1) const;
  double ChisqPenaltyCorrelated   (bool all_syst = true) const;
  double ChisqPenaltyUncorrelated (bool all_syst = true) const;
  double ChisqPenaltyEigenDecomp  (bool all_syst = true) const;

  // Get & tweak flags	
  bool IsTweaked           (void) const { return fIsTweaked;         } ///< BANFF params tweaked?
  bool IsFluxWghtUpdated   (void) const { return fIsFluxWghtUpdated; } ///< NuFlux object updated flux in response?
  void MarkFluxWghtUpdated (void)       { fIsFluxWghtUpdated = true; }

  // Print
  void PrintSystematicVariations (ostream & stream = std::cout) const;
   
  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

private:

  void Init      (void);
  void CleanUp   (void);
  int  Param2Idx   (conv::PhysParamId_t p) const; ///< convert BANFF PhysParamId_t to a covariance matrix or vector index
  int  UnParam2Idx (conv::PhysParamId_t p) const; ///< convert UnBANFF PhysParamId_t to a covariance matrix or vector index
  int  EDParam2Idx (conv::PhysParamId_t p) const; ///< convert EDBANFF PhysParamId_t to a covariance matrix or vector index
  void UpdateCorrelatedCholesky(void);    ///< use the currently stored uncorrelated values to update the current correlated values, using the cholesky method
  void UpdateCorrelatedEigenDecomp(void); ///< use the currently stored uncorrelated values to update the current correlated values, using the eigendecomposition method
  void CholeskyDecomposeMatrix (void);    ///< Cholesky decompose the covariance matrix

  TVectorD *  fBestFitValues;     ///< BANFF best fit values
  TVectorD *  fCurrentValues;     ///< Current values of BANFF parameters 
  TVectorD *  fUncorrelatedCurrentValues; ///< Uncorrelated Cholesky transformed values of BANFF current values
  TMatrixD *  fCovMtrx;           ///< BANFF covariance matrix
  TMatrixD *  fCovMtrxInv;        ///< Inverse BANFF covariance matrix
  TMatrixD *  fCovMtrxChoL;       ///< Triangular matrix obtained from Cholesky decomposition of BANFF covariance matrix
  bool        fIsTweaked;         ///< Were BANF params tweaked?
  bool        fIsFluxWghtUpdated; ///< Were Flux weights updated as a result of the last BANFF param variation
      
  TVectorD    * fEigenValues;           ///< Eigenvalues
  TMatrixD    * fEigenVectors;          ///< Eigenvectors
  TMatrixD    * fEigenVectorsT;         ///< Transpose of eigenvectors
  TMatrixDSym * fEigenCovMtrx;          ///< EigenDecomposition covariance matrix
  TVectorD    * fEigenBestFitValues;    ///< Bestfit values of BANFF parameters in eigendecomposition space
  TVectorD    * fEigenCurrentValues;    ///< Current values of BANFF parameters in eigendecomposition space

ClassDef(BANFF, 1)
};

} //t2kana

#endif
