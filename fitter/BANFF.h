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
#include <TMatrixTSym.h>
#include <TXMLEngine.h>

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
   
  // Get BANFF best fit values, covariance matrix and cholesky decomposition of covariance matrix
  const TVectorD & BestFitValues     (void) const { return *fBestFitValues; }
  const TMatrixD & CovMatrix         (void) const { return *fCovMtrx;       }
  const TMatrixD & CovMatrixCholesky (void) const { return *fCovMtrxChoL;   }

  // Systematic variations (wrt the best-fit BANFF params)
  void Tweak (bool force_input_twkv = false, TVectorD * input_twkv = 0, bool correlated = true);

  // A special systematic variation which cancels the effect of BANFF tuning and reverts back to the nominal MC
  void RestoreNominalMC(void);

  // Get & tweak flags	
  bool IsTweaked           (void) const { return fIsTweaked;         } ///< BANFF params tweaked?
  bool IsFluxWghtUpdated   (void) const { return fIsFluxWghtUpdated; } ///< NuFlux object updated flux in response?
  void MarkFluxWghtUpdated (void)       { fIsFluxWghtUpdated = true; }

  // Print
  void PrintSystematicVariations (void) const;
   
  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

private:

  void Init      (void);
  void CleanUp   (void);
  int  Param2Idx (conv::PhysParamId_t p) const; ///< convert BANFF PhysParamId_t to a covariance matrix or vector index
  void DoTweak   (TVectorD * input_twkv); ///< apply input tweaking dials

  TVectorD *  fBestFitValues;     ///< BANFF best fit values
  TVectorD *  fCurrentValues;     ///< Current values of BANFF parameters (best_fit * twk_dial)
  TMatrixD *  fCovMtrx;           ///< BANFF covariance matrix
  TMatrixD *  fCovMtrxChoL;       ///< Triangular matrix obtained from Cholesky decomposition of BANFF covariance matrix
  bool        fIsTweaked;         ///< Were BANF params tweaked?
  bool        fIsFluxWghtUpdated; ///< Were Flux weights updated as a result of the last BANFF param variation
      
ClassDef(BANFF, 1)
};

} //t2kana

#endif
