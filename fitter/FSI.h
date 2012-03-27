//____________________________________________________________________________
/*!
   
\class    t2kana::FSI
      
\brief    Final state interaction (FSI) info used in the oscillation fit
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 21, 2012

*/
//____________________________________________________________________________
   
#ifndef _t2kana_fsi_h_
#define _t2kana_fsi_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TH1D.h>
#include <TH2D.h>

#include "Conventions.h"

using t2kana::conv::kns;
using t2kana::conv::knm;

namespace t2kana {

class FSI: public TObject
{
public:

  FSI();
 ~FSI();

  // Calculate weight in response to systematic parameter variations of FSI parameters
  double SystWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Tweak FSI params (either force the input params or generate correlated params using the Cholesky method)
  void Tweak (bool force_inputs = false, TVectorD * input_params = 0);

  // Get current value of systematic params handled by this class
  double CurrentValue (conv::PhysParamId_t p) const;

  // Get input FSI covariance matrix and its Cholesky decomposition
  const TMatrixD & CovMatrix         (void) const { return *fCovMtrx;       }
  const TMatrixD & CovMatrixCholesky (void) const { return *fCovMtrxChoL;   }

  // Penalty term included in the chisq if FSI nuisance terms are included in the fit       
  double ChisqPenalty(void) const;

  // Print
  void PrintSystematicVariations (void) const;

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  //
  void BuildFracMCTemplates(TH2D* input_mctempl_array[kns][knm]);
   
private:

  void Init      (void);
  void CleanUp   (void);
  int  Param2Idx (conv::PhysParamId_t p) const; ///< convert FSI PhysParamId_t to a covariance matrix or vector index

  //
  // Private data members
  //

  TVectorD *  fCurrentValues;    ///< Current values of FSI params
  TVectorD *  fNominalValues;    ///< Nominal values of FSI params
  TMatrixD *  fCovMtrx;          ///< FSI covariance matrix
  TMatrixD *  fCovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix
  TH1D *      fFracTun[kns];     ///<

ClassDef(FSI, 1);
};

} // t2kana

#endif
