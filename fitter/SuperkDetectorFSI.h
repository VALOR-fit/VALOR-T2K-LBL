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

#include "valor/FitCore/CorrelatedSystematics.h"

#include "Enumerations.h"
#include "Conventions.h"
#include "Utilities.h"

namespace valor {
  
class FitParamList;
class PhysParamList;

class SuperkDetectorFSI: public CorrelatedSystematics
{
public:

  SuperkDetectorFSI();
  virtual ~SuperkDetectorFSI();
  
  // What weight to apply
  double SystWght(
    conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double /*Etrue*/, std::string option = "") const;
    
  double EnergyScale(conv::SuperkSample_t /*s*/ = conv::kSKsUndefined) const;
      
  // Apply the given FSI parameter to the given sample, mode, Ereco?
  bool ApplySyst            (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco) const;
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco) const;

private:

  // Check that SKFSISystConv.h is setup correctly
  void CheckListOfSKDetFSISyst(int nparams, const conv::PhysParamId_t * params, const char * vname, bool allow_gaps);
  void CheckERange            (int nparams, const conv::PhysParamId_t * params, const double * erange, const char * vname);
  void CheckSKMode            (int nparams, const conv::PhysParamId_t * params, const bool * modes, const char * vname, bool fatal_no_effect);
  
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
