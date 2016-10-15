
//____________________________________________________________________________
/*!

\class    valor::BANFF

\brief    Class to hold analysis inputs from BANFF (Beam and ND280 Flux
          Extrapolation Task Force). Now a CorrelatedSystematics() subclass.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
          Steve Dennis <s.r.dennis \at University of Liverpool>
          University of Liverpool

\created  February 09, 2012

\copyright (c) 2010-2016, VaLOR group.

*/
//____________________________________________________________________________

#ifndef _valor_banff_h_
#define _valor_banff_h_

#include "valor/FitCore/CorrelatedSystematics.h"

namespace valor {

class BANFF: public CorrelatedSystematics
{
public:
  BANFF();   
 ~BANFF();
  void TweakSampleConstrainedOnly  (bool force_inputs = false, TVectorD * input_vec = 0, bool unconstrainted=false);
  void TweakSampleFluxOnly  (bool force_inputs = false, TVectorD * input_vec = 0);
  void TweakSampleXSecOnly  (bool force_inputs = false, TVectorD * input_vec = 0);
  void TweakSampleXSecNumuNueOnly     (bool force_inputs = false, TVectorD * input_vec = 0);
  void TweakSampleMEC  (bool force_inputs = false, TVectorD * input_vec = 0);
  void TweakSampleMECBar  (bool force_inputs = false, TVectorD * input_vec = 0);
  void TweakSampleNCGamma  (bool force_inputs = false, TVectorD * input_vec = 0);

  
 // The old interface is gone. It never needed to exist anyway.
 // Here's what to replace your calls with:
 // Read*               -> Read
 // ChisqPenalty*       -> ChisqPenalty
 // Print*              -> PrintSystematicVariations
 // *CurrentValue*      -> CurrentValues
 // *BestFitValue*      -> InitialValue
 // *Tweak*             -> Tweak
 // Update/Split/Merge* -> [removed]
 // Restore*Nominal     -> Tweak("untuned")
 // Restore*Tuned       -> Tweak(1,0)
 // Asterisks represent the plethora of useless variants we had in the old interface.
 
private:
  
ClassDef(BANFF, 2)
};

} //valor

#endif
