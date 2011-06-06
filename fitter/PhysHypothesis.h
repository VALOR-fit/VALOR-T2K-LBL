//____________________________________________________________________________
/*!
   
\class    t2kana::PhysHypothesis
      
\brief    An hypothesis for the fitted SuperK distributions (oscillation,
          decay, decoherence)
   
\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups
  
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>   

\created  June 06, 2011

*/
//____________________________________________________________________________
   
#ifndef _t2kana_phys_hypothesis_h_
#define _t2kana_phys_hypothesis_h_

#include "Conventions.h"

namespace t2kana {

class PhysHypothesis
{
public:
   PhysHypothesis();
  ~PhysHypothesis();

  void EnableOscillations   (bool on_off);
  void EnableDecay          (bool on_off);
  void EnableDecoherence    (bool on_off);
  void UseOscillationParams (double sinsq2theta23, double dm23sq); 
  void UseDecayParams       (double sinsq2theta23, double alpha ); 
  void UseDecoherenceParams (double sinsq2theta23, double musq  ); 

  bool   Oscillations      (void) const; ///<
  bool   Decay             (void) const; ///<
  bool   Decoherence       (void) const; ///<
  double CurrSinsq2Theta23 (void) const; ///<
  double CurrDm23sq        (void) const; ///<
  double CurrAlphaDecay    (void) const; ///<
  double CurrMusqDecoh     (void) const; ///<

private:

  void Init(void);

  bool    fHaveNuOsc;
  bool    fHaveNuDecay;
  bool    fHaveNuDecoh;
  double  fCurrSinsq2Theta23;
  double  fCurrDm23sq;
  double  fCurrAlphaDecay;
  double  fCurrMusqDecoh;  
};

} // t2kana

#endif
