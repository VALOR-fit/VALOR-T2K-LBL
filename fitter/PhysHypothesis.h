//____________________________________________________________________________
/*!
   
\class    t2kana::PhysHypothesis
      
\brief    A physics hypothesis (oscillations, decay, decoherence, ...)
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  June 06, 2011

*/
//____________________________________________________________________________
   
#ifndef _t2kana_phys_hypothesis_h_
#define _t2kana_phys_hypothesis_h_

#include <map>

#include "Conventions.h"
#include "FitParamList.h"

using std::map;

namespace t2kana {

class PhysHypothesis
{
public:
   PhysHypothesis(conv::PhysHypothesis_t h = conv::kHypUndefined);
   PhysHypothesis(const FitParamList & fp);
  ~PhysHypothesis();

   void Init (conv::PhysHypothesis_t h);
   void Copy (const PhysHypothesis & ph);

   conv::PhysHypothesis_t Hypothesis (void) const { return fHypothesis; }

   double CurrValue    (conv::PhysParamId_t p) const;
   void   SetCurrValue (conv::PhysParamId_t p, double v);

private:

   conv::PhysHypothesis_t fHypothesis;
   map<int,double>        fParams;
};

} // t2kana

#endif
