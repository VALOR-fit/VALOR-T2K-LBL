//____________________________________________________________________________
/*!
   
\class     valor::T2KPhysicalParameterLimits
      
\brief     A set of physics parameter limits.
   
\author   Steve Dennis <stephen.dennis \at warwick.ac.uk>
          University of Warwick, Rutherford Appleton Laboratory

\created  February 2015

\copyright (c) 2010-2015, The VALOR group.

*/
//____________________________________________________________________________
   
#ifndef _valor_t2k_physical_limits_h_
#define _valor_t2k_physical_limits_h_

#include <iostream>
#include <string>
#include <map>

#include <TNamed.h>

#include "valor/Conventions/PhysParamId.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/PhysParamList.h"

namespace valor {

class T2KPhysicalParameterLimits
{
public:

  static T2KPhysicalParameterLimits * Instance();
  
  double MinValue    (conv::PhysParamId_t p) const;
  void   SetMinValue (conv::PhysParamId_t p, double v);
  double MaxValue    (conv::PhysParamId_t p) const;
  void   SetMaxValue (conv::PhysParamId_t p, double v);
  
  double GetNearestPhysicalValue(conv::PhysParamId_t p, double v) const;
  
  bool IsWithinLimits(conv::PhysParamId_t p, double v) const;
  bool AllWithinLimits(const FitParamList & fparm) const;
  bool AllWithinLimits(const PhysParamList & fparm) const;

  const map<int,double> & Minima(void) const;
  const map<int,double> & Maxima(void) const;
  
  void SetDefaults(void);
  void SetLimits  (FitParamList * fpl);

private:

  std::map<int,double> fParamMinima;
  std::map<int,double> fParamMaxima;
  
  
  T2KPhysicalParameterLimits();
  T2KPhysicalParameterLimits(const T2KPhysicalParameterLimits &);
 ~T2KPhysicalParameterLimits();

  static T2KPhysicalParameterLimits * fInstance;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (T2KPhysicalParameterLimits::fInstance !=0) {
            delete T2KPhysicalParameterLimits::fInstance;
            T2KPhysicalParameterLimits::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
  
};

} // valor

#endif
