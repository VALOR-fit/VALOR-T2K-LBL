//____________________________________________________________________________
/*!

\class     t2kana::ParameterStorage

\brief     Singleton to keep global best fit params to throw around

\author   Steve Dennis <stephen.dennis \at warwick.ac.uk>
          University of Warwick ; STFC, Rutherford Appleton Laboratory
          
\created  May 01 2014

Changes since 2014 v0.0.1:

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_paramstore_
#define _t2kana_paramstore_

#include <iostream>

#include "Conventions.h"
#include "PhysParamId.h"

using std::ostream;

namespace t2kana {

class ParameterStorage {

public:

  static ParameterStorage * Instance();
  
  double GetGlobalBestFit     (conv::PhysParamId_t p) const;
  double GetGlobalFitError    (conv::PhysParamId_t p) const;
  double GetInitialValue      (conv::PhysParamId_t p) const;
  double GetRandomGaus        (conv::PhysParamId_t p) const;
  
  void SetGlobalBestFit     (conv::PhysParamId_t p, double value);
  void SetGlobalFitError    (conv::PhysParamId_t p, double value);
  void SetInitialValue      (conv::PhysParamId_t p, double value);
  
  void SetDefaults(void);
  
  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const ParameterStorage & obj);
  
  static double BadResult(void) { return kBadResult;};
  
private:

  double fGlobalBestFit       [conv::knppmax];
  double fGlobalFitError      [conv::knppmax];
  double fInitialValue        [conv::knppmax];

  ParameterStorage();
  ParameterStorage(const ParameterStorage &);
 ~ParameterStorage();

  static ParameterStorage * fInstance;
  
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (ParameterStorage::fInstance !=0) {
            delete ParameterStorage::fInstance;
            ParameterStorage::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
  
  static const double kBadResult = -999.;
  
  
  //ClassDef(ParameterStorage,1);
};

ostream & operator<< (ostream & stream, const ParameterStorage & plst);

}      // t2kana namespace
#endif // __t2kana_paramstore_
