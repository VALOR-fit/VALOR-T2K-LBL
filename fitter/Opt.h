//____________________________________________________________________________
/*!

\class    t2kana::Opt

\brief    Singleton to keep global analysis options

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef __t2kana_opt__
#define __t2kana_opt__

#include <iostream>

#include "Conventions.h"

using std::ostream;

namespace t2kana {

class Opt {

public:
  static Opt * Instance();

  // Mass hierarchy 
  conv::NuMassHierarchy_t gMassHierarchyAssumed;
 
  // Choice of mixing parameter (sinsq2theta, sinsqtheta, ...)
  // Optimum choice depends on fit
  // Use sinsq2theta for the numu analysis, sinsqtheta for the global nue+numu analysis etc
  conv::NuOscMixParam_t gOscMixingParamChoice;

  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const Opt & opt);

private:

  Opt();
  Opt(const Opt & opt);
 ~Opt();

  static Opt * fInstance;

  void SetDefaults(void);

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (Opt::fInstance !=0) {
            delete Opt::fInstance;
            Opt::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // t2kana namespace
#endif // __t2kana_opt__
