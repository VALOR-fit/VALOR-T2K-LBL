//____________________________________________________________________________
/*!

\class     t2kana::Opt

\brief     Singleton to keep global analysis options

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory

\created   October 20, 2010

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_opt_
#define _t2kana_opt_

#include <iostream>
#include <TNamed.h>

#include "Conventions.h"

using std::ostream;

namespace t2kana {

class Opt {

public:
  class Persistable; // a lightweight class to extract state in a persistable named object
  static Opt * Instance();

  // Mass hierarchy 
  conv::NuMassHierarchy_t gMassHierarchyAssumed;
 
  // Choice of mixing parameter (sinsq2theta, sinsqtheta, ...)
  // Optimum choice depends on fit
  // Use sinsq2theta for the numu analysis, sinsqtheta for the global nue+numu analysis etc
  conv::NuOscMixParam_t gOscMixingParamChoice;

  // Choice of mixing parameter octant for theta23 (1, 2)
  // Only used when in sinsq2theta 3-flavour mode
  conv::NuOscOctant_t gOscOctantChoice;

  // Choice of independent squared mass splittings, eg                                                                                                                           
  // - dm^2 = Dm^2_{12} and Dm^2 = |Dm^2_{23}|, or                                                                                                                                 
  // - dm^2 = Dm^2_{12} and Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2| (Fogli,Lisi)                                                                                                  
  conv::NuOscDmsqScheme_t gOscDmsqScheme;

  // Choice of deltaCP, use it as purely phase or its sin
  conv::NuOscDeltaCP_t gOscDeltaCP;


  // Choice of oscillation analysis type (combined, numu disapp only, nue app only)
  conv::OscAnaType_t gOscAnaTypeChoice;

  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const Opt & opt);

  double GetGlobalBestFit     (conv::PhysParamId_t p) const;
  double GetGlobalFitError    (conv::PhysParamId_t p) const;
  double GetRandomGaus        (conv::PhysParamId_t p) const;
  
  void SetGlobalBestFit     (conv::PhysParamId_t p, double value);
  void SetGlobalFitError    (conv::PhysParamId_t p, double value);
  
  void SetDefaultParameterValues(void);

  static double BadResult(void) { return kBadResult;};

  // Check for any parameters being undefined
  bool Undefined(void);

  // 
  void WriteTo  (Persistable & optp);
  void ReadFrom (Persistable & optp);

  // nested lightweight class to extract state in a persistable named object
  class Persistable: public TNamed {
   public:
     Persistable() : 
     TNamed() 
     {
       mh    = conv::kNmhNormal;
       mixp  = conv::kNuMixSinsq2Theta;
       oct   = conv::kNuOctantOne;
       dmsqs = conv::kNuDmsq12andAbs23;
       dcp   = conv::kNuDeltaCPPhase;
     }
     Persistable(const Persistable & p) : 
     TNamed() 
     {
       mh    = p.mh;
       mixp  = p.mixp;
       oct   = p.oct;
       dmsqs = p.dmsqs;
       dcp   = p.dcp;
     }
    ~Persistable() { }
     conv::NuMassHierarchy_t mh;
     conv::NuOscMixParam_t   mixp;
     conv::NuOscOctant_t     oct;
     conv::NuOscDmsqScheme_t dmsqs;
     conv::NuOscDeltaCP_t    dcp;
   ClassDef(Persistable,1);
  };
private:

  Opt();
  Opt(const Opt &);
 ~Opt();

  static Opt * fInstance;

  void SetDefaults(void);

  double fGlobalBestFit       [conv::knppmax];
  double fGlobalFitError      [conv::knppmax];
  
  static const double kBadResult = -999;

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

ostream & operator<< (ostream & stream, const Opt & plst);

}      // t2kana namespace
#endif // __t2kana_opt__
