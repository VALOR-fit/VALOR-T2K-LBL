//____________________________________________________________________________
/*!

\class     valor::T2KOpt

\brief     Singleton to keep global analysis options

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory

\created   October 20, 2010

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _valor_t2kopt_
#define _valor_t2kopt_

#include <iostream>
#include <TNamed.h>

#include "Conventions.h"
#include "Conventions/Opt.h"

using std::ostream;

namespace valor {

class T2KOpt {

public:
  static T2KOpt * Instance();


  // Choice of oscillation analysis type (combined, numu disapp only, nue app only)
  conv::OscAnaType_t gOscAnaTypeChoice;

  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const T2KOpt & opt);

  double GetGlobalBestFit     (conv::PhysParamId_t p) const;
  double GetGlobalFitError    (conv::PhysParamId_t p) const;
  double GetRandomGaus        (conv::PhysParamId_t p) const;
  
  void SetGlobalBestFit     (conv::PhysParamId_t p, double value);
  void SetGlobalFitError    (conv::PhysParamId_t p, double value);
  
  void SetDefaultParameterValues(void);

  static double BadResult(void) { return kBadResult;};

  // Check for any parameters being undefined
  bool Undefined(void);

private:

  T2KOpt();
  T2KOpt(const T2KOpt &);
 ~T2KOpt();

  static T2KOpt * fInstance;

  void SetDefaults(void);

  double fGlobalBestFit       [conv::knppmax];
  double fGlobalFitError      [conv::knppmax];
  
  static const double kBadResult = -999;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (T2KOpt::fInstance !=0) {
            delete T2KOpt::fInstance;
            T2KOpt::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

ostream & operator<< (ostream & stream, const T2KOpt & plst);

}      // valor namespace
#endif // __valor_opt__
