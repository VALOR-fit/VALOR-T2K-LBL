//____________________________________________________________________________
/*!

\class     valor::T2KOpt

\brief     Singleton to keep global analysis options

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory

\created   October 20, 2010

\copyright (c) 2010-2014, The VALOR group.

21-7-2014 Davide
Added samples for multiring forward horn current
 
*/
//____________________________________________________________________________

#ifndef _valor_t2kopt_
#define _valor_t2kopt_

#include <iostream>
#include <TNamed.h>

#include "Conventions.h"
#include "Conventions/Opt.h"

namespace valor {

class T2KOpt {

public:
  static T2KOpt * Instance();


  // Choice of oscillation analysis type (combined, numu disapp only, nue app only)
  conv::OscAnaType_t gOscAnaTypeChoice;

  // Paths in xml
  
  string fConfigFileName;
  void ReadXmlPath(string config_file);

  string fPathBanffCovMatr;
  string fPathBanffParam;
  string fPathRHCBanffCovMatr;
  string fPathRHCBanffParam;

  string fPathEscale1RmuErr;
  string fPathEscale1ReErr;
  string fPathSKDetFSICovMatr;
  string fPathRHCEscale1RmuErr;
  string fPathRHCEscale1ReErr;
  string fPathRHCSKDetFSICovMatr;

  string fPathEscaleMultiRmuErr;
  string fPathEscaleMultiReErr;
  string fPathMultiRingSKDetFSICovMatr;

  //

  void Print(std::ostream & stream) const;
  friend std::ostream & operator << (std::ostream & stream, const T2KOpt & opt);

  // Check for any parameters being undefined
  bool Undefined(void);
  
  // These parameters are filled based on the condition of the compiler flags.
  // Altering them does not, at this point, change anything.
  bool fFit_1Rmu    ;
  bool fFit_1Re     ;
  bool fFit_RHC_1Rmu;
  bool fFit_RHC_1Re ;
  bool fFit_MultiRmu;
  bool fFit_MultiRe ;
  
private:

  T2KOpt();
  T2KOpt(const T2KOpt &);
 ~T2KOpt();

  static T2KOpt * fInstance;

  void SetDefaults(void);

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

std::ostream & operator<< (std::ostream & stream, const T2KOpt & plst);

}      // valor namespace
#endif // __valor_opt__
