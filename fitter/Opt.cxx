//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @ 18/2/13 - LE                                                                                                                               
 Added back the selection of dmsq schema,from branch 2010_v3: choice between usual Dm12sq + Dm23sq and Fogli&Lisi schema

*/
//________________________________________________________________________________________

#include <string>
#include <cstdlib>

#include "Opt.h"
#include "Utilities.h"

using std::endl;
using std::string;

using namespace t2kana;
using namespace t2kana::conv;

//____________________________________________________________________________
namespace t2kana {
  ostream & operator<<(ostream & stream, const Opt & opt)
  {
    opt.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
Opt * Opt::fInstance = 0;
//____________________________________________________________________________
Opt::Opt()
{
  this->SetDefaults();
  fInstance = 0;
}
//____________________________________________________________________________
Opt::~Opt()
{

}
//____________________________________________________________________________
Opt * Opt::Instance()
{
  if(fInstance == 0) {
    static Opt::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new Opt;
  }
  return fInstance;
}
//____________________________________________________________________________
void Opt::SetDefaults(void) 
{
  gMassHierarchyAssumed = kNmhNormal;
  gOscMixingParamChoice = kNuMixSinsq2Theta;
  gOscDmsqScheme        = kNuDmsq12andAbs23;

  bool fit_1Rmu = false;
  bool fit_1Re  = false;

  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsFCFV1RingMuLike) fit_1Rmu = true;
    if(s == kSKsFCFV1RingELike ) fit_1Re  = true;
  }

  if(fit_1Rmu && fit_1Re) {
    gOscAnaTypeChoice = kAnaCombined;
  }
  else
  if(fit_1Rmu && !fit_1Re) {
    gOscAnaTypeChoice = kAnaNuMuDisappOnly;
  }
  else
  if(!fit_1Rmu && fit_1Re) {
    gOscAnaTypeChoice = kAnaNuEAppOnly;
  }
  else {
    gOscAnaTypeChoice = kAnaUndefined;
  }

}
//____________________________________________________________________________
void Opt::Print(ostream & stream) const
{
  stream << "[Global analysis options]" 
         << endl;
  stream << " - Hierarchy: " 
         << util::AsString(gMassHierarchyAssumed) 
         << endl;
  stream << " - Mixing param choice: " 
         << util::AsString(gOscMixingParamChoice) 
         << endl;
  stream << " - Independent squared mass splittings: "
         << util::AsString(gOscDmsqScheme)
         << endl;
  stream << " - Oscillation analysis choice: " 
         << util::AsString(gOscAnaTypeChoice) 
         << endl;
}
//____________________________________________________________________________

