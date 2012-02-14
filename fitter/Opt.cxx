//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

#include <string>

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
  gOscAnaTypeChoice     = kAnaCombined;
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
  stream << " - Oscillation analysis  choice: " 
         << util::AsString(gOscAnaTypeChoice) 
         << endl;
}
//____________________________________________________________________________

