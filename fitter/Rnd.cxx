//________________________________________________________________________________________
/*

*/
//________________________________________________________________________________________

#include "Rnd.h"
#include "Build.h"

using std::endl;

using namespace t2kana;

//____________________________________________________________________________
namespace t2kana {
  ostream & operator<<(ostream & stream, const Rnd & rnd)
  {
    rnd.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
Rnd * Rnd::fInstance = 0;
//____________________________________________________________________________
Rnd::Rnd()
{
  this->SetSeeds(); // default seed: 0
  fInstance = 0;
}
//____________________________________________________________________________
Rnd::~Rnd()
{

}
//____________________________________________________________________________
Rnd * Rnd::Instance()
{
  if(fInstance == 0) {
    static Rnd::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new Rnd;
  }
  return fInstance;
}
//____________________________________________________________________________
void Rnd::SetSeeds(
  UInt_t toymcstat, UInt_t toymcsyst, UInt_t toymcosc, UInt_t fitter, UInt_t other)
{
  fSeedToyMCstat = toymcstat;
  fSeedToyMCsyst = toymcsyst;
  fSeedToyMCosc  = toymcosc ;
  fSeedFitter    = fitter   ;
  fSeedOther     = other    ;
  
#ifndef __SEED_TRANDOM_ZERO__
  unsigned int randseeds[5] = {};
  // If any seeds are zero, read from /dev/urandom to seed them.
  // /dev/random might be slow on batch machines
  if (! (toymcstat && toymcsyst && toymcosc && fitter && other) ) {
    FILE * devrandom = fopen("/dev/urandom","r");
    fread(randseeds, sizeof(unsigned int), 5, devrandom);
    fclose(devrandom);
  }
    
  fSeedToyMCstat = ( toymcstat ? toymcstat : randseeds[0] );
  fSeedToyMCsyst = ( toymcsyst ? toymcsyst : randseeds[1] );
  fSeedToyMCosc  = ( toymcosc  ? toymcosc  : randseeds[2] );
  fSeedFitter    = ( fitter    ? fitter    : randseeds[3] );
  fSeedOther     = ( other     ? other     : randseeds[4] );
#endif

  fRndToyMCstat . SetSeed(fSeedToyMCstat);
  fRndToyMCsyst . SetSeed(fSeedToyMCsyst);
  fRndToyMCosc  . SetSeed(fSeedToyMCosc);
  fRndFitter    . SetSeed(fSeedFitter);
  fRndOther     . SetSeed(fSeedOther);
  
}
//____________________________________________________________________________
void Rnd::Print(ostream & stream) const
{
  stream << "[Random number generator seed numbers]" 
         << endl;

  stream << " - Toy MC (Stat) : ";
  if (fSeedToyMCstat != 0) {
    stream << fRndToyMCstat.GetSeed();
  } else {
    stream << " ** randomized ** (" << fRndToyMCstat.GetSeed() << ")";
  }
  stream << endl;

  stream << " - Toy MC (Syst) : ";
  if (fSeedToyMCsyst != 0) {
    stream << fRndToyMCsyst.GetSeed();
  } else {
    stream << " ** randomized ** (" << fRndToyMCsyst.GetSeed() << ")";
  }
  stream << endl;

  stream << " - Toy MC (Osc) : ";
  if (fSeedToyMCosc != 0) {
    stream << fRndToyMCosc.GetSeed();
  } else {
    stream << " ** randomized ** (" << fRndToyMCosc.GetSeed() << ")";
  }
  stream << endl;

  stream << " - Fitter : ";
  if (fSeedFitter != 0) {
    stream << fRndFitter.GetSeed();
  } else {
    stream << " ** randomized ** (" << fRndFitter.GetSeed() << ")";
  }
  stream << endl;

  stream << " - Other : ";
  if (fSeedOther != 0) {
    stream << fRndOther.GetSeed();
  } else {
    stream << " ** randomized ** (" << fRndOther.GetSeed() << ")";
  }
  stream << endl;
}
//____________________________________________________________________________

