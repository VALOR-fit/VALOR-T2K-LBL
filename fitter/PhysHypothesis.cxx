//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 06/06/2011 - CA
   PhysHypothesis class added so as to simplify SuperMC
 @ 21/06/2011 - CA
   Add PhysHypothesis::Copy(const PhysHypothesis & ph).
 @ 05/07/2011 - CA
   The class was restructured, using the conv::PhysHypothesis_t enumeration and an STL
   map to store the actual params. The interface for setting/getting param values matched 
   the one in FitParamList.
 @ 11/07/2011 - CA
   Implemented kHypOsc3F case.
 @ 18/07/2011 - CA
   Bug fix: Reverse order of kHypOsc3F and kHypNuMuDisapp2F cases (kHypOsc3F now first) 
   in PhysHypothesis(const FitParamList &) ctor so that kHypOsc3F case can be picked.
 @ 27/07/2011 - CA
   Added Print(ostream &) and `<<' operator.
 @ 14/09/2011 - CA
   Added PhysHypothesis::PhysHypothesis(string h) constructor so as to simplify some of
   the test macros.

*/
//________________________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__

#include <iomanip>

#include "PhysHypothesis.h"
#include "Utilities.h"

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________

namespace t2kana {
  ostream & operator << (ostream & stream, const PhysHypothesis & ph)
  {
    ph.Print(stream);
    return stream;
  }
}
//________________________________________________________________________________________
PhysHypothesis::PhysHypothesis(PhysHypothesis_t h)
{
  this->Init(h);
}
//________________________________________________________________________________________
PhysHypothesis::PhysHypothesis(string sh)
{
  PhysHypothesis_t h = util::FromString(sh);
  this->Init(h);
}
//________________________________________________________________________________________
PhysHypothesis::PhysHypothesis(const FitParamList & fp)
{
  PhysHypothesis_t hyp = fp.PhysicsHypothesis();

  // neutrino oscillations (3-flavour framework)
  if(hyp == kHypOsc3F)
  {
    this->Init(kHypOsc3F);
    this->SetCurrValue(kPhypDm12sq,        fp.CurrValue(kPhypDm12sq));
    this->SetCurrValue(kPhypDm23sq,        fp.CurrValue(kPhypDm23sq));
    this->SetCurrValue(kPhypSinsq2Theta12, fp.CurrValue(kPhypSinsq2Theta12));
    this->SetCurrValue(kPhypSinsq2Theta13, fp.CurrValue(kPhypSinsq2Theta13));
    this->SetCurrValue(kPhypSinsq2Theta23, fp.CurrValue(kPhypSinsq2Theta23));
    this->SetCurrValue(kPhypDeltaCP,       fp.CurrValue(kPhypDeltaCP));
  }

  // numu disappearance, 2-flavour neutrino oscillation framework 
  else if(hyp == kHypNuMuDisapp2F) 
  {
    this->Init(kHypNuMuDisapp2F);
    this->SetCurrValue(kPhypDm23sq,        fp.CurrValue(kPhypDm23sq));
    this->SetCurrValue(kPhypSinsq2Theta23, fp.CurrValue(kPhypSinsq2Theta23));
  }

  // neutrino decay
  else if(hyp == kHypDecay) 
  {
    this->Init(kHypDecay);
    this->SetCurrValue(kPhypAlphaDecay,       fp.CurrValue(kPhypAlphaDecay));
    this->SetCurrValue(kPhypSinsq2ThetaDecay, fp.CurrValue(kPhypSinsq2ThetaDecay));
  }

  // neutrino decoherence
  else if(hyp == kHypDecoh) 
  {
    this->Init(kHypDecoh);
    this->SetCurrValue(kPhypMusqDecoh,          fp.CurrValue(kPhypMusqDecoh));
    this->SetCurrValue(kPhypSinsq2ThetaDecoh,   fp.CurrValue(kPhypSinsq2ThetaDecoh));
  }

  // undefined
  else {  
    this->Init(kHypUndefined);
  }

#ifdef __debug_level_2__
  cout << " PhysHypothesis ctor > "
       << "Input physics hypothesis = " << util::AsString(fHypothesis) << endl;
#endif
}
//________________________________________________________________________________________
PhysHypothesis::~PhysHypothesis()
{
  fHypothesis = kHypUndefined;
  fParams.clear();
}
//________________________________________________________________________________________
double PhysHypothesis::CurrValue(PhysParamId_t p) const
{
  map<int,double>::const_iterator it = fParams.find(p);
  if(it == fParams.end()) {
     exit(1);
  } 
  return it->second;
}
//________________________________________________________________________________________
void PhysHypothesis::SetCurrValue(PhysParamId_t p, double v)
{
  fParams[p] = v;
}
//________________________________________________________________________________________
void PhysHypothesis::Copy(const PhysHypothesis & ph)
{
  this->Init(ph.fHypothesis);
  map<int,double>::const_iterator it = ph.fParams.begin();
  for( ; it != ph.fParams.end(); ++it) {
     fParams[it->first] = it->second;
  }
}
//________________________________________________________________________________________
void PhysHypothesis::Init(PhysHypothesis_t h)
{
  fHypothesis = h;

  fParams.clear();

  // numu disappearance, 2-flavour neutrino oscillation framework 
  if(fHypothesis == kHypNuMuDisapp2F) {
    fParams[kPhypSinsq2Theta23]    = 0;
    fParams[kPhypDm23sq]           = 0;
  }
  // neutrino oscillations (3-flavour framework)
  if(fHypothesis == kHypOsc3F) {
    fParams[kPhypSinsq2Theta12]    = 0;
    fParams[kPhypSinsq2Theta13]    = 0;
    fParams[kPhypSinsq2Theta23]    = 0;
    fParams[kPhypDeltaCP]          = 0;
    fParams[kPhypDm23sq]           = 0;
    fParams[kPhypDm12sq]           = 0;
  }
  // neutrino decay 
  if(fHypothesis == kHypDecay) {
    fParams[kPhypSinsq2ThetaDecay] = 0;
    fParams[kPhypAlphaDecay]       = 0;
  }
  // neutrino decoherence 
  if(fHypothesis == kHypDecoh) {
    fParams[kPhypSinsq2ThetaDecoh] = 0;
    fParams[kPhypMusqDecoh]        = 0;
  }

}
//________________________________________________________________________________________
void PhysHypothesis::Print(ostream & stream) const
{
  stream << endl;
  stream << "................................................................." << endl;
  stream << "Hypothesis: \n" 
         << util::AsString(fHypothesis) << endl;
  stream << "Parameters: " << endl;
  map<int,double>::const_iterator it = fParams.begin();
  for( ; it != fParams.end(); ++it) {
     PhysParamId_t p = (PhysParamId_t) it->first;
     double        v = it->second;
     stream << " - "
            << setw(15)  << util::AsString(p) << " = "
            << setprecision(8) << setw(12) << v
            << endl;
  }
  stream << "................................................................." << endl;
}
//________________________________________________________________________________________

