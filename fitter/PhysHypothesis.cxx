//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

//#define __show_init_mesg_info__

#include <iomanip>
#include <cstdlib>

#include "PhysParamList.h"
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
  PhysHypothesis_t h = util::PhysHypFromString(sh);
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
    this->SetCurrValue(kPhypDm12sq,   fp.CurrValue(kPhypDm12sq));
    this->SetCurrValue(kPhypDm23sq,   fp.CurrValue(kPhypDm23sq));
    this->SetCurrValue(kPhypMixing12, fp.CurrValue(kPhypMixing12));
    this->SetCurrValue(kPhypMixing13, fp.CurrValue(kPhypMixing13));
    this->SetCurrValue(kPhypMixing23, fp.CurrValue(kPhypMixing23));
    this->SetCurrValue(kPhypDeltaCP,  fp.CurrValue(kPhypDeltaCP));
  }

  // numu disappearance, 2-flavour neutrino oscillation framework 
  else if(hyp == kHypNuMuDisapp2F) 
  {
    this->Init(kHypNuMuDisapp2F);
    this->SetCurrValue(kPhypDm23sq,   fp.CurrValue(kPhypDm23sq));
    this->SetCurrValue(kPhypMixing23, fp.CurrValue(kPhypMixing23));
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
  
  // Steriles
  else if(hyp == kHypOscSterile4F || hyp == kHypOscSterile5F ) {
    this->Init(hyp);
    for (int i = 0 ; i < fp.Size() ; i++) {
      PhysParamId_t p = fp.WhichParam(i);
      this->SetCurrValue(p, fp.CurrValue(p));
    }
  }

  // undefined
  else {  
    this->Init(kHypUndefined);
  }

#ifdef __show_init_mesg_info__
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
    fParams[kPhypMixing23] = 0;
    fParams[kPhypDm23sq]   = 0;
  }
  // neutrino oscillations (3-flavour framework)
  if(fHypothesis == kHypOsc3F) {
    fParams[kPhypMixing12]  = 0;
    fParams[kPhypMixing13]  = 0;
    fParams[kPhypMixing23]  = 0;
    fParams[kPhypDeltaCP]   = 0;
    fParams[kPhypDm23sq]    = 0;
    fParams[kPhypDm12sq]    = 0;
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
PhysParamList PhysHypothesis::AsPhysParamList(void) const
{
  PhysParamList pl;
  map<int, double>::const_iterator iter = fParams.begin();
  for( ; iter != fParams.end(); ++iter)
  {
    PhysParamId_t p = (PhysParamId_t)iter->first;
    double v = iter->second;
    pl.SetValue(p,v);
  }
  return pl;
}
//________________________________________________________________________________________

