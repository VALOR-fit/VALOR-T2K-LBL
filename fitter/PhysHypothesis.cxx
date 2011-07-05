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
*/
//________________________________________________________________________________________

#include "PhysHypothesis.h"

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
PhysHypothesis::PhysHypothesis(PhysHypothesis_t h)
{
  this->Init(h);
}
//________________________________________________________________________________________
PhysHypothesis::PhysHypothesis(const FitParamList & fp)
{
  if( fp.Exists(kPhypDm23sq) && fp.Exists(kPhypSinsq2Theta23) ) 
  {
    this->Init(kHypNuMuDisapp2F);
    this->SetCurrValue(kPhypDm23sq,        fp.CurrValue(kPhypDm23sq));
    this->SetCurrValue(kPhypSinsq2Theta23, fp.CurrValue(kPhypSinsq2Theta23));
  }
  else
  if( fp.Exists(kPhypAlphaDecay) && fp.Exists(kPhypSinsq2Theta23) )
  {
    this->Init(kHypDecay);
    this->SetCurrValue(kPhypAlphaDecay,    fp.CurrValue(kPhypAlphaDecay));
    this->SetCurrValue(kPhypSinsq2Theta23, fp.CurrValue(kPhypSinsq2Theta23));
  }
  else
  if( fp.Exists(kPhypMusqDecoherence) && fp.Exists(kPhypSinsq2Theta23) )
  {
    this->Init(kHypDecoh);
    this->SetCurrValue(kPhypMusqDecoherence, fp.CurrValue(kPhypMusqDecoherence));
    this->SetCurrValue(kPhypSinsq2Theta23,   fp.CurrValue(kPhypSinsq2Theta23));
  }
  else {  
    this->Init(kHypUndefined);
  }
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

  // neutrino decay 
  if(fHypothesis == kHypDecay) {
    fParams[kPhypSinsq2Theta23]   = 0;
    fParams[kPhypAlphaDecay]      = 0;
  }
  // neutrino decoherence 
  if(fHypothesis == kHypDecoh) {
    fParams[kPhypSinsq2Theta23]   = 0;
    fParams[kPhypMusqDecoherence] = 0;
  }
  // numu disappearance, 2-flavour neutrino oscillation framework 
  if(fHypothesis == kHypNuMuDisapp2F) {
    fParams[kPhypSinsq2Theta23] = 0;
    fParams[kPhypDm23sq]        = 0;
  }
  // neutrino oscillations (3-flavour framework)
  if(fHypothesis == kHypOsc3F) {
  }

}
//________________________________________________________________________________________

