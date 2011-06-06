//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 06/06/2011 - CA
   PhysHypothesis class added so as to simplify SuperMC

*/
//________________________________________________________________________________________

#include "PhysHypothesis.h"

using namespace t2kana;

//________________________________________________________________________________________
PhysHypothesis::PhysHypothesis()
{
  this->Init();
}
//________________________________________________________________________________________
PhysHypothesis::~PhysHypothesis()
{

}
//________________________________________________________________________________________
void PhysHypothesis::EnableOscillations(bool on_off) 
{ 
  fHaveNuOsc   = on_off; 
}
//________________________________________________________________________________________
void PhysHypothesis::EnableDecay(bool on_off) 
{ 
  fHaveNuDecay = on_off; 
}
//________________________________________________________________________________________
void PhysHypothesis::EnableDecoherence(bool on_off) 
{ 
  fHaveNuDecoh = on_off; 
}
//________________________________________________________________________________________
void PhysHypothesis::UseOscillationParams(double sinsq2theta23, double dm23sq) 
{ 
  fCurrSinsq2Theta23 = sinsq2theta23; 
  fCurrDm23sq        = dm23sq; 
}
//________________________________________________________________________________________
void PhysHypothesis::UseDecayParams(double sinsq2theta23, double alpha ) 
{ 
  fCurrSinsq2Theta23 = sinsq2theta23;
  fCurrAlphaDecay    = alpha;  
}
//________________________________________________________________________________________
void PhysHypothesis::UseDecoherenceParams(double sinsq2theta23, double musq) 
{ 
  fCurrSinsq2Theta23 = sinsq2theta23; 
  fCurrMusqDecoh     = musq;   
}
//________________________________________________________________________________________
bool PhysHypothesis::Oscillations(void) const
{
  return fHaveNuOsc;
}
//________________________________________________________________________________________
bool PhysHypothesis::Decay(void) const
{
  return fHaveNuDecay;
}
//________________________________________________________________________________________
bool PhysHypothesis::Decoherence(void) const
{
  return fHaveNuDecoh;
}
//________________________________________________________________________________________
double PhysHypothesis::CurrSinsq2Theta23(void) const
{
  return fCurrSinsq2Theta23;
}
//________________________________________________________________________________________
double PhysHypothesis::CurrDm23sq(void) const
{
  return fCurrDm23sq;
}
//________________________________________________________________________________________
double PhysHypothesis::CurrAlphaDecay(void) const
{
  return fCurrAlphaDecay;
}
//________________________________________________________________________________________
double PhysHypothesis::CurrMusqDecoh(void) const
{
  return fCurrMusqDecoh;
}
//________________________________________________________________________________________
void PhysHypothesis::Init(void)
{
  fHaveNuOsc         = true;
  fHaveNuDecay       = false;
  fHaveNuDecoh       = false;
  fCurrSinsq2Theta23 = 0.;
  fCurrDm23sq        = 0.;
  fCurrAlphaDecay    = 0.;
  fCurrMusqDecoh     = 0.;  
}
//________________________________________________________________________________________

