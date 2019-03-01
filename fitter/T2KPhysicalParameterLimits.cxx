//________________________________________________________________________________________
/*

Copyright (c) 2010-2015, The VALOR group. All rights reserved.

ChangeLog:
 
23-06-2016 Chris
  Added method to set the physical limits of the parameters in the FitParamList,
  limits are now abstracted out from the framework

*/
//________________________________________________________________________________________

#include <string>
#include <algorithm>
#include <iostream>

#include <TROOT.h>

#include "T2KPhysicalParameterLimits.h"
#include "Utilities.h"
#include "Enumerations.h"
#include "valor/Messenger/Messenger.h"

using namespace valor;
using namespace valor::conv;

using std::endl;
using std::cout;
using std::cerr;
using std::map;

//#define show_debug_limits

//____________________________________________________________________________
T2KPhysicalParameterLimits * T2KPhysicalParameterLimits::fInstance = 0;
//____________________________________________________________________________
valor::T2KPhysicalParameterLimits::T2KPhysicalParameterLimits()
{
  this->SetDefaults();
  fInstance = 0;
}
//____________________________________________________________________________
valor::T2KPhysicalParameterLimits::~T2KPhysicalParameterLimits()
{
  return;
}
//____________________________________________________________________________
T2KPhysicalParameterLimits * T2KPhysicalParameterLimits::Instance()
{
  if(fInstance == 0) {
    static T2KPhysicalParameterLimits::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new T2KPhysicalParameterLimits;
  }
  return fInstance;
}
//____________________________________________________________________________
double  valor::T2KPhysicalParameterLimits::MinValue    (conv::PhysParamId_t p) const
{
  map<int,double>::const_iterator it = fParamMinima.find(p);
  if(it == fParamMinima.end()) {
     cerr << "T2KPhysicalParameterLimits >> Parameter [" 
          << util::RawParamAsString(p) << "] was not set!" << endl;
     return -99999;
  } 
  return it->second;
}
//____________________________________________________________________________
void    valor::T2KPhysicalParameterLimits::SetMinValue (conv::PhysParamId_t p, double v)
{
  fParamMinima[p] = v;
}
//____________________________________________________________________________
double  valor::T2KPhysicalParameterLimits::MaxValue    (conv::PhysParamId_t p) const
{
  map<int,double>::const_iterator it = fParamMaxima.find(p);
  if(it == fParamMaxima.end()) {
     cerr << "T2KPhysicalParameterLimits >> Parameter [" 
          << util::RawParamAsString(p) << "] was not set!" << endl;
     return -99999;
  } 
  return it->second;
}
//____________________________________________________________________________
void    valor::T2KPhysicalParameterLimits::SetMaxValue (conv::PhysParamId_t p, double v)
{
  fParamMaxima[p] = v;
}
//____________________________________________________________________________
bool valor::T2KPhysicalParameterLimits::IsWithinLimits(conv::PhysParamId_t p, double v) const
{
  if (v >= this->MinValue(p) && v <= this->MaxValue(p)) {
    return true;
  }
#ifdef show_debug_limits
  std::cout<<"T2KPhysicalParameterLimits/IsWithinLimits >> False for param "<<util::AsString(p)<<" with value "<<v<<std::endl;
#endif
  return false;
}
//____________________________________________________________________________
bool valor::T2KPhysicalParameterLimits::AllWithinLimits(const FitParamList & fparm) const
{
  for (int iparam = 0 ; iparam < fparm.Size(); iparam++) {
    PhysParamId_t p = fparm.WhichParam(iparam);
    double value = fparm.CurrValueByIndex(iparam);
    if (not this->IsWithinLimits(p,value)) {
      return false;
    }
  }
  return true;
}
//____________________________________________________________________________
bool valor::T2KPhysicalParameterLimits::AllWithinLimits(const PhysParamList & pparm) const
{
  for (int iparam = 0 ; iparam < conv::knppmax; iparam++) {
    PhysParamId_t p = PhysParamId_t(iparam);
    if (not pparm.Contains(p)) continue;
    double value = pparm.Value(p);
    if (not this->IsWithinLimits(p,value)) {
      return false;
    }
  }
  return true;
}
//____________________________________________________________________________
const map<int,double> & valor::T2KPhysicalParameterLimits::Minima(void) const
{
  return fParamMinima;
}
//____________________________________________________________________________
const map<int,double> & valor::T2KPhysicalParameterLimits::Maxima(void) const
{
  return fParamMaxima;
}
//____________________________________________________________________________
double valor::T2KPhysicalParameterLimits::GetNearestPhysicalValue(conv::PhysParamId_t p, double v) const
{
  if (this->IsWithinLimits(p,v)) {return v;}
  if (v < this->MinValue(p))     {return this->MinValue(p);}
  return this->MaxValue(p);
}
//____________________________________________________________________________
void valor::T2KPhysicalParameterLimits::SetDefaults(void)
{
  for (int ip = 0 ; ip < knppmax; ip++) {
    PhysParamId_t p = (PhysParamId_t)ip;
    fParamMinima[p] = -9999;
    fParamMaxima[p] =  9999;
    
    if (util::IsNormalisation(p)) {fParamMinima[p] = 0;}
    //if (util::IsEnergyScale  (p)) {fParamMinima[p] = 0;}
    if (util::IsMixingAngle  (p)) {fParamMinima[p] = 0; fParamMaxima[p] = 1;}
  }
  
  // Parameter limits from Asher's email (27th Jan 2015)
  
  // MA
  fParamMinima[conv::kNuisBANFFMAQE    ] = 0;
  fParamMinima[conv::kRHCNuisBANFFMAQE ] = 0;
  
  fParamMinima[conv::kNuisBANFFMARES   ] = 0;
  fParamMinima[conv::kRHCNuisBANFFMARES] = 0;
  
  // PFO
  fParamMinima[conv::kNuisBANFFFermiMom    ] = 0.89;
  fParamMaxima[conv::kNuisBANFFFermiMom    ] = 1.22;
  
  fParamMinima[conv::kRHCNuisBANFFFermiMom ] = 0.89;
  fParamMaxima[conv::kRHCNuisBANFFFermiMom ] = 1.22;
  
  // isospin
  fParamMinima[conv::kNuisBANFFBgRES    ] = 0.;
  fParamMinima[conv::kRHCNuisBANFFBgRES ] = 0.;
  
  // ca5
  fParamMinima[conv::kNuisBANFFCA5    ] = 0.;
  fParamMinima[conv::kRHCNuisBANFFCA5 ] = 0.;

  // BeRPA
  fParamMinima[conv::kNuisBANFFBeRPA_A] = 0.;
  fParamMinima[conv::kRHCNuisBANFFBeRPA_A] = 0.;
  fParamMinima[conv::kNuisBANFFBeRPA_B] = 0.;
  fParamMinima[conv::kRHCNuisBANFFBeRPA_B] = 0.;
  fParamMinima[conv::kNuisBANFFBeRPA_D] = 0.;
  fParamMinima[conv::kRHCNuisBANFFBeRPA_D] = 0.;
  fParamMinima[conv::kNuisBANFFBeRPA_E] = 0.;
  fParamMinima[conv::kRHCNuisBANFFBeRPA_E] = 0.;
  fParamMinima[conv::kNuisBANFFBeRPA_U] = 0.;
  fParamMinima[conv::kRHCNuisBANFFBeRPA_U] = 0.;

  // Shape 2p2h
  fParamMinima[conv::kNuisBANFFShape2p2h] = -1.;
  fParamMinima[conv::kRHCNuisBANFFShape2p2h] = -1.;
  fParamMaxima[conv::kNuisBANFFShape2p2h] = 1.;
  fParamMaxima[conv::kRHCNuisBANFFShape2p2h] = 1.;

  // Eb
  fParamMinima[conv::kNuisBANFFEb] = -3.;
  fParamMinima[conv::kRHCNuisBANFFEb] = -3.;
  fParamMaxima[conv::kNuisBANFFEb] = 10.;
  fParamMaxima[conv::kRHCNuisBANFFEb] = 10.;

#ifdef ADDBANFF00
  fParamMinima[conv::kNuisBANFFNew00] = 0.;
#endif //ADDBANFF00
#ifdef ADDBANFF01
  fParamMinima[conv::kNuisBANFFNew01] = 0.;
#endif //ADDBANFF01
#ifdef ADDBANFF02
  fParamMinima[conv::kNuisBANFFNew02] = 0.;
#endif //ADDBANFF02
#ifdef ADDBANFF03
  fParamMinima[conv::kNuisBANFFNew03] = 0.;
#endif //ADDBANFF03
#ifdef ADDBANFF04
  fParamMinima[conv::kNuisBANFFNew04] = 0.;
#endif //ADDBANFF04
#ifdef ADDBANFF05
  fParamMinima[conv::kNuisBANFFNew05] = 0.;
#endif //ADDBANFF05
#ifdef ADDBANFF06
  fParamMinima[conv::kNuisBANFFNew06] = 0.;
#endif //ADDBANFF06
#ifdef ADDBANFF07
  fParamMinima[conv::kNuisBANFFNew07] = 0.;
#endif //ADDBANFF07
#ifdef ADDBANFF08
  fParamMinima[conv::kNuisBANFFNew08] = 0.;
#endif //ADDBANFF08
#ifdef ADDBANFF09
  fParamMinima[conv::kNuisBANFFNew09] = 0.;
#endif //ADDBANFF09

  return;
}
//____________________________________________________________________________
void valor::T2KPhysicalParameterLimits::SetLimits(FitParamList * fpl)
{
  for (int ip = 0 ; ip < knppmax; ip++) {
    PhysParamId_t p = (PhysParamId_t)ip;
    if (fpl->Exists(p)){
      map<int,double>::const_iterator it_min = fParamMinima.find(p);
      map<int,double>::const_iterator it_max = fParamMaxima.find(p);
      if (util::IsNormalisation(p)) {
        // FPL set max valuex
        LOG("T2KPhysicalParameterLimits",pDEBUG) << "Setting maxima and minima in FitParamList for parameter: \n" << util::AsString(p);
        if (fpl->MinValue(p) < 0)
        fpl->SetMinValue(p,it_min->second);
        // Rely on the upper limits being set correctly
        // fpl->SetMaxValue(p,it_max->second);
      }
    }
    //    if (util::IsMixingAngle  (p)) {}
  }

}
