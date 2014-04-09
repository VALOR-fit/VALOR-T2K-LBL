//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <TMath.h>

#include "UtlOsc.h"

#include "ParameterValues.h"
#include "Conventions.h"
#include "UtlNumStat.h"
#include "Opt.h"

using namespace t2kana;
using namespace t2kana::conv;

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setfill;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::ofstream;
using std::ostringstream;

double t2kana::util::NuMuSurvProbOsc2Flavour(
  double mixing23, double dm23sq, double Et, bool unphysical_osc_prob = false)
{
  if(Et <= 0.) return 0.;

  double sinsq2theta23 = 0.;

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
   case(kNuMixSinsq2Theta):
       // mixing23 parameter is sinsq2theta23
       sinsq2theta23 = mixing23;
       break;
   case(kNuMixSinsqTheta):
       // mixing23 parameter is sinsqtheta23
       // convert to sinsq2theta23
       sinsq2theta23 = 4*mixing23*(1-mixing23);
       break;
   default:
       return 0.;
  }

  double osc   = TMath::Sin(1.267*dm23sq*kBaseline/Et);
  double posc  = sinsq2theta23 * TMath::Power(osc,2);
  double psurv = 1-posc;

  if(!unphysical_osc_prob) {
    psurv = TMath::Min(1.,psurv);
    psurv = TMath::Max(0.,psurv);
  }

  return psurv;
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbOsc2Flavour(
  double mixing23, double dm23sq, double Et, t2kana::conv::SuperkMode_t m,
  bool unphysical_osc_prob = false)
{
// Note: Using only the numu disappearance analysis modes
//
  if(m == kSKmNuMuNC1piC      || 
     m == kSKmNuMuNC1pi0      ||
     m == kSKmNuMuNCcoh       ||
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuNCoth_e     ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNC1pi0   ||
     m == kSKmNuMuBarNCcoh    ||
     m == kSKmNuMuBarNCoth_mu ||
     m == kSKmNuMuBarNCoth_e  ||
     m == kSKmNuECCQE         ||
     m == kSKmNuECC1pi        ||
     m == kSKmNuECCcoh        ||
     m == kSKmNuECCoth     ||
     //     m == kSKmNuECCoth      ||
     m == kSKmNuENC1piC       ||
     m == kSKmNuENC1pi0       ||
     m == kSKmNuENCcoh        ||
     m == kSKmNuENCoth_mu     ||
     m == kSKmNuENCoth_e      ||
     m == kSKmNuEBarCCQE      ||
     m == kSKmNuEBarCC1pi     ||
     m == kSKmNuEBarCCcoh     ||
     //     m == kSKmNuEBarCCoth  ||
     m == kSKmNuEBarCCoth   ||
     m == kSKmNuEBarNC1piC    ||
     m == kSKmNuEBarNC1pi0    ||
     m == kSKmNuEBarNCcoh     ||
     m == kSKmNuEBarNCoth_mu  ||
     m == kSKmNuEBarNCoth_e   )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     //     m == kSKmNuMuCCoth    ||
     m == kSKmNuMuCCoth     ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     //     m == kSKmNuMuBarCCoth || 
     m == kSKmNuMuBarCCoth  ||
     m == kSKmMEC             )
  {
    return NuMuSurvProbOsc2Flavour(mixing23, dm23sq, Et, unphysical_osc_prob);
  }
  else 
  if(
     m == kSKmOscNuMuCCQE     ||
     m == kSKmOscNuMuCC1pi    ||
     m == kSKmOscNuMuCCcoh    ||
     //     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuMuCCoth  ||
     m == kSKmOscNuMuBarCCQE     ||
     m == kSKmOscNuMuBarCC1pi    ||
     m == kSKmOscNuMuBarCCcoh    ||
     //     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuMuBarCCoth  ||
     m == kSKmOscNuEBarCCQE     ||
     m == kSKmOscNuEBarCC1pi    ||
     m == kSKmOscNuEBarCCcoh    ||
     //     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuEBarCCoth 
)
  {
    return 0;
  }
  else {
    cerr << "NuMuSurvProbOsc2D > Can't handle SK mode = " << m << endl;
    exit(1);
  }
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbDecay (
   double sinsq2theta, double alpha, double Et)
{
// Note: In oscillation formula (Dm^2 * L / 4E) in natural units becomes 
// 1.267 * Dm^2 * L / E where Dm^2 in eV2, L in km, E in GeV. 
// Apply the same conversion in the decoherence survival prob formula (Eq. 13) 
// given in Barger et al., astro-ph/9810121v1.
//
// a L / 2E (natural units) --> 2 * 1.267 a (eV^2) L(km) /E (GeV)
//

  double cossq2theta = TMath::Sqrt( TMath::Max(0.,1.-sinsq2theta) );
  double cos2theta   = TMath::Sqrt(cossq2theta);
  double sinsqtheta  = (1.-cos2theta)/2.;
  double cossqtheta  = (1.+cos2theta)/2.;
  double expaL       = TMath::Exp(-2*1.267*alpha*kBaseline/Et);

  double psurv = TMath::Power(sinsqtheta + cossqtheta*expaL, 2.);
  return psurv;
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbDecay (
  double sinsq2theta, double alpha, double Et, 
  t2kana::conv::SuperkMode_t m)
{
// Note: Using only the numu disappearance analysis modes
//
  if(m == kSKmNuECCQE        ||
     m == kSKmNuECC1pi       ||
     m == kSKmNuECCcoh       ||
     m == kSKmNuECCoth    ||
     m == kSKmNuENC1piC      ||
     m == kSKmNuENCoth_mu    ||
     m == kSKmNuEBarCCQE     ||
     m == kSKmNuEBarCC1pi    ||
     m == kSKmNuEBarCCcoh    ||
     m == kSKmNuEBarCCoth ||
     m == kSKmNuEBarNC1piC   ||
     m == kSKmNuEBarNCoth_mu )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     m == kSKmNuMuCCoth    ||
     m == kSKmNuMuNC1piC      ||
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     m == kSKmNuMuBarCCoth ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNCoth_mu ||
     m == kSKmMEC             )
  {
    return NuMuSurvProbDecay(sinsq2theta, alpha, Et);
  }
  else 
  if(
m == kSKmOscNuMuCCQE     ||
     m == kSKmOscNuMuCC1pi    ||
     m == kSKmOscNuMuCCcoh    ||
//     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuMuCCoth  ||
     m == kSKmOscNuMuBarCCQE     ||
     m == kSKmOscNuMuBarCC1pi    ||
     m == kSKmOscNuMuBarCCcoh    ||
//     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuMuBarCCoth  ||
     m == kSKmOscNuEBarCCQE     ||
     m == kSKmOscNuEBarCC1pi    ||
     m == kSKmOscNuEBarCCcoh    ||
//     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuEBarCCoth
     )
  {
    return 0;
  }
  else {
    cerr << "NuMuSurvProbDecay > Can't handle SK mode = " << m << endl;
    exit(1);
  }
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbDecoh (
   double sinsq2theta, double musq, double Et)
{
//
// Note: In oscillation formula (Dm^2 * L / 4E) in natural units becomes 
// 1.267 * Dm^2 * L / E where Dm^2 in eV2, L in km, E in GeV. 
// Apply the same conversion in the decoherence survival prob formula (Eq. 5) 
// given in Fogli et al., PRD67,093006(2003).
//
// \mu^2 L / 2 E (natural units) --> 2 * 1.267 \mu^2 (eV^2) L(km) /E (GeV)
//

  double expmusqLE = TMath::Exp(-2*1.267*musq*kBaseline/Et);
  double psurv = 1. - 0.5*sinsq2theta * ( 1 - expmusqLE );
  return psurv;
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbDecoh (
  double sinsq2theta, double musq, double Et, 
  t2kana::conv::SuperkMode_t m)
{
// Note: Using only the numu disappearance analysis modes
//
  if(m == kSKmNuECCQE        ||
     m == kSKmNuECC1pi       ||
     m == kSKmNuECCcoh       ||
     m == kSKmNuECCoth    ||
     m == kSKmNuENC1piC      ||
     m == kSKmNuENCoth_mu    ||
     m == kSKmNuEBarCCQE     ||
     m == kSKmNuEBarCC1pi    ||
     m == kSKmNuEBarCCcoh    ||
     m == kSKmNuEBarCCoth ||
     m == kSKmNuEBarNC1piC   ||
     m == kSKmNuEBarNCoth_mu )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     m == kSKmNuMuCCoth    ||
     m == kSKmNuMuNC1piC      ||
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     m == kSKmNuMuBarCCoth ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNCoth_mu ||
     m == kSKmMEC             )
   {
    return NuMuSurvProbDecoh(sinsq2theta, musq, Et);
  }
  else 
  if(
     m == kSKmOscNuMuCCQE     ||
     m == kSKmOscNuMuCC1pi    ||
     m == kSKmOscNuMuCCcoh    ||
//     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuMuCCoth  ||
     m == kSKmOscNuMuBarCCQE     ||
     m == kSKmOscNuMuBarCC1pi    ||
     m == kSKmOscNuMuBarCCcoh    ||
//     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuMuBarCCoth  ||
     m == kSKmOscNuEBarCCQE     ||
     m == kSKmOscNuEBarCC1pi    ||
     m == kSKmOscNuEBarCCcoh    ||
//     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuEBarCCoth
     )
  {
    return 0;
  }
  else {
    cerr << "NuMuSurvProbDecoh > Can't handle SK mode = " << m << endl;
    exit(1);
  }
}
//____________________________________________________________________________
