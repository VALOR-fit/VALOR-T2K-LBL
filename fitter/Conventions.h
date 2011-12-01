//____________________________________________________________________________
/*!

\class    t2kana::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since v1.0.0:

 @ 01-20/07/11 - CA
   Added nuisance and 12-, 13-sector osc params in PhysParamId_t.
   Added PhysHypothesis_t, OscChannelId_t and NuMassHierarchy_t.
   Added some constants (GF,hbar,...) used in calculations and added values 
   and uncertainties of oscillation params from a global analysis.
 @ 28/09/11 - CA
   Added kMassHierarchyAssumed constant to specify the mass hierarchy to be
   used in the t2kana library build. Used separate PhysParamId_t for the
   decoherence and decay sinsq2theta parameter.

*/
//____________________________________________________________________________

#ifndef _t2kana_conv_h_
#define _t2kana_conv_h_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace t2kana {
namespace conv   {

// confidence level
typedef enum EConfLevel {
   kCLUndefined = -1, 
   kCL68, 
   kCL90, 
   kCL95, 
   kCL99
} ConfLevel_t; 

// enumeration of event samples 
static const int kNSuperkSamples  = 2;  // number of samples included in analysis
static const int kns = kNSuperkSamples;
typedef enum ESuperkSample {
   kSKsUndefined = -1, 
   kSKsFCFV1RingMuLike, // fully-contained fiducial-volume single-ring mu-like
   kSKsFCFV1RingELike,  // fully-contained fiducial-volume single-ring e-like
} SuperkSample_t; 

// enumeration of (true) reaction modes that may be considered separately within each SuperK sample
static const int kNSuperkModes = 30; // number of reaction modes considered separately within each SuperK sample
static const int knm = kNSuperkModes; 
typedef enum ESuperkMode {
   kSKmUndefined = -1, 
   kSKmNuMuCCQE,         // numu CC QE
   kSKmNuMuCC1pi,        // numu CC 1pi
   kSKmNuMuCCcoh,        // numu CC coherent
   kSKmNuMuCCoth_mu,     // numu CC other, nu_mu disappearance analysis definition
   kSKmNuMuCCoth_e,      // numu CC other, nu_e appearance analysis definition
   kSKmNuMuNC,           // numu NC (all)
   kSKmNuMuNC1pi0,       // numu NC 1pi0
   kSKmNuMuNCcoh,        // numu NC coherent
   kSKmNuMuNCoth_e,      // numu NC other, nu_e appearance analysis definition
   kSKmNuMuBarCCQE,      // numubar CCQE
   kSKmNuMuBarCC1pi,     // numubar CC 1pi
   kSKmNuMuBarCCcoh,     // numubar CC coh
   kSKmNuMuBarCCoth_mu,  // numubar CC other, nu_mu disappearance analysis definition
   kSKmNuMuBarCCoth_e,   // numubar CC other, nu_e appearance analysis definition
   kSKmNuMuBarNC,        // numubar NC (all)
   kSKmNuMuBarNC1pi0,    // numubar NC 1pi0
   kSKmNuMuBarNCcoh,     // numubar NC coherent
   kSKmNuMuBarNCoth_e,   // numubar NC other, nu_e appearance analysis definition
   kSKmNuECCQE,          // nue CC QE
   kSKmNuECC1pi,         // nue CC 1pi
   kSKmNuECCcoh,         // nue CC coherent
   kSKmNuECCoth_mu,      // nue CC other, nu_mu disappearance analysis definition
   kSKmNuECCoth_e,       // nue CC other, nu_e appearance analysis definition
   kSKmNuENC,            // nue NC (all)
   kSKmNuENC1pi0,        // nue NC 1pi0
   kSKmNuENCcoh,         // nue NC coherent
   kSKmNuENCoth_e,       // nue NC other, nu_e appearance analysis definition
   kSKmSigNuECCQE,       // sig nue
   kSKmSigNuECC1pi,      // sig nue
   kSKmSigNuECCoth       // sig nue
} SuperkMode_t; 

//
// Define which SK samples are considered
//

static const SuperkSample_t kSKSamples[kns] = {
    kSKsFCFV1RingMuLike,
    kSKsFCFV1RingELike
};

//
// Define which SK modes are considered for each SK sample included in the analysis
//

static const SuperkMode_t kSKReactionModes[kns][knm] = {

  // > 1 muon-like ring (numu disappearance)
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi, 
    kSKmNuMuCCoth_mu,                      
    kSKmNuMuNC,          
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarNC,       
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCoth_mu,                                
    kSKmNuENC,                                                               
    kSKmSigNuECCQE, 
    kSKmSigNuECC1pi, 
    kSKmSigNuECCoth,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  },

  // > 1 electron-like ring (nue appearance)
  {
    kSKmNuMuCCQE, 
    kSKmNuMuCC1pi, 
    kSKmNuMuCCcoh, 
    kSKmNuMuCCoth_e, 
    kSKmNuMuNC1pi0, 
    kSKmNuMuNCcoh, 
    kSKmNuMuNCoth_e,                            
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    kSKmNuMuBarCCoth_e,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    kSKmNuECCoth_e,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCoth,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  }
};

// (true) reaction modes considered separately in nd280 inputs
static const int kNND280Modes = 21; // number of reaction modes considered separately within the ND280 sample
static const int knmnd = kNND280Modes; 
typedef enum END280Mode {
   kNDmUndefined = -1, 
   kNDmNuMuCCQE, 
   kNDmNuMuCCRES, 
   kNDmNuMuCCCOH, 
   kNDmNuMuCCOth, 
   kNDmNuMuNCRES, 
   kNDmNuMuNCCOH, 
   kNDmNuMuNCOth, 
   kNDmNuMuBarCCQE, 
   kNDmNuMuBarCCRES, 
   kNDmNuMuBarCCCOH, 
   kNDmNuMuBarCCOth, 
   kNDmNuMuBarNCRES, 
   kNDmNuMuBarNCCOH, 
   kNDmNuMuBarNCOth, 
   kNDmNuECCQE, 
   kNDmNuECCRES, 
   kNDmNuECCCOH, 
   kNDmNuECCOth, 
   kNDmNuENCRES, 
   kNDmNuENCCOH, 
   kNDmNuENCOth
} ND280Mode_t; 

// fit parameters
typedef enum EPhysParamId {   
   kPhypUndefined = -1, 
   kPhypSinsq2Theta23,       // mixing angle sin^2(2*\theta_{23}) 
   kPhypSinsq2Theta13,       // mixing angle sin^2(2*\theta_{13}) 
   kPhypSinsq2Theta12,       // mixing angle sin^2(2*\theta_{12}) 
   kPhypDeltaCP,             // CP-invariance violating phase 
   kPhypDm23sq,              // modulus of squared-mass splitting \Deltam_{23}^{2} (always positive independent of mass hierarchy) 
   kPhypDm12sq,              // squared-mass splitting \Deltam_{12}^{2} 
   kPhypAlphaDecay,          //
   kPhypSinsq2ThetaDecay,    // 
   kPhypMusqDecoh,           //
   kPhypSinsq2ThetaDecoh,    // 
   kNuisSuperkEffCCQEoth,    //
   kNuisSuperkEffCCQErc0u,   //
   kNuisSuperkEffCCQErc1u,   //
   kNuisSuperkEffCCQErc2u,   //
   kNuisSuperkEffCCnonQE,    //
   kNuisSuperkEffNC,         //
   kNuisSuperkEffNueCC,      //
   kNuisSuperkEnergyScale,   //
   kNuisNd280Eff,            //
   kNuisFluxNorm,            //
   kNuisFluxShapeNuMu0u,     //
   kNuisFluxShapeNuMu1u,     //
   kNuisFluxShapeNuMu2u,     //
   kNuisFluxShapeNuMu3u,     //
   kNuisFluxShapeNuMu4u,     //
   kNuisFluxShapeNuMu5u,     //
   kNuisFluxShapeNuMu6u,     //
   kNuisFluxShapeNuMu7u,     //
   kNuisFluxShapeNuMu8u,     //
   kNuisFluxShapeNuMu9u,     //
   kNuisFluxShapeNuMu10u,    //
   kNuisFluxShapeNuMu11u,    //
   kNuisFluxShapeNuMu12u,    //
   kNuisFluxShapeNuMu13u,    //
   kNuisFluxShapeNuMu14u,    //
   kNuisFluxShapeNuMu15u,    //
   kNuisFluxShapeNuMu16u,    //
   kNuisFluxShapeNuMu17u,    //
   kNuisFluxShapeNuMu18u,    //
   kNuisXSecCCQE,            //  
   kNuisXSecCC1piToCCQE,     //
   kNuisXSecCCothToCCQE,     //
   kNuisXSecNCToCCQE,        //
   kNuisXSecNuclEffect       //
} PhysParamId_t;

// physics hypotheses
typedef enum EPhysHypothesis {   
   kHypUndefined = -1,
   kHypDecay,             // neutrino decay 
   kHypDecoh,             // neutrino decoherence 
   kHypNuMuDisapp2F,      // muon-neutrino disappearance, simple 2-flavour oscillation framework 
   kHypOsc3F              // neutrino oscillations, 3-flavour framework 
} PhysHypothesis_t;

// oscillation channels
typedef enum EOscChannelId {   
   kOscUndefined = -1, 
   kOscNuE2NuE,
   kOscNuE2NuMu,
   kOscNuE2NuTau,
   kOscNuMu2NuE,
   kOscNuMu2NuMu,
   kOscNuMu2NuTau,
   kOscNuTau2NuE,
   kOscNuTau2NuMu,
   kOscNuTau2NuTau,
   kOscNuEBar2NuEBar,
   kOscNuEBar2NuMuBar,
   kOscNuEBar2NuTauBar,
   kOscNuMuBar2NuEBar,
   kOscNuMuBar2NuMuBar,
   kOscNuMuBar2NuTauBar,
   kOscNuTauBar2NuEBar,
   kOscNuTauBar2NuMuBar,
   kOscNuTauBar2NuTauBar
} OscChannelId_t;

// neutrino mass hierarchy
typedef enum ENuMassHierarchy {   
   kNmhUndefined = -1, 
   kNmhNormal,
   kNmhInverted
} NuMassHierarchy_t;

static const NuMassHierarchy_t kMassHierarchyAssumed = kNmhNormal; // mass hierarchy used in T2KAna build

static const int kPdgNuE      =  12;
static const int kPdgNuMu     =  14;
static const int kPdgNuTau    =  16;
static const int kPdgNuEBar   = -12;
static const int kPdgNuMuBar  = -14;
static const int kPdgNuTauBar = -16;

// Some constants
static const double kASmallNum      = 1E-10;                         // 
static const double kBaseline       = 295.335;                       // Tokai->Kamioka (km)
static const double kOscDimensFactor= 1.26693281;                    // Dimensional factor (GeV/eV^2/km) appearing in netrino oscillation formula
static const double kCrustDensity   = 2.7;                           // Average density of Earth's continental crust (gr / cm^3)
static const double kpi             = 3.141592653589793;             // pi
static const double khbar_SI        = 6.582119e-16;                  // hbar (eV s)
static const double kc_SI           = 2.9979246e10;                  // c (cm/s)
static const double khbarc_SI       = khbar_SI*kc_SI;                // hbar*c (eV * cm)
static const double khbarc3_SI      = khbarc_SI*khbarc_SI*khbarc_SI; // (hbar*c)^3 (eV^3 * cm^3)
static const double kGF             = 1.16637e-5;                    // Fermi coupling constant, (GeV^-2, natural units)
static const double kGF_SI          = kGF*khbarc3_SI*1.0e-18;        // Fermi coupling constant, (eV cm^3, SI units), 1e-18 is used for GeV^-2 -> eV^-2
static const double kNA             = 6.022142e23;                   // Avogadro number (mol^-1)

// Oscillation parameters, recent global fit values & 1 sigma errors
static const double kGlobFitSinsq2Theta23    = 1.00;
static const double kGlobFitSinsq2Theta13    = 0.10;
static const double kGlobFitSinsq2Theta12    = 0.86;
static const double kGlobFitDm23sq           = 2.35E-3;
static const double kGlobFitDm12sq           = 7.58E-5;
static const double kGlobFitErrSinsq2Theta23 = 0.15; // ~15%
static const double kGlobFitErrSinsq2Theta12 = 0.06; // ~6%
static const double kGlobFitErrDm23sq        = 0.05; // ~5%
static const double kGlobFitErrDm12sq        = 0.03; // ~3%

//
// ------------------------------------------------------------------------
// Binning
// ------------------------------------------------------------------------
//

//
// Variable Etrue binning scheme used by the Beam WG.
// Used in certain Beam WG inputs (eg numu flux covariance matrix) 
//
static const int kNEtrueBinsBeamWG = 19;
static const double kEtrueBinEdgesBeamWG[kNEtrueBinsBeamWG+1] = 
{
   0.000, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800,  1.000, 
   1.200, 1.500, 2.000, 2.500, 3.000, 3.500, 4.000, 5.000, 7.000, 10.000
};

//
// Binning for SK MC templates
//

static const double kEtrueMin    = 0;
static const double kEtrueMax    = 12;
static const double kEtrueBinSz  = 0.05;
static const int    kNEtrueBins  = (int) ((kEtrueMax - kEtrueMin)/kEtrueBinSz);

//
// Ereco binning for SK 1Rmu sample
//

// 200 50-MeV bins from 0 to 10 GeV
static const int kNErecoBins1Rmu = 200;
static const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
   1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950, 
   2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450, 
   2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950, 
   3.000,  3.050,  3.100,  3.150,  3.200,  3.250,  3.300,  3.350,  3.400,  3.450, 
   3.500,  3.550,  3.600,  3.650,  3.700,  3.750,  3.800,  3.850,  3.900,  3.950, 
   4.000,  4.050,  4.100,  4.150,  4.200,  4.250,  4.300,  4.350,  4.400,  4.450, 
   4.500,  4.550,  4.600,  4.650,  4.700,  4.750,  4.800,  4.850,  4.900,  4.950, 
   5.000,  5.050,  5.100,  5.150,  5.200,  5.250,  5.300,  5.350,  5.400,  5.450, 
   5.500,  5.550,  5.600,  5.650,  5.700,  5.750,  5.800,  5.850,  5.900,  5.950, 
   6.000,  6.050,  6.100,  6.150,  6.200,  6.250,  6.300,  6.350,  6.400,  6.450, 
   6.500,  6.550,  6.600,  6.650,  6.700,  6.750,  6.800,  6.850,  6.900,  6.950, 
   7.000,  7.050,  7.100,  7.150,  7.200,  7.250,  7.300,  7.350,  7.400,  7.450, 
   7.500,  7.550,  7.600,  7.650,  7.700,  7.750,  7.800,  7.850,  7.900,  7.950, 
   8.000,  8.050,  8.100,  8.150,  8.200,  8.250,  8.300,  8.350,  8.400,  8.450, 
   8.500,  8.550,  8.600,  8.650,  8.700,  8.750,  8.800,  8.850,  8.900,  8.950, 
   9.000,  9.050,  9.100,  9.150,  9.200,  9.250,  9.300,  9.350,  9.400,  9.450, 
   9.500,  9.550,  9.600,  9.650,  9.700,  9.750,  9.800,  9.850,  9.900,  9.950, 
  10.000
};
/*
// Fewer bins at higher energies for improved speed
// -   50 MeV bins in 0-3 GeV
// -  200 MeV bins in 3-4 GeV
// -  500 MeV bins in 4-6 GeV
// - 1000 MeV bins > 6 GeV
static const int kNErecoBins1Rmu = 73;
static const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
   1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950, 
   2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450, 
   2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950, 
   3.000,  3.200,  3.400,  3.600,  3.800,
   4.000,  4.500,  5.000,  5.500,  
   6.000,  7.000,  8.000,  9.000, 10.000
};
*/

//
// Ereco binning for SK 1Re sample
//

// 1 bin from 0 to 10 GeV (counting expt)
static const int kNErecoBins1Re = 1;
static const double kErecoBinEdges1Re[kNErecoBins1Re+1] = 
{ 
   0., 10. 
};
/*
// Same binning as the standard 1Rmu one: 200 50-MeV bins from 0 to 10 GeV
static const int kNErecoBins1Re = 200;
static const double kErecoBinEdges1Re[kNErecoBins1Re+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
   1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950, 
   2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450, 
   2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950, 
   3.000,  3.050,  3.100,  3.150,  3.200,  3.250,  3.300,  3.350,  3.400,  3.450, 
   3.500,  3.550,  3.600,  3.650,  3.700,  3.750,  3.800,  3.850,  3.900,  3.950, 
   4.000,  4.050,  4.100,  4.150,  4.200,  4.250,  4.300,  4.350,  4.400,  4.450, 
   4.500,  4.550,  4.600,  4.650,  4.700,  4.750,  4.800,  4.850,  4.900,  4.950, 
   5.000,  5.050,  5.100,  5.150,  5.200,  5.250,  5.300,  5.350,  5.400,  5.450, 
   5.500,  5.550,  5.600,  5.650,  5.700,  5.750,  5.800,  5.850,  5.900,  5.950, 
   6.000,  6.050,  6.100,  6.150,  6.200,  6.250,  6.300,  6.350,  6.400,  6.450, 
   6.500,  6.550,  6.600,  6.650,  6.700,  6.750,  6.800,  6.850,  6.900,  6.950, 
   7.000,  7.050,  7.100,  7.150,  7.200,  7.250,  7.300,  7.350,  7.400,  7.450, 
   7.500,  7.550,  7.600,  7.650,  7.700,  7.750,  7.800,  7.850,  7.900,  7.950, 
   8.000,  8.050,  8.100,  8.150,  8.200,  8.250,  8.300,  8.350,  8.400,  8.450, 
   8.500,  8.550,  8.600,  8.650,  8.700,  8.750,  8.800,  8.850,  8.900,  8.950, 
   9.000,  9.050,  9.100,  9.150,  9.200,  9.250,  9.300,  9.350,  9.400,  9.450, 
   9.500,  9.550,  9.600,  9.650,  9.700,  9.750,  9.800,  9.850,  9.900,  9.950, 
  10.000
};
*/


} // conv
} // t2kana
  
#endif 


