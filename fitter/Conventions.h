//____________________________________________________________________________
/*!

\class    t2kana::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since 2011_v0.9.0:

*/
//____________________________________________________________________________

#ifndef _t2kana_conv_h_
#define _t2kana_conv_h_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace t2kana {
namespace conv   {

//
// ------------------------------------------------------------------------
// Constants & commonly used codes
// ------------------------------------------------------------------------
//

// PDG Codes
static const int kPdgNuE      =  12;
static const int kPdgNuMu     =  14;
static const int kPdgNuTau    =  16;
static const int kPdgNuEBar   = -12;
static const int kPdgNuMuBar  = -14;
static const int kPdgNuTauBar = -16;

//
// Oscillation parameters, recent global fit values & 1 sigma errors and
// typical initialization and range values for fitter
//

// Dm^2 global fit results
static const double kGlobFitDm23sq           = 2.35E-3;
static const double kGlobFitDm12sq           = 7.58E-5;
static const double kGlobFitErrDm23sq        = 0.05; // ~5%
static const double kGlobFitErrDm12sq        = 0.03; // ~3%
// deltacp global fit results
static const double kGlobFitDeltaCP          = 0.0;
static const double kGlobFitErrDeltaCP       = 0.10; // ~10%
// sinsq2theta global fit results
static const double kGlobFitSinsq2Theta23    = 1.00;
static const double kGlobFitSinsq2Theta13    = 0.10;
static const double kGlobFitSinsq2Theta12    = 0.86;
static const double kGlobFitErrSinsq2Theta23 = 0.15; // ~15%
static const double kGlobFitErrSinsq2Theta12 = 0.06; // ~6%
static const double kGlobFitErrSinsq2Theta13 = 0.10; // ~10%
// sinsqtheta global fit results
static const double kGlobFitSinsqTheta23    = 0.420;
static const double kGlobFitSinsqTheta13    = 0.025;
static const double kGlobFitSinsqTheta12    = 0.312;
static const double kGlobFitErrSinsqTheta23 = 0.261;  // ~26% (1 sigma range)
static const double kGlobFitErrSinsqTheta12 = 0.110;  // ~11% (1 sigma range)
static const double kGlobFitErrSinsqTheta13 = 0.560;  // ~56% (1 sigma range)

// Dm^2 fitter init values and range
static const double kInitFitDm23sq           = kGlobFitDm23sq;
static const double kMinFitDm23sq            = 1.0E-3;
static const double kMaxFitDm23sq            = 4.0E-3;
static const double kInitFitDm12sq           = kGlobFitDm12sq;
static const double kMinFitDm12sq            = 6.0E-5;
static const double kMaxFitDm12sq            = 9.0E-5;
// sinsq2theta fitter init values and range
static const double kInitFitSinsq2Theta23    = 0.90;
static const double kMinFitSinsq2Theta23     = 0.60;
static const double kMaxFitSinsq2Theta23     = 1.00;
static const double kInitFitSinsq2Theta13    = kGlobFitSinsq2Theta13;
static const double kMinFitSinsq2Theta13     = 0.00;
static const double kMaxFitSinsq2Theta13     = 0.40;
static const double kInitFitSinsq2Theta12    = kGlobFitSinsq2Theta12;
static const double kMinFitSinsq2Theta12     = 0.70;
static const double kMaxFitSinsq2Theta12     = 1.00;
// sinsqtheta fitter init values and range
static const double kInitFitSinsqTheta23     = 0.90;
static const double kMinFitSinsqTheta23      = 0.60;
static const double kMaxFitSinsqTheta23      = 1.00;
static const double kInitFitSinsqTheta13     = kGlobFitSinsqTheta13;
static const double kMinFitSinsqTheta13      = 0.00;
static const double kMaxFitSinsqTheta13      = 0.40;
static const double kInitFitSinsqTheta12     = kGlobFitSinsqTheta12;
static const double kMinFitSinsqTheta12      = 0.70;
static const double kMaxFitSinsqTheta12      = 1.00;

//
// Other constants
//
static const double kMaxAllowedWght = 1E+3;                          // Complain if any computed weights exceeds this value
static const double kASmallNum      = 1E-10;                         // 
static const double kOneMeV         = 0.001;                         // 
static const double kBaseline       = 295.335;                       // Tokai->Kamioka (km)
static const double kOscDimensFactor= 1.26693281;                    // Dimensional factor (GeV/eV^2/km) appearing in netrino oscillation formula
static const double kCrustDensity   = 2.6;                           // Average density of Earth's continental crust (gr / cm^3)
static const double kpi             = 3.141592653589793;             // pi
static const double khbar_SI        = 6.582119e-16;                  // hbar (eV s)
static const double kc_SI           = 2.9979246e10;                  // c (cm/s)
static const double khbarc_SI       = khbar_SI*kc_SI;                // hbar*c (eV * cm)
static const double khbarc3_SI      = khbarc_SI*khbarc_SI*khbarc_SI; // (hbar*c)^3 (eV^3 * cm^3)
static const double kGF             = 1.16637e-5;                    // Fermi coupling constant, (GeV^-2, natural units)
static const double kGF_SI          = kGF*khbarc3_SI*1.0e-18;        // Fermi coupling constant, (eV cm^3, SI units), 1e-18 is used for GeV^-2 -> eV^-2
static const double kNA             = 6.022142e23;                   // Avogadro number (mol^-1)

//
// ------------------------------------------------------------------------
// Enumerations
// ------------------------------------------------------------------------
//

// Indices for array with default graphic attributes (see below)
static const int kNGraphicAttributes = 5; 
typedef enum EGraphicAttributeIdx
{   
   kGrAttUndefined  = -1, 
   kGrAttFillStyle  =  0, 
   kGrAttFillColour =  1, 
   kGrAttLineStyle  =  2, 
   kGrAttLineWidth  =  3, 
   kGrAttLineColour =  4 
} 
GraphicAttributeIdx_t;

// Confidence level
typedef enum EConfLevel {
   kCLUndefined = -1, 
   kCL68, 
   kCL90, 
   kCL95, 
   kCL99
} 
ConfLevel_t; 

// Physics hypotheses
typedef enum EPhysHypothesis 
{   
   kHypUndefined = -1,
   kHypDecay,             // neutrino decay 
   kHypDecoh,             // neutrino decoherence 
   kHypNuMuDisapp2F,      // muon-neutrino disappearance, simple 2-flavour oscillation framework 
   kHypOsc3F              // neutrino oscillations, 3-flavour framework 
} 
PhysHypothesis_t;

// Combined or single analysis 
typedef enum EOscAnaType 
{   
   kAnaUndefined = -1, 
   kAnaCombined,
   kAnaNuMuDisappOnly,
   kAnaNuEAppOnly
} 
OscAnaType_t;

// Neutrino mass hierarchy
typedef enum ENuMassHierarchy 
{   
   kNmhUndefined = -1, 
   kNmhNormal,
   kNmhInverted
} 
NuMassHierarchy_t;

// Choices of oscillation mixing parameter
typedef enum ENuOscMixParam
{   
   kNuMixUndefined = -1, 
   kNuMixSinsq2Theta,      // sin^2(2*theta)
   kNuMixSinsqTheta        // sin^2(theta), preferred choice in globa nue appearance + numu disappearance fit
} 
NuOscMixParam_t;

// Oscillation channels
typedef enum EOscChannelId 
{   
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
} 
OscChannelId_t;

// Enumeration of event samples that may be included in the analysis
// Actual choice for current analysis is given in `Configuration' section.
static const int kNSuperkSamples  = 2;  // number of samples 
static const int kns = kNSuperkSamples; // as above - shorter name
typedef enum ESuperkSample 
{
   kSKsUndefined = -1, 
   kSKsFCFV1RingMuLike, // fully-contained fiducial-volume single-ring mu-like
   kSKsFCFV1RingELike,  // fully-contained fiducial-volume single-ring e-like
} 
SuperkSample_t; 

// Enumeration of (true) reaction modes that may be considered separately within each SuperK sample
// Actual choice for current analysis is given in `Configuration' section.
static const int kNSuperkModes = 45;  // number of reaction modes kept track of within each sample
static const int knm = kNSuperkModes; // as above - shorter name
typedef enum ESuperkMode 
{
   kSKmUndefined = -1, 
   kSKmNuMuCCQE,         // numu CC QE
   kSKmNuMuCC1pi,        // numu CC 1pi
   kSKmNuMuCCcoh,        // numu CC coherent
   kSKmNuMuCCoth_mu,     // numu CC other, nu_mu disappearance analysis definition
   kSKmNuMuCCoth_e,      // numu CC other, nu_e appearance analysis definition
   kSKmNuMuNC1piC,       // numu NC 1 charged pi
   kSKmNuMuNC1pi0,       // numu NC 1pi0
   kSKmNuMuNCcoh,        // numu NC coherent
   kSKmNuMuNCoth_mu,     // numu NC other, nu_mu appearance analysis definition
   kSKmNuMuNCoth_e,      // numu NC other, nu_e appearance analysis definition
   kSKmNuMuBarCCQE,      // numubar CCQE
   kSKmNuMuBarCC1pi,     // numubar CC 1pi
   kSKmNuMuBarCCcoh,     // numubar CC coh
   kSKmNuMuBarCCoth_mu,  // numubar CC other, nu_mu disappearance analysis definition
   kSKmNuMuBarCCoth_e,   // numubar CC other, nu_e appearance analysis definition
   kSKmNuMuBarNC1piC,    // numubar NC 1 charged pi
   kSKmNuMuBarNC1pi0,    // numubar NC 1pi0
   kSKmNuMuBarNCcoh,     // numubar NC coherent
   kSKmNuMuBarNCoth_mu,  // numubar NC other, nu_mu appearance analysis definition
   kSKmNuMuBarNCoth_e,   // numubar NC other, nu_e appearance analysis definition
   kSKmNuECCQE,          // nue CC QE
   kSKmNuECC1pi,         // nue CC 1pi
   kSKmNuECCcoh,         // nue CC coherent
   kSKmNuECCoth_mu,      // nue CC other, nu_mu disappearance analysis definition
   kSKmNuECCoth_e,       // nue CC other, nu_e appearance analysis definition
   kSKmNuENC1piC,        // nue NC 1 charged pi
   kSKmNuENC1pi0,        // nue NC 1pi0
   kSKmNuENCcoh,         // nue NC coherent
   kSKmNuENCoth_mu,      // nue NC other, nu_mu appearance analysis definition
   kSKmNuENCoth_e,       // nue NC other, nu_e appearance analysis definition
   kSKmNuEBarCCQE,       // nuebar CC QE
   kSKmNuEBarCC1pi,      // nuebar CC 1pi
   kSKmNuEBarCCcoh,      // nuebar CC coherent
   kSKmNuEBarCCoth_mu,   // nuebar CC other, nu_mu disappearance analysis definition
   kSKmNuEBarCCoth_e,    // nuebar CC other, nu_e appearance analysis definition
   kSKmNuEBarNC1piC,     // nuebar NC 1 charged pi
   kSKmNuEBarNC1pi0,     // nuebar NC 1pi0
   kSKmNuEBarNCcoh,      // nuebar NC coherent
   kSKmNuEBarNCoth_mu,   // nuebar NC other, nu_mu appearance analysis definition
   kSKmNuEBarNCoth_e,    // nuebar NC other, nu_e appearance analysis definition
   kSKmSigNuECCQE,       // sig nue CCQE
   kSKmSigNuECC1pi,      // sig nue CC 1pi
   kSKmSigNuECCcoh,      // sig nue coh
   kSKmSigNuECCoth_mu,   // sig nue other, nu_mu disappearance analysis definition
   kSKmSigNuECCoth_e     // sig nue other, nu_e appearance analysis definition
} 
SuperkMode_t; 

// Type of SK MC template
typedef enum ESKMCTemplateType
{   
   kSKMCtUndefined = -1, 
   kSKMCtOri,   // nominal MC
   kSKMCtNom,   // nominal MC + tuning (ext) weights 
   kSKMCtTun,   // nominal MC + tuning (ext) weights + tuning (T2K) weights
   kSKMCtSys,   // nominal MC + tuning (ext) weights + tuning (T2K) weights + systematic variations
   kSKMCtFin,   // nominal MC + tuning (ext) weights + tuning (T2K) weights + systematic variations + physics hypothesis 
} 
SKMCTemplateType_t;

// Enumeration of flux histograms types held by NuFlux
typedef enum EFluxType
{   
   kFlxUndefined      = -1, 
   kFlxNominal        =  0,  // flux used for MC production
   kFlxExtTuning      =  1,  // flux used for MC production + tuning wghts based on external data (NA61,...)
   kFlxFinalTuning    =  2,  // flux used for MC production + tuning wghts based on external data (NA61,...) + tuning wghts based on T2K data (ND280 data fit by BANFF,...)
   kFlxCurrent        =  3   // current flux from systematic param variation of kFlxFinalTuning flux
} 
FluxType_t;

// Fit parameters
static const int knppmax = 100; // maximum number of physics params listed below (for booking large enough param arrays)
typedef enum EPhysParamId 
{   
   kPhypUndefined = -1, 
   kPhypMixing23,              // 23-mixing parameter (either sin^2(2*\theta_{23}) or sin^2(\theta_{23}) depending on config choices)
   kPhypMixing13,              // 13-mixing parameter (either sin^2(2*\theta_{13}) or sin^2(\theta_{13}) depending on config choices)
   kPhypMixing12,              // 12-mixing parameter (either sin^2(2*\theta_{12}) or sin^2(\theta_{12}) depending on config choices)
   kPhypDeltaCP,               // CP-invariance violating phase 
   kPhypDm23sq,                // modulus of squared-mass splitting \Deltam_{23}^{2} (always positive independent of mass hierarchy) 
   kPhypDm12sq,                // squared-mass splitting \Deltam_{12}^{2} 
   kPhypAlphaDecay,            //
   kPhypSinsq2ThetaDecay,      // 
   kPhypMusqDecoh,             //
   kPhypSinsq2ThetaDecoh,      //

   /*
    * List of 2011 analysis nuisance parameters 
    * -----------------------------------------
   */

   /*
    * SuperK detector errors (energy scale and efficiencies)
   */   
   kNuisSuperkEnergyScale,     // SK energy scale
   kNuisSuperkEff1RmuCCQEoth,  // SK 1Rmu CCQE efficiency (other)
   kNuisSuperkEff1RmuCCQErc0u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 1)
   kNuisSuperkEff1RmuCCQErc1u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 2)
   kNuisSuperkEff1RmuCCQErc2u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 3)
   kNuisSuperkEff1RmuCCnonQE,  // SK 1Rmu CCnonQE efficiency
   kNuisSuperkEff1RmuNC,       // SK 1Rmu NC efficiency
   kNuisSuperkEff1RmuNueCC,    // SK 1Rmu nue CC efficiency

   /*
    * 27 (flux and cross-section) corellated parameters used in the BANFF fit
   */   
   kNuisBANFF00,               // # 00, SK     numu    flux,  E =  0.0 -  0.4 GeV
   kNuisBANFF01,               // # 01, SK     numu    flux,  E =  0.4 -  0.5 GeV
   kNuisBANFF02,               // # 02, SK     numu    flux,  E =  0.5 -  0.6 GeV
   kNuisBANFF03,               // # 03, SK     numu    flux,  E =  0.6 -  0.7 GeV
   kNuisBANFF04,               // # 04, SK     numu    flux,  E =  0.7 -  1.0 GeV
   kNuisBANFF05,               // # 05, SK     numu    flux,  E =  1.0 -  1.5 GeV
   kNuisBANFF06,               // # 06, SK     numu    flux,  E =  1.5 -  2.5 GeV
   kNuisBANFF07,               // # 07, SK     numu    flux,  E =  2.5 -  3.5 GeV
   kNuisBANFF08,               // # 08, SK     numu    flux,  E =  3.5 -  5.0 GeV
   kNuisBANFF09,               // # 09, SK     numu    flux,  E =  5.0 -  7.0 GeV
   kNuisBANFF10,               // # 10, SK     numu    flux,  E =  7.0 - 30.0 GeV
   kNuisBANFF11,               // # 11, SK     numubar flux,  E =  0.0 -  1.5 GeV
   kNuisBANFF12,               // # 12, SK     numubar flux,  E =  1.5 - 30.0 GeV
   kNuisBANFF13,               // # 13, SK     nue     flux,  E =  0.0 -  0.5 GeV
   kNuisBANFF14,               // # 14, SK     nue     flux,  E =  0.5 -  0.7 GeV
   kNuisBANFF15,               // # 15, SK     nue     flux,  E =  0.7 -  0.8 GeV
   kNuisBANFF16,               // # 16, SK     nue     flux,  E =  0.8 -  1.5 GeV
   kNuisBANFF17,               // # 17, SK     nue     flux,  E =  1.5 -  2.5 GeV
   kNuisBANFF18,               // # 18, SK     nue     flux,  E =  2.5 -  4.0 GeV
   kNuisBANFF19,               // # 19, SK     nue     flux,  E =  4.0 - 30.0 GeV
   kNuisBANFF20,               // # 20, SK     nuebar  flux,  E =  0.0 -  2.5 GeV
   kNuisBANFF21,               // # 21, SK     nuebar  flux,  E =  2.5 - 30.0 GeV
   kNuisBANFF22,               // # 22, M_{A}  QE
   kNuisBANFF23,               // # 23, M_{A}  RES
   kNuisBANFF24,               // # 24, CCQE   normalization, 0 - 1.5 GeV
   kNuisBANFF25,               // # 25, CC1pi  normalization, 0 - 2.5 GeV
   kNuisBANFF26,               // # 26, NC1pi0 normalization

   /*
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.
   */   
   //systematics which use splines
   kNuisXSecFermiMom,          // Fermi momentum
   kNuisXSecBindE,             // binding energy
   kNuisXSecSpecFunc,          // spectral function on/off
   kNuisXSecWshape,            // W-shape
   kNuisXSecShapeCCoth,        // CC other shape
   //systematics that are normalisations
   kNuisXSecNormCCQEE1,        // CCQE norm 1.5 - 3.5 GeV
   kNuisXSecNormCCQEE2,        // CCQE norm > 3.5 GeV
   kNuisXSecNormCC1piE1,       // CC1pi norm > 2.5 GeV   
   kNuisXSecNormCCcoh,         // CC coherent norm
   kNuisXSecNormNCcoh,         // NC coherent norm
   kNuisXSecNormNCoth,         // NC other norm
   kNuisXSecNormNuEToNuMu,     // nue / numu cross section

   /*
    * FSI parameters.
    * Assuming here a 10x10 covariance matrix will be made available.
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)
    */
    kNuisFSI00,                // # 00, ??
    kNuisFSI01,                // # 01, ??
    kNuisFSI02,                // # 02, ??
    kNuisFSI03,                // # 03, ??
    kNuisFSI04,                // # 04, ??
    kNuisFSI05,                // # 05, ??
    kNuisFSI06,                // # 06, ??
    kNuisFSI07,                // # 07, ??
    kNuisFSI08,                // # 08, ??
    kNuisFSI09                 // # 09, ??
} 
PhysParamId_t;


//
// ------------------------------------------------------------------------
// ** NOTE: IMPORTANT CHOICES MADE HERE **
//
// Choice of samples included in the analysis 
// Choice of modes (per sample) we are going to keep track of.
//
// Note that some of the arrays defined further below (such as formatting 
// choices or binning choices for systematic parameter weight splines) 
// depend on kSKReactionModes. Make sure changes are made consistently.
// ------------------------------------------------------------------------
//

//
// Define which SK samples are considered
//

static const SuperkSample_t kSKSamples[kns] = {
    kSKsFCFV1RingMuLike,
    kSKsUndefined            //<-- add kSKsFCFV1RingELike for the combined fit
};

//
// Define which SK modes are considered for each SK sample included in the analysis
//

static const SuperkMode_t kSKReactionModes[kns][knm] = {

  // > 1 muon-like ring (numu disappearance)
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCCoth_mu,                      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,          
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNCoth_mu,
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECCoth_mu,                                
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth_mu,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
    kSKmSigNuECCQE,
    kSKmSigNuECC1pi,
    kSKmSigNuECCcoh,
    kSKmSigNuECCoth_mu,
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
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth_e,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth_e,
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCoth_e,
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
  }
};

//
// ------------------------------------------------------------------------
//
// Info on systematic parameters 
//
// ------------------------------------------------------------------------
//

//
// Number of nuisance params parameters for 1Rmu, 1Re, 1Rmu+1Re analysis
// 1Rmu : 57 nuisance params
//   - 1 SK E scale 
//   - 7 SK efficiencies (3 CCQE-RC, CCQE-Oth, CCnonQE, NC, nue CC) 
//   - 27 BANFF
//   - 12 XSec
//   - 10 FSI
// 1Re :
//   -
// Joint :
//   -
//
static const int kNumOfNuisanceParams1Rmu  = 57;  
static const int kNumOfNuisanceParams1Re   = 0;  // not implemented for 2012a yet
static const int kNumOfNuisanceParamsJoint = 0;  // not implemented for 2012a yet

//
// Number of BANFF fit parameters
//
static const int kNumOfParamsInBANFFfit = 27;  

//
// ** Flux systematics
//

// Mark which BANFF parameters are SuperK flux params
static const bool kParamInBANFFfitIsSuperkFluxParam[kNumOfParamsInBANFFfit] =
{
   true,                 // kNuisBANFF00
   true,                 // kNuisBANFF01
   true,                 // kNuisBANFF02
   true,                 // kNuisBANFF03
   true,                 // kNuisBANFF04
   true,                 // kNuisBANFF05
   true,                 // kNuisBANFF06
   true,                 // kNuisBANFF07
   true,                 // kNuisBANFF08
   true,                 // kNuisBANFF09
   true,                 // kNuisBANFF10
   true,                 // kNuisBANFF11
   true,                 // kNuisBANFF12
   true,                 // kNuisBANFF13
   true,                 // kNuisBANFF14
   true,                 // kNuisBANFF15
   true,                 // kNuisBANFF16
   true,                 // kNuisBANFF17
   true,                 // kNuisBANFF18
   true,                 // kNuisBANFF19
   true,                 // kNuisBANFF20
   true,                 // kNuisBANFF21
   false,                // kNuisBANFF22
   false,                // kNuisBANFF23
   false,                // kNuisBANFF24
   false,                // kNuisBANFF25
   false                 // kNuisBANFF26
};

// flux-related BANFF parameters --> flux neutrino PDG code
static const int kNuPdgCodeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit] =
{
   kPdgNuMu,             // kNuisBANFF00
   kPdgNuMu,             // kNuisBANFF01
   kPdgNuMu,             // kNuisBANFF02
   kPdgNuMu,             // kNuisBANFF03
   kPdgNuMu,             // kNuisBANFF04
   kPdgNuMu,             // kNuisBANFF05
   kPdgNuMu,             // kNuisBANFF06
   kPdgNuMu,             // kNuisBANFF07
   kPdgNuMu,             // kNuisBANFF08
   kPdgNuMu,             // kNuisBANFF09
   kPdgNuMu,             // kNuisBANFF10
   kPdgNuMuBar,          // kNuisBANFF11
   kPdgNuMuBar,          // kNuisBANFF12
   kPdgNuE,              // kNuisBANFF13
   kPdgNuE,              // kNuisBANFF14
   kPdgNuE,              // kNuisBANFF15
   kPdgNuE,              // kNuisBANFF16
   kPdgNuE,              // kNuisBANFF17
   kPdgNuE,              // kNuisBANFF18
   kPdgNuE,              // kNuisBANFF19
   kPdgNuEBar,           // kNuisBANFF20
   kPdgNuEBar,           // kNuisBANFF21
   0,                    // kNuisBANFF22, not a flux param
   0,                    // kNuisBANFF23, not a flux param
   0,                    // kNuisBANFF24, not a flux param
   0,                    // kNuisBANFF25, not a flux param
   0                     // kNuisBANFF26, not a flux param
};

// Etrue energy range for flux-related BANFF parameters
static const double kEtrueRangeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit][2] =
{
   { 0.0,  0.4 },        // kNuisBANFF00
   { 0.4,  0.5 },        // kNuisBANFF01
   { 0.5,  0.6 },        // kNuisBANFF02
   { 0.6,  0.7 },        // kNuisBANFF03
   { 0.7,  1.0 },        // kNuisBANFF04
   { 1.0,  1.5 },        // kNuisBANFF05
   { 1.5,  2.5 },        // kNuisBANFF06
   { 2.5,  3.5 },        // kNuisBANFF07
   { 3.5,  5.0 },        // kNuisBANFF08
   { 5.0,  7.0 },        // kNuisBANFF09
   { 7.0, 30.0 },        // kNuisBANFF10
   { 0.0,  1.5 },        // kNuisBANFF11
   { 1.5, 30.0 },        // kNuisBANFF12
   { 0.0,  0.5 },        // kNuisBANFF13
   { 0.5,  0.7 },        // kNuisBANFF14
   { 0.7,  0.8 },        // kNuisBANFF15
   { 0.8,  1.5 },        // kNuisBANFF16
   { 1.5,  2.5 },        // kNuisBANFF17
   { 2.5,  4.0 },        // kNuisBANFF18
   { 4.0, 30.0 },        // kNuisBANFF19
   { 0.0,  2.5 },        // kNuisBANFF20
   { 2.5, 30.0 },        // kNuisBANFF21
   {  -1,   -1 },        // kNuisBANFF22, not a flux param
   {  -1,   -1 },        // kNuisBANFF23, not a flux param
   {  -1,   -1 },        // kNuisBANFF24, not a flux param
   {  -1,   -1 },        // kNuisBANFF25, not a flux param
   {  -1,   -1 }         // kNuisBANFF26, not a flux param
};

//
// ** Cross-section systematics
//

// Nominal values of cross-section related params 
static const double kNomMAQE  = 1.21;
static const double kNomMARES = 1.21;

// Number of cross-section related systematics 
static const int kNumOfXSecSyst = 17;

// The following is the list of the above kNumOfXSecSyst parameters
static const PhysParamId_t kListOfXSecSyst[kNumOfXSecSyst] = 
{
   kNuisBANFF22,               // # 22, M_{A}  QE
   kNuisBANFF23,               // # 23, M_{A}  RES
   kNuisBANFF24,               // # 24, CCQE   normalization, 0 - 1.5 GeV
   kNuisBANFF25,               // # 25, CC1pi  normalization, 0 - 2.5 GeV
   kNuisBANFF26,               // # 26, NC1pi0 normalization
   kNuisXSecFermiMom,          // Fermi momentum
   kNuisXSecBindE,             // binding energy
   kNuisXSecSpecFunc,          // spectral function on/off
   kNuisXSecWshape,            // W-shape
   kNuisXSecShapeCCoth,        // CC other shape
   kNuisXSecNormCCQEE1,        // CCQE norm 1.5 - 3.5 GeV
   kNuisXSecNormCCQEE2,        // CCQE norm > 3.5 GeV
   kNuisXSecNormCC1piE1,       // CC1pi norm > 2.5 GeV   
   kNuisXSecNormCCcoh,         // CC coherent norm
   kNuisXSecNormNCcoh,         // NC coherent norm
   kNuisXSecNormNCoth,         // NC other norm
   kNuisXSecNormNuEToNuMu      // nue / numu cross section
};

// The following list specifies which cross-section related systematic has its
// effect on MC templates pre-computed and stored in cubic splines.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const bool kXSecSystHasWghtSplines[kNumOfXSecSyst] = 
{
   true,                      // kNuisBANFF22
   true,                      // kNuisBANFF23
   false,                     // kNuisBANFF24
   false,                     // kNuisBANFF25
   false,                     // kNuisBANFF26
   true,                      // Fermi momentum
   true,                      // binding energy
   true,                      // spectral function on/off
   true,                      // W-shape
   true,                      // CC other shape
   false,                     // CCQE norm 1.5 - 3.5 GeV
   false,                     // CCQE norm > 3.5 GeV
   false,                     // CC1pi norm > 2.5 GeV   
   false,                     // CC coherent norm
   false,                     // NC coherent norm
   false,                     // NC other norm
   false                      // nue / numu cross section
};

// The following list specifies which cross-section related systematic 
// parameters were included in the BANFF fit.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const bool kXSecSystIncludedInBANFF[kNumOfXSecSyst] = 
{
   true,                      // kNuisBANFF22
   true,                      // kNuisBANFF23
   true,                      // kNuisBANFF24
   true,                      // kNuisBANFF25
   true,                      // kNuisBANFF26
   false,                     // Fermi momentum
   false,                     // binding energy
   false,                     // spectral function on/off
   false,                     // W-shape
   false,                     // CC other shape
   false,                     // CCQE norm 1.5 - 3.5 GeV
   false,                     // CCQE norm > 3.5 GeV
   false,                     // CC1pi norm > 2.5 GeV   
   false,                     // CC coherent norm
   false,                     // NC coherent norm
   false,                     // NC other norm
   false                      // nue / numu cross section
};

// The following list specifies which cross-section related systematic is
// calculated for a var.-size Ereco,Etrue bin and then needs to be distributed
// unequally between the MC template bins (eg in proportion to number of events in each bin)
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const bool kXSecSystWghtUnequallySharedBetweenBins[kNumOfXSecSyst] = 
{
   true,                      // kNuisBANFF22
   true,                      // kNuisBANFF23
   false,                     // kNuisBANFF24
   false,                     // kNuisBANFF25
   false,                     // kNuisBANFF26
   true,                      // Fermi momentum
   true,                      // binding energy
   true,                      // spectral function on/off
   true,                      // W-shape
   true,                      // CC other shape
   false,                     // CCQE norm 1.5 - 3.5 GeV
   false,                     // CCQE norm > 3.5 GeV
   false,                     // CC1pi norm > 2.5 GeV   
   false,                     // CC coherent norm
   false,                     // NC coherent norm
   false,                     // NC other norm
   false                      // nue / numu cross section
};

// Ereco and Etrue cuts for the application of cross-section systematics.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const double kXSecSystErecoCut[kNumOfXSecSyst][2] =
{
   { 0.0, 100.0 },            // kNuisBANFF22
   { 0.0, 100.0 },            // kNuisBANFF23
   { 0.0, 100.0 },            // kNuisBANFF24
   { 0.0, 100.0 },            // kNuisBANFF25
   { 0.0, 100.0 },            // kNuisBANFF26
   { 0.0, 100.0 },            // Fermi momentum
   { 0.0, 100.0 },            // binding energy
   { 0.0, 100.0 },            // spectral function on/off
   { 0.0, 100.0 },            // W-shape
   { 0.0, 100.0 },            // CC other shape
   { 0.0, 100.0 },            // CCQE norm 1.5 - 3.5 GeV
   { 0.0, 100.0 },            // CCQE norm > 3.5 GeV
   { 0.0, 100.0 },            // CC1pi norm > 2.5 GeV   
   { 0.0, 100.0 },            // CC coherent norm
   { 0.0, 100.0 },            // NC coherent norm
   { 0.0, 100.0 },            // NC other norm
   { 0.0, 100.0 }             // nue / numu cross section
};
static const double kXSecSystEtrueCut[kNumOfXSecSyst][2] =
{
   { 0.0, 100.0 },            // kNuisBANFF22
   { 0.0, 100.0 },            // kNuisBANFF23
   { 0.0,   1.5 },            // kNuisBANFF24
   { 0.0,   2.5 },            // kNuisBANFF25
   { 0.0, 100.0 },            // kNuisBANFF26
   { 0.0, 100.0 },            // Fermi momentum
   { 0.0, 100.0 },            // binding energy
   { 0.0, 100.0 },            // spectral function on/off
   { 0.0, 100.0 },            // W-shape
   { 0.0, 100.0 },            // CC other shape
   { 1.5,   3.5 },            // CCQE norm 1.5 - 3.5 GeV
   { 3.5, 100.0 },            // CCQE norm > 3.5 GeV
   { 2.5, 100.0 },            // CC1pi norm > 2.5 GeV   
   { 0.0, 100.0 },            // CC coherent norm
   { 0.0, 100.0 },            // NC coherent norm
   { 0.0, 100.0 },            // NC other norm
   { 0.0, 100.0 }             // nue / numu cross section
};

// The following list describles in which MC templates we are supposed to apply the
// weights as a result of tweaking any of the above systematic parameters.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

  // > 1 muon-like ring 
  {
    // banff22  banff23  banff24  banff25  banff26  pF    Eb   spectral   W-shape  ccothshape
    // ccqeE2   ccqeE3   cc1piE2  cccoh    nccoh  ncoth    nue/numu
    {  1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNCoth_mu
    {  1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarNCoth_mu
    {  1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, // kSKmNuENCoth_mu
    {  1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, // kSKmNuEBarNCoth_mu
    {  1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1 }, // kSKmSigNueCCQE        
    {  0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1 }, // kSKmSigNueCC1pi
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1 }, // kSKmSigNueCCcoh
    {  0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmSigNueCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } // kSKmUndefined
  },

  // > 1 electron-like ring 
  {
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmSigNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
};

//
// The weight splines are pre-computed for the following Ereco bins.
// For the time being, the same Ereco binning is used for all samples and all modes.
//

static const int kXsNErecoBins = 28;

static const double kXsErecoBinEdges[kXsNErecoBins+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  
   3.000, 10.000, 30.000
};

//
// The weight splines are pre-computed in the following Etrue bins (for each Ereco bin)
// The Etrue bins can, in principle, be different for each sample/mode combination.
// For the time-being, the Etrue bins are the BANFF ones (depend on the neutrino code).
// Note that the order of modes below should match the one defined in kSKReactionModes array.
//

static const int kXsNEtrueBins = 11;

static const double kXsEtrueBinEdges[kns][knm][kXsNEtrueBins+1] = {

  // > 1 muon-like ring 
  {
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCQE         
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCC1pi
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCcoh
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCoth_mu                      
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNC1piC
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNCoth_mu
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCQE      
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCC1pi
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCcoh
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCoth_mu                  
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNC1piC
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNCoth_mu
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCQE 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECC1pi 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCcoh
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCoth_mu                                
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENC1piC
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENCoth_mu
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCQE
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCC1pi
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCcoth_mu
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCoth_mu
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNC1piC
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNCoth_mu
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCQE 
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECC1pi
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCcoh
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCoth_mu
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }  // kSKmUndefined
  },

  // > 1 electron-like ring 
  {
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCQE 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCC1pi 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCcoh 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCoth_e 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNC1pi0 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNCcoh 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNCoth_e                            
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCQE 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCC1pi 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCcoh 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCoth_e  
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNC1pi0 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNCcoh 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNCoth_e
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCQE 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECC1pi 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCcoh 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCoth_e
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENC1pi0  
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENCcoh 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENCoth_e
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCQE
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCC1pi
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCcoh
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCoth_e
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNC1pi0
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNCcoh
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNCoth_e
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCQE  
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECC1pi  
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCoth_e
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }  // kSKmUndefined
  }
};

// Number of cross-section weight splines
static const int kNumOfXSecSystWghtSpl = kNSuperkSamples * kNSuperkModes * kNumOfXSecSyst * kXsNErecoBins * kXsNEtrueBins;
static const int knwspl = kNumOfXSecSystWghtSpl; // as above - shorter name

//
// ** FSI systematics
//

// Number of FSI systematic params
// NOTE: COMPLETELY HYPOTHETICAL FOR THE TIME-BEING
static const int kNumOfFSIParams = 10;

// Ereco energy range for FSI systematic params
// NOTE: COMPLETELY HYPOTHETICAL FOR THE TIME-BEING
static const double kFSIParamsErecoBinEdges[kNumOfFSIParams+1] =
{
//    [00]  [01]  [02]  [03]  [04]  [05]  [06]  [07]  [08]  [09]
   0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0 
};


//
// ------------------------------------------------------------------------
// Formatting etc
// ------------------------------------------------------------------------
//

//
// Default graphic attributes for output spectra. Each line refers to a mode of a given sample.
// List of samples and modes should match the one given in kSKReactionModes[kns][knm] above.
// The 5 listed attributes for each sample are: Fill style, Fill color, Line style, Line width, Line color.
//
static const int kSKHistoAttributes[kns][knm][kNGraphicAttributes] = {

  // > 1 muon-like ring 
  {
    {   -1,  kWhite,    1,  1,  kBlack}, // kSKmNuMuCCQE
    { 3544,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCC1pi
    { 3244,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCcoh
    { 3006,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCoth_mu
    { 3554,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC1piC
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNCoth_mu
    { 3001,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCQE
    { 3345,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC1pi
    { 3244,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCcoh
    { 3354,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCoth_mu
    { 3544,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC1piC
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNCoth_mu
    { 3001,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCQE
    { 3345,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC1pi
    { 3244,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCcoh
    { 3354,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCoth_mu
    { 3544,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC1pic
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENCoth_mu
    { 3001,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCQE
    { 3345,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC1pi
    { 3244,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCcoh
    { 3354,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCoth_mu
    { 3544,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC1pic
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNCoth_mu
    { 3001,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECCQE
    { 3345,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECC1pi
    { 3244,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECCcoh
    { 3354,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECCoth_mu
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }  // kSKmUndefined
  },

  // > 1 electron-like ring 
  {
    {   -1,    kWhite,  1,  1, kBlack }, // kSKmNuMuCCQE
    { 3544,    kBlack,  1,  1, kBlack }, // kSKmNuMuCC1pi
    { 3244,    kGreen,  1,  1, kBlack }, // kSKmNuMuCCcoh
    { 3244,     kBlue,  1,  1, kBlack }, // kSKmNuMuCCoth_e
    { 3006,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1pi0
    { 3006,    kRed-1,  1,  1, kBlack }, // kSKmNuMuNCcoh
    { 3006,    kRed-2,  1,  1, kBlack }, // kSKmNuMuNCoth_e
    { 3001,  kAzure-3,  1,  1, kBlack }, // kSKmNuMuBarCCQE
    { 3345,  kAzure-4,  1,  1, kBlack }, // kSKmNuMuBarCC1pi
    { 3354,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCCcoh
    { 3354,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarCCoth_e
    { 3006,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-9,  1,  1, kBlack }, // kSKmNuMuBarNCcoh
    { 3006,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarNCoth_e
    { 3001,  kGreen+1,  1,  1, kBlack }, // kSKmNuECCQE
    { 3345,  kGreen+2,  1,  1, kBlack }, // kSKmNuECC1pi
    { 3354,  kGreen+2,  1,  1, kBlack }, // kSKmNuECCcoh
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuECCoth_e
    { 3006,  kGreen+5,  1,  1, kBlack }, // kSKmNuENC1pi0
    { 3006,  kGreen+7,  1,  1, kBlack }, // kSKmNuENCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuENCoth_e
    { 3001,  kBlue+1,   1,  1, kBlack }, // kSKmNuEBarCCQE
    { 3345,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCC1pi
    { 3354,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCCcoh
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuEBarCCoth_e
    { 3006,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarNC1pi0
    { 3006,  kBlue+7,   1,  1, kBlack }, // kSKmNuEBarNCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuEBarNCoth_e
    { 3001,  kGreen+8,  1,  1, kBlack }, // kSKmSigNuECCQE
    { 3345, kGreen+10,  1,  1, kBlack }, // kSKmSigNuECC1pi
    { 3354,      kRed,  1,  1, kBlack }, // kSKmSigNuECCoth_e
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }  // kSKmUndefined

  }
};

//
// Define the order in which TH1D spectra are added in a THStack.
// The order is chosen to increase visibility of smaller components.
//

static const SuperkMode_t kSKReactionModesOrderInStack[kns][knm] = {

  // > 1 muon-like ring (numu disappearance)
  {
    kSKmSigNuECCcoh, 
    kSKmSigNuECC1pi,
    kSKmSigNuECCoth_mu,
    kSKmSigNuECCQE, 
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCoth_mu,                                
    kSKmNuEBarCCQE, 
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmNuECCcoh,
    kSKmNuECC1pi,
    kSKmNuECCoth_mu,                                
    kSKmNuECCQE, 
    kSKmNuMuBarNC1piC, 
    kSKmNuMuBarNCoth_mu,
    kSKmNuMuBarCCcoh, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarCCQE,      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,
    kSKmNuMuCCoth_mu,                      
    kSKmNuMuCCcoh, 
    kSKmNuMuCC1pi, 
    kSKmNuMuCCQE,         
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
    kSKmNuEBarCCQE, 
    kSKmNuEBarCC1pi, 
    kSKmNuEBarCCcoh, 
    kSKmNuEBarCCoth_e,
    kSKmNuEBarNC1pi0,  
    kSKmNuEBarNCcoh, 
    kSKmNuEBarNCoth_e,
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    kSKmNuECCoth_e,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCoth_e,
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
  }
};

//
// ------------------------------------------------------------------------
// Binning 
// ------------------------------------------------------------------------
//

//
// Etrue binning for SK MC templates
//

static const double kEtrueMin    = 0;
static const double kEtrueMax    = 30;
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
// 50 bins from 0 to 10 GeV (same as the 2010a NIWG input for SKFSI)
static const int kNErecoBins1Re = 6;
static const double kErecoBinEdges1Re[kNErecoBins1Re+1] = 
{
  0.000,  0.200,  0.400, 0.600, 0.800,  1.000,  1.250
};
/*
// 1 bin from 0 to 10 GeV (counting expt)
//WARNING this gives wrong errors for anything that is Ereco dependant (SKFSI)
static const int kNErecoBins1Re = 1;
static const double kErecoBinEdges1Re[kNErecoBins1Re+1] = 
{ 
   0., 1.250 
};
*/
/*
// Same binning as the standard 1Rmu one: 200 50-MeV bins from 0 to 10 GeV
static const int kNErecoBins1Re = 25;
static const double kErecoBinEdges1Re[kNErecoBins1Re+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250
};
*/

} // conv
} // t2kana
  
#endif 


