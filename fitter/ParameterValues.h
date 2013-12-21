//____________________________________________________________________________
/*!

\class    t2kana::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since 2011_v0.9.0:

17-12-2013  Nick  Split t2kana::conv namespace into several smaller files
*/
//____________________________________________________________________________

#ifndef parameter_values_h
#define parameter_values_h

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
// Updated to 2012 PDG + reactors
static const double kGlobFitDm23sq           = 2.4E-3;
static const double kGlobFitDmsqFL           = 2.36E-3;
static const double kGlobFitDm12sq           = 7.5E-5;
static const double kGlobFitErrDm23sq        = 0.04; // ~4%
static const double kGlobFitErrDm12sq        = 0.027; // ~2.7%
// deltacp global fit results
static const double kGlobFitDeltaCP          = 0.0;
static const double kGlobFitErrDeltaCP       = 0.10; // ~10%
// sinsq2theta global fit results
static const double kGlobFitSinsq2Theta23    = 1.00;
static const double kGlobFitSinsq2Theta13    = 0.098;
static const double kGlobFitSinsq2Theta12    = 0.8495; //fogli2011
//static const double kGlobFitSinsq2Theta12    = 0.857; //pdg2012
static const double kGlobFitErrSinsq2Theta23 = 0.05;  // ~5%
static const double kGlobFitErrSinsq2Theta12 = 0.031; // ~3.1% fogli2011
//static const double kGlobFitErrSinsq2Theta12 = 0.028; // ~2.8% pdg2012
static const double kGlobFitErrSinsq2Theta13 = 0.133; // ~13.3%
// sinsqtheta global fit results
static const double kGlobFitSinsqTheta23    = 0.50;
static const double kGlobFitSinsqTheta13    = 0.0251;
static const double kGlobFitSinsqTheta12    = 0.306; //fogli2011
//static const double kGlobFitSinsqTheta12    = 0.311; //pdg2012
static const double kGlobFitErrSinsqTheta23 = 0.06;  // ~26% (1 sigma range)
static const double kGlobFitErrSinsqTheta12 = 0.059;  // ~5.9% fogli2011
//static const double kGlobFitErrSinsqTheta12 = 0.053;  // ~5.3% pdg2012
static const double kGlobFitErrSinsqTheta13 = 0.139;  // ~13.9%

// projected future errors
//atmos params unknown
static const double kFutureAbsErrSinsqTheta23  = 0.05;//red line, page 36, PAC report final
static const double kFutureAbsErrSinsq2Theta23 = 0.19;//conversion of above (seems huge!)
static const double kFutureAbsErrDm23sq        = 0.045E-3;//red line, page 36, PAC report final
//mix13 using current daya bay systematic error
static const double kFutureAbsErrSinsqTheta13  = 0.0;
static const double kFutureAbsErrSinsq2Theta13 = 0.005;
//solar are known
static const double kFutureAbsErrSinsqTheta12  = 0.0;//if you make these non zero, also change FitInputs.cxx and add the relevant penalty term
static const double kFutureAbsErrSinsq2Theta12 = 0.0;
static const double kFutureAbsErrDm12sq        = 0.0;



// Dm^2 fitter init values and range and step
static const double kInitFitDm23sq           = kGlobFitDm23sq;
static const double kMinFitDm23sq            = 1.0E-3;
static const double kMaxFitDm23sq            = 6.0E-3;
static const double kStepDm23sq              = 1.0E-4;
static const double kInitFitDm12sq           = kGlobFitDm12sq;
static const double kMinFitDm12sq            = 6.0E-5;
static const double kMaxFitDm12sq            = 9.0E-5;
static const double kStepDm12sq              = 2.0E-6;
// sinsq2theta fitter init values and range
static const double kInitFitSinsq2Theta23    = 0.90;
static const double kMinFitSinsq2Theta23     = 0.00;
static const double kMaxFitSinsq2Theta23     = 1.00;
static const double kStepSinsq2Theta23       = 5.0E-3;
static const double kInitFitSinsq2Theta13    = kGlobFitSinsq2Theta13;
static const double kMinFitSinsq2Theta13     = 0.00;
static const double kMaxFitSinsq2Theta13     = 1.00;
static const double kStepSinsq2Theta13       = 1.0E-3;
static const double kInitFitSinsq2Theta12    = kGlobFitSinsq2Theta12;
static const double kMinFitSinsq2Theta12     = 0.00;
static const double kMaxFitSinsq2Theta12     = 1.00;
static const double kStepSinsq2Theta12       = 1.0E-2;
// sinsqtheta fitter init values and range
static const double kInitFitSinsqTheta23     = 0.60;
static const double kMinFitSinsqTheta23      = 0.00;
static const double kMaxFitSinsqTheta23      = 1.00;
static const double kStepSinsqTheta23        = 5.0E-3;
static const double kInitFitSinsqTheta13     = kGlobFitSinsqTheta13;
static const double kMinFitSinsqTheta13      = 0.00;
static const double kMaxFitSinsqTheta13      = 1.00;
static const double kStepSinsqTheta13        = 1.0E-3;
static const double kInitFitSinsqTheta12     = kGlobFitSinsqTheta12;
static const double kMinFitSinsqTheta12      = 0.00;
static const double kMaxFitSinsqTheta12      = 1.00;
static const double kStepSinsqTheta12        = 1.0E-2;

// fitter chisq surface ranges
static const double kMinChisqSurfaceDm23sq            = 1.5E-3;
static const double kMaxChisqSurfaceDm23sq            = 3.3E-3;
static const double kMinChisqSurfaceSinsq2Theta23     = 0.80;
static const double kMaxChisqSurfaceSinsq2Theta23     = 1.00;


//
// Other constants
//
static const double kMaxAllowedWght = 2E+3;                          // Complain if any computed weights exceeds this value
static const double kSetNegEventsTo = 1E-8;                          // Set any bins with nevent < 0 to this value
static const double kASmallNum      = 1E-10;                         // 
static const double kOneMeV         = 0.001;                         // 
static const double kBaseline       = 295;                           // Tokai->Kamioka (km)
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

} // conv
} // t2kana
  
#endif 

