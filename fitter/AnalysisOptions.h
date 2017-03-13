
//____________________________________________________________________________
/*!

\class    valor::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since 2011_v0.9.0:

17-12-2013  Nick  
Split valor::conv namespace into several smaller files

*/
//____________________________________________________________________________

#ifndef ana_options_h
#define ana_options_h

#include "Build.h"

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#define NUMOFTOYSTOMARGINALIZE 5 // Ntoys

#define SYST_NSIGMA   3   // number of sigma away from nominal nuisance parameters are allowed to be
#define SYST_STEP     0.1 // step size (in units of sigma)

// NEW 2D TEMPLATES
#ifndef __USE_2D_TEMPLATES__ 
#define __USE_2D_TEMPLATES__ // Use 2D templates
#endif

//#define __ERECO_VS_Q2_TEMPLATES__ // Ereco Vs Q2 parametrization
#ifndef __ERECO_VS_Q2_TEMPLATES__ 
#define __ERECO_VS_THETA_TEMPLATES__ // Ereco Vs THETA parametrization
#endif

//

#ifndef __FULLY_CORRELATED__ 
//#define __FULLY_CORRELATED__ // fully correlated rhc-fhc samples
#endif

// XSEC (either niwg2015 or niwg2016)
#define __USE_NIWG2017__
#ifndef __USE_NIWG2017__
#define __USE_NIWG2015__
#endif

// NEW SKDETFSI (either 2013 or 2015)
//#define __USE_SKDETFSI2015__
#ifndef __USE_SKDETFSI2015__
#define __USE_SKDETFSI2013__
#endif


#ifndef __BANFF_UNCORRELATED__
//#define __BANFF_UNCORRELATED__
#endif

#ifndef SWAP_RHC_NUMU_NUMUBAR_ERRORS
//#define SWAP_RHC_NUMU_NUMUBAR_ERRORS
#endif

#ifndef NEW_RHC_BANFF
//#define NEW_RHC_BANFF //1.06* fhc error
#endif

/// NEW MULTIRING 2
#ifndef __FULLY_CORRELATED_MULTIRING__ // multiring-fhc samples
//#define __FULLY_CORRELATED_MULTIRING__
#endif
///

#ifndef HKHACK
//#define HKHACK //use signal/background only systematics
#endif

#ifndef HK
#define HK //allow mix23 & dm32sq to be nuisance for current osc params
#endif

#ifndef BANFF_BINS
#ifdef __NUMUONLY__
#define BANFF_BINS 1 //numu disapp without 2p2h // keep fixed SCC parameters
#else
#ifdef __NUEONLY__
#define BANFF_BINS 3 //nue app
#else
#define BANFF_BINS 1 //joint
#endif
#endif
#endif

#ifndef FSI_BINS
#ifdef __NUMUONLY__
#define FSI_BINS 1 //numu disapp analysis, 6 bins //Keep same as joint
#else
#ifdef __NUEONLY__
#define FSI_BINS 3 //nue app analysis, 12 bins
#else
#define FSI_BINS 1 //joint analysis, 18 bins
#endif
#endif
#endif

#ifndef MR_FSI_BINS
#define MR_FSI_BINS 3
#endif

#ifndef ONEOVERMAQE
//#define ONEOVERMAQE
// Turn on if using a BANFF matrix with 1overMAQE, rather than MAQE. 
// Changes output name in Utilities string methods, and well as relative error in NuXSec tweak methods.
#endif 

#ifndef ETRUE_BINS_1RMU
//#define ETRUE_BINS_1RMU 1 // 50 MeV binning 0-30 GeV (600 bins)
#define ETRUE_BINS_1RMU 2 // DEFAULT 50 MeV 0-0.3 GeV, 25 MeV 0.3-1 GeV, 50 MeV 1-3 GeV, 100 MeV 3-3.5 GeV, 3.5-4 GeV, 4-5 GeV, 5-7 GeV, 7-10 GeV, 10-30 GeV (84 bins)
//#define ETRUE_BINS_1RMU 3 // same binning as BANFF numu flux binning (11 bins)
//#define ETRUE_BINS_1RMU 4 // 1 bin
#endif

#ifndef ETRUE_BINS_1RE
//#define ETRUE_BINS_1RE 1 // 50 MeV binning 0-3 GeV , 1GeV/bin 3-4 GeV,  10-30 GeV (69 bins)
//#define ETRUE_BINS_1RE 3 //  5 MeV binning 0-1 GeV, 50 MeV binning 1-30 GeV (780 bins)
//#define ETRUE_BINS_1RE 4 //  SPECIAL BINNING for testing -> 11 bins (numu) flux binning in BANFF
#define ETRUE_BINS_1RE 2 // DEFAULT same as disap. 50 MeV 0-0.3 GeV, 25 MeV 0.3-1 GeV, 50 MeV 1-3 GeV, 100 MeV 3-3.5 GeV, 3.5-4 GeV, 4-5 GeV, 5-7 GeV, 7-10 GeV, 10-30 GeV (84 bins)
#endif

#ifndef ERECO_XSEC_BINS_1RMU
#define ERECO_XSEC_BINS_1RMU 1 // DEFAULT 0-0.3 GeV, 50 MeV 0.3-1 GeV, 100 MeV 1-3 GeV, 3-3.5 GeV, 3.5-4 GeV, 4-5 GeV, 5-7 GeV, 7-10 GeV, 10-30 GeV (41 bins)
//#define ERECO_XSEC_BINS_1RMU 2 // 50 MeV binning 0-30 GeV (600 bins)
#endif

#ifndef ERECO_XSEC_BINS_1RE
//#define ERECO_XSEC_BINS_1RE 1 // DEFAULT 0-0.3 GeV, 50 MeV 0.3-1.25 GeV (20 bins)
#define ERECO_XSEC_BINS_1RE 2 // DEFAULT 50 MeV 0.0-1.25 GeV (25 bins)
#endif

#ifndef SPLINE_TYPE
// NOTE1: binning is identical; just produced by different people
// NOTE2: don't need to change config xml file
//#define SPLINE_TYPE -1 // Special mode turned on in which splines are not read, and weights involving splines are set to 1
#define SPLINE_TYPE 1 // DEFAULT Use VaLOR splines
#endif

#ifndef ERECO_BINS_1RMU
//#define ERECO_BINS_1RMU 1 // 50 MeV binning 0-30 GeV
#define ERECO_BINS_1RMU 2 // DEFAULT 50 MeV bins in 0-3 GeV, 250 MeV bins in 3-4 GeV, 500 MeV bins in 4-6 GeV, 1000 MeV bins 6-10 GeV, 1 bin 10-30 GeV (73 bins)
//#define ERECO_BINS_1RMU 3 // 41 bins, same as ERECO_XSEC_BINS_1RMU 1
//#define ERECO_BINS_1RMU 4 //Same as Linda's binning scheme
//#define ERECO_BINS_1RMU 5 //39 bins, added for T2K Run 5-6 RHC analysis
#endif

#ifndef ERECO_BINS_1RE
//#define ERECO_BINS_1RE 1 // DEFAULT 200 MeV bins 0-1 GeV, 1 bin from 1-1.25 GeV (6 bins)
//#define ERECO_BINS_1RE 2 // 1 bin (0-1.25 GeV)
#define ERECO_BINS_1RE 3 // 50 MeV bins 0-1.25 GeV (25 bins)
//#define ERECO_BINS_1RE 4 //Same as Linda's binning scheme
#endif




// NEW 2D TEMPLATES

#ifndef Q2_XSEC_BINS_1RMU
//#define Q2_XSEC_BINS_1RMU 1 // Q^2 parametrization (19 bins) 
#define Q2_XSEC_BINS_1RMU 2 // 1 single bin (1D template) 
//#define Q2_XSEC_BINS_1RMU 3 // Qd^2 parametrization, finer close to 1 (19 bins) 
//#define Q2_XSEC_BINS_1RMU 4 // Inelasticity parametrization, finer close to 1 (19 bins) 
//#define Q2_XSEC_BINS_1RMU 5 // XSec term parametrization, (19 bins) 
//#define Q2_XSEC_BINS_1RMU 6 // XSec term parametrization, (92 bins) 
//#define Q2_XSEC_BINS_1RMU 7 // XSec term parametrization, Optimized (21 bins) 
//#define Q2_XSEC_BINS_1RMU 8 // Q^2 parametrization, Optimized (21 bins) 
#endif

#ifndef Q2_XSEC_BINS_1RE
//#define Q2_XSEC_BINS_1RE 1 // Q^2 parametrization (19 bins) 
//#define Q2_XSEC_BINS_1RE 2 // 1 single bin (1D template)
//#define Q2_XSEC_BINS_1RE 3 // Qd^2 parametrization, finer close to 1 (19 bins) 
//#define Q2_XSEC_BINS_1RE 4 // Inelasticity parametrization, finer close to 1 (19 bins) 
//#define Q2_XSEC_BINS_1RE 5 // XSec term parametrization, (19 bins) 
//#define Q2_XSEC_BINS_1RE 6 // XSec term parametrization, (89 bins)  
#define Q2_XSEC_BINS_1RE 7 // XSec term parametrization, Optimized (21 bins)  
//#define Q2_XSEC_BINS_1RE 8 // Q^2 parametrization, Optimized (21 bins) 
#endif

#ifndef Q2_BINS_1RMU
//#define Q2_BINS_1RMU 1 // Q^2 (19 bins)
#define Q2_BINS_1RMU 2 // 1 single bin (1D template)
//#define Q2_BINS_1RMU 3 // Qd^2 parametrization, finer close to 1 (19 bins) 
//#define Q2_BINS_1RMU 4 // Inelasticity parametrization, finer close to 1 (19 bins) 
//#define Q2_BINS_1RMU 5 // XSec term parametrization, (19 bins) 
//#define Q2_BINS_1RMU 6 // XSec term parametrization, (92 bins)   
//#define Q2_BINS_1RMU 7 // XSec term parametrization, Optimized (21 bins)   
//#define Q2_BINS_1RMU 8 // Q^2 parametrization, Optimized (21 bins)
#endif

#ifndef Q2_BINS_1RE
//#define Q2_BINS_1RE 1 // Q^2 (19 bins)
//#define Q2_BINS_1RE 2 // 1 single bin (1D template)
//#define Q2_BINS_1RE 3 // Qd^2 parametrization, finer close to 1 (19 bins) 
//#define Q2_BINS_1RE 4 // Inelasticity parametrization, finer close to 1 (19 bins) 
//#define Q2_BINS_1RE 5 // XSec term parametrization, (19 bins) 
//#define Q2_BINS_1RE 6 // XSec term parametrization, (89 bins)
#define Q2_BINS_1RE 7 // XSec term parametrization, Optimized (21 bins)
//#define Q2_BINS_1RE 8 // Q^2 parametrization, Optimized (21 bins)
#endif


#ifndef THETA_XSEC_BINS_1RMU
//#define THETA_XSEC_BINS_1RMU 1 // DEFAULT (15 bins) 
#define THETA_XSEC_BINS_1RMU 2 // 1 single bin 
#endif

#ifndef THETA_XSEC_BINS_1RE
//#define THETA_XSEC_BINS_1RE 1 // DEFAULT (15 bins) 
#define THETA_XSEC_BINS_1RE 2 // 1 single bin 
#endif

#ifndef THETA_XSEC_BINS_1RE_MULTIRING
//#define THETA_XSEC_BINS_1RE_MULTIRING 1 // DEFAULT (15 bins) 
#define THETA_XSEC_BINS_1RE_MULTIRING 2 // 1 single bin 
#endif

#ifndef THETA_BINS_1RMU
//#define THETA_BINS_1RMU 1 // DEFAULT (15 bins)
#define THETA_BINS_1RMU 2 // 1 single bin
#endif

#ifndef THETA_BINS_1RE
//#define THETA_BINS_1RE 1 // DEFAULT (15 bins)
#define THETA_BINS_1RE 2 // 1 single bin
#endif

#ifndef THETA_BINS_1RE_MULTIRING
//#define THETA_BINS_1RE_MULTIRING 1 // DEFAULT (15 bins)
#define THETA_BINS_1RE_MULTIRING 2 // 1 single bin
#endif

// Check to make sure that the samples are on in MultiRing FHC Running

#if defined(__MULTIRING_FHC_RUNNING__)
#if !(defined(__MULTIRING_FHC_NUMU__) || defined(__MULTIRING_FHC_NUE__))
#define __MULTIRING_FHC_NUMU__
#define __MULTIRING_FHC_NUE__
#endif
#endif

#endif //ana_options_h

