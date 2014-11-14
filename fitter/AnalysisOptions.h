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

#define SYST_NSIGMA   3   // number of sigma away from nominal nuisance parameters are allowed to be
#define SYST_STEP     0.1 // step size (in units of sigma)

#ifndef __FULLY_CORRELATED__  // rhc-fhc samples
#define __FULLY_CORRELATED__
#endif

// NEW XSEC (either niwg2013 or niwg2015)
#define __USE_NIWG2015__
#ifndef __USE_NIWG2015__
#define __USE_NIWG2013__
#endif

// NEW SKDETFSI (either 2013 or 2015)
#define __USE_SKDETFSI2015__
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
#define BANFF_BINS 2 //numu disapp without MEC
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
#define FSI_BINS 2 //numu disapp analysis, 6 bins
#else
#ifdef __NUEONLY__
#define FSI_BINS 3 //nue app analysis, 12 bins
#else
#define FSI_BINS 1 //joint analysis, 18 bins
#endif
#endif
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
//#define ERECO_BINS_1RMU 4 //Same as Linda's binning scheme
//#define ERECO_BINS_1RMU 3 // 41 bins, same as ERECO_XSEC_BINS_1RMU 1
#endif

#ifndef ERECO_BINS_1RE
//#define ERECO_BINS_1RE 1 // DEFAULT 200 MeV bins 0-1 GeV, 1 bin from 1-1.25 GeV (6 bins)
//#define ERECO_BINS_1RE 2 // 1 bin (0-1.25 GeV)
#define ERECO_BINS_1RE 3 // 50 MeV bins 0-1.25 GeV (25 bins)
//#define ERECO_BINS_1RE 4 //Same as Linda's binning scheme
#endif

#endif

