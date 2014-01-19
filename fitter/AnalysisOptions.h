//____________________________________________________________________________
/*!

\class    t2kana::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since 2011_v0.9.0:

17-12-2013  Nick  
Split t2kana::conv namespace into several smaller files

*/
//____________________________________________________________________________

#ifndef ana_options_h
#define ana_options_h

#include "Build.h"

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#ifndef __FULLY_CORRELATED__
#define __FULLY_CORRELATED__
#endif

#ifndef SKFSITOGETHER
#define SKFSITOGETHER
// If turned on, errors for SK detector efficiencies and FSI are combined
// fits are faster in this way. 
#endif

#ifndef HKHACK
//#define HKHACK //use signal/background only systematics
#endif

#ifndef HK
#define HK //allow mix23 & dm32sq to be nuisance for current osc params
#endif

#ifndef BANFF_BINS
#ifdef __NUMUONLY__
#define BANFF_BINS 2 //numu disapp without MEC
//#define BANFF_BINS 4 //numu disapp analysis, 25 bins (includes MEC systematics)
#else
#ifdef __NUEONLY__
#define BANFF_BINS 3 //nue app
#else
#define BANFF_BINS 1 //joint
#endif
#endif
#endif


#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
//#ifndef __NUEONLY__
//#define __USE_OSC_NUMU_MC_TEMPLATES__
// If turned on, oscillated numu (from nue) templates will be added (to be weighted with P(nue->numu)).
// Turn off for compatibility with 2012-era analysis.
//#endif
#endif


#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
//#ifndef __NUEONLY__
//#define __USE_OSC_NUMUBAR_MC_TEMPLATES__
// If turned on, oscillated numu (from nue) templates will be added (to be weighted with P(nue->numu)).
// Turn off for compatibility with 2012-era analysis.
//#endif
#endif


#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
//#ifndef __NUEONLY__
//#define __USE_OSC_NUEBAR_MC_TEMPLATES__
// If turned on, oscillated numu (from nue) templates will be added (to be weighted with P(nue->numu)).
// Turn off for compatibility with 2012-era analysis.
//#endif
#endif

//#ifndef __USE_EXTRA_MC_TEMPLATES__
//#ifndef __NUEONLY__
//#define __USE_EXTRA_MC_TEMPLATES__
// If turned on, oscillated numu (from nue) templates will be added (to be weighted with P(nue->numu)).
// Turn off for compatibility with 2012-era analysis.
//#endif
//#endif


#ifndef USE_MEC
//#define USE_MEC //if this is turned off, kSKmMEC is not included in kSKReactionModes
#endif

#ifndef FSI_BINS
#ifdef __NUMUONLY__
#define FSI_BINS 2 //joint analysis, 18 bins
#else
#ifdef __NUEONLY__
#define FSI_BINS 3 //numu disapp analysis, 6 bins
#else
#define FSI_BINS 1 //numu disapp analysis, 6 bins
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
#define ERECO_XSEC_BINS_1RE 2 // DEFAULT 0-0.3 GeV, 50 MeV 0.3-1.25 GeV (20 bins)
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
#endif

#ifndef ERECO_BINS_1RE
//#define ERECO_BINS_1RE 1 // DEFAULT 200 MeV bins 0-1 GeV, 1 bin from 1-1.25 GeV (6 bins)
//#define ERECO_BINS_1RE 2 // 1 bin (0-1.25 GeV)
#define ERECO_BINS_1RE 3 // 50 MeV bins 0-1.25 GeV (25 bins)
#endif

#endif

