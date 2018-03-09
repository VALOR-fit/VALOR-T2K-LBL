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

  31-12-2013  Nick
  Added extra arrays for reverse horn current modes

  21-7-2014 Davide
  Added samples for multiring forward horn current

*/
//____________________________________________________________________________

#ifndef sk_fsi_syst_conv_h
#define sk_fsi_syst_conv_h

#include "AnalysisOptions.h"
#include "Enumerations.h"
#include "SamplesAndModes.h"
#

namespace valor {
  namespace conv   {

    //
    // ** SKDET+FSI systematics
    //

    // Number of SKDET+FSI systematic params

    //
    // NB:
    //
    // We don't need kNumOfFHCpRHCSKDetFSIParams* because each skdetfsi parameter
    // is definitely independent from the others. Indeed they belong to different data samples.
    // It's different from the BANFF case, where the xsec parameters are common to
    // both the FHC and RHC modes
    //

    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
#if FSI_BINS == 1
    static const int kNumOfSKDetFSIParams     = 18;
    static const int kNumOfSKDetFSIParams1Rmu = 6;
    static const int kNumOfSKDetFSIParams1Re  = 12;
    static const int kNumOfRHCSKDetFSIParams     = 18;
    static const int kNumOfRHCSKDetFSIParams1Rmu = 6;
    static const int kNumOfRHCSKDetFSIParams1Re  = 12;
#elif FSI_BINS == 2
    static const int kNumOfSKDetFSIParams     = 6;
    static const int kNumOfSKDetFSIParams1Rmu = 6;
    static const int kNumOfSKDetFSIParams1Re  = 0;
    static const int kNumOfRHCSKDetFSIParams     = 6;
    static const int kNumOfRHCSKDetFSIParams1Rmu = 6;
    static const int kNumOfRHCSKDetFSIParams1Re  = 0;
#elif FSI_BINS == 3
    static const int kNumOfSKDetFSIParams     = 12;
    static const int kNumOfSKDetFSIParams1Rmu = 0;
    static const int kNumOfSKDetFSIParams1Re  = 12;
    static const int kNumOfRHCSKDetFSIParams     = 12;
    static const int kNumOfRHCSKDetFSIParams1Rmu = 0;
    static const int kNumOfRHCSKDetFSIParams1Re  = 12;
#endif
    
#if MR_FSI_BINS == 1
    static const int kNumOfSKDetFSIParamsMultiRing = 18;
    static const int kNumOfSKDetFSIParamsMultiRmu  = 6;
    static const int kNumOfSKDetFSIParamsMultiRe   = 12;
    
#elif MR_FSI_BINS == 2
    static const int kNumOfSKDetFSIParamsMultiRing = 6;
    static const int kNumOfSKDetFSIParamsMultiRmu  = 6;
    static const int kNumOfSKDetFSIParamsMultiRe   = 0;
    
// Have 12 parameters but empty covariances in the lowest energy bins
// 4 parameters with 3 bins each -> 12 - 4    
#elif MR_FSI_BINS == 3 
    static const int kNumOfSKDetFSIParamsMultiRing = 8; //12;
    static const int kNumOfSKDetFSIParamsMultiRmu  = 0;
    static const int kNumOfSKDetFSIParamsMultiRe   = 8; //12;    
#endif
    
#endif

    // NEW FHC+RHC JOINT
    static const int kNumOfFHCpRHCSKDetFSIParams     = kNumOfSKDetFSIParams + kNumOfRHCSKDetFSIParams; // 36
    static const int kNumOfFHCpRHCSKDetFSIParams1Rmu = kNumOfSKDetFSIParams1Rmu + kNumOfRHCSKDetFSIParams1Rmu; //12
    static const int kNumOfFHCpRHCSKDetFSIParams1Re  = kNumOfSKDetFSIParams1Re + kNumOfRHCSKDetFSIParams1Re; //24
    //

    // MultiRing + JOINT
    static const int kNumOfFHCpRHCSKDetFSIParamsMultiRing = kNumOfSKDetFSIParams + kNumOfRHCSKDetFSIParams + kNumOfSKDetFSIParamsMultiRing; // 52
    static const int kNumOfFHCpRHCSKDetFSIParamsMultiRmu  = kNumOfSKDetFSIParams + kNumOfRHCSKDetFSIParams + kNumOfSKDetFSIParamsMultiRmu; // 42
    static const int kNumOfFHCpRHCSKDetFSIParamsMultiRe   = kNumOfSKDetFSIParams + kNumOfRHCSKDetFSIParams + kNumOfSKDetFSIParamsMultiRe; // 48
  

    // The following is the list of the above kNumOfFSIParams
    static const PhysParamId_t kListOfSKDetFSIParams[kNumOfSKDetFSIParams] =
      {
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
#if FSI_BINS != 3
        kNuisSKDetFSI000,
        kNuisSKDetFSI001,
        kNuisSKDetFSI002,
        kNuisSKDetFSI003,
        kNuisSKDetFSI004,
        kNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
        kNuisSKDetFSI006,
        kNuisSKDetFSI007,
        kNuisSKDetFSI008,
        kNuisSKDetFSI009,
        kNuisSKDetFSI010,
        kNuisSKDetFSI011,
        kNuisSKDetFSI012,
        kNuisSKDetFSI013,
        kNuisSKDetFSI014,
        kNuisSKDetFSI015,
        kNuisSKDetFSI016,
        kNuisSKDetFSI017
#endif
#elif defined(__USE_SKDETFSI2015__)
#if FSI_BINS != 3
        kNuisSKDetFSI000,
        kNuisSKDetFSI001,
        kNuisSKDetFSI002,
        kNuisSKDetFSI003,
        kNuisSKDetFSI004,
        kNuisSKDetFSI005,
        kNuisSKDetFSI006,
        kNuisSKDetFSI007,
        kNuisSKDetFSI008,
        kNuisSKDetFSI009,
        kNuisSKDetFSI010,
        kNuisSKDetFSI011,
#endif
#if FSI_BINS != 2
        kNuisSKDetFSI012,
        kNuisSKDetFSI013,
        kNuisSKDetFSI014,
        kNuisSKDetFSI015,
        kNuisSKDetFSI016,
        kNuisSKDetFSI017,
        kNuisSKDetFSI018,
        kNuisSKDetFSI019,
        kNuisSKDetFSI020,
        kNuisSKDetFSI021,
        kNuisSKDetFSI022,
        kNuisSKDetFSI023
#endif
#endif
      };

    static const PhysParamId_t kListOfRHCSKDetFSIParams[kNumOfRHCSKDetFSIParams] =
      {
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
#if FSI_BINS != 3
        kRHCNuisSKDetFSI000,
        kRHCNuisSKDetFSI001,
        kRHCNuisSKDetFSI002,
        kRHCNuisSKDetFSI003,
        kRHCNuisSKDetFSI004,
        kRHCNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
        kRHCNuisSKDetFSI006,
        kRHCNuisSKDetFSI007,
        kRHCNuisSKDetFSI008,
        kRHCNuisSKDetFSI009,
        kRHCNuisSKDetFSI010,
        kRHCNuisSKDetFSI011,
        kRHCNuisSKDetFSI012,
        kRHCNuisSKDetFSI013,
        kRHCNuisSKDetFSI014,
        kRHCNuisSKDetFSI015,
        kRHCNuisSKDetFSI016,
        kRHCNuisSKDetFSI017
#endif

#elif defined(__USE_SKDETFSI2015__)
#if FSI_BINS != 3
        kRHCNuisSKDetFSI000,
        kRHCNuisSKDetFSI001,
        kRHCNuisSKDetFSI002,
        kRHCNuisSKDetFSI003,
        kRHCNuisSKDetFSI004,
        kRHCNuisSKDetFSI005,
        kRHCNuisSKDetFSI006,
        kRHCNuisSKDetFSI007,
        kRHCNuisSKDetFSI008,
        kRHCNuisSKDetFSI009,
        kRHCNuisSKDetFSI010,
        kRHCNuisSKDetFSI011,
#endif
#if FSI_BINS != 2
        kRHCNuisSKDetFSI012,
        kRHCNuisSKDetFSI013,
        kRHCNuisSKDetFSI014,
        kRHCNuisSKDetFSI015,
        kRHCNuisSKDetFSI016,
        kRHCNuisSKDetFSI017,
        kRHCNuisSKDetFSI018,
        kRHCNuisSKDetFSI019,
        kRHCNuisSKDetFSI020,
        kRHCNuisSKDetFSI021,
        kRHCNuisSKDetFSI022,
        kRHCNuisSKDetFSI023
#endif
#endif
      };


    // NEW FHC+RHC JOINT
    static const PhysParamId_t kListOfFHCpRHCSKDetFSIParams[kNumOfFHCpRHCSKDetFSIParams] =
      {
#if defined(__USE_SKDETFSI2013__)
#if FSI_BINS != 3
        kNuisSKDetFSI000,
        kNuisSKDetFSI001,
        kNuisSKDetFSI002,
        kNuisSKDetFSI003,
        kNuisSKDetFSI004,
        kNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
        kNuisSKDetFSI006,
        kNuisSKDetFSI007,
        kNuisSKDetFSI008,
        kNuisSKDetFSI009,
        kNuisSKDetFSI010,
        kNuisSKDetFSI011,
        kNuisSKDetFSI012,
        kNuisSKDetFSI013,
        kNuisSKDetFSI014,
        kNuisSKDetFSI015,
        kNuisSKDetFSI016,
        kNuisSKDetFSI017,
#endif
#elif defined(__USE_SKDETFSI2015__)
#if FSI_BINS != 3
        kNuisSKDetFSI000,
        kNuisSKDetFSI001,
        kNuisSKDetFSI002,
        kNuisSKDetFSI003,
        kNuisSKDetFSI004,
        kNuisSKDetFSI005,
        kNuisSKDetFSI006,
        kNuisSKDetFSI007,
        kNuisSKDetFSI008,
        kNuisSKDetFSI009,
        kNuisSKDetFSI010,
        kNuisSKDetFSI011,
#endif
#if FSI_BINS != 2
        kNuisSKDetFSI012,
        kNuisSKDetFSI013,
        kNuisSKDetFSI014,
        kNuisSKDetFSI015,
        kNuisSKDetFSI016,
        kNuisSKDetFSI017,
        kNuisSKDetFSI018,
        kNuisSKDetFSI019,
        kNuisSKDetFSI020,
        kNuisSKDetFSI021,
        kNuisSKDetFSI022,
        kNuisSKDetFSI023,
#endif
#endif

#if defined(__USE_SKDETFSI2013__)
#if FSI_BINS != 3
        kRHCNuisSKDetFSI000,
        kRHCNuisSKDetFSI001,
        kRHCNuisSKDetFSI002,
        kRHCNuisSKDetFSI003,
        kRHCNuisSKDetFSI004,
        kRHCNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
        kRHCNuisSKDetFSI006,
        kRHCNuisSKDetFSI007,
        kRHCNuisSKDetFSI008,
        kRHCNuisSKDetFSI009,
        kRHCNuisSKDetFSI010,
        kRHCNuisSKDetFSI011,
        kRHCNuisSKDetFSI012,
        kRHCNuisSKDetFSI013,
        kRHCNuisSKDetFSI014,
        kRHCNuisSKDetFSI015,
        kRHCNuisSKDetFSI016,
        kRHCNuisSKDetFSI017
#endif

#elif defined(__USE_SKDETFSI2015__)
#if FSI_BINS != 3
        kRHCNuisSKDetFSI000,
        kRHCNuisSKDetFSI001,
        kRHCNuisSKDetFSI002,
        kRHCNuisSKDetFSI003,
        kRHCNuisSKDetFSI004,
        kRHCNuisSKDetFSI005,
        kRHCNuisSKDetFSI006,
        kRHCNuisSKDetFSI007,
        kRHCNuisSKDetFSI008,
        kRHCNuisSKDetFSI009,
        kRHCNuisSKDetFSI010,
        kRHCNuisSKDetFSI011,
#endif
#if FSI_BINS != 2
        kRHCNuisSKDetFSI012,
        kRHCNuisSKDetFSI013,
        kRHCNuisSKDetFSI014,
        kRHCNuisSKDetFSI015,
        kRHCNuisSKDetFSI016,
        kRHCNuisSKDetFSI017,
        kRHCNuisSKDetFSI018,
        kRHCNuisSKDetFSI019,
        kRHCNuisSKDetFSI020,
        kRHCNuisSKDetFSI021,
        kRHCNuisSKDetFSI022,
        kRHCNuisSKDetFSI023
#endif
#endif
      };
    //



    /// NEW MULTIRING 2
    static const PhysParamId_t kListOfMultiRingSKDetFSIParams[kNumOfSKDetFSIParamsMultiRing] =
      {
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
#if MR_FSI_BINS != 3
        kMultiRingNuisSKDetFSI000,
        kMultiRingNuisSKDetFSI001,
        kMultiRingNuisSKDetFSI002,
        kMultiRingNuisSKDetFSI003,
        kMultiRingNuisSKDetFSI004,
        kMultiRingNuisSKDetFSI005,
#endif
#if MR_FSI_BINS != 2
        // "Masked" the parameters here as the lowest energy bins have no entries
        //kMultiRingNuisSKDetFSI006,
        kMultiRingNuisSKDetFSI007,
        kMultiRingNuisSKDetFSI008,
        //kMultiRingNuisSKDetFSI009,
        kMultiRingNuisSKDetFSI010,
        kMultiRingNuisSKDetFSI011,
        //kMultiRingNuisSKDetFSI012,
        kMultiRingNuisSKDetFSI013,
        kMultiRingNuisSKDetFSI014,
        //kMultiRingNuisSKDetFSI015,
        kMultiRingNuisSKDetFSI016,
        kMultiRingNuisSKDetFSI017
#endif
#endif
      };
    ////

    static const double kSKDetFSIParamsErecoCut[kNumOfSKDetFSIParams][2] =
      {
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
#if FSI_BINS != 3
        { 0.0,   0.4 }, // FSI1Rm000
        { 0.4,   1.1 }, // FSI1Rm001
        { 1.1,  30.0 }, // FSI1Rm002
        { 0.0,  30.0 }, // FSI1Rm003
        { 0.0,  30.0 }, // FSI1Rm004
        { 0.0,  30.0 }, // FSI1Rm005
#endif
#if FSI_BINS != 2
        { 0.0,  0.35 }, // FSI1Re000
        { 0.35,  0.8 }, // FSI1Re001
        { 0.8,  1.25 }, // FSI1Re002
        { 0.0,  0.35 }, // FSI1Re003
        { 0.35,  0.8 }, // FSI1Re004
        { 0.8,  1.25 }, // FSI1Re005
        { 0.0,  0.35 }, // FSI1Re006
        { 0.35,  0.8 }, // FSI1Re007
        { 0.8,  1.25 }, // FSI1Re008
        { 0.0,  0.35 }, // FSI1Re009
        { 0.35,  0.8 }, // FSI1Re010
        { 0.8,  1.25 }  // FSI1Re011
#endif

#elif defined(__USE_SKDETFSI2015__)
#if FSI_BINS != 3
        { 0.0,   0.4 }, // FSI1Rm000 - numu CCQE
        { 0.4,   1.1 }, // FSI1Rm001 - numu CCQE
        { 1.1,  30.0 }, // FSI1Rm002 - numu CCQE
        { 0.0,   0.4 }, // FSI1Rm003 - numu CC
        { 0.4,   1.1 }, // FSI1Rm004 - numu CC
        { 1.1,  30.0 }, // FSI1Rm005 - numu CC
        { 0.0,   0.4 }, // FSI1Rm006 - nue CC
        { 0.4,   1.1 }, // FSI1Rm007 - nue CC
        { 1.1,  30.0 }, // FSI1Rm008 - nue CC
        { 0.0,   0.4 }, // FSI1Rm009 - NC
        { 0.4,   1.1 }, // FSI1Rm010 - NC
        { 1.1,  30.0 }, // FSI1Rm011 - NC
#endif
#if FSI_BINS != 2
        { 0.0,  0.35 }, // FSI1Re000 - osc nue CCQE
        { 0.35,  0.8 }, // FSI1Re001 - osc nue CCQE
        { 0.8,  1.25 }, // FSI1Re002 - osc nue CCQE
        { 0.0,  0.35 }, // FSI1Re003 - numu CC
        { 0.35,  0.8 }, // FSI1Re004 - numu CC
        { 0.8,  1.25 }, // FSI1Re005 - numu CC
        { 0.0,  0.35 }, // FSI1Re006 - nue CC
        { 0.35,  0.8 }, // FSI1Re007 - nue CC
        { 0.8,  1.25 }, // FSI1Re008 - nue CC
        { 0.0,  0.35 }, // FSI1Re009 - NC
        { 0.35,  0.8 }, // FSI1Re010 - NC
        { 0.8,  1.25 }  // FSI1Re011 - NC
#endif
#endif
      };

    static const double kRHCSKDetFSIParamsErecoCut[kNumOfRHCSKDetFSIParams][2] =
      {
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
#if FSI_BINS != 3
        { 0.0,   0.4 }, // FSI1Rm000
        { 0.4,   1.1 }, // FSI1Rm001
        { 1.1,  30.0 }, // FSI1Rm002
        { 0.0,  30.0 }, // FSI1Rm003
        { 0.0,  30.0 }, // FSI1Rm004
        { 0.0,  30.0 }, // FSI1Rm005
#endif
#if FSI_BINS != 2
        { 0.0,  0.35 }, // FSI1Re000
        { 0.35,  0.8 }, // FSI1Re001
        { 0.8,  1.25 }, // FSI1Re002
        { 0.0,  0.35 }, // FSI1Re003
        { 0.35,  0.8 }, // FSI1Re004
        { 0.8,  1.25 }, // FSI1Re005
        { 0.0,  0.35 }, // FSI1Re006
        { 0.35,  0.8 }, // FSI1Re007
        { 0.8,  1.25 }, // FSI1Re008
        { 0.0,  0.35 }, // FSI1Re009
        { 0.35,  0.8 }, // FSI1Re010
        { 0.8,  1.25 }  // FSI1Re011
#endif

#elif defined(__USE_SKDETFSI2015__)
#if FSI_BINS != 3
        { 0.0,   0.4 }, // FSI1Rm000 - numu CCQE
        { 0.4,   1.1 }, // FSI1Rm001 - numu CCQE
        { 1.1,  30.0 }, // FSI1Rm002 - numu CCQE
        { 0.0,   0.4 }, // FSI1Rm003 - numu CC
        { 0.4,   1.1 }, // FSI1Rm004 - numu CC
        { 1.1,  30.0 }, // FSI1Rm005 - numu CC
        { 0.0,   0.4 }, // FSI1Rm006 - nue CC
        { 0.4,   1.1 }, // FSI1Rm007 - nue CC
        { 1.1,  30.0 }, // FSI1Rm008 - nue CC
        { 0.0,   0.4 }, // FSI1Rm009 - NC
        { 0.4,   1.1 }, // FSI1Rm010 - NC
        { 1.1,  30.0 }, // FSI1Rm011 - NC
#endif
#if FSI_BINS != 2
        { 0.0,  0.35 }, // FSI1Re000 - osc nue CCQE
        { 0.35,  0.8 }, // FSI1Re001 - osc nue CCQE
        { 0.8,  1.25 }, // FSI1Re002 - osc nue CCQE
        { 0.0,  0.35 }, // FSI1Re003 - numu CC
        { 0.35,  0.8 }, // FSI1Re004 - numu CC
        { 0.8,  1.25 }, // FSI1Re005 - numu CC
        { 0.0,  0.35 }, // FSI1Re006 - nue CC
        { 0.35,  0.8 }, // FSI1Re007 - nue CC
        { 0.8,  1.25 }, // FSI1Re008 - nue CC
        { 0.0,  0.35 }, // FSI1Re009 - NC
        { 0.35,  0.8 }, // FSI1Re010 - NC
        { 0.8,  1.25 }  // FSI1Re011 - NC
#endif
#endif
      };

    /// NEW MULTIRING 2
    static const double kMultiRingSKDetFSIParamsErecoCut[kNumOfSKDetFSIParamsMultiRing][2] =
      {
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
#if MR_FSI_BINS != 3
        { 0.0,   0.4 }, // FSI000
        { 0.4,   1.1 }, // FSI001
        { 1.1,  30.0 }, // FSI002
        { 0.0,  30.0 }, // FSI003
        { 0.0,  30.0 }, // FSI004
        { 0.0,  30.0 }, // FSI005
#endif
#if MR_FSI_BINS != 2
        // Extended second bin to the threshold for Delta resonance production ~340 MeV
        // lowest energy bins to allow for the low stats event
        // Edge is at 0.3 GeV to correspond with Ereco binning
        //{ 0.0,  0.35 }, // FSI1Re000
        //{ 0.35,  0.8 }, // FSI1Re001
        { 0.3,  0.8 }, // FSI1Re001
        { 0.8,  1.25 }, // FSI1Re002
        //{ 0.0,  0.35 }, // FSI1Re003
        //{ 0.35,  0.8 }, // FSI1Re004
        { 0.3,  0.8 }, // FSI1Re004
        { 0.8,  1.25 }, // FSI1Re005
        //{ 0.0,  0.35 }, // FSI1Re006
        //{ 0.35,  0.8 }, // FSI1Re007
        { 0.3,  0.8 }, // FSI1Re007
        { 0.8,  1.25 }, // FSI1Re008
        //{ 0.0,  0.35 }, // FSI1Re009
        //{ 0.35,  0.8 }, // FSI1Re010
        { 0.3,  0.8 }, // FSI1Re010
        { 0.8,  1.25 }  // FSI1Re011
#endif
#endif
      };
    ///


    /// NEW MULTIRING 2
#if FSI_BINS == 1 && MR_FSI_BINS == 3  // Joint fit

    //Todo: make this work if kNumOfSKDetFSIParams is different from kNumOfRHCSKDetFSIParams and kNumOfSKDetFSIParamsMultiRing

    static const bool kSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParams] = {

      // > 1 muon-like ring
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // 1Re:  oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        //
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1Gamma

#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1Gamma

#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1Gamma


        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1Gamma


        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif

      }//1Rmu
      ,
      // > 1 electron-like ring
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1Gamma

#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1Gamma

#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1Gamma


        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1Gamma

        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif

      }//1Re

    };



    static const bool kRHCSKDetFSIParamToSKMode[2][knm][kNumOfRHCSKDetFSIParams] = {

      // > 1 muon-like ring
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // 1Re:  oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        //
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1Gamma

#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1Gamma
#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1Gamma


        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1Gamma


        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif

      }//1Rmu
      ,
      // > 1 electron-like ring
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1Gamma
#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1Gamma

#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1Gamma

        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1Gamma

        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif
      }//1Re
    };


    static const bool kMultiRingSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParamsMultiRing] = {

      // > 1 muon-like ring MultiRing copy
      {
        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // 1Re:  oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
	// Removed lowest energy bins for the three categories above
        // Due to zeros (no events) in the lowest energy bins for the
        // covariance matrix
        //
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCQE
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC1pi
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCcoh
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuNC1piC
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuNCcoh
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuNCoth
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuNC1Gamma

#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarNCoth
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarNC1Gamma

#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuECCQE
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuECC1pi
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuECCcoh
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuECC2p2h
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuECCoth
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuENC1piC
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuENC1pi0
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuENCcoh
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuENCoth
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuENC1Gamma


        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuECCQE
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuECC1pi
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuECCcoh
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuECCoth

        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCQE
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCoth
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuEBarNCoth
        { 0, 0, 0, 1,   0, 0, 0, 0 }, // kSKmNuEBarNC1Gamma


        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        { 0, 0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif

      }//MultiRmu (1Rmu copy)
      ,
      // > 1 electron-like ring MultiRing copy
      //
      // oscnueCC x 2 / numuCC x 2 / nueCC x 2 / NC x 2
      {
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuCCQE
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuCC1pi
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuCCcoh
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuCCoth
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuNC1piC
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuNC1pi0
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuNCcoh
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuNCoth
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuNC1Gamma
#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuBarNC1piC
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuBarNC1pi0
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuBarNCcoh
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuBarNCoth
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuMuBarNC1Gamma
#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        { 0, 0, 1, 1,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuECCQE
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuECC1pi
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuECCcoh
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuECC2p2h
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuECCoth
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuENC1piC
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuENC1pi0
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuENCcoh
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuENCoth
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuENC1Gamma

        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuECCQE
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuECC1pi
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuECCcoh
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuECCoth

        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuEBarCCQE
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuEBarCC1pi
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuEBarCCcoh
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuEBarCC2p2h
        { 0, 0, 0, 0,   1, 1, 0, 0 }, // kSKmNuEBarCCoth
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuEBarNC1piC
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuEBarNC1pi0
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuEBarNCcoh
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuEBarNCoth
        { 0, 0, 0, 0,   0, 0, 1, 1 }, // kSKmNuEBarNC1Gamma

        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        { 1, 1, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
        { 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#  endif
      }//MultiRe (1Re copy)
    };

#elif FSI_BINS == 1 // Joint fit

    //Todo: make this work if kNumOfSKDetFSIParams is different from kNumOfRHCSKDetFSIParams and kNumOfSKDetFSIParamsMultiRing

    static const bool kSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParams] = {

      // > 1 muon-like ring
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // 1Re:  oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        //
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1Gamma

#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1Gamma

#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1Gamma


        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1Gamma


        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif

      }//1Rmu
      ,
      // > 1 electron-like ring
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1Gamma

#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1Gamma

#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1Gamma


        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1Gamma

        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif

      }//1Re

    };



    static const bool kRHCSKDetFSIParamToSKMode[2][knm][kNumOfRHCSKDetFSIParams] = {

      // > 1 muon-like ring
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // 1Re:  oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        //
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1Gamma

#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1Gamma
#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1Gamma


        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1Gamma


        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif

      }//1Rmu
      ,
      // > 1 electron-like ring
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1Gamma
#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1Gamma

#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1Gamma

        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1Gamma

        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif
      }//1Re
    };


    static const bool kMultiRingSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParamsMultiRing] = {

      // > 1 muon-like ring MultiRing copy
      {

        // NEW SKDETFSI

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // 1Re:  oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        //
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif	
#  endif
       	
	
      }//MultiRmu (1Rmu copy)
      ,
      // > 1 electron-like ring MultiRing copy 
      {

#  if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1Gamma
#    ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1Gamma        
#    ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#    endif
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1Gamma        

        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBar1Gamma

        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#    ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#    ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#    endif
#  endif

      }//MultiRe (1Re copy)
};
    
#elif FSI_BINS == 2 // numu disappearance

    //Todo: make this work if kNumOfSKDetFSIParams is different from kNumOfRHCSKDetFSIParams

    static const bool kSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParams] = {

      // > 1 muon-like ring
      {

        // NEW SKDETFSI

#if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#endif

      }//1Rmu
      ,
      // > 1 electron-like ring
      {

        // NEW SKDETFSI

#if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#endif
      }//1Re

    };


    static const bool kRHCSKDetFSIParamToSKMode[2][knm][kNumOfRHCSKDetFSIParams] = {

      // > 1 muon-like ring RHC
      {

        // NEW SKDETFSI

#if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#endif

      }//1RmuRHC
      ,
      // > 1 electron-like ring RHC
      {

        // NEW SKDETFSI

#if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#endif

      }//1ReRHC
    };


    static const bool kMultiRingSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParamsMultiRing] = {

      // > 1 muon-like ring RHC MultiRing copy
      {

        // NEW SKDETFSI

#if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#endif

      }//MultiRmu (1Rmu copy)
      ,
      // > 1 electron-like ring MultiRing copy
      {

        // NEW SKDETFSI

#if defined(__USE_SKDETFSI2013__)
        // 2013
        //
        // 1Rmu: CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
        //
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#endif
      }//MultiRe (1Re copy)

    };


#elif FSI_BINS == 3 // nue appearance

    // SAME PARAMETRIZATION FOR 2013 AND 2015 SKDETFSI PARAMETRIZATIONS

    //Todo: make this work if kNumOfSKDetFSIParams is different from kNumOfRHCSKDetFSIParams

    static const bool kSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParams] = {

      // > 1 muon-like ring
      //
      // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
      {
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

      }//1Rmu
      ,
      // > 1 electron-like ring
      //
      // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
      {
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth

        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth

        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
      }//1Re
    };


    static const bool kRHCSKDetFSIParamToSKMode[2][knm][kNumOfRHCSKDetFSIParams] = {
      // > 1 muon-like ring
      //
      // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
      {
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

      }//1Rmu
      ,
      // > 1 electron-like ring
      //
      // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
      {
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth

        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth

        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
      }//1Re
    };

    static const bool kMultiRingSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParamsMultiRing] = {
      // > 1 muon-like ring MultiRing copy
      //
      // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
      {
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth

        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

      }//MultiRmu (1Rmu copy)
      ,
      // > 1 electron-like ring MultiRing copy
      //
      // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
      {
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1Gamma        
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC2p2h
        {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC2p2h
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1Gamma        

        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC2p2h
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC2p2h
        {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1piC
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1Gamma        

        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC2p2h
        {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
        {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
      }//MultiRe (1Re copy)
    };

#endif


    //
    // Number of nuisance params parameters for 1Rmu, 1Re, 1Rmu+1Re analysis
    // 1Rmu : 57 nuisance params
    //   - 1 SK E scale
    //   - 7 SK efficiencies (3 CCQE-RC, CCQE-Oth, CCnonQE, NC, nue CC)
    //   - 31 BANFF flux
    //   - 8 BANFF XSec
    //   - 6 XSec parameters treated with splines
    //   - 5 XSec normalisations
    //   - 15 FSI
    // 1Re :
    //   - 25 BANFF flux
    //   -  8 BANFF xsec
    // Joint :
    //   - 25 BANFF flux
    //   -  8 BANFF xsec
    //

    // TODO
    // It doesn't support RHC systematic parameters if not fully correlated w/ FHC sample --> to fix it
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)

#if defined(__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__) && !defined(__RHC_RUNNING__)  // fhc all + skdetfsi multiring (no escale)
#  if defined(__MULTIRING_FHC_NUMU__) && ! defined(__MULTIRING_FHC_NUE__)
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRmu;  // Tot # par FHC+MultiRmu
#  elif defined (__MULTIRING_FHC_NUE__) && ! defined(__MULTIRING_FHC_NUMU__)
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRe;  // Tot # par FHC+MultiRmu
#  else
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRing;  // Tot # par FHC+MultiRing
#  endif
    // NEW FLUX
#elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__)// fhc+rhc (fully correlated)
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams + 1/*energy scale*/  + kNumOfUncXSecSyst ;
    //    
    
    // NEW FHC+RHC JOINT and multiring
#elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && !defined(__FULLY_CORRELATED__) // fhc+rhc (w/ correlations)
#  if defined(__MULTIRING_FHC_NUMU__) && ! defined(__MULTIRING_FHC_NUE__)
    static const int kNumOfNuisanceParams = kNumOfFHCpRHCParamsInBANFFfit + kNumOfFHCpRHCSKDetFSIParams + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRmu;  // Tot # par FHC+MultiRmu
#  elif defined (__MULTIRING_FHC_NUE__) && ! defined(__MULTIRING_FHC_NUMU__)
    static const int kNumOfNuisanceParams = kNumOfFHCpRHCParamsInBANFFfit + kNumOfFHCpRHCSKDetFSIParams + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRe;  // Tot # par FHC+MultiRmu
#  else
    static const int kNumOfNuisanceParams = kNumOfFHCpRHCParamsInBANFFfit + kNumOfFHCpRHCSKDetFSIParams + 1/*energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRing;
#  endif
    
      //

#elif defined(__FHC_RUNNING__)  // fhc only
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams + 1/*energy scale*/  + kNumOfUncXSecSyst ;

    // NEW FLUX
#elif defined(__RHC_RUNNING__) // rhc only
    static const int kNumOfNuisanceParams = kNumOfRHCParamsInBANFFfit + kNumOfRHCSKDetFSIParams + 1/*energy scale*/  + kNumOfUncXSecSyst ; // uncorr xsec is identical to fhc-only
    //
    
#endif

#else
#ifdef __MULTIRING_FHC_RUNNING__ // multiring fhc only
#  if defined(__MULTIRING_FHC_NUMU__) && ! defined(__MULTIRING_FHC_NUE__)
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRmu;  // Tot # par FHC+MultiRmu
#  elif defined (__MULTIRING_FHC_NUE__) && ! defined(__MULTIRING_FHC_NUMU__)
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRe;  // Tot # par FHC+MultiRmu
#  else
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRing;  // Tot # par FHC+MultiRing
#  endif

#else
    cerr << "No samples have been chosen for the analysis!!!" << endl;
    exit(1);
#endif
#endif // sample analysis option

  } // conv
} // valor

#endif

