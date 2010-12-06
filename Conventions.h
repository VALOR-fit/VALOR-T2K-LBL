//____________________________________________________________________________
/*!

\class    t2kana::yr1numu::conv

\brief    Conventions

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>
          
\created  Oct 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_year1_numu_analysis_conv_h_
#define _t2kana_year1_numu_analysis_conv_h_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace t2kana {
namespace yr1numu {
namespace conv {

// confidence level
typedef enum EConfLevel {
   kCLUndefined = -1, 
   kCL68, 
   kCL90, 
   kCL95, 
   kCL99
} ConfLevel_t; 

// samples included in analysis
typedef enum ESuperkSample {
   kSKsUndefined = -1, 
   kSKsFCFV1RingMuLike, 
   kSKsFCFVmRingMuLike
} SuperkSample_t; 

// (true) reaction modes considered separately within each SuperK sample
typedef enum ESuperkMode {
   kSKmUndefined = -1, 
   kSKmNuMuCCQE, 
   kSKmNuMuCC1pi, 
   kSKmNuMuCCoth, 
   kSKmNC,    
   kSKmNuMuBar, 
   kSKmNuE
} SuperkMode_t; 

// (true) reaction modes considered separately in nd280 inputs
typedef enum END280Mode {
   kNDmUndefined = -1, 
   kNDmNuMuCCQE, 
   kNDmNuMuCCRES, 
   kNDmNuMuCCCOH, 
   kNDmNuMuCCOth, 
   kNDmNuMuNCEL, 
   kNDmNuMuNCRES, 
   kNDmNuMuNCCOH, 
   kNDmNuMuNCOth, 
   kNDmNuMuBarCCQE, 
   kNDmNuMuBarCCRES, 
   kNDmNuMuBarCCCOH, 
   kNDmNuMuBarCCOth, 
   kNDmNuMuBarNCEL, 
   kNDmNuMuBarNCRES, 
   kNDmNuMuBarNCCOH, 
   kNDmNuMuBarNCOth, 
   kNDmNuECCQE, 
   kNDmNuECCRES, 
   kNDmNuECCCOH, 
   kNDmNuECCOth, 
   kNDmNuENCEL, 
   kNDmNuENCRES, 
   kNDmNuENCCOH, 
   kNDmNuENCOth
} ND280Mode_t; 

// fit parameters
typedef enum EPhysParamId {   
   kPhypUndefined = -1, 
   kPhypSinsq2Theta23, 
   kPhypDm23sq,
   kPhypNorm,
} PhysParamId_t;


static const int kPdgNuE      =  12;
static const int kPdgNuMu     =  14;
static const int kPdgNuTau    =  16;
static const int kPdgNuEBar   = -12;
static const int kPdgNuMuBar  = -14;
static const int kPdgNuTauBar = -16;

static const int kNSuperkSamples  = 1;  // number of samples included in analysis
static const int kNSuperkModes    = 6;  // number of reaction modes considered separately within each SuperK sample
static const int kNND280Modes     = 24; // number of reaction modes considered separately within the ND280 sample
//shorter names
static const int kns   = kNSuperkSamples;
static const int knm   = kNSuperkModes; 
static const int knmnd = kNND280Modes; 

static const double kASmallNum   = 1E-10;
static const double kBaseline    = 295;        // Tokai->Kamioka (km)

static const double kEtrueMin    = 0;
static const double kEtrueMax    = 10;
static const double kEtrueBinSz  = 0.05;
static const double kNEtrueBins  = (int) ((kEtrueMax - kEtrueMin)/kEtrueBinSz);

// fixed bins: 50 MeV from 0 to 4 GeV
static const int kNvErecoBins = 80;
static const double kNvErecoBinEdges[kNvErecoBins+1] = {
   0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 
   0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 
   1.000, 1.050, 1.100, 1.150, 1.200, 1.250, 1.300, 1.350, 1.400, 1.450, 
   1.500, 1.550, 1.600, 1.650, 1.700, 1.750, 1.800, 1.850, 1.900, 1.950, 
   2.000, 2.050, 2.100, 2.150, 2.200, 2.250, 2.300, 2.350, 2.400, 2.450, 
   2.500, 2.550, 2.600, 2.650, 2.700, 2.750, 2.800, 2.850, 2.900, 2.950, 
   3.000, 3.050, 3.100, 3.150, 3.200, 3.250, 3.300, 3.350, 3.400, 3.450, 
   3.500, 3.550, 3.600, 3.650, 3.700, 3.750, 3.800, 3.850, 3.900, 3.950, 
   4.000
};

// variable bins
/*
static const int kNvErecoBins = 7;
static const double kNvErecoBinEdges[kNvErecoBins+1] = {
   0.000, 0.300, 0.450, 0.600, 0.750, 1.000, 1.500, 4.000 
};
*/

} // conv
} // yr1numu
} // t2kana
  
#endif 


