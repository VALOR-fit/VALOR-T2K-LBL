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
   kSKmNuMuNonCCQE, 
   kSKmNC    
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


static const int kNSuperkSamples  = 2;  // number of samples included in analysis
static const int kNSuperkModes    = 3;  // number of reaction modes considered separately within each SuperK sample
static const int kNND280Modes     = 24; // number of reaction modes considered separately within the ND280 sample
//shorter names
static const int kns   = kNSuperkSamples;
static const int knm   = kNSuperkModes; 
static const int knmnd = kNND280Modes; 

static const double kASmallNum = 1E-10;
static const double kBaseline  = 295;        // Tokai->Kamioka (km)
static const double kNA        = 6.023E+23;  // Avogadro's number

} // conv
} // yr1numu
} // t2kana
  
#endif 


