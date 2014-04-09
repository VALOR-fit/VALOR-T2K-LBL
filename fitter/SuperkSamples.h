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

31-12-2013  Nick
Add samples and modes for reverse horn current running

*/
//____________________________________________________________________________

#ifndef __skmodes_samples_h__
#define __skmodes_samples_h__

#include "Enumerations.h"
#include "AnalysisOptions.h"

namespace t2kana {
namespace conv   {

// Enumeration of event samples that may be included in the analysis
// Actual choice for current analysis is given in `Configuration' section.
static const int kNSuperkSamples  = 4;  // number of samples 
static const int kns = kNSuperkSamples; // as above - shorter name
typedef enum ESuperkSample 
{
   kSKsUndefined = -1, 
   kSKsFCFV1RingMuLike, // fully-contained fiducial-volume single-ring mu-like
   kSKsFCFV1RingELike,  // fully-contained fiducial-volume single-ring e-like
   kSKsRHCFCFV1RingMuLike, //fully-contained fiducial-volume single-ring mu-like with reverse horn current 
   kSKsRHCFCFV1RingELike  //fully-contained fiducial-volume single-ring e-like with reverse horn current 
} 
SuperkSample_t;

} // conv

namespace util {
  // Samples
  bool IsRHCSample        (conv::SuperkSample_t s);
  bool IsFHCSample        (conv::SuperkSample_t s);
  bool Is1RMuSample       (conv::SuperkSample_t s);
  bool Is1RESample        (conv::SuperkSample_t s);
  
  // String stuff
  conv::SuperkSample_t  SKSampleFromString        (string h);
  string AsString           (conv::SuperkSample_t     s);
  string AsStringCompact    (conv::SuperkSample_t     s);
  string AsStringTex        (conv::SuperkSample_t     s);
  string AsStringRooTex     (conv::SuperkSample_t     s);
  
} // util
} // t2kana
  
#endif 

