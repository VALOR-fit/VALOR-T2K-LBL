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
Add samples and modes for reverse horn current running

21-7-2014 Davide
Added samples for multiring forward horn current

*/
//____________________________________________________________________________

#ifndef _sksamples_h_
#define _sksamples_h_

#include "Enumerations.h"
#include "AnalysisOptions.h"

namespace valor {
namespace conv   {

// Enumeration of event samples that may be included in the analysis
// Actual choice for current analysis is given in `Configuration' section.
static const int kNSuperkSamples  = 6;  // number of samples 
static const int kns = kNSuperkSamples; // as above - shorter name
typedef enum ESuperkSample 
{
   kSKsUndefined = -1, 
   kSKsFCFV1RingMuLike, // fully-contained fiducial-volume single-ring mu-like
   kSKsFCFV1RingELike,  // fully-contained fiducial-volume single-ring e-like
   kSKsRHCFCFV1RingMuLike, //fully-contained fiducial-volume single-ring mu-like with reverse horn current 
   kSKsRHCFCFV1RingELike,  //fully-contained fiducial-volume single-ring e-like with reverse horn current 
   kSKsFCFVMultiRingMuLike, // fully-contained fiducial-volume multi-ring mu-like
   kSKsFCFVMultiRingELike  // fully-contained fiducial-volume multi-ring e-like
}
SuperkSample_t;

} // conv

namespace util {
  // Samples
  bool IsRHCSample        (conv::SuperkSample_t s);
  bool IsFHCSample        (conv::SuperkSample_t s);
  bool IsMultiRingSample  (conv::SuperkSample_t s);
  bool Is1RMuSample       (conv::SuperkSample_t s);
  bool Is1ReSample        (conv::SuperkSample_t s);
  bool IsMultiRMuSample   (conv::SuperkSample_t s);
  bool IsMultiReSample    (conv::SuperkSample_t s);
  
  string GetXMLMCTreeName  (conv::SuperkSample_t     s);
  string GetXMLDatasetName (conv::SuperkSample_t     s);
  
  // String stuff
  conv::SuperkSample_t  SKSampleFromString        (string h);
  string AsString           (conv::SuperkSample_t     s);
  string AsStringCompact    (conv::SuperkSample_t     s);
  string AsStringTex        (conv::SuperkSample_t     s);
  string AsStringRooTex     (conv::SuperkSample_t     s);
  // Aliases
  string SampleAsString       (conv::SuperkSample_t   s);
  string SampleAsStringCompact(conv::SuperkSample_t   s);
  string SampleAsStringTex    (conv::SuperkSample_t   s);
  string SampleAsStringRooTex (conv::SuperkSample_t   s);

  // Get conv::SuperkSample_t enum in array 
  conv::SuperkSample_t GetSKSample (int s);
  
} // util
} // valor
  
#endif 

