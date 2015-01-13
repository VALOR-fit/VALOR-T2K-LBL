//____________________________________________________________________________
/*!

\class    valor skmodes

\brief    Superk Modes for VALOR

\author   Steve Dennis <stephen.dennis \at warwick.ac.uk>
          University of Warwick ; STFC, Rutherford Appleton Laboratory
          
\created  April 08 2014

Changes since 2014 v0.0.1:


*/
//____________________________________________________________________________

#ifndef _skmodes_modes_h_
#define _skmodes_modes_h_

#include "Enumerations.h"
#include "AnalysisOptions.h"

namespace valor {
  
namespace conv   {
  
// Enumeration of (true) reaction modes that may be considered separately within each SuperK sample
// Actual choice for current analysis is given in `Configuration' section.
static const int kNSuperkModes = 56;  // number of reaction modes kept track of within each sample
static const int knm = kNSuperkModes; // as above - shorter name
typedef enum ESuperkMode 
{
   kSKmUndefined = -1, 
   kSKmNuMuCCQE,         // numu CC QE
   kSKmNuMuCC1pi,        // numu CC 1pi
   kSKmNuMuCCcoh,        // numu CC coherent
   kSKmNuMuCCMEC,        // numu CC MEC
   kSKmNuMuCCoth,        // numu CC other
   kSKmNuMuNC1piC,       // numu NC 1 charged pi
   kSKmNuMuNC1pi0,       // numu NC 1pi0
   kSKmNuMuNCcoh,        // numu NC coherent
   kSKmNuMuNCoth,        // numu NC other
   kSKmOscNuMuCCQE,      // oscillated numu CCQE
   kSKmOscNuMuCC1pi,     // oscillated numu CC 1pi
   kSKmOscNuMuCCcoh,     // oscillated numu coh
   kSKmOscNuMuCCMEC,     // oscillated numu CC MEC
   kSKmOscNuMuCCoth,     // oscillated numu CC other
   kSKmNuMuBarCCQE,      // numubar CCQE
   kSKmNuMuBarCC1pi,     // numubar CC 1pi
   kSKmNuMuBarCCcoh,     // numubar CC coh
   kSKmNuMuBarCCMEC,     // numubar CC MEC
   kSKmNuMuBarCCoth,     // numubar CC other
   kSKmNuMuBarNC1piC,    // numubar NC 1 charged pi
   kSKmNuMuBarNC1pi0,    // numubar NC 1pi0
   kSKmNuMuBarNCcoh,     // numubar NC coherent
   kSKmNuMuBarNCoth,     // numubar NC other
   kSKmOscNuMuBarCCQE,   // oscillated numubar CCQE
   kSKmOscNuMuBarCC1pi,  // oscillated numubar CC 1pi
   kSKmOscNuMuBarCCcoh,  // oscillated numubar CC coh
   kSKmOscNuMuBarCCMEC,  // oscillated numubar CC MEC
   kSKmOscNuMuBarCCoth,  // oscillated numubar CC other
   kSKmNuECCQE,          // nue CC QE
   kSKmNuECC1pi,         // nue CC 1pi
   kSKmNuECCcoh,         // nue CC coherent
   kSKmNuECCMEC,         // nue CC MEC
   kSKmNuECCoth,         // nue CC other
   kSKmNuENC1piC,        // nue NC 1 charged pi
   kSKmNuENC1pi0,        // nue NC 1pi0
   kSKmNuENCcoh,         // nue NC coherent
   kSKmNuENCoth,         // nue NC other, nu_mu appearance analysis definition
   kSKmOscNuEBarCCQE,    // oscillated nuebar CC QE
   kSKmOscNuEBarCC1pi,   // oscillated nuebar CC 1pi
   kSKmOscNuEBarCCcoh,   // oscillated nuebar CC coherent
   kSKmOscNuEBarCCMEC,   // oscillated nuebar CC MEC
   kSKmOscNuEBarCCoth,   // oscillated nuebar CC other
   kSKmNuEBarCCQE,       // nuebar CC QE
   kSKmNuEBarCC1pi,      // nuebar CC 1pi
   kSKmNuEBarCCcoh,      // nuebar CC coherent
   kSKmNuEBarCCMEC,      // nuebar CC MEC
   kSKmNuEBarCCoth,      // nuebar CC other
   kSKmNuEBarNC1piC,     // nuebar NC 1 charged pi
   kSKmNuEBarNC1pi0,     // nuebar NC 1pi0
   kSKmNuEBarNCcoh,      // nuebar NC coherent
   kSKmNuEBarNCoth,      // nuebar NC other
   kSKmOscNuECCQE,       // oscillated nue CCQE
   kSKmOscNuECC1pi,      // oscillated nue CC 1pi
   kSKmOscNuECCcoh,      // oscillated nue CC coh
   kSKmOscNuECCMEC,      // oscillated nue CC MEC
   kSKmOscNuECCoth,      // oscillated nue CC other
} 
SuperkMode_t; 
} // conv

namespace util {
  // Modes
  int SuperkMode2FluxNuPdg(conv::SuperkMode_t m);
  int SuperkMode2SKMCPdg  (conv::SuperkMode_t m);

  bool IsNuMu             (conv::SuperkMode_t  m);
  bool IsOscNuMu          (conv::SuperkMode_t  m);
  bool IsNuMuBar          (conv::SuperkMode_t  m);
  bool IsOscNuMuBar       (conv::SuperkMode_t  m);
  bool IsNuE              (conv::SuperkMode_t  m);
  bool IsNuEBar           (conv::SuperkMode_t  m);
  bool IsOscNuEBar        (conv::SuperkMode_t  m);
  bool IsOscNuE           (conv::SuperkMode_t  m);
  bool IsCC               (conv::SuperkMode_t  m);
  bool IsNC               (conv::SuperkMode_t  m);
  bool IsMEC              (conv::SuperkMode_t  m);
  
  // String stuff
  conv::SuperkMode_t  SuperkModeFromString      (string h);
  string AsString           (conv::SuperkMode_t       m);
  string AsStringCompact    (conv::SuperkMode_t       m);
  string AsStringTex        (conv::SuperkMode_t       m);
  string AsStringRooTex     (conv::SuperkMode_t       m);
  // Aliases
  string ModeAsString           (conv::SuperkMode_t       m);
  string ModeAsStringCompact    (conv::SuperkMode_t       m);
  string ModeAsStringTex        (conv::SuperkMode_t       m);
  string ModeAsStringRooTex     (conv::SuperkMode_t       m);
  
} // util
} // valor
  
#endif 

