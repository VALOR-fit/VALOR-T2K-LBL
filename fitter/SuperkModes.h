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
static const int kNSuperkModes = 53;  // number of reaction modes kept track of within each sample
static const int knm = kNSuperkModes; // as above - shorter name
typedef enum ESuperkMode 
{
   kSKmUndefined = -1, 
   kSKmNuMuCCQE,         // numu CC QE
   kSKmNuMuCC1pi,        // numu CC 1pi
   kSKmNuMuCCcoh,        // numu CC coherent
   kSKmNuMuCCoth,        // numu CC other
   kSKmNuMuNC1piC,       // numu NC 1 charged pi
   kSKmNuMuNC1pi0,       // numu NC 1pi0
   kSKmNuMuNCcoh,        // numu NC coherent
   kSKmNuMuNCoth_mu,     // numu NC other, nu_mu appearance analysis definition
   kSKmNuMuNCoth_e,      // numu NC other, nu_e appearance analysis definition
   kSKmOscNuMuCCQE,      // oscillated numu CCQE
   kSKmOscNuMuCC1pi,     // oscillated numu CC 1pi
   kSKmOscNuMuCCcoh,     // oscillated numu coh
   kSKmOscNuMuCCoth,     // oscillated numu CC other
   kSKmNuMuBarCCQE,      // numubar CCQE
   kSKmNuMuBarCC1pi,     // numubar CC 1pi
   kSKmNuMuBarCCcoh,     // numubar CC coh
   kSKmNuMuBarCCoth,     // numubar CC other
   kSKmNuMuBarNC1piC,    // numubar NC 1 charged pi
   kSKmNuMuBarNC1pi0,    // numubar NC 1pi0
   kSKmNuMuBarNCcoh,     // numubar NC coherent
   kSKmNuMuBarNCoth_mu,  // numubar NC other, nu_mu appearance analysis definition
   kSKmNuMuBarNCoth_e,   // numubar NC other, nu_e appearance analysis definition
   kSKmOscNuMuBarCCQE,   // oscillated numubar CCQE
   kSKmOscNuMuBarCC1pi,  // oscillated numubar CC 1pi
   kSKmOscNuMuBarCCcoh,  // oscillated numubar CC coh
   kSKmOscNuMuBarCCoth,  // oscillated numubar CC other
   kSKmNuECCQE,          // nue CC QE
   kSKmNuECC1pi,         // nue CC 1pi
   kSKmNuECCcoh,         // nue CC coherent
   kSKmNuECCoth,         // nue CC other
   kSKmNuENC1piC,        // nue NC 1 charged pi
   kSKmNuENC1pi0,        // nue NC 1pi0
   kSKmNuENCcoh,         // nue NC coherent
   kSKmNuENCoth_mu,      // nue NC other, nu_mu appearance analysis definition
   kSKmNuENCoth_e,       // nue NC other, nu_e appearance analysis definition
   kSKmOscNuEBarCCQE,    // oscillated nuebar CC QE
   kSKmOscNuEBarCC1pi,   // oscillated nuebar CC 1pi
   kSKmOscNuEBarCCcoh,   // oscillated nuebar CC coherent
   kSKmOscNuEBarCCoth,   // oscillated nuebar CC other
   kSKmNuEBarCCQE,       // nuebar CC QE
   kSKmNuEBarCC1pi,      // nuebar CC 1pi
   kSKmNuEBarCCcoh,      // nuebar CC coherent
   kSKmNuEBarCCoth,      // nuebar CC other
   kSKmNuEBarNC1piC,     // nuebar NC 1 charged pi
   kSKmNuEBarNC1pi0,     // nuebar NC 1pi0
   kSKmNuEBarNCcoh,      // nuebar NC coherent
   kSKmNuEBarNCoth_mu,   // nuebar NC other, nu_mu appearance analysis definition
   kSKmNuEBarNCoth_e,    // nuebar NC other, nu_e appearance analysis definition
   kSKmOscNuECCQE,       // oscillated nue CCQE
   kSKmOscNuECC1pi,      // oscillated nue CC 1pi
   kSKmOscNuECCcoh,      // oscillated nue CC coh
   kSKmOscNuECCoth,      // oscillated nue CC other
   kSKmMEC               // numu CC MEC
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
  
} // util
} // valor
  
#endif 

