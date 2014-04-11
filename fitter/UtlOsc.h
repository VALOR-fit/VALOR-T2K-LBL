//____________________________________________________________________________
/*!

\class     t2kana::utils

\brief     Oscillation calculation utilities

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _valor_osc_util_h_
#define _valor_osc_util_h_

#include "Enumerations.h"

namespace t2kana {
namespace util {

  //
  // Simple physics hypotheses (effective 2-flavour oscillation, decay and decoherence).
  // Code for 3-flavour oscillations incl. matter effects is in NuOsc3Flavour.cxx
  //

  double NuMuSurvProbOsc2Flavour (double mixing23, double dm23sq, double Et, bool unphysical_osc_prob);
  double NuMuSurvProbOsc2Flavour (double mixing23, double dm23sq, double Et, conv::SuperkMode_t m, bool unphysical_osc_prob);
  double NuMuSurvProbDecay       (double sinsq2thetadecay, double alpha, double Et);
  double NuMuSurvProbDecay       (double sinsq2thetadecay, double alpha, double Et, conv::SuperkMode_t m);
  double NuMuSurvProbDecoh       (double sinsq2thetadecoh, double musq, double Et);
  double NuMuSurvProbDecoh       (double sinsq2thetadecoh, double musq, double Et, conv::SuperkMode_t m);

} //util
} //t2kana

#endif
