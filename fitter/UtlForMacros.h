//____________________________________________________________________________
/*!

  \class     valor::utils

  \brief     Utilities for macros

  \author    Tom Dealtry <t.dealtry \at lancaster.ac.uk>
             Lancaster University
          
  \created   Aug 12, 2015

  \copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _valor_macros_utils_h_
#define _valor_macros_utils_h_

#include "OscParams.h"

#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#include <string>

namespace valor {
  class SuperkExpectation;
  class SuperkObservation;
  class PhysHypothesis;
  class FitMarginalization;

namespace util {
  
  // Get the string related to the systematics which are being tweaked
  std::string SysThrowIdToString(unsigned int tweak_systematics);
  // Tweak the relevant systematics
  void PerformSysThrow(unsigned int tweak_systematics, SuperkExpectation * skexp, bool print = true);
  //TODO should we change the uint to an enum SysTweakId_t?

  void PerformMargFit(int argc, char ** argv, FitMarginalization *fitmarg, SuperkExpectation *skexp, SuperkObservation *skobs_fakedata, bool print = true); // take inputs from command line
  void MargFit_1D_deltaCP(FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);

  // Construct an absolute file path, relative to $T2KOSCANA
  const char * FormFilePath(const char * file, const char * dir);

  // Set the mass hierarchy in Opt and SKMC (for both bar & non-bar) using either
  // - the current value in Opt
  // - a prior
  void SetMassHierarchy(SuperkExpectation * skexp);

  // Helper function. Get a TH2D from a file & check its valid
  TH2D * GetTH2DFromTFile(TFile * f, const char * name);
  
  conv::PhysParamId_t GetFreePar(string freeparname);

} // util
} // valor

#endif
