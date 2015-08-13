//____________________________________________________________________________
/*!

  \class     valor::utils

  \brief     Utilities for macros

  \author    Tom Dealtry <t.dealtry \at lancaster.ac.uk>
  University of Lancaster
          
  \created   Aug 12, 2015

  \copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _valor_macros_utils_h_
#define _valor_macros_utils_h_

#include "Conventions/NuMassHierarchy.h"

class TH2D;

#include <string>

namespace valor {
  class SuperkExpectation;
  class PhysHypothesis;

namespace util {

  struct OscParams {
    //NH
    double mixing12_nh;
    double mixing13_nh;
    double mixing23_nh;
    double dcp_nh;
    double dm12sq_nh;
    double dm23sq_nh;
    //IH
    double mixing12_ih;
    double mixing13_ih;
    double mixing23_ih;
    double dcp_ih;
    double dm12sq_ih;
    double dm23sq_ih;
    //Other
    double beta_bar;
    double mh_nh_prior;
    int    iMH;
    conv::NuMassHierarchy_t MH;
    //Randomize?
    //NH
    bool randomize_mixing12_nh;
    bool randomize_mixing13_nh;
    bool randomize_mixing23_nh;
    bool randomize_dcp_nh;
    bool randomize_dm12sq_nh;
    bool randomize_dm23sq_nh;
    //IH
    bool randomize_mixing12_ih;
    bool randomize_mixing13_ih;
    bool randomize_mixing23_ih;
    bool randomize_dcp_ih;
    bool randomize_dm12sq_ih;
    bool randomize_dm23sq_ih;
    //Other
    bool randomize_beta_bar;
    bool marginalize_mh;
  };
  
  // Get the string related to the systematics which are being tweaked
  std::string SysThrowIdToString(unsigned int tweak_systematics);
  // Tweak the relevant systematics
  void PerformSysThrow(unsigned int tweak_systematics, SuperkExpectation * skexp,
		       bool print = true);
  //TODO should we change the uint to an enum SysTweakId_t?

  // Construct an absolute file path, relative to $T2KOSCANA
  const char * FormFilePath(const char * file, const char * dir);

  // Setup the entire Opt class in a single line
  void SetGlobalOptions(const char * dmsqscheme, const char * hierarchy,
			const char * mixparamchoice, const char * mixparamoctant,
			const char * cpvchoice);

  // Check ALL oscillation parameters for -99 (set to global best fits) and -999 (randomize) in a single function
  //  Creates an OscParams struct that should be used for everything else in the macro
  OscParams SetOscParamOptions(double input_mixing12_nh, double input_mixing13_nh, double input_mixing23_nh,
			       double input_deltacp_nh,  double input_dm12sq_nh,   double input_dm23sq_nh,
			       double input_beta_bar,
			       double input_mixing12_ih, double input_mixing13_ih, double input_mixing23_ih,
			       double input_deltacp_ih,  double input_dm12sq_ih,   double input_dm23sq_ih,
			       double mh_nh_prior,       bool marginalize_MH);

  // Randomize ALL oscillation parameters (if requested) in a single function
  //  Input values passed by reference, so can be updated
  void RandomizeOscParams(OscParams & inputs);
  
  // Update all the values in the PhysHypothesis, for any 2/3 flavour hypothesis type, in a single function
  //  PhysHypothesis passed by reference, so can be updated
  void SetOscParams(PhysHypothesis& ph, OscParams inputs, bool print = true);

  // Set the mass hierarchy in Opt and SKMC (for both bar & non-bar) using either
  // - the current value in Opt
  // - a prior
  void SetMassHierarchy(SuperkExpectation * skexp, OscParams & inputs);

  // Get the prior histograms from the file
  //  and create the 'values used in toys' histograms
  void GetPriorHists(const char * file_priors,
		     TH2D ** hPrior_th13_dcp, TH2D ** hPrior_23,
		     TH2D ** hToy_th13_dcp, TH2D ** hToy_23);
 
} // util
} // valor
  
#endif
