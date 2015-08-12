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

#include "SuperkExpectation.h"
#include "FitCore/PhysHypothesis.h"

namespace valor {
namespace util {

  // Get the string related to the systematics which are being tweaked
  string SysThrowIdToString(unsigned int tweak_systematics);
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
  //  Everything passed by reference, so all variables can be updated
  void SetOscParamOptions(double& input_mixing12_nh, double& input_mixing13_nh, double& input_mixing23_nh,
			  double& input_deltacp_nh,  double& input_dm12sq_nh,   double& input_dm23sq_nh,
			  bool& randomize_mixing12_nh, bool& randomize_mixing13_nh, bool& randomize_mixing23_nh,
			  bool& randomize_deltacp_nh,  bool& randomize_dm12sq_nh,   bool& randomize_dm23sq_nh,
			  double& input_beta_bar, bool& randomize_beta_bar,
			  double& input_mixing12_ih, double& input_mixing13_ih, double& input_mixing23_ih,
			  double& input_deltacp_ih,  double& input_dm12sq_ih,   double& input_dm23sq_ih,
			  bool& randomize_mixing12_ih, bool& randomize_mixing13_ih, bool& randomize_mixing23_ih,
			  bool& randomize_deltacp_ih,  bool& randomize_dm12sq_ih,   bool& randomize_dm23sq_ih
			  );

    // Randomize ALL oscillation parameters (if requested) in a single function
    //  Input values passed by reference, so can be updated
    void RandomizeOscParams(double& input_mixing12_nh, double& input_mixing13_nh, double& input_mixing23_nh,
			    double& input_deltacp_nh,  double& input_dm12sq_nh,   double& input_dm23sq_nh,
			    bool randomize_mixing12_nh, bool randomize_mixing13_nh, bool randomize_mixing23_nh,
			    bool randomize_deltacp_nh,  bool randomize_dm12sq_nh,   bool randomize_dm23sq_nh,
			    double& input_beta_bar, bool randomize_beta_bar,
			    double& input_mixing12_ih, double& input_mixing13_ih, double& input_mixing23_ih,
			    double& input_deltacp_ih,  double& input_dm12sq_ih,   double& input_dm23sq_ih,
			    bool randomize_mixing12_ih, bool randomize_mixing13_ih, bool randomize_mixing23_ih,
			    bool randomize_deltacp_ih,  bool randomize_dm12sq_ih,   bool randomize_dm23sq_ih
			    );
  
    // Update all the values in the PhysHypothesis, for any 2/3 flavour hypothesis type, in a single function
    //  PhysHypothesis passed by reference, so can be updated
    void SetOscParams(PhysHypothesis& ph,
		      double input_mixing12_used, double input_mixing13_used, double input_mixing23_used,
		      double input_deltacp_used,  double input_dm12sq_used,   double input_dm23sq_used,
		      double input_beta_bar_used,
		      bool print = true
		      );
    
    
} // util
} // valor
  
#endif
