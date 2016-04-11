//~ #include "SuperkObservation.h"
#include "SuperkExpectation.h"
//~ #include "T2KPhysicalParameterLimits.h"

#include "Utilities.h"
//~ #include "Binning.h"
//~ #include "AnalysisOptions.h"

#include "Build.h"
//~ #include "UtlForMacrosT2XK.h"

//~ #include "Utilities/UtlForMacros.h"
//~ #include "Utilities/OscParams.h"
//~ #include "Utilities/RunOpt.h"
#include "Utilities/Env.h"
//~ #include "Conventions/PhysHypothesisId.h"

//~ #include "TROOT.h"
//~ #include "TRandom3.h"
#include "TSystem.h"
//~ #include "TStopwatch.h"
#include "TString.h"

//~ #include <sstream>
//~ #include <fstream>
//~ #include <iomanip>

using namespace valor;

int main(/*int argc, char ** argv*/)
{
  Env::Instance()->AutoSet("t2xk");
    
  const char * cfile = "analysis_configuration.xml";
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), cfile ));
  
  SuperkExpectation skexp;
  skexp.Read(config.Data());
}


