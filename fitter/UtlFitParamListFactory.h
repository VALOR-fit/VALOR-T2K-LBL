//____________________________________________________________________________
/*!

\class     valor::util::FitParamListFactory

\brief     List of parameters included in the fit

\author    Steve Dennis

\created   May 07, 2014

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _fit_param_list_factory_h_
#define _fit_param_list_factory_h_

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "Build.h"

#include "FitCore/FitParamList.h"

using std::string;
using std::map;

namespace valor {
  
  namespace util
  {
    FitParamList * NamedParameterList (string name);
    void AppendAllSystematicsTo(FitParamList * fp, int start);
    
    void  AppendAllUncorrelatedXSecSystematics       (FitParamList * fp, int start);
    void  AppendAllSuperkDetectorFSISystematics      (FitParamList * fp, int start);
    void  AppendAllCorrelatedBANFFSystematics        (FitParamList * fp, int start);

#ifdef __NUMUONLY__
    void  AppendUncorrelatedXSecSystematicsDisap     (FitParamList * fp, int start);
    void  AppendCorrelatedBANFFSystematicsDisap      (FitParamList * fp, int start);
    void  AppendSuperkDetectorFSISystematicsDisap    (FitParamList * fp, int start);
#endif

#ifdef __NUEONLY__
    void  AppendUncorrelatedXSecSystematicsApp       (FitParamList * fp, int start);
    void  AppendCorrelatedBANFFSystematicsApp        (FitParamList * fp, int start);
    void  AppendSuperkDetectorFSISystematicsApp      (FitParamList * fp, int start);
#endif
  
  }
} // valor
  
#endif


