
// For testing the Parameter and ParameterPool class

// Francis Bench


#include "valor/Conventions/ParameterPool.h"
#include "valor/Conventions/Parameter.h"

void test_FitParamNew
(
  const char * param_file             = "params/params_test.xml",
  const char * default_osc_param_file = "params/params_test_default_osc_params.xml"
)
{
  using std::cout;
  using std::endl;
  using namespace valor;
  using namespace valor::conv;

  // Analysis configuration
  TString param_file_full             = Form("%s/config/%s", gSystem->Getenv("VALOR_EXPT"), param_file);
  TString default_osc_param_file_full = Form("%s/config/%s", gSystem->Getenv("VALOR_EXPT"), default_osc_param_file);

  ParameterPool::MutableInstance()->ReadXML(default_osc_param_file_full.Data());
  std::cout << "\n\nAfter reading default osc param settings:" << std::endl;
  ParameterPool::Instance()->Print(std::cout);

  ParameterPool::MutableInstance()->ReadXML(param_file_full.Data());
  ParameterPool::Instance()->CheckForCollisions();
  std::cout << "\n\nAfter reading the experiment-specific syst and osc params:" << std::endl;
  ParameterPool::Instance()->Print(std::cout);

}
