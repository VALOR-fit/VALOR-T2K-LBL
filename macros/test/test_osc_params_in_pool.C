
// For testing the Parameter and ParameterPool class

// Francis Bench

void test_osc_params_in_pool
(
  const char * config_file = "analysis/t2k_analysis_configuration_oa2018_run1-9d_v1_postfit.xml"
)
{
  using std::cout;
  using std::endl;
  using namespace valor;
  using namespace valor::conv;

  // Analysis configuration
  TString config_file_full = Form("%s/config/%s", gSystem->Getenv("VALOR_EXPT"), config_file);
  CompoundFitElements* elements = new CompoundFitElements();
  elements->Read(config_file_full.Data(), "nomc nosplines");

  ParameterPool::Instance()->CheckForCollisions();
  std::cout << "\n\nAfter reading the experiment-specific syst and osc params:" << std::endl;
  ParameterPool::Instance()->Print(std::cout);

  FitParamList * fp = util::GetFitParameterList("3f_osc", "");
  fp->Print();

}
