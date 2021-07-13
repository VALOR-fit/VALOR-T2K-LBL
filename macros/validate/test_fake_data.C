//
//______________________________________________________________________________________________________________________
//
// Francis Bench
// 2020-11-25
//______________________________________________________________________________________________________________________
//

using namespace valor;
namespace vutil = valor::util;

//______________________________________________________________________________________________________________________
// Main
void test_fake_data(
  const char * config_file        = "analysis/t2k_analysis_configuration_oa2020_run1-10_v1_postfit.xml",
  const char * fakedata_filename  = "/scratch/fbench/t2k_run10/fc_toys/run1-10_FC_fake_data_normal_cpv_reactor_0_1_40000.bin",
  int          ntoys              = 40000,
  bool         use_simple_binning = 0
)
{
  std::cout << "*** Starting" << std::endl;
  CompoundFitElements* elements = new CompoundFitElements();
  elements->Read(util::CompleteConfigPath(config_file));
  FitModelABC* model = elements->FitModel();
  FitDataABC* data = elements->FitData();
  std::cout << "*** Set up CFM and CFD" << std::endl;

  const char* read_opt_data = "data_vf";
  if (use_simple_binning) {
    read_opt_data = "data_vf+SB";
  }

  // First generate a dummy fake data toy to be overwritten (so that the data arrays can be initialised)
  FitParamList* fpl_void = util::GetFitParameterList("void", "");
  bool with_stat_fluct = true;
  data->Generate(model, *fpl_void, with_stat_fluct);
  std::cout << "*** Set up dummy fake data" << std::endl;

  // Fit to data read from an input binary toy file
  util::ToyStream* toy_stream = data->InitToyStream(fakedata_filename, "read", read_opt_data);
  std::cout << "*** Initialised ToyStream" << std::endl;

  // Test reading each toy
  for(int i = 0; i < ntoys; ++i) {
    data->Read(toy_stream);
    if (i % 100 == 0) std::cout << "*** Read " << i << "th toy" << std::endl;
  }
  std::cout << "*** Successfully read all toys in sequence" << std::endl;

  // Test jumping around the file to specific fake data toys
  toy_stream->SetPositionToToy(0);
  std::cout << "*** Set position to the 1st toy" << std::endl;
  data->Read(toy_stream);
  std::cout << "*** Successfully read the 1st toy" << std::endl;
  toy_stream->SetPositionToToy(ntoys-1);
  std::cout << "*** Set position to the last toy" << std::endl;
  data->Read(toy_stream);
  std::cout << "*** Successfully read the last toy" << std::endl;
  toy_stream->SetPositionToToy((int)(ntoys/2));
  std::cout << "*** Set position to the halfway toy" << std::endl;
  data->Read(toy_stream);
  std::cout << "*** Successfully read the halfway toy" << std::endl;

  std::cout << "*** Done!" << std::endl;
}

//______________________________________________________________________________________________________________________