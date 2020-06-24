//
// Read data from a config file and write the data PDFs in both VALORArray and ROOT format
// to a file ready for plotting etc.
//
// Francis Bench <F.Bench \at liverpool.ac.uk>
// December 2019


#include <vector>
#include <utility>

//___________________________________________________________________________
void write_data_pdfs
(
 const char* output_tag     = "SK_data_Run1-9d_valor_PDFs",
 const char* config_file    = ""
 )
{
  CompoundFitElements* elements = new CompoundFitElements();
  elements->Read(util::CompleteConfigPath(config_file));
  FitDataABC* data = elements->FitData();

  TFile* fout = new TFile(output_tag, "RECREATE");

  TObjArray * array = new TObjArray();
  data->Write(array, "valor", "all+SB");

  std::cout << "Writing objects:";
  array->Print();

  array->Write();

  fout->Close();
  delete array;
  delete fout;

  std::cout << "Done!" << std::endl;
}
//___________________________________________________________________________
