//
// ------------------------------------------------------------------------------------------------
// Generate / store mock data
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void write_out_mock_data(
  double sinsq2theta23=0.90, double dm23sq=2.4E-3, const char * filename="mock_data_tight.root")
{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  bool fluctuate = false;

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->CreateMockData(*skexp, *fparm, fluctuate);

  TH1D* mock_data = skobs->AsTH1D(conv::kSKsFCFV1RingMuLike);

  TFile f(filename, "recreate");
  mock_data->Write("mock_data");
  f.Close();
}

