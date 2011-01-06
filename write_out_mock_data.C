//
// ------------------------------------------------------------------------------------------------
// Generate / store mock data.
// The data can be stored in a ROOT file. By default Ereco values get printed as a comma-separated
// list so that can be easily copied & pasted in inputs/data.xml
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void write_out_mock_data(
  double sinsq2theta23=1.00, double dm23sq=2.4E-3, const char * root_filename=0)
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

  gSystem->Load("libT2KDisapp2010a");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::numu10a;

  TRandom rand(0);

  bool fluctuate = true;

  // Oscillation params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::numu10a::conv::kPhypSinsq2Theta23, sinsq2theta23);
  fparm->Update(t2kana::numu10a::conv::kPhypDm23sq,        dm23sq);

  // Get superk spectrum PDF & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // Generate mock data
  SuperkObservation * skobs = new SuperkObservation;
  skobs->CreateMockData(*skexp, *fparm, fluctuate);

  TH1D* mock_data = skobs->ErecoSpectrum(conv::kSKsFCFV1RingMuLike);

  // Print-out Ereco values as a comma separated list.
  // Assuming uniform distribution within each bin, randomize values within each bin.
  bool is_first=true;
  for(int ibin = 1; ibin <= mock_data->GetNbinsX(); ibin++) {
    double Nb   = mock_data->GetBinContent (ibin); // number of events in this Ereco bin
  //double Ec   = mock_data->GetXaxis()->GetBinCenter  (ibin); // Ereco bin centre
    double Emin = mock_data->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
    double Emax = mock_data->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
    for(int iev = 0; iev < Nb; iev++) {
      double E = Emin + rand.Uniform() * (Emax-Emin);
      string delimeter = ",";
      if(is_first) {
	delimeter = "";
        is_first = false;
      }
      cout << delimeter << E;
    }//iev
  }//ibin
  cout << endl;

  // If a ROOT filename was specified, save the ROOT histogram
  if(root_filename) {
    TFile f(filename, "recreate");
    mock_data->Write("mock_data");
    f.Close();
  }

}

