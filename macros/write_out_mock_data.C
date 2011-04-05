//
// ------------------------------------------------------------------------------------------------
// Generate / store mock data.
// The data can be stored in an XML (formatted like inputs/data.xml) or a ROOT file. 
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void write_out_mock_data(
  double sinsq2theta23=1.00, 
  double dm23sq=2.4E-3, 
  bool include_systematics=false, 
  double nexperiments=1,
  const char * xml_filename=0,
  const char * root_filename=0)
{
  using namespace std;
  using namespace t2kana;

  TRandom rand(0);

  bool fluctuate = true;

  TFile * f = 0;
  if(root_filename) {
    f = new TFile(root_filename, "recreate");
  }

  ofstream fxml;
  if(xml_filename) {
    fxml.open(xml_filename);
    fxml << "<?xml version=\"1.0\"?> " << endl;
    fxml << "<analysis_inputs>       " << endl;
    fxml << "  <superk>              " << endl;
    fxml << "    <observation>       " << endl;
    fxml << "       <one_mu_ring>    " << endl;
    fxml << "         <ereco_values> " << endl;
  }

  // Oscillation params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::conv::kPhypSinsq2Theta23, sinsq2theta23);
  fparm->Update(t2kana::conv::kPhypDm23sq,        dm23sq);

  // Get superk spectrum PDF & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  for(int iexperiment=0; iexperiment < nexperiments; iexperiment++) {

     cout << " **** toy experiment : " << iexperiment << endl;

     if(include_systematics) {
        // tweak all sources of systematic err (gaussian err)
        skexp->TweakAllSystematics();
     }

     // Generate mock data
     SuperkObservation * skobs = new SuperkObservation;
     skobs->CreateMockData(*skexp, *fparm, fluctuate);

     TH1D* mock_data = skobs->ErecoSpectrum(conv::kSKsFCFV1RingMuLike);
     mock_data->SetDirectory(0);

     // Print-out Ereco values as a comma separated list.
     // Assuming uniform distribution within each bin, randomize values within each bin.
     if(xml_filename) {
       fxml << "           " << Form("<data_set_%d>",iexperiment) << endl;
     }
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
            delimeter = "\t\t";
            is_first = false;
         }
         cout << delimeter << E;
         if(xml_filename) {
            fxml << delimeter << E;           
         }
       }//iev
     }//ibin
     cout << endl;
     if(xml_filename) {
       fxml << endl;
       fxml << "           " << Form("</data_set_%d>",iexperiment) << endl;
     }

     // If a ROOT filename was specified, save the ROOT histogram too
     if(f) {
       f->cd();
       mock_data->Write(Form("mock_data_%d",iexperiment));
     }

  }//iexperiment

  if(f) {
    f->Close();
    delete f;
  }

  if(xml_filename) {
    fxml << "         </ereco_values> " << endl;
    fxml << "       </one_mu_ring>    " << endl;
    fxml << "    </observation>       " << endl;
    fxml << "  </superk>              " << endl;
    fxml << "</analysis_inputs>       ";// << endl;
    fxml.close();
  }

}

