//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction()
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

  // oscillation params
  double input_sinsq2theta23 = 0.00;
  double input_dm23sq        = 0;

  // histogram booking
  TH1D * hErecoTot = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
  hErecoTot->SetDirectory(0);
  hErecoTot->SetLineWidth(2);
  
  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // oscillation params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
    double Ereco = hErecoTot->GetBinCenter(ibin);
    double Nexp  = skexp.Generate(Ereco, *fparm);
    cout << "NumExpSK (Ereco = " << Ereco << " / total) = " << Nexp << endl;
    hErecoTot->SetBinContent(ibin,Nexp);
  }

  // draw
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  hErecoTot->Draw();
  c->Update();
  double E1 = 0.75; // GeV
  double E2 = 1.00; // GeV
  double E3 = 2.00; // GeV
  double ntot = hErecoTot->Integral();
  cout << "Total number of events: " << ntot << endl;
  cout << "Total number of events (< " << E1 << " GeV): " << hErecoTot->Integral(1,hErecoTot->FindBin(E1)) << endl;
  cout << "Total number of events (< " << E2 << " GeV): " << hErecoTot->Integral(1,hErecoTot->FindBin(E2)) << endl;
  cout << "Total number of events (<  "<< E3 << " GeV): " << hErecoTot->Integral(1,hErecoTot->FindBin(E3)) << endl;
}
