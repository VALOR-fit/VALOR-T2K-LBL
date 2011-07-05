//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction
(
  const char * config_file = "analysis_configuration.xml"
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // oscillation params
  const int ncases = 2;
  double input_sinsq2theta23[ncases] = { 0.00, 1.00    };
  double input_dm23sq[ncases]        = { 0.00, 2.32E-3 };

  double line_width[ncases] = {      2,       2 };
  double line_style[ncases] = { kSolid, kDashed };

  // histogram booking
  TH1D * hErecoTot[ncases];

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // physics hypothesis
  PhysHypothesis ph(t2kana::conv::kHypNuMuDisapp2F);

  for(int icase = 0; icase < ncases; icase++) {
    ph.SetCurrValue(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
    ph.SetCurrValue(t2kana::conv::kPhypDm23sq,        input_dm23sq);
    skexp.Generate(ph);
    hErecoTot[icase] = new TH1D(*skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike));
    hErecoTot[icase]->SetDirectory(0);
    hErecoTot[icase]->SetLineWidth(line_width[icase]);
    hErecoTot[icase]->SetLineStyle(line_style[icase]);
  }//icase

  // draw
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TLegend * legend = new TLegend(0.5,0.1,0.9,0.5);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot[0]->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot[0]->GetYaxis()->SetTitle("number of events in bin");
  for(int icase = 0; icase < ncases; icase++) {
     const char * opt   = (icase == 0) ? "" : "same";
     const char * label = (icase == 0) ? "No oscillation" : 
       Form("#Delta,^{2}_{23}=%.5f eV^{2}/c^{4}, sin^{2}2#theta_{23}=%.1f",input_dm23sq[icase],input_sinsq2theta23[icase]);
     hErecoTot[icase]->Draw(opt);
     legend->AddEntry(hErecoTot[icase], label, "L");
     double E1 = 0.75; // GeV
     double E2 = 1.00; // GeV
     double E3 = 2.00; // GeV
     double ntot = hErecoTot[icase]->Integral();
     cout << "*** " << label << endl;
     cout << "  Total number of events: " << ntot << endl;
     cout << "  Total number of events (< " << E1 << " GeV): " << hErecoTot[icase]->Integral(1,hErecoTot[icase]->FindBin(E1)) << endl;
     cout << "  Total number of events (< " << E2 << " GeV): " << hErecoTot[icase]->Integral(1,hErecoTot[icase]->FindBin(E2)) << endl;
     cout << "  Total number of events (<  "<< E3 << " GeV): " << hErecoTot[icase]->Integral(1,hErecoTot[icase]->FindBin(E3)) << endl;
  }
  legend->Draw();
  c->Update();
}
