//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_contributions()
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
  double input_sinsq2theta23 = 1.00;
  double input_dm23sq        = 2.4E-3;

  const unsigned int nm = 6;
  conv::SuperkMode_t modes [nm] = {
     conv::kSKmNuMuCCQE, 
     conv::kSKmNuMuCC1pi,
     conv::kSKmNuMuCCoth, 
     conv::kSKmNC, 
     conv::kSKmNuMuBar, 
     conv::kSKmNuE
  };

  // histogram booking
  TH1D * hErecoTot = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
  hErecoTot->SetDirectory(0);
  TH1D * hErecoMode[nm]; //< contribution from mode
  for(int im = 0; im < nm; im++) {
    hErecoMode[im] = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
    hErecoMode[im]->SetDirectory(0);
  }

  hErecoTot                       -> SetLineWidth ( 2        );
  hErecoMode[conv::kSKmNuMuCC1pi] -> SetFillStyle ( 3354     );
  hErecoMode[conv::kSKmNuMuCCoth] -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmNC]        -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuMuBar]   -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmNuE]       -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmNuMuCC1pi] -> SetFillColor ( kBlue    );
  hErecoMode[conv::kSKmNuMuCCoth] -> SetFillColor ( kMagenta );
  hErecoMode[conv::kSKmNC]        -> SetFillColor ( kRed     );
  hErecoMode[conv::kSKmNuMuBar]   -> SetFillColor ( kYellow  );
  hErecoMode[conv::kSKmNuE]       -> SetFillColor ( kGreen   );
  hErecoMode[conv::kSKmNuMuCCQE]  -> SetLineColor ( kBlack   );
  hErecoMode[conv::kSKmNuMuCC1pi] -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCoth] -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNC]        -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCQE]  -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBar]   -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuE]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNC]        -> SetLineStyle ( kDashed  );
  

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  // generate predictions
  for(int im = 0; im < nm; im++) {
    for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
       double Ereco = hErecoMode[im]->GetBinCenter(ibin);
       double Nexp  = skexp.Generate(Ereco, *fparm, conv::kSKsFCFV1RingMuLike, modes[im]);
       cout << "NumExpSK (Ereco = " << Ereco << " / " << util::AsString(modes[im]) << ") = " << Nexp << endl;
       hErecoMode[im]->SetBinContent(ibin,Nexp);
    }
  }
  for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
    double Ereco = hErecoTot->GetBinCenter(ibin);
    double Nexp  = skexp.Generate(Ereco, *fparm);
    cout << "NumExpSK (Ereco = " << Ereco << " / total) = " << Nexp << endl;
    hErecoTot->SetBinContent(ibin,Nexp);
  }

  // draw
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot->Draw();
  hErecoTot->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot->GetYaxis()->SetTitle("Number of events");
  THStack * stack = new THStack("stack","spectrum contributions");
  stack->Add(hErecoMode[conv::kSKmNuMuCCQE]);
  stack->Add(hErecoMode[conv::kSKmNuMuCC1pi]);
  stack->Add(hErecoMode[conv::kSKmNuMuCCoth]);
  stack->Add(hErecoMode[conv::kSKmNC]);
  stack->Add(hErecoMode[conv::kSKmNuMuBar]);
  stack->Add(hErecoMode[conv::kSKmNuE]);
  stack->Draw("same");
  legend->AddEntry(hErecoTot,                       "total",                   "L" );
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCQE],  "#nu_{#mu} CCQE",          "L" );
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCC1pi], "#nu_{#mu} CC1#pi",        "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCoth], "#nu_{#mu} CC other",      "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNC],        "#nu_{#mu}/#nu_{#tau} NC", "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBar],   "#bar{#nu_{#mu}}",         "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuE],       "#nu_{e}",                 "LF");
  legend->Draw();
  c->Update();
  double E1 = 0.75; // GeV
  double E2 = 1.00; // GeV
  double E3 = 2.00; // GeV
  double ntot = hErecoTot->Integral();
  cout << "Total number of events: " << ntot << endl;
  cout << "Total number of events (< " << E1 << " GeV): " << hErecoTot->Integral(1,hErecoTot->FindBin(E1)) << endl;
  cout << "Total number of events (< " << E2 << " GeV): " << hErecoTot->Integral(1,hErecoTot->FindBin(E2)) << endl;
  cout << "Total number of events (<  "<< E3 << " GeV): " << hErecoTot->Integral(1,hErecoTot->FindBin(E3)) << endl;
  for(int im = 0; im < nm; im++) {
     double nmode = hErecoMode[im]->Integral();
     cout << "Contrib from mode: " << util::AsString(modes[im]) << " : " << nmode 
          << " events (" << 100*nmode/ntot << " %)" << endl;
  }
}
