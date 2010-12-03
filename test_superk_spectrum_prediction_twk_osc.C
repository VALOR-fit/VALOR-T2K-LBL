//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various oscillation hypotheses
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_osc()
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
/*
  const int ncases = 5;
  double sinsq2theta23[ncases] = { 0., 1.0,    1.0,    1.0,    1.0  };
  double dm23sq       [ncases] = { 0., 2.0E-3, 2.4E-3, 3.3E-3, 8E-3 };
*/
  const int ncases = 5;
  double sinsq2theta23[ncases] = { 0., 1.0,    0.8,    0.6,    0.4    };
  double dm23sq       [ncases] = { 0., 2.4E-3, 2.4E-3, 2.4E-3, 2.4E-3 };

  // histogram drawing options
  int hwidth [ncases] = { 2,      1,       1,       1,       1        };
  int hcolor [ncases] = { kBlack, kBlue,   kGreen,  kRed,    kBlack   };
  int hstyle [ncases] = { kBlack, kDashed, kDashed, kDashed, kDashed  };

  // histogram booking
  TH1D * hEreco[ncases];
  TH1D * hErecoFracDeficit[ncases];
  for(int icase = 0; icase < ncases; icase++) {
    hEreco[icase] = new TH1D("","",conv::kNErecoBins,conv::kErecoMin,conv::kErecoMax);
    hEreco[icase]->SetDirectory(0);
    hEreco[icase]->SetLineWidth(hwidth[icase]);
    hEreco[icase]->SetLineColor(hcolor[icase]);
    hEreco[icase]->SetLineStyle(hstyle[icase]);
    hErecoFracDeficit[icase] = new TH1D("","",conv::kNErecoBins,conv::kErecoMin,conv::kErecoMax);
    hErecoFracDeficit[icase]->SetDirectory(0);
    hErecoFracDeficit[icase]->SetLineWidth(hwidth[icase]);
    hErecoFracDeficit[icase]->SetLineColor(hcolor[icase]);
    hErecoFracDeficit[icase]->SetLineStyle(hstyle[icase]);
  }

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // generate predictions
  for(int icase = 0; icase < ncases; icase++) {
    fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, sinsq2theta23[icase]);
    fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        dm23sq[icase]);
    for(int ibin=1; ibin <= conv::kNErecoBins; ibin++) {
       double Ereco = hEreco[icase]->GetBinCenter(ibin);
       double Nexp  = skexp.Generate(Ereco, *fparm);
       cout << "NumExpSK (Ereco = " << Ereco << ") = " << Nexp << endl;
       hEreco[icase]->SetBinContent(ibin,Nexp);
    }
  }
  // calculate fractional deficit
  for(int icase = 1; icase < ncases; icase++) {
    for(int ibin=1; ibin <= conv::kNErecoBins; ibin++) {
       double Nunosc = hEreco[0]    ->GetBinContent(ibin); 
       double Nosc   = hEreco[icase]->GetBinContent(ibin);
       if(Nunosc>0) {
          double frac_deficit = (Nunosc-Nosc)/Nunosc;
          hErecoFracDeficit[icase]->SetBinContent(ibin,frac_deficit);
       }
    }
  }

  // draw
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TPad * upp = new TPad("upp","",0.05,0.35,0.95,0.97);
  TPad * low = new TPad("low","",0.05,0.02,0.95,0.30);
  upp->Draw();
  low->Draw();
  upp->cd();
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  for(int icase = 0; icase < ncases; icase++) {
     const char * opt   = (icase==0) ? ""             : "same";
     const char * label = (icase==0) ? "Unoscillated" : 
        Form("#Deltam^{2} = %.5f eV^{2}/c^{4}, sin^{2}2#theta = %.2f", dm23sq[icase], sinsq2theta23[icase]);
     if(icase==0) {
       hEreco[icase]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       hEreco[icase]->GetYaxis()->SetTitle("Number of events");
     }
     hEreco[icase]->Draw(opt);
     legend->AddEntry(hEreco[icase], label, "L");
     double E1 = 0.75; // GeV
     double E2 = 1.00; // GeV
     double E3 = 2.00; // GeV
     cout << " *** " << label << " : " << endl;
     cout << "NumEvents = " << hEreco[icase]->Integral() << endl;
     cout << "NumEvents (< " << E1 << " GeV): " << hEreco[icase]->Integral(1,hEreco[icase]->FindBin(E1)) << endl;
     cout << "NumEvents (< " << E2 << " GeV): " << hEreco[icase]->Integral(1,hEreco[icase]->FindBin(E2)) << endl;
     cout << "NumEvents (<  "<< E3 << " GeV): " << hEreco[icase]->Integral(1,hEreco[icase]->FindBin(E3)) << endl;

  }
  legend->Draw();
  low->cd();
  TH1F * hframe = (TH1F*) low->DrawFrame(0., 0., hEreco[0]->GetXaxis()->GetXmax(), 1);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Fractional deficit");
  hframe->Draw();
  for(int icase = 1; icase < ncases; icase++) {
     hErecoFracDeficit[icase]->Draw("same");
  }
  c->Update();
}
