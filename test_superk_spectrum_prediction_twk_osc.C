//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various oscillation hypotheses
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_osc(bool showdata=false)
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

  // oscillation params
  const int ncases = 5;
  double sinsq2theta23[ncases] = { 0., 1.0,    1.0,    1.0,    1.0  };
  double dm23sq       [ncases] = { 0., 2.0E-3, 2.4E-3, 3.3E-3, 8E-3 };
/*
  const int ncases = 5;
  double sinsq2theta23[ncases] = { 0., 1.0,    0.8,    0.6,    0.4    };
  double dm23sq       [ncases] = { 0., 2.4E-3, 2.4E-3, 2.4E-3, 2.4E-3 };
*/

  // histogram drawing options
  int hwidth [ncases] = { 2,      1,       1,       1,       1        };
  int hcolor [ncases] = { kBlack, kBlue,   kGreen,  kRed,    kBlack   };
  int hstyle [ncases] = { kBlack, kDashed, kDashed, kDashed, kDashed  };

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // calculate spectra
  TH1D * hEreco[ncases];
  for(int icase = 0; icase < ncases; icase++) {
    // update osc params
    fparm->Update(t2kana::numu10a::conv::kPhypSinsq2Theta23, sinsq2theta23[icase]);
    fparm->Update(t2kana::numu10a::conv::kPhypDm23sq,        dm23sq[icase]);
    // generate predictions  
    skexp.Generate(*fparm);
    // store spectrum
    hEreco[icase] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
    hEreco[icase]->SetDirectory(0);
  }
  // calculate fractional deficit
  TH1D * hErecoFracDeficit[ncases];
  for(int icase = 0; icase < ncases; icase++) {
    hErecoFracDeficit[icase] = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
    hErecoFracDeficit[icase]->SetDirectory(0);
    for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
       double Nunosc = hEreco[0]    ->GetBinContent(ibin); 
       double Nosc   = hEreco[icase]->GetBinContent(ibin);
       if(Nunosc>0) {
          double frac_deficit = (Nunosc-Nosc)/Nunosc;
          hErecoFracDeficit[icase]->SetBinContent(ibin,frac_deficit);
       }
    }
  }

  //
  // draw
  //
  for(int icase = 0; icase < ncases; icase++) {
    hEreco[icase]->SetLineWidth(hwidth[icase]);
    hEreco[icase]->SetLineColor(hcolor[icase]);
    hEreco[icase]->SetLineStyle(hstyle[icase]);
    hErecoFracDeficit[icase]->SetLineWidth(hwidth[icase]);
    hErecoFracDeficit[icase]->SetLineColor(hcolor[icase]);
    hErecoFracDeficit[icase]->SetLineStyle(hstyle[icase]);
  }
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,700);
  TPad * upp = new TPad("upp","",0.05,0.45,0.95,0.97);
  TPad * low = new TPad("low","",0.05,0.03,0.95,0.40);
  low->SetBottomMargin(.2);
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

     cout << " *** " << label << " : " << endl;

     cout << "NumEvents = " 
          <<  setw(10) << setfill(' ') << setprecision(6) << hEreco[icase]->Integral() << endl;

     const int nemax = 7;
     double emax[nemax] = { 0.50, 0.75, 1.00, 1.50, 2.00, 4.00, 10.00 };

     for (int iemax = 0; iemax < nemax; iemax++) {
          cout << "NumEvents (< " 
               << setw(6)  << setfill(' ') << setprecision(3) 
               << emax[iemax] << " GeV): " 
               << setw(10) << setfill(' ') << setprecision(6)
               << hEreco[icase]->Integral(1,hEreco[icase]->FindBin(emax[iemax])) << endl;
     }
     if(showdata) {
       SuperkObservation * skobs = new SuperkObservation;
       skobs->Read("inputs/analysis_configuration.xml");
       TH1D * data = new TH1D(*skobs->ErecoSpectrum(conv::kSKsFCFV1RingMuLike));
       data->SetMarkerStyle(20);
       data->SetMarkerSize(1);
       data->Draw("PERRSAME");
     }

  }
  legend->Draw();
  low->cd();
  TH1F * hframe = (TH1F*) low->DrawFrame(0., 0., hEreco[0]->GetXaxis()->GetXmax(), 1);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Fractional deficit");
  hframe->GetYaxis()->SetLabelSize(0.08);
  hframe->GetYaxis()->SetTitleSize(0.08);
  hframe->GetYaxis()->SetTitleOffset(0.75);
  hframe->GetXaxis()->SetLabelSize(0.07);
  hframe->GetXaxis()->SetTitleSize(0.08);
  hframe->Draw();
  for(int icase = 1; icase < ncases; icase++) {
     hErecoFracDeficit[icase]->Draw("same");
  }
  c->Update();
}
