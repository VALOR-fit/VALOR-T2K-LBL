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

  // energy range and binning
  const double Emin  = 0;
  const double Emax  = 4;
  const double dE    = 0.05;
  const int    nbins = (int) ((Emax - Emin)/dE);

  // oscillation params
  const int ncases = 4;
  double sinsq2theta23[ncases] = { 0., 1.0,    1.0,     1.0    };
  double dm23sq       [ncases] = { 0., 2.0E-3, 2.35E-3, 3.3E-3 };

  // histogram drawing options
  int hwidth [ncases] = { 2,      1,       1,       1       };
  int hcolor [ncases] = { kBlack, kBlue,   kGreen,  kRed    };
  int hstyle [ncases] = { kBlack, kDashed, kDashed, kDashed };

  // histogram booking
  TH1D * hEreco[ncases];
  for(int icase = 0; icase < ncases; icase++) {
    hEreco[icase] = new TH1D("","",nbins,Emin,Emax);
    hEreco[icase]->SetDirectory(0);
    hEreco[icase]->SetLineWidth(hwidth[icase]);
    hEreco[icase]->SetLineColor(hcolor[icase]);
    hEreco[icase]->SetLineStyle(hstyle[icase]);
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
    for(int ibin=1; ibin<=nbins; ibin++) {
       double Ereco = hEreco[icase]->GetBinCenter(ibin);
       double Nexp  = skexp.Generate(Ereco, *fparm);
       cout << "NumExpSK (Ereco = " << Ereco << ") = " << Nexp << endl;
       hEreco[icase]->SetBinContent(ibin,Nexp);
    }
  }

  // draw
  gStyle->SetOptStat(0);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  for(int icase = 0; icase < ncases; icase++) {
     const char * opt   = (icase==0) ? ""             : "same";
     const char * label = (icase==0) ? "Unoscillated" : Form("#Deltam^{2} = %.5f eV^{2}/c^{4}, sin^{2}2#theta = %.2f", dm23sq[icase], sinsq2theta23[icase]);
     if(icase==0) {
       hEreco[icase]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       hEreco[icase]->GetYaxis()->SetTitle("Number of events");
     }
     hEreco[icase]->Draw(opt);
     legend->AddEntry(hEreco[icase], label, "L");

     cout << " *** " << label << " : NumEvents = " << hEreco[icase]->Integral() << endl;
  }
  legend->Draw();
  
}
