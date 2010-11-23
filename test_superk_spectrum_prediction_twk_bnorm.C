//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of SuperK
// cross section tweaking dials.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_bnorm()
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

  const int ntwkdials = 5;
  double twkdial[ntwkdials] = { -1.0, -0.5, 0.0, +0.5, +1.0 };

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2,       1,       1      };
  int hcolor [ntwkdials] = { kBlue,   kBlue,   kBlack,  kRed,    kRed   };
  int hstyle [ntwkdials] = { kDotted, kDashed, kSolid,  kDashed, kDotted};

  // histogram booking
  TH1D * hEreco[ntwkdials];
  for(int idial = 0; idial < ntwkdials; idial++) {
     hEreco[idial] = new TH1D("","",nbins,Emin,Emax);
     hEreco[idial]->SetDirectory(0);
     hEreco[idial]->SetLineWidth(hwidth[idial]);
     hEreco[idial]->SetLineColor(hcolor[idial]);
     hEreco[idial]->SetLineStyle(hstyle[idial]);
  }

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, 0);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        0);

  // generate predictions
  for(int idial = 0; idial < ntwkdials; idial++) {
     skexp.TweakBeamNorm( twkdial[idial] );
     for(int ibin=1; ibin<=nbins; ibin++) {
        double Ereco = hEreco[idial]->GetBinCenter(ibin);
        double Nexp  = skexp.Generate(Ereco, *fparm);
        cout << "NumExpSK (Ereco = " << Ereco << ") = " << Nexp << endl;
        hEreco[idial]->SetBinContent(ibin,Nexp);
     }//ibin
  }//idial

  // draw
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  for(int idial = 0; idial < ntwkdials; idial++) {
    const char * opt   = (idial==0) ? "" : "same";
    const char * label = Form("Beam normalization tweaking dial = %.2f", twkdial[idial]);
    if(idial==0) {
      hEreco[idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hEreco[idial]->GetYaxis()->SetTitle("Number of events");
    }
    hEreco[idial]->Draw(opt);
    legend->AddEntry(hEreco[idial], label, "L");
    cout << " *** " << label << " : NumEvents = " << hEreco[idial]->Integral() << endl;
  }//idial
  legend->Draw();
  c->Update();
}
