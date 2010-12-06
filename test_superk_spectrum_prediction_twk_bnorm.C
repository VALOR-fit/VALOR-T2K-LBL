//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various flux normalization
// tweaking dials.
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

  const int ntwkdials = 3;
  double twkdial[ntwkdials] = { -1.0, 0.0, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2      };
  int hcolor [ntwkdials] = { kBlue,   kBlack,  kRed   };
  int hstyle [ntwkdials] = { kDashed, kSolid, kDashed };

  // histogram booking
  TH1D * hEreco   [ntwkdials];
  TH1D * hErecoErr[ntwkdials];
  for(int idial = 0; idial < ntwkdials; idial++) {
     hEreco   [idial] = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
     hEreco   [idial]->SetDirectory(0);
     hEreco   [idial]->SetLineWidth(hwidth[idial]);
     hEreco   [idial]->SetLineColor(hcolor[idial]);
     hEreco   [idial]->SetLineStyle(hstyle[idial]);
     hErecoErr[idial] = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
     hErecoErr[idial]->SetDirectory(0);
     hErecoErr[idial]->SetLineWidth(hwidth[idial]);
     hErecoErr[idial]->SetLineColor(hcolor[idial]);
     hErecoErr[idial]->SetLineStyle(hstyle[idial]);
  }

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, 1.);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        2.4E-3);

  // generate predictions
  for(int idial = 0; idial < ntwkdials; idial++) {
     cout << " ** Beam norm tweaking dial = " << twkdial[idial] << endl;
     skexp.GetBeamInputs()->TweakNorm( twkdial[idial] );
     for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
        cout << "ibin = " << ibin << endl;
        double Ereco = hEreco[idial]->GetBinCenter(ibin);
        cout << "Ereco = " << Ereco << endl;
        double Nexp  = skexp.Generate(Ereco, *fparm);
        cout << "NumExpSK (Ereco = " << Ereco << ") = " << Nexp << endl;
        hEreco[idial]->SetBinContent(ibin,Nexp);
     }//ibin
  }//idial

  // get fractional err
  double maxerrpos = -1.;
  double maxerrneg = +1.;
  for(int idial = 0; idial < ntwkdials; idial++) {
       for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
         double ncurr  = hEreco[idial]->GetBinContent(ibin);
         double ndeflt = hEreco[deflt]->GetBinContent(ibin);
         if(ndeflt>0) {
           double err = (ncurr-ndeflt)/ndeflt;
           maxerrpos = TMath::Max(maxerrpos, err);
           maxerrneg = TMath::Min(maxerrneg, err);
           hErecoErr[idial]->SetBinContent(ibin, (ncurr-ndeflt)/ndeflt);
         }
       }//ibin
  }//idial

  // draw
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  c->Divide(1,2);
  c->cd(1);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetHeader("Beam normalization");
  for(int idial = 0; idial < ntwkdials; idial++) {
    const char * opt   = (idial==0) ? "" : "same";
    const char * label = Form("Tweaking dial = %.2f", twkdial[idial]);
    if(idial==0) {
      hEreco[idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hEreco[idial]->GetYaxis()->SetTitle("Number of events");
    }
    hEreco[idial]->Draw(opt);
    legend->AddEntry(hEreco[idial], label, "L");
    cout << " *** " << label << " : NumEvents = " << hEreco[idial]->Integral() << endl;
  }//idial
  legend->Draw();
  c->cd(2);
  TH1F * hframe = (TH1F*) c->GetPad(2)->DrawFrame(
       0., 1.2*maxerrneg, hErecoErr[deflt]->GetXaxis()->GetXmax(), 1.2*maxerrpos);
    hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hframe->GetYaxis()->SetTitle("Fractional variation");
    hframe->Draw();
    for(int idial = 0; idial < ntwkdials; idial++) {
      if(idial==deflt) continue;
      hErecoErr[idial]->Draw("same");
  }
  c->Update();
}
