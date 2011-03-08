//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of SuperK
// energy scale variations.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_skescale
(
  double input_sinsq2theta23 = 1.0,
  double input_dm23sq        = 2.4E-3
)
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

  const int ntwkdials = 3;
  double twkdial[ntwkdials] = { -1.0, 0.0, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2      };
  int hcolor [ntwkdials] = { kBlue,   kBlack,  kRed   };
  int hstyle [ntwkdials] = { kDashed, kSolid, kDashed };

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::numu10a::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::numu10a::conv::kPhypDm23sq,        input_dm23sq);
  
  // calc spectra
  TH1D * hEreco[ntwkdials];
  for(int idial = 0; idial < ntwkdials; idial++) {
     // tweak spectrum generation model
     skexp.GetSuperkInputs()->TweakEnergyScale( twkdial[idial] );
     // generate predictions
     skexp.Generate(*fparm);
     // store spectrum
     hEreco[idial] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     hEreco[idial]->SetDirectory(0);
  }//idial
  // get fractional err
  double maxerrpos = -1.;
  double maxerrneg = +1.;
  TH1D * hErecoErr[ntwkdials];
  for(int idial = 0; idial < ntwkdials; idial++) {
    hErecoErr[idial] = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
    hErecoErr[idial]->SetDirectory(0);
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

  //
  // draw
  //
  for(int idial = 0; idial < ntwkdials; idial++) {
     hEreco   [idial]->SetLineWidth(hwidth[idial]);
     hEreco   [idial]->SetLineColor(hcolor[idial]);
     hEreco   [idial]->SetLineStyle(hstyle[idial]);
     hErecoErr[idial]->SetLineWidth(hwidth[idial]);
     hErecoErr[idial]->SetLineColor(hcolor[idial]);
     hErecoErr[idial]->SetLineStyle(hstyle[idial]);
  }
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  c->Divide(1,2);
  c->cd(1);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetHeader("SuperK energy scale");
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

  // write-out in file
  TFile f(Form("sk_ereco_spectrum-tweak_sk_energy_scale-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",input_sinsq2theta23,input_dm23sq),"recreate");
  for(int idial = 0; idial < ntwkdials; idial++) {
     const char * dial = "";
     if( TMath::Abs(twkdial[idial]-1.0) < conv::kASmallNum ) dial="plus_one_sigma";
     if( TMath::Abs(twkdial[idial]-0.5) < conv::kASmallNum ) dial="plus_half_sigma";
     if( TMath::Abs(twkdial[idial]    ) < conv::kASmallNum ) dial="nominal";
     if( TMath::Abs(twkdial[idial]+1.0) < conv::kASmallNum ) dial="minus_one_sigma";
     if( TMath::Abs(twkdial[idial]+0.5) < conv::kASmallNum ) dial="minus_half_sigma";
     const char * title = Form("sk_ereco_spectrum_twk_sk_energy_scale_%s",dial);
     hEreco[idial]->Write(title);
  }
  f.Close();

}
