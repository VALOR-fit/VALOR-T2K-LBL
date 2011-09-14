//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of SuperK
// energy scale variations.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_sk_ereco_spectrum_twk_sk_escale
(
  const char * hypothesis              = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_sinsq2theta12     = 0.86,                           /* */
  double       input_sinsq2theta13     = 0.10,                           /* */
  double       input_sinsq2theta23     = 1.00,                           /* */
  double       input_deltacp           = 0.00,                           /* */
  double       input_dm12sq            = 7.58E-5,                        /* */
  double       input_dm23sq            = 2.32E-3,                        /* */
  const char * config_file             = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  const int ntwkdials = 3;
  double twkdial[ntwkdials] = { -1.0, 0.0, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2      };
  int hcolor [ntwkdials] = { kBlue,   kBlack,  kRed   };
  int hstyle [ntwkdials] = { kDashed, kSolid, kDashed };

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypSinsq2Theta12,  input_sinsq2theta12);
    ph.SetCurrValue(conv::kPhypSinsq2Theta13,  input_sinsq2theta13);
    ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
    ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  }
  cout << "Physics hypothesis set: " << endl;
  cout << ph;
  
  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation skexp;
  skexp.Read(config);
  
  // calc spectra
  TH1D * hEreco[ntwkdials];
  for(int idial = 0; idial < ntwkdials; idial++) {
     // tweak spectrum generation model
     skexp.GetSuperkMC()->TweakEnergyScale( twkdial[idial] );
     // generate predictions
     skexp.Generate(ph);
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
  string oscparams = util::AsString(hypothesis,
     input_sinsq2theta12,input_sinsq2theta13,input_sinsq2theta23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_sk_escale-%s.root",oscparams.c_str()),"recreate");
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
