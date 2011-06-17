//
// ------------------------------------------------------------------------------------------------
// Calculate SuperK Ereco spectrum error envelope. Tweaking all systematics.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_err_envelope
(
  double       input_sinsq2theta23 = 1.0,
  double       input_dm23sq        = 2.32E-3,
  const char * config_file         = "analysis_configuration.xml"
)
{
  gROOT->Reset();

  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  const int nexp = 2500; // number of toy MC experiments

  // output histograms and graphs
  TH1D * hEreco [nexp];        //
  TH1D * hErecoFracErrP  = 0;  //
  TH1D * hErecoFracErrM = 0;   //
  TH1D * hErecoDef = 0;        //
  TH1D * hErecoMax = 0;        //
  TH1D * hErecoMin = 0;        //

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

  // generate nexp random spectra
  for(int iexp = 0; iexp < nexp; iexp++) {
     // tweak spectrum generation model
     if(iexp%50==0) {
        cout << "* Current toy MC experiment: " << iexp << endl;
     }
     if(iexp>0) {
       // use 0^th entry to store the untweaked SK reco energy spectrum
       skexp.TweakAllSystematics();
     }
     // generate predictions
     skexp.Generate(*fparm);
     // store spectrum
     hEreco[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     hEreco[iexp]->SetDirectory(0);
  }

  // store nominal
  hErecoDef = new TH1D(*hEreco[0]);
  hErecoDef->SetDirectory(0);

  // find min/max and calculate fractional deviation from nominal
  hErecoFracErrP = new  TH1D(*hEreco[0]);
  hErecoFracErrM = new  TH1D(*hEreco[0]);
  hErecoMax      = new TH1D(*hEreco[0]);
  hErecoMin      = new TH1D(*hEreco[0]);
  hErecoFracErrP -> SetDirectory(0);
  hErecoFracErrM -> SetDirectory(0);
  hErecoMax      -> SetDirectory(0);
  hErecoMin      -> SetDirectory(0);
  
  for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
    double nom   = hErecoDef->GetBinContent(ibin);
    double sumsq = 0;
    for(int iexp = 1; iexp < nexp; iexp++) {
       double N = hEreco[iexp]->GetBinContent(ibin);
       sumsq += TMath::Power(N-nom,2);	
    }//iexp
    double sig  = TMath::Sqrt(sumsq/(nexp-1));
    double min  = nom-sig;
    double max  = nom+sig; 
    double errp = 0;
    double errm = 0;
    if(nom>0){
      errp = (max-nom)/nom;
      errm = (min-nom)/nom;
    }
    hErecoFracErrP -> SetBinContent (ibin,  errp);
    hErecoFracErrM -> SetBinContent (ibin,  errm);
    hErecoMax      -> SetBinContent (ibin,  max);
    hErecoMin      -> SetBinContent (ibin,  min);
  }//ibin
  hErecoFracErrP->GetYaxis()->SetRangeUser(-1,1);
  hErecoFracErrM->GetYaxis()->SetRangeUser(-1,1);

  // draw
  hErecoDef->SetLineStyle(kDashed);
  hErecoMin->SetLineStyle(kSolid);
  hErecoMax->SetLineStyle(kSolid);
  hErecoMin->SetLineWidth(2);
  hErecoMax->SetLineWidth(2);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->Divide(1,2);
  c->cd(1);
  hErecoMax  -> Draw();
  hErecoMin  -> Draw("same");
  hErecoDef  -> Draw("same");
  hErecoMax->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoMax->GetYaxis()->SetTitle("Number of events");
  c->cd(2);
  hErecoFracErrP->SetLineWidth(2);
  hErecoFracErrM->SetLineWidth(2);
  hErecoFracErrP -> Draw();
  hErecoFracErrM -> Draw("same");
  hErecoFracErrP->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoFracErrP->GetYaxis()->SetTitle("Fractional error");
  c->Update();

  // write-out in file
  TFile f(Form("sk_ereco_spectrum-err_envelope-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",input_sinsq2theta23,input_dm23sq),"recreate");
  hErecoFracErrP -> Write("spectrum_frac_err_plus");
  hErecoFracErrM -> Write("spectrum_frac_err_minus");
  hErecoMax      -> Write("spectrum_error_envelope_upper_boundary");
  hErecoMin      -> Write("spectrum_error_envelope_lower_boundary");
  hErecoDef      -> Write("spectrum_nominal");
  f.Close();
}
