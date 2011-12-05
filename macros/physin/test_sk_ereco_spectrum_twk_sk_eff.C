//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of SuperK
// nuisance parameter variations that affect the SuperK efficiency.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_sk_ereco_spectrum_twk_sk_eff
(
  const char * hypothesis         = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12     = 0.86,                           /* */
  double       input_mixing13     = 0.10,                           /* */
  double       input_mixing23     = 1.00,                           /* */
  double       input_deltacp      = 0.00,                           /* */
  double       input_dm12sq       = 7.58E-5,                        /* */
  double       input_dm23sq       = 2.32E-3,                        /* */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;
	
  const int nexperiments = 100;

  const int ntwkdials = 5;
  double twkdial[ntwkdials] = { -1.0, -0.5, 0.0, +0.5, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  const int ntwkmodes = 8;

  string systematic[ntwkmodes] = { 
     "sk_ccqe_other_eff", 
     "sk_ccqe_ring_counting_eff_0", /* change RC efficiency < 0.4 GeV,        ignore correlations with other bins (for validation purposes only) */
     "sk_ccqe_ring_counting_eff_1", /* change RC efficiency in 0.4 - 1.1 GeV, ignore correlations with other bins (for validation purposes only) */
     "sk_ccqe_ring_counting_eff_2", /* change RC efficiency > 1.1 GeV,        ignore correlations with other bins (for validation purposes only) */
     "sk_ccqe_ring_counting_eff",   /* change all RC effiencies + include correlations */
     "sk_ccnonqe_eff", 
     "sk_nc_eff", 
     "sk_nuecc_eff" 
  };
  bool generate_toy_expts[ntwkmodes] = { 
     false,    /* sk_ccqe_other_eff              */
     false,    /* sk_ccqe_ring_counting_eff_0    */
     false,    /* sk_ccqe_ring_counting_eff_1    */ 
     false,    /* sk_ccqe_ring_counting_eff_2    */
     true,     /* sk_ccqe_ring_counting_eff      */
     false,    /* sk_ccnonqe_eff                 */
     false,    /* sk_nc_eff                      */
     false     /* sk_nuecc_eff                   */
  };
  string legend_label[ntwkmodes] = { 
    "CCQE (other)",
    "CCQE (RC0 only)",
    "CCQE (RC1 only)",
    "CCQE (RC2 only)",
    "CCQE (RC)",
    "CCnonQE",
    "NC",
    "#nu_{e} CC"
  };

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2,       1,       1       };
  int hcolor [ntwkdials] = { kBlue,   kBlue,   kBlack,  kRed,    kRed    };
  int hstyle [ntwkdials] = { kDotted, kDashed, kSolid,  kDashed, kDotted };

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }
  cout << "Physics hypothesis set: " << endl;
  cout << ph;
  
  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation skexp;
  skexp.Read(config);


  cout << "\n\nCalculating effects of SuperK efficiency uncertainties\n" << endl;

  // calc spectra
  double maxnev[ntwkmodes] = {-1.};
  TH1D * hEreco [ntwkmodes][ntwkdials];
  for(int imode = 0; imode < ntwkmodes; imode++) {
     cout << "*Calculating effect of " << systematic[imode] << " systematic..." << endl;
     if(!generate_toy_expts[imode]) 
     {
        for(int idial = 0; idial < ntwkdials; idial++) {
           // reset previously set tweaks
           skexp.ResetSystematicVariations();
           // tweak spectrum generation model
           skexp.TweakSystematic(systematic[imode], twkdial[idial]);
           // generate predictions  
           skexp.Generate(ph);
           // store spectrum
           hEreco[imode][idial] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
           hEreco[imode][idial]->SetDirectory(0);
        }//idial
     } 
     else {
        cout << "Generating toy MC... " << endl;
        TH1D * hErecoToyMC[nexperiments]; 
        for(int iexp = 0; iexp < nexperiments; iexp++) {
           if(iexp%100==0) {
               cout << "Generated " << iexp << " toy MC experiments so far... " << endl;
           }
           // reset previously set tweaks
           skexp.ResetSystematicVariations();
           // tweak spectrum generation model
           if(iexp>0) {
               skexp.TweakSystematic(systematic[imode]);
           }
           // generate predictions  
           skexp.Generate(ph);
           // store spectrum
           hErecoToyMC[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
           hErecoToyMC[iexp]->SetDirectory(0);
        }//iexp
        // Initialise Ereco distributions for all twk dials to nominal distribution
        for(int idial = 0; idial < ntwkdials; idial++) {
           hEreco[imode][idial] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
        } // idial
        // Caluclate 1-sigma for each bin 
        for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
           double nom   = hErecoToyMC[0]->GetBinContent(ibin);
           double sumsq = 0;
           for(int iexp = 1; iexp < nexperiments; iexp++) {
              double N = hErecoToyMC[iexp]->GetBinContent(ibin);
              sumsq += TMath::Power(N-nom,2);	
           }//iexp
           double sig = TMath::Sqrt(sumsq/(nexperiments-1));   
           for(int idial = 0; idial < ntwkdials; idial++) {
              double nevtwk = nom + twkdial[idial] * sig;
              hEreco[imode][idial]->SetBinContent(ibin, nevtwk);
           } //idial          
        }//ibin
        //delete [] hErecoToyMC;
     }//generate toy MC?
  }//imode
  // get fractional err
  double maxerrpos[ntwkmodes] = {-1.};
  double maxerrneg[ntwkmodes] = {+1.};
  TH1D * hErecoErr[ntwkmodes][ntwkdials];
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
      hErecoErr[imode][idial] = new TH1D("","",conv::kNErecoBins1Rmu,conv::kErecoBinEdges1Rmu);
      hErecoErr[imode][idial]->SetDirectory(0);
      for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
         double ncurr  = hEreco[imode][idial]->GetBinContent(ibin);
         double ndeflt = hEreco[imode][deflt]->GetBinContent(ibin);
         if(ndeflt>0) {
            double err = (ncurr-ndeflt)/ndeflt;
            maxerrpos[imode] = TMath::Max(maxerrpos[imode], err);
            maxerrneg[imode] = TMath::Min(maxerrneg[imode], err);
            maxnev[imode]    = TMath::Max(maxnev[imode], ncurr);
            hErecoErr[imode][idial]->SetBinContent(ibin, (ncurr-ndeflt)/ndeflt);
         }
      }//ibin
    }//idial
  }//imode

  // draw
  for(int imode = 0; imode < ntwkmodes; imode++) {
     for(int idial = 0; idial < ntwkdials; idial++) {
       hEreco   [imode][idial]->SetLineWidth(hwidth[idial]);
       hEreco   [imode][idial]->SetLineColor(hcolor[idial]);
       hEreco   [imode][idial]->SetLineStyle(hstyle[idial]);
       hErecoErr[imode][idial]->SetLineWidth(hwidth[idial]);
       hErecoErr[imode][idial]->SetLineColor(hcolor[idial]);
       hErecoErr[imode][idial]->SetLineStyle(hstyle[idial]);
     }
  }
  TCanvas * c[ntwkmodes];
  gStyle->SetOptStat(0);
  for(int imode = 0; imode < ntwkmodes; imode++) {
    c[imode] = new TCanvas(Form("c%d",imode),"",20,20,900,500);
    c[imode]->Divide(1,2);
    c[imode]->cd(1);
    TLegend * legend1 = new TLegend(0.3, 0.7, 0.85, 0.9);
    legend1->SetFillColor(0);
    legend1->SetBorderSize(0);
    legend1->SetHeader(Form("SuperK %s efficiency", legend_label[imode].c_str()));
    cout << "*** " << Form("SuperK %s efficiency", legend_label[imode].c_str()) << endl;
    TH1F * hframe1 = (TH1F*) c[imode]->GetPad(1)->DrawFrame(
       0., 0., hErecoErr[imode][deflt]->GetXaxis()->GetXmax(), 1.15*maxnev[imode]);
    hframe1->Draw();
    hframe1->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hframe1->GetYaxis()->SetTitle("Number of events");
    for(int idial = 0; idial < ntwkdials; idial++) {
      const char * label = Form("tweak dial = %.2f", twkdial[idial]);
      hEreco[imode][idial]->Draw("same");
      legend1->AddEntry(hEreco[imode][idial], label, "L");
      cout << " * " << label << " : NumEvents = " << hEreco[imode][idial]->Integral() << endl;
    }//idial
    legend1->Draw();
    c[imode]->cd(2);
    TH1F * hframe2 = (TH1F*) c[imode]->GetPad(2)->DrawFrame(
       0., 1.2*maxerrneg[imode], hErecoErr[imode][deflt]->GetXaxis()->GetXmax(), 1.2*maxerrpos[imode]);
    hframe2->GetXaxis()->SetTitle("E_{reco} (GeV)");                                
    hframe2->GetYaxis()->SetTitle("Fractional variation");
    hframe2->Draw();
    for(int idial = 0; idial < ntwkdials; idial++) {
       if(idial==deflt) continue;
       hErecoErr[imode][idial]->Draw("same");
     }
     c[imode]->Update();
  }//imode

  // write-out in file
  string oscparams = util::AsString(hypothesis,
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_sk_eff-%s.root",oscparams.c_str()),"recreate");
  for(int imode = 0; imode < ntwkmodes; imode++) {
     for(int idial = 0; idial < ntwkdials; idial++) {
       const char * dial = "";
       if( TMath::Abs(twkdial[idial]-1.0) < conv::kASmallNum ) dial="plus_one_sigma";
       if( TMath::Abs(twkdial[idial]-0.5) < conv::kASmallNum ) dial="plus_half_sigma";
       if( TMath::Abs(twkdial[idial]    ) < conv::kASmallNum ) dial="nominal";
       if( TMath::Abs(twkdial[idial]+1.0) < conv::kASmallNum ) dial="minus_one_sigma";
       if( TMath::Abs(twkdial[idial]+0.5) < conv::kASmallNum ) dial="minus_half_sigma";
       const char * title = Form("sk_ereco_spectrum_twk_%s_%s",systematic[imode].c_str(),dial);
       hEreco[imode][idial]->Write(title);
     }
  }
  f.Close();
}
