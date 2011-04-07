//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of 
// cross section tweaking dials.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_xsec
(
  double       input_sinsq2theta23 = 1.0,
  double       input_dm23sq        = 2.4E-3,
  const char * config_file         = "analysis_configuration.xml"
)
{
  using namespace std;
  using namespace t2kana;
  
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  const int ntwkmodes = 5;
  string twklabel[ntwkmodes] = {
   	"CCQE cross section",
	"CC1#pi/CCQE cross section ratio",
        "CCother/CCQE cross section ratio",
        "NC/CCQE cross section ratio",
        "Nuclear effect"
  };
  const int ntwkdials = 5;
  double twkdial[ntwkdials] = { -1.0, -0.5, 0.0, +0.5, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2,       1,       1      };
  int hcolor [ntwkdials] = { kBlue,   kBlue,   kBlack,  kRed,    kRed   };
  int hstyle [ntwkdials] = { kDotted, kDashed, kSolid,  kDashed, kDotted};

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

  // calculate spectra
  TH1D * hEreco   [ntwkmodes][ntwkdials];
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
       // reset previously set tweaks
       skexp.ResetSystematicVariations(); 
       // tweak spectrum generation model
       if      (imode==0) skexp.GetNuXSec()->TweakCCQEXSec            ( twkdial[idial] );
       else if (imode==1) skexp.GetNuXSec()->TweakCC1piToCCQEXSecRatio( twkdial[idial] );
       else if (imode==2) skexp.GetNuXSec()->TweakCCOthToCCQEXSecRatio( twkdial[idial] );
       else if (imode==3) skexp.GetNuXSec()->TweakNCToCCQEXSecRatio   ( twkdial[idial] );
       else if (imode==4) skexp.GetNuXSec()->TweakNuclEffect          ( twkdial[idial] );
       // generate predictions  
       skexp.Generate(*fparm);
       // store spectrum
       hEreco[imode][idial] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
       hEreco[imode][idial]->SetDirectory(0);
    }//idial
  }//imode
  // get fractional err
  TH1D * hErecoErr[ntwkmodes][ntwkdials];
  double maxerrpos = -1.;
  double maxerrneg = +1.;
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
       hErecoErr[imode][idial] = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
       hErecoErr[imode][idial]->SetDirectory(0);
       for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
         double ncurr  = hEreco[imode][idial]->GetBinContent(ibin);
         double ndeflt = hEreco[imode][deflt]->GetBinContent(ibin);
         if(ndeflt>0) {
           double err = (ncurr-ndeflt)/ndeflt;
           maxerrpos = TMath::Max(maxerrpos, err);
           maxerrneg = TMath::Min(maxerrneg, err);          
           hErecoErr[imode][idial]->SetBinContent(ibin, (ncurr-ndeflt)/ndeflt);
         }
       }//ibin
    }//idial
  }//imode

  //
  // draw
  //
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
    c[imode] = new TCanvas(Form("c%d",imode),"",20,20,500,500);
    c[imode]->Divide(1,2);
    c[imode]->cd(1);
    TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetHeader(twklabel[imode].c_str());
    cout << " *** " << twklabel[imode] << endl;
    for(int idial = 0; idial < ntwkdials; idial++) {
       const char * opt   = (idial==0) ? "" : "same";
       const char * label = Form("tweak dial = %.2f", twkdial[idial]);
       if(idial==0) {
        hEreco[imode][idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
        hEreco[imode][idial]->GetYaxis()->SetTitle("Number of events");
       }
       hEreco[imode][idial]->Draw(opt);
       legend->AddEntry(hEreco[imode][idial], label, "L");
       cout << " * " << label << " : NumEvents = " << hEreco[imode][idial]->Integral() << endl;
    }//idial
    legend->Draw();
    c[imode]->cd(2);
    TH1F * hframe = (TH1F*) c[imode]->GetPad(2)->DrawFrame(
       0., 1.2*maxerrneg, hErecoErr[imode][deflt]->GetXaxis()->GetXmax(), 1.2*maxerrpos);
    hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hframe->GetYaxis()->SetTitle("Fractional variation");
    hframe->Draw();
    for(int idial = 0; idial < ntwkdials; idial++) {
      if(idial==deflt) continue;
      hErecoErr[imode][idial]->Draw("same");
    }
    c[imode]->Update();
  }//imode

  // write-out in file
  TFile f(Form("sk_ereco_spectrum-tweak_cross_sections-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",input_sinsq2theta23,input_dm23sq),"recreate");
  for(int imode = 0; imode < ntwkmodes; imode++) {
     for(int idial = 0; idial < ntwkdials; idial++) {
       const char * syst = "";
       if(imode==0) syst="ccqe_xsec";
       if(imode==1) syst="cc1pi_to_ccqe_xsec_ratio";
       if(imode==2) syst="ccoth_to_ccqe_xsec_ratio";
       if(imode==3) syst="nc_to_ccqe_xsec_ratio";
       if(imode==4) syst="nuclear_effect";
       const char * dial = "";
       if( TMath::Abs(twkdial[idial]-1.0) < conv::kASmallNum ) dial="plus_one_sigma";
       if( TMath::Abs(twkdial[idial]-0.5) < conv::kASmallNum ) dial="plus_half_sigma";
       if( TMath::Abs(twkdial[idial]    ) < conv::kASmallNum ) dial="nominal";
       if( TMath::Abs(twkdial[idial]+1.0) < conv::kASmallNum ) dial="minus_one_sigma";
       if( TMath::Abs(twkdial[idial]+0.5) < conv::kASmallNum ) dial="minus_half_sigma";
       const char * title = Form("sk_ereco_spectrum_twk_%s_%s",syst,dial);
       hEreco[imode][idial]->Write(title);
     }
  }
  f.Close();

}
