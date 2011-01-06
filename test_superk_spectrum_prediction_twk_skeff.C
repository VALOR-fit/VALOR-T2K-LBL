//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of SuperK
// nuisance parameter variations that affect the SuperK efficiency.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_skeff()
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
  double input_sinsq2theta23 = 1.;
  double input_dm23sq        = 2.4E-3;

  const int ntwkdials = 5;
  double twkdial[ntwkdials] = { -1.0, -0.5, 0.0, +0.5, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  const int ntwkmodes = 4;

  string systematic[ntwkmodes] = { "sk_ccqe_eff", "sk_ccnonqe_eff", "sk_nc_eff", "sk_nuecc_eff" };

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2,       1,       1      };
  int hcolor [ntwkdials] = { kBlue,   kBlue,   kBlack,  kRed,    kRed   };
  int hstyle [ntwkdials] = { kDotted, kDashed, kSolid,  kDashed, kDotted};

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::numu10a::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::numu10a::conv::kPhypDm23sq,        input_dm23sq);

  // calc spectra
  double maxnev[ntwkmodes] = {-1.};
  TH1D * hEreco [ntwkmodes][ntwkdials];
  for(int imode = 0; imode < ntwkmodes; imode++) {
     for(int idial = 0; idial < ntwkdials; idial++) {
        // reset previously set tweaks
        skexp.ResetSystematicVariations();
        // tweak spectrum generation model
        skexp.TweakSystematic(systematic[imode], twkdial[idial]);
        // generate predictions  
        skexp.Generate(*fparm);
        // store spectrum
        hEreco[imode][idial] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
        hEreco[imode][idial]->SetDirectory(0);
     }//idial
  }//imode
  // get fractional err
  double maxerrpos[ntwkmodes] = {-1.};
  double maxerrneg[ntwkmodes] = {+1.};
  TH1D * hErecoErr[ntwkmodes][ntwkdials];
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
      hErecoErr[imode][idial] = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
      hErecoErr[imode][idial]->SetDirectory(0);
      for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
         double ncurr  = hEreco[imode][idial]->GetBinContent(ibin);
         double ndeflt = hEreco[imode][deflt]->GetBinContent(ibin);
         if(ndeflt>0) {
            double err = (ncurr-ndeflt)/ndeflt;
            maxerrpos[imode] = TMath::Max(maxerrpos[imode], err);
            maxerrneg[imode] = TMath::Min(maxerrneg[imode], err);
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
    const char * label = "";
    if(imode==0) label = "CCQE";
    if(imode==1) label = "CCnonQE";
    if(imode==2) label = "NC";
    if(imode==3) label = "#nu_{e} CC";
    legend1->SetHeader(Form("SuperK %s efficiency", label));
    cout << "*** " << Form("SuperK %s efficiency", label) << endl;
    TH1F * hframe1 = (TH1F*) c[imode]->GetPad(1)->DrawFrame(
       0., 0., hErecoErr[imode][deflt]->GetXaxis()->GetXmax(), maxnev[imode]);
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
  TFile f(Form("sk_ereco_spectrum-tweak_sk_efficiency-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",input_sinsq2theta23,input_dm23sq),"recreate");
  for(int imode = 0; imode < ntwkmodes; imode++) {
     for(int idial = 0; idial < ntwkdials; idial++) {
       const char * syst = "";
       if(imode==0) syst="ccqe_efficiency";
       if(imode==1) syst="ccnonqe_efficiency";
       if(imode==2) syst="nc_efficiency";
       if(imode==3) syst="nuecc_efficiency";
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
