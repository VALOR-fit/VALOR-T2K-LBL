//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of 
// cross section tweaking dials.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_sig()
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

  // vary cross sections at far only or at near+far (see cross section err cancelation)
  bool correlate_near_far_sigtwk = true;

  // energy range and binning
  const double Emin  = 0;
  const double Emax  = 4;
  const double dE    = 0.05;
  const int    nbins = (int) ((Emax - Emin)/dE);

  const int ntwkmodes = 3;
  t2kana::yr1numu::conv::SuperkMode_t twkmodes[ntwkmodes] = {
    t2kana::yr1numu::conv::kSKmNuMuCCQE, 
    t2kana::yr1numu::conv::kSKmNuMuNonCCQE, 
    t2kana::yr1numu::conv::kSKmNC    
  };

  const int ntwkdials = 5;
  double twkdial[ntwkdials] = { -1.0, -0.5, 0.0, +0.5, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2,       1,       1      };
  int hcolor [ntwkdials] = { kBlue,   kBlue,   kBlack,  kRed,    kRed   };
  int hstyle [ntwkdials] = { kDotted, kDashed, kSolid,  kDashed, kDotted};

  // histogram booking
  TH1D * hEreco   [ntwkmodes][ntwkdials];
  TH1D * hErecoErr[ntwkmodes][ntwkdials];
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
       hEreco   [imode][idial] = new TH1D("","",nbins,Emin,Emax);
       hEreco   [imode][idial]->SetDirectory(0);
       hEreco   [imode][idial]->SetLineWidth(hwidth[idial]);
       hEreco   [imode][idial]->SetLineColor(hcolor[idial]);
       hEreco   [imode][idial]->SetLineStyle(hstyle[idial]);
       hErecoErr[imode][idial] = new TH1D("","",nbins,Emin,Emax);
       hErecoErr[imode][idial]->SetDirectory(0);
       hErecoErr[imode][idial]->SetLineWidth(hwidth[idial]);
       hErecoErr[imode][idial]->SetLineColor(hcolor[idial]);
       hErecoErr[imode][idial]->SetLineStyle(hstyle[idial]);
    }
  }

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, 0);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        0);

  // generate predictions
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
       skexp.GetGeneratorInputs().TweakSigmaModeSK(twkmodes[imode], twkdial[idial], correlate_near_far_sigtwk);
       for(int ibin=1; ibin<=nbins; ibin++) {
          double Ereco = hEreco[imode][idial]->GetBinCenter(ibin);
          double Nexp  = skexp.Generate(Ereco, *fparm);
          cout << "NumExpSK (Ereco = " << Ereco << ") = " << Nexp << endl;
          hEreco[imode][idial]->SetBinContent(ibin,Nexp);
       }//ibin
    }//idial
  }//imode
  // get fractional err
  double maxerrpos = -1.;
  double maxerrneg = +1.;
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
       for(int ibin=1; ibin<=nbins; ibin++) {
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

  // draw
  TCanvas * c[ntwkmodes];
  gStyle->SetOptStat(0);
  for(int imode = 0; imode < ntwkmodes; imode++) {
    c[imode] = new TCanvas(Form("c%d",imode),"",20,20,500,500);
    c[imode]->Divide(1,2);
    c[imode]->cd(1);
    TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    for(int idial = 0; idial < ntwkdials; idial++) {
       const char * opt   = (idial==0) ? "" : "same";
       const char * label = 
         Form("%s cross section tweaking dial = %.2f", 
             util::AsString(twkmodes[imode]).c_str(), twkdial[idial]);
       if(idial==0) {
        hEreco[imode][idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
        hEreco[imode][idial]->GetYaxis()->SetTitle("Number of events");
       }
       hEreco[imode][idial]->Draw(opt);
       legend->AddEntry(hEreco[imode][idial], label, "L");
       cout << " *** " << label << " : NumEvents = " << hEreco[imode][idial]->Integral() << endl;
    }//idial
    legend->Draw();
    c[imode]->cd(2);
    TH1F * hframe = (TH1F*) c[imode]->DrawFrame(
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
}
