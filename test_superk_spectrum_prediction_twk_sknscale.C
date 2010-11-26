//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of SuperK
// nuisance parameter variations that affect the overall normalization such as fiducial mass,
// or "reconstruction efficiency" (in SK language: effect of Reduction, OD cut, vis. energy cut)
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_norm()
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
  const double Emax  = 10;
  const double dE    = 0.05;
  const int    nbins = (int) ((Emax - Emin)/dE);

  const int ntwkdials = 5;
  double twkdial[ntwkdials] = { -1.0, -0.5, 0.0, +0.5, +1.0 };

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2,       1,       1      };
  int hcolor [ntwkdials] = { kBlue,   kBlue,   kBlack,  kRed,    kRed   };
  int hstyle [ntwkdials] = { kDotted, kDashed, kSolid,  kDashed, kDotted};

  // histogram booking
  TH1D * hEreco_fidmass[ntwkdials];
  TH1D * hEreco_recoeff[ntwkdials];
  for(int idial = 0; idial < ntwkdials; idial++) {
     hEreco_fidmass[idial] = new TH1D("","",nbins,Emin,Emax);
     hEreco_fidmass[idial]->SetDirectory(0);
     hEreco_fidmass[idial]->SetLineWidth(hwidth[idial]);
     hEreco_fidmass[idial]->SetLineColor(hcolor[idial]);
     hEreco_fidmass[idial]->SetLineStyle(hstyle[idial]);
     hEreco_recoeff[idial] = new TH1D("","",nbins,Emin,Emax);
     hEreco_recoeff[idial]->SetDirectory(0);
     hEreco_recoeff[idial]->SetLineWidth(hwidth[idial]);
     hEreco_recoeff[idial]->SetLineColor(hcolor[idial]);
     hEreco_recoeff[idial]->SetLineStyle(hstyle[idial]);
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
     skexp.GetSuperkInputs()->TweakFiducialMass( twkdial[idial] );
     for(int ibin=1; ibin<=nbins; ibin++) {
        double Ereco = hEreco_fidmass[idial]->GetBinCenter(ibin);
        double Nexp  = skexp.Generate(Ereco, *fparm);
        cout << "NumExpSK (Ereco = " << Ereco << ", fiducial_twk = " 
             << twkdial[idial] << ") = " << Nexp << endl;
        hEreco_fidmass[idial]->SetBinContent(ibin,Nexp);
     }//ibin
  }//idial
  //reset
  skexp.GetSuperkInputs()->TweakFiducialMass(0.);

  for(int idial = 0; idial < ntwkdials; idial++) {
     skexp.GetSuperkInputs()->TweakRecoEfficiency( twkdial[idial] );
     for(int ibin=1; ibin<=nbins; ibin++) {
        double Ereco = hEreco_recoeff[idial]->GetBinCenter(ibin);
        double Nexp  = skexp.Generate(Ereco, *fparm);
        cout << "NumExpSK (Ereco = " << Ereco << ", recoeff_twk = " 
             << twkdial[idial] << ") = " << Nexp << endl;
        hEreco_recoeff[idial]->SetBinContent(ibin,Nexp);
     }//ibin
  }//idial
  //reset
  skexp.GetSuperkInputs()->TweakRecoEfficiency(0.);

  // draw
  TCanvas * c = new TCanvas("c","",20,20,900,500);
  c->Divide(1,2);
  gStyle->SetOptStat(0);
  c->cd(1);
  TLegend * legend1 = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);
  for(int idial = 0; idial < ntwkdials; idial++) {
    const char * opt   = (idial==0) ? "" : "same";
    const char * label = Form("SuperK fiducial tweaking dial = %.2f", twkdial[idial]);
    if(idial==0) {
      hEreco_fidmass[idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hEreco_fidmass[idial]->GetYaxis()->SetTitle("Number of events");
    }
    hEreco_fidmass[idial]->Draw(opt);
    legend1->AddEntry(hEreco_fidmass[idial], label, "L");
    cout << " *** " << label << " : NumEvents = " << hEreco_fidmass[idial]->Integral() << endl;
  }//idial
  legend1->Draw();
  c->cd(2);
  TLegend * legend2 = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);
  for(int idial = 0; idial < ntwkdials; idial++) {
    const char * opt   = (idial==0) ? "" : "same";
    const char * label = Form("SuperK reco efficiency tweaking dial = %.2f", twkdial[idial]);
    if(idial==0) {
      hEreco_recoeff[idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hEreco_recoeff[idial]->GetYaxis()->SetTitle("Number of events");
    }
    hEreco_recoeff[idial]->Draw(opt);
    legend2->AddEntry(hEreco_recoeff[idial], label, "L");
    cout << " *** " << label << " : NumEvents = " << hEreco_recoeff[idial]->Integral() << endl;
  }//idial
  legend2->Draw();
  c->Update();
}
