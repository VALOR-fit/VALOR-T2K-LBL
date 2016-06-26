
//example of macro illustrating how to superimpose two histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"


//
// It takes expected histogram from test_sk_ereco_spectrum_nominal_mainmode_nuebarapp.C output
//

using namespace std;
using namespace valor;
using namespace valor::conv;

void plot_bestfit_on_top_of_data_and_prediction
(
 const bool draw_preliminary = false,

 const char * sk_sample            = "RHC1Re",                         /* options: "1Rmu", "1Re", "MultiRmu" */

 // RUN 1-7
 
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_NoOsc_1Rmu.root",
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_NoOsc_1Re.root",
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_NoOsc_RHC1Rmu.root",
 const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_NoOsc_RHC1Re.root",
 
 const char * datafile = "data_run1to7/SKdata_Run17.root",

 // RUN 1-6
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/files/spectra_NoOsc_RHC1Rmu.root",
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/files/spectra_NoOsc_RHC1Re_2d.root",
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/files/spectra_NoOsc_1Re_2d.root",
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/files/spectra_NoOsc_1Rmu.root",
 //const char * datafile = "data_run1to6/SKdata_Run16.root", 
 )
{ 
  const int NBestFitFiles = 1; 
  const char * bestfitfile[NBestFitFiles] = {
    //"../TECH_NOTE_14_3_16/BESTFIT/RUN16/spectrum_bestfit_fitmarg_joint_sinsqth23Vsdcp_react_nh_dcp0.root",
    //"../TECH_NOTE_14_3_16/BESTFIT/RUN16/spectrum_bestfit_fitmarg_joint_sinsqth23Vsdcp_react_nh_dcpPi.root",

    "../TECH_NOTE_14_3_16/BESTFIT/RUN17/spectrum_bestfit_fitmarg_joint_sinsqth23Vsdcp_react_nh_Run17.root",
  };
  
  const char * bestfitlatex[NBestFitFiles] = {
    "Best-fit"
    //"#delta_{CP} = 0",
    //"#delta_{CP} = #pi"
  };


  // Get SK sample
  
  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  
  // Get Best-Fit
  
  TFile *fbestfit[NBestFitFiles];
  TH1D *hBestFit[NBestFitFiles];
  for(int ifile=0;ifile<NBestFitFiles;ifile++){
    fbestfit[ifile] = new TFile(bestfitfile[ifile],"READ");
    if(s==kSKsFCFV1RingELike)          hBestFit[ifile] = (TH1D*) fbestfit[ifile]->Get("h_bestfit_marg_1Re");
    else if(s==kSKsRHCFCFV1RingELike)  hBestFit[ifile] = (TH1D*) fbestfit[ifile]->Get("h_bestfit_marg_RHC1Re");
    else if(s==kSKsFCFV1RingMuLike)    hBestFit[ifile] = (TH1D*) fbestfit[ifile]->Get("h_bestfit_marg_1Rmu");
    else if(s==kSKsRHCFCFV1RingMuLike) hBestFit[ifile] = (TH1D*) fbestfit[ifile]->Get("h_bestfit_marg_RHC1Rmu");
    else{
      cerr << endl;
      cerr << "Sample " << util::AsStringCompact(s) << " is not defined" << endl;
      cerr << endl;
      exit(1);
    }
  }

  // Get Expectation
  
  TFile *fexp = new TFile(expfile,"READ");
  TH1D *hErecoExp = (TH1D*) fexp->Get("hErecoTot");
  TH1D *hQ2Exp = (TH1D*) fexp->Get("hQ2Tot");

  // Get Data
    
  TString datafile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), datafile));
  cout << endl;
  cout << "Reading data from: " << datafile2 << endl; 
  cout << endl;

  SuperkObservation * skobs_fakedata = new SuperkObservation;
  skobs_fakedata->Read(datafile2.Data(),0,false);

  TH2D *hErecoVsQ2Data = (TH2D*)skobs_fakedata->TwoDSpectrum(s);
  TH1D *hQ2Data = (TH1D*)hErecoVsQ2Data->ProjectionX()->Clone("hQ2Data");
  TH1D *hErecoData = (TH1D*)hErecoVsQ2Data->ProjectionY()->Clone("hErecoData");
  
  

  // Print histograms

  cout << endl;
  for(int ifile=0;ifile<NBestFitFiles;ifile++){
    cout << "File: " << bestfitfile[ifile] << endl;	
    for(int ibin=0;ibin<hBestFit[ifile]->GetNbinsX();ibin++){
      cout << hBestFit[ifile]->GetBinContent(ibin+1) << ", ";
    }
    cout << endl;
  }

  // Print # of events
  
  cout << endl;
  cout << "Tot # of events (expected) = " << hErecoExp->Integral() << endl;
  for(int ifile=0;ifile<NBestFitFiles;ifile++){  
    cout << "File: " << bestfitfile[ifile] << endl;
    cout << "Tot # of events (best-fit) = " << hBestFit[ifile]->Integral() << endl;
  }
  cout << "Tot # of events (observed) = " << hErecoVsQ2Data->Integral() << endl;
  cout << endl;

 



  // Draw 
  
  hErecoVsQ2Data->SetMarkerColor(kGray);
  //hErecoVsQ2Data->SetMarkerStyle(1);
  hErecoVsQ2Data->SetMarkerSize(3);

  TCanvas *Ereco_ExpVsData = new TCanvas("Ereco_ExpVsData","Ereco_ExpVsData");
  hErecoExp->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoExp->GetYaxis()->SetTitle("Number of events");  
  hErecoExp->SetLineColor(kBlue);
  hErecoExp->SetLineWidth(2);
  hErecoExp->DrawCopy("histo");  
  for(int ifile=0;ifile<NBestFitFiles;ifile++){
    hBestFit[ifile]->SetLineColor(2+ifile);
    hBestFit[ifile]->SetLineWidth(2);
    hBestFit[ifile]->Draw("same");
  }
  hErecoData->DrawCopy("pe same");

  TLegend *leg_Ereco = new TLegend(0.6, 0.4, 0.9, 0.7);
  leg_Ereco->SetFillColor(0);
  leg_Ereco->SetBorderSize(0);
  leg_Ereco->AddEntry(hErecoExp,"No oscillation","L");
  leg_Ereco->AddEntry(hErecoData,"Data","P");
  for(int ifile=0;ifile<NBestFitFiles;ifile++){
    leg_Ereco->AddEntry(hBestFit[ifile],bestfitlatex[ifile],"L");
  }
  leg_Ereco->Draw();
    

  TLatex *latPreliminary = new TLatex(hBestFit[0]->GetXaxis()->GetBinLowEdge(5),hBestFit[0]->GetMaximum(),"PRELIMINARY");
  if(draw_preliminary) latPreliminary->Draw("same");


  // Make checks if the # of best-fit files is > 1

  if(NBestFitFiles>1){
    
    cout << endl;

    for(int ifile=0;ifile<NBestFitFiles;ifile++){
    
      int bin_0p6GeV = hBestFit[ifile]->FindBin(0.6);
      double TotEvents = hBestFit[ifile]->Integral();
      double TotEvents_m0p6GeV = hBestFit[ifile]->Integral(1,bin_0p6GeV);
      double FracEvents_m0p6GeV = TotEvents_m0p6GeV / TotEvents;
      double MeanEnergy = hBestFit[ifile]->GetMean();

      cout << bestfitlatex[ifile] << ":" << endl;
      cout << " - Fraction of events below 0.6 GeV = " << FracEvents_m0p6GeV << endl;
      cout << " - Mean Energy = " << MeanEnergy << endl;
    }

    int bin_0p6GeV_data = hErecoData->FindBin(0.6);
    double TotEvents_data = hErecoData->Integral(); 
    double TotEvents_m0p6GeV_data = hErecoData->Integral(1,bin_0p6GeV_data); 
    double FracEvents_m0p6GeV_data = TotEvents_m0p6GeV_data / TotEvents_data;
    double MeanEnergy_data = hErecoData->GetMean();
    cout << "Data:" << endl;
    cout << " - Fraction of events below 0.6 GeV = " << FracEvents_m0p6GeV_data << endl;
    cout << " - Mean Energy = " << MeanEnergy_data << endl;
  }

  return;
}
