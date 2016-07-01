
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

void plot_ereco_on_top_of_many_prediction_1Donly
(
 const char * sk_sample            = "1Re",                         /* options: "1Rmu", "1Re", "MultiRmu" */
 
 const char * rootfile = "data_run1to7/SKdata_Run17.root",
 //const char * rootfile = "data_run1to6/SKdata_Run16.root",
 //const char * rootfile = "../outputs/test/SKdata_Run16.root",

 const bool drawdata = true,

 float myscale = 1.,
 )
{  
  const int Nfiles = 2;
  const char *expfile[Nfiles] = {
    //"../TECH_NOTE_14_3_16/SPECTRA/files/spectra_RHC1Rmu.root",
    //"../TECH_NOTE_14_3_16/SPECTRA/files/spectra_RHC1Re_2d.root",
    //"spectra_prova_RHC1Rmu.root"

    //"../TECH_NOTE_14_3_16/SPECTRA/files/spectrum_dcp0_1Re.root",
    //"../TECH_NOTE_14_3_16/SPECTRA/files/spectrum_dcpPi_1Re.root",
    //"../TECH_NOTE_14_3_16/SPECTRA/files/spectrum_dcp0_RHC1Re.root",
    //"../TECH_NOTE_14_3_16/SPECTRA/files/spectrum_dcpPi_RHC1Re.root",

    // "../TECH_NOTE_14_3_16/SPECTRA/files/spectrum_pseudobf_dcp0_1Re.root",
    // "../TECH_NOTE_14_3_16/SPECTRA/files/spectrum_pseudobf_dcpPi_1Re.root"
    //"../TECH_NOTE_14_3_16/SPECTRA/files/spectrum_pseudobf_dcp0_RHC1Re.root",
    //"../TECH_NOTE_14_3_16/SPECTRA/files/spectrum_pseudobf_dcpPi_RHC1Re.root",
    
    "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_Asimov_A_dcp0_1Re.root",
    "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_Asimov_A_dcpPi_1Re.root",    
    // "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_Asimov_A_dcp0_RHC1Re.root",
    // "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_Asimov_A_dcpPi_RHC1Re.root",
  };
  
  const char *legname[Nfiles] = {
    //"Normal Hierarchy",
    //"Inverted Hierarchy"
    "#delta_{CP} = 0",
    "#delta_{CP} = #pi"
  };


  
  // Get Data
    
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));
  cout << endl;
  cout << "Reading data from: " << rootfile2 << endl; 
  cout << endl;

  SuperkObservation * skobs_fakedata = new SuperkObservation;
  skobs_fakedata->Read(rootfile2.Data(),0,false);

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  TH2D *hErecoVsQ2Data = (TH2D*)skobs_fakedata->TwoDSpectrum(s);

  TH1D *hQ2Data = (TH1D*)hErecoVsQ2Data->ProjectionX()->Clone("hQ2Data");
  TH1D *hErecoData = (TH1D*)hErecoVsQ2Data->ProjectionY()->Clone("hErecoData");


  

  // Get Expectation
  
  TFile *fexp[Nfiles];
  TH1D *hErecoExp[Nfiles];
  TH1D *hQ2Exp[Nfiles];

  TCanvas *cErecoCompare = new TCanvas("cErecoCompare","cErecoCompare");
  
  TLegend *leg_Ereco = new TLegend(0.6, 0.4, 0.9, 0.7);
  leg_Ereco->SetFillColor(0);
  leg_Ereco->SetBorderSize(0);
  leg_Ereco->SetTextSize(0.06);
 
  for(int ifile=0;ifile<Nfiles;ifile++){
    fexp[ifile] = new TFile(expfile[ifile],"READ");
    hErecoExp[ifile] = (TH1D*) fexp[ifile]->Get("hErecoTot");
    hQ2Exp[ifile] = (TH1D*) fexp[ifile]->Get("hQ2Tot");
    
    hErecoExp[ifile]->SetLineWidth(2);
    hErecoExp[ifile]->SetLineColor(ifile+1);

    hErecoExp[ifile]->GetYaxis()->SetTitle("# of entries");

    if(ifile==0) hErecoExp[ifile]->Draw();
    else         hErecoExp[ifile]->Draw("same");

    leg_Ereco->AddEntry(hErecoExp[ifile],legname[ifile],"L");
  }

  if(drawdata){
    hErecoData->Draw("pe same");
    leg_Ereco->AddEntry(hErecoData,"Data","P");
  }
  
  leg_Ereco->Draw();


  cout << endl;
  for(int ifile=0;ifile<Nfiles;ifile++){
    cout << "File: " << legname[ifile] << endl;
    cout << " - # of events = " << hErecoExp[ifile]->Integral() << endl;
    
    int bin_0p6GeV = hErecoExp[ifile]->FindBin(0.6);
    cout << " - # of events (<0.6 GeV) = " << hErecoExp[ifile]->Integral(1,bin_0p6GeV) << endl;
    cout << " - Frac of events (<0.6 GeV) = " << hErecoExp[ifile]->Integral(1,bin_0p6GeV) / hErecoExp[ifile]->Integral() << endl;
  }
  cout << "Data # of events = " << hErecoData->Integral() << endl;
  cout << " - # of events (<0.6 GeV) = " << hErecoData->Integral(1,bin_0p6GeV) << endl;
  cout << " - Frac of events (<0.6 GeV) = " << hErecoData->Integral(1,bin_0p6GeV) / hErecoData->Integral() << endl;
  cout << " - Mean E = " << hErecoData->GetMean() << endl;
  cout << endl;
  
  

  return;
}
