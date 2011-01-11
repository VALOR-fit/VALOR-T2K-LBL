//
// -----------------------------------------------------------------------------
// Plot contour with / without systematics.
//
// Uses the output of run_oscillation_contour_generator.C
// and averages chisq surfaces generated over many fitter runs.
//
// Costas Andreopoulos
// -----------------------------------------------------------------------------
//

#include <iostream>

#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TLine.h>

using namespace std;

void plot_contour(string input_filename_template="./fit-run_%d.root")
{
 const int nruns = 50;

 TH2D * chisq_surf[nruns];       /// average chisq surface per run (calculated in run_oscillation_contour_generator.C)
 TH2D * chisq_surf_syst    = 0;  /// average chisq surface for all runs
 TH2D * chisq_surf_no_syst = 0;  /// chisq surface for fit without systematics

 double bestfit_dm23sq        = 0;
 double bestfit_sinsq2theta23 = 0;

 int ngoodruns = 0;

 for(int irun = 0; irun < nruns; irun++) {

   cout << "Run: " << irun << endl;
   chisq_surf[irun] = 0;

   const char * filename = Form(input_filename_template.c_str(),irun);

   bool file_exists = ! (gSystem->AccessPathName(filename));
   if(!file_exists) continue;

   TFile * f = new TFile(filename, "read");
   if(!f) {
	cout << "Null file??" << endl;
	continue;
   }
   TH2D * h = (TH2D*)f->Get("chisqsurf_avg");
   if(!h) {
	cout << "Null histo??" << endl;
	continue;
   }

   // get avg chisq surfafe from this run
   chisq_surf[irun] = new TH2D(*h);
   for(int i = 1; i <= chisq_surf[irun]->GetNbinsX(); i++) {
      for(int j = 1; j <= chisq_surf[irun]->GetNbinsY(); j++) {
        chisq_surf[irun]->SetBinContent(i,j,h->GetBinContent(i,j));
      }//i
   }//j

   // get no syst contour
   if(irun == 0) {
      TH2D * h2 = (TH2D*)f->Get("chisqsurf_0");
      chisq_surf_no_syst = new TH2D(*h2);
      for(int i = 1; i <= chisq_surf_no_syst->GetNbinsX(); i++) {
         for(int j = 1; j <= chisq_surf_no_syst->GetNbinsY(); j++) {
            chisq_surf_no_syst->SetBinContent(i,j,h2->GetBinContent(i,j));
         }//i
      }//j
   }

   cout << "OK" << endl;

   ngoodruns++;
 }//irun

 cout << "Got all histos" << endl;

 // average evenrything
 chisq_surf_syst = new TH2D(*chisq_surf_no_syst);

 cout << "Creating avg chisq surface" << endl;

 double chisq_avg_min = 99999999;
 for(int i = 1; i <= chisq_surf_syst->GetNbinsX(); i++) {
    for(int j = 1; j <= chisq_surf_syst->GetNbinsY(); j++) {
       double chisq_avg = 0;
        cout << "bin (" << i << ", " << j << ")" << endl;
        for(int irun = 0; irun < nruns; irun++) {
              cout << "Run: " << irun << endl;
	      if(!chisq_surf[irun]) continue;
              assert(chisq_surf[irun]);
              double chisq = chisq_surf[irun]->GetBinContent(i,j);
              chisq_avg += chisq;   
              cout << "OK" << endl;
          }
          cout << "Done" << endl;
          chisq_avg /= ngoodruns;
          chisq_surf_syst->SetBinContent(i,j,chisq_avg);
          chisq_avg_min = TMath::Min(chisq_avg_min, chisq_avg);
      }//j
   }//i     
  for(int i = 1; i <= chisq_surf_syst->GetNbinsX(); i++) {
      for(int j = 1; j <= chisq_surf_syst->GetNbinsY(); j++) {
          chisq_surf_syst->SetBinContent(i,j, chisq_surf_syst->GetBinContent(i,j) - chisq_avg_min);
      }
  }

  double chisq_min = 99999999;
  for(int i = 1; i <= chisq_surf_no_syst->GetNbinsX(); i++) {
      for(int j = 1; j <= chisq_surf_no_syst->GetNbinsY(); j++) {
          double chisq = chisq_surf_no_syst->GetBinContent(i,j);
          if(chisq < chisq_min) {
               chisq_min = chisq;
               bestfit_sinsq2theta23 = chisq_surf_no_syst->GetXaxis()->GetBinCenter(i);
               bestfit_dm23sq        = chisq_surf_no_syst->GetYaxis()->GetBinCenter(j);
          }
      }
  }
  for(int i = 1; i <= chisq_surf_no_syst->GetNbinsX(); i++) {
      for(int j = 1; j <= chisq_surf_no_syst->GetNbinsY(); j++) {
          chisq_surf_no_syst->SetBinContent(i,j, chisq_surf_no_syst->GetBinContent(i,j) - chisq_min);
      }
  }


  const int nlev = 2;
  double levels[nlev] = { 2.30, 4.61 };

  chisq_surf_syst   ->SetContour(nlev,levels);
  chisq_surf_no_syst->SetContour(nlev,levels);
  chisq_surf_no_syst->SetLineColor(kRed);
 
  chisq_surf_syst->Draw("cont2");
  chisq_surf_no_syst->Draw("cont2same");

  TLegend * legend = new TLegend(0.2,0.2,0.6,0.6);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->AddEntry(chisq_surf_syst,    "stat+syst", "L");
  legend->AddEntry(chisq_surf_no_syst, "stat",      "L");
  legend->Draw();


  TMarker * best_fit_point =
          new TMarker(bestfit_sinsq2theta23, bestfit_dm23sq, 29);
  best_fit_point->SetMarkerSize(1.4);
  best_fit_point->Draw();
    
  TString chain_files = input_filename_template;
  chain_files.ReplaceAll("%d","*");
  TChain * chain = new TChain("sysnt");
  chain->Add(chain_files.Data());
  int nfits = chain->GetEntries();
  if(nfits==0) {
    delete chain;
  } else {
    TH1D * hdm23sq        = new TH1D("hdm23sq",        "", 1000, 1E-4,1E-2);
    TH1D * hsinsq2theta23 = new TH1D("hsinsq2theta23", "",  100, 0., 1.);
          
    chain->Draw("dm23sq>>hdm23sq", "", "goff");
    chain->Draw("sinsq2theta23>>hsinsq2theta23", "", "goff");
  
    double dm23sq_err        = hdm23sq       ->GetRMS();
    double sinsq2theta23_err = hsinsq2theta23->GetRMS();
       
    cout << "delta (Deltam2) = " <<  dm23sq_err << endl;
    cout << "delta (Singsq2Theta) = " << sinsq2theta23_err << endl;

    TLine * line1 = new TLine(bestfit_sinsq2theta23-sinsq2theta23_err/2,
                              bestfit_dm23sq,
                              bestfit_sinsq2theta23+sinsq2theta23_err/2,   
                              bestfit_dm23sq);
    TLine * line2 = new TLine(bestfit_sinsq2theta23,
                              bestfit_dm23sq-dm23sq_err/2,
                              bestfit_sinsq2theta23,
                              bestfit_dm23sq+dm23sq_err/2);
    line1->Draw();
    line2->Draw();
  }
          

}
