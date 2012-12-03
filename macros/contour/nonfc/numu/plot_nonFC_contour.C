//-----------------------------------------------------------------------------------------------------
// Plot nonFC 68%, 90%, 95%, 99% CL regions.
// Uses the output of "makeChiSqSurface_nonFC.C" 
// Inputs:
//  - label : 
//          label to be prepended before each legend entry for own confidence regions
//          eg "T2K 1.431E+20 POT"
//  - smooth (boolean):
//          to smooth or not to smooth the \Delta\chi^{2} histograms.
//  - show_contours_with_systematics    : (boolean) Show contours with stat+syst err
//  - show_contours_without_systematics : (boolean) Show contours with syst err only
//  - show_external_contours            : (boolean) Show MINOS/K2K/SK contours
//  - show_68cl : (boolean) Show the 68% CL contour?
//  - show_90cl : (boolean) Show the 90% CL contour?
//  - show_95cl : (boolean) Show the 95% CL contour?
//  - show_99cl : (boolean) Show the 99% CL contour?
//
//@ 11/01/2012 Jelena
//-----------------------------------------------------------------------------------------------------

//
void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style);
void set_contour_style(TGraph*contour, Int_t color, Int_t width, Int_t style);

void pl_nonFC_contour(
            const char * label = 0,
            bool smooth=true,
            bool show_contours_without_systematics=true,
            bool show_contours_with_systematics=true,
            bool show_external_contours=false,
            bool show_68cl=true,
            bool show_90cl=true,
            bool show_95cl=false,
            bool show_99cl=false)
{
      gROOT->Reset();

      // Define dm23sq, sinsq2theta23 range and step  
      double sinsq2theta23_min   = 0.80;
      double sinsq2theta23_max   = 1.0;
      double sinsq2theta23_step  = 0.01;
      double dm23sq_min          = 2.0E-3;
      double dm23sq_max          = 3.4E-3;
      double dm23sq_step         = 1.0E-4;

      // Define grid so that min and max are at the bin centre 
      // and not at the low/high edge of the first/last bin 
      const int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
      const int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );

      double sinsq2theta23_grid_min = sinsq2theta23_min - 0.5*sinsq2theta23_step;
      double sinsq2theta23_grid_max = sinsq2theta23_max + 0.5*sinsq2theta23_step;
      double dm23sq_grid_min        = dm23sq_min        - 0.5*dm23sq_step;
      double dm23sq_grid_max        = dm23sq_max        + 0.5*dm23sq_step;
      
      // Book some histograms (chiSqSurface & confidence region) 
      //stat
      TH2D *hCL68Stat = new TH2D("hCL68Stat", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      TH2D *hCL90Stat = new TH2D("hCL90Stat", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      TH2D *hCL95Stat = new TH2D("hCL95Stat", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      TH2D *hCL99Stat = new TH2D("hCL99Stat", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      //syst
      TH2D *hCL68StatSyst = new TH2D("hCL68StatSyst", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      TH2D *hCL90StatSyst = new TH2D("hCL90StatSyst", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      TH2D *hCL95StatSyst = new TH2D("hCL95StatSyst", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      TH2D *hCL99StatSyst = new TH2D("hCL99StatSyst", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      //==========================

      TFile *dataFile = TFile::Open("ChisqSurface.root", "READ");
      
      dataFile->cd();
      TH2D * hstat = (TH2D*)dataFile->Get("hChi2SurfaceStat");
      TH2D * hsyst = (TH2D*)dataFile->Get("hChi2SurfaceSyst");
      
      TH2D * hstat_min = (TH2D*)dataFile->Get("hChi2SurfaceStat_min");
      TH2D * hsyst_min = (TH2D*)dataFile->Get("hChi2SurfaceSyst_min");
      
      
      
      TH2D * hChisqSurfaceStat = new TH2D(*hstat);
      TH2D * hChisqSurfaceSyst = new TH2D(*hsyst);

      hChisqSurfaceStat->Add(hChi2SurfaceStat_min, -1);
      hChisqSurfaceSyst->Add(hChi2SurfaceStat_min, -1);

//-----------------------------------------------
      //Plot
      for(int sinsq2theta23_bin = 1; sinsq2theta23_bin <= sinsq2theta23_nbins; sinsq2theta23_bin++) {
            for(int dm23sq_bin = 1; dm23sq_bin <= dm23sq_nbins; dm23sq_bin++) {

                  double sinsq2theta23 = hCL68Stat->GetXaxis()->GetBinCenter(sinsq2theta23_bin);
                  double dm23sq        = hCL68Stat->GetYaxis()->GetBinCenter(dm23sq_bin);

                  double dchisq_data_stat=hChisqSurfaceStat->GetBinContent(hChisqSurfaceStat->GetXaxis()->FindBin(sinsq2theta23),hChisqSurfaceStat->GetYaxis()->FindBin(dm23sq));
                  double dchisq_data_sys=hChisqSurfaceSyst->GetBinContent(hChisqSurfaceSyst->GetXaxis()->FindBin(sinsq2theta23),hChisqSurfaceSyst->GetYaxis()->FindBin(dm23sq));

                  double dchisq68_crit = 2.30;
                  double dchisq90_crit = 4.61;
                  double dchisq95_crit = 5.99;
                  double dchisq99_crit = 9.21;
                  //stat
                  hCL68Stat->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data_stat - dchisq68_crit);
                  hCL90Stat->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data_stat - dchisq90_crit);
                  hCL95Stat->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data_stat - dchisq95_crit);
                  hCL99Stat->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data_stat - dchisq99_crit);
                  
                  //syst
                  hCL68StatSyst->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data_sys - dchisq68_crit);
                  hCL90StatSyst->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data_sys - dchisq90_crit);
                  hCL95StatSyst->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data_sys - dchisq95_crit);
                  hCL99StatSyst->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data_sys - dchisq99_crit);

            }//dm23sq_bin
      }//sinsq2theta23_bin

      // Plot contours

      TGaxis::SetMaxDigits(2);

      TCanvas *ConfidenceContour = new TCanvas("ConfidenceContour", "ConfidenceContour", 800, 600);
      ConfidenceContour->Divide(1);
      ConfidenceContour->cd(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);

      TH2D * hframe = (TH2D*) gPad->DrawFrame(sinsq2theta23_min, dm23sq_min, sinsq2theta23_max, dm23sq_max);
      hframe->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
      hframe->GetYaxis()->SetTitle("|#Deltam_{23}^{2}| (eV^{2}/c^{4})");
      hframe->GetXaxis()->SetTitleOffset(1.2);
      hframe->GetYaxis()->SetTitleOffset(1.4);
      hframe->GetXaxis()->CenterTitle();
      hframe->GetYaxis()->CenterTitle();
      hframe->Draw();

      set_contour_style(hCL68Stat, kRed,   2, 2);//colour, line width, line style
      set_contour_style(hCL90Stat, kBlue,  2, 2);
      set_contour_style(hCL95Stat, kGreen, 2, 2);
      set_contour_style(hCL99Stat, kBlack, 2, 2);

      set_contour_style(hCL68StatSyst, kRed,   2, 1);
      set_contour_style(hCL90StatSyst, kBlue,  2, 1);
      set_contour_style(hCL95StatSyst, kGreen, 2, 1);
      set_contour_style(hCL99StatSyst, kBlack, 2, 1);
      
      //smoothing
      if(smooth){
            hCL99Stat -> Smooth();  
            hCL95Stat -> Smooth();  
            hCL90Stat -> Smooth();  
            hCL68Stat -> Smooth();  
            hCL99StatSyst -> Smooth();
            hCL95StatSyst -> Smooth();  
            hCL90StatSyst -> Smooth();  
            hCL68StatSyst -> Smooth();  
      }

      //plot the contours
      if(show_contours_without_systematics) {     
            if(show_99cl) { hCL99Stat -> Draw("cont3same");}  
            if(show_95cl) { hCL95Stat -> Draw("cont3same");}  
            if(show_90cl) { hCL90Stat -> Draw("cont3same");}  
            if(show_68cl) { hCL68Stat -> Draw("cont3same");}  
      }
      if(show_contours_with_systematics){ 
            if(show_99cl) { hCL99StatSyst -> Draw("cont3same");}
            if(show_95cl) { hCL95StatSyst -> Draw("cont3same");}  
            if(show_90cl) { hCL90StatSyst -> Draw("cont3same");}  
            if(show_68cl) { hCL68StatSyst -> Draw("cont3same");}  
      }

      
      if (show_external_contours){
            TGraph * minos2011_68cl= 0;
            TGraph * minos2011_90cl= 0;
            TGraph * skle_90cl     = 0;
            TGraph * sk_90cl       = 0;
            if(plot) {
                  const char * extcontourfile = Form("%s/inputs/external/contours/numu_disappearance.root", gSystem->Getenv("T2KOSCANA"));
                  TFile f(extcontourfile,"read");
                  minos2011_68cl = new TGraph(* (TGraph*)f.Get("minos2011_68cl") );
                  minos2011_90cl = new TGraph(* (TGraph*)f.Get("minos2011_90cl") );
                  skle_90cl      = new TGraph(* (TGraph*)f.Get("skle_90cl")      );
                  sk_90cl        = new TGraph(* (TGraph*)f.Get("sk_90cl")        );
                  set_contour_style(minos2011_68cl, kBlack,     2, 3);
                  set_contour_style(minos2011_90cl, kBlack,     2, 1);
                  set_contour_style(skle_90cl,      kMagenta-3, 2, 1);
                  set_contour_style(sk_90cl,        kGreen+2,   2, 1);

                  if(show_90cl) {
                        minos2011_90cl->Draw("l");
                        skle_90cl->Draw("l");
                        sk_90cl->Draw("l");
                  }
                  if(show_68cl) {
                        minos2011_68cl->Draw("l");
                  }
                  f.Close();
            }
      }

      // build legend labels
      TString pl_68cl_statonly = "Plain 68% CL (stat)";
      TString pl_90cl_statonly = "Plain 90% CL (stat)";
      TString pl_95cl_statonly = "Plain 95% CL (stat)";
      TString pl_99cl_statonly = "Plain 99% CL (stat)";

      TString pl_68cl = "Plain 68% CL (stat + syst)";
      TString pl_90cl = "Plain 90% CL (stat + syst)";
      TString pl_95cl = "Plain 95% CL (stat + syst)";
      TString pl_99cl = "Plain 99% CL (stat + syst)";

      if(label) {
            pl_68cl_statonly = TString(label) + TString(", ") + pl_68cl_statonly;
            pl_90cl_statonly = TString(label) + TString(", ") + pl_90cl_statonly;
            pl_95cl_statonly = TString(label) + TString(", ") + pl_95cl_statonly;
            pl_99cl_statonly = TString(label) + TString(", ") + pl_99cl_statonly;

            pl_68cl = TString(label) + TString(", ") + pl_68cl;
            pl_90cl = TString(label) + TString(", ") + pl_90cl;
            pl_95cl = TString(label) + TString(", ") + pl_95cl;
            pl_99cl = TString(label) + TString(", ") + pl_99cl;

      }

      TLegend * legend = new TLegend(0.17,0.72,0.36,0.88);
      legend->SetFillColor(0);
      legend->SetTextSize(0.03);
      legend->SetBorderSize(0);

      if(show_contours_without_systematics) {
            if(show_99cl) { legend->AddEntry(hCL99Stat, pl_99cl_statonly.Data(), "L");  }
            if(show_95cl) { legend->AddEntry(hCL95Stat, pl_95cl_statonly.Data(), "L");  }
            if(show_90cl) { legend->AddEntry(hCL90Stat, pl_90cl_statonly.Data(), "L");  }
            if(show_68cl) { legend->AddEntry(hCL68Stat, pl_68cl_statonly.Data(), "L");  }
      }
      if(show_contours_with_systematics) {
            if(show_99cl) { legend->AddEntry(hCL99StatSyst, pl_99cl.Data(), "L");  }
            if(show_95cl) { legend->AddEntry(hCL95StatSyst, pl_95cl.Data(), "L");  }
            if(show_90cl) { legend->AddEntry(hCL90StatSyst, pl_90cl.Data(), "L");  }
            if(show_68cl) { legend->AddEntry(hCL68StatSyst, pl_68cl.Data(), "L");  }
      }
      if(show_external_contours) {
            if(show_90cl){
                  legend->AddEntry(sk_90cl,        "Super-K, 90% CL");
                  legend->AddEntry(skle_90cl,      "Super-K L/E, 90% CL");
                  legend->AddEntry(minos2011_90cl, "MINOS 2011 (7.25E+20 POT), 90% CL");
            }
            if(show_68cl) {legend->AddEntry(minos2011_68cl, "MINOS 2011 (7.25E+20 POT), 68% CL");}
      }

      legend->Draw();

      //Save the plot;
      ConfidenceContour->SaveAs("Contour_nonFC.pdf");
      
}
//=============================================

void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style)
{
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
  contour->SetContour(1);//number of contour levels
  contour->SetContourLevel(0,0.0);//SetContourLevel(int level, double value)
  contour->SetStats(0);
}

void set_contour_style(TGraph *contour, Int_t color, Int_t width, Int_t style)
{
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
}

