//
//-----------------------------------------------------------------------------------------------------
//
// This macro allows superimposing FC contours for multiple fits allowing for different critical 
// \Delta\chi^{2} values for each fit (eg fit to different Run periods etc).
//
// C.Andreopoulos
//
//-----------------------------------------------------------------------------------------------------
//

void set_contour_style (TH2D *contour, Int_t color, Int_t width, Int_t style );
void set_contour_style (TGraph*contour, Int_t color, Int_t width, Int_t style);

void plot_feldman_cousins_contour_nfits_ddchisqcrit
(
  /* Number of fits in file below */
  int nfits,                                           

  /* Text file pointing to results allowing to plot the contour for each fit.
     The input file contains 4 lines per fit.
       - a label for the fit (to appear on the plot legend)
       - chisq_file (with full path)  
       - delta_chisq_crit_file_statonly (with full path)
       - delta_chisq_crit_file_statsyst (with full path)  
     Lines starting with # will be ignored so they can be used for adding comments for each fit.
     A typical input text file might look like the following (more fits can be added):
     -----------
     # Official 2012 disappearance fit for the Run1+2+3 dataset
     T2K Run 1+2+3 (3.010E+20 POT)
     /opt/ppd/t2k/data/2012oa/fits/run123fit.root
     /opt/ppd/t2k/data/2012oa/fc/dchisq_crit_statonly.root
     /opt/ppd/t2k/data/2012oa/fc/dchisq_crit_statsyst.root
     # Official 201a disappearance fit for the Run1+2 dataset
     T2K Run 1+2 (1.45E+20 POT)
     /opt/ppd/t2k/data/2010oa/fits/run12fit.root
     /opt/ppd/t2k/data/2010oa/fc/dchisq_crit_statonly.root
     /opt/ppd/t2k/data/2010oa/fc/dchisq_crit_statsyst.root
     -----------
  */
  const char * fit_result_list,    

  /* \Delta\chi^{2}_{critical} surface smoothing factor */
  int nsmoothings = 1,     

  /* superipose contours from MINOS and SK? */
  bool show_external_contours = false,

  /* */
  bool show_contours_with_systematics    = true,
  bool show_contours_without_systematics = false,

  /* Which CLs to show? More than one could be shown but since multiple expts are superimposed
     and, potentially, contours both with and without systematics may be shown,
     the plots looks best if only one CL is shown (typically 90% CL).
     If more than one CLs is shown, the color scheme needs to be tweaked.
     Currently, contours for each fit are shown with their own color.
     Contours with systematics are shown with a solid line, contours without systematics are
     shown with a dashed line. But, the same color is used for all CLs.
  */
  bool show_68cl = false, 
  bool show_90cl = true, 
  bool show_95cl = false,
  bool show_99cl = false,

  /* show best-fit points? */
  bool show_bestfit_point = true
)
{
  gROOT->Reset();

  //
  // Define dm23sq, mixing23 range and step  
  //
  double mixing23_min   = 0.75;
  double mixing23_max   = 1.00;
  double mixing23_step  = 0.01;
  double dm23sq_min     = 1.5E-3;
  double dm23sq_max     = 4.0E-3;
  double dm23sq_step    = 1.0E-4;

  //
  // Define grid so that min and max are at the bin centre 
  // and not at the low/high edge of the first/last bin                                                                                                                               
  //
  const int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
  const int dm23sq_nbins   = 1 + (int) ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  
  double mixing23_grid_min = mixing23_min - 0.5*mixing23_step;
  double mixing23_grid_max = mixing23_max + 0.5*mixing23_step;
  double dm23sq_grid_min   = dm23sq_min   - 0.5*dm23sq_step;
  double dm23sq_grid_max   = dm23sq_max   + 0.5*dm23sq_step;

  // Max number of fits expected
  const int nmax = 5; // if you increase, assign color below
  if(nfits > nmax) {
     cerr << "Increase the max number of input fits expected" << endl;
  }

  // Colors used for each expt (note: Colors used for external MINOS+SK contours may not used)
  int color[nmax] = {
    kBlue, kRed, kGreen, kYellow-2, kViolet+2
  };

  //
  // Book arrays with info displayed for each fit
  //

  string label[nmax];

  TH2D *hCL68StatOnly[nmax];
  TH2D *hCL90StatOnly[nmax];
  TH2D *hCL95StatOnly[nmax];
  TH2D *hCL99StatOnly[nmax];

  TH2D *hCL68StatSyst[nmax];
  TH2D *hCL90StatSyst[nmax];
  TH2D *hCL95StatSyst[nmax];
  TH2D *hCL99StatSyst[nmax];

  TMarker * bestfit_point [nmax];


  for(int ifit=0; ifit<nfits; ifit++) {
     hCL68StatOnly[ifit]  = new TH2D( Form("hCL68StatOnly_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL90StatOnly[ifit]  = new TH2D( Form("hCL90StatOnly_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL95StatOnly[ifit]  = new TH2D( Form("hCL95StatOnly_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL99StatOnly[ifit]  = new TH2D( Form("hCL99StatOnly_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);

     hCL68StatSyst[ifit]  = new TH2D( Form("hCL68StatSyst_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL90StatSyst[ifit]  = new TH2D( Form("hCL90StatSyst_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL95StatSyst[ifit]  = new TH2D( Form("hCL95StatSyst_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL99StatSyst[ifit] = new TH2D( Form("hCL99StatSyst_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  }


  //
  // Read \Delta\chi^{2} surface from oscillation fit and \Delta\chi^{2}_{crit} surfaces from FC runs
  //
  const int buffsize = 2048; 
  ifstream fitinp(fit_result_list, ios::in);

  for(int ifit=0; ifit<nfits; ifit++) {

    cout << "*** Fit: " << ifit << endl;

    string fit_label                  = "";
    string fit_result_file            = "";
    string dchisq_crit_statonly_file  = "";
    string dchisq_crit_statsyst_file  = "";
    
    while(fitinp.peek() == '#') {
      fitinp.ignore(buffsize, '\n');
    }

    getline(fitinp, fit_label,                 '\n');
    getline(fitinp, fit_result_file,           '\n');
    getline(fitinp, dchisq_crit_statonly_file, '\n');
    getline(fitinp, dchisq_crit_statsyst_file, '\n');

    label[ifit] = fit_label;
    cout << "Fit label: " << label[ifit] << endl;

    cout << "Reading fit result from : " << fit_result_file << endl;
    double chisq_min        = 9999999.;
    double mixing23_bestfit = -1.;
    double dm23sq_bestfit   = -1.;
    TFile *chisqFile = TFile::Open(fit_result_file.c_str(), "READ");
    chisqFile->cd();
    TH2D * h2 = (TH2D*)chisqFile->Get("chisq_surface");
    assert(h2);
    TH2D * fChisqSurface = new TH2D(*h2);
    for(int i = 1; i <= fChisqSurface->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface->GetNbinsY(); j++) {
        double chisq = fChisqSurface->GetBinContent(i,j);
          if(chisq < chisq_min) {
            chisq_min = chisq;
            mixing23_bestfit = fChisqSurface->GetXaxis()->GetBinCenter(i);
            dm23sq_bestfit   = fChisqSurface->GetYaxis()->GetBinCenter(j);
          }
      }//j loop
    }//i loop
    for(int i = 1; i <= fChisqSurface->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface->GetNbinsY(); j++) {
        double chisq = fChisqSurface->GetBinContent(i,j);
        fChisqSurface->SetBinContent(i,j, TMath::Max(0.,chisq-chisq_min));
      }//j loop
    }//i loop
    cout << "chi\2_{min} = " << chisq_min << endl;
    cout << "dm23sq_{bestfit} = " << dm23sq_bestfit << " eV^2" << endl;
    cout << "sinsq2theta_{bestfit} = " << mixing23_bestfit << endl;

    cout << "Reading stat-only \delta\chi^2_{crit} results from : " << dchisq_crit_statonly_file << endl;
    TFile *chisqCritStatOnlyFile = TFile::Open(dchisq_crit_statonly_file.c_str(), "READ");
    chisqCritStatOnlyFile->cd();
    TH2D * hDeltaChisqCrit68StatOnly = new TH2D(*hDeltaChi2Surface68);
    TH2D * hDeltaChisqCrit90StatOnly = new TH2D(*hDeltaChi2Surface90);
    TH2D * hDeltaChisqCrit95StatOnly = new TH2D(*hDeltaChi2Surface95);
    TH2D * hDeltaChisqCrit99StatOnly = new TH2D(*hDeltaChi2Surface99);

    cout << "Reading stat+syst \delta\chi^2_{crit} results from : " << dchisq_crit_statsyst_file << endl;
    TFile *chisqCritStatSystFile = TFile::Open(dchisq_crit_statsyst_file.c_str(), "READ");
    chisqCritStatSystFile->cd();
    TH2D * hDeltaChisqCrit68StatSyst = new TH2D(*hDeltaChi2Surface68);
    TH2D * hDeltaChisqCrit90StatSyst = new TH2D(*hDeltaChi2Surface90);
    TH2D * hDeltaChisqCrit95StatSyst = new TH2D(*hDeltaChi2Surface95);
    TH2D * hDeltaChisqCrit99StatSyst = new TH2D(*hDeltaChi2Surface99);
    
    cout << "Smoothing \delta\chi^2_{crit} results" << endl;
    for(int ism = 0; ism < nsmoothings; ism++) {
       hDeltaChisqCrit68StatOnly -> Smooth();
       hDeltaChisqCrit90StatOnly -> Smooth();
       hDeltaChisqCrit95StatOnly -> Smooth();
       hDeltaChisqCrit99StatOnly -> Smooth();
       hDeltaChisqCrit68StatSyst -> Smooth();
       hDeltaChisqCrit90StatSyst -> Smooth();
       hDeltaChisqCrit95StatSyst -> Smooth();
       hDeltaChisqCrit99StatSyst -> Smooth();
    }

    cout << "Computing FC confidence regions..." << endl;
    for(int mixing23_bin = 1; mixing23_bin <= mixing23_nbins; mixing23_bin++) {
        for(int dm23sq_bin = 1; dm23sq_bin <= dm23sq_nbins; dm23sq_bin++) {

          double mixing23 = hCL68StatOnly[ifit]->GetXaxis()->GetBinCenter(mixing23_bin);
          double dm23sq   = hCL68StatOnly[ifit]->GetYaxis()->GetBinCenter(dm23sq_bin);

          double dchisq_data = fChisqSurface->GetBinContent(  
             fChisqSurface->GetXaxis()->FindBin(mixing23), fChisqSurface->GetYaxis()->FindBin(dm23sq) );

          double dchisq68_crit_statonly = hDeltaChisqCrit68StatOnly->GetBinContent( 
               hDeltaChisqCrit68StatOnly->GetXaxis()->FindBin(mixing23), 
               hDeltaChisqCrit68StatOnly->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statonly = hDeltaChisqCrit90StatOnly->GetBinContent( 
               hDeltaChisqCrit90StatOnly->GetXaxis()->FindBin(mixing23), 
               hDeltaChisqCrit90StatOnly->GetYaxis()->FindBin(dm23sq) );
          double dchisq95_crit_statonly = hDeltaChisqCrit95StatOnly->GetBinContent( 
               hDeltaChisqCrit95StatOnly->GetXaxis()->FindBin(mixing23), 
               hDeltaChisqCrit95StatOnly->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statonly = hDeltaChisqCrit99StatOnly->GetBinContent( 
               hDeltaChisqCrit99StatOnly->GetXaxis()->FindBin(mixing23), 
               hDeltaChisqCrit99StatOnly->GetYaxis()->FindBin(dm23sq) );
          hCL68StatOnly[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statonly);
          hCL90StatOnly[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statonly);
          hCL95StatOnly[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statonly);
          hCL99StatOnly[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statonly);

          double dchisq68_crit_statsyst = hDeltaChisqCrit68StatSyst->GetBinContent( 
               hDeltaChisqCrit68StatSyst->GetXaxis()->FindBin(mixing23), 
               hDeltaChisqCrit68StatSyst->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statsyst = hDeltaChisqCrit90StatSyst->GetBinContent( 
               hDeltaChisqCrit90StatSyst->GetXaxis()->FindBin(mixing23), 
               hDeltaChisqCrit90StatSyst->GetYaxis()->FindBin(dm23sq) );
          double dchisq95_crit_statsyst = hDeltaChisqCrit95StatSyst->GetBinContent( 
               hDeltaChisqCrit95StatSyst->GetXaxis()->FindBin(mixing23), 
               hDeltaChisqCrit95StatSyst->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statsyst = hDeltaChisqCrit99StatSyst->GetBinContent( 
               hDeltaChisqCrit99StatSyst->GetXaxis()->FindBin(mixing23), 
               hDeltaChisqCrit99StatSyst->GetYaxis()->FindBin(dm23sq) );
          hCL68StatSyst[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statsyst);
          hCL90StatSyst[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statsyst);
          hCL95StatSyst[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statsyst);
          hCL99StatSyst[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statsyst);

        }//dm23sq_bin
    }//mixing23_bin

    bestfit_point [ifit] = new TMarker(mixing23_bestfit, dm23sq_bestfit, 29);

  }//ifit


  //
  // Plot contours
  //

  TGaxis::SetMaxDigits(2);

  TCanvas *ConfidenceContour = new TCanvas("ConfidenceContour", "ConfidenceContour", 800, 600);
  ConfidenceContour->Divide(1);
  ConfidenceContour->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  TH2D * hframe = (TH2D*) gPad->DrawFrame(mixing23_min, dm23sq_min, mixing23_max, dm23sq_max);
  hframe->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  hframe->GetYaxis()->SetTitle("|#Deltam_{23}^{2}| (eV^{2}/c^{4})");
//hframe->GetXaxis()->SetLimits(0.,1.002);
//hframe->GetYaxis()->SetLimits(0.001,0.004);
  hframe->GetXaxis()->SetTitleOffset(1.2);
  hframe->GetYaxis()->SetTitleOffset(1.4);
  hframe->GetXaxis()->CenterTitle();
  hframe->GetYaxis()->CenterTitle();
  hframe->Draw();

  for(int ifit=0; ifit<nfits; ifit++) {

     set_contour_style(hCL68StatSyst[ifit],  color[ifit], 2, kSolid);
     set_contour_style(hCL90StatSyst[ifit],  color[ifit], 2, kSolid);
     set_contour_style(hCL95StatSyst[ifit],  color[ifit], 2, kSolid);
     set_contour_style(hCL99StatSyst[ifit],  color[ifit], 2, kSolid);
     set_contour_style(hCL68StatOnly[ifit],  color[ifit], 2, kDashed);
     set_contour_style(hCL90StatOnly[ifit],  color[ifit], 2, kDashed);
     set_contour_style(hCL95StatOnly[ifit],  color[ifit], 2, kDashed);
     set_contour_style(hCL99StatOnly[ifit],  color[ifit], 2, kDashed);

     bestfit_point [ifit] -> SetMarkerSize(1.8);
     bestfit_point [ifit] -> SetMarkerColor(color[ifit]);

     if(show_99cl) {
       if (show_contours_with_systematics   ) { hCL99StatSyst[ifit] -> Draw("cont3same");  }
       if (show_contours_without_systematics) { hCL99StatOnly[ifit] -> Draw("cont3same");  }
     }
     if(show_95cl) {
       if (show_contours_with_systematics   ) { hCL95StatSyst[ifit] -> Draw("cont3same");  }
       if (show_contours_without_systematics) { hCL95StatOnly[ifit] -> Draw("cont3same");  }
     }
     if(show_90cl) {
       if (show_contours_with_systematics   ) { hCL90StatSyst[ifit] -> Draw("cont3same");  }
       if (show_contours_without_systematics) { hCL90StatOnly[ifit] -> Draw("cont3same");  }
     }
     if(show_68cl) {
       if (show_contours_with_systematics   ) { hCL68StatSyst[ifit] -> Draw("cont3same");  }
       if (show_contours_without_systematics) { hCL68StatOnly[ifit] -> Draw("cont3same");  }
     }

     if(show_bestfit_point) {
       bestfit_point [ifit] -> Draw("same");
     }
  }//ifit

  // Plot external contours
  TGraph * minos2012_90cl = 0;
  TGraph * skle2012_90cl  = 0;
  TGraph * skzen2012_90cl = 0;
  if(show_external_contours) {
     const char * extcontourfile = Form("%s/inputs/contours/contours_nont2k.root", gSystem->Getenv("T2KOSCANA"));
     TFile f(extcontourfile,"read");
     minos2012_90cl = new TGraph(* (TGraph*)f.Get("minos2012_combined_90CL") );
     skle2012_90cl  = new TGraph(* (TGraph*)f.Get("sk2012_le2nu_90CL")       );
     skzen2012_90cl = new TGraph(* (TGraph*)f.Get("sk2012_zenith3nu_90CL")   );
     set_contour_style(minos2012_90cl, kBlack,     1, kSolid);
     set_contour_style(skle2012_90cl,  kMagenta-3, 1, kSolid);
     set_contour_style(skzen2012_90cl, kGreen+2,   1, kSolid);
     if(show_90cl) {
       minos2012_90cl->Draw("l");
       skle2012_90cl ->Draw("l");
       skzen2012_90cl->Draw("l");
     }
     f.Close();
  }

  //
  // Add legend
  //

  TLegend * legend = new TLegend(0.17,0.72,0.76,0.9);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  const char * pct = "%";
  for(int ifit=0; ifit<nfits; ifit++) {
     if(show_99cl) {
       if (show_contours_with_systematics   ) { legend -> AddEntry(hCL99StatSyst[ifit], Form("%s, FC 99\% CL (stat+syst)",label[ifit].c_str())); }
       if (show_contours_without_systematics) { legend -> AddEntry(hCL99StatOnly[ifit], Form("%s, FC 99\% CL (stat)",     label[ifit].c_str())); }
     }
     if(show_95cl) {
       if (show_contours_with_systematics   ) { legend -> AddEntry(hCL95StatSyst[ifit], Form("%s, FC 95\% CL (stat+syst)",label[ifit].c_str())); }
       if (show_contours_without_systematics) { legend -> AddEntry(hCL95StatOnly[ifit], Form("%s, FC 95\% CL (stat)",     label[ifit].c_str())); }
     }
     if(show_90cl) {
       if (show_contours_with_systematics   ) { legend -> AddEntry(hCL90StatSyst[ifit], Form("%s, FC 90%s CL (stat+syst)",label[ifit].c_str(),pct)); }
       if (show_contours_without_systematics) { legend -> AddEntry(hCL90StatOnly[ifit], Form("%s, FC 90%s CL (stat)",     label[ifit].c_str(),pct)); }
     }
     if(show_68cl) {
       if (show_contours_with_systematics   ) { legend -> AddEntry(hCL68StatSyst[ifit], Form("%s, FC 68%s CL (stat+syst)",label[ifit].c_str(),pct)); }
       if (show_contours_without_systematics) { legend -> AddEntry(hCL68StatOnly[ifit], Form("%s, FC 68%s CL (stat)",     label[ifit].c_str(),pct)); }
     }
  }//ifit
  if(show_external_contours) {
     if(show_90cl) {
       legend->AddEntry(minos2012_90cl, "MINOS 2012 (10.71#times10^{20} POT + 37.9 kton*yr), 90% CL");
       legend->AddEntry(skzen2012_90cl, "SK1-4 zenith angle, 90% CL");
       legend->AddEntry(skle2012_90cl,  "SK1-4 L/E, 90% CL");
     }
  }
  legend->SetBorderSize(0);
  legend->Draw();

}

void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style)
{
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
  contour->SetContour(1);
  contour->SetContourLevel(0,0.0);
  contour->SetStats(0);
}

void set_contour_style(TGraph *contour, Int_t color, Int_t width, Int_t style)
{
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
}
     
   
