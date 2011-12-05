//
//-----------------------------------------------------------------------------------------------------
// A variation of plot_feldman_cousins_contour.C
//
// Allows drawing FC contours for an arbitrary number of fits.
// The same critical values of \Delta\chi^{2} are assumed for all fits 
// (see plot_feldman_cousins_contour_3fits for a macro that doesn't make this assumption).
// Usefull, for example, in order to superimpose contours from fits to N toy-MC experiments
// corresponding to the same exposure but generated with different oscillation parameters.
//
// C.Andreopoulos
//
// ChangeLog:
// @ 29/03/2011, Costas 
//   Added macro during Run1+2 sensitivity studies.
//
//-----------------------------------------------------------------------------------------------------
//

void set_contour_style(TH2D*  contour, Int_t color, Int_t width, Int_t style);
void set_contour_style(TGraph*contour, Int_t color, Int_t width, Int_t style);

void plot_feldman_cousins_contour_nfits
(
  const char * label = 0,                                 /* A label to be added in legend  */
  const char * delta_chisq_critical_file_mc_statonly,     /* Values of \Delta\chi^{2}_{critical} for each \Deltam^{2}, sin^{2}2\theta value for 1st fit (stat err only) */
  const char * delta_chisq_critical_file_mc_statsyst,     /* Values of \Delta\chi^{2}_{critical} for each \Deltam^{2}, sin^{2}2\theta value for 1st fit (stat+syst err) */
  const char * fit_result_list,                           /* Text file pointing to fit results, 1 line per fit - 3 entries per line: `input_mixing23 input_dm23sq  path/to/the/output/fit/file.root' */
  int  nfits,                                             /* Number of gits in above file */
  int  nsmoothings                       = 1,             /* */
  bool show_constant_dchisq_contours     = false,         /* */
  bool show_contours_with_systematics    = true,          /* */
  bool show_contours_without_systematics = true,          /* */
  bool show_external_contours            = false,         /* */
  bool show_68cl                         = true,          /* */
  bool show_90cl                         = true,          /* */
  bool show_95cl                         = false,         /* */
  bool show_99cl                         = false          /* */
)
{
  gROOT->Reset();

  //
  // Define dm23sq, mixing23 range and step  
  //
  double mixing23_min   = 0.3; // 0.0;
  double mixing23_max   = 1.0;
  double mixing23_step  = 0.01;
  double dm23sq_min          = 1.0E-3; //1E-3;
  double dm23sq_max          = 5.0E-3; //6E-3;
  double dm23sq_step         = 1E-4;

  //
  // Define grid so that min and max are at the bin centre 
  // and not at the low/high edge of the first/last bin                                                                                                                               
  //
  const int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
  const int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );
  
  double mixing23_grid_min = mixing23_min - 0.5*mixing23_step;
  double mixing23_grid_max = mixing23_max + 0.5*mixing23_step;
  double dm23sq_grid_min        = dm23sq_min        - 0.5*dm23sq_step;
  double dm23sq_grid_max        = dm23sq_max        + 0.5*dm23sq_step;

  //
  // Book confidence region histograms
  //

  const int nmax = 100;

  TH2D *hCL68StatOnlyNonFC[nmax];
  TH2D *hCL90StatOnlyNonFC[nmax]; 
  TH2D *hCL95StatOnlyNonFC[nmax]; 
  TH2D *hCL99StatOnlyNonFC[nmax]; 

  TH2D *hCL68StatOnly[nmax]; 
  TH2D *hCL90StatOnly[nmax]; 
  TH2D *hCL95StatOnly[nmax]; 
  TH2D *hCL99StatOnly[nmax]; 

  TH2D *hCL68StatSyst[nmax]; 
  TH2D *hCL90StatSyst[nmax]; 
  TH2D *hCL95StatSyst[nmax]; 
  TH2D *hCL99StatSyst[nmax]; 

  TMarker * input_point   [nmax];
  TMarker * bestfit_point [nmax];

  for(int ifit=0; ifit<nfits; ifit++) {
     hCL68StatOnlyNonFC[ifit]  = new TH2D( Form("hCL68StatOnlyNonFC_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL90StatOnlyNonFC[ifit]  = new TH2D( Form("hCL90StatOnlyNonFC_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL95StatOnlyNonFC[ifit]  = new TH2D( Form("hCL95StatOnlyNonFC_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
     hCL99StatOnlyNonFC[ifit]  = new TH2D( Form("hCL99StatOnlyNonFC_%d",ifit), "",
            mixing23_nbins, mixing23_grid_min, mixing23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);

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
  //
  //

  TFile *mcFileStatOnly = TFile::Open(delta_chisq_critical_file_mc_statonly, "READ");
  mcFileStatOnly->cd();
  TH2D * hDeltaChisqCrit68StatOnly = new TH2D(*hDeltaChi2Surface68);	
  TH2D * hDeltaChisqCrit90StatOnly = new TH2D(*hDeltaChi2Surface90);
  TH2D * hDeltaChisqCrit95StatOnly = new TH2D(*hDeltaChi2Surface95);	
  TH2D * hDeltaChisqCrit99StatOnly = new TH2D(*hDeltaChi2Surface99);

  TFile *mcFileStatSyst = TFile::Open(delta_chisq_critical_file_mc_statsyst, "READ");
  mcFileStatSyst->cd();
  TH2D * hDeltaChisqCrit68StatSyst = new TH2D(*hDeltaChi2Surface68);	
  TH2D * hDeltaChisqCrit90StatSyst = new TH2D(*hDeltaChi2Surface90);
  TH2D * hDeltaChisqCrit95StatSyst = new TH2D(*hDeltaChi2Surface95);	
  TH2D * hDeltaChisqCrit99StatSyst = new TH2D(*hDeltaChi2Surface99);

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

  //
  // Read \Delta\chi^{2} surface from oscillation fit
  //

  ifstream fitinp(fit_result_list, ios::in);

  for(int ifit=0; ifit<nfits; ifit++) {

    double chisq_min             = 999999.0;
    double mixing23_bestfit = -1;
    double dm23sq_bestfit        = -1;
    double mixing23_input   = -1;
    double dm23sq_input          = -1;    
    string fit_result_file       = "";

    fitinp >> mixing23_input >>  dm23sq_input >> fit_result_file;

    cout << "Reading fit result from : " << fit_result_file << endl;

    TFile *dataFile = TFile::Open(fit_result_file.c_str(), "READ");
    dataFile->cd();
    TH2D * h2 = (TH2D*)dataFile->Get("chisq_surface");
    TH2D * fChisqSurface = new TH2D(*h2);

    for(int i = 1; i <= fChisqSurface->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface->GetNbinsY(); j++) {
        double chisq = fChisqSurface->GetBinContent(i,j);
          if(chisq < chisq_min) {
            chisq_min = chisq;
            mixing23_bestfit = fChisqSurface->GetXaxis()->GetBinCenter(i);
	    dm23sq_bestfit        = fChisqSurface->GetYaxis()->GetBinCenter(j);
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
    cout << "dm23sq_{bestfit} = " << dm23sq_bestfit << endl;
    cout << "sinsq2theta_{bestfit} = " << mixing23_bestfit << endl;

    for(int mixing23_bin = 1; mixing23_bin <= mixing23_nbins; mixing23_bin++) {
        for(int dm23sq_bin = 1; dm23sq_bin <= dm23sq_nbins; dm23sq_bin++) {

          double mixing23 = hCL68StatOnly[ifit]->GetXaxis()->GetBinCenter(mixing23_bin);
          double dm23sq        = hCL68StatOnly[ifit]->GetYaxis()->GetBinCenter(dm23sq_bin);

          double dchisq_data = fChisqSurface->GetBinContent(  fChisqSurface->GetXaxis()->FindBin(mixing23), fChisqSurface->GetYaxis()->FindBin(dm23sq) );

          double dchisq68_crit_statonly_nonfc = 2.30;
          double dchisq90_crit_statonly_nonfc = 4.61;
          double dchisq95_crit_statonly_nonfc = 5.99;
          double dchisq99_crit_statonly_nonfc = 9.21;
	  hCL68StatOnlyNonFC[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statonly_nonfc);
	  hCL90StatOnlyNonFC[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statonly_nonfc);
	  hCL95StatOnlyNonFC[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statonly_nonfc);
	  hCL99StatOnlyNonFC[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statonly_nonfc);

          double dchisq68_crit_statonly = hDeltaChisqCrit68StatOnly->GetBinContent( hDeltaChisqCrit68StatOnly->GetXaxis()->FindBin(mixing23), hDeltaChisqCrit68StatOnly->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statonly = hDeltaChisqCrit90StatOnly->GetBinContent( hDeltaChisqCrit90StatOnly->GetXaxis()->FindBin(mixing23), hDeltaChisqCrit90StatOnly->GetYaxis()->FindBin(dm23sq) ); 
          double dchisq95_crit_statonly = hDeltaChisqCrit95StatOnly->GetBinContent( hDeltaChisqCrit95StatOnly->GetXaxis()->FindBin(mixing23), hDeltaChisqCrit95StatOnly->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statonly = hDeltaChisqCrit99StatOnly->GetBinContent( hDeltaChisqCrit99StatOnly->GetXaxis()->FindBin(mixing23), hDeltaChisqCrit99StatOnly->GetYaxis()->FindBin(dm23sq) );
          hCL68StatOnly[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statonly);
          hCL90StatOnly[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statonly);
          hCL95StatOnly[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statonly);
          hCL99StatOnly[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statonly);

          double dchisq68_crit_statsyst = hDeltaChisqCrit68StatSyst->GetBinContent( hDeltaChisqCrit68StatSyst->GetXaxis()->FindBin(mixing23), hDeltaChisqCrit68StatSyst->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statsyst = hDeltaChisqCrit90StatSyst->GetBinContent( hDeltaChisqCrit90StatSyst->GetXaxis()->FindBin(mixing23), hDeltaChisqCrit90StatSyst->GetYaxis()->FindBin(dm23sq) ); 
          double dchisq95_crit_statsyst = hDeltaChisqCrit95StatSyst->GetBinContent( hDeltaChisqCrit95StatSyst->GetXaxis()->FindBin(mixing23), hDeltaChisqCrit95StatSyst->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statsyst = hDeltaChisqCrit99StatSyst->GetBinContent( hDeltaChisqCrit99StatSyst->GetXaxis()->FindBin(mixing23), hDeltaChisqCrit99StatSyst->GetYaxis()->FindBin(dm23sq) );
          hCL68StatSyst[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statsyst);
          hCL90StatSyst[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statsyst);
          hCL95StatSyst[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statsyst);
          hCL99StatSyst[ifit]->SetBinContent(mixing23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statsyst);

        }//dm23sq_bin
    }//mixing23_bin

   input_point   [ifit] = new TMarker(mixing23_input,   dm23sq_input,   29);
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

     set_contour_style(hCL68StatOnlyNonFC[ifit], kRed,   2, 3);
     set_contour_style(hCL90StatOnlyNonFC[ifit], kBlue,  2, 3);
     set_contour_style(hCL95StatOnlyNonFC[ifit], kGreen, 2, 3);
     set_contour_style(hCL99StatOnlyNonFC[ifit], kBlack, 2, 3);
     set_contour_style(hCL68StatOnly[ifit],      kRed,   2, 2);
     set_contour_style(hCL90StatOnly[ifit],      kBlue,  2, 2);
     set_contour_style(hCL95StatOnly[ifit],      kGreen, 2, 2);
     set_contour_style(hCL99StatOnly[ifit],      kBlack, 2, 2);
     set_contour_style(hCL68StatSyst[ifit],      kRed,   2, 1);
     set_contour_style(hCL90StatSyst[ifit],      kBlue,  2, 1);
     set_contour_style(hCL95StatSyst[ifit],      kGreen, 2, 1);
     set_contour_style(hCL99StatSyst[ifit],      kBlack, 2, 1);

     input_point   [ifit] -> SetMarkerSize(1.8);
     input_point   [ifit] -> SetMarkerColor(kRed);
     bestfit_point [ifit] -> SetMarkerSize(1.8);
     bestfit_point [ifit] -> SetMarkerColor(kBlack);

     if(show_99cl) {
       if (show_contours_with_systematics   ) { hCL99StatSyst[ifit]      -> Draw("cont3same");  }
       if (show_contours_without_systematics) { hCL99StatOnly[ifit]      -> Draw("cont3same");  }
       if (show_constant_dchisq_contours    ) { hCL99StatOnlyNonFC[ifit] -> Draw("cont3same");  }
     }
     if(show_95cl) {
       if (show_contours_with_systematics   ) { hCL95StatSyst[ifit]      -> Draw("cont3same");  }
       if (show_contours_without_systematics) { hCL95StatOnly[ifit]      -> Draw("cont3same");  }
       if (show_constant_dchisq_contours    ) { hCL95StatOnlyNonFC[ifit] -> Draw("cont3same");  }
     }
     if(show_90cl) {
       if (show_contours_with_systematics   ) { hCL90StatSyst[ifit]      -> Draw("cont3same");  }
       if (show_contours_without_systematics) { hCL90StatOnly[ifit]      -> Draw("cont3same");  }
       if (show_constant_dchisq_contours    ) { hCL90StatOnlyNonFC[ifit] -> Draw("cont3same");  }
     }
     if(show_68cl) {   
       if (show_contours_with_systematics   ) { hCL68StatSyst[ifit]      -> Draw("cont3same");  }
       if (show_contours_without_systematics) { hCL68StatOnly[ifit]      -> Draw("cont3same");  }
       if (show_constant_dchisq_contours    ) { hCL68StatOnlyNonFC[ifit] -> Draw("cont3same");  }
     }

     //input_point   [ifit] -> Draw("same");
     //bestfit_point [ifit] -> Draw("same");

  }//ifit


  // Plot external contours
  TGraph * minos2011_68cl = 0;
  TGraph * minos2011_90cl = 0;
  TGraph * skle_90cl      = 0;
  TGraph * sk_90cl        = 0;
  if(show_external_contours) {
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
     if(show_68cl) {
       minos2011_68cl->Draw("l");
     }
     if(show_90cl) {
       minos2011_90cl->Draw("l");
       skle_90cl->Draw("l");
       sk_90cl->Draw("l");
     }
     f.Close();
  }
    

  // build legend labels
  TString fc_68cl_statsyst = "Feldman-Cousins 68% CL (stat+syst)";
  TString fc_90cl_statsyst = "Feldman-Cousins 90% CL (stat+syst)";
  TString fc_95cl_statsyst = "Feldman-Cousins 95% CL (stat+syst)";
  TString fc_99cl_statsyst = "Feldman-Cousins 99% CL (stat+syst)";
  TString fc_68cl_statonly = "Feldman-Cousins 68% CL (stat)";
  TString fc_90cl_statonly = "Feldman-Cousins 90% CL (stat)";
  TString fc_95cl_statonly = "Feldman-Cousins 95% CL (stat)";
  TString fc_99cl_statonly = "Feldman-Cousins 99% CL (stat)";
  TString pl_68cl_statonly = "Plain (const. #Delta#chi^{2}_{crit.}(=2.30)) 68% CL (stat)";
  TString pl_90cl_statonly = "Plain (const. #Delta#chi^{2}_{crit.}(=4.61)) 90% CL (stat)";
  TString pl_95cl_statonly = "Plain (const. #Delta#chi^{2}_{crit.}(=5.99)) 95% CL (stat)";
  TString pl_99cl_statonly = "Plain (const. #Delta#chi^{2}_{crit.}(=9.21)) 99% CL (stat)";
    
  if(label) {
    fc_68cl_statsyst = TString(label) + TString(", ") + fc_68cl_statsyst;
    fc_90cl_statsyst = TString(label) + TString(", ") + fc_90cl_statsyst;
    fc_95cl_statsyst = TString(label) + TString(", ") + fc_95cl_statsyst;
    fc_99cl_statsyst = TString(label) + TString(", ") + fc_99cl_statsyst;
    fc_68cl_statonly = TString(label) + TString(", ") + fc_68cl_statonly;
    fc_90cl_statonly = TString(label) + TString(", ") + fc_90cl_statonly;
    fc_95cl_statonly = TString(label) + TString(", ") + fc_95cl_statonly;
    fc_99cl_statonly = TString(label) + TString(", ") + fc_99cl_statonly;
    pl_68cl_statonly = TString(label) + TString(", ") + pl_68cl_statonly;
    pl_90cl_statonly = TString(label) + TString(", ") + pl_90cl_statonly;     
    pl_95cl_statonly = TString(label) + TString(", ") + pl_95cl_statonly;
    pl_99cl_statonly = TString(label) + TString(", ") + pl_99cl_statonly;
  }  


  TLegend * legend = new TLegend(0.17,0.72,0.36,0.9);
  legend->SetFillColor(0);

  if(show_contours_with_systematics) {
    if(show_99cl) { legend->AddEntry(hCL99StatSyst[0], fc_99cl_statsyst.Data(), "L");  }
    if(show_95cl) { legend->AddEntry(hCL95StatSyst[0], fc_95cl_statsyst.Data(), "L");  }
    if(show_90cl) { legend->AddEntry(hCL90StatSyst[0], fc_90cl_statsyst.Data(), "L");  }
    if(show_68cl) { legend->AddEntry(hCL68StatSyst[0], fc_68cl_statsyst.Data(), "L");  }
  }
  if(show_contours_without_systematics) {
    if(show_99cl) { legend->AddEntry(hCL99StatOnly[0], fc_99cl_statonly.Data(), "L");  }
    if(show_95cl) { legend->AddEntry(hCL95StatOnly[0], fc_95cl_statonly.Data(), "L");  }
    if(show_90cl) { legend->AddEntry(hCL90StatOnly[0], fc_90cl_statonly.Data(), "L");  }
    if(show_68cl) { legend->AddEntry(hCL68StatOnly[0], fc_68cl_statonly.Data(), "L");  }
  }
  if(show_constant_dchisq_contours) {
    if(show_99cl) { legend->AddEntry(hCL99StatOnlyNonFC[0], pl_68cl_statonly.Data(), "L");  }
    if(show_95cl) { legend->AddEntry(hCL95StatOnlyNonFC[0], pl_68cl_statonly.Data(), "L");  }
    if(show_90cl) { legend->AddEntry(hCL90StatOnlyNonFC[0], pl_68cl_statonly.Data(), "L");  }
    if(show_68cl) { legend->AddEntry(hCL68StatOnlyNonFC[0], pl_68cl_statonly.Data(), "L");  }
  }
  if(show_external_contours) {
    if(show_90cl) {
      legend->AddEntry(minos2011_90cl, "MINOS 2011 (7.25E+20 POT), 90% CL");
    }
    if(show_68cl) {
      legend->AddEntry(minos2011_68cl, "MINOS 2011 (7.25E+20 POT), 68% CL");
    }
    if(show_90cl) {
      legend->AddEntry(sk_90cl,   "Super-K, 90% CL");
      legend->AddEntry(skle_90cl, "Super-K L/E, 90% CL");
    }
  }
  legend->SetBorderSize(0);
  legend->Draw();

//ConfidenceContour->SaveAs("ConfidenceContours.png");
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

