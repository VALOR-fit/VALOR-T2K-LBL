//
//-----------------------------------------------------------------------------------------------------
// An extension of plot_feldman_cousins_contour.C
//
// Allows drawing FC contours for up to 3 different fits allowing for different critical 
// \Delta\chi^{2} values for each fit.
// This macro can be used, for example, in order to superimpose fits from different T2K datasets.
//
// C.Andreopoulos
//
// ChangeLog:
// @ 30/05/2011, Costas
//   Added macro during Run1+2 analysis.            
//
//-----------------------------------------------------------------------------------------------------
//

void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style);
void set_contour_style(TGraph*contour, Int_t color, Int_t width, Int_t style);

void plot_feldman_cousins_contour_3fits
(
  const char * label_1                                 = 0,     /* A label for 1st fit  */
  const char * delta_chisq_critical_file_mc_statonly_1 = 0,     /* Values of \Delta\chi^{2}_{critical} for each \Deltam^{2}, sin^{2}2\theta value for 1st fit (stat err only) */
  const char * delta_chisq_critical_file_mc_statsyst_1 = 0,     /* Values of \Delta\chi^{2}_{critical} for each \Deltam^{2}, sin^{2}2\theta value for 1st fit (stat+syst err) */
  const char * chisq_file_data_1                       = 0,     /* The \chi^{2} surface for 1st fit */
  const char * label_2                                 = 0,     /* ... as above but for the 2nd fit */
  const char * delta_chisq_critical_file_mc_statonly_2 = 0,     /* ... */
  const char * delta_chisq_critical_file_mc_statsyst_2 = 0,     /* ... */
  const char * chisq_file_data_2                       = 0,     /* ... */
  const char * label_3                                 = 0,     /* ... as above but for the 3nd fit */
  const char * delta_chisq_critical_file_mc_statonly_3 = 0,     /* ... */
  const char * delta_chisq_critical_file_mc_statsyst_3 = 0,     /* ... */
  const char * chisq_file_data_3                       = 0,     /* ... */
  int          nsmoothings                             = 1,     /* \Delta\chi^{2}_{critical} surface smoothing factor */
  bool         show_external_contours                  = false,
  bool         show_statsyst                           = true,
  bool         show_statonly                           = false
)
{
  gROOT->Reset();

  //
  // Define dm23sq, sinsq2theta23 range and step  
  //
  double sinsq2theta23_min   = 0.75;
  double sinsq2theta23_max   = 1.00;
  double sinsq2theta23_step  = 0.01;
  double dm23sq_min          = 1.5E-3;
  double dm23sq_max          = 4.0E-3;
  double dm23sq_step         = 1.0E-4;

  //
  // Define grid so that min and max are at the bin centre 
  // and not at the low/high edge of the first/last bin                                                                                                                               
  //
  const int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
  const int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );
  
  double sinsq2theta23_grid_min = sinsq2theta23_min - 0.5*sinsq2theta23_step;
  double sinsq2theta23_grid_max = sinsq2theta23_max + 0.5*sinsq2theta23_step;
  double dm23sq_grid_min        = dm23sq_min        - 0.5*dm23sq_step;
  double dm23sq_grid_max        = dm23sq_max        + 0.5*dm23sq_step;

  bool have1 = (label_1!=0);
  bool have2 = (label_2!=0);
  bool have3 = (label_3!=0);

  //
  // Book confidence region histograms
  //

  TH2D *hCL68StatOnly1 = new TH2D("hCL68StatOnly1", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL90StatOnly1 = new TH2D("hCL90StatOnly1", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL95StatOnly1 = new TH2D("hCL95StatOnly1", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL99StatOnly1 = new TH2D("hCL99StatOnly1", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);

  TH2D *hCL68StatSyst1 = new TH2D("hCL68StatSyst1", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL90StatSyst1 = new TH2D("hCL90StatSyst1", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL95StatSyst1 = new TH2D("hCL95StatSyst1", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL99StatSyst1 = new TH2D("hCL99StatSyst1", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);


  TH2D *hCL68StatOnly2 = new TH2D("hCL68StatOnly2", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL90StatOnly2 = new TH2D("hCL90StatOnly2", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL95StatOnly2 = new TH2D("hCL95StatOnly2", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL99StatOnly2 = new TH2D("hCL99StatOnly2", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);

  TH2D *hCL68StatSyst2 = new TH2D("hCL68StatSyst2", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL90StatSyst2 = new TH2D("hCL90StatSyst2", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL95StatSyst2 = new TH2D("hCL95StatSyst2", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL99StatSyst2 = new TH2D("hCL99StatSyst2", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);


  TH2D *hCL68StatOnly3 = new TH2D("hCL68StatOnly3", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL90StatOnly3 = new TH2D("hCL90StatOnly3", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL95StatOnly3 = new TH2D("hCL95StatOnly3", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL99StatOnly3 = new TH2D("hCL99StatOnly3", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);

  TH2D *hCL68StatSyst3 = new TH2D("hCL68StatSyst3", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL90StatSyst3 = new TH2D("hCL90StatSyst3", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL95StatSyst3 = new TH2D("hCL95StatSyst3", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
  TH2D *hCL99StatSyst3 = new TH2D("hCL99StatSyst3", "",
            sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);


  //
  // Get \Delta\chi^2_{critical}
  //

  // \Delta\chi^2_{critical} for 1st fit

  TH2D * hDeltaChisqCrit68StatOnly1 = 0;
  TH2D * hDeltaChisqCrit90StatOnly1 = 0;
  TH2D * hDeltaChisqCrit95StatOnly1 = 0;
  TH2D * hDeltaChisqCrit99StatOnly1 = 0;
  TH2D * hDeltaChisqCrit68StatSyst1 = 0;
  TH2D * hDeltaChisqCrit90StatSyst1 = 0;
  TH2D * hDeltaChisqCrit95StatSyst1 = 0;
  TH2D * hDeltaChisqCrit99StatSyst1 = 0;

  if(have1) {
    TFile *mcFileStatOnly1 = TFile::Open(delta_chisq_critical_file_mc_statonly_1, "READ");
    mcFileStatOnly1->cd();
    TH2D * hDeltaChisqCrit68StatOnly1 = new TH2D(*hDeltaChi2Surface68);	
    TH2D * hDeltaChisqCrit90StatOnly1 = new TH2D(*hDeltaChi2Surface90);
    TH2D * hDeltaChisqCrit95StatOnly1 = new TH2D(*hDeltaChi2Surface95);	
    TH2D * hDeltaChisqCrit99StatOnly1 = new TH2D(*hDeltaChi2Surface99);

    TFile *mcFileStatSyst1 = TFile::Open(delta_chisq_critical_file_mc_statsyst_1, "READ");
    mcFileStatSyst1->cd();
    TH2D * hDeltaChisqCrit68StatSyst1 = new TH2D(*hDeltaChi2Surface68);	
    TH2D * hDeltaChisqCrit90StatSyst1 = new TH2D(*hDeltaChi2Surface90);
    TH2D * hDeltaChisqCrit95StatSyst1 = new TH2D(*hDeltaChi2Surface95);	
    TH2D * hDeltaChisqCrit99StatSyst1 = new TH2D(*hDeltaChi2Surface99);

    for(int ism = 0; ism < nsmoothings; ism++) {
      hDeltaChisqCrit68StatOnly1 -> Smooth();
      hDeltaChisqCrit90StatOnly1 -> Smooth();
      hDeltaChisqCrit95StatOnly1 -> Smooth();
      hDeltaChisqCrit99StatOnly1 -> Smooth();
      hDeltaChisqCrit68StatSyst1 -> Smooth();
      hDeltaChisqCrit90StatSyst1 -> Smooth();
      hDeltaChisqCrit95StatSyst1 -> Smooth();
      hDeltaChisqCrit99StatSyst1 -> Smooth();
    }
  }

  // \Delta\chi^2_{critical} for 2nd fit

  TH2D * hDeltaChisqCrit68StatOnly2 = 0;
  TH2D * hDeltaChisqCrit90StatOnly2 = 0;
  TH2D * hDeltaChisqCrit95StatOnly2 = 0;
  TH2D * hDeltaChisqCrit99StatOnly2 = 0;
  TH2D * hDeltaChisqCrit68StatSyst2 = 0;
  TH2D * hDeltaChisqCrit90StatSyst2 = 0;
  TH2D * hDeltaChisqCrit95StatSyst2 = 0;
  TH2D * hDeltaChisqCrit99StatSyst2 = 0;

  if(have2) {
    TFile *mcFileStatOnly2 = TFile::Open(delta_chisq_critical_file_mc_statonly_2, "READ");
    mcFileStatOnly2->cd();
    hDeltaChisqCrit68StatOnly2 = new TH2D(*hDeltaChi2Surface68);	
    hDeltaChisqCrit90StatOnly2 = new TH2D(*hDeltaChi2Surface90);
    hDeltaChisqCrit95StatOnly2 = new TH2D(*hDeltaChi2Surface95);	
    hDeltaChisqCrit99StatOnly2 = new TH2D(*hDeltaChi2Surface99);

    TFile *mcFileStatSyst2 = TFile::Open(delta_chisq_critical_file_mc_statsyst_2, "READ");
    mcFileStatSyst2->cd();
    hDeltaChisqCrit68StatSyst2 = new TH2D(*hDeltaChi2Surface68);	
    hDeltaChisqCrit90StatSyst2 = new TH2D(*hDeltaChi2Surface90);
    hDeltaChisqCrit95StatSyst2 = new TH2D(*hDeltaChi2Surface95);	
    hDeltaChisqCrit99StatSyst2 = new TH2D(*hDeltaChi2Surface99);

    for(int ism = 0; ism < nsmoothings; ism++) {
      hDeltaChisqCrit68StatOnly2 -> Smooth();
      hDeltaChisqCrit90StatOnly2 -> Smooth();
      hDeltaChisqCrit95StatOnly2 -> Smooth();
      hDeltaChisqCrit99StatOnly2 -> Smooth();
      hDeltaChisqCrit68StatSyst2 -> Smooth();
      hDeltaChisqCrit90StatSyst2 -> Smooth();
      hDeltaChisqCrit95StatSyst2 -> Smooth();
      hDeltaChisqCrit99StatSyst2 -> Smooth();
    }
  }

  // \Delta\chi^2_{critical} for 3nd fit

  TH2D * hDeltaChisqCrit68StatOnly3 = 0;
  TH2D * hDeltaChisqCrit90StatOnly3 = 0;
  TH2D * hDeltaChisqCrit95StatOnly3 = 0;
  TH2D * hDeltaChisqCrit99StatOnly3 = 0;
  TH2D * hDeltaChisqCrit68StatSyst3 = 0;
  TH2D * hDeltaChisqCrit90StatSyst3 = 0;
  TH2D * hDeltaChisqCrit95StatSyst3 = 0;
  TH2D * hDeltaChisqCrit99StatSyst3 = 0;

  if(have3) {
    TFile *mcFileStatOnly3 = TFile::Open(delta_chisq_critical_file_mc_statonly_3, "READ");
    mcFileStatOnly3->cd();
    hDeltaChisqCrit68StatOnly3 = new TH2D(*hDeltaChi2Surface68);	
    hDeltaChisqCrit90StatOnly3 = new TH2D(*hDeltaChi2Surface90);
    hDeltaChisqCrit95StatOnly3 = new TH2D(*hDeltaChi2Surface95);	
    hDeltaChisqCrit99StatOnly3 = new TH2D(*hDeltaChi2Surface99);

    TFile *mcFileStatSyst3 = TFile::Open(delta_chisq_critical_file_mc_statsyst_2, "READ");
    mcFileStatSyst3->cd();
    hDeltaChisqCrit68StatSyst3 = new TH2D(*hDeltaChi2Surface68);	
    hDeltaChisqCrit90StatSyst3 = new TH2D(*hDeltaChi2Surface90);
    hDeltaChisqCrit95StatSyst3 = new TH2D(*hDeltaChi2Surface95);	
    hDeltaChisqCrit99StatSyst3 = new TH2D(*hDeltaChi2Surface99);

    for(int ism = 0; ism < nsmoothings; ism++) {
      hDeltaChisqCrit68StatOnly3 -> Smooth();
      hDeltaChisqCrit90StatOnly3 -> Smooth();
      hDeltaChisqCrit95StatOnly3 -> Smooth();
      hDeltaChisqCrit99StatOnly3 -> Smooth();
      hDeltaChisqCrit68StatSyst3 -> Smooth();
      hDeltaChisqCrit90StatSyst3 -> Smooth();
      hDeltaChisqCrit95StatSyst3 -> Smooth();
      hDeltaChisqCrit99StatSyst3 -> Smooth();
    }
  }

  //
  // Read \Delta\chi^{2} surface from oscillation fit
  //

  double sinsq2theta23_bestfit1 = -1;
  double dm23sq_bestfit1        = -1;

  if(have1) {
    TFile *dataFile1 = TFile::Open(chisq_file_data_1, "READ");
    dataFile1->cd();
    TH2D * h21 = (TH2D*)dataFile1->Get("chisq_surface");
    TH2D * fChisqSurface1 = new TH2D(*h21);
    double chisq_min1 = 999999.0;
    for(int i = 1; i <= fChisqSurface1->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface1->GetNbinsY(); j++) {
        double chisq = fChisqSurface1->GetBinContent(i,j);
        if(chisq < chisq_min1) {
           chisq_min1             = chisq;
           sinsq2theta23_bestfit1 = fChisqSurface1->GetXaxis()->GetBinCenter(i);
	   dm23sq_bestfit1        = fChisqSurface1->GetYaxis()->GetBinCenter(j);
        }
      }//j loop
    }//i loop
    for(int i = 1; i <= fChisqSurface1->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface1->GetNbinsY(); j++) {
        double chisq = fChisqSurface1->GetBinContent(i,j);
        fChisqSurface1->SetBinContent(i,j, TMath::Max(0.,chisq-chisq_min1));
      }//j loop
    }//i loop
    cout << "chi\2_{min} = " << chisq_min1 << endl;
    cout << "dm23sq_{bestfit} = " << dm23sq_bestfit1 << endl;
    cout << "sinsq2theta_{bestfit} = " << sinsq2theta23_bestfit1 << endl;

    for(int sinsq2theta23_bin = 1; sinsq2theta23_bin <= sinsq2theta23_nbins; sinsq2theta23_bin++) {
        for(int dm23sq_bin = 1; dm23sq_bin <= dm23sq_nbins; dm23sq_bin++) {

          double sinsq2theta23 = hCL68StatOnly1->GetXaxis()->GetBinCenter(sinsq2theta23_bin);
          double dm23sq        = hCL68StatOnly1->GetYaxis()->GetBinCenter(dm23sq_bin);
          double dchisq_data = fChisqSurface1->GetBinContent(  fChisqSurface1->GetXaxis()->FindBin(sinsq2theta23), fChisqSurface1->GetYaxis()->FindBin(dm23sq) );

          double dchisq68_crit_statonly = hDeltaChisqCrit68StatOnly1->GetBinContent( hDeltaChisqCrit68StatOnly1->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit68StatOnly1->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statonly = hDeltaChisqCrit90StatOnly1->GetBinContent( hDeltaChisqCrit90StatOnly1->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit90StatOnly1->GetYaxis()->FindBin(dm23sq) ); 
          double dchisq95_crit_statonly = hDeltaChisqCrit95StatOnly1->GetBinContent( hDeltaChisqCrit95StatOnly1->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit95StatOnly1->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statonly = hDeltaChisqCrit99StatOnly1->GetBinContent( hDeltaChisqCrit99StatOnly1->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit99StatOnly1->GetYaxis()->FindBin(dm23sq) );
          hCL68StatOnly1->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statonly);
	  hCL90StatOnly1->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statonly);
	  hCL95StatOnly1->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statonly);
	  hCL99StatOnly1->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statonly);

          double dchisq68_crit_statsyst = hDeltaChisqCrit68StatSyst1->GetBinContent( hDeltaChisqCrit68StatSyst1->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit68StatSyst1->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statsyst = hDeltaChisqCrit90StatSyst1->GetBinContent( hDeltaChisqCrit90StatSyst1->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit90StatSyst1->GetYaxis()->FindBin(dm23sq) ); 
          double dchisq95_crit_statsyst = hDeltaChisqCrit95StatSyst1->GetBinContent( hDeltaChisqCrit95StatSyst1->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit95StatSyst1->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statsyst = hDeltaChisqCrit99StatSyst1->GetBinContent( hDeltaChisqCrit99StatSyst1->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit99StatSyst1->GetYaxis()->FindBin(dm23sq) );
	  hCL68StatSyst1->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statsyst);
	  hCL90StatSyst1->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statsyst);
	  hCL95StatSyst1->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statsyst);
	  hCL99StatSyst1->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statsyst);

        }//dm23sq_bin
    }//sinsq2theta23_bin
  }

  double sinsq2theta23_bestfit2 = -1;
  double dm23sq_bestfit2        = -1;

  if(have2) {
    TFile *dataFile2 = TFile::Open(chisq_file_data_2, "READ");
    dataFile2->cd();
    TH2D * h22 = (TH2D*)dataFile2->Get("chisq_surface");
    TH2D * fChisqSurface2 = new TH2D(*h22);
    double chisq_min2 = 999999.0;
    for(int i = 1; i <= fChisqSurface2->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface2->GetNbinsY(); j++) {
        double chisq = fChisqSurface2->GetBinContent(i,j);
        if(chisq < chisq_min2) {
           chisq_min2 = chisq;
           sinsq2theta23_bestfit2 = fChisqSurface2->GetXaxis()->GetBinCenter(i);
	   dm23sq_bestfit2        = fChisqSurface2->GetYaxis()->GetBinCenter(j);
        }
      }//j loop
    }//i loop
    for(int i = 1; i <= fChisqSurface2->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface2->GetNbinsY(); j++) {
        double chisq = fChisqSurface2->GetBinContent(i,j);
        fChisqSurface2->SetBinContent(i,j, TMath::Max(0.,chisq-chisq_min2));
      }//j loop
    }//i loop
    cout << "chi\2_{min} = " << chisq_min2 << endl;
    cout << "dm23sq_{bestfit} = " << dm23sq_bestfit2 << endl;
    cout << "sinsq2theta_{bestfit} = " << sinsq2theta23_bestfit2 << endl;

    for(int sinsq2theta23_bin = 1; sinsq2theta23_bin <= sinsq2theta23_nbins; sinsq2theta23_bin++) {
        for(int dm23sq_bin = 1; dm23sq_bin <= dm23sq_nbins; dm23sq_bin++) {

          double sinsq2theta23 = hCL68StatOnly2->GetXaxis()->GetBinCenter(sinsq2theta23_bin);
          double dm23sq        = hCL68StatOnly2->GetYaxis()->GetBinCenter(dm23sq_bin);
          double dchisq_data   = fChisqSurface2->GetBinContent(  fChisqSurface2->GetXaxis()->FindBin(sinsq2theta23), fChisqSurface2->GetYaxis()->FindBin(dm23sq) );

          double dchisq68_crit_statonly = hDeltaChisqCrit68StatOnly2->GetBinContent( hDeltaChisqCrit68StatOnly2->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit68StatOnly2->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statonly = hDeltaChisqCrit90StatOnly2->GetBinContent( hDeltaChisqCrit90StatOnly2->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit90StatOnly2->GetYaxis()->FindBin(dm23sq) ); 
          double dchisq95_crit_statonly = hDeltaChisqCrit95StatOnly2->GetBinContent( hDeltaChisqCrit95StatOnly2->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit95StatOnly2->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statonly = hDeltaChisqCrit99StatOnly2->GetBinContent( hDeltaChisqCrit99StatOnly2->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit99StatOnly2->GetYaxis()->FindBin(dm23sq) );
	  hCL68StatOnly2->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statonly);
	  hCL90StatOnly2->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statonly);
	  hCL95StatOnly2->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statonly);
	  hCL99StatOnly2->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statonly);

          double dchisq68_crit_statsyst = hDeltaChisqCrit68StatSyst2->GetBinContent( hDeltaChisqCrit68StatSyst2->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit68StatSyst2->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statsyst = hDeltaChisqCrit90StatSyst2->GetBinContent( hDeltaChisqCrit90StatSyst2->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit90StatSyst2->GetYaxis()->FindBin(dm23sq) ); 
          double dchisq95_crit_statsyst = hDeltaChisqCrit95StatSyst2->GetBinContent( hDeltaChisqCrit95StatSyst2->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit95StatSyst2->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statsyst = hDeltaChisqCrit99StatSyst2->GetBinContent( hDeltaChisqCrit99StatSyst2->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit99StatSyst2->GetYaxis()->FindBin(dm23sq) );
	  hCL68StatSyst2->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statsyst);
	  hCL90StatSyst2->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statsyst);
	  hCL95StatSyst2->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statsyst);
	  hCL99StatSyst2->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statsyst);

        }//dm23sq_bin
    }//sinsq2theta23_bin
  }//2nd

  double sinsq2theta23_bestfit3 = -1;
  double dm23sq_bestfit3        = -1;

  if(have3) {
    TFile *dataFile3 = TFile::Open(chisq_file_data_3, "READ");
    dataFile3->cd();
    TH2D * h23 = (TH2D*)dataFile3->Get("chisq_surface");
    TH2D * fChisqSurface3 = new TH2D(*h23);
    double chisq_min3 = 999999.0;
    for(int i = 1; i <= fChisqSurface3->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface3->GetNbinsY(); j++) {
        double chisq = fChisqSurface3->GetBinContent(i,j);
        if(chisq < chisq_min3) {
           chisq_min3 = chisq;
           sinsq2theta23_bestfit3 = fChisqSurface3->GetXaxis()->GetBinCenter(i);
	   dm23sq_bestfit3        = fChisqSurface3->GetYaxis()->GetBinCenter(j);
        }
      }//j loop
    }//i loop
    for(int i = 1; i <= fChisqSurface3->GetNbinsX(); i++) {
      for(int j = 1; j <= fChisqSurface3->GetNbinsY(); j++) {
        double chisq = fChisqSurface3->GetBinContent(i,j);
        fChisqSurface3->SetBinContent(i,j, TMath::Max(0.,chisq-chisq_min3));
      }//j loop
    }//i loop
    cout << "chi\2_{min} = " << chisq_min3 << endl;
    cout << "dm23sq_{bestfit} = " << dm23sq_bestfit3 << endl;
    cout << "sinsq2theta_{bestfit} = " << sinsq2theta23_bestfit3 << endl;

    for(int sinsq2theta23_bin = 1; sinsq2theta23_bin <= sinsq2theta23_nbins; sinsq2theta23_bin++) {
        for(int dm23sq_bin = 1; dm23sq_bin <= dm23sq_nbins; dm23sq_bin++) {

          double sinsq2theta23 = hCL68StatOnly3->GetXaxis()->GetBinCenter(sinsq2theta23_bin);
          double dm23sq        = hCL68StatOnly3->GetYaxis()->GetBinCenter(dm23sq_bin);
          double dchisq_data   = fChisqSurface3->GetBinContent(  fChisqSurface3->GetXaxis()->FindBin(sinsq2theta23), fChisqSurface3->GetYaxis()->FindBin(dm23sq) );

          double dchisq68_crit_statonly = hDeltaChisqCrit68StatOnly3->GetBinContent( hDeltaChisqCrit68StatOnly3->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit68StatOnly3->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statonly = hDeltaChisqCrit90StatOnly3->GetBinContent( hDeltaChisqCrit90StatOnly3->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit90StatOnly3->GetYaxis()->FindBin(dm23sq) ); 
          double dchisq95_crit_statonly = hDeltaChisqCrit95StatOnly3->GetBinContent( hDeltaChisqCrit95StatOnly3->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit95StatOnly3->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statonly = hDeltaChisqCrit99StatOnly3->GetBinContent( hDeltaChisqCrit99StatOnly3->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit99StatOnly3->GetYaxis()->FindBin(dm23sq) );
	  hCL68StatOnly3->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statonly);
	  hCL90StatOnly3->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statonly);
	  hCL95StatOnly3->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statonly);
	  hCL99StatOnly3->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statonly);

          double dchisq68_crit_statsyst = hDeltaChisqCrit68StatSyst3->GetBinContent( hDeltaChisqCrit68StatSyst3->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit68StatSyst3->GetYaxis()->FindBin(dm23sq) );
          double dchisq90_crit_statsyst = hDeltaChisqCrit90StatSyst3->GetBinContent( hDeltaChisqCrit90StatSyst3->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit90StatSyst3->GetYaxis()->FindBin(dm23sq) ); 
          double dchisq95_crit_statsyst = hDeltaChisqCrit95StatSyst3->GetBinContent( hDeltaChisqCrit95StatSyst3->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit95StatSyst3->GetYaxis()->FindBin(dm23sq) );
          double dchisq99_crit_statsyst = hDeltaChisqCrit99StatSyst3->GetBinContent( hDeltaChisqCrit99StatSyst3->GetXaxis()->FindBin(sinsq2theta23), hDeltaChisqCrit99StatSyst3->GetYaxis()->FindBin(dm23sq) );
	  hCL68StatSyst3->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq68_crit_statsyst);
	  hCL90StatSyst3->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq90_crit_statsyst);
	  hCL95StatSyst3->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq95_crit_statsyst);
	  hCL99StatSyst3->SetBinContent(sinsq2theta23_bin, dm23sq_bin, dchisq_data - dchisq99_crit_statsyst);

        }//dm23sq_bin
    }//sinsq2theta23_bin
  }


  //
  // Plot contours
  //

  TGaxis::SetMaxDigits(2);

  TCanvas *ConfidenceContour = new TCanvas("ConfidenceContour", "ConfidenceContour", 800, 600);
  ConfidenceContour->Divide(1);
  ConfidenceContour->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  TH2D * hframe = (TH2D*) gPad->DrawFrame(sinsq2theta23_min, dm23sq_min, sinsq2theta23_max, dm23sq_max);
  hframe->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  hframe->GetYaxis()->SetTitle("|#Deltam_{23}^{2}| (eV^{2}/c^{4})");
//hframe->GetXaxis()->SetLimits(0.,1.002);
//hframe->GetYaxis()->SetLimits(0.001,0.004);
  hframe->GetXaxis()->SetTitleOffset(1.2);
  hframe->GetYaxis()->SetTitleOffset(1.4);
  hframe->GetXaxis()->CenterTitle();
  hframe->GetYaxis()->CenterTitle();
  hframe->Draw();

  if(have1) {
    set_contour_style(hCL90StatOnly1, kBlue,     2, 2);
    set_contour_style(hCL90StatSyst1, kBlue,     2, 1);
    if(show_statonly) hCL90StatOnly1->Draw("cont3same");
    if(show_statsyst) hCL90StatSyst1->Draw("cont3same");
  }
  if(have2) {
    set_contour_style(hCL90StatOnly2, kRed,   2, 2);
    set_contour_style(hCL90StatSyst2, kRed,   2, 1);
    if(show_statonly) hCL90StatOnly2->Draw("cont3same");
    if(show_statsyst) hCL90StatSyst2->Draw("cont3same");
  }
  if(have3) {
    set_contour_style(hCL90StatOnly3, kGreen,    2, 2);
    set_contour_style(hCL90StatSyst3, kGreen,    2, 1);
    if(show_statonly) hCL90StatOnly3->Draw("cont3same");
    if(show_statsyst) hCL90StatSyst3->Draw("cont3same");
  }

  // Plot external contours
  TGraph * minos2011_68cl = 0;
  TGraph * minos2011_90cl = 0;
  TGraph * skle_90cl      = 0;
  TGraph * sk_90cl        = 0;
  double sk_90cl_sinsq2theta[27] = {1.0,0.992857,0.991747,0.985714,0.978571,0.975489,0.971429,0.967109,0.964286,0.961218,0.957143,0.956349,
                                    0.952912,0.950173,0.950867,0.95322,0.956397,0.957143,0.962394,0.964286,0.970629,0.971429,0.978571,
                                    0.985714,0.986626,0.992857,1.0};
  double sk_90cl_dmsq[27] = {0.00252342,0.00251499,0.00251215,0.00249424,0.00246071,0.00244087,0.0024114,0.00237162,0.00234248,0.00230433,
                             0.00225104,0.00223895,0.00217543,0.00211371,0.00205374,0.00199547,0.00193885,0.00193125,0.00188384,0.00187013,
                             0.0018304,0.00182693,0.00180104,0.00178044,0.00177846,0.00176761,0.00175943};
  double skle_90cl_sinsq2theta[32] = {1.0,0.994963,0.988889,0.983704,0.979259,0.974815,0.971111,0.966815,0.962815,0.958519,0.954815,
                                      0.950963,0.949259,0.947037,0.944593,0.942963,0.94237,0.943333,0.945556,0.948519,0.952,0.956148,
                                      0.960593,0.964815,0.969482,0.974445,0.979111,0.983556,0.988519,0.993333,0.996667,1.00015};
  double skle_90cl_dmsq[32] = {0.00191758,0.00191758,0.0019271,0.00192948,0.00194146,0.0019511,0.00196929,0.0019754,0.00198766,0.00200248,
                               0.0020224,0.00204379,0.00206286,0.00208596,0.00211718,0.00216086,0.00218372,0.00222464,0.00227054,0.00231309,
                               0.00234916,0.00239022,0.00242299,0.00245013,0.00248679,0.00251776,0.00254597,0.00256177,0.00258087,0.00260655,
                               0.00261301,0.00262273};
  if(show_external_contours) {
     const char * extcontourfile = Form("%s/inputs/external/contours/numu_disappearance.root", gSystem->Getenv("T2KOSCANA"));
     TFile f(extcontourfile,"read");
     minos2011_68cl = new TGraph(* (TGraph*)f.Get("minos2011_68cl") );
     minos2011_90cl = new TGraph(* (TGraph*)f.Get("minos2011_90cl") );
     sk_90cl        = new TGraph(27, sk_90cl_sinsq2theta, sk_90cl_dmsq);
     skle_90cl      = new TGraph(32, skle_90cl_sinsq2theta,skle_90cl_dmsq);
     set_contour_style(minos2011_68cl, kBlack,     2, 3);
     set_contour_style(minos2011_90cl, kBlack,     2, 1);
     set_contour_style(skle_90cl,      kMagenta-3, 2, 1);
     set_contour_style(sk_90cl,        kGreen+2,   2, 1);
//     if(show_68cl) {
//       minos2011_68cl->Draw("l");
//     }
//     if(show_90cl) {
       minos2011_90cl->Draw("l");
//       sk_90cl->Draw("l");
//       skle_90cl->Draw("l");
//     }
     f.Close();
  }

  TLegend * legend = new TLegend(0.17,0.72,0.36,0.9);
  legend->SetFillColor(0);
  if(have1) {
    TString fc_90cl_statonly = TString(label_1) + TString(", Feldman-Cousins 90% CL (stat)");
    TString fc_90cl_statsyst = TString(label_1) + TString(", Feldman-Cousins 90% CL (stat+syst)");
    if(show_statonly) legend->AddEntry(hCL90StatOnly1, fc_90cl_statonly.Data(), "L"); 
    if(show_statsyst) legend->AddEntry(hCL90StatSyst1, fc_90cl_statsyst.Data(), "L");
  }
  if(have2) {
    TString fc_90cl_statonly = TString(label_2) + TString(", Feldman-Cousins 90% CL (stat)");
    TString fc_90cl_statsyst = TString(label_2) + TString(", Feldman-Cousins 90% CL (stat+syst)");
    if(show_statonly) legend->AddEntry(hCL90StatOnly2, fc_90cl_statonly.Data(), "L");
    if(show_statsyst) legend->AddEntry(hCL90StatSyst2, fc_90cl_statsyst.Data(), "L");
  }
  if(have3) {
    TString fc_90cl_statonly = TString(label_3) + TString(", Feldman-Cousins 90% CL (stat)");
    TString fc_90cl_statsyst = TString(label_3) + TString(", Feldman-Cousins 90% CL (stat+syst)");
    if(show_statonly) legend->AddEntry(hCL90StatOnly3, fc_90cl_statonly.Data(), "L");
    if(show_statsyst) legend->AddEntry(hCL90StatSyst3, fc_90cl_statsyst.Data(), "L");
  }

  if(show_external_contours) {
  //legend->AddEntry(minos2011_68cl, "MINOS 2011 (7.25E+20 POT), 68% CL", "L");
    legend->AddEntry(minos2011_90cl, "MINOS 2011 (7.25E+20 POT), 90% CL", "L");
  //legend->AddEntry(sk_90cl,        "SuperK Zenith, 90% CL", "L");
  //legend->AddEntry(skle_90cl,      "SuperK L/E, 90% CL", "L");
  }
  legend->SetBorderSize(0);
  legend->Draw();

  if(have1) {
    TMarker * m1 = new TMarker(sinsq2theta23_bestfit1, dm23sq_bestfit1, 29);
    m1->SetMarkerColor(kBlue);
    m1->SetMarkerSize(1.4);
    m1->SetMarkerStyle(29);
    m1->Draw();
  }
  if(have2) {
    TMarker * m2 = new TMarker(sinsq2theta23_bestfit2, dm23sq_bestfit2, 29);
    m2->SetMarkerColor(kRed);
    m2->SetMarkerSize(1.4);
    m2->SetMarkerStyle(29);
    m2->Draw();
  }
  if(have3) {
    TMarker * m3 = new TMarker(sinsq2theta23_bestfit3, dm23sq_bestfit3, 29);
    m3->SetMarkerColor(kGreen);
    m3->SetMarkerSize(1.4);
    m3->SetMarkerStyle(29);
    m3->Draw();
  }
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
     
   
