//
// ------------------------------------------------------------------------------------------------
// Read in N T2K legacy fake data toys and plot the mean total PDF spectra
// Optionally project to 1D and plot RMS stat+syst error envelopes.
// In both cases compare the legacy spectra to input SDK mean and RMS error histograms
//
// Francis Bench
// Sept. 2019
// ------------------------------------------------------------------------------------------------
//

double get_max_in_range_1D(TH1D* hist, double xmax);
double get_nonzero_min_in_range_1D(TH1D* hist, double xmax);
bool sample_is_mu_like(int isample);
void make_palette();
void draw_title(int isample, double min_x, double max_x, double min_y, double max_y);
void draw_frame(TCanvas* canv, bool is_1D, std::string title_x, std::string title_y, std::string title_z, double min_x, double max_x, double min_y, double max_y);
void set_empty_bins_to_one(TH1D* hist);
void set_negative_bins_to_zero(TH1D* hist);

int n_samples = 5;
vector<std::string> legacy_sample_names = {
  "1Rmu",
  "RHC1Rmu",
  "1Re",
  "RHC1Re",
  "MultiRe"
};

vector<std::string> sdk_sample_names = {
  "superk__JPARC_FHC__1rmu",
  "superk__JPARC_RHC__1rmu",
  "superk__JPARC_FHC__1re",
  "superk__JPARC_RHC__1re",
  "superk__JPARC_FHC__multire"
};

//___________________________________________________________________________
void compare_toy_mean_spectra_to_legacy(
  const char * filename_legacy,
  const char * filename_sdk,
  int ntoys = 10000,
  bool plot_1D_RMS_error = true,           // If true, projects all spectra to 1D Ereco and plots RMS error bands
  const char * out_tag = "compare_mean_spectra_legacy_SDK",
  const char * image_save_format = "pdf",
  int canv_size_x = 700,
  int canv_size_y = 700,
  double kinevar_x_max_1D = 3.,      // Max value of the kinematic variable drawn on the x-axis (1D plots)
  double kinevar_x_max_2D = 1.25,    // Max value of the kinematic variable drawn on the x-axis (2D plots)
  double kinevar_y_max_2D = 180.,    // Max value of the kinematic variable drawn on the y-axis (2D plots)
  const char * title_x = "#nu Reconstructed Energy (GeV)",
  const char * title_y = "#theta (degrees)"                 // Used for 2D plots only
)
{
  using std::cout;
  using std::endl;
  using std::vector;
  using namespace valor;
  using namespace valor::conv;
  namespace vutil = valor::util;

  std::string save_name = ((std::string)out_tag);
  if(plot_1D_RMS_error) {
    save_name += "_error_band";
  }

  vector<TH1*> mean_hists(n_samples, nullptr);

  // Vector which we'll use to store the content of each bin for each toy in order to generate error bands
  vector<vector<double>> toy_bin_content;   // 1st index = sample; 2nd index = bin; content = running sum of bin content^2

  // Read the legacy toys from a tree
  TFile* file_legacy = new TFile(filename_legacy, "read");
  TTree* tree_legacy = (TTree*) file_legacy->Get("fitnt");
  TH2D* hist_legacy_1Rmu    = nullptr;
  TH2D* hist_legacy_RHC1Rmu = nullptr;
  TH2D* hist_legacy_1Re     = nullptr;
  TH2D* hist_legacy_RHC1Re  = nullptr;
  TH2D* hist_legacy_MultiRe = nullptr;
  vector<TH2D**> legacy_hists = {&hist_legacy_1Rmu, &hist_legacy_RHC1Rmu, &hist_legacy_1Re, &hist_legacy_RHC1Re, &hist_legacy_MultiRe};
  tree_legacy->SetBranchAddress("hErecoToy_1Rmu"   , &hist_legacy_1Rmu);
  tree_legacy->SetBranchAddress("hErecoToy_RHC1Rmu", &hist_legacy_RHC1Rmu);
  tree_legacy->SetBranchAddress("hErecoToy_1Re"    , &hist_legacy_1Re);
  tree_legacy->SetBranchAddress("hErecoToy_RHC1Re" , &hist_legacy_RHC1Re);
  tree_legacy->SetBranchAddress("hErecoToy_MultiRe", &hist_legacy_MultiRe);
  if(ntoys > tree_legacy->GetEntries()) {
    std::cerr << "ERROR: The input legacy tree contains " << tree_legacy->GetEntries()
              << "toys, but you tried to read " << ntoys << " toys." << endl;
    vutil::Die(-1);
  }

  for(int itoy = 0; itoy < ntoys; ++itoy) {
    if(itoy % 1000 == 0) { cout << "Reading legacy toy #" << itoy << endl; }
    tree_legacy->GetEntry(itoy);

    for(int isample = 0; isample < n_samples; ++isample) {
      TH2D* legacy_hist_2D = nullptr;
      TH1D* legacy_hist_1D = nullptr;
      assert( legacy_hists.at(isample));
      assert(*legacy_hists.at(isample));
      bool is_1D = ( sample_is_mu_like(isample) || plot_1D_RMS_error );
      if(is_1D) {
        TString new_hame = TString::Format("%s_%d", legacy_sample_names.at(isample).c_str(), itoy);
        legacy_hist_1D = (TH1D*) (*legacy_hists.at(isample))->ProjectionY(new_hame);  // In the legacy analysis, thetareco is the x-axis, instead of Ereco
      }
      else {
        legacy_hist_2D = (TH2D*) (*legacy_hists.at(isample))->Clone();
        vutil::SwapAxes(legacy_hist_2D);  // In the legacy analysis, thetareco is the x-axis, instead of Ereco
      }

      if(itoy == 0) {
        // Initialise the vector of mean histograms
        if(is_1D) {
          assert(legacy_hist_1D);
          mean_hists.at(isample) = (TH1*) legacy_hist_1D->Clone(legacy_sample_names.at(isample).c_str());

          // Initialise a vector which we'll use to store the content of each bin for each toy in order to generate error bands
          toy_bin_content.push_back( vector<double>( legacy_hist_1D->GetNbinsX(), 0. ) );
        }
        else {
          assert(legacy_hist_2D);
          mean_hists.at(isample) = (TH1*) legacy_hist_2D->Clone(legacy_sample_names.at(isample).c_str());
        }
        mean_hists.at(isample)->SetDirectory(0);
      }
      else {
        // Add the spectrum for this toy
        assert(mean_hists.at(isample));
        if(is_1D) {
          mean_hists.at(isample)->Add(legacy_hist_1D);
        }
        else {
          assert(legacy_hist_2D);
          mean_hists.at(isample)->Add(legacy_hist_2D);
        }
      }

      if(is_1D) {
        for(int ibin_1D = 0; ibin_1D < legacy_hist_1D->GetNbinsX(); ++ibin_1D) {
          toy_bin_content[isample][ibin_1D] += TMath::Power(legacy_hist_1D->GetBinContent(ibin_1D +1), 2.);  // +1 because the 1st hist bin is an underflow
        }
      }

      // Convert sum to mean
      if(itoy == ntoys -1) { mean_hists.at(isample)->Scale(1.0/((double)ntoys)); }

      delete legacy_hist_2D;
      delete legacy_hist_1D;
    } // isample
  } // itoy
  file_legacy->Close();
  delete file_legacy;
  cout << "Finished reading legacy fakedata toys. Will read in the SDK PDF hists" << endl;

  // Read in the SDK hists
  TFile* file_sdk = new TFile(filename_sdk, "read");
  vector<TH1*>  sdk_mean_hists;
  vector<TH1D*> sdk_up_err_hists;
  vector<TH1D*> sdk_down_err_hists;
  for(int isample = 0; isample < n_samples; ++isample) {
    TString hname = TString::Format("%s_%s", sdk_sample_names.at(isample).c_str(), "mean_PDF");
    sdk_mean_hists.push_back((TH1*) file_sdk->Get(hname));
    assert(*--sdk_mean_hists.end());

    if(plot_1D_RMS_error) {
      hname = TString::Format("%s_%s", sdk_sample_names.at(isample).c_str(), "mean_PDF_err_up");
      sdk_up_err_hists.push_back((TH1D*) file_sdk->Get(hname));
      assert(*--sdk_up_err_hists.end());

      hname = TString::Format("%s_%s", sdk_sample_names.at(isample).c_str(), "mean_PDF_err_down");
      sdk_down_err_hists.push_back((TH1D*) file_sdk->Get(hname));
      assert(*--sdk_down_err_hists.end());
    }
  }
  cout << "Finished reading SDK PDF hists. Will Start plotting" << endl;

  // Draw and save plots of the mean spectra to the output file
  for(int isample = 0; isample < n_samples; ++isample) {
    bool is_1D = ( sample_is_mu_like(isample) || plot_1D_RMS_error );

    TCanvas* canv = new TCanvas(sdk_sample_names.at(isample).c_str(), sdk_sample_names.at(isample).c_str(), canv_size_x, canv_size_y);


    // Draw and save hists
    double max_y = 0.;
    double min_x = 0.;
    double max_x = 0.;
    if(is_1D) {
      canv->cd();
      TH1D* mean_hist     = (TH1D*) mean_hists.at(isample)->Clone();
      TH1D* err_up_hist   = (TH1D*) mean_hist->Clone();
      TH1D* err_down_hist = (TH1D*) mean_hist->Clone();
      TLegend* leg = new TLegend(0.70, 0.80, 0.90, 0.95);
      err_up_hist  ->Reset();
      err_down_hist->Reset();

      if(plot_1D_RMS_error) {
        // Calculate the RMS error for each bin
        for(int ibin = 0; ibin < mean_hist->GetNbinsX(); ++ibin) {
          double n = (double) ntoys;
          double mean = mean_hist->GetBinContent(ibin +1);
          double rms_sq = toy_bin_content[isample][ibin] / n;
          double err = TMath::Sqrt( (n/(n-1.)) * (rms_sq - TMath::Power(mean,2)) );   // Sample standard deviation

          err_up_hist  ->SetBinContent(ibin +1, mean + err);  // +1 because the 1st hist bin is an underflow
          err_down_hist->SetBinContent(ibin +1, mean - err);  // +1 because the 1st hist bin is an underflow
        }
        set_negative_bins_to_zero((TH1D*) sdk_mean_hists    .at(isample));
        set_negative_bins_to_zero((TH1D*) sdk_up_err_hists  .at(isample));
        set_negative_bins_to_zero((TH1D*) sdk_down_err_hists.at(isample));
      }

      // Draw within the limits given at the command line
      max_y = TMath::Max(get_max_in_range_1D(mean_hist  , kinevar_x_max_1D),
                         get_max_in_range_1D(err_up_hist, kinevar_x_max_1D));
      if(plot_1D_RMS_error) {
        max_y = TMath::Max(max_y, get_max_in_range_1D(sdk_up_err_hists.at(isample), kinevar_x_max_1D));
      }
      max_y *= 1.2;
      min_x = mean_hist->GetXaxis()->GetBinLowEdge(1);
      if(sample_is_mu_like(isample)) { max_x = kinevar_x_max_1D; }
      else                           { max_x = kinevar_x_max_2D; }

      draw_frame(canv, true, title_x, "Number of Events", "", 0, max_x, 0, max_y);

      mean_hist->SetLineColor(kBlue);
      if(plot_1D_RMS_error) {
        sdk_up_err_hists  .at(isample)->SetFillColor(kRed);
        sdk_up_err_hists  .at(isample)->SetFillStyle(1001);
        sdk_up_err_hists  .at(isample)->SetLineWidth(1);
        sdk_up_err_hists  .at(isample)->DrawClone("hist,same");
        sdk_down_err_hists.at(isample)->SetFillColor(10);
        sdk_down_err_hists.at(isample)->SetFillStyle(1001);
        sdk_down_err_hists.at(isample)->SetLineWidth(1);
        sdk_down_err_hists.at(isample)->DrawClone("hist,same");

        err_up_hist  ->SetLineColor(kBlue);
        err_up_hist  ->SetLineStyle(2);
        err_up_hist  ->SetLineWidth(2);
        err_up_hist  ->DrawClone("hist,same");
        err_down_hist->SetLineColor(kBlue);
        err_down_hist->SetLineStyle(2);
        err_down_hist->SetLineWidth(2);
        err_down_hist->DrawClone("hist,same");

        leg->AddEntry(sdk_mean_hists  .at(isample), "SDK: Mean"        , "l");
        leg->AddEntry(sdk_up_err_hists.at(isample), "SDK: RMS Error"   , "f");
        leg->AddEntry(mean_hist                   , "Legacy: Mean"     , "l");
        leg->AddEntry(err_up_hist                 , "Legacy: RMS Error", "l");
        leg->DrawClone();
      }

      mean_hist->DrawClone("hist,same");
      sdk_mean_hists.at(isample)->DrawClone("hist,same");
      draw_title(isample, min_x, max_x, 0, max_y);

      // Save canvas to the output file
      std::string canv_save_name = save_name + "." + ((std::string)image_save_format);
      if(isample == 0            && n_samples > 1) { canv_save_name += "("; }
      canv->SaveAs(canv_save_name.c_str());
      delete canv;

      // Draw ratio plots (SDK mean / legacy mean, SDK +err / legacy +err, etc.)
      if(plot_1D_RMS_error) {
        TString canv_ratio_name = TString::Format("%s_%s", sdk_sample_names.at(isample).c_str(), "ratio");
        TCanvas* canv_ratio = new TCanvas(canv_ratio_name, canv_ratio_name, canv_size_x, canv_size_y);
        TH1D* mean_ratio_hist     = (TH1D*) sdk_mean_hists    .at(isample)->Clone();
        TH1D* up_err_ratio_hist   = (TH1D*) sdk_up_err_hists  .at(isample)->Clone();
        TH1D* down_err_ratio_hist = (TH1D*) sdk_down_err_hists.at(isample)->Clone();
        TH1D* RMS_ratio_hist      = (TH1D*) up_err_ratio_hist->Clone();
        RMS_ratio_hist->Add(mean_ratio_hist, -1.);
        TH1D* RMS_hist_legacy = (TH1D*) err_up_hist->Clone();
        RMS_hist_legacy->Add(mean_hist, -1.);
        RMS_ratio_hist     ->Divide(RMS_hist_legacy);
        mean_ratio_hist    ->Divide(mean_hist);
        up_err_ratio_hist  ->Divide(err_up_hist);
        down_err_ratio_hist->Divide(err_down_hist);
        set_empty_bins_to_one(mean_ratio_hist);
        set_empty_bins_to_one(up_err_ratio_hist);
        set_empty_bins_to_one(down_err_ratio_hist);
        set_empty_bins_to_one(RMS_ratio_hist);

        up_err_ratio_hist  ->SetLineColor(kBlue);
        down_err_ratio_hist->SetLineColor(kRed);
        RMS_ratio_hist     ->SetLineColor(kGreen);
        up_err_ratio_hist  ->SetLineStyle(2);
        down_err_ratio_hist->SetLineStyle(3);
        RMS_ratio_hist     ->SetLineStyle(4);
        up_err_ratio_hist  ->SetLineWidth(2);
        down_err_ratio_hist->SetLineWidth(2);
        RMS_ratio_hist     ->SetLineWidth(2);
        up_err_ratio_hist  ->SetFillStyle(0);
        down_err_ratio_hist->SetFillStyle(0);
        RMS_ratio_hist     ->SetFillStyle(0);

        // Draw a faint grey line at ratio = 1, to guide the eye
        TLine* line = new TLine(0, 1., max_x, 1.);
        line->SetLineColor(15);

        max_y = TMath::Max(get_max_in_range_1D(mean_ratio_hist  , max_x),
                           get_max_in_range_1D(up_err_ratio_hist, max_x));
        max_y = TMath::Max(max_y, get_max_in_range_1D(RMS_ratio_hist, max_x));
        max_y *= 1.03;
        //double min_y = 0.8 * TMath::Min(get_nonzero_min_in_range_1D(mean_ratio_hist    , max_x),
        //                                get_nonzero_min_in_range_1D(down_err_ratio_hist, max_x));
        double min_y = 0.97 * get_nonzero_min_in_range_1D(mean_ratio_hist, max_x);

        draw_frame(canv_ratio, true, title_x, "SDK / Legacy", "", 0, max_x, min_y, max_y);
        line->DrawClone("same");
        mean_ratio_hist    ->DrawClone("hist,same");
        up_err_ratio_hist  ->DrawClone("hist,same");
        //down_err_ratio_hist->DrawClone("hist,same");  // I've decided to not draw this, as the number of events is usually close to 0, so the ratio is all over the place (a small absolute change is a much larger fractional change to the ratio close to 0 than further away)
        RMS_ratio_hist     ->DrawClone("hist,same");

        //TLegend* leg_ratio = new TLegend(0.70, 0.80, 0.90, 0.95);
        TLegend* leg_ratio = new TLegend(0.70, 0.16, 0.90, 0.31);
        leg_ratio->AddEntry(mean_ratio_hist    , "Mean"      , "l");
        leg_ratio->AddEntry(up_err_ratio_hist  , "Mean + RMS", "l");
        //leg_ratio->AddEntry(down_err_ratio_hist, "Mean - RMS", "l");
        leg_ratio->AddEntry(RMS_ratio_hist     , "RMS", "l");
        leg_ratio->DrawClone();
        draw_title(isample, min_x, max_x, min_y, max_y);

        if(isample == n_samples -1 && n_samples > 1) { canv_save_name += ")"; }
        canv_ratio->SaveAs(canv_save_name.c_str());

        delete mean_ratio_hist;
        delete up_err_ratio_hist;
        delete down_err_ratio_hist;
        delete RMS_ratio_hist;
        delete leg_ratio;
      }

      delete mean_hist;
      delete err_up_hist;
      delete err_down_hist;
      delete leg;
    } // is 1D
    else {
      // Draw the fractional difference between the SDK and legacy spectra
      canv->cd();
      TH2D* mean_hist = (TH2D*) mean_hists.at(isample)->Clone();
      TH2D* mean_hist_sdk = (TH2D*) sdk_mean_hists.at(isample)->Clone();
      TH2D* diff_hist = (TH2D*) mean_hist->Clone();
      assert(mean_hist);
      assert(mean_hist_sdk);

      for(int ibin_x = 0; ibin_x < mean_hist->GetNbinsX(); ++ibin_x) {
        for(int ibin_y = 0; ibin_y < mean_hist->GetNbinsY(); ++ibin_y) {
          double val_legacy = mean_hist    ->GetBinContent(ibin_x, ibin_y);
          double val_sdk    = mean_hist_sdk->GetBinContent(ibin_x, ibin_y);
          double frac_diff = 0.;
          if(val_sdk != 0) {
            frac_diff = (val_sdk - val_legacy) / val_sdk;
          }
          else {
            frac_diff = -val_legacy;
          }
          diff_hist->SetBinContent(ibin_x, ibin_y, frac_diff);
        }
      }

      max_y = kinevar_y_max_2D;
      min_x = diff_hist->GetXaxis()->GetBinLowEdge(1);
      max_x = kinevar_x_max_2D;
      double max_abs_z = TMath::Max(diff_hist->GetMaximum(), TMath::Abs(diff_hist->GetMinimum()));

      draw_frame(canv, false, title_x, title_y, "(SDK - legacy) / SDK", 0, max_x, 0, max_y);

      make_palette();
      diff_hist->GetZaxis()->SetRangeUser(-max_abs_z, max_abs_z);
      diff_hist->DrawClone("colz,same");
      draw_title(isample, min_x, max_x, 0, max_y);

      // Save canvas to the output file
      std::string canv_save_name = save_name + "." + ((std::string)image_save_format);
      if(isample == 0            && n_samples > 1) { canv_save_name += "("; }
      if(isample == n_samples -1 && n_samples > 1) { canv_save_name += ")"; }
      canv->SaveAs(canv_save_name.c_str());

      delete canv;
      delete mean_hist;
      delete mean_hist_sdk;
      delete diff_hist;
    } // is 2D
  } // isample

  file_sdk->Close();
  delete file_sdk;

}

//___________________________________________________________________________
double get_max_in_range_1D(TH1D* hist, double xmax)
{
  double max_n = 0.;
  for(int xbin = 1; xbin <= hist->GetNbinsX(); xbin++) {
    if(hist->GetBinCenter(xbin) < xmax ) {
      max_n = TMath::Max(max_n, hist->GetBinContent(xbin));
    } //<xmax
  } //xbin
  return max_n;
}
//___________________________________________________________________________
double get_nonzero_min_in_range_1D(TH1D* hist, double xmax)
{
  double min_n = 999.;
  for(int xbin = 1; xbin <= hist->GetNbinsX(); xbin++) {
    if(hist->GetBinCenter(xbin) < xmax && hist->GetBinContent(xbin) != 0) {
      min_n = TMath::Min(min_n, hist->GetBinContent(xbin));
    } //<xmax
  } //xbin
  return min_n;
}
//___________________________________________________________________________
bool sample_is_mu_like(int isample)
{
  if(isample < 0 || isample > 4){
    std::cerr << "ERROR: legacy sample index [" << isample << "] out of range" << endl;
  }
  else {
    if(isample < 2) { return true;  }
    else            { return false; }
  }
  return false;
}
//___________________________________________________________________________
void make_palette()
{
  const int npoints = 3;
  Double_t r[npoints] = {0.0, 1.0, 1.0};
  Double_t g[npoints] = {0.0, 1.0, 0.0};
  Double_t b[npoints] = {1.0, 1.0, 0.0};
  Double_t s[npoints] = {0.0, 0.5, 1.0};
  Int_t ncontours = 255;
  TColor::CreateGradientColorTable(npoints, s, r, g, b, ncontours);
  gStyle->SetNumberContours(ncontours);
}
//___________________________________________________________________________
void draw_title(int isample, double min_x, double max_x, double min_y, double max_y)
{
  TString title_str = TString::Format("%s", sdk_sample_names.at(isample).c_str());
  TPaveText * title = new TPaveText(min_x * 1.05, max_y, max_x, max_y + (max_y-min_y) * 0.07 );
  title->SetShadowColor(0);
  title->SetFillColor(0);
  title->SetFillStyle(0);
  title->SetLineColor(0);
  title->AddText(title_str);
  title->SetTextSize(0.035);
  title->Draw("NB,same");
}
//___________________________________________________________________________
void draw_frame(TCanvas* canv, bool is_1D, std::string title_x, std::string title_y, std::string title_z, double min_x, double max_x, double min_y, double max_y)
{
  TVirtualPad * this_pad = canv->cd(1);
  if(is_1D) {
    this_pad->SetLeftMargin(0.2);
  }
  else {
    this_pad->SetLeftMargin(0.16);
    this_pad->SetRightMargin(0.12);
  }
  TH1F * pdf_frame = (TH1F*)this_pad->DrawFrame(min_x, min_y, max_x, max_y);
  pdf_frame->SetDirectory(0);
  pdf_frame->GetXaxis()->SetTitle(title_x.c_str());
  pdf_frame->GetYaxis()->SetTitle(title_y.c_str());
  if(!is_1D) {
    pdf_frame->GetZaxis()->SetTitle(title_z.c_str());
  }
  pdf_frame->Draw();
}
//___________________________________________________________________________
void set_empty_bins_to_one(TH1D* hist)
{
  for(int ibin = 1; ibin < hist->GetNbinsX() +1; ++ibin) {
    if(hist->GetBinContent(ibin) == 0) {
      hist->SetBinContent(ibin, 1);
    }
  }
}
//___________________________________________________________________________
void set_negative_bins_to_zero(TH1D* hist)
{
  for(int ibin = 1; ibin < hist->GetNbinsX() +1; ++ibin) {
    if(hist->GetBinContent(ibin) < 0) {
      hist->SetBinContent(ibin, 0);
    }
  }
}
//___________________________________________________________________________
