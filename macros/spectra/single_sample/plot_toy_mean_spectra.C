//
// ------------------------------------------------------------------------------------------------
// Read in N fake data toys and plot the mean total PDF spectra
// Optionally project to 1D and plot standard deviation stat+syst error envelopes
//
// Francis Bench
// Sept. 2019
// ------------------------------------------------------------------------------------------------
//

double get_max_in_range_1D(TH1D* hist, double xmax);

//___________________________________________________________________________
void plot_toy_mean_spectra(
  const char * config_file,
  const char * filename,
  int ntoys = 10000,
  bool plot_1D_std_error = true,              // If true, projects all spectra to 1D (according to project_to_kinevar) and plots standard deviation error bands
  const char * project_to_kinevar = "Ereco",  // Only used if plot_1D_std_error = true. See KineSpace.cxx
  const char * out_tag = "mean_spectra",
  const bool by_mode = false,                 // If true, plots the mean model PDF for each mode. As this uses the model, statistical fluctuations are not applied.
  const char * image_save_format = "pdf",
  bool stat_fluct = true,
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

  gStyle->SetPalette(56);   // Palette 56 is "inverted dark-body radiator"

  const char * write_opt       = "data_rf+slices";  // TODO: have the option to use simple binning instead!
  const char * write_opt_model = "pdfcs_rf_pdft_rf";// TODO: have the option to use simple binning instead!
  const char * gen_opt         = "";                // TODO: have the option to use simple binning instead!
  const char * read_opt        = "pdft";            // TODO: have the option to use simple binning instead!
  const char * tag             = "fakedata";
  if(by_mode) {
    read_opt = "pdft_pdfc";
  }

  if(by_mode && !plot_1D_std_error) {
    cout << "ERROR: The by_mode option must only be true if plot_1D_std_error is also true (due to unimplemented write options in GSSFM::Write())" << endl;
    vutil::Die(-1);
  }

  std::string save_name = ((std::string)out_tag);
  if(plot_1D_std_error) {
    save_name += "_error_band";
  }
  if(by_mode) {
    save_name += "_by_mode";
  }

  CompoundFitElements* elements = new CompoundFitElements();
  elements->Read(vutil::CompleteConfigPath(config_file), "nomc nosplines");
  FitModelABC* model = elements->FitModel();
  FitDataABC * data  = elements->FitData();
  CompoundFitModel * cm = dynamic_cast<CompoundFitModel*>(model);
  CompoundFitData  * cd = dynamic_cast<CompoundFitData*> (data);

  cout << "Reading fake data from file [" << filename << "]" << endl;
  util::ToyStream* toy_stream = model->InitToyStream(filename, "read", read_opt);

  vector<TH1*> mean_hists(cd->GetNSubData(), nullptr);
  vector<vector<TH1*>> mean_hists_mode(cd->GetNSubData(), vector<TH1*>(conv::kNReacModes, nullptr));

  // Vector which we'll use to store the content of each bin for each toy in order to generate error bands
  vector<vector<double>> toy_bin_content;   // 1st index = subdata; 2nd index = bin; content = running sum of bin content^2
  vector<vector<vector<double>>> toy_bin_content_mode;   // 1st index = subdata; 2nd index = mode; 3rd index = bin; content = running sum of bin content^2

  int n_reac_modes_to_iterate = 1;
  if(by_mode) { n_reac_modes_to_iterate = conv::kNReacModes; }

  for(int itoy = 0; itoy < ntoys; ++itoy) {
    if(itoy % 1000 == 0) { cout << "Reading toy #" << itoy << endl; }

    // Read a fakedata toy from the input file and write hists to an array
    model->Read(toy_stream);
    TObjArray * data_array = new TObjArray();
    if(!by_mode) {
      model->Write(data_array, tag, write_opt);
    }
    else {
      model->Write(data_array, tag, write_opt_model);
    }

    // Extract the histogram for each subdata and add to the mean histogram
    for(int isubdata = 0; isubdata < cd->GetNSubData(); ++isubdata) {
      GenericSingleSampleFitData * subdata  = dynamic_cast<GenericSingleSampleFitData*>(cd->GetSubData(isubdata));
      GenericSingleSampleFitModel* submodel = nullptr;
      if(by_mode) { submodel = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(isubdata)); }

      for(int im = 0; im < n_reac_modes_to_iterate; ++im) {
        conv::ReacMode_t m = conv::kmUndefined;
        if(by_mode) {
          m = vutil::ReacModeFromInt(im);
          if(m == conv::kmUndefined) { continue; }
          int mode_idx = vutil::ArrIdx::Instance()->ReacModeIdx(submodel->GetBeam(), submodel->GetSample(), m);
          if(mode_idx < 0) { continue; }
        }

        TString ob_name;
        if(!by_mode) {
          if(plot_1D_std_error) { ob_name = TString::Format("%s_data_ob_rf_slice_%s_%s", tag, subdata->Name().c_str(), project_to_kinevar); }
          else {                  ob_name = TString::Format("%s_data_ob_rf_%s", tag, subdata->Name().c_str()); }
        }
        else {
          ob_name = TString::Format("%s_pdfcs_ob_rf_%s_%s_%s", tag, submodel->Name().c_str(), project_to_kinevar, vutil::AsString(m).c_str());
        }

        TH1 * hist_ob = (TH1*) data_array->FindObject(ob_name.Data());
        if(!hist_ob) {
          cerr << "ERROR: could not find hist named [" << ob_name << "] in the array with contents:" << endl;
          data_array->Print();
          vutil::Die(-1);
        }

        if(itoy == 0) {
          // Initialise the vector of mean histograms
          mean_hists.at(isubdata) = (TH1*) hist_ob->Clone(subdata->Name().c_str());
          if(by_mode) {
            TString hname_mean( TString::Format("%s_%s", submodel->Name().c_str(), vutil::AsString(m).c_str()) );
            mean_hists_mode.at(isubdata).at(im) = (TH1*) hist_ob->Clone(hname_mean);
          }

          // Initialise a vector which we'll use to store the content of each bin for each toy in order to generate error bands
          if(plot_1D_std_error) {
            toy_bin_content.push_back( vector<double>( hist_ob->GetNbinsX(), 0. ) );
            if(by_mode) {
              toy_bin_content_mode.push_back( vector<vector<double>>( n_reac_modes_to_iterate, vector<double>( hist_ob->GetNbinsX(), 0. ) ) );
            }
          }
        }
        else {
          // Add the spectrum for this toy
          mean_hists.at(isubdata)->Add(hist_ob);
          if(by_mode) {
            mean_hists_mode.at(isubdata).at(im)->Add(hist_ob);
          }
        }

        if(plot_1D_std_error) {
          for(int ibin_1D = 0; ibin_1D < hist_ob->GetNbinsX(); ++ibin_1D) {
            toy_bin_content[isubdata][ibin_1D] += TMath::Power(hist_ob->GetBinContent(ibin_1D +1), 2.);  // +1 because the 1st hist bin is an underflow
            if(by_mode) {
              toy_bin_content_mode[isubdata][im][ibin_1D] += TMath::Power(hist_ob->GetBinContent(ibin_1D + 1), 2.);  // +1 because the 1st hist bin is an underflow
            }
          }
        }

        // Convert sum to mean
        if(itoy == ntoys -1) {
          mean_hists.at(isubdata)->Scale(1.0/((double)ntoys));
          if(by_mode) {
            mean_hists_mode.at(isubdata).at(im)->Scale(1.0/((double)ntoys));
          }
        }
      } // im
    } // isubdata

    // Free up memory
    for(int i = 0; i < data_array->GetEntries(); ++i) {
      delete data_array->At(i);
    }
    delete data_array;
  } // itoy
  delete toy_stream;
  cout << "Finished reading fakedata toys. Will start plotting" << endl;

  // Draw and save plots of the mean spectra to the output file
  std::string root_save_name = save_name + ".root";
  TFile* root_file_out = new TFile(root_save_name.c_str(), "recreate");
  bool first_saved_plot = true;
  for(int isubdata = 0; isubdata < cd->GetNSubData(); ++isubdata) {
    // Get the kinematic variables used for binning to decide whether to plot in 1D or 2D
    GenericSingleSampleFitData * subdata = dynamic_cast<GenericSingleSampleFitData*>(cd->GetSubData(isubdata));
    GenericSingleSampleFitModel* submodel = nullptr;
    if(by_mode) { submodel = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(isubdata)); }
    const SingleSampleKineBinning* sskb = nullptr;
    if(!by_mode) {  sskb = subdata->GetBinning(); }
    else { sskb = submodel->GetBinning(); }
    const VALORBinning * reco_opt = sskb->RecoOptBinning();
    const vector<KineVar_t> & reco_vars = reco_opt->KineVars();
    if(!(reco_vars.size() == 1 || reco_vars.size() == 2)) {
      std::cerr << "ERROR: The number of reco dimensions to plot was != 1 or 2." << endl;
      vutil::Die(-1);
    }
    const bool is_1D = ( reco_vars.size() == 1 || plot_1D_std_error || by_mode );

    // Get the index of the last used reaction mode, so we know when to close the output PDF file
    int last_mode = 0;
    if(by_mode) {
      for(int im = 0; im < n_reac_modes_to_iterate; ++im) {
        conv::ReacMode_t m = conv::kmUndefined;
        m = vutil::ReacModeFromInt(im);
        if(m == conv::kmUndefined) { continue; }
        int mode_idx = vutil::ArrIdx::Instance()->ReacModeIdx(submodel->GetBeam(), submodel->GetSample(), m);
        if(mode_idx < 0) { continue; }
        last_mode = im;
      }
    }

    // Draw and save hists
    double max_y = 0.;
    double min_x = 0.;
    double max_x = 0.;
    if(is_1D) {
      for(int im = 0; im < n_reac_modes_to_iterate; ++im) {
        conv::ReacMode_t m = conv::kmUndefined;
        if(by_mode) {
          m = vutil::ReacModeFromInt(im);
          if(m == conv::kmUndefined) { continue; }
          int mode_idx = vutil::ArrIdx::Instance()->ReacModeIdx(submodel->GetBeam(), submodel->GetSample(), m);
          if(mode_idx < 0) { continue; }
        }

        TString canv_name(subdata->Name().c_str());
        if(by_mode) {
          canv_name = TString::Format("%s_%s", submodel->Name().c_str(), vutil::AsString(m).c_str());
        }
        TCanvas* canv = new TCanvas(canv_name, canv_name, canv_size_x, canv_size_y);
        canv->cd();

        TH1D* mean_hist = nullptr;
        if(!by_mode) { mean_hist     = (TH1D*) mean_hists     .at(isubdata)       ->Clone(); }
        else         { mean_hist     = (TH1D*) mean_hists_mode.at(isubdata).at(im)->Clone(); }
        TH1D* err_up_hist   = (TH1D*) mean_hist->Clone();
        TH1D* err_down_hist = (TH1D*) mean_hist->Clone();

        TLegend* leg = new TLegend(0.70, 0.80, 0.90, 0.95);
        err_up_hist  ->Reset();
        err_down_hist->Reset();

        if(plot_1D_std_error) {
          // Calculate the std error for each bin
          for(int ibin = 0; ibin < mean_hist->GetNbinsX(); ++ibin) {
            double rms_sq = -1.;
            double n = (double) ntoys;
            double mean = mean_hist->GetBinContent(ibin +1);
            if(!by_mode) {
              rms_sq = toy_bin_content[isubdata][ibin] / n;
            }
            else {
              rms_sq = toy_bin_content_mode[isubdata][im][ibin] / n;
            }
            double err = TMath::Sqrt( (n/(n-1.)) * (rms_sq - TMath::Power(mean,2)) );   // Sample standard deviation

            err_up_hist  ->SetBinContent(ibin +1, mean + err);  // +1 because the 1st hist bin is an underflow
            err_down_hist->SetBinContent(ibin +1, mean - err);  // +1 because the 1st hist bin is an underflow
          }
        }

        // Draw within the limits given at the command line
        max_y = 1.2 * TMath::Max(get_max_in_range_1D(mean_hist  , kinevar_x_max_1D),
                                 get_max_in_range_1D(err_up_hist, kinevar_x_max_1D));
        min_x = mean_hist->GetXaxis()->GetBinLowEdge(1);
        if(reco_vars.size() == 1) { max_x = kinevar_x_max_1D; }
        else { max_x = kinevar_x_max_2D; }
        TVirtualPad * this_pad = canv->cd(1);
        this_pad->SetLeftMargin(0.2);
        TH1F * pdf_frame = (TH1F*)this_pad->DrawFrame(0, 0, max_x, max_y);
        pdf_frame->GetXaxis()->SetTitle(title_x);
        pdf_frame->GetYaxis()->SetTitle("Number of Events");
        pdf_frame->Draw();

        if(plot_1D_std_error) {
          err_up_hist  ->SetFillColor(kRed);
          err_up_hist  ->SetFillStyle(1001);
          err_up_hist  ->SetLineWidth(1);
          err_up_hist  ->DrawClone("hist,same");
          err_down_hist->SetFillColor(10);
          err_down_hist->SetFillStyle(1001);
          err_down_hist->SetLineWidth(1);
          err_down_hist->DrawClone("hist,same");

          leg->AddEntry(mean_hist  , "Mean"     , "l");
          leg->AddEntry(err_up_hist, "#pm 1 #sigma", "f");
          leg->DrawClone();
        }

        mean_hist->DrawClone("hist,same");

        root_file_out->cd();
        if(!by_mode) {
          root_file_out  ->WriteObject(mean_hist    , TString::Format("%s_%s", subdata->Name().c_str(), "mean_PDF").Data());
          if(plot_1D_std_error) {
            root_file_out->WriteObject(err_up_hist  , TString::Format("%s_%s", subdata->Name().c_str(), "mean_PDF_err_up").Data());
            root_file_out->WriteObject(err_down_hist, TString::Format("%s_%s", subdata->Name().c_str(), "mean_PDF_err_down").Data());
          }
        }
        else {
          root_file_out  ->WriteObject(mean_hist    , TString::Format("%s_%s_%s", submodel->Name().c_str(), vutil::AsString(m).c_str(), "mean_PDF").Data());
          if(plot_1D_std_error) {
            root_file_out->WriteObject(err_up_hist  , TString::Format("%s_%s_%s", submodel->Name().c_str(), vutil::AsString(m).c_str(), "mean_PDF_err_up").Data());
            root_file_out->WriteObject(err_down_hist, TString::Format("%s_%s_%s", submodel->Name().c_str(), vutil::AsString(m).c_str(), "mean_PDF_err_down").Data());
          }
        }

        // Draw title
        TString title_str = TString::Format("%s", subdata->Name().c_str());
        if(by_mode) { title_str = TString::Format("%s_%s", submodel->Name().c_str(), vutil::AsString(m).c_str()); }
        TPaveText * title = new TPaveText(min_x * 1.05, max_y, max_x, max_y * 1.07 );
        title->SetShadowColor(0);
        title->SetFillColor(0);
        title->SetFillStyle(0);
        title->SetLineColor(0);
        title->AddText(title_str);
        title->SetTextSize(0.035);
        title->Draw("NB,same");

        // Save canvas to the output file
        std::string canv_save_name = save_name + "." + ((std::string)image_save_format);
        if(first_saved_plot && (cd->GetNSubData() > 1 || n_reac_modes_to_iterate > 1)) {
          canv_save_name += "(";
          first_saved_plot = false;
        }
        if(isubdata == cd->GetNSubData() -1 && im == last_mode &&
           (cd->GetNSubData() > 1 || n_reac_modes_to_iterate > 1)) {
          canv_save_name += ")";
        }
        canv->SaveAs(canv_save_name.c_str());

        delete mean_hist;
        delete err_up_hist;
        delete err_down_hist;
        delete leg;
        delete title;
        delete canv;
      } // im
    } // is 1D
    else {
      TString canv_name(subdata->Name().c_str());
      TCanvas* canv = new TCanvas(canv_name, canv_name, canv_size_x, canv_size_y);
      canv->cd();
      TH2D* mean_hist = (TH2D*) mean_hists.at(isubdata)->Clone();

      max_y = kinevar_y_max_2D;
      min_x = mean_hist->GetXaxis()->GetBinLowEdge(1);
      max_x = kinevar_x_max_2D;
      TVirtualPad *this_pad = canv->cd(1);
      this_pad->SetLeftMargin(0.16);
      this_pad->SetRightMargin(0.12);
      TH1F *pdf_frame = (TH1F *) this_pad->DrawFrame(0, 0, kinevar_x_max_2D, kinevar_y_max_2D);
      pdf_frame->GetXaxis()->SetTitle(title_x);
      pdf_frame->GetYaxis()->SetTitle(title_y);
      pdf_frame->GetZaxis()->SetTitle("Number of Events");
      pdf_frame->Draw();
      mean_hist->DrawClone("colz,same");

      root_file_out->cd();
      root_file_out->WriteObject(mean_hist, TString::Format("%s_%s", subdata->Name().c_str(), "mean_PDF").Data());

      // Draw title
      TString title_str = TString::Format("%s", subdata->Name().c_str());
      TPaveText * title = new TPaveText(min_x * 1.05, max_y, max_x, max_y * 1.07 );
      title->SetShadowColor(0);
      title->SetFillColor(0);
      title->SetFillStyle(0);
      title->SetLineColor(0);
      title->AddText(title_str);
      title->SetTextSize(0.035);
      title->Draw("NB,same");

      // Save canvas to the output file
      std::string canv_save_name = save_name + "." + ((std::string)image_save_format);
      if(first_saved_plot && cd->GetNSubData() > 1) {
        canv_save_name += "(";
        first_saved_plot = false;
      }
      if(isubdata == cd->GetNSubData() -1 && cd->GetNSubData() > 1) { canv_save_name += ")"; }
      canv->SaveAs(canv_save_name.c_str());

      delete mean_hist;
      delete title;
      delete canv;
    } // is 2D
  } // isubdata

  root_file_out->Close();
  delete root_file_out;

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