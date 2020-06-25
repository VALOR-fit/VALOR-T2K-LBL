//
// ------------------------------------------------------------------------------------------------
// Simple program to plot DUNE/ND PDFs and MC templates written out by write_nd_pdf.C
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

std::string kinevar = "Ereco"; // make this configurable later

void rj_plot_pdf(
  const char * config_file,
  const char * filename,
  const char * tag = "valor",
  double       per_range = -999., // Show events per what ereco (GeV). If negative, show per bin 
  TString plotting_format = "", // If empty or 'stacked': plot stacked interaction breakdown (valor default). If 'total': plot total signal + bg.
  const char * data_filename = "", // if we want to plot a dataset on top, this is that file
  const char * spectra_filename = "", // Save the stacked plot, total hist and data plot (if included) to a TFile with this name
  const char * data_tag = "valor",
  TString   option = "",
  const char * out_tag = "valor",
  const char * image_save_format = "",
  int canv_size_x = 700,
  int canv_size_y = 700,
  double enu_reco_max = 3.
)
{
  using namespace std;
  using namespace valor;
  namespace vutil = valor::util;

  const ExperimentSetup &es = *ExperimentSetup::Instance();

  const int n_pdf_slices = 10; //
  
  const int max_plots_per_canv = 12;

  int pdf_slice_color[n_pdf_slices] = 
  {
     kRed+1,   kOrange+7,   kYellow, kGreen+2,  kAzure+1, kViolet-1
  };
  
  bool write_image(image_save_format[0]);

  TFile * fin = new TFile(filename);
  assert(fin); assert(fin->IsOpen());
  
  bool plot_data(data_filename[0]);
  TFile * fdata = 0;
  if (plot_data) {
    fdata = new TFile(data_filename);
    assert(fdata); assert(fdata->IsOpen());
  }

  // If we have chosen to save the plots, get the filename
  bool save_spectra(spectra_filename[0]);
  TFile *temp_file = 0;
  if(save_spectra){
    std::cout << " Chosen to save plots " << std::endl;
    temp_file = new TFile(spectra_filename,"RECREATE");
    temp_file->Close();
    delete temp_file;
  }
  else{
    std::cout << " Chosen not to save plots " << std::endl;
    delete temp_file;
  }

  const char * binning = (option.Contains("OB")?"ob":"sb");

  // 0 == total, 1 = stacked
  int format = (plotting_format.Contains("total")?0:1);
  
  // Couple of vectors to keep all our objects
  vector<TCanvas*> all_canvases;
  vector<TH1*> all_hists;
  
  // First build a FitElements to tell us which samples we're plotting
  std::string config_path = vutil::CompleteConfigPath(config_file);
  CompoundFitElements elems;
  elems.Read(config_path.c_str(),"nomc nosplines");
  FitModelABC * modelabc = elems.FitModel();
  CompoundFitModel & model = *dynamic_cast<CompoundFitModel*>(modelabc);
  
  // A vector storing the list of submodels to plot on each canvas
  // Outer vector is per canvas, inner vector is plots on said canvas
  vector<vector<int> > canvas_content;
  canvas_content.push_back(vector<int>());
  
  conv::Detector_t last_det  = conv::kdUndefined;
  conv::BeamConf_t last_beam = conv::kbUndefined;
  
  // A buffer to write our event rates too
  std::ostringstream event_rates;
  
  // First loop over our selections to work out how many canvases we're drawing
  // We'll do a new canvas every time the detector changes, the beam current changes or we go above max_plots_per_canv
  for (unsigned int isubmodel = 0 ; isubmodel < model.GetNSubModels() ; isubmodel++) {
    FitModelABC * submodelabc = model.GetSubModel(isubmodel);
    GenericSingleSampleFitModel & submodel = *dynamic_cast<GenericSingleSampleFitModel*>(submodelabc);
    conv::Detector_t d = submodel.GetDetector();
    conv::BeamConf_t b = submodel.GetBeam();
    if (d == conv::kdUndefined) {continue;}
    if (b == conv::kbUndefined) {continue;}
    
    // If we don't match the last detector/beam start a new TCanvas
    if ((d!=last_det || b!= last_beam) && isubmodel) {
      canvas_content.push_back(vector<int>());
    }
    
    last_det  = d;
    last_beam = b;
    
    // If our canvas is full, start a new one
    if (canvas_content.back().size() == max_plots_per_canv) {
      canvas_content.push_back(vector<int>());
    }
    
    // Append the submodel number to the canvas' list
    canvas_content.back().push_back(isubmodel);
  }
  
  // Loop over the canvases we'll be drawing
  for (unsigned int icanvas = 0 ; icanvas < canvas_content.size() ; icanvas++) {
    // Calculate number of vertical and horizontal plots we need
    int plots_on_this_canv = canvas_content[icanvas].size();
    int nrows = 1;
    int ncols = 1;
    while (nrows * ncols < plots_on_this_canv) {
      if (canv_size_x > canv_size_y) {
        if (ncols < nrows) {ncols++;}
        else {nrows++;}
      }
      else {
        if (nrows < ncols) {nrows++;}
        else {ncols++;}
      }
    }
    // Get the first submodel so we can get beam/detector info
    unsigned int isubmodel = canvas_content[icanvas][0];
    FitModelABC * submodelabc = model.GetSubModel(isubmodel);
    GenericSingleSampleFitModel & submodel = *dynamic_cast<GenericSingleSampleFitModel*>(submodelabc);
    conv::Detector_t d = submodel.GetDetector();
    conv::BeamConf_t b = submodel.GetBeam();
    conv::Sample_t   s = submodel.GetSample();
    
    // Make a canvas for the mode plots
    TString c_name_contrib  = TString::Format("c_pdf_contrib_%s_%s_%u",
      vutil::AsString(d).c_str(),vutil::AsString(b).c_str(),icanvas);
    TString c_title_contrib = TString::Format("PDF (by mode) %s %s %u",
      vutil::AsString(d).c_str(),vutil::AsString(b).c_str(),icanvas);
    TCanvas * c_pdf_contrib = new TCanvas(c_name_contrib.Data(),c_title_contrib.Data(),20,20,canv_size_x,canv_size_y);
    c_pdf_contrib->Divide(nrows,ncols);
    all_canvases.push_back(c_pdf_contrib);
    
    // Actually read the PDFs
    for (unsigned int iplot = 0 ; iplot < canvas_content[icanvas].size() ; iplot++) {
      isubmodel = canvas_content[icanvas][iplot];
      submodelabc = model.GetSubModel(isubmodel);
      submodel = *dynamic_cast<GenericSingleSampleFitModel*>(submodelabc);
      d = submodel.GetDetector();
      b = submodel.GetBeam();
      s = submodel.GetSample();
      
      std::string model_name = submodel.Name();
      
      int isample = ExperimentSetup::Instance()->SampleIndex(b, s);
      int ibeam   = ExperimentSetup::Instance()->BeamConfIndex(b);
      
      { // Divided by mode
        THStack * contrib_stack = new THStack("stack","spectrum contributions");
        
        const size_t NReacModesCoarse = es.ReactionModesCoarse()[ibeam][isample].size();
        const size_t NReacModes       = es.ReactionModes()[ibeam][isample].size();
        std::cout << " Number of reaction modes : " << NReacModes << std::endl;

        std::vector<TH1D *> pdf_contrib             (NReacModes      ,nullptr);
        std::vector<TH1D *> pdf_contrib_width       (NReacModes      ,nullptr);
        std::vector<TH1D *> pdf_contrib_coarse      (NReacModesCoarse,nullptr);
        std::vector<TH1D *> pdf_contrib_coarse_width(NReacModesCoarse,nullptr);
        for(int imode = 0; imode < NReacModes; imode++) {
          conv::ReacMode_t m = es.ReactionModes()[ibeam][isample][imode];
          if(m == conv::kmUndefined) continue;
          cout << "Taking mode: " << vutil::AsString(m) <<" for "<<model_name<< endl;
          cout << "ibeam : " << ibeam << ", isample: " << isample << ", imode: " << imode << endl;
          TString contrib_name = TString::Format("%s_%s_%s_%s_%s_%s_%s",
            tag, "pdfcs",binning,"rf", model_name.c_str(),kinevar.c_str(), vutil::AsString(m).c_str());
          TObject * tmp_obj = fin->Get(contrib_name.Data());
          pdf_contrib[imode] = (TH1D*) tmp_obj;
          if(pdf_contrib[imode]) {
            conv::ReacModeCoarse_t mc = conv::kmcUndefined;
            bool found = false;
            for (int imcp = 0 ; imcp < conv::kMaxCoarseReacOptions ; imcp++) {
              mc = conv::kCoarseReacPriority[m][imcp];;
              if(mc == conv::kmcUndefined) continue; 
              for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
                if(mc == es.ReactionModesCoarse()[ibeam][isample][imodec]) {
                    found = true;
                    cout << "Mode: " << vutil::AsString(m) << " is mapped to plot type: " << vutil::AsString(mc) << endl;
                    break;
                }
              }
              if (found) break;
            }
            if(!found) {
              cout << "WTF! Did not find coarse reaction mode for mode: " << vutil::AsString(m) 
                   << " beam: " << vutil::AsString(b) << " / sample: " << vutil::AsString(s) << endl;
              cout << "Options were: "<<std::endl;
              for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
                cout<<"\t"<<vutil::AsString(es.ReactionModesCoarse()[ibeam][isample][imodec])<<std::endl;
              }
              exit(1);
            }                  
            for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
              if(mc != es.ReactionModesCoarse()[ibeam][isample][imodec]) continue;
              
              // Do we renormalise the bin-width before adding it?
              if (per_range > 0) {
                TH1D * width_norm = new TH1D(*pdf_contrib[imode]);
                for (int ibin = 1 ; ibin <= width_norm->GetXaxis()->GetNbins(); ibin++) {
                  double old_bin_content = width_norm->GetBinContent(ibin);
                  double bin_width = width_norm->GetXaxis()->GetBinWidth(ibin);
                  double new_bin_content = old_bin_content * per_range / bin_width;
                  width_norm->SetBinContent(ibin,new_bin_content);
                }
                if(pdf_contrib_coarse_width[imodec] == 0) {
                   pdf_contrib_coarse_width[imodec] = new TH1D(*width_norm);
                }
                else {
                   pdf_contrib_coarse_width[imodec]->Add(width_norm);
                }
                pdf_contrib_width[imode] = width_norm;
              }// per width
              if(pdf_contrib_coarse[imodec] == 0) {
                 pdf_contrib_coarse[imodec] = new TH1D(*pdf_contrib[imode]);
              }
              else {
                 pdf_contrib_coarse[imodec]->Add(pdf_contrib[imode]);
              }
            }//imodec
          }//!null
          else {
            cerr << "Null [" << contrib_name << "]" << endl;
            cerr << "Giving up." << endl;
            cerr << "This is usually caused by running this script with a different config file to write_pdf"<<std::endl;
            exit(-1);
          }
        }//imode
        
        // Scan for max in desired Erange
        double max_n = 0.;
        for(int ebin = 1; ebin <= pdf_contrib[0]->GetNbinsX(); ebin++) {
           if(pdf_contrib[0]->GetBinCenter(ebin) < enu_reco_max ) {
             double sum = 0;
             for(int imode = 0; imode < NReacModes; imode++) {
                TH1D * hist = (per_range > 0 ? pdf_contrib_width[imode] : pdf_contrib[imode]);
                if (hist) {
                  double val = hist->GetBinContent(ebin);
                  sum += val;
                }
             }
             max_n = TMath::Max(max_n, sum);           
           }//<emax
        }//ebin
        
        TH1D * data_hist = 0;
        float data_max = -99999.;
        // If we need data, get that histogram now
        if (plot_data) {
          TString data_name = TString::Format("%s_%s_%s_%s_%s_%s",
              data_tag, "pdfts",binning,"rf", model_name.c_str(),kinevar.c_str());
          TObject * tmp = fdata->Get(data_name.Data());
          std::cout<<data_name<<"\t"<<tmp<<std::endl;
          TH1D * htmp = (TH1D*)tmp;
          data_hist = new TH1D(*htmp);

          // and if necessary normalise to bin width
          if (per_range > 0) {
            for (int ibin = 1 ; ibin <= data_hist->GetXaxis()->GetNbins(); ibin++) {
              double old_bin_content = data_hist->GetBinContent(ibin);
              double bin_width = data_hist->GetXaxis()->GetBinWidth(ibin);
              double new_bin_content = old_bin_content * per_range / bin_width;
              data_hist->SetBinContent(ibin,new_bin_content);
            }
          }

          data_hist->SetLineStyle(0);
          data_hist->SetMarkerColor(kRed);
          data_hist->SetMarkerStyle(kFullTriangleDown);
          data_max = data_hist->GetMaximum();
        }
        // If totals rather than by mode, get that 
        TH1D * total_hist = 0;
        float total_max = -99999.;
        if (format == 0) {
          TString total_name = TString::Format("%s_%s_%s_%s_%s_%s",
              tag, "pdfts",binning,"rf", model_name.c_str(),kinevar.c_str());
          TObject * tmp = fin->Get(total_name.Data());
          std::cout<<total_name<<"\t"<<tmp<<std::endl;
          TH1D * htmp = (TH1D*)tmp;
          total_hist = new TH1D(*htmp);

          // and if necessary normalise to bin width
          if (per_range > 0) {
            for (int ibin = 1 ; ibin <= total_hist->GetXaxis()->GetNbins(); ibin++) {
              double old_bin_content = total_hist->GetBinContent(ibin);
              double bin_width = total_hist->GetXaxis()->GetBinWidth(ibin);
              double new_bin_content = old_bin_content * per_range / bin_width;
              total_hist->SetBinContent(ibin,new_bin_content);
            }
          }

          total_hist->SetLineStyle(2);
          total_hist->SetLineWidth(2);
          total_hist->SetLineColor(8);
          total_hist->GetMaximum();
        }

        // Start drawing
        TLegend * contrib_legend = new TLegend(0.7,0.3,0.9,0.8);  
        contrib_legend->SetFillColor(0);
        contrib_legend->SetLineColor(0);
        TH1D * cc_hist    = nullptr;
        TH1D * nc_hist    = nullptr;
      
        float mode_max = 0.;
        for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
          conv::ReacModeCoarse_t mc = es.ReactionModesCoarse()[ibeam][isample][imodec];
          TString mode_str = vutil::AsString(mc);
          //cout<<"Drawing: "<<vutil::AsString(b)<<" "<<vutil::AsString(s)<<"\t"<<imodec<<"\t"<<vutil::AsString(mc)<<"\t"<<pdf_contrib_coarse[imodec]<<std::endl;
          if(mc == conv::kmcUndefined) continue;
          TH1D * relevant_hist = per_range > 0. ? pdf_contrib_coarse_width[imodec] : pdf_contrib_coarse[imodec];
          if(per_range > .0){
            if(!cc_hist && !(mode_str.Contains("NC")))
              cc_hist = (TH1D*)pdf_contrib_coarse_width[imodec]->Clone("cc_hist");
            else if(!nc_hist && (mode_str.Contains("NC")))
              nc_hist = (TH1D*)pdf_contrib_coarse_width[imodec]->Clone("nc_hist");
          }
          else {
            if(!cc_hist && !(mode_str.Contains("NC")))
              cc_hist = (TH1D*)pdf_contrib_coarse[imodec]->Clone("cc_hist");
            else if(!nc_hist && (mode_str.Contains("NC")))
              nc_hist = (TH1D*)pdf_contrib_coarse[imodec]->Clone("nc_hist");
          }
          relevant_hist->SetLineWidth(1);
          relevant_hist->SetLineColor(es.PlotStyle()[ibeam][isample][imodec][0]);
          relevant_hist->SetFillColor(es.PlotStyle()[ibeam][isample][imodec][1]);
          relevant_hist->SetFillStyle(es.PlotStyle()[ibeam][isample][imodec][2]);
          contrib_stack->Add(relevant_hist);
          contrib_legend->AddEntry(relevant_hist,Form("%s",vutil::AsStringTex(mc,true).c_str()), "LF");

          if(cc_hist && !(mode_str.Contains("NC")))
            cc_hist->Add(relevant_hist);
          else if(nc_hist && (mode_str.Contains("NC")))
            nc_hist->Add(relevant_hist);

          mode_max += relevant_hist->GetMaximum();
        }//imodec
        
        float max_y = std::max(mode_max, std::max(data_max, total_max)) * 1.1;

        if(format == 0) total_hist->GetYaxis()->SetRangeUser(0,max_y);
        
        TVirtualPad * this_pad = c_pdf_contrib->cd(iplot+1); // Pad for this plot
        this_pad->SetLeftMargin(0.2);
        TH1F * contrib_frame = (TH1F*)this_pad->DrawFrame(0,0,enu_reco_max,1.2*mode_max);
        contrib_frame->Draw();
        TString contrib_title_str = TString::Format("%s / %s / %s",
          vutil::AsString(d).c_str(), vutil::AsString(b).c_str(), vutil::AsStringTeX(s,true).c_str());
        TLatex * contrib_title = new TLatex(0.3,1.07*mode_max,contrib_title_str.Data());
        contrib_title->Draw(); 
        if(vutil::IsCC(s)) {  
           contrib_frame->GetXaxis()->SetTitle("Reco. E_{#nu} (GeV)");
        } else {
           contrib_frame->GetXaxis()->SetTitle("Visible Energy (GeV)");
        }
        if (per_range > 0.) {
          TString title = TString::Format("Events / %lg GeV",per_range);
          contrib_frame->GetYaxis()->SetTitle(title.Data());
        }

        else {
          contrib_frame->GetYaxis()->SetTitle("Events /bin");
        }
        if(format == 1){
          contrib_stack->Draw("same");
          if (plot_data) {data_hist->Draw("e same");}
          contrib_legend->Draw("same");
        }
        else{
          total_hist->Draw("same");
          if (plot_data) {data_hist->Draw("e same");}
        }
       
        double tot_event_rate = 0.;
        std::ostringstream sample_event_rates;
        // Write the event rates to the output stream
        for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
          conv::ReacModeCoarse_t mc = es.ReactionModesCoarse()[ibeam][isample][imodec];
          if (mc == conv::kmcUndefined) {continue;}
          TH1D * tmphist = (TH1D*)pdf_contrib_coarse[imodec]->Clone();
          double rate = tmphist->Integral();
          sample_event_rates<<"    "<<vutil::AsString(mc)<<": "<<rate<<std::endl;
          tot_event_rate += rate;
          delete tmphist;
        }
        event_rates<<vutil::AsString(b)<<" "<<vutil::AsString(d)<<" "<<vutil::AsString(s)<<" (total) : "<<tot_event_rate<<std::endl;
        event_rates<<sample_event_rates.str()<<std::endl; 
        
        // If we have chosen to save the spectra, save both the stacked distribution and the total
        // If we have also included the data/oscillated spectra, save that too
        if(save_spectra){
          TFile spectra_plots(spectra_filename,"UPDATE");
          contrib_stack->SetName((vutil::AsString(d)+"_stacked").c_str());
          contrib_legend->SetName((vutil::AsString(d)+"_legend").c_str());
          total_hist->SetName((vutil::AsString(d)+"_integrated").c_str());
          cc_hist->SetName((vutil::AsString(d)+"_cc").c_str());
          nc_hist->SetName((vutil::AsString(d)+"_nc").c_str());
          contrib_stack->Write();
          contrib_legend->Write();
          total_hist->Write();
          cc_hist->Write();
          nc_hist->Write();
          if(plot_data){
            data_hist->SetName((vutil::AsString(d)+"_oscillated").c_str());
            data_hist->Write();
          }
          spectra_plots.Write();
          spectra_plots.Close();
        } // save spectra
      } // imodec
    } //iplot
    
    // Update the canvases
    //~ c_pdf_tot    ->Update();
    //~ c_pdf_slice  ->Update();
    c_pdf_contrib->Update();
    
    // And write to files if needed
    if (write_image) {
      //~ TString pdf_tot_filename = TString::Format("%s_spectrum_%s_%s_%u_2D.%s",
          //~ out_tag, vutil::AsString(d).c_str(), vutil::AsString(b).c_str(),icanvas,image_save_format);
      //~ c_pdf_tot->SaveAs(pdf_tot_filename.Data());
      
      //~ TString slice_filename = TString::Format("%s_spectrum_%s_%s_%u_y.%s",
          //~ out_tag, vutil::AsString(d).c_str(), vutil::AsString(b).c_str(),icanvas,image_save_format);
      //~ c_pdf_slice->Print(slice_filename.Data());
      
      TString contrib_filename = TString::Format("%s_spectrum_%s_%s_%u_modes.%s",
        out_tag, vutil::AsString(d).c_str(), vutil::AsString(b).c_str(),icanvas, image_save_format);
      c_pdf_contrib->Print(contrib_filename.Data());
    } //write_image
      
  } // icanvas
  
  std::cout<<event_rates.str();
    
  return;
}
