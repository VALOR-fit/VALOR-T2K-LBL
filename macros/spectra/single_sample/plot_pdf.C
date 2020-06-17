//
// ------------------------------------------------------------------------------------------------
// Simple program to plot DUNE/ND PDFs and MC templates written out by write_nd_pdf.C
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

bool t2k               = false;   // make this configurable later
bool one_plot_per_canv = false;   // make this configurable later
double title_text_size = 0.05;
std::string caveat_text_str = "T2K Run1-9 Preliminary";

void plot_pdf(
  const char * config_file,
  const char * filename,
  const char * tag               = "valor",
  const char * kinevar_c         = "Ereco", // Plot slices in this kinematic variable
  double       per_range         = -999.,   // Show events per what ereco (GeV). If negative, show per bin
  const char * data_filename     = "",      // if we want to plot a dataset on top, this is that file
  const char * data_tag          = "valor",
  TString      option            = "",
  const char * out_tag           = "valor",
  const char * image_save_format = "",
  int          canv_size_x       = 700,
  int          canv_size_y       = 700,
  double       kinevar_x_max_1D  = 3.,      // Max value of the kinematic variable drawn on the x-axis (1D plots)
  double       kinevar_x_max_2D  = 1.25,    // Max value of the kinematic variable drawn on the x-axis (2D plots)
  double       kinevar_y_max_2D  = 180.,    // Max value of the kinematic variable drawn on the y-axis (2D plots)
  bool         draw_caveat       = false,   // Draws a text box in the top right containing the contents of caveat_text_str
  const char * title_x           = "#nu Reconstructed Energy (GeV)",  // Used for 2D plots only
  const char * title_y           = "#theta (degrees)"                 // Used for 2D plots only
)
{
  using namespace std;
  using namespace valor;
  namespace vutil = valor::util;

  if(!vutil::plotting::SetupCustomPalette("valor_inverted_dark_body_radiator")) {
    std::cout << "WARNING: the custom color palette could not be set up. Reverting to default" << std::endl;
  }

  const int n_pdf_slices = 10; //
  
  const int max_plots_per_canv = 12;

  int pdf_slice_color[n_pdf_slices] = 
  {
     kRed+1,   kOrange+7,   kYellow, kGreen+2,  kAzure+1, kViolet-1
  };

  std::string kinevar(kinevar_c);
  std::cout << "Plotting " << kinevar << " slices" << std::endl;
  
  bool write_image(image_save_format[0]);

  TFile * fin = new TFile(filename);
  assert(fin); assert(fin->IsOpen());
  
  bool plot_data(data_filename[0]);
  TFile * fdata = 0;
  if (plot_data) {
    fdata = new TFile(data_filename);
    assert(fdata);
  }

  const char * binning = (option.Contains("OB")?"ob":"sb");
  
  
  // Couple of vectors to keep all our objects
  vector<TCanvas*> all_canvases;
  vector<TH1*> all_hists;
  
  // First build a FitElements to tell us which samples we're plotting
  std::string config_path = vutil::CompleteConfigPath(config_file);
  CompoundFitElements elems;
  elems.Read(config_path.c_str(),"nomc nosplines");
  FitModelABC * modelabc = elems.FitModel();
  FitDataABC  * dataabc  = elems.FitData();
  CompoundFitModel & model = *dynamic_cast<CompoundFitModel*>(modelabc);
  CompoundFitData  & data  = *dynamic_cast<CompoundFitData*> (dataabc);

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
    if ((((d!=last_det || b!= last_beam) && isubmodel) || one_plot_per_canv) && canvas_content.back().size() != 0) {
      canvas_content.push_back(vector<int>());
    }
    
    last_det  = d;
    last_beam = b;
    
    // If our canvas is full, start a new one
    if ((canvas_content.back().size() == max_plots_per_canv) && !one_plot_per_canv) {
      canvas_content.push_back(vector<int>());
    }
    
    // Append the submodel number to the canvas' list
    canvas_content.back().push_back(isubmodel);
  }
  
  // Loop over the canvases we'll be drawing
  for (unsigned int icanvas = 0 ; icanvas < canvas_content.size() ; icanvas++) {
    bool draw_pdft_2D = true;
    bool draw_pdfc    = true;

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
    FitDataABC  * subdataabc  = data .GetSubData (isubmodel);
    GenericSingleSampleFitModel & submodel = *dynamic_cast<GenericSingleSampleFitModel*>(submodelabc);
    GenericSingleSampleFitData  & subdata  = *dynamic_cast<GenericSingleSampleFitData*> (subdataabc);
    conv::Detector_t d = submodel.GetDetector();
    conv::BeamConf_t b = submodel.GetBeam();
    conv::Sample_t   s = submodel.GetSample();
    const ModeList* ml = submodel.GetModeList();
    
    // make a canvas for the 2d plots
    TString c_name_pdf_tot  = TString::Format("c_pdf_tot_%s_%s_%u",
      vutil::AsString(d).c_str(),vutil::AsString(b).c_str(),icanvas);
    TString c_title_pdf_tot = TString::Format("Total PDF %s %s %u",
    vutil::AsString(d).c_str(),vutil::AsString(b).c_str(),icanvas);
    TCanvas * c_pdf_tot = new TCanvas(c_name_pdf_tot.Data(),c_title_pdf_tot.Data(),20,20,canv_size_x,canv_size_y);
    c_pdf_tot->Divide(nrows,ncols);
    all_canvases.push_back(c_pdf_tot);
    
    //~ // Make a canvas for the 1D y-sliced plots
    //~ TString c_name_pdf_slice  = TString::Format("c_pdf_slice_%s_%s_%u",
      //~ vutil::AsString(d).c_str(),vutil::AsString(b).c_str(),icanvas);
    //~ TString c_title_pdf_slice = TString::Format("PDF (by y_reco) %s %s %u",
      //~ vutil::AsString(d).c_str(),vutil::AsString(b).c_str(),icanvas);
    //~ TCanvas * c_pdf_slice = new TCanvas(c_name_pdf_slice.Data(),c_title_pdf_slice.Data(),20,20,canv_size_x,canv_size_y);
    //~ c_pdf_slice->Divide(nrows,ncols);
    //~ all_canvases.push_back(c_pdf_slice);
    
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
      subdataabc  = data .GetSubData (isubmodel);
      submodel = *dynamic_cast<GenericSingleSampleFitModel*>(submodelabc);
      subdata  = *dynamic_cast<GenericSingleSampleFitData*> (subdataabc);
      d = submodel.GetDetector();
      b = submodel.GetBeam();
      s = submodel.GetSample();

      if(t2k) {
        if(kinevar == "Ereco") {
          if(vutil::SampleIsMuLike(s))  { kinevar_x_max_1D = 3.00; }
          else                          { kinevar_x_max_1D = 1.25; }
        }
        else if(kinevar == "Thetareco") { kinevar_x_max_1D = 180.; }
      }
      
      std::string model_name = submodel.Name();

      const SingleSampleKineBinning* sskb = submodel.GetBinning();
      conv::KineSpace_t ks = conv::kspUndefined;
      if(binning == std::string("ob")) {
        const VALORBinning & reco_opt = *(sskb->RecoOptBinning());
        ks = vutil::KineSpaceFromKineVars(reco_opt.KineVars());
      }
      else {
        const VALORBinning & reco_smp = *(sskb->RecoSmpBinning());
        ks = vutil::KineSpaceFromKineVars(reco_smp.KineVars());
      }

      if(vutil::NDim(ks) == 2) { // 2D PDF
        TString pdf_tot_name = TString::Format("%s_%s_%s_%s_%s_%s", tag, "pdft", binning, "rf", model_name.c_str(),
                                               vutil::AsString(ks).c_str());

        TObject *tmp_obj = fin->Get(pdf_tot_name.Data());
        if(!tmp_obj) {
          std::cout << "WARNING: Hist named [" << pdf_tot_name.Data() << "] not found. Skipping 2D PDFT draw."
                    << std::endl;
          draw_pdft_2D = false;
        } else {
          TH2D *pdf_tot = new TH2D(*((TH2D *) tmp_obj));
          all_hists.push_back(pdf_tot);

          TGraph * data_graph_inner = 0;
          TGraph * data_graph_outer = 0;
          // If we need data, get that TGraph now (TGraph because data should be unbinned on 2D plots)
          if (plot_data) {
            TString data_name = TString::Format("%s_%s_%s_%s_%s",
                                                data_tag, "data", vutil::AsString(d).c_str(),
                                                vutil::AsString(b).c_str(), vutil::AsString(s).c_str());
            TObject * tmp = fdata->Get(data_name.Data());
            std::cout << "Reading 2D data TGraph [" << data_name << "]\t pointer = " << tmp << std::endl;
            TGraph * gtmp = (TGraph*)tmp;
            data_graph_outer = new TGraph(*gtmp);
            data_graph_inner = new TGraph(*gtmp);

            // Make sure that the PDF kinevars match those of the TGraph
            conv::KineVar_t kv_pdf_x = conv::kvUndefined;
            conv::KineVar_t kv_pdf_y = conv::kvUndefined;
            if(binning == std::string("ob")) {
              const VALORBinning & reco_opt = *(sskb->RecoOptBinning());
              kv_pdf_x = reco_opt.KineVars().at(0);
              kv_pdf_y = reco_opt.KineVars().at(1);
            }
            else {
              const VALORBinning & reco_smp = *(sskb->RecoSmpBinning());
              kv_pdf_x = reco_smp.KineVars().at(0);
              kv_pdf_y = reco_smp.KineVars().at(1);
            }
            if(data_graph_outer->GetXaxis()->GetTitle() != vutil::AsString(kv_pdf_x)) {
              std::cerr << "ERROR: Data TGraph [" << data_name << "] x-axis kinevar ["
                        << data_graph_outer->GetXaxis()->GetTitle() << "] does not match that of the PDF ["
                        << vutil::AsString(kv_pdf_x) << "]" << std::endl;
              vutil::Die(-1);
            }
            if(data_graph_outer->GetYaxis()->GetTitle() != vutil::AsString(kv_pdf_y)) {
              std::cerr << "ERROR: Data TGraph [" << data_name << "] y-axis kinevar ["
                        << data_graph_outer->GetYaxis()->GetTitle() << "] does not match that of the PDF ["
                        << vutil::AsString(kv_pdf_y) << "]" << std::endl;
              vutil::Die(-1);
            }

            // Draw the markers as a circle surrounded by a ring of a different colour.
            // This helps the markers to be more visible to colourblind people
            data_graph_outer->SetMarkerColor(kAzure - 3);
            data_graph_outer->SetMarkerStyle(kFullCircle);
            data_graph_outer->SetMarkerSize(2);
            data_graph_inner->SetMarkerColor(10);
            data_graph_inner->SetMarkerStyle(kFullCircle);
            data_graph_inner->SetMarkerSize(1);
          } // plot data?

          // Draw the 2D plots
          TVirtualPad *this_pad = c_pdf_tot->cd(iplot + 1); // Pad for this plot
          this_pad->SetLeftMargin(0.12);  //0.16
          this_pad->SetRightMargin(0.16); //0.12
          TH1F *pdf_frame = (TH1F *) this_pad->DrawFrame(0, 0, kinevar_x_max_2D, kinevar_y_max_2D);
          pdf_frame->GetXaxis()->SetTitle(title_x);
          pdf_frame->GetYaxis()->SetTitle(title_y);
          pdf_frame->Draw();

          TLatex * pdf_title = nullptr;
          if(!one_plot_per_canv) {
            TString pdf_title_str = TString::Format("%s / %s / %s",
                                                    vutil::AsString(d).c_str(), vutil::AsString(b).c_str(),
                                                    vutil::AsStringTeX(s,true).c_str());
            double pdf_title_xlow = 0.0;
            TLatex * pdf_title = new TLatex(pdf_title_xlow, 1.01*kinevar_y_max_2D, pdf_title_str.Data());
            pdf_title->SetTextSize(title_text_size);
            pdf_title->Draw();
          }

          TLatex * caveat = nullptr;
          if(draw_caveat) {
            double tx = pdf_tot->GetXaxis()->GetBinLowEdge(pdf_tot->GetXaxis()->GetLast());
            caveat = new TLatex(0.5*tx, 1.01*kinevar_y_max_2D, caveat_text_str.c_str());
            caveat->SetTextSize(title_text_size);
            caveat->Draw();
          }

          pdf_tot->GetXaxis()->SetTitle(title_x);
          pdf_tot->GetYaxis()->SetTitle(title_y);
          pdf_tot->GetZaxis()->SetTitle("Number of Events");
          pdf_tot->GetZaxis()->SetLabelSize(pdf_frame->GetXaxis()->GetLabelSize()); // The x and y axis labels and titles are set somewhere when running write_pdf, so just copy the styles here
          pdf_tot->GetZaxis()->SetLabelOffset(pdf_frame->GetXaxis()->GetLabelOffset());
          pdf_tot->GetZaxis()->SetLabelFont(pdf_frame->GetXaxis()->GetTitleFont());
          pdf_tot->GetZaxis()->SetTitleSize(pdf_frame->GetXaxis()->GetTitleSize());
          pdf_tot->GetZaxis()->SetTitleOffset(0.85);
          pdf_tot->GetZaxis()->SetTitleFont(pdf_frame->GetXaxis()->GetTitleFont());
          pdf_tot->Draw("colzsame");

          // Draw the data markers
          if (plot_data) {
            data_graph_outer->Draw("P0,same");
            data_graph_inner->Draw("P0,same");
          }
        }
      }
      

      //~ { // Sliced in Y
        //~ THStack * slice_stack = new THStack("stack","spectrum contributions");
        //~ const int n_y_reco_bins_sb = vutil::NumYrecoBinsSb(s);
        //~ const double * y_reco_bin_edges_sb = vutil::YrecoBinEdgesSb(s);
        //~ vector<TH1*> pdf_slice;
        //~ for(int iy = 0; iy < n_y_reco_bins_sb; iy++) {
          //~ TString slice_name = TString::Format("%s_%s_%s_%s_%s_%s_y%d",
            //~ tag, "pdfts",binning,"rf", vutil::AsString(b).c_str(), vutil::AsString(s).c_str(), iy);
            //~ TObject * tmp_obj = fin->Get( slice_name );
            //~ pdf_slice.push_back((TH1D*)tmp_obj);
        //~ }//iy
        
        //~ double max_n = 0.;
        //~ // Scan for max in desired Erange
        //~ for(int ebin = 1; ebin <= pdf_slice[0]->GetNbinsX(); ebin++) {
           //~ if(pdf_slice[0]->GetBinCenter(ebin) < kinevar_x_max_1D ) {
             //~ double sum = 0;
             //~ for(int iy = 0; iy < n_y_reco_bins_sb; iy++) {
                //~ double val = pdf_slice[iy]->GetBinContent(ebin);
                //~ sum += val;
             //~ }
             //~ max_n = TMath::Max(max_n, sum);           
            //~ }//<emax
        //~ }//ebin
        
        //~ // Make slice legend
        //~ TLegend * slice_legend = new TLegend(0.70,0.40,0.88,0.80);
        //~ slice_legend->SetFillColor(0);
        //~ for(int iy = 0; iy < n_y_reco_bins_sb; iy++) {
          //~ pdf_slice[iy]->SetFillColor(pdf_slice_color[iy]); 
          //~ pdf_slice[iy]->SetLineWidth(1); 
          //~ slice_stack->Add(pdf_slice[iy]);
          //~ slice_legend->AddEntry(pdf_slice[iy],Form("y = [%.2f,%.2f]", y_reco_bin_edges_sb[iy], y_reco_bin_edges_sb[iy+1]),"LF");
        //~ }
        
        //~ // Draw the slice plots
        //~ TVirtualPad * this_pad = c_pdf_slice->cd(iplot+1); // Pad for this plot
        //~ TH1F * slice_frame = (TH1F*)this_pad->DrawFrame(0,0,kinevar_x_max_1D,1.2*max_n);
        //~ slice_frame->Draw();
        //~ TString slice_title_str = TString::Format("%s / %s", vutil::AsString(b).c_str(), vutil::AsStringTeX(s,true).c_str());
        //~ TLatex * slice_title = new TLatex(1,1.07*max_n,slice_title_str.Data());
        //~ slice_title->Draw(); 
        //~ if(vutil::IsCC(s)) {
          //~ slice_frame->GetXaxis()->SetTitle("Reco. E_{#nu} (GeV)");
        //~ } else {
          //~ slice_frame->GetXaxis()->SetTitle("Visible Energy (GeV)");
        //~ }
        //~ slice_frame->GetYaxis()->SetTitle("Events /bin");
        //~ slice_stack->Draw("same");
        //~ slice_legend->Draw();
      //~ } // end slice in y
      
      
      { // Divided by mode
        THStack * contrib_stack = new THStack("stack","spectrum contributions");

        const std::vector<conv::ReacMode_t>       modes        = ml->Modes();
        const std::vector<conv::ReacModeCoarse_t> coarse_modes = ml->CoarseModes();
        const size_t NReacModesCoarse = coarse_modes.size();
        const size_t NReacModes       = modes.size();

        std::vector<TH1D *> pdf_contrib             (NReacModes      ,nullptr);
        std::vector<TH1D *> pdf_contrib_width       (NReacModes      ,nullptr);
        std::vector<TH1D *> pdf_contrib_coarse      (NReacModesCoarse,nullptr);
        std::vector<TH1D *> pdf_contrib_coarse_width(NReacModesCoarse,nullptr);
        for(int imode = 0; imode < NReacModes; imode++) {
          conv::ReacMode_t m = modes.at(imode);
          cout << "Taking mode: " << vutil::AsString(m) <<" for "<<model_name<< endl;
          cout << "beam : " << vutil::AsString(b) << ", sample: " << vutil::AsString(s) << ", imode: " << imode << endl;
          TString contrib_name = TString::Format("%s_%s_%s_%s_%s_%s_%s",
            tag, "pdfcs",binning,"rf", model_name.c_str(),kinevar.c_str(), vutil::AsString(m).c_str());
          cout << "beam : " << vutil::AsString(b) << ", sample: " << vutil::AsString(s) << ", imode: " << imode << endl;
          TObject * tmp_obj = fin->Get(contrib_name.Data());
          pdf_contrib[imode] = (TH1D*) tmp_obj;
          if(pdf_contrib[imode]) {
            if(t2k && kinevar == "Ereco") {
              // For T2K, the binning used in the FHC1Rmu and RHC1Rmu plots are coarser than the fit binning.
              // As problems are caused when filling the MC templates when the simple binning has fewer bins than the opt
              // binning, rebin the relevant hists here instead of changing the binning.
              if(s == conv::ks1Rmu && b == conv::kbJPARC_FHC) {
                pdf_contrib[imode]->Rebin(2);
              }
              else if(s == conv::ks1Rmu && b == conv::kbJPARC_RHC) {
                pdf_contrib[imode]->Rebin(3);
              }
            }
            conv::ReacModeCoarse_t mc = conv::kmcUndefined;
            bool found = false;
            for (int imcp = 0 ; imcp < conv::kMaxCoarseReacOptions ; imcp++) {
              mc = conv::kCoarseReacPriority[m][imcp];;
              if(mc == conv::kmcUndefined) continue; 
              for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
                if(mc == coarse_modes.at(imodec)) {
                    found = true;
                    cout << "Mode: " << vutil::AsString(m) << " is mapped to plot type: " << vutil::AsString(mc) << endl;
                    break;
                }
              }
              if (found) break;
            }
            if(!found) {
              cout << "WTF! Did not find coarse reaction mode for mode: " << vutil::AsString(m) 
                   << " / beam: " << vutil::AsString(b) << " / sample: " << vutil::AsString(s) << endl;
              cout << "Options were: "<<std::endl;
              for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
                cout<<"\t"<<vutil::AsString(coarse_modes.at(imodec))<<std::endl;
              }
              exit(1);
            }                  
            for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
              if(mc != coarse_modes.at(imodec)) continue;
              
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
            cerr << "WARNING: Null [" << contrib_name << "]. Skipping PDFC draw." << endl;
            draw_pdfc = false;
            break;
          }
        }//imode

        if (draw_pdfc) {
          // Scan for max in desired Erange
          double max_n = 0.;
          for(int ebin = 1; ebin <= pdf_contrib[0]->GetNbinsX(); ebin++) {
             if(pdf_contrib[0]->GetBinCenter(ebin) < kinevar_x_max_1D ) {
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
          // If we need data, get that histogram now (as 1D data should be binned)
          if (plot_data) {
            // First try reading data from a stored TGraph
            TGraph * data_graph = 0;
            TString data_name = TString::Format("%s_%s_%s_%s_%s",
                                                data_tag, "data", vutil::AsString(d).c_str(),
                                                vutil::AsString(b).c_str(), vutil::AsString(s).c_str());
            TObject * tmp = fdata->Get(data_name.Data());
            std::cout << "Reading 2D data TGraph [" << data_name << "]\t pointer = " << tmp << std::endl;
            TGraph * gtmp = (TGraph*)tmp;
            data_graph = new TGraph(*gtmp);

            if(data_graph) {
              // Copy the binning of the hists in the stack
              data_hist = new TH1D(*pdf_contrib_coarse[0]);
              data_hist->Reset();

              // Fill the data hist from the TGraph
              for(int ipoint = 0; ipoint < data_graph->GetN(); ++ipoint) {
                double kine_val = -999, x = -999, y = -999;
                data_graph->GetPoint(ipoint, x, y);
                if     (std::string(data_graph->GetXaxis()->GetTitle()) == kinevar) { kine_val = x; }
                else if(std::string(data_graph->GetYaxis()->GetTitle()) == kinevar) { kine_val = y; }
                else {
                  std::cerr << "ERROR: Could not read 1D data from input Tgraph ["
                            << data_name.Data() << "]" << std::endl;
                  vutil::Die(-1);
                }
                data_hist->Fill(kine_val);
              }

              if(!data_hist) {
                std::cerr << "ERROR: Could not fill 1D data hist from the TGraph ["
                          << data_name.Data() << "]" << std::endl;
                vutil::Die(-1);
              }
            } // found data graph
            else {
              // Look for a 1D data hist instead
              data_name = TString::Format("%s_%s_%s_%s_%s_%s",
                                          data_tag, "data", binning, "rf_slice", subdata.Name().c_str(), kinevar.c_str());
              TObject * tmp = fdata->Get(data_name.Data());
              std::cout << "Reading data array [" << data_name << "]\t pointer = " << tmp << std::endl;
              TH1D * htmp = (TH1D*)tmp;
              data_hist = new TH1D(*htmp);

              if(!data_hist) {
                std::cerr << "ERROR: Could not read 1D data hist directly ["
                          << data_name.Data() << "]" << std::endl;
                vutil::Die(-1);
              }
            } // did not find data graph

            // and if necessary normalise to bin width
            if (per_range > 0) {
              for (int ibin = 1 ; ibin <= data_hist->GetXaxis()->GetNbins(); ibin++) {
                double old_bin_content = data_hist->GetBinContent(ibin);
                double bin_width = data_hist->GetXaxis()->GetBinWidth(ibin);
                double new_bin_content = old_bin_content * per_range / bin_width;
                data_hist->SetBinContent(ibin,new_bin_content);
              }
            }

            if(t2k) {
              data_hist->SetLineStyle(0);
              data_hist->SetLineColor(kBlack);
              data_hist->SetMarkerColor(kBlack);
              data_hist->SetMarkerStyle(kFullCircle);
              data_hist->SetMarkerSize(1.5);
              data_hist->SetBarWidth(0.4);
            }
            else {
              data_hist->SetLineStyle(0);
              data_hist->SetMarkerColor(kRed);
              data_hist->SetMarkerStyle(kFullTriangleDown);
            }

            // Scan the data + error for max in the desired range
            for(int ebin = 1; ebin <= data_hist->GetNbinsX(); ebin++) {
              if(data_hist->GetBinCenter(ebin) < kinevar_x_max_1D ) {
                double val = data_hist->GetBinContent(ebin) +  data_hist->GetBinError(ebin);
                max_n = TMath::Max(max_n, val);
              }//<emax
            }//ebin
          } // plot data?

          // Start drawing
          TLegend * contrib_legend = new TLegend(0.7,0.55,0.9,0.95);
          contrib_legend->SetFillColor(0);
          for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
              conv::ReacModeCoarse_t mc = coarse_modes.at(imodec);
              //~ cout<<"Drawing: "<<vutil::AsString(b)<<" "<<vutil::AsString(s)<<"\t"<<imodec<<"\t"<<vutil::AsString(mc)<<"\t"<<pdf_contrib_coarse[imodec]<<std::endl;
              TH1D * relevant_hist = per_range > 0. ? pdf_contrib_coarse_width[imodec] : pdf_contrib_coarse[imodec];
              relevant_hist->SetLineWidth(1);
              relevant_hist->SetLineColor(ml->GetCoarseModeLineCol(mc));
              relevant_hist->SetFillColor(ml->GetCoarseModeFillCol(mc));
              relevant_hist->SetFillStyle(ml->GetCoarseModeFillSty(mc));
              contrib_stack->Add(relevant_hist);
              contrib_legend->AddEntry(relevant_hist,Form("%s",vutil::AsStringTex(mc,true).c_str()), "F");
          }//imodec

          double max_y_1D = 1.1 * max_n;
          TVirtualPad * this_pad = c_pdf_contrib->cd(iplot+1); // Pad for this plot
          this_pad->SetLeftMargin(0.12);
          TH1F * contrib_frame = (TH1F*)this_pad->DrawFrame(0,0,kinevar_x_max_1D,max_y_1D);
          contrib_frame->Draw();

          TString contrib_title_str = TString::Format("%s / %s / %s",
            vutil::AsString(d).c_str(), vutil::AsString(b).c_str(), vutil::AsStringTeX(s,true).c_str());
          double contrib_title_xlow = 0.0;
          TLatex * contrib_title = new TLatex(contrib_title_xlow, 1.01*max_y_1D, contrib_title_str.Data());
          contrib_title->SetTextSize(title_text_size);
          if(!one_plot_per_canv) { contrib_title->Draw(); }

          TLatex * caveat = nullptr;
          if(draw_caveat) {
            double tx = contrib_stack->GetXaxis()->GetBinLowEdge(contrib_stack->GetXaxis()->GetLast());
            caveat = new TLatex(0.5*tx, 1.01*max_y_1D, caveat_text_str.c_str());
            caveat->SetTextSize(title_text_size);
            caveat->Draw();
          }

          // Set the x-axis title based on the kinevar we've projected to
          if(kinevar == "Ereco") {
            if(vutil::IsCC(s)) {
              if(t2k) { contrib_frame->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)"); }
              else    { contrib_frame->GetXaxis()->SetTitle("Reco. E_{#nu} (GeV)"); }
            } else {
               contrib_frame->GetXaxis()->SetTitle("Visible Energy (GeV)");
            }
          }
          else if(kinevar == "Thetareco") {
            contrib_frame->GetXaxis()->SetTitle("#theta (degrees)");
          }
          else {
            std::cout << "WARNING: did not recogise kinevar name [" << kinevar
                      << "], so could not set the x-axis title" << std::endl;
          }

          if (per_range > 0.) {
            TString title = TString::Format("Events /%lg GeV",per_range);
            contrib_frame->GetYaxis()->SetTitle(title.Data());
          }
          else {
            if(t2k) { contrib_frame->GetYaxis()->SetTitle("Number of Events"); }
            else    { contrib_frame->GetYaxis()->SetTitle("Events /bin"); }
          }
          contrib_stack->Draw("same");
          if (plot_data) {data_hist->Draw("e same");}
          gPad->RedrawAxis();  // Draw the axes again, so that the axis tick marks are not obscured by the stack
          contrib_legend->SetTextFont(contrib_frame->GetXaxis()->GetTitleFont());
          contrib_legend->SetTextSize(0.75 * contrib_frame->GetXaxis()->GetTitleSize());
          contrib_legend->Draw();

          double tot_event_rate = 0.;
          std::ostringstream sample_event_rates;
          // Write the event rates to the output stream
          for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
            conv::ReacModeCoarse_t mc = coarse_modes.at(imodec);
            TH1D * tmphist = (TH1D*)pdf_contrib_coarse[imodec]->Clone();
            double rate = tmphist->Integral();
            sample_event_rates<<"    "<<vutil::AsString(mc)<<": "<<rate<<std::endl;
            tot_event_rate += rate;
            delete tmphist;
          }
          event_rates<<vutil::AsString(b)<<" "<<vutil::AsString(d)<<" "<<vutil::AsString(s)<<" (total) : "<<tot_event_rate<<std::endl;
          event_rates<<sample_event_rates.str()<<std::endl;

          // And I guess compare them to the total spectra as a sanity check.
          TString total_name = TString::Format("%s_%s_%s_%s_%s",
            tag, "pdft",binning,"vf", model_name.c_str());
          VALORArray * tot_arr = dynamic_cast<VALORArray*>(fin->Get(total_name.Data()));
          std::cerr<<total_name.Data()<<std::endl;
          assert(tot_arr);
          double tot_event_rate_from_varray = tot_arr->Integral();
          if (std::fabs(tot_event_rate_from_varray / tot_event_rate - 1) > conv::kASmallNum) {
            event_rates<<"\033[1;31m"<<"Warning - stack total is "<<tot_event_rate
                       <<" while array total is "<<tot_event_rate_from_varray
                       <<" - something has gone wrong (unless the hist has been rebinned"
                       << " [like when using the \'t2k\' option], in which case a slight discrepancy is expected)!\n\n"
                       <<"\033[0m"<<std::endl;
          }
        } // draw_pdfc
        
      } // imodec
      
    
    } //iplot
    
    // Update the canvases
    //~ c_pdf_tot    ->Update();
    //~ c_pdf_slice  ->Update();
    c_pdf_contrib->Update();
    
    // And write to files if needed
    if (write_image) {
      if(draw_pdft_2D) {
        TString pdf_tot_filename("");
        if(one_plot_per_canv) {
          pdf_tot_filename = TString::Format("%s_spectrum_%s_2D.%s",
                                             out_tag, submodel.Name().c_str(), image_save_format);
        }
        else {
          pdf_tot_filename = TString::Format("%s_spectrum_%s_%s_%u_2D.%s",
                                             out_tag, vutil::AsString(d).c_str(), vutil::AsString(b).c_str(),
                                             icanvas, image_save_format);
        }

        c_pdf_tot->SaveAs(pdf_tot_filename.Data());
      }
      
      //~ TString slice_filename = TString::Format("%s_spectrum_%s_%s_%u_y.%s",
          //~ out_tag, vutil::AsString(d).c_str(), vutil::AsString(b).c_str(),icanvas,image_save_format);
      //~ c_pdf_slice->Print(slice_filename.Data());

      if(draw_pdfc) {
        TString contrib_filename("");
        if(one_plot_per_canv) {
          contrib_filename = TString::Format("%s_spectrum_%s_%s_modes.%s",
                                             out_tag, submodel.Name().c_str(), kinevar.c_str(), image_save_format);
        }
        else {
          contrib_filename = TString::Format("%s_spectrum_%s_%s_%s_%u_modes.%s",
                                             out_tag, vutil::AsString(d).c_str(), vutil::AsString(b).c_str(),
                                             kinevar.c_str(), icanvas, image_save_format);
        }
        c_pdf_contrib->Print(contrib_filename.Data());
      }
    } //write_image
      
  } // icanvas
  
  std::cout<<event_rates.str();

  std::cout << "Done! (a following segfault is expected...)" << std::endl;
    
  return;
}
