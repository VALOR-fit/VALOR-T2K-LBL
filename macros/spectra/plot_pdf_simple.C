//
//______________________________________________________________________________________________________________________
// Simple program to plot T2K PDFs written out by valor-write_pdf.
// Adapted for T2K from plot_pdf.C
//
// Francis Bench
// 2020-09-24
//______________________________________________________________________________________________________________________
//

// TODO: Option of drawing all plots of the same sample with the same y or z axis

using namespace valor;
namespace vutil = valor::util;

//______________________________________________________________________________________________________________________
// Hard-coded global variables
double      _title_text_size        = 0.04;
std::string _caveat_text_str        = "T2K Run1-10 Preliminary";
bool        _mark_empty_bins_2d     = true;
double      _erec_max_mulike        = 3.;
double      _erec_max_elike         = 1.25;
double      _theta_max_elike        = 180.;
std::string _title_x                = "#nu Reconstructed Energy (GeV)";
std::string _title_y                = "#theta (degrees)";
double      _error_bar_significance = 0.6827;
std::string _data_label             = "T2K Data";

//______________________________________________________________________________________________________________________
// Global variables assigned from command-line input
bool        _first_draw      = true;
TFile *     _fin             = nullptr;
std::string _output_filename = "";
int         _canv_size_x     = 0;
int         _canv_size_y     = 0;
bool        _draw_caveat     = false;
bool        _rebin_mulike    = false;
std::string _binning         = "";
std::string _data_filename   = "";
std::string _data_tag        = "";
std::string _tag             = "";
bool        _thesis_format   = false;

//______________________________________________________________________________________________________________________
// Prototypes
void                 plot_PDF_totals              (CompoundFitElements* elems);
void                 plot_PDF_contribs_2D         (CompoundFitElements* elems);
void                 plot_PDF_projections         (CompoundFitElements* elems);
void                 draw_single_PDF_2D           (GenericSingleSampleFitModel& submodel, conv::ReacMode_t m = conv::kmUndefined);
void                 draw_single_PDF_contrib_1D   (GenericSingleSampleFitModel& submodel, GenericSingleSampleFitData& subdata, conv::KineVar_t kv);
conv::KineSpace_t    GetKineSpace                 (GenericSingleSampleFitModel& submodel);
TGraphAsymmErrors*   GetDataGraph1D               (GenericSingleSampleFitData& subdata, TH1D* binning_hist, conv::KineVar_t kv, conv::Detector_t d, conv::BeamConf_t b, conv::Sample_t s);
std::vector<TGraph*> GetDataGraphs2D              (GenericSingleSampleFitModel& submodel, conv::Detector_t d, conv::BeamConf_t b, conv::Sample_t s);
void                 save_to_file                 (TCanvas* canv, std::string bookmark_text = "", bool close_file = false);
TLatex*              make_caveat                  (TH1* hist_to_draw, double y_max);
TLatex*              make_title                   (TString title_text, double ymax);
std::string          make_title_text              (conv::BeamConf_t b, conv::Sample_t s, conv::ReacMode_t m = conv::kmUndefined, bool no_tex = false, bool no_spaces = false);
std::string          get_beam_current_name_simple (conv::BeamConf_t b);                    // TODO: move to BeamConf.h
std::string          get_sample_name_simple       (conv::Sample_t s, bool no_tex = false); // TODO: move to Sample.h
void                 rebin                        (TH1D* hist, conv::BeamConf_t b, conv::Sample_t s, conv::KineVar_t kv);
double               get_max_in_xrange_by_mode    (double xmax, std::vector<TH1D*> hists_by_mode, const int NReacModes);
double               get_max_in_xrange            (double xmax, TH1D* hist, bool with_error);
double               get_max_in_xrange_asymm_err  (double xmax, TGraphAsymmErrors* graph);
void                 insert_section_break         (std::string title_text);
double               get_kinevar_max              (conv::Sample_t s, conv::KineVar_t kv);
TGraph*              get_empty_bin_markers        (TH2D* pdf);

//______________________________________________________________________________________________________________________
// Main
void plot_pdf_simple(
  const char * config_file,
  const char * output_filename    = "spectra_all_sb_t2k_run10_postfit_asimova",
  const char * input_filename     = "write_all_sb_t2k_run10_postfit_aimova.root",
  const char * in_tag             = "valor", // The tag used when writing the hists in valor-write_pdf
  const char * data_filename      = "",      // if we want to plot a dataset on top, this is that file
  const char * data_tag           = "valor",
  TString      binning            = "",      // 'SB' (simple binning - recommended) or 'OB' (optimised binning)
  int          canv_size_x        = 700,     // In pixels
  int          canv_size_y        = 700,
  bool         draw_caveat        = true,    // Draws a text box in the top right containing the contents of _caveat_text_str
  bool         rebin_mulike       = true,    // rebin the mu-like sample 1D plots with coarser binning, for cosmetic reasons
  bool         mark_empty_bins_2d = false,   // Draws a cross in each empty 2D PDF bin
  bool         thesis_format      = false
)
{
  if(!vutil::plotting::SetupCustomPalette("valor_inverted_dark_body_radiator")) {
    std::cout << "WARNING: the custom color palette could not be set up. Reverting to default" << std::endl;
  }

  // Assign global variables from command line input
  _fin = new TFile(input_filename, "READ");
  assert(_fin); assert(_fin->IsOpen());
  _output_filename    = std::string(output_filename);
  _canv_size_x        = canv_size_x;
  _canv_size_y        = canv_size_y;
  _draw_caveat        = draw_caveat;
  _binning            = std::string(binning.Contains("OB")?"ob":"sb");
  _data_filename      = data_filename;
  _data_tag           = data_tag;
  _tag                = in_tag;
  _mark_empty_bins_2d = mark_empty_bins_2d;
  _rebin_mulike       = rebin_mulike;
  _thesis_format      = thesis_format;

  // First build a FitElements to tell us which samples we're plotting
  std::string config_path = vutil::CompleteConfigPath(config_file);
  CompoundFitElements* elems = new CompoundFitElements();
  elems->Read(config_path.c_str(),"nomc nosplines");

  // Start plotting
  insert_section_break("Total Spectra");
  plot_PDF_totals(elems);

  insert_section_break("#splitline{1D Projections}{of the 2D Spectra}");
  plot_PDF_projections(elems);

  insert_section_break("#splitline{Contributions (by-mode)}{to the 2D Spectra}");
  plot_PDF_contribs_2D(elems);

  save_to_file(nullptr, "", true);  // Close the PDF
  _fin->Close();
  return;
}

//______________________________________________________________________________________________________________________
void plot_PDF_totals(CompoundFitElements* elems)
{
  FitModelABC * modelabc = elems->FitModel();
  FitDataABC  * dataabc  = elems->FitData();
  CompoundFitModel & model = *dynamic_cast<CompoundFitModel*>(modelabc);
  CompoundFitData  & data  = *dynamic_cast<CompoundFitData*> (dataabc);
  
  for (unsigned int isubmodel = 0 ; isubmodel < model.GetNSubModels() ; isubmodel++) {
    // Get sub-model/data
    FitModelABC *submodelabc = model.GetSubModel(isubmodel);
    FitDataABC  *subdataabc  = data .GetSubData (isubmodel);
    GenericSingleSampleFitModel &submodel = *dynamic_cast<GenericSingleSampleFitModel *>(submodelabc);
    GenericSingleSampleFitData  &subdata  = *dynamic_cast<GenericSingleSampleFitData  *>(subdataabc);

    conv::KineSpace_t ks = GetKineSpace(submodel);
    std::vector<conv::KineVar_t> kvs = vutil::KineVar(ks);

    if(vutil::NDim(ks) == 1) {
      draw_single_PDF_contrib_1D(submodel, subdata, kvs.at(0));
    }
    else if(vutil::NDim(ks) == 2) { // 2D PDF
      draw_single_PDF_2D(submodel);
    }
    else {
      std::cerr << "ERROR: Unexpected number of kinematic dimensions to plot!" << std::endl;
      vutil::Die(1);
    }
  } // isubmodel
}

//______________________________________________________________________________________________________________________
void plot_PDF_projections(CompoundFitElements* elems)
{
  FitModelABC * modelabc = elems->FitModel();
  FitDataABC  * dataabc  = elems->FitData();
  CompoundFitModel & model = *dynamic_cast<CompoundFitModel*>(modelabc);
  CompoundFitData  & data  = *dynamic_cast<CompoundFitData*> (dataabc);

  for (unsigned int isubmodel = 0 ; isubmodel < model.GetNSubModels() ; isubmodel++) {
    // Get sub-model/data
    FitModelABC *submodelabc = model.GetSubModel(isubmodel);
    FitDataABC  *subdataabc  = data .GetSubData (isubmodel);
    GenericSingleSampleFitModel &submodel = *dynamic_cast<GenericSingleSampleFitModel *>(submodelabc);
    GenericSingleSampleFitData  &subdata  = *dynamic_cast<GenericSingleSampleFitData  *>(subdataabc);

    conv::KineSpace_t ks = GetKineSpace(submodel);
    std::vector<conv::KineVar_t> kvs = vutil::KineVar(ks);

    if(vutil::NDim(ks) == 2) { // 2D PDF
      for(auto kv : kvs) {
        draw_single_PDF_contrib_1D(submodel, subdata, kv);
      }
    }
  } // isubmodel
}

//______________________________________________________________________________________________________________________
void plot_PDF_contribs_2D(CompoundFitElements* elems)
{
  FitModelABC * modelabc = elems->FitModel();
  CompoundFitModel & model = *dynamic_cast<CompoundFitModel*>(modelabc);
  
  for (unsigned int isubmodel = 0 ; isubmodel < model.GetNSubModels() ; isubmodel++) {
    // Get sub-model/data
    FitModelABC *submodelabc = model.GetSubModel(isubmodel);
    GenericSingleSampleFitModel &submodel = *dynamic_cast<GenericSingleSampleFitModel *>(submodelabc);

    conv::KineSpace_t ks = GetKineSpace(submodel);
    if(vutil::NDim(ks) != 2) {  continue; }
    
    const ModeList* ml = submodel.GetModeList();
    const std::vector<conv::ReacMode_t> modes = ml->Modes();
    const size_t NReacModes = modes.size();
    for(int imode = 0; imode < NReacModes; imode++) {
      conv::ReacMode_t m = modes.at(imode);
      draw_single_PDF_2D(submodel, m);
    }
  } // isubmodel
}

//______________________________________________________________________________________________________________________
conv::KineSpace_t GetKineSpace(GenericSingleSampleFitModel& submodel)
{
  const SingleSampleKineBinning* sskb = submodel.GetBinning();
  conv::KineSpace_t ks = conv::kspUndefined;
  if(_binning == "ob") {
    const VALORBinning & reco_opt = *(sskb->RecoOptBinning());
    ks = vutil::KineSpaceFromKineVars(reco_opt.KineVars());
  }
  else {
    const VALORBinning & reco_smp = *(sskb->RecoSmpBinning());
    ks = vutil::KineSpaceFromKineVars(reco_smp.KineVars());
  }
  return ks;
}

//______________________________________________________________________________________________________________________
std::vector<TGraph*> GetDataGraphs2D(GenericSingleSampleFitModel& submodel, conv::Detector_t d, conv::BeamConf_t b, conv::Sample_t s)
{
  std::vector<TGraph*> data_graphs(2, nullptr);

  bool plot_data(_data_filename[0] && !_mark_empty_bins_2d);
  TFile * fdata = 0;
  if (plot_data) {
    fdata = new TFile(_data_filename.c_str());
    if(!fdata->IsOpen()) {
      std::cout << "ERROR: Could not open data file \'" << _data_filename << "\'" << std::endl;
      vutil::Die(-1);
    }
  }
  else {
    return data_graphs;
  }

  TString data_name = TString::Format("%s_%s_%s_%s_%s",
                                      _data_tag.c_str(), "data", vutil::AsString(d).c_str(),
                                      vutil::AsString(b).c_str(), vutil::AsString(s).c_str());
  TObject * tmp = fdata->Get(data_name.Data());
  TGraph * gtmp = (TGraph*)tmp;
  TGraph * data_graph_outer = new TGraph(*gtmp);
  TGraph * data_graph_inner = new TGraph(*gtmp);

  // Make sure that the PDF kinevars match those of the TGraph
  conv::KineVar_t kv_pdf_x = conv::kvUndefined;
  conv::KineVar_t kv_pdf_y = conv::kvUndefined;
  const SingleSampleKineBinning* sskb = submodel.GetBinning();
  if(_binning == "ob") {
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

  data_graphs[0] = data_graph_inner;
  data_graphs[1] = data_graph_outer;
  return data_graphs;
}

//______________________________________________________________________________________________________________________
TGraphAsymmErrors* GetDataGraph1D(GenericSingleSampleFitData& subdata, TH1D* binning_hist, conv::KineVar_t kv, conv::Detector_t d, conv::BeamConf_t b, conv::Sample_t s)
{
  TGraphAsymmErrors * data_graph_draw = nullptr;
  TH1D * data_hist = nullptr;

  bool plot_data(_data_filename[0]);
  TFile * fdata = 0;
  if (plot_data) {
    fdata = new TFile(_data_filename.c_str());
    if(!fdata->IsOpen()) {
      std::cout << "ERROR: Could not open data file \'" << _data_filename << "\'" << std::endl;
      vutil::Die(-1);
    }
  }
  else {
    return data_graph_draw;
  }

  // First try reading data from a stored TGraph
  data_graph_draw = new	TGraphAsymmErrors();
  TGraph * data_graph = 0;
  TString data_name = TString::Format("%s_%s_%s_%s_%s",
                                      _data_tag.c_str(), "data", vutil::AsString(d).c_str(),
                                      vutil::AsString(b).c_str(), vutil::AsString(s).c_str());
  TObject * tmp = fdata->Get(data_name.Data());
  std::cout << "Reading 2D data TGraph [" << data_name << "]\t pointer = " << tmp << std::endl;
  if(tmp) {
    TGraph * gtmp = (TGraph*)tmp;
    data_graph = new TGraph(*gtmp);
  }

  if(data_graph) {
    data_hist = new TH1D(*binning_hist);
    data_hist->Reset();

    // Fill the data hist from the TGraph
    for(int ipoint = 0; ipoint < data_graph->GetN(); ++ipoint) {
      double kine_val = -999, x = -999, y = -999;
      data_graph->GetPoint(ipoint, x, y);
      if     (std::string(data_graph->GetXaxis()->GetTitle()) == vutil::AsString(kv)) { kine_val = x; }
      else if(std::string(data_graph->GetYaxis()->GetTitle()) == vutil::AsString(kv)) { kine_val = y; }
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
                                _data_tag.c_str(), "data", _binning.c_str(), "rf_slice",
                                subdata.Name().c_str(), vutil::AsString(kv).c_str());
    TObject * tmp = fdata->Get(data_name.Data());
    std::cout << "Reading data array [" << data_name << "]\t pointer = " << tmp << std::endl;
    if(tmp) {
      TH1D * htmp = (TH1D*)tmp;
      data_hist = new TH1D(*htmp);
    }

    if(!data_hist) {
      std::cerr << "ERROR: Could not read 1D data hist directly ["
                << data_name.Data() << "]" << std::endl;
      vutil::Die(-1);
    }
  } // did not find data graph

  // Convert the data hist to a TGraphAsymmErrors so we can set the data error bars appropriately.
  // Setting sqrt(N) as the statistical error on a data point is wrong, as the data has no stat error.
  // Instead calculate an error band showing how far away the prediction can be while still keeping
  // the data point within its Poisson stat errors.
  const double alpha = 1 - _error_bar_significance;
  for (int ibin = 1 ; ibin <= data_hist->GetXaxis()->GetNbins(); ibin++) {
    double center  = data_hist->GetXaxis()->GetBinCenter(ibin);
    double content = data_hist->GetBinContent(ibin);
    double width   = data_hist->GetXaxis()->GetBinWidth(ibin);

    double err_ylow  = (content==0) ? content : content - (ROOT::Math::gamma_quantile(alpha/2, content, 1.));
    double err_yhigh = ROOT::Math::gamma_quantile_c(alpha/2, content+1, 1) - content;
    double err_x     = width*0.5;

    data_graph_draw->SetPoint(ibin-1, center, content);
    data_graph_draw->SetPointError(ibin-1, err_x, err_x, err_ylow, err_yhigh);
  }
  
  data_graph_draw->SetLineStyle(0);
  data_graph_draw->SetLineColor(kBlack);
  data_graph_draw->SetMarkerColor(kBlack);
  data_graph_draw->SetMarkerStyle(kFullCircle);
  data_graph_draw->SetMarkerSize(1.5);
  //data_graph_draw->SetBarWidth(0.4);
  
  return data_graph_draw;
}

//______________________________________________________________________________________________________________________
// Draw a 2D PDF as a coloured hist
// If no mode is supplied, it draws the total PDF.
void draw_single_PDF_2D(GenericSingleSampleFitModel& submodel, conv::ReacMode_t m)
{
  conv::Detector_t d = submodel.GetDetector();
  conv::BeamConf_t b = submodel.GetBeam();
  conv::Sample_t   s = submodel.GetSample();
  std::string model_name = submodel.Name();
  conv::KineSpace_t ks = GetKineSpace(submodel);
  if(d == conv::kdUndefined) { return; }
  if(b == conv::kbUndefined) { return; }
  if(s == conv::ksUndefined) { return; }

  // Setup canvas
  TString canv_name  = TString::Format("canv_%s_%s_%s"  , vutil::AsString(d).c_str(), vutil::AsString(b).c_str(), vutil::AsString(s).c_str());
  TString canv_title = TString::Format("PDF 2D %s %s %s", vutil::AsString(d).c_str(), vutil::AsString(b).c_str(), vutil::AsString(s).c_str());
  std::string format = "pdft";
  if(m != conv::kmUndefined) {
    canv_name  = TString::Format("%s_%s", canv_name .Data(), vutil::AsString(m).c_str());
    canv_title = TString::Format("%s_%s", canv_title.Data(), vutil::AsString(m).c_str());
    format = "pdfc";
  }

  TCanvas * canv = new TCanvas(canv_name.Data(), canv_title.Data(), 20, 20, _canv_size_x, _canv_size_y);

  // Get the object to draw
  TString read_obj_name = TString::Format("%s_%s_%s_%s_%s_%s", _tag.c_str(), format.c_str(), _binning.c_str(), "rf",
                                          model_name.c_str(), vutil::AsString(ks).c_str());
  if(m != conv::kmUndefined) {
    read_obj_name = TString::Format("%s_%s", read_obj_name.Data(), vutil::AsString(m).c_str());
  }
  TObject *tmp_obj = _fin->Get(read_obj_name.Data());
  if(!tmp_obj) {
    std::cout << "WARNING: Hist named [" << read_obj_name.Data() << "] not found. Skipping 2D PDFT draw."
              << std::endl;
    return;
  }
  TH2D *pdf = new TH2D(*((TH2D *) tmp_obj));

  // Draw the 2D plots
  canv->SetLeftMargin(0.12);  //0.16
  canv->SetRightMargin(0.16); //0.12
  TH1F *pdf_frame = (TH1F *) canv->DrawFrame(0, 0, get_kinevar_max(s, conv::kvEreco), get_kinevar_max(s, conv::kvThetareco));
  pdf_frame->GetXaxis()->SetTitle(_title_x.c_str());
  pdf_frame->GetYaxis()->SetTitle(_title_y.c_str());
  pdf_frame->GetYaxis()->SetTitleOffset(0.95);
  pdf_frame->GetZaxis()->SetTitleOffset(0.90);

  std::string title_text = make_title_text(b, s, m);
  TLatex * title = make_title(title_text.c_str(), get_kinevar_max(s, conv::kvThetareco));
  title->Draw();

  TLatex* caveat = nullptr;
  if(_draw_caveat) {
    caveat = make_caveat((TH1*) pdf, get_kinevar_max(s, conv::kvThetareco));
    caveat->Draw();
  }

  pdf->GetZaxis()->SetTitle("Number of Events");
  pdf->GetZaxis()->SetLabelSize  (pdf_frame->GetXaxis()->GetLabelSize()); // The x and y axis labels and titles are set somewhere when running write_pdf, so just copy the styles here
  pdf->GetZaxis()->SetLabelOffset(pdf_frame->GetXaxis()->GetLabelOffset());
  pdf->GetZaxis()->SetLabelFont  (pdf_frame->GetXaxis()->GetTitleFont());
  pdf->GetZaxis()->SetTitleSize  (pdf_frame->GetXaxis()->GetTitleSize());
  pdf->GetZaxis()->SetTitleFont  (pdf_frame->GetXaxis()->GetTitleFont());
  pdf->GetZaxis()->SetTitleOffset(0.90);
  pdf->Draw("colzsame");

  // Draw the data markers as a TGraph (because data should be unbinned on 2D plots)
  std::vector<TGraph*> data_graphs = GetDataGraphs2D(submodel, d, b, s);
  assert(data_graphs.size() == 2);
  if (data_graphs[0] && data_graphs[1]) {
    data_graphs[1]->Draw("P0,same");
    data_graphs[0]->Draw("P0,same");
  }

  // Optionally mark empty bins in the 2D spectra
  TGraph* empty_markers = get_empty_bin_markers(pdf);
  if(empty_markers) {
    empty_markers->Draw("P0,same");
  }

  // Redraw the axes over the 2D hist, as the hist usually obscures the tick marks.
  pdf_frame->Draw("AXIS,SAME");

  std::string descriptor = "Total";
  if(m != conv::kmUndefined) {
    descriptor = "Contrib. 2D";
  }
  std::string bookmark_text = descriptor + " " + make_title_text(b, s, m, true, false);
  save_to_file(canv, bookmark_text);

  // Clean up objects we made with 'new'
  if(title)          { delete title;          }
  if(caveat)         { delete caveat;         }
  if(pdf)            { delete pdf;            }
  if(data_graphs[0]) { delete data_graphs[0]; }
  if(data_graphs[1]) { delete data_graphs[1]; }
  if(canv)           { delete canv;           }
}

//______________________________________________________________________________________________________________________
// Draw a 1D PDF as a stack of contributing coarse reaction modes
void draw_single_PDF_contrib_1D(GenericSingleSampleFitModel& submodel, GenericSingleSampleFitData& subdata, conv::KineVar_t kv)
{
  conv::Detector_t d = submodel.GetDetector();
  conv::BeamConf_t b = submodel.GetBeam();
  conv::Sample_t   s = submodel.GetSample();
  if(d == conv::kdUndefined) { return; }
  if(b == conv::kbUndefined) { return; }
  if(s == conv::ksUndefined) { return; }
  std::string model_name = submodel.Name();

  // Setup canvas
  TString canv_name  = TString::Format("c_pdf_contrib_1d_%s_%s_%s", vutil::AsString(d).c_str(),
                                            vutil::AsString(b).c_str(), vutil::AsString(s).c_str());
  TString canv_title = TString::Format("PDF Contributions %s %s %s", vutil::AsString(d).c_str(),
                                            vutil::AsString(b).c_str(), vutil::AsString(s).c_str());
  TCanvas * canv = new TCanvas(canv_name.Data(), canv_title.Data(), 20, 20, _canv_size_x, _canv_size_y);

  // Get the list and number of modes
  const ModeList* ml = submodel.GetModeList();
  const std::vector<conv::ReacMode_t>       modes        = ml->Modes();
  const std::vector<conv::ReacModeCoarse_t> coarse_modes = ml->CoarseModes();
  const size_t NReacModesCoarse = coarse_modes.size();
  const size_t NReacModes       = modes.size();

  std::vector<TH1D *> pdf_contrib        (NReacModes      , nullptr);
  std::vector<TH1D *> pdf_contrib_coarse (NReacModesCoarse, nullptr);
  bool draw_pdfc = true;
  for(int imode = 0; imode < NReacModes; imode++) {
    conv::ReacMode_t m = modes.at(imode);
    TString read_obj_name = TString::Format("%s_%s_%s_%s_%s_%s_%s", _tag.c_str(), "pdfcs", _binning.c_str(), "rf",
                                            model_name.c_str(), vutil::AsString(kv).c_str(), vutil::AsString(m).c_str());

    TObject * tmp_obj = _fin->Get(read_obj_name.Data());
    pdf_contrib[imode] = (TH1D*) tmp_obj;
    if(!pdf_contrib[imode]) {
      std::cerr << "WARNING: Null [" << read_obj_name << "]. Skipping PDFC draw for beam : "
           << vutil::AsString(b) << ", sample: " << vutil::AsString(s) << std::endl;
      draw_pdfc = false;
      break;
    }

    rebin(pdf_contrib[imode], b, s, kv);

    // Decide which coarse reaction mode this mode lies in, and add the corresponding mode PDF to the coarse mode stack
    conv::ReacModeCoarse_t mc = conv::kmcUndefined;
    bool found = false;
    for (int imcp = 0 ; imcp < conv::kMaxCoarseReacOptions ; imcp++) {
      mc = conv::kCoarseReacPriority[m][imcp];;
      if(mc == conv::kmcUndefined) continue;
      for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
        if(mc == coarse_modes.at(imodec)) {
          found = true;
          break;
        }
      }
      if (found) break;
    }
    if(!found) {
      std::cout << "WTF! Did not find coarse reaction mode for mode: " << vutil::AsString(m)
                << " / beam: " << vutil::AsString(b) << " / sample: " << vutil::AsString(s) << std::endl;
      std::cout << "Options were: "<<std::endl;
      for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
        std::cout<<"\t"<<vutil::AsString(coarse_modes.at(imodec))<<std::endl;
      }
      exit(1);
    }

    for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
      if(mc != coarse_modes.at(imodec)) continue;
      if(pdf_contrib_coarse[imodec] == 0) { pdf_contrib_coarse[imodec] = new TH1D(*pdf_contrib[imode]); }
      else                                { pdf_contrib_coarse[imodec]       ->Add(pdf_contrib[imode]); }
    }//imodec

  }//imode

  if(!draw_pdfc) { return; }

  double max_n = get_max_in_xrange_by_mode(get_kinevar_max(s, kv), pdf_contrib, NReacModes);

  TGraphAsymmErrors* data_graph = GetDataGraph1D(subdata, pdf_contrib[0], kv, d, b, s);
  if(data_graph) {
    // Scan the data + error for max in the desired range
    double data_max_n = get_max_in_xrange_asymm_err(get_kinevar_max(s, kv), data_graph);
    max_n = TMath::Max(max_n, data_max_n);
  }

  // Start drawing
  THStack * contrib_stack = new THStack("stack","spectrum contributions");
  TLegend * contrib_legend = new TLegend(0.75,0.55,0.95,0.95);
  contrib_legend->SetFillColor(0);
  for(int imodec = 0; imodec < NReacModesCoarse; imodec++) {
    conv::ReacModeCoarse_t mc = coarse_modes.at(imodec);
    pdf_contrib_coarse[imodec]->SetLineWidth(1);
    pdf_contrib_coarse[imodec]->SetLineColor(ml->GetCoarseModeLineCol(mc));
    pdf_contrib_coarse[imodec]->SetFillColor(ml->GetCoarseModeFillCol(mc));
    pdf_contrib_coarse[imodec]->SetFillStyle(ml->GetCoarseModeFillSty(mc));
    contrib_stack->Add(pdf_contrib_coarse[imodec]);

    std::string coarse_name = vutil::AsStringTex(mc, true);
    if(_thesis_format) {
      if     (mc == conv::kmcOscNuECC)    { coarse_name  =  "#nu_{#mu} #rightarrow #nu_{e} CC"; }
      else if(mc == conv::kmcOscNuEBarCC) { coarse_name  =  "#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e} CC"; }
    }

    contrib_legend->AddEntry(pdf_contrib_coarse[imodec],Form("%s",coarse_name.c_str()), "F");
  }//imodec
  if(data_graph) { contrib_legend->AddEntry(data_graph, _data_label.c_str(), "LP"); }

  double max_y_1D = 1.1 * max_n;
  canv->SetLeftMargin(0.12);
  TH1F * contrib_frame = (TH1F*)canv->DrawFrame(0, 0, get_kinevar_max(s, kv), max_y_1D);
  contrib_frame->Draw();

  std::string title_text = make_title_text(b, s);
  TLatex * title = make_title(title_text.c_str(), max_y_1D);
  title->Draw();

  TLatex* caveat = nullptr;
  if(_draw_caveat) {
    caveat = make_caveat((TH1*) contrib_stack->GetHistogram(), max_y_1D);
    caveat->Draw();
  }

  // Set the x-axis title based on the kinevar we've projected to
  if(vutil::AsString(kv) == "Ereco") {
    contrib_frame->GetXaxis()->SetTitle("#nu Reconstructed Energy (GeV)");
  }
  else if(vutil::AsString(kv) == "Thetareco") {
    contrib_frame->GetXaxis()->SetTitle("#theta (degrees)");
  }
  else {
    std::cout << "WARNING: did not recogise kinevar name [" << vutil::AsString(kv)
              << "], so could not set the x-axis title" << std::endl;
  }
  contrib_frame->GetYaxis()->SetTitle("Number of Events");

  contrib_stack->Draw("same");
  if (data_graph) {data_graph->Draw("P,()");}
  gPad->RedrawAxis();  // Draw the axes again, so that the axis tick marks are not obscured by the stack
  contrib_legend->SetTextFont(contrib_frame->GetXaxis()->GetTitleFont());
  contrib_legend->SetTextSize(0.75 * contrib_frame->GetXaxis()->GetTitleSize());
  contrib_legend->Draw();


  std::string bookmark_text = "Contrib. 1D " + make_title_text(b, s, conv::kmUndefined, true, false);
  save_to_file(canv, bookmark_text);

  // Clean up objects we made with 'new'
  if(title)          { delete title;          }
  if(caveat)         { delete caveat;         }
  if(contrib_stack)  { delete contrib_stack;  }
  if(contrib_legend) { delete contrib_legend; }
  if(data_graph)     { delete data_graph;     }
  if(canv)           { delete canv;           }
  for(auto h : pdf_contrib)        { delete h; }
  for(auto h : pdf_contrib_coarse) { delete h; }

}

//______________________________________________________________________________________________________________________
void save_to_file(TCanvas* canv, std::string title_text, bool close_file)
{
  // If we don't know in advance when we're going to close the file (or are too lazy to work it out),
  // this lets us just write a blank page at the end and close the file.
  if(!canv) {
    close_file = true;
    canv = new TCanvas("blank", "blank", _canv_size_x, _canv_size_y);
  }

  // Allows us to write plots to multiple pages of the output PDF
  std::string file_suffix = "";
  if(_first_draw && !close_file) {
    _first_draw = false;
    file_suffix = "(";
  }
  else if(close_file) {
    file_suffix = ")";
  }

  std::string output_filename_use = _output_filename + ".pdf" + file_suffix;

  // Optionally write the plot title as the page bookmark (useful for navigating large documents)
  if(title_text == "") {
    title_text = "-";
  }
  std::string bookmark_name = "Title:" + title_text;

  canv->Print(output_filename_use.c_str(), bookmark_name.c_str());

  if(close_file) {
    std::cout << "All plots written to file \'" << _output_filename + ".pdf" << "\'. Closing File." << std::endl;
  }
}

//______________________________________________________________________________________________________________________
void insert_section_break(std::string title_text)
{
  TCanvas* canv = new TCanvas("break", "break", _canv_size_x, _canv_size_y);

  TLatex * title = new TLatex(0.15, 0.5, title_text.c_str());
  title->SetTextSize(0.1);
  title->Draw();

  save_to_file(canv, "-");
}

//______________________________________________________________________________________________________________________
TLatex* make_caveat(TH1* hist_to_draw, double y_max)
{
  double tx = hist_to_draw->GetXaxis()->GetBinLowEdge(hist_to_draw->GetXaxis()->GetLast());
  TLatex * caveat = new TLatex(0.65*tx, 1.02*y_max, _caveat_text_str.c_str());
  caveat->SetTextSize(_title_text_size);
  return caveat;
}

//______________________________________________________________________________________________________________________
TLatex* make_title(TString title_text, double ymax)
{
  double title_xlow = 0.0;
  TLatex * title = new TLatex(title_xlow, 1.02*ymax, title_text.Data());
  title->SetTextSize(_title_text_size);
  return title;
}

//______________________________________________________________________________________________________________________
std::string make_title_text(conv::BeamConf_t b, conv::Sample_t s, conv::ReacMode_t m, bool no_tex, bool no_spaces)
{
  std::string title_text = get_beam_current_name_simple(b) + " " + get_sample_name_simple(s, b, no_tex);
  if(m != conv::kmUndefined) {
    std::string padding = (_thesis_format) ? "; " : " ";
    if(no_tex) { title_text += padding + vutil::AsString   (m); }
    else       { title_text += padding + vutil::AsStringTeX(m); }
  }
  if(no_spaces) {
    std::replace( title_text.begin(), title_text.end(), ' ', '_');
  }
  return title_text;
}

//______________________________________________________________________________________________________________________
std::string get_beam_current_name_simple(conv::BeamConf_t b)
{
  std::string beam_str = "";
  if     (b == conv::kbJPARC_FHC) { beam_str = (_thesis_format) ? "#nu-mode"       : "FHC"; }
  else if(b == conv::kbJPARC_RHC) { beam_str = (_thesis_format) ? "#bar{#nu}-mode" : "RHC"; }
  else {
    std::cerr << "ERROR: Unrecogised beam \'" << vutil::AsString(b) << "\'" << std::endl;
    vutil::Die(1);
  }
  return beam_str;
}

//______________________________________________________________________________________________________________________
std::string get_sample_name_simple(conv::Sample_t s, conv::BeamConf_t b, bool no_tex)
{
  std::string sample_str = "";
  if(!_thesis_format) {
    if     (s == conv::ks1Rmu   ) { sample_str = "1R#mu";     }
    else if(s == conv::ks1Re    ) { sample_str = "1Re";       }
    else if(s == conv::ks1Re1dcy) { sample_str = "1Re 1d.e."; }
    else {
      std::cerr << "ERROR: Unrecogised sample \'" << vutil::AsString(s) << "\'" << std::endl;
      vutil::Die(1);
    }
  }
  else {
    bool rhc = (b == conv::kbJPARC_RHC);
    if     (s == conv::ks1Rmu   ) { sample_str = (rhc) ? "#bar{#nu}_{#mu}CCQE-like" : "#nu_{#mu}CCQE-like"; }
    else if(s == conv::ks1Re    ) { sample_str = (rhc) ? "#bar{#nu}_{e}CCQE-like"   : "#nu_{e}CCQE-like";   }
    else if(s == conv::ksMultiRe) { sample_str = "#nu_{e}CC1#pi^{+}-like"; }
    else {
      std::cerr << "ERROR: Unrecogised sample \'" << vutil::AsString(s) << "\'" << std::endl;
      vutil::Die(1);
    }
  }

  if(no_tex) {
    sample_str.erase(std::remove(sample_str.begin(), sample_str.end(), '#'), sample_str.end());
  }

  return sample_str;
}

//______________________________________________________________________________________________________________________
void rebin(TH1D* hist, conv::BeamConf_t b, conv::Sample_t s, conv::KineVar_t kv)
{
  if(_rebin_mulike && vutil::AsString(kv) == "Ereco") {
    // For T2K, the binning used in the FHC1Rmu and RHC1Rmu plots are coarser than the fit binning.
    // As problems are caused when filling the MC templates when the simple binning has fewer bins than the opt
    // binning, rebin the relevant hists here instead of changing the simple binning.
    if(s == conv::ks1Rmu && b == conv::kbJPARC_FHC) {
      hist->Rebin(2);
    }
    else if(s == conv::ks1Rmu && b == conv::kbJPARC_RHC) {
      hist->Rebin(3);
    }
  }
}

//______________________________________________________________________________________________________________________
// Scan for the max bin content in the desired x-axis range and summed across reaction modes
double get_max_in_xrange_by_mode(double xmax, std::vector<TH1D*> hists_by_mode, const int NReacModes)
{
  double max_n = 0.;
  for(int ibin = 1; ibin <= hists_by_mode[0]->GetNbinsX(); ibin++) {
    if(hists_by_mode[0]->GetBinCenter(ibin) < xmax ) {
      double sum = 0;
      for(int imode = 0; imode < NReacModes; imode++) {
        if (hists_by_mode[imode]) {
          double val = hists_by_mode[imode]->GetBinContent(ibin);
          sum += val;
        }
      }
      max_n = TMath::Max(max_n, sum);
    }//<emax
  }//ibin
  return max_n;
}

//______________________________________________________________________________________________________________________
// Scan for the max bin content in the desired x-axis range
double get_max_in_xrange(double xmax, TH1D* hist, bool with_error)
{
  double max_n = 0.;
  for(int ibin = 1; ibin <= hist->GetNbinsX(); ibin++) {
    if(hist->GetBinCenter(ibin) < xmax ) {
      double val = hist->GetBinContent(ibin);
      if(with_error) { val += hist->GetBinError(ibin); }
      max_n = TMath::Max(max_n, val);
    }//<emax
  }//ibin
  return max_n;
}

//______________________________________________________________________________________________________________________
// Scan for the max bin content in the desired x-axis range
double get_max_in_xrange_asymm_err(double xmax, TGraphAsymmErrors* graph)
{
  double max_n = 0.;
  for(int ipt = 0; ipt < graph->GetN(); ipt++) {
    double X = 0;
    double Y = 0;
    graph->GetPoint(ipt, X, Y);
    if(X < xmax) {
      double val = Y + graph->GetErrorYhigh(ipt);
      max_n = TMath::Max(max_n, val);
    }//<emax
  }//ipt
  return max_n;
}

//______________________________________________________________________________________________________________________
double get_kinevar_max(conv::Sample_t s, conv::KineVar_t kv)
{
  double max = -999;
  if(kv == conv::kvEreco) {
    if(s == conv::ks1Rmu) { max = _erec_max_mulike; }
    else                  { max = _erec_max_elike;  }
  }
  else if(kv == conv::kvThetareco) {
    max = _theta_max_elike;
  }
  else {
    std::cerr << "ERROR: Unrecogised kinevar \'" << vutil::AsString(kv) << "\'" << std::endl;
    vutil::Die(1);
  }
  return max;
}

//______________________________________________________________________________________________________________________
TGraph* get_empty_bin_markers(TH2D* pdf)
{
  if(!_mark_empty_bins_2d) { return nullptr; }
  std::vector<double> vec_empty_x;
  std::vector<double> vec_empty_y;
  int n_empty_bins = 0;
  for(int ix = 1; ix < pdf->GetNbinsX() +1; ++ix) {
    for(int iy = 1; iy < pdf->GetNbinsY() +1; ++iy) {
      if(TMath::Abs(pdf->GetBinContent(ix,iy)) < 1e-10) {
        n_empty_bins++;
        vec_empty_x.push_back(pdf->GetXaxis()->GetBinCenter(ix));
        vec_empty_y.push_back(pdf->GetYaxis()->GetBinCenter(iy));
      }
    }
  }

  TVectorD empty_x(n_empty_bins, vec_empty_x.data());
  TVectorD empty_y(n_empty_bins, vec_empty_y.data());
  TGraph* empty_markers = new TGraph(empty_x, empty_y);

  empty_markers->SetMarkerColor(kAzure - 3);
  empty_markers->SetMarkerStyle(kFullCrossX);
  empty_markers->SetMarkerSize(2);

  return empty_markers;
}

//______________________________________________________________________________________________________________________