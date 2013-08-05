//plot an infinite number of contours
#include <algorithm>

void contour_plots_nfits(
	   const char * contour_filelist,    // filename for the list of contours.
	   /* Should be in format "mixing_in\tdm_in\t68%_contour_name\t90%_contour_name\tlegend_tag\tpath\to\file.root" */
	   const int    ncontours,           // number of contours to draw
	   const char * output_tag,
	   bool         plot_68 = false,     // plot 68% contours?
	   bool         plot_90 = true,      // plot 90% contours?
	   bool         plot_bf = true,      // plot best-fit markers?
	   bool         plot_in = false,     // plot input-point markers?
	   double       mixing23_min = 0.9,  // edges of the canvas
	   double       mixing23_max = 1.0,
	   double       dm32sq_min   = 0.002,
	   double       dm32sq_max   = 0.003,
	   const char * mixingparam  = "sinsq2theta", /* for the x-axis */
	   double       leg_font_size = 0.05,
	   bool         cl_and_bf_in_leg = false
)
{
  using namespace std;
  using namespace t2kana;

  // param space
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.16);
  // dm32sq, mixing23 range and step
  double mixing23_step  = 0.005;
  double dm32sq_step    = 5E-5;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int) ((1    - 0   ) / mixing23_step);
  int dm32sq_nbins   = 1 + (int) ((6E-3 - 1E-3) / dm32sq_step  );
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       0 - 0.5*mixing23_step, 
       1 + 0.5*mixing23_step,
       dm32sq_nbins, 
       1E-3 - 0.5*dm32sq_step, 
       6E-3 + 0.5*dm32sq_step);
  osc_param_space->GetXaxis()->SetRangeUser(mixing23_min, mixing23_max);
  osc_param_space->GetYaxis()->SetRangeUser(dm32sq_min,   dm32sq_max);
  if(mixingparam == "sinsq2theta")
    osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  else if(mixingparam == "sinsqtheta")
    osc_param_space->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  //osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  // book histograms etc
  const int nmax = 100;
  assert(ncontours < nmax);
  TH2D    * h68          [nmax]; 
  TH2D    * h90          [nmax]; 
  TMarker * bestfit_point[nmax];
  TMarker * input_point  [nmax];
  //Legend
  TLegend * leg = new TLegend(0.17,0.68,0.53,0.93);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->SetTextSize(leg_font_size);
  TLegend * legbf = new TLegend(0.17,0.18,0.53,0.3);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->SetTextSize(leg_font_size);

  // Read in contours
  ifstream fitinp(contour_filelist, ios::in);
  vector<int> skip;
  for(int ifit=0; ifit<ncontours; ifit++) {

    double chisq_min        = 999999.0;
    double mixing23_bestfit = -1;
    double dm32sq_bestfit   = -1;
    double mixing23_input   = -1;
    double dm32sq_input     = -1;
    string fit_result_file  = "";
    string h68_histname     = "";
    string h90_histname     = "";
    string tag              = "";

    fitinp >> mixing23_input >>  dm32sq_input >> h68_histname >> h90_histname >> tag >> fit_result_file;

    cout << "Reading contours " << h68_histname << " and " << h90_histname << " from : " << fit_result_file << endl;

    TFile * dataFile = TFile::Open(fit_result_file.c_str(), "READ");
    if(!dataFile) {
      skip.push_back(ifit);
      continue;
    }
    dataFile->cd();
    TH2D * htemp68 = 0, *htemp90 = 0;
    dataFile->GetObject(h68_histname.c_str(), htemp68);
    dataFile->GetObject(h90_histname.c_str(), htemp90);
    h68[ifit] = new TH2D(*htemp68);
    h90[ifit] = new TH2D(*htemp90);
    if(!h68[ifit] || !h90[ifit]) {
      skip.push_back(ifit);
      continue;
    }

    //swap '$' with ' ' in tag
    for (int i = 0; i < tag.length(); ++i) {
      if (tag[i] == '$')
	tag[i] = ' ';
    }

    //find best-fit point
    FitParamList * fparm = 0;
    dataFile->GetObject("fparm",fparm);
    //get it from the histogram
    if(!fparm) {
      cout << "Retreiving best-fit parameters from the chisq surface" << endl;
      for(int i = 1; i <= h68[ifit]->GetNbinsX(); i++) {
	for(int j = 1; j <= h68[ifit]->GetNbinsY(); j++) {
	  double chisq = h68[ifit]->GetBinContent(i,j);
	  if(chisq < chisq_min) {
	    chisq_min        = chisq;
	    mixing23_bestfit = h68[ifit]->GetXaxis()->GetBinCenter(i);
	    dm32sq_bestfit   = h68[ifit]->GetYaxis()->GetBinCenter(j);
	  }
	}//j loop
      }//i loop
    }
    //get it from the fit param list
    else {
      cout << "Retreiving best-fit parameters from the FitParamList" << endl;
      mixing23_bestfit = fparm->BestFitValue(conv::kPhypMixing23);
      dm32sq_bestfit   = fparm->BestFitValue(conv::kPhypDm23sq);
      chisq_min        = fparm->GetBestFitChisq();
    }
    cout << "chi\2_{min} = " << chisq_min << endl;
    cout << "dm32sq_{bestfit} = " << dm32sq_bestfit << endl;
    cout << "sinsq2theta_{bestfit} = " << mixing23_bestfit << endl;

    input_point   [ifit] = new TMarker(mixing23_input,   dm32sq_input,   29);
    bestfit_point [ifit] = new TMarker(mixing23_bestfit, dm32sq_bestfit, 8);

    //style
    input_point  [ifit]->SetMarkerSize(0.5);
    if(ncontours <= 3) {
      if(ifit == 0) {
	bestfit_point[ifit]->SetMarkerSize(1.8);
	bestfit_point[ifit]->SetMarkerStyle(29);
      }
      else if(ifit == 1) {
	bestfit_point[ifit]->SetMarkerSize(1.8);
	bestfit_point[ifit]->SetMarkerStyle(24);
      }
      else if(ifit == 2) {
	bestfit_point[ifit]->SetMarkerSize(1.8);
	bestfit_point[ifit]->SetMarkerStyle(25);
      }
    }
    else
      bestfit_point[ifit]->SetMarkerSize(0.5);
    h68[ifit]->SetLineWidth(2);
    h90[ifit]->SetLineWidth(2);
    h68[ifit]->SetLineStyle(2);
    h90[ifit]->SetLineStyle(1);
    //color
    Color_t icolor = kBlack;
    switch(ifit%10) {
    case 0:
      icolor = kBlack; break;
    case 1:
      icolor = kRed; break;
    case 2:
      icolor = kBlue; break;
    case 3:
      icolor = kGreen+2; break;
    case 4:
      icolor = kOrange+7; break;
    case 5:
      icolor = kYellow+2; break;
    case 6:
      icolor = kMagenta+1; break;
    case 7:
      icolor = kCyan+1; break;
    case 8:
      icolor = kGray; break;
    case 9:
      icolor = kViolet+5; break;
    }//switch(ifit%10)
    bestfit_point[ifit]->SetMarkerColor(icolor);
    input_point  [ifit]->SetMarkerColor(icolor);
    h68          [ifit]->SetLineColor(icolor);
    h90          [ifit]->SetLineColor(icolor);
    const char * percent = "%";
    if(cl_and_bf_in_leg) {
      leg->AddEntry(h68[ifit], Form("%s 68%s CL",tag.c_str(),percent), "L");
      leg->AddEntry(h90[ifit], Form("%s 90%s CL",tag.c_str(),percent), "L");
      legbf->AddEntry(bestfit_point[ifit], Form("%s best fit",tag.c_str()), "P");
    }
    else
      leg->AddEntry(h90[ifit], tag.c_str(), "L");
  }//ifit

  //Plot
  osc_param_space->Draw();
  for(int ifit=0; ifit<ncontours; ifit++) {
    //check for null file/spline
    bool skip_this = false;
    for(int i = 0; i<skip.size(); i++) {
      if(skip.at(i) == ifit) {
	skip_this = true;
	break;
      }
    }
    if(skip_this)
      continue;
    if(plot_68)
      h68[ifit]->Draw("cont3same");
    if(plot_90)
      h90[ifit]->Draw("cont3same");
    if(plot_bf)
      bestfit_point[ifit]->Draw();
    if(plot_in)
      input_point[ifit]->Draw();
  }
  leg->Draw();
  if(cl_and_bf_in_leg) {
    legbf->Draw();
  }

  c->SaveAs(Form("contours_%s.pdf",output_tag));
  c->SaveAs(Form("contours_%s.eps",output_tag));
  c->SaveAs(Form("contours_%s.C",  output_tag));

  if(mixingparam == "sinsq2theta") {
    double linewidth_ext = 2;
    //get externals
    TFile fminos("$T2KOSCANA/inputs/contours/MINOS_NumuDisappearance_LogLikelihoodSurface_May2013.root");
    TH2D * ll_minos = 0;
    fminos.GetObject("LogLikelihoodSurface",ll_minos);
    ll_minos->SetContour(1);
    ll_minos->SetContourLevel(0,4.61/2);
    ll_minos->SetLineStyle(3);
    ll_minos->SetLineWidth(linewidth_ext);
    TFile fskle("$T2KOSCANA/inputs/contours/contours_skloe.root");
    TGraph * g_skle = 0;
    fskle.GetObject("sk3f_loe_90cl",g_skle);
    g_skle->SetLineColor(kGreen+2);
    g_skle->SetLineStyle(3);
    g_skle->SetLineWidth(linewidth_ext);
    TFile fskz("$T2KOSCANA/inputs/contours/contours_skzen.root");
    TGraph * g_skzen = 0;
    fskz.GetObject("sk3f_23vars_overlay_normal_q13_fixed",g_skzen);
    g_skzen->SetLineColor(kBlue);
    g_skzen->SetLineStyle(3);
    g_skzen->SetLineWidth(linewidth_ext);
    TFile ft2k_run123_octant1("$T2KOSCANA/2012results/NEWRESULTS/datafit/dchisq_2all.root");
    TH2D * h_t2k_o1 = 0;
    ft2k_run123_octant1.GetObject("hdchisq",h_t2k_o1);
    // TH2D * h_t2k_o1_68 = new TH2D(*h_t2k_o1);
    // h_t2k_o1_68->SetContour(1);
    // h_t2k_o1_68->SetContourLevel(0,2.30);
    // h_t2k_o1_68->SetLineStyle(2);
    // h_t2k_o1_68->SetLineWidth(linewidth_ext);
    TH2D * h_t2k_o1_90 = new TH2D(*h_t2k_o1);
    h_t2k_o1_90->SetContour(1);
    h_t2k_o1_90->SetContourLevel(0,4.61);
    h_t2k_o1_90->SetLineStyle(3);
    h_t2k_o1_90->SetLineWidth(linewidth_ext);
    TFile ft2k_run123_octant2("$T2KOSCANA/2012results/DATA/octant2_23only/combination/fitresults_3.root");
    TH2D * h_t2k_o2 = 0;
    ft2k_run123_octant2.GetObject("hdchisq",h_t2k_o2);
    // TH2D * h_t2k_o2_68 = new TH2D(*h_t2k_o2);
    // h_t2k_o2_68->SetContour(1);
    // h_t2k_o2_68->SetContourLevel(0,2.30);
    // h_t2k_o2_68->SetLineStyle(2);
    // h_t2k_o2_68->SetLineWidth(linewidth_ext);
    TH2D * h_t2k_o2_90 = new TH2D(*h_t2k_o2);
    h_t2k_o2_90->SetContour(1);
    h_t2k_o2_90->SetContourLevel(0,4.61);
    h_t2k_o2_90->SetLineStyle(3);
    h_t2k_o2_90->SetLineWidth(linewidth_ext);
    //style
    h_t2k_o1_90->SetLineColor(kBlack);
    h_t2k_o2_90->SetLineColor(kRed);
    ll_minos->SetLineColor(kBlue);
    g_skle->SetLineColor(kGreen+2);
    g_skzen->SetLineColor(kMagenta-4);
    //legend
    g_skle->SetName("g_skle");
    g_skzen->SetName("g_skzen");
    TLegend * legext = new TLegend(0.15,0.2,0.5,0.4);
    legext->SetFillColor(0);
    legext->SetFillStyle(0);
    legext->SetBorderSize(0);
    legext->SetTextFont(132);
    legext->SetTextSize(leg_font_size);
    legext->AddEntry(ll_minos, "MINOS 2013 2#nu", "l");
    legext->AddEntry(g_skle, "SuperK L/E 2012 2#nu", "l");
    legext->AddEntry(g_skzen, "SuperK zenith 2012 3#nu", "l");
    legext->AddEntry(h_t2k_o1_90, "T2K Run 1-3 3#nu #theta_{23} #leq #pi/4","l");
    legext->AddEntry(h_t2k_o2_90, "T2K Run 1-3 3#nu #theta_{23} #geq #pi/4","l");
    //draw
    osc_param_space->GetXaxis()->SetRangeUser(0.8,1);
    osc_param_space->GetYaxis()->SetRangeUser(2E-3,3.5E-3);
    osc_param_space->Draw();
    h_t2k_o1_90->Draw("cont3same");
    h_t2k_o2_90->Draw("cont3same");
    ll_minos->Draw("cont3same");
    g_skle->Draw("same");
    g_skzen->Draw("same");
    for(int ifit=0; ifit<ncontours; ifit++) {
      //check for null file/spline
      bool skip_this = false;
      for(int i = 0; i<skip.size(); i++) {
	if(skip.at(i) == ifit) {
	  skip_this = true;
	  break;
	}
      }
      if(skip_this)
	continue;
      h90[ifit]->Draw("cont3same");
      if(plot_bf)
	bestfit_point[ifit]->Draw();
      if(plot_in)
	input_point[ifit]->Draw();
    }
    leg->Draw();
    legext->Draw();

    c->SaveAs(Form("contours_%s_external.pdf",output_tag));
  }//sinsq2theta externals

}
