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
	   const char * mixingparam  = "sinsq2theta" /* for the x-axis */
)
{
  using namespace std;

  // param space
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.16);
  // dm32sq, mixing23 range and step
  double mixing23_step  = 0.005;
  double dm32sq_step    = 5E-5;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
  int dm32sq_nbins   = 1 + (int) ((dm32sq_max   - dm32sq_min  ) / dm32sq_step  );
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm32sq_nbins, 
       dm32sq_min - 0.5*dm32sq_step, 
       dm32sq_max + 0.5*dm32sq_step);
  if(mixingparam == "sinsq2theta")
    osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  else if(mixingparam == "sinsqtheta")
    osc_param_space->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
  //osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  osc_param_space->GetYaxis()->SetTitle("#Deltam^{2} (eV^{2}/c^{4})");
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
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);

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

    fitinp >> mixing23_input >>  dm32sq_input >> mixing23_bestfit >> dm32sq_bestfit >> h68_histname >> h90_histname >> tag >> fit_result_file;

    if(tag.find("pddeff") != string::npos)
      fit_result_file = "";

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

    //find best-fit point
    bool mixing23_min_from_grid = true;
    if(mixing23_bestfit > 0)
      mixing23_min_from_grid = false;
    bool dm32sq_min_from_grid = true;
    if(mixing23_bestfit > 0)
      dm32sq_min_from_grid = false;
    for(int i = 1; i <= h68[ifit]->GetNbinsX(); i++) {
      for(int j = 1; j <= h68[ifit]->GetNbinsY(); j++) {
        double chisq = h68[ifit]->GetBinContent(i,j);
	if(chisq < chisq_min) {
	  chisq_min        = chisq;
	  if(mixing23_min_from_grid)
	    mixing23_bestfit = h68[ifit]->GetXaxis()->GetBinCenter(i);
	  if(dm32sq_min_from_grid)
	    dm32sq_bestfit   = h68[ifit]->GetYaxis()->GetBinCenter(j);
	}
      }//j loop
    }//i loop
    cout << "chi\2_{min} = " << chisq_min << endl;
    cout << "dm32sq_{bestfit} = " << dm32sq_bestfit << endl;
    cout << "sinsq2theta_{bestfit} = " << mixing23_bestfit << endl;

    input_point   [ifit] = new TMarker(mixing23_input,   dm32sq_input,   29);
    bestfit_point [ifit] = new TMarker(mixing23_bestfit, dm32sq_bestfit, 8);

    //style
    input_point  [ifit]->SetMarkerSize(0.5);
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
      icolor = kYellow-7; break;
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

  c->SaveAs(Form("contours_%s.pdf",output_tag));

  /*
  //get externals
  TFile fext("$T2KOSCANA/inputs/contours/contours_nont2k.root");
  TGraph * g_minos = 0;
  fext.GetObject("minos2012_combined_90CL",g_minos);
  g_minos->SetLineColor(kBlack);
  TGraph * g_skle = 0;
  fext.GetObject("sk2012_le2nu_90CL",g_skle);
  g_skle->SetLineColor(kGreen+2);
  TGraph * g_skzen = 0;
  fext.GetObject("sk2012_zenith3nu_90CL",g_skzen);
  g_skzen->SetLineColor(kBlue);
  //names
  g_minos->SetName("g_minos");
  g_skle->SetName("g_skle");
  g_skzen->SetName("g_skzen");

  g_minos->Draw("l");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry("g_minos", "MINOS 2012 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->Draw();
  */
}
