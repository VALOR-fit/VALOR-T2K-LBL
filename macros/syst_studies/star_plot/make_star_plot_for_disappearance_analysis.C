//
// ------------------------------------------------------------------------------------------------
// Use the combined outputs from test_bestfit_point_shift_due_to_syst_twk.C calls in order to
// construct the `star plot' showing the effect of systematics on best-fit params.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//      

//
// formatting for best-fit point
//
const int kTruePointMarkerStyle  = 20;
const int kTruePointMarkerSize   =  1;
const int kTruePointMarkerColour =  kBlack;

//
// formatting for ellipse showing `cut' for considering a systematic in the fit
//
const int    kCutLineStyle = kDashed;
const int    kCutFillStyle = 0; // transparent
const double kCutLineWidth = 1;

//
// formatting for star-plot arrows
//
const double kSystShiftArrowSize = 0.02;
const double kSystShiftArrowLineWidth = 1;

const int kNSyst = 33; // maximum number of systematics (with unique formatting style)

int kArrowColour   [kNSyst] = {
   kBlack, kBlue, kAzure+10, kCyan-2, kGreen, kGreen+2, kYellow+2, kOrange+7, kRed+1, kMagenta, kMagenta+2,
   kBlack, kBlue, kAzure+10, kCyan-2, kGreen, kGreen+2, kYellow+2, kOrange+7, kRed+1, kMagenta, kMagenta+2,
   kBlack, kBlue, kAzure+10, kCyan-2, kGreen, kGreen+2, kYellow+2, kOrange+7, kRed+1, kMagenta, kMagenta+2
};
int kArrowLineStyle[kNSyst] = {
   kSolid,  kSolid,  kSolid,  kSolid,  kSolid,  kSolid,  kSolid,  kSolid,  kSolid,  kSolid,  kSolid, 
   kDashed, kDashed, kDashed, kDashed, kDashed, kDashed, kDashed, kDashed, kDashed, kDashed, kDashed,
   kDotted, kDotted, kDotted, kDotted, kDotted, kDotted, kDotted, kDotted, kDotted, kDotted, kDotted
};

void make_star_plot_for_disappearance_analysis
(
  string       tweaked_toy_mc_expt_fit_data, /* merged outputs of test_bestfit_point_shift_due_to_syst_twk.C calls */
  double       true_mixing23,
  double       true_dm23sq,
  double       stat_av_mixing23,
  double       stat_av_dm23sq
)
{
  using namespace std;
  using namespace t2kana;

  gStyle->SetOptStat(0);

  if(true_dm23sq   < -1) true_dm23sq   = conv::kGlobFitDm23sq;
  if(true_mixing23 < -1) true_mixing23 = conv::kGlobFitSinsq2Theta23;

  //
  // Plot osc param frame
  //

  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.5;
  double mixing23_max   = 1.;
  double mixing23_step  = 0.001;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 4E-3;
  double dm23sq_step    = 1E-5;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + TMath::CeilNint ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + TMath::CeilNint ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  // param space
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
  osc_param_space->GetXaxis()->SetTitle("sin^{2}2#theta_{23}"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);
  osc_param_space->GetXaxis()->CenterTitle();
  osc_param_space->GetYaxis()->CenterTitle();
  TCanvas * c = new TCanvas("c","T2K numu disappearance analysis star plot",20,20,1000,1000);
  c->SetLeftMargin(0.2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  osc_param_space->GetYaxis()->SetTitleOffset(1.75);
  osc_param_space->Draw();

  // Add legend
  TLegend * legend = new TLegend(0.80, 0.16, 0.95, 0.95);
  legend->SetFillColor(0);
  legend->SetBorderSize(1);
  legend->SetHeader("Systematic");
  

  //
  // Add true/input oscillation parameter point
  //

  TMarker * true_point = new TMarker(true_mixing23,true_dm23sq,kTruePointMarkerStyle);
  true_point->SetMarkerSize (kTruePointMarkerSize);
  true_point->SetMarkerColor(kTruePointMarkerColour);
  true_point->Draw();

  //
  double kRmixing23 = TMath::Abs(stat_av_mixing23);
  double kRdm23sq   = TMath::Abs(stat_av_dm23sq);
  TEllipse * cut = new TEllipse(true_mixing23,true_dm23sq, kRmixing23 / 4, kRdm23sq / 4);
  cut->SetLineWidth(kCutLineWidth);
  cut->SetLineStyle(kCutLineStyle);
  cut->SetFillStyle(kCutFillStyle);
  cut->Draw();

  //
  // Read tweaked toy MC fit results and assemble the star...
  //

  cout << "Reading tweaked toy MC fit data from: " << tweaked_toy_mc_expt_fit_data << endl;

  ifstream fit_data(tweaked_toy_mc_expt_fit_data.c_str());
  if (!fit_data.good() ) {
      cerr << "** Can't open data file: " << tweaked_toy_mc_expt_fit_data << endl;
      return;
  }

  int isyst=0;

  while(1) {

      // skip header lines staring with #
      if(fit_data.peek() == '#') {
         fit_data.ignore(1000, '\n');

       } else {
         string systematic =  "";
         double bestfit_mixing12_p1sig  = 0;
         double bestfit_mixing13_p1sig  = 0.;
         double bestfit_mixing23_p1sig  = 0.;
         double bestfit_deltacp_p1sig   = 0.;
         double bestfit_dm12sq_p1sig    = 0.;
         double bestfit_dm23sq_p1sig    = 0.;
         double bestfit_mixing12_m1sig  = 0.;
         double bestfit_mixing13_m1sig  = 0.;
         double bestfit_mixing23_m1sig  = 0.;
         double bestfit_deltacp_m1sig   = 0.;
         double bestfit_dm12sq_m1sig    = 0.;
         double bestfit_dm23sq_m1sig    = 0.;

         fit_data >> systematic 
                  >> bestfit_mixing12_p1sig 
                  >> bestfit_mixing13_p1sig 
                  >> bestfit_mixing23_p1sig 
                  >> bestfit_deltacp_p1sig  
                  >> bestfit_dm12sq_p1sig   
                  >> bestfit_dm23sq_p1sig   
                  >> bestfit_mixing12_m1sig 
                  >> bestfit_mixing13_m1sig 
                  >> bestfit_mixing23_m1sig 
                  >> bestfit_deltacp_m1sig  
                  >> bestfit_dm12sq_m1sig   
                  >> bestfit_dm23sq_m1sig;

         if(fit_data.eof()) break;            
         cout << "Read best-fit osc. param shifts for syst: " << systematic << endl;

         // plot arrows
         TArrow * arr_p1sig = new TArrow(true_mixing23,true_dm23sq,bestfit_mixing23_p1sig,bestfit_dm23sq_p1sig,kSystShiftArrowSize,"|>");
         arr_p1sig->SetLineWidth(kSystShiftArrowLineWidth);
         arr_p1sig->SetLineColor(kArrowColour[isyst]);
         arr_p1sig->SetFillColor(kArrowColour[isyst]);
         arr_p1sig->SetLineStyle(kArrowLineStyle[isyst]);
         arr_p1sig->Draw();
         TArrow * arr_m1sig = new TArrow(true_mixing23,true_dm23sq,bestfit_mixing23_m1sig,bestfit_dm23sq_m1sig,kSystShiftArrowSize,"|>");
         arr_m1sig->SetLineWidth(kSystShiftArrowLineWidth);
         arr_m1sig->SetLineColor(kArrowColour[isyst]);
         arr_m1sig->SetFillColor(kArrowColour[isyst]);
         arr_m1sig->SetLineStyle(kArrowLineStyle[isyst]);
	 if(systematic != "xsec_sf")
	   arr_m1sig->Draw();

         // add legend entry
	 if(systematic.find("osc_") != string::npos)
	   systematic.erase(systematic.begin(), systematic.begin()+4);
	 if(systematic == "nominal")
	   legend->AddEntry(arr_p1sig, "nominal", "L");
	 else
	   legend->AddEntry(arr_p1sig, util::AsStringPaperRooTex(util::PhysParamIdFromString(systematic)).c_str(), "L");

         isyst++;
         // reset and re-use formatting styles, shouldn't have anywhere near that many systematics
         // in order to be able to see anything on the plot
         if(isyst==kNSyst) isyst=0; 

       }//not a comment line
  }//while(1)

  fit_data.close();

  if(TMath::Abs(1.0-true_mixing23) < conv::kASmallNum) {
    legend->SetX1(0.2);
    legend->SetX2(0.35);
  }
  legend->Draw();

  //Save a couple of zoom levels of the canvas automatically
  vector<string> token = util::Tokens (tweaked_toy_mc_expt_fit_data, "-");
  string tokened = token[token.size()-1];
  vector<string> groupname = util::Tokens(tokened, ".");
  string group = groupname[0];

  c->SaveAs(Form("%s-%f-%f_big.eps",group.c_str(),true_mixing23,true_dm23sq));

  if(TMath::Abs(0.9-true_mixing23) < conv::kASmallNum && TMath::Abs(2.4E-3-true_dm23sq) < conv::kASmallNum) {
    osc_param_space->GetXaxis()->SetRangeUser(0.878,0.916);
    osc_param_space->GetYaxis()->SetRangeUser(0.00227,0.00253);
    c->SaveAs(Form("%s-%f-%f.eps",group.c_str(),true_mixing23,true_dm23sq));

    osc_param_space->GetXaxis()->SetRangeUser(0.895,0.905);
    osc_param_space->GetYaxis()->SetRangeUser(0.002375,0.002425);
    c->SaveAs(Form("%s-%f-%f_zoom.eps",group.c_str(),true_mixing23,true_dm23sq));

    osc_param_space->GetXaxis()->SetRangeUser(0.898,0.902);
    osc_param_space->GetYaxis()->SetRangeUser(0.002395,0.002405);
    c->SaveAs(Form("%s-%f-%f_zoom2.eps",group.c_str(),true_mixing23,true_dm23sq));
  }

  else if(TMath::Abs(1.0-true_mixing23) < conv::kASmallNum && TMath::Abs(2.4E-3-true_dm23sq) < conv::kASmallNum) {
    osc_param_space->GetXaxis()->SetRangeUser(0.982,1.);
    osc_param_space->GetYaxis()->SetRangeUser(0.00227,0.00253);
    c->SaveAs(Form("%s-%f-%f.eps",group.c_str(),true_mixing23,true_dm23sq));

    osc_param_space->GetXaxis()->SetRangeUser(0.995,1.);
    osc_param_space->GetYaxis()->SetRangeUser(0.002375,0.002425);
    c->SaveAs(Form("%s-%f-%f_zoom.eps",group.c_str(),true_mixing23,true_dm23sq));

    osc_param_space->GetXaxis()->SetRangeUser(0.998,1.);
    osc_param_space->GetYaxis()->SetRangeUser(0.002395,0.002405);
    c->SaveAs(Form("%s-%f-%f_zoom2.eps",group.c_str(),true_mixing23,true_dm23sq));
  }
  
  //cout << osc_param_space->GetXaxis()->GetNdivisions() << endl;
  osc_param_space->GetXaxis()->SetNdivisions(507,true);
  osc_param_space->GetXaxis()->SetRangeUser(true_mixing23 - (1.1* stat_av_mixing23), true_mixing23 + (1.1* stat_av_mixing23));
  osc_param_space->GetYaxis()->SetRangeUser(true_dm23sq - (1.1* stat_av_dm23sq), true_dm23sq + (1.1* stat_av_dm23sq));
  c->SaveAs(Form("%s-%f-%f_testbig.eps",group.c_str(),true_mixing23,true_dm23sq));

  osc_param_space->GetXaxis()->SetRangeUser(true_mixing23 - (0.6* stat_av_mixing23), true_mixing23 + (0.6* stat_av_mixing23));
  osc_param_space->GetYaxis()->SetRangeUser(true_dm23sq - (0.6* stat_av_dm23sq), true_dm23sq + (0.6* stat_av_dm23sq));
  c->SaveAs(Form("%s-%f-%f_testzoom.eps",group.c_str(),true_mixing23,true_dm23sq));
  
  cout << "Done! " << endl;
}
