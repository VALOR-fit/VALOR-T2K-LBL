//
// ------------------------------------------------------------------------------------------------
// Plot 3-flavour oscillation probabilities :
//  - P (numu     -> numu   )
//  - P (numubar  -> numubar)
//  - P (nue      -> nue    )
//  - P (nuebar   -> nuebar )
//  - P (numu     -> nue    )
//
// Probabilities are plotted for the input oscillated params.
// Defaults are used for unset oscillation params.
// 
// Examples:
// % root2k -l -b -q '$T2KOSCANA/macros/misc/osc/plot_3f_osc_prob.C(0.957,0.00239)'
// % root2k -l -b -q '$T2KOSCANA/macros/misc/osc/plot_3f_osc_prob.C(1.000,0.00244)'
//
// ------------------------------------------------------------------------------------------------
//

void plot_3f_osc_prob
(
  double       mixing23   = -99,                            /* */
  double       dm23sq     = -99,                            /* */
  double       mixing12   = -99,                            /* */
  double       mixing13   = -99,                            /* */
  double       deltacp    = -99,                            /* */
  double       dm12sq     = -99                             /* */
)
{
  using namespace std;
  using namespace valor;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  if(mixing13 < -1) mixing13 = conv::kGlobFitSinsq2Theta13;
  if(mixing23 < -1) mixing23 = conv::kGlobFitSinsq2Theta23;
  if(mixing12 < -1) mixing12 = conv::kGlobFitSinsq2Theta12;
  if(dm23sq < -1)   dm23sq   = conv::kGlobFitDm23sq;
  if(dm12sq < -1)   dm12sq   = conv::kGlobFitDm12sq;
  if(deltacp < -90) deltacp  = conv::kGlobFitDeltaCP;

  //outfile
  string oscparams = util::AsString(conv::kHypOsc3F,
	  mixing12,mixing13,mixing23,deltacp,dm12sq,dm23sq);
  TFile out(Form("osc_probs-%s.root",oscparams.c_str()),"RECREATE");

  const int nbins = 44;
  const double binedges[nbins+1] = {
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,
   0.300,  0.325,  0.350,  0.375,  0.400,  0.425,  0.450,  0.475,  0.500,  0.525,
   0.550,  0.575,  0.600,  0.625,  0.650,  0.675,  0.700,  0.725,  0.750,  0.775,
   0.800,  0.825,  0.850,  0.875,  0.900,  0.925,  0.950,  0.975,
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
   1.500
  };

  //axis, for the binning
  TAxis * axis = new TAxis(nbins,binedges);
  //to store the bin centers
  double binning[nbins];

  NuOsc3Flavour  fNuOsc3Flavour;
  //parameters assumed to be sin^2(2theta)
  fNuOsc3Flavour.SetOscParameters(mixing12,mixing13,mixing23,dm12sq,dm23sq,deltacp);
  
  //store the osc probs
  double p_numu_numu      [nbins];
  double p_numubar_numubar[nbins];
  double p_nue_nue        [nbins];
  double p_nuebar_nuebar  [nbins];
  double p_numu_nue       [nbins];

  // loop over Etrue
  for(int it = 1; it <= nbins; it++) {
    // calculate oscillation probabilities (incl. matter effects) for the current Etrue value
    double Etrue = axis->GetBinCenter(it);
    binning[it-1] = Etrue;

    p_numu_numu      [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuMu,    Etrue);
    p_numubar_numubar[it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuMuBar, Etrue);
    p_nue_nue        [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuE,     conv::kPdgNuE,     Etrue);
    p_nuebar_nuebar  [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuEBar,  conv::kPdgNuEBar,  Etrue);
    p_numu_nue       [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuE,     Etrue);

    /*
    cout << " >>  E = " << Etrue << " GeV" << endl;
    cout << "     P (numu     -> numu   ) = " << p_numu_numu      [it-1] << endl;
    cout << "     P (numubar  -> numubar) = " << p_numubar_numubar[it-1] << endl;
    cout << "     P (nue      -> nue    ) = " << p_nue_nue        [it-1] << endl;
    cout << "     P (nuebar   -> nuebar ) = " << p_nuebar_nuebar  [it-1] << endl;
    cout << "     P (numu     -> nue    ) = " << p_numu_nue       [it-1] << endl;
    */
  }


  // outputs
  TGraph * g_numu_numu       = new TGraph(nbins, binning, p_numu_numu);
  TGraph * g_numubar_numubar = new TGraph(nbins, binning, p_numubar_numubar);
  TGraph * g_nue_nue         = new TGraph(nbins, binning, p_nue_nue);
  TGraph * g_nuebar_nuebar   = new TGraph(nbins, binning, p_nuebar_nuebar);
  TGraph * g_numu_nue        = new TGraph(nbins, binning, p_numu_nue);

  //Draw
  new TCanvas(); g_numu_numu->Draw("LPA");
  new TCanvas(); g_numubar_numubar->Draw("LPA");
  new TCanvas(); g_nue_nue->Draw("LPA");
  new TCanvas(); g_nuebar_nuebar->Draw("LPA");
  new TCanvas(); g_numu_nue->Draw("LPA");

  //Save
  g_numu_numu      ->Write("g_numu_numu");
  g_numubar_numubar->Write("g_numubar_numubar");
  g_nue_nue        ->Write("g_nue_nue");
  g_nuebar_nuebar  ->Write("g_nuebar_nuebar");
  g_numu_nue       ->Write("g_numu_nue");
}
