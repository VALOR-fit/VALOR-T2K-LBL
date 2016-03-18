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
 // double       mixing23   = -99,                            /* */
 // double       dm23sq     = -99,                            /* */
 // double       mixing12   = -99,                            /* */
 // double       mixing13   = -99,                            /* */
 // double       deltacp    = -99,                            /* */
 // double       dm12sq     = -99,                             /* */

 // // sinsq2theta - T2K reactor
 // double       mixing23   = 0.997,                            /* */
 // double       dm23sq     = 0.00251,                            /* */
 // double       mixing12   = 0.8495,                            /* */
 // double       mixing13   = 0.0967, //0.0025, //0.0967,                            /* */
 // double       deltacp    = -1.55,                            /* */
 // double       dm12sq     = 7.5E-5,                             /* */

 // // sinsq2theta T2K only
 // double       mixing23   = 9.97696000000000027e-01,            /* */
 // double       dm23sq     = 0.00251,                            /* */
 // double       mixing12   = 8.49455999999999989e-01,            /* */
 // double       mixing13   = 1.61676640000000038e-01,            /* */
 // double       deltacp    = 0, //1.91,                            /* */
 // double       dm12sq     = 7.5E-5,                             /* */

 // sinsqtheta - T2K reactor
 double       mixing23   = 0.527,                            /* */
 double       dm23sq     = 0.00251,                            /* */
 double       mixing12   = 0.,                            /* */
 double       mixing13   = 0., //0.000001, //0.0248,                            /* */
 double       deltacp    = -1.55,                            /* */
 double       dm12sq     = 7.53E-5,                             /* */

 // // sinsqtheta T2K only
 // double       mixing23   = 0.524,                            /* */
 // double       dm23sq     = 0.00251,                            /* */
 // double       mixing12   = 0.306,                            /* */
 // double       mixing13   = 0.0422,                            /* */
 // double       deltacp    = 1.91,                            /* */
 // double       dm12sq     = 7.5E-5,                             /* */
 
 const char * dmsqscheme     = "abs23",      /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy      = "normal",     /* normal or inverted */
 const char * mixparamchoice = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant = "two"         /* one==kNuOctantOne two==kNuOctantTwo */
 )
{
  using namespace std;
  using namespace valor;

  
  double p_osc[3][3];
  double p_osc_bar[3][3];


  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  Opt::Instance()->SetDefaultParameterValues();

  
  if(Opt::Instance()->gOscMixingParamChoice==conv::kNuMixSinsqTheta){
    if(mixing13 < -1) mixing13 = conv::kGlobFitSinsqTheta13;
    if(mixing23 < -1) mixing23 = conv::kGlobFitSinsqTheta23;
    if(mixing12 < -1) mixing12 = conv::kGlobFitSinsqTheta12;
  }
  else if(Opt::Instance()->gOscMixingParamChoice==conv::kNuMixSinsq2Theta){
    if(mixing13 < -1) mixing13 = conv::kGlobFitSinsq2Theta13;
    if(mixing23 < -1) mixing23 = conv::kGlobFitSinsq2Theta23;
    if(mixing12 < -1) mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  if(dm23sq < -1)   dm23sq   = conv::kGlobFitDm23sq;
  if(dm12sq < -1)   dm12sq   = conv::kGlobFitDm12sq;
  if(deltacp < -90) deltacp  = conv::kGlobFitDeltaCP;

  //outfile
  string oscparams = util::AsString(conv::kHypOsc3F,
	  mixing12,mixing13,mixing23,deltacp,dm12sq,dm23sq);
  TFile out(Form("osc_probs-%s.root",oscparams.c_str()),"RECREATE");

  // const int nbins = 44;
  // const double binedges[nbins+1] = {
  //  0.000,  0.050,  0.100,  0.150,  0.200,  0.250,
  //  0.300,  0.325,  0.350,  0.375,  0.400,  0.425,  0.450,  0.475,  0.500,  0.525,
  //  0.550,  0.575,  0.600,  0.625,  0.650,  0.675,  0.700,  0.725,  0.750,  0.775,
  //  0.800,  0.825,  0.850,  0.875,  0.900,  0.925,  0.950,  0.975,
  //  1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
  //  1.500
  // };

  const int nbins = 100;
  double binedges[nbins+1];
  for(int ibin=0;ibin<(nbins+1);ibin++){
    binedges[ibin] = (double)ibin / ((double)nbins);
    //cout << binedges[ibin] << endl;
  }

  //axis, for the binning
  TAxis * axis = new TAxis(nbins,binedges);
  //to store the bin centers
  double binning[nbins];

  NuOsc3Flavour  fNuOsc3Flavour;
  NuOsc3Flavour  fNuOsc3FlavourBar;
  
  //parameters assumed to be sin^2(2theta)
  if(Opt::Instance()->gOscMixingParamChoice==conv::kNuMixSinsq2Theta)
    fNuOsc3Flavour.SetOscParameters(mixing12,mixing13,mixing23,dm12sq,dm23sq,deltacp);
  else if(Opt::Instance()->gOscMixingParamChoice==conv::kNuMixSinsqTheta)
    fNuOsc3Flavour.SetOscParametersSinSqTAll(mixing12,mixing13,mixing23,dm12sq,dm23sq,deltacp);

  //parameters assumed to be sin^2(2theta)
  if(Opt::Instance()->gOscMixingParamChoice==conv::kNuMixSinsq2Theta)
    fNuOsc3FlavourBar.SetOscParameters(mixing12,mixing13,mixing23,dm12sq,dm23sq,deltacp);
  else if(Opt::Instance()->gOscMixingParamChoice==conv::kNuMixSinsqTheta)
    fNuOsc3FlavourBar.SetOscParametersSinSqTAll(mixing12,mixing13,mixing23,dm12sq,dm23sq,deltacp);
   
  //store the osc probs
  double p_numu_numu      [nbins];
  double p_numubar_numubar[nbins];
  double p_nue_nue        [nbins];
  double p_nuebar_nuebar  [nbins];
  double p_numu_nue       [nbins];
  double p_nue_numu       [nbins];
  double p_numubar_nuebar [nbins];
  double p_nuebar_numubar [nbins];

  double p_numubar_nutaubar [nbins];
  double p_nutaubar_nuebar [nbins];
  double p_numubar_nutaubar_nuebar [nbins];

  double p_nuebar_nutaubar [nbins];
  double p_nutaubar_nuebar [nbins];
  double p_nuebar_nutaubar_nuebar [nbins];


  // loop over Etrue
  for(int it = 1; it <= nbins; it++) {
    // calculate oscillation probabilities (incl. matter effects) for the current Etrue value
    double Etrue = axis->GetBinCenter(it);
    binning[it-1] = Etrue;


    
    
    // fNuOsc3Flavour.CalcMatterOscProbAll(false, Etrue, p_osc    );
    // fNuOsc3FlavourBar.CalcMatterOscProbAll(true , Etrue, p_osc_bar);   
    // p_numu_numu       [it-1] = p_osc    [1][1];
    // p_nue_nue         [it-1] = p_osc    [0][0];
    // p_numu_nue        [it-1] = p_osc    [1][0];
    // p_nue_numu        [it-1] = p_osc    [0][1];
    // p_numubar_numubar [it-1] = p_osc_bar[1][1];
    // p_nuebar_nuebar   [it-1] = p_osc_bar[0][0];
    // p_numubar_nuebar  [it-1] = p_osc_bar[1][0];
    // p_nuebar_numubar  [it-1] = p_osc_bar[0][1];
   

    p_numu_numu      [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuMu,    Etrue);
    p_nue_nue        [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuE,     conv::kPdgNuE,     Etrue);
    p_numu_nue       [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuE,     Etrue);
    p_nue_numu       [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuE,    conv::kPdgNuMu,     Etrue);
    p_numubar_numubar[it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuMuBar, Etrue);
    p_nuebar_nuebar  [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuEBar,  conv::kPdgNuEBar,  Etrue);
    p_numubar_nuebar [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuEBar,  Etrue);
    p_nuebar_numubar [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuEBar, conv::kPdgNuMuBar,  Etrue);

    p_numubar_nutaubar [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuTauBar,  Etrue);
    p_nuebar_nutaubar [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuEBar, conv::kPdgNuTauBar,  Etrue);
    p_nutaubar_nuebar [it-1] = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuTauBar, conv::kPdgNuEBar,  Etrue);
 
    p_numubar_nutaubar_nuebar [it-1] = p_numubar_nutaubar [it-1] * p_nutaubar_nuebar [it-1];
    p_nuebar_nutaubar_nuebar [it-1] = p_nuebar_nutaubar [it-1] * p_nutaubar_nuebar [it-1];
   


    // if( fabs(Etrue-0.6)<0.001 ){
    //   cout << endl;
    //   cout << " >>  E = " << Etrue << " GeV" << endl;
    //   cout << "     P (numu     -> numu   ) = " << p_numu_numu      [it-1] << endl;
    //   cout << "     P (numubar  -> numubar) = " << p_numubar_numubar[it-1] << endl;
    //   cout << "     P (nue      -> nue    ) = " << p_nue_nue        [it-1] << endl;
    //   cout << "     P (nuebar   -> nuebar ) = " << p_nuebar_nuebar  [it-1] << endl;
    //   cout << "     P (numu     -> nue    ) = " << p_numu_nue       [it-1] << endl;
    //   cout << "     P (nue      -> numu   ) = " << p_nue_numu       [it-1] << endl;
    //   cout << "     P (numubar  -> nuebar ) = " << p_numubar_nuebar [it-1] << endl;
    //   cout << "     P (nuebar   -> numubar) = " << p_nuebar_numubar [it-1] << endl;
    //   cout << "     P (numubar  -> nutaubar) = " << p_numubar_nutaubar [it-1] << endl;
    //   cout << "     P (nutaubar -> nuebar) = "   << p_nutaubar_nuebar [it-1] << endl;
    //   cout << "     P (numubar -> nutaubar -> nuebar) = "   << p_numubar_nutaubar_nuebar [it-1] << endl;      
    //   cout << endl;
    // }

    
    cout << "Etrue = " << Etrue 
      << " --> p_numu_nue = " << p_numu_nue [it-1] 
      //<< " --> p_numubar_nuebar = " << p_numubar_nuebar [it-1] 
      //<< " --> p_numubar_nutaubar_nuebar = " << p_numubar_nutaubar_nuebar [it-1] 
	 << endl;

  }


  // outputs
  TGraph * g_numu_numu       = new TGraph(nbins, binning, p_numu_numu);
  TGraph * g_numubar_numubar = new TGraph(nbins, binning, p_numubar_numubar);
  TGraph * g_nue_nue         = new TGraph(nbins, binning, p_nue_nue);
  TGraph * g_nuebar_nuebar   = new TGraph(nbins, binning, p_nuebar_nuebar);
  TGraph * g_numu_nue        = new TGraph(nbins, binning, p_numu_nue);
  TGraph * g_numubar_nuebar  = new TGraph(nbins, binning, p_numubar_nuebar);
  TGraph * g_numubar_nutaubar  = new TGraph(nbins, binning, p_numubar_nutaubar);
  TGraph * g_nuebar_nutaubar  = new TGraph(nbins, binning, p_nuebar_nutaubar);
  TGraph * g_nutaubar_nuebar  = new TGraph(nbins, binning, p_nutaubar_nuebar);
  TGraph * g_numubar_nutaubar_nuebar  = new TGraph(nbins, binning, p_numubar_nutaubar_nuebar);
  TGraph * g_nuebar_nutaubar_nuebar  = new TGraph(nbins, binning, p_nuebar_nutaubar_nuebar);

  //Draw
  TCanvas *c_numu_numu = new TCanvas("c_numu_numu","c_numu_numu"); 
  //g_numu_numu->Draw("LPA");
  g_numu_numu->Draw("LA");
  TCanvas *c_numubar_numubar = new TCanvas("c_numubar_numubar","c_numubar_numubar");
  //g_numubar_numubar->Draw("LPA");
  g_numubar_numubar->Draw("LA");
  TCanvas *c_nue_nue = new TCanvas("c_nue_nue","c_nue_nue"); 
  //g_nue_nue->Draw("LPA");
  g_nue_nue->Draw("LA");
  TCanvas *c_nuebar_nuebar = new TCanvas("c_nuebar_nuebar","c_nuebar_nuebar"); 
  //g_nuebar_nuebar->Draw("LPA");
  g_nuebar_nuebar->Draw("LA");
  TCanvas *c_numu_nue = new TCanvas("c_numu_nue","c_numu_nue"); 
  //g_numu_nue->Draw("LPA");
  g_numu_nue->Draw("LA");
  TCanvas *c_numubar_nuebar = new TCanvas("c_numubar_nuebar","c_numubar_nuebar"); 
  //g_numubar_nuebar->Draw("LPA");
  g_numubar_nuebar->Draw("LA");


  // TCanvas *c_numubar_nutaubar = new TCanvas("c_numubar_nutaubar","c_numubar_nutaubar"); 
  // //g_numubar_nutaubar->Draw("LPA");
  // g_numubar_nutaubar->Draw("LA"); 
  // TCanvas *c_nuebar_nutaubar = new TCanvas("c_nuebar_nutaubar","c_nuebar_nutaubar"); 
  // //g_nuebar_nutaubar->Draw("LPA");
  // g_nuebar_nutaubar->Draw("LA"); 
  // TCanvas *c_nutaubar_nuebar = new TCanvas("c_nutaubar_numubar","c_nutaubar_nuebar"); 
  // //g_nutaubar_nuebar->Draw("LPA");
  // g_nutaubar_nuebar->Draw("LA");
  TCanvas *c_numubar_nutaubar_nuebar = new TCanvas("c_numubar_nutaubar_nuebar","c_numubar_nutaubar_nuebar"); 
  //g_numubar_nutaubar_nuebar->Draw("LPA");
  g_numubar_nutaubar_nuebar->Draw("LA"); 
  TCanvas *c_nuebar_nutaubar_nuebar = new TCanvas("c_nuebar_nutaubar_nuebar","c_nuebar_nutaubar_nuebar"); 
  //g_nuebar_nutaubar_nuebar->Draw("LPA");
  g_nuebar_nutaubar_nuebar->Draw("LA"); 


  //Save
  g_numu_numu      ->Write("g_numu_numu");
  g_numubar_numubar->Write("g_numubar_numubar");
  g_nue_nue        ->Write("g_nue_nue");
  g_nuebar_nuebar  ->Write("g_nuebar_nuebar");
  g_numu_nue       ->Write("g_numu_nue");
  g_numubar_nuebar ->Write("g_numubar_nuebar");
  g_numubar_nutaubar->Write("g_numubar_nutaubar");
  g_nuebar_nutaubar->Write("g_nuebar_nutaubar");
  g_nutaubar_nuebar->Write("g_nutaubar_nuebar");
  g_numubar_nutaubar_nuebar->Write("g_numubar_nutaubar_nuebar");
  g_nuebar_nutaubar_nuebar->Write("g_nuebar_nutaubar_nuebar");

  cout << endl;
  cout << "Used set of oscillation parameters:" << endl;

  if(Opt::Instance()->gOscMixingParamChoice==conv::kNuMixSinsqTheta){
    cout << "sinsqth23 = " << mixing23 << endl;
    cout << "dm23sq = " << dm23sq << endl;
    cout << "sinsqth12 = " << mixing12 << endl;
    cout << "dm12sq = " << dm12sq << endl;
    cout << "sinsqth13 = " << mixing13 << endl;
    cout << "deltacp = " << deltacp << endl;
    cout << endl;
  }
  else if(Opt::Instance()->gOscMixingParamChoice==conv::kNuMixSinsq2Theta){
    cout << "sinsq2th23 = " << mixing23 << endl;
    cout << "dm23sq = " << dm23sq << endl;
    cout << "sinsq2th12 = " << mixing12 << endl;
    cout << "dm12sq = " << dm12sq << endl;
    cout << "sinsq2th13 = " << mixing13 << endl;
    cout << "deltacp = " << deltacp << endl;
    cout << endl;
  }

}
