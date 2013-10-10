#include "TFile.h"
#include "TGraphPainter.h"
#include "TH2D.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <stdio.h>
#include <iostream>

void plot_contour_comparisons
(

)
{
  using std::cout;

  //VaLOR data results
  TFile f_dc_28nob("2012results/DATA/20_2flavour/sensitivity_nobound.root");
  TFile f_dc_28b("2012results/DATA/20_2flavour/sensitivity.root");
  TFile f_dc_20("2012results/DATA/bigger_surface/20_from28/fitresults_3.root");
  TFile f_dc_123("2012results/DATA/bigger_surface/28/fitresults_3.root");
  TFile f_dc_12("2012results/DATA/bigger_surface/28/fitresults_4.root");
  TFile f_dc_3("2012results/DATA/bigger_surface/28/fitresults_5.root");
  TFile f_dc_28_fine("2012results/DATA/fine_grain_28/fitresults_3.root");
  //TFile f_dcall_123("2012results/DATA/2all/fitresults_3.root");
  TFile f_dcall_123("2012results/NEWRESULTS/datafit/dchisq_2all.root");
  TFile f_dcall_12("2012results/DATA/2all/fitresults_4.root");
  TFile f_dcall_3("2012results/DATA/2all/fitresults_5.root");
  TFile f_dcall_28_fine("2012results/DATA/fine_grain_2all/fitresults_3.root");
  //TFile f_fc_123("2012results/FC/fc_contours.root");
  TFile f_fc_123("2012results/NEWRESULTS/FINALFCCONTOURS/fc_contours.root");
  //  TFile f_fcall_123("2012results/FC/fc_contours_2all.root");
  //TFile f_fcall_123("fc_contours_2all.root");
  TFile f_fcall_123("2012results/NEWRESULTS/FINALFCCONTOURS/fc_oldsurf_newres_5E5_nosub.root");
  TFile f_dc_kit("2012results/DATA/KIT/run1-3_constchi2_v2_fine.root");
  //TFile f_fc_kit("2012results/DATA/KIT/run1-3_fc.root");
  //TFile f_fc_kit("of_fc_valor.root");
  TFile f_fc_kit("2012results/DATA/KIT/of_fc_kit.root");
  TFile f_crit("2012results/FC/delta_chisq_critical_Mark2.root");
  TFile f_crit_all("2012results/FC/delta_chisq_critical_2all.root");
  TFile f_crit_kit("2012results/FC/KIT/run1-3_criticalchi2.root");
  TFile f_sens("2012results/SENSITIVITY/t2k_28/sensitivity.root");
  TFile f_dcall_banff13v0("2012results/DATA/banff13v0/combination/fitresults_3.root");
  TFile f_dcall_banff13v3("2012results/DATA/banff13v3/combination/fitresults_3.root");
  TFile f_dcall_octant2only23("2012results/DATA/octant2_23only/combination/fitresults_3.root");
  TFile f_dcall_octant2only23_hack("2012results/DATA/octant2_23only/combination/final.root");
  TFile f_fcall_octant2only23("2012results/FC/fc_contours_2all_octant2.root");

  //Fill in TMarker's with the best fit point (is the same for 2+N, 2+8)
  TMarker bestfit(1.000,0.002442,29);
  bestfit.SetMarkerSize(1.8);
  bestfit.SetMarkerColor(kBlack);
  TMarker bestfit_2all(1.000,0.002443,29);
  bestfit_2all.SetMarkerSize(1.8);
  bestfit_2all.SetMarkerColor(kBlack);
  TMarker bestfit_2all_octant2(0.999,0.002441,29);
  bestfit_2all.SetMarkerSize(1.8);
  bestfit_2all.SetMarkerColor(kBlack);
  TMarker bestfit_2all_from_sinsqtheta(0.999,0.00244,29);
  bestfit_2all_from_sinsqtheta.SetMarkerSize(1.8);
  bestfit_2all_from_sinsqtheta.SetMarkerColor(kBlack);
  TMarker xbestfit_2all_from_sinsqtheta(0.999,2.44,29);
  xbestfit_2all_from_sinsqtheta.SetMarkerSize(1.8);
  xbestfit_2all_from_sinsqtheta.SetMarkerColor(kBlack);
  TMarker bestfit_2all_2nu(1.000,0.002475,29);
  bestfit_2all_2nu.SetMarkerSize(1.8);
  bestfit_2all_2nu.SetMarkerColor(kBlack);

  double linewidth_val = 3;

  //20 const dchisq
  TH2D * h_dc_20 = 0;
  f_dc_20.GetObject("hdchisq",h_dc_20);
  TH2D * h_dc_20_68 = new TH2D(*h_dc_20);
  TH2D * h_dc_20_90 = new TH2D(*h_dc_20);
  h_dc_20_68->SetContour(1);
  h_dc_20_68->SetContourLevel(0,2.30);
  h_dc_20_68->SetLineStyle(2);
  h_dc_20_68->SetLineWidth(linewidth_val);
  h_dc_20_90->SetContour(1);
  h_dc_20_90->SetContourLevel(0,4.61);
  h_dc_20_90->SetLineStyle(1);
  h_dc_20_90->SetLineWidth(linewidth_val);

  //28nob const dchisq
  TH2D * h_dc_28nob_90 = 0, * h_dc_28nob_68 = 0;
  f_dc_28nob.GetObject("dchisq_for_cont1_90",h_dc_28nob_90);
  f_dc_28nob.GetObject("dchisq_for_cont1_68",h_dc_28nob_68);
  h_dc_28nob_68->SetContour(1);
  h_dc_28nob_68->SetContourLevel(0,2.30);
  h_dc_28nob_68->SetLineStyle(2);
  h_dc_28nob_68->SetLineWidth(linewidth_val);
  h_dc_28nob_90->SetContour(1);
  h_dc_28nob_90->SetContourLevel(0,4.61);
  h_dc_28nob_90->SetLineStyle(1);
  h_dc_28nob_90->SetLineWidth(linewidth_val);

  //28b const dchisq
  TH2D * h_dc_28b_90 = 0, * h_dc_28b_68 = 0;
  f_dc_28b.GetObject("dchisq_for_cont1_90",h_dc_28b_90);
  f_dc_28b.GetObject("dchisq_for_cont1_68",h_dc_28b_68);
  h_dc_28b_68->SetContour(1);
  h_dc_28b_68->SetContourLevel(0,2.30);
  h_dc_28b_68->SetLineStyle(2);
  h_dc_28b_68->SetLineWidth(linewidth_val);
  h_dc_28b_90->SetContour(1);
  h_dc_28b_90->SetContourLevel(0,4.61);
  h_dc_28b_90->SetLineStyle(1);
  h_dc_28b_90->SetLineWidth(linewidth_val);

  //28 const dchisq
  TH2D * h_dc_123 = 0;
  f_dc_123.GetObject("hdchisq",h_dc_123);
  TH2D * h_dc_123_68 = new TH2D(*h_dc_123);
  TH2D * h_dc_123_90 = new TH2D(*h_dc_123);
  h_dc_123_68->SetContour(1);
  h_dc_123_68->SetContourLevel(0,2.30);
  h_dc_123_68->SetLineStyle(2);
  h_dc_123_68->SetLineWidth(linewidth_val);
  h_dc_123_90->SetContour(1);
  h_dc_123_90->SetContourLevel(0,4.61);
  h_dc_123_90->SetLineStyle(1);
  h_dc_123_90->SetLineWidth(linewidth_val);

  TH2D * h_dc_12 = 0;
  f_dc_12.GetObject("hdchisq",h_dc_12);
  TH2D * h_dc_12_68 = new TH2D(*h_dc_12);
  TH2D * h_dc_12_90 = new TH2D(*h_dc_12);
  h_dc_12_68->SetContour(1);
  h_dc_12_68->SetContourLevel(0,2.30);
  h_dc_12_68->SetLineStyle(2);
  h_dc_12_68->SetLineWidth(linewidth_val);
  h_dc_12_90->SetContour(1);
  h_dc_12_90->SetContourLevel(0,4.61);
  h_dc_12_90->SetLineStyle(1);
  h_dc_12_90->SetLineWidth(linewidth_val);

  TH2D * h_dc_3 = 0;
  f_dc_3.GetObject("hdchisq",h_dc_3);
  TH2D * h_dc_3_68 = new TH2D(*h_dc_3);
  TH2D * h_dc_3_90 = new TH2D(*h_dc_3);
  h_dc_3_68->SetContour(1);
  h_dc_3_68->SetContourLevel(0,2.30);
  h_dc_3_68->SetLineStyle(2);
  h_dc_3_68->SetLineWidth(linewidth_val);
  h_dc_3_90->SetContour(1);
  h_dc_3_90->SetContourLevel(0,4.61);
  h_dc_3_90->SetLineStyle(1);
  h_dc_3_90->SetLineWidth(linewidth_val);

  //28 const dchisq (fine grain)
  TH2D * h_dc_28_fine = 0;
  f_dc_28_fine.GetObject("hdchisq",h_dc_28_fine);
  TH2D * h_dc_28_fine_68 = new TH2D(*h_dc_28_fine);
  TH2D * h_dc_28_fine_90 = new TH2D(*h_dc_28_fine);
  h_dc_28_fine_68->SetContour(1);
  h_dc_28_fine_68->SetContourLevel(0,2.30);
  h_dc_28_fine_68->SetLineStyle(2);
  h_dc_28_fine_68->SetLineWidth(linewidth_val);
  h_dc_28_fine_90->SetContour(1);
  h_dc_28_fine_90->SetContourLevel(0,4.61);
  h_dc_28_fine_90->SetLineStyle(1);
  h_dc_28_fine_90->SetLineWidth(linewidth_val);

  //28 FC
  TH2D * h_fc_123_68 = 0;
  TH2D * h_fc_123_90 = 0;
  f_fc_123.GetObject("hCL68StatSyst",h_fc_123_68);
  f_fc_123.GetObject("hCL90StatSyst",h_fc_123_90);
  h_fc_123_68->SetContour(1);
  h_fc_123_68->SetContourLevel(0,0.0);
  h_fc_123_68->SetLineStyle(2);
  h_fc_123_68->SetLineWidth(linewidth_val);
  h_fc_123_90->SetContour(1);
  h_fc_123_90->SetContourLevel(0,0.0);
  h_fc_123_90->SetLineStyle(1);
  h_fc_123_90->SetLineWidth(linewidth_val);

  //2all FC
  TH2D * h_fcall_123_68 = 0;
  TH2D * h_fcall_123_90 = 0;
  f_fcall_123.GetObject("hCL68StatSyst",h_fcall_123_68);
  f_fcall_123.GetObject("hCL90StatSyst",h_fcall_123_90);
  h_fcall_123_68->SetContour(1);
  h_fcall_123_68->SetContourLevel(0,0.0);
  h_fcall_123_68->SetLineStyle(2);
  h_fcall_123_68->SetLineWidth(linewidth_val);
  h_fcall_123_90->SetContour(1);
  h_fcall_123_90->SetContourLevel(0,0.0);
  h_fcall_123_90->SetLineStyle(1);
  h_fcall_123_90->SetLineWidth(linewidth_val);

  //2all FC 2nd octant
  TH2D * h_fcall_octant2only23_68 = 0;
  TH2D * h_fcall_octant2only23_90 = 0;
  f_fcall_octant2only23.GetObject("hCL68StatSyst",h_fcall_octant2only23_68);
  f_fcall_octant2only23.GetObject("hCL90StatSyst",h_fcall_octant2only23_90);
  h_fcall_octant2only23_68->SetContour(1);
  h_fcall_octant2only23_68->SetContourLevel(0,0.0);
  h_fcall_octant2only23_68->SetLineStyle(2);
  h_fcall_octant2only23_68->SetLineWidth(linewidth_val);
  h_fcall_octant2only23_90->SetContour(1);
  h_fcall_octant2only23_90->SetContourLevel(0,0.0);
  h_fcall_octant2only23_90->SetLineStyle(1);
  h_fcall_octant2only23_90->SetLineWidth(linewidth_val);

  //2all const dchisq
  TH2D * h_dcall_123 = 0;
  f_dcall_123.GetObject("hdchisq",h_dcall_123);
  TH2D * h_dcall_123_68 = new TH2D(*h_dcall_123);
  TH2D * h_dcall_123_90 = new TH2D(*h_dcall_123);
  h_dcall_123_68->SetContour(1);
  h_dcall_123_68->SetContourLevel(0,2.30);
  h_dcall_123_68->SetLineStyle(2);
  h_dcall_123_68->SetLineWidth(linewidth_val);
  h_dcall_123_68->SetLineWidth(linewidth_val);
  h_dcall_123_90->SetContour(1);
  h_dcall_123_90->SetContourLevel(0,4.61);
  h_dcall_123_90->SetLineStyle(1);
  h_dcall_123_90->SetLineWidth(linewidth_val);
  h_dcall_123_90->SetLineWidth(linewidth_val);

  TH2D * h_dcall_12 = 0;
  f_dcall_12.GetObject("hdchisq",h_dcall_12);
  TH2D * h_dcall_12_68 = new TH2D(*h_dcall_12);
  TH2D * h_dcall_12_90 = new TH2D(*h_dcall_12);
  h_dcall_12_68->SetContour(1);
  h_dcall_12_68->SetContourLevel(0,2.30);
  h_dcall_12_68->SetLineStyle(2);
  h_dcall_12_68->SetLineWidth(linewidth_val);
  h_dcall_12_90->SetContour(1);
  h_dcall_12_90->SetContourLevel(0,4.61);
  h_dcall_12_90->SetLineStyle(1);
  h_dcall_12_90->SetLineWidth(linewidth_val);

  TH2D * h_dcall_3 = 0;
  f_dcall_3.GetObject("hdchisq",h_dcall_3);
  TH2D * h_dcall_3_68 = new TH2D(*h_dcall_3);
  TH2D * h_dcall_3_90 = new TH2D(*h_dcall_3);
  h_dcall_3_68->SetContour(1);
  h_dcall_3_68->SetContourLevel(0,2.30);
  h_dcall_3_68->SetLineStyle(2);
  h_dcall_3_68->SetLineWidth(linewidth_val);
  h_dcall_3_90->SetContour(1);
  h_dcall_3_90->SetContourLevel(0,4.61);
  h_dcall_3_90->SetLineStyle(1);
  h_dcall_3_90->SetLineWidth(linewidth_val);

  //2all with 2013 banff
  TH2D * h_dcall_banff13v0 = 0;
  f_dcall_banff13v0.GetObject("hdchisq",h_dcall_banff13v0);
  TH2D * h_dcall_banff13v0_68 = new TH2D(*h_dcall_banff13v0);
  TH2D * h_dcall_banff13v0_90 = new TH2D(*h_dcall_banff13v0);
  h_dcall_banff13v0_68->SetContour(1);
  h_dcall_banff13v0_68->SetContourLevel(0,2.30);
  h_dcall_banff13v0_68->SetLineStyle(2);
  h_dcall_banff13v0_68->SetLineWidth(linewidth_val);
  h_dcall_banff13v0_90->SetContour(1);
  h_dcall_banff13v0_90->SetContourLevel(0,4.61);
  h_dcall_banff13v0_90->SetLineStyle(1);
  h_dcall_banff13v0_90->SetLineWidth(linewidth_val);

  TH2D * h_dcall_banff13v3 = 0;
  f_dcall_banff13v3.GetObject("hdchisq",h_dcall_banff13v3);
  TH2D * h_dcall_banff13v3_68 = new TH2D(*h_dcall_banff13v3);
  TH2D * h_dcall_banff13v3_90 = new TH2D(*h_dcall_banff13v3);
  h_dcall_banff13v3_68->SetContour(1);
  h_dcall_banff13v3_68->SetContourLevel(0,2.30);
  h_dcall_banff13v3_68->SetLineStyle(2);
  h_dcall_banff13v3_68->SetLineWidth(linewidth_val);
  h_dcall_banff13v3_90->SetContour(1);
  h_dcall_banff13v3_90->SetContourLevel(0,4.61);
  h_dcall_banff13v3_90->SetLineStyle(1);
  h_dcall_banff13v3_90->SetLineWidth(linewidth_val);

  //OCTANT 2
  TH2D * h_dcall_octant2only23 = 0;
  f_dcall_octant2only23.GetObject("hdchisq",h_dcall_octant2only23);
  TH2D * h_dcall_octant2only23_68 = new TH2D(*h_dcall_octant2only23);
  TH2D * h_dcall_octant2only23_90 = new TH2D(*h_dcall_octant2only23);
  h_dcall_octant2only23_68->SetContour(1);
  h_dcall_octant2only23_68->SetContourLevel(0,2.30);
  //h_dcall_octant2only23_68->SetContourLevel(1,2.30);
  h_dcall_octant2only23_68->SetLineStyle(2);
  h_dcall_octant2only23_68->SetLineWidth(linewidth_val);
  h_dcall_octant2only23_90->SetContour(1);
  h_dcall_octant2only23_90->SetContourLevel(0,4.61);
  //h_dcall_octant2only23_90->SetContourLevel(1,4.61);
  h_dcall_octant2only23_90->SetLineStyle(1);
  h_dcall_octant2only23_90->SetLineWidth(linewidth_val);
  cout << h_dcall_octant2only23_90->GetContourLevel(0) << "\t" << h_dcall_octant2only23_90->GetContourLevel(1) << endl;
  cout << h_dcall_octant2only23_68->GetContourLevel(0) << "\t" << h_dcall_octant2only23_68->GetContourLevel(1) << endl;
  /*
  double octant2only23_chisq_bf = 56.03;
  TH2D * h_dcall_octant2only23 = 0;
  f_dcall_octant2only23_hack.GetObject("chisqwithsyst_3",h_dcall_octant2only23);
  TH2D * h_dcall_octant2only23_68 = new TH2D(*h_dcall_octant2only23);
  TH2D * h_dcall_octant2only23_90 = new TH2D(*h_dcall_octant2only23);
  h_dcall_octant2only23_68->SetContour(1);
  h_dcall_octant2only23_68->SetContourLevel(0,2.30+octant2only23_chisq_bf);
  h_dcall_octant2only23_68->SetLineStyle(2);
  h_dcall_octant2only23_68->SetLineWidth(linewidth_val);
  h_dcall_octant2only23_90->SetContour(1);
  h_dcall_octant2only23_90->SetContourLevel(0,4.61+octant2only23_chisq_bf);
  h_dcall_octant2only23_90->SetLineStyle(1);
  h_dcall_octant2only23_90->SetLineWidth(linewidth_val);
  */
  h_dcall_octant2only23_90->SetLineWidth(2);
  h_dcall_octant2only23_68->SetLineWidth(2);

  //KIT DC
  TH2D * h_dc_kit_68 = 0;
  TH2D * h_dc_kit_90 = 0;
  f_dc_kit.GetObject("h_cont_68cl",h_dc_kit_68);
  f_dc_kit.GetObject("h_cont_90cl",h_dc_kit_90);
  //h_dc_kit_68->SetContour(1);
  //h_dc_kit_68->SetContourLevel(0,0.0);
  h_dc_kit_68->SetLineStyle(2);
  h_dc_kit_68->SetLineWidth(linewidth_val);
  //h_dc_kit_90->SetContour(1);
  //h_dc_kit_90->SetContourLevel(0,0.0);
  h_dc_kit_90->SetLineStyle(1);
  h_dc_kit_90->SetLineWidth(linewidth_val);
  
  //KIT FC
  TH2D * h_fc_kit_68 = 0;
  TH2D * h_fc_kit_90 = 0;
  f_fc_kit.GetObject("h_fccont_68",h_fc_kit_68);
  f_fc_kit.GetObject("h_fccont_90",h_fc_kit_90);
  if(!h_fc_kit_68 || !h_fc_kit_90) {
    f_fc_kit.GetObject("hfccontvalor",h_fc_kit_68);
    f_fc_kit.GetObject("hfccontvalor",h_fc_kit_90);
    /*
    TH2D * h_fc_kit_temp = 0;
    f_fc_kit.GetObject("hfccontvalor",h_fc_kit_temp);
    h_fc_kit_68 = new TH2D(*h_fc_kit_temp);
    //h_fc_kit_68->SetContour(1);
    //h_fc_kit_68->SetContourLevel(0,2.30);
    h_fc_kit_90 = new TH2D(*h_fc_kit_temp);
    //h_fc_kit_90->SetContour(1);
    //h_fc_kit_90->SetContourLevel(0,4.61);
    */
  }
  if(!h_fc_kit_68 || !h_fc_kit_90) {
    f_fc_kit.GetObject("hfccontkit68",h_fc_kit_68);
    f_fc_kit.GetObject("hfccontkit90",h_fc_kit_90);
  }
  h_fc_kit_68->SetLineStyle(2);
  h_fc_kit_68->SetLineWidth(linewidth_val);
  h_fc_kit_90->SetLineStyle(1);
  h_fc_kit_90->SetLineWidth(linewidth_val);



  //KIT BF
  /*
    TGraph * kitbf = 0;
    f_dc_kit.GetObject("g_bstosc",kitbf);
    double kitx = 0, kity = 0;
    kitbf->GetPoint(0,kitx,kity);
    TMarker bestfit_kit(kitx,kity,29);
  */
  TMarker bestfit_kit(1.000,0.002444,29);
  bestfit_kit.SetMarkerSize(1.8);
  bestfit_kit.SetMarkerColor(kGreen+2);

  //externals + olds
  TFile fminos("$T2KOSCANA/inputs/contours/MINOS_NumuDisappearance_LogLikelihoodSurface_May2013.root");
  TH2D * ll_minos = 0;
  fminos.GetObject("LogLikelihoodSurface",ll_minos);
  ll_minos->SetContour(1);
  ll_minos->SetContourLevel(0,4.61/2);
  ll_minos->SetLineStyle(1);
  TFile fskle("$T2KOSCANA/inputs/contours/contours_skloe.root");
  TGraph * g_skle = 0;
  fskle.GetObject("sk3f_loe_90cl",g_skle);
  g_skle->SetLineColor(kGreen+2);
  g_skle->SetLineStyle(1);
  TFile fskz("$T2KOSCANA/inputs/contours/contours_skzen.root");
  TGraph * g_skzen = 0;
  fskz.GetObject("sk3f_23vars_overlay_normal_q13_fixed",g_skzen);
  g_skzen->SetLineColor(kBlue);
  g_skzen->SetLineStyle(1);
  TFile fextern("$T2KOSCANA/inputs/contours/contours_nont2k.root");
  TGraph * g_ice = 0;
  fextern.GetObject("icecube20102011_90CL",g_ice);
  g_ice->SetLineColor(kBlue);
  g_ice->SetLineStyle(1);

  TFile fold("$T2KOSCANA/inputs/contours/contours_oldt2k.root");
  TGraph * g_run12_68_statsyst = 0;
  fold.GetObject("t2k2010_run1p2_68CL_statsyst",g_run12_68_statsyst);
  TGraph * g_run12_90_statsyst = 0;
  fold.GetObject("t2k2010_run1p2_90CL_statsyst",g_run12_90_statsyst);
  TGraph * g_run12_68_stat = 0;
  fold.GetObject("t2k2010_run1p2_68CL_stat",g_run12_68_stat);
  TGraph * g_run12_90_stat = 0;
  fold.GetObject("t2k2010_run1p2_90CL_stat",g_run12_90_stat);

  //so we can put the graphs in the legends
  g_skle->SetName("g_skle");
  g_skzen->SetName("g_skzen");
  g_ice->SetName("g_ice");
  g_run12_68_stat->SetName("g_run12_68_stat");
  g_run12_90_stat->SetName("g_run12_90_stat");
  g_run12_68_statsyst->SetName("g_run12_68_statsyst");
  g_run12_90_statsyst->SetName("g_run12_90_statsyst");
  //width
  double linewidth_ext = 1.5;
  g_skle->SetLineWidth(linewidth_ext);
  g_skzen->SetLineWidth(linewidth_ext);
  g_ice->SetLineWidth(linewidth_ext);
  g_run12_68_stat->SetLineWidth(linewidth_ext);
  g_run12_90_stat->SetLineWidth(linewidth_ext);
  g_run12_68_statsyst->SetLineWidth(linewidth_ext);
  g_run12_90_statsyst->SetLineWidth(linewidth_ext);
  ll_minos->SetLineWidth(linewidth_ext);


  // param space
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.16);

  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.6;
  double mixing23_max   = 1.5;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 6E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + TMath::CeilNint ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + TMath::CeilNint ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
  osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  //osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  TLegend * leg = new TLegend(0.17,0.65,0.43,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);

  //colors
  h_dc_20_68->SetLineColor(kBlue+1);
  h_dc_20_90->SetLineColor(kBlue+1);
  h_dc_28_fine_68->SetLineColor(kRed);
  h_dc_28_fine_90->SetLineColor(kRed);
  h_dc_12_68->SetLineColor(kBlue);
  h_dc_12_90->SetLineColor(kBlue);
  h_dc_3_68->SetLineColor(kGreen+2);
  h_dc_3_90->SetLineColor(kGreen+2);
  h_dcall_123_68->SetLineColor(kBlack);
  h_dcall_123_90->SetLineColor(kBlack);
  h_dcall_12_68->SetLineColor(kBlue+3);
  h_dcall_12_90->SetLineColor(kBlue+3);
  h_dcall_3_68->SetLineColor(kGreen+4);
  h_dcall_3_90->SetLineColor(kGreen+4);
  h_fc_123_68->SetLineColor(kBlack);
  h_fc_123_90->SetLineColor(kBlack);
  h_fcall_123_68->SetLineColor(kGreen+2);
  h_fcall_123_90->SetLineColor(kGreen+2);
  g_run12_68_statsyst->SetLineColor(kRed);
  g_run12_90_statsyst->SetLineColor(kRed);
  h_dc_kit_68->SetLineColor(kGreen+2);
  h_dc_kit_90->SetLineColor(kGreen+2);
  h_fc_kit_68->SetLineColor(kGreen+2);
  h_fc_kit_90->SetLineColor(kGreen+2);

  //get rid of the crap at the top of the FC contour
  h_fc_123_68->GetYaxis()->SetRangeUser(0,0.33);
  h_fc_123_90->GetYaxis()->SetRangeUser(0,0.33);

  //separate run markers (with the same color as the dc plots)
  TMarker bestfit123(1.000,0.002442,29);
  bestfit123.SetMarkerSize(1.8);
  bestfit123.SetMarkerColor(kRed);
  TMarker bestfit12(0.986,0.002589,29);
  bestfit12.SetMarkerSize(1.8);
  bestfit12.SetMarkerColor(kBlue);
  TMarker bestfit3(1.000,0.002385,29);
  bestfit3.SetMarkerSize(1.8);
  bestfit3.SetMarkerColor(kGreen+2);
  TMarker bestfit12_2all(0.986,0.002588,29);
  bestfit12_2all.SetMarkerSize(1.8);
  bestfit12_2all.SetMarkerColor(kBlue);
  TMarker bestfit3_2all(1.000,0.002383,29);
  bestfit3_2all.SetMarkerSize(1.8);
  bestfit3_2all.SetMarkerColor(kGreen+2);

  //Run 1-3 const dchisq vs FC (2+8)
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_dc_28_fine_68->SetLineColor(kRed);
  h_dc_28_fine_90->SetLineColor(kRed);
  h_fc_123_68->SetLineColor(kBlack);
  h_fc_123_90->SetLineColor(kBlack);
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_28_fine_68, "const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_fc_123_68, "FC 68% CL","l");
  leg->AddEntry(h_fc_123_90, "FC 90% CL","l");
  leg->Draw();
  bestfit.SetMarkerStyle(29);
  bestfit.SetMarkerColor(kBlack);
  bestfit.Draw();
  c->SaveAs("compare_FC_constchisq.eps");
  c->SaveAs("compare_FC_constchisq.C");

  //Run 1-3 const dchisq vs FC (2+all)
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_dcall_123_68->SetLineColor(kRed);
  h_dcall_123_90->SetLineColor(kRed);
  h_fcall_123_68->SetLineColor(kBlack);
  h_fcall_123_90->SetLineColor(kBlack);
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dcall_123_68, "const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_fcall_123_68, "FC 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "FC 90% CL","l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.Draw();
  c->SaveAs("compare_FC_constchisq_2N.eps");

  //Run 1-3 const dchisq vs FC (2+all 2nd octant)
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_dcall_octant2only23_68->SetLineColor(kRed);
  h_dcall_octant2only23_90->SetLineColor(kRed);
  h_fcall_octant2only23_68->SetLineColor(kBlack);
  h_fcall_octant2only23_90->SetLineColor(kBlack);
  h_dcall_octant2only23_68->Draw("cont3same");
  h_dcall_octant2only23_90->Draw("cont3same");
  h_fcall_octant2only23_68->Draw("cont3same");
  h_fcall_octant2only23_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dcall_octant2only23_68, "const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_octant2only23_90, "const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_fcall_octant2only23_68, "FC 68% CL","l");
  leg->AddEntry(h_fcall_octant2only23_90, "FC 90% CL","l");
  leg->Draw();
  bestfit_2all_octant2.SetMarkerStyle(29);
  bestfit_2all_octant2.SetMarkerColor(kBlack);
  bestfit_2all_octant2.Draw();
  c->SaveAs("compare_FC_constchisq_2NO2.eps");
  c->SaveAs("compare_FC_constchisq_2NO2.pdf");

  //Run 1-3 const dchisq banff 12v8 vs 13v0
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_dcall_123_68->SetLineColor(kRed);
  h_dcall_123_90->SetLineColor(kRed);
  h_dcall_banff13v0_68->SetLineColor(kBlack);
  h_dcall_banff13v0_90->SetLineColor(kBlack);
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  h_dcall_banff13v0_68->Draw("cont3same");
  h_dcall_banff13v0_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dcall_123_68, "BANFF12v8 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "BANFF12v8 90% CL","l");
  leg->AddEntry(h_dcall_banff13v0_68, "BANFF13v0 68% CL","l");
  leg->AddEntry(h_dcall_banff13v0_90, "BANFF13v0 90% CL","l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.Draw();
  c->SaveAs("compare_banff_12v8_13v0.pdf");

  //Run 1-3 const dchisq banff 12v8 vs 13v3
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_dcall_123_68->SetLineColor(kRed);
  h_dcall_123_90->SetLineColor(kRed);
  h_dcall_banff13v3_68->SetLineColor(kBlack);
  h_dcall_banff13v3_90->SetLineColor(kBlack);
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  h_dcall_banff13v3_68->Draw("cont3same");
  h_dcall_banff13v3_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dcall_123_68, "BANFF12v8 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "BANFF12v8 90% CL","l");
  leg->AddEntry(h_dcall_banff13v3_68, "BANFF13v3 68% CL","l");
  leg->AddEntry(h_dcall_banff13v3_90, "BANFF13v3 90% CL","l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.Draw();
  c->SaveAs("compare_banff_12v8_13v3.pdf");

  //Run 1-3 const dchisq octant 1 vs octant 2 (23 only)
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.0033);
  osc_param_space->Draw();  
  h_dcall_123_68->SetLineColor(kBlack);
  h_dcall_123_90->SetLineColor(kBlack);
  h_dcall_octant2only23_68->SetLineColor(kRed);
  h_dcall_octant2only23_90->SetLineColor(kRed);
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  h_dcall_octant2only23_68->Draw("cont3same");
  h_dcall_octant2only23_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dcall_123_68, "Octant 1 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "Octant 1 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dcall_octant2only23_68, "Octant 2 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_octant2only23_90, "Octant 2 const #Delta#chi^{2} 90% CL","l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.Draw();
  bestfit_2all_octant2.SetMarkerStyle(25);
  bestfit_2all_octant2.SetMarkerColor(kRed);
  bestfit_2all_octant2.SetMarkerSize(1.8);
  bestfit_2all_octant2.Draw();
  c->SaveAs("compare_octants.pdf");
  c->SaveAs("compare_octants.eps");

  //Run 1-3 FC octant 1 vs octant 2 (23 only)
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.0033);
  osc_param_space->Draw();  
  h_fcall_123_68->SetLineColor(kBlack);
  h_fcall_123_90->SetLineColor(kBlack);
  h_fcall_octant2only23_68->SetLineColor(kRed);
  h_fcall_octant2only23_90->SetLineColor(kRed);
  h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  h_fcall_octant2only23_68->Draw("cont3same");
  h_fcall_octant2only23_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_fcall_123_68, "Octant 1 FC 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "Octant 1 FC 90% CL","l");
  leg->AddEntry(h_fcall_octant2only23_68, "Octant 2 FC 68% CL","l");
  leg->AddEntry(h_fcall_octant2only23_90, "Octant 2 FC 90% CL","l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.Draw();
  bestfit_2all_octant2.SetMarkerStyle(25);
  bestfit_2all_octant2.SetMarkerColor(kRed);
  bestfit_2all_octant2.SetMarkerSize(1.8);
  bestfit_2all_octant2.Draw();
  c->SaveAs("compare_octants_FC.pdf");
  c->SaveAs("compare_octants_FC.eps");

  //Run 1-3 const dchisq KIT vs VaLOR
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  h_dc_kit_68->Draw("cont3same");
  h_dc_kit_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_28_fine_68, "VaLOR const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "VaLOR const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dc_kit_68, "KIT const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_kit_90, "KIT const #Delta#chi^{2} 90% CL","l");
  leg->Draw();
  //bestfit.Draw();
  c->SaveAs("compare_KIT_constchisq.eps");

  //Run 1-3 FC KIT vs VaLOR
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  h_fc_kit_68->Draw("cont3same");
  h_fc_kit_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_fc_123_68, "VaLOR FC 68% CL","l");
  leg->AddEntry(h_fc_123_90, "VaLOR FC 90% CL","l");
  leg->AddEntry(h_fc_kit_68, "KIT FC 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "KIT FC 90% CL","l");
  leg->Draw();
  //bestfit.Draw();
  c->SaveAs("compare_KIT_FC.eps");
  leg->Clear();
  c->SaveAs("compare_KIT_FC.pdf");

  //Run 1-3 const dchisq vs FC + externals (2+8)
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0035);
  osc_param_space->Draw();  
  //h_dc_28_fine_68->SetLineColor(kRed);
  h_dc_28_fine_90->SetLineColor(kRed);
  //h_fc_123_68->SetLineColor(kBlack);
  h_fc_123_90->SetLineColor(kBlack);
  //h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  //h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  leg->Clear();
  //leg->AddEntry(h_dc_28_fine_68, "const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "T2K const #Delta#chi^{2} 3#nu 90% CL","l");
  //leg->AddEntry(h_fc_123_68, "FC 68% CL","l");
  leg->AddEntry(h_fc_123_90, "T2K FC 3#nu 90% CL","l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->Draw();
  bestfit.SetMarkerStyle(29);
  bestfit.SetMarkerColor(kBlack);
  bestfit.Draw();
  c->SaveAs("compare_FC_constchisq_external.eps");

  //Run 1-3 const dchisq vs FC + externals (2+N)
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0035);
  osc_param_space->Draw();  
  //h_dcall_123_68->SetLineColor(kRed);
  h_dcall_123_90->SetLineColor(kRed);
  //h_fcall_123_68->SetLineColor(kBlack);
  h_fcall_123_90->SetLineColor(kBlack);
  //h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  //h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  leg->Clear();
  //leg->AddEntry(h_dcall_123_68, "const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "T2K const #Delta#chi^{2} 3#nu 90% CL","l");
  //leg->AddEntry(h_fcall_123_68, "FC 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K FC 3#nu 90% CL","l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.Draw();
  c->SaveAs("compare_FC_constchisq_external_2N.eps");

  //Run 1-3 const dchisq vs FC + externals (2+N octant2)
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0035);
  osc_param_space->Draw();  
  //h_dcall_octant2only23_68->SetLineColor(kRed);
  h_dcall_octant2only23_90->SetLineColor(kRed);
  //h_fcall_octant2only23_68->SetLineColor(kBlack);
  h_fcall_octant2only23_90->SetLineColor(kBlack);
  //h_dcall_octant2only23_68->Draw("cont3same");
  h_dcall_octant2only23_90->Draw("cont3same");
  //h_fcall_octant2only23_68->Draw("cont3same");
  h_fcall_octant2only23_90->Draw("cont3same");
  leg->Clear();
  //leg->AddEntry(h_dcall_octant2only23_68, "const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_octant2only23_90, "T2K const #Delta#chi^{2} 3#nu 90% CL","l");
  //leg->AddEntry(h_fcall_octant2only23_68, "FC 68% CL","l");
  leg->AddEntry(h_fcall_octant2only23_90, "T2K FC 3#nu 90% CL","l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->Draw();
  bestfit_2all_octant2.SetMarkerStyle(29);
  bestfit_2all_octant2.SetMarkerColor(kBlack);
  bestfit_2all_octant2.Draw();
  c->SaveAs("compare_FC_constchisq_external_2NO2.eps");

  //Run 1-3 FC octant 1 vs octant2 + externals (2+N)
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0035);
  osc_param_space->Draw();  
  h_fcall_octant2only23_90->SetLineColor(kRed);
  h_fcall_123_90->SetLineColor(kBlack);
  h_fcall_octant2only23_90->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_fcall_123_90, "T2K FC 3#nu (octant 1) 90% CL","l");
  leg->AddEntry(h_fcall_octant2only23_90, "T2K FC 3#nu (octant 2) 90% CL","l");
  ll_minos->SetLineColor(kBlack);
  g_skzen->SetLineColor(kBlue);
  g_skle->SetLineColor(kGreen+2);
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.Draw();
  bestfit_2all_octant2.SetMarkerStyle(25);
  bestfit_2all_octant2.SetMarkerColor(kRed);
  bestfit_2all_octant2.Draw();
  c->SaveAs("compare_FC_octants_external.eps");
  c->SaveAs("compare_FC_octants_external.pdf");


  //run1+2 2012 const dchisq vs 2010 FC
  osc_param_space->GetXaxis()->SetRangeUser(0.75,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0035);
  osc_param_space->Draw();
  h_dc_12_68->SetLineColor(kBlack);
  h_dc_12_90->SetLineColor(kBlack);
  h_dc_12_68->Draw("cont3same");
  h_dc_12_90->Draw("cont3same");
  g_run12_68_statsyst->Draw("l");
  g_run12_90_statsyst->Draw("l");
  leg->Clear();
  leg->AddEntry(h_dc_12_68, "2012 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_12_90, "2012 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry("g_run12_68_statsyst", "2010 FC 68% CL", "l");
  leg->AddEntry("g_run12_90_statsyst", "2010 FC 90% CL", "l");
  leg->Draw();
  bestfit12.SetMarkerStyle(29);
  bestfit12.SetMarkerColor(kBlack);
  bestfit12.Draw();
  TMarker bestfit12_2010(0.970,0.002600,29);
  bestfit123.SetMarkerSize(1.8);
  bestfit123.SetMarkerColor(kRed);
  c->SaveAs("compare_run12_2012_2010.eps");

  //Run 1-3 28, 2all const dchisq
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->Draw();  
  h_dc_28_fine_68->SetLineColor(kBlack);
  h_dc_28_fine_90->SetLineColor(kBlack);
  h_dcall_123_68->SetLineColor(kRed);
  h_dcall_123_90->SetLineColor(kRed);
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_28_fine_68, "2+8 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "2+8 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dcall_123_68, "2+N const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "2+N const #Delta#chi^{2} 90% CL","l");
  leg->Draw();
  bestfit.SetMarkerStyle(29);
  bestfit.SetMarkerColor(kBlack);
  bestfit.Draw();
  bestfit_2all.SetMarkerColor(kRed);
  bestfit_2all.SetMarkerStyle(24);
  bestfit_2all.Draw();
  c->SaveAs("compare_2all_28.eps");

  //Run 1-3 28, 2all const dchisq with 2+0
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->Draw();
  h_dc_20_68->SetLineColor(kBlue);
  h_dc_20_90->SetLineColor(kBlue);
  h_dc_28_fine_68->SetLineColor(kBlack);
  h_dc_28_fine_90->SetLineColor(kBlack);
  h_dcall_123_68->SetLineColor(kRed);
  h_dcall_123_90->SetLineColor(kRed);
  h_dc_20_68->Draw("cont3same");
  h_dc_20_90->Draw("cont3same");
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_20_68, "2+0 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_20_90, "2+0 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dc_28_fine_68, "2+8 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "2+8 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dcall_123_68, "2+N const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "2+N const #Delta#chi^{2} 90% CL","l");
  leg->Draw();
  bestfit.SetMarkerStyle(29);
  bestfit.SetMarkerColor(kBlack);
  bestfit.Draw();
  bestfit_2all.SetMarkerStyle(24);
  bestfit_2all.SetMarkerColor(kRed);
  bestfit_2all.Draw();
  TMarker bestfit20(1.000,0.002457,25);
  bestfit20.SetMarkerSize(1.8);
  bestfit20.SetMarkerColor(kBlue);
  bestfit20.Draw();
  c->SaveAs("compare_2all_28_20.eps");

  //compare 2+8 with/without boundary
  osc_param_space->GetXaxis()->SetRangeUser(0.8,1.2);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0033);
  osc_param_space->Draw();
  h_dc_28nob_68->SetLineColor(kRed);
  h_dc_28nob_90->SetLineColor(kRed);
  h_dc_28nob_68->Draw("cont3same");
  h_dc_28nob_90->Draw("cont3same");
  TMarker bestfit28nob(1.025,0.00247,29);
  bestfit28nob.SetMarkerSize(1.8);
  bestfit28nob.SetMarkerColor(kRed);
  bestfit28nob.Draw();
  c->SaveAs("compare_28nob_only.pdf");
  h_dc_28b_68->SetLineColor(kBlack);
  h_dc_28b_90->SetLineColor(kBlack);
  h_dc_28b_68->GetXaxis()->SetRangeUser(0,1);
  h_dc_28b_90->GetXaxis()->SetRangeUser(0,1);
  h_dc_28b_68->Draw("cont3same");
  h_dc_28b_90->Draw("cont3same");
  TMarker bestfit28b(1.000,0.00247,29);
  bestfit28b.SetMarkerSize(1.8);
  bestfit28b.SetMarkerColor(kBlack);
  bestfit28b.Draw();
  leg->Clear();
  leg->AddEntry(h_dc_28b_68, "2#nu w/ boundary 68% CL","l");
  leg->AddEntry(h_dc_28b_90, "2#nu w/ boundary 90% CL","l");
  leg->AddEntry(h_dc_28nob_68, "2#nu w/o boundary 68% CL","l");
  leg->AddEntry(h_dc_28nob_90, "2#nu w/o boundary 90% CL","l");
  leg->Draw();
  c->SaveAs("compare_28_boundaries.pdf");


  //Run 1-3 28, 2all FC
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->Draw();  
  h_fc_123_68->SetLineColor(kBlack);
  h_fc_123_90->SetLineColor(kBlack);
  h_fcall_123_68->SetLineColor(kRed);
  h_fcall_123_90->SetLineColor(kRed);
  h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_fc_123_68, "2+8 FC 68% CL","l");
  leg->AddEntry(h_fc_123_90, "2+8 FC 90% CL","l");
  leg->AddEntry(h_fcall_123_68, "2+N FC 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "2+N FC 90% CL","l");
  leg->Draw();
  bestfit.SetMarkerStyle(29);
  bestfit.SetMarkerColor(kBlack);
  bestfit.Draw();
  bestfit_2all.SetMarkerStyle(24);
  bestfit_2all.SetMarkerColor(kRed);
  bestfit_2all.Draw();
  c->SaveAs("compare_2all_28_FC.eps");

  //Run 1-3 2all FC with 1D limits
  osc_param_space->GetXaxis()->SetRangeUser(0.92,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.003);
  osc_param_space->Draw();  
  h_fcall_123_68->SetLineColor(kBlack);
  h_fcall_123_90->SetLineColor(kBlack);
  h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  TLine l_limit_mix    (0.95,0.00200,0.95,0.00300);
  TLine l_limit_dm_low (0.92,0.00229,1.00,0.00229);
  TLine l_limit_dm_high(0.92,0.00261,1.00,0.00261);
  l_limit_mix.SetLineColor(kRed);
  l_limit_dm_low.SetLineColor(kRed);
  l_limit_dm_high.SetLineColor(kRed);
  l_limit_mix.Draw();
  l_limit_dm_low.Draw();
  l_limit_dm_high.Draw();
  c->SaveAs("2D_contour_with_1D_limits.pdf");


  //compare 2+8, 2+N const dchisq surfaces
  TH2D * h_dc_dchisqratio1_282N = new TH2D(*h_dcall_123);
  TH2D * h_dc_dchisqratio2_282N = new TH2D(*h_dcall_123);
  for(int i = 1; i <= h_dcall_123->GetNbinsX(); i++) {
    int ibin = h_dc_123->GetXaxis()->FindBin(h_dcall_123->GetXaxis()->GetBinCenter(i));
    for(int j = 1; j <= h_dcall_123->GetNbinsY(); j++) {
      int jbin = h_dc_123->GetYaxis()->FindBin(h_dcall_123->GetYaxis()->GetBinCenter(j));
      double n28 = h_dc_123->GetBinContent(ibin, jbin);
      double n2n = h_dcall_123->GetBinContent(i,j);
      double frac1 = (n28-n2n)/n2n;
      double frac2 = (n2n-n28)/((n2n+n28)/2);
      //cout << i << "\t" << j << "\t" << n28 << "\t" << n2n << "\t" << frac1 << "\t" << frac2 << endl;
      h_dc_dchisqratio1_282N->SetBinContent(i,j, frac1);
      h_dc_dchisqratio2_282N->SetBinContent(i,j, frac2);
    }
  }
  osc_param_space->GetXaxis()->SetRangeUser(0.8,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.003);
  osc_param_space->Draw();  
  ll_minos->SetLineColor(kBlack);
  h_dc_dchisqratio1_282N->SetContour(255);
  h_dc_dchisqratio2_282N->SetContour(255);
  c->SetRightMargin(0.13);
  h_dc_dchisqratio1_282N->Draw("colzsame");
  ll_minos->Draw("cont3same");
  c->SaveAs("compare_2all_28_dchisqdiff_8MnOn.eps");
  osc_param_space->Draw();  
  h_dc_dchisqratio2_282N->Draw("colzsame");
  ll_minos->Draw("cont3same");
  c->SaveAs("compare_2all_28_dchisqdiff_nM8OnP8O2.eps");
  //draw again, with t2k contour
  Color_t fccol = h_fc_123_90->GetLineColor();
  h_fc_123_90->SetLineColor(kBlack);
  osc_param_space->Draw();
  h_dc_dchisqratio1_282N->Draw("colzsame");
  h_fc_123_90->Draw("cont3same");
  c->SaveAs("compare_2all_28_dchisqdiff_8MnOn_t2kcont.eps");
  osc_param_space->Draw();  
  h_dc_dchisqratio2_282N->Draw("colzsame");
  h_fc_123_90->Draw("cont3same");
  c->SaveAs("compare_2all_28_dchisqdiff_nM8OnP8O2_t2kcont.eps");
  //zoom
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.0028);
  h_dc_dchisqratio1_282N->GetXaxis()->SetRangeUser(0.85,1.0);
  h_dc_dchisqratio1_282N->GetYaxis()->SetRangeUser(0.002,0.0028);
  h_dc_dchisqratio2_282N->GetXaxis()->SetRangeUser(0.85,1.0);
  h_dc_dchisqratio2_282N->GetYaxis()->SetRangeUser(0.002,0.0028);
  osc_param_space->Draw();  
  h_dc_dchisqratio1_282N->Draw("colzsame");
  ll_minos->Draw("cont3same");
  c->SaveAs("compare_2all_28_dchisqdiff_8MnOn_zoom.eps");
  osc_param_space->Draw();  
  h_dc_dchisqratio2_282N->Draw("colzsame");
  ll_minos->Draw("cont3same");
  c->SaveAs("compare_2all_28_dchisqdiff_nM8OnP8O2_zoom.eps");
  //zoom into t2k contour region
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.003);
  h_dc_dchisqratio1_282N->GetXaxis()->SetRangeUser(0.9,1.0);
  h_dc_dchisqratio1_282N->GetYaxis()->SetRangeUser(0.002,0.003);
  h_dc_dchisqratio2_282N->GetXaxis()->SetRangeUser(0.9,1.0);
  h_dc_dchisqratio2_282N->GetYaxis()->SetRangeUser(0.002,0.003);
  //draw again, with t2k contour
  osc_param_space->Draw();
  h_dc_dchisqratio1_282N->Draw("colzsame");
  h_fc_123_90->Draw("cont3same");
  c->SaveAs("compare_2all_28_dchisqdiff_8MnOn_zoom_t2kcont.eps");
  osc_param_space->Draw();  
  h_dc_dchisqratio2_282N->Draw("colzsame");
  h_fc_123_90->Draw("cont3same");
  c->SaveAs("compare_2all_28_dchisqdiff_nM8OnP8O2_zoom_t2kcont.eps");
  h_fc_123_90->SetLineColor(fccol);

  //compare dchisq surfaces for different octants
  TH2D * h_dc_dchisqratio1_O1O2 = new TH2D(*h_dcall_octant2only23);
  TH2D * h_dc_dchisqratio2_O1O2 = new TH2D(*h_dcall_octant2only23);
  for(int i = 1; i <= h_dcall_octant2only23->GetNbinsX(); i++) {
    int ibin = h_dcall_123->GetXaxis()->FindBin(h_dcall_octant2only23->GetXaxis()->GetBinCenter(i));
    for(int j = 1; j <= h_dcall_octant2only23->GetNbinsY(); j++) {
      int jbin = h_dcall_123->GetYaxis()->FindBin(h_dcall_octant2only23->GetYaxis()->GetBinCenter(j));
      double no1 = h_dcall_123->GetBinContent(ibin, jbin);
      double no2 = h_dcall_octant2only23->GetBinContent(i,j);
      double frac1 = (no2-no1)/no1;
      double frac2 = (no1-no2)/((no1+no2)/2);
      /*
      cout << i << "\t" << j << "\t" << no1 << "\t" << no2 << "\t" << frac1 << "\t" << frac2 << endl;
      cout << "\t" << ibin << "," << jbin << "\t" 
	   << h_dcall_123->GetXaxis()->GetBinCenter(ibin)        << "," << h_dcall_123->GetYaxis()->GetBinCenter(jbin) << endl;
      cout << "\t" << i    << "," << j    << "\t"
	   << h_dcall_octant2only23->GetXaxis()->GetBinCenter(i) << "," << h_dcall_octant2only23->GetYaxis()->GetBinCenter(j) << endl;
      */
      h_dc_dchisqratio1_O1O2->SetBinContent(i,j, frac1);
      h_dc_dchisqratio2_O1O2->SetBinContent(i,j, frac2);
      /*
	if(TMath::Abs(h_dcall_octant2only23->GetXaxis()->GetBinCenter(i) - 0.95) < 1E-8
	&& TMath::Abs(h_dcall_octant2only23->GetYaxis()->GetBinCenter(j) - 0.0025) < 1E-8) {
	cout << endl << endl << no1 << "\t" << no2 << "\t" << frac2 << endl << endl;
	exit(1);
      }
      */
    }
  }
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.003);
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0021,0.0029);
  osc_param_space->Draw();  
  ll_minos->SetLineColor(kBlack);
  h_dc_dchisqratio1_O1O2->SetContour(255);
  h_dc_dchisqratio2_O1O2->SetContour(255);
  c->SetRightMargin(0.13);
  h_dc_dchisqratio1_O1O2->Draw("colzsame");
  ll_minos->Draw("cont3same");
  c->SaveAs("compare_octants_dchisqdiff_8MnOn.eps");
  osc_param_space->Draw();  
  h_dc_dchisqratio2_O1O2->Draw("colzsame");
  ll_minos->Draw("cont3same");
  c->SaveAs("compare_octants_dchisqdiff_nM8OnP8O2.eps");
  //draw again, with t2k contour
  Color_t fccol = h_fc_123_90->GetLineColor();
  h_fc_123_90->SetLineColor(kBlack);
  osc_param_space->Draw();
  h_dc_dchisqratio1_O1O2->Draw("colzsame");
  h_fc_123_90->Draw("cont3same");
  c->SaveAs("compare_octants_dchisqdiff_8MnOn_t2kcont.eps");
  osc_param_space->Draw();  
  h_dc_dchisqratio2_O1O2->Draw("colzsame");
  h_fc_123_90->Draw("cont3same");
  c->SaveAs("compare_octants_dchisqdiff_nM8OnP8O2_t2kcont.eps");
  //zoom
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.0028);
  h_dc_dchisqratio1_O1O2->GetXaxis()->SetRangeUser(0.85,1.0);
  h_dc_dchisqratio1_O1O2->GetYaxis()->SetRangeUser(0.002,0.0028);
  h_dc_dchisqratio2_O1O2->GetXaxis()->SetRangeUser(0.85,1.0);
  h_dc_dchisqratio2_O1O2->GetYaxis()->SetRangeUser(0.002,0.0028);
  osc_param_space->Draw();  
  h_dc_dchisqratio1_O1O2->Draw("colzsame");
  ll_minos->Draw("cont3same");
  c->SaveAs("compare_octants_dchisqdiff_8MnOn_zoom.eps");
  osc_param_space->Draw();  
  h_dc_dchisqratio2_O1O2->Draw("colzsame");
  ll_minos->Draw("cont3same");
  c->SaveAs("compare_octants_dchisqdiff_nM8OnP8O2_zoom.eps");
  //zoom into t2k contour region
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.003);
  h_dc_dchisqratio1_O1O2->GetXaxis()->SetRangeUser(0.9,1.0);
  h_dc_dchisqratio1_O1O2->GetYaxis()->SetRangeUser(0.002,0.003);
  h_dc_dchisqratio2_O1O2->GetXaxis()->SetRangeUser(0.9,1.0);
  h_dc_dchisqratio2_O1O2->GetYaxis()->SetRangeUser(0.002,0.003);
  //draw again, with t2k contour
  osc_param_space->Draw();
  h_dc_dchisqratio1_O1O2->Draw("colzsame");
  h_fc_123_90->Draw("cont3same");
  c->SaveAs("compare_octants_dchisqdiff_8MnOn_zoom_t2kcont.eps");
  osc_param_space->Draw();  
  h_dc_dchisqratio2_O1O2->Draw("colzsame");
  h_fc_123_90->Draw("cont3same");
  c->SaveAs("compare_octants_dchisqdiff_nM8OnP8O2_zoom_t2kcont.eps");
  h_fc_123_90->SetLineColor(fccol);

  //compare KIT/VaLOR dchisq_crit surfaces
  TH2D * h_crit_68 = 0;
  f_crit.GetObject("hDeltaChi2Surface68",h_crit_68);
  TH2D * h_crit_68_kit = 0;
  f_crit_kit.GetObject("hcchi2_68",h_crit_68_kit);
  TH2D * h_crit_ratio = new TH2D(*h_crit_68_kit);
  h_crit_ratio->Clear();
  for(int i = 1; i <= h_crit_ratio->GetNbinsX(); i++) {
    int ibin = h_crit_68->GetXaxis()->FindBin(h_crit_ratio->GetXaxis()->GetBinCenter(i));
    //cout << i << endl;
    for(int j = 1; j <= h_crit_ratio->GetNbinsY(); j++) {
      h_crit_ratio->SetBinContent(i,j,0);
      int jbin = h_crit_68->GetYaxis()->FindBin(h_crit_ratio->GetYaxis()->GetBinCenter(j));
      double val = h_crit_68->GetBinContent(ibin, jbin);
      double kit = h_crit_68_kit->GetBinContent(i,j);
      double frac = 2*(val-kit)/(kit+val);
      //cout << i << "\t" << j << "\t" << val << "\t" << kit << "\t" << frac << endl;
      h_crit_ratio->SetBinContent(i,j,frac);
    }
  }
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.003);
  osc_param_space->Draw("");
  h_crit_ratio->SetContour(255);
  h_crit_ratio->Draw("colzsame");
  c->SaveAs("compare_fc_crit.eps");

  
  h_crit_68->GetZaxis()->SetRangeUser(1,3);
  h_crit_68_kit->GetZaxis()->SetRangeUser(1,3);

  osc_param_space->Draw("");
  h_crit_68->Draw("colzsame");
  c->SaveAs("valcrit.eps");
  osc_param_space->Draw("");
  h_crit_68_kit->Draw("colzsame");
  c->SaveAs("kitcrit.eps");
  

  //compare 2+8/2+N dchisq_crit surfaces
  TH2D * h_crit_68_all = 0;
  f_crit_all.GetObject("hDeltaChi2Surface68",h_crit_68_all);
  delete h_crit_ratio;
  h_crit_ratio = new TH2D(*h_crit_68_all);
  h_crit_ratio->Clear();
  for(int i = 1; i <= h_crit_ratio->GetNbinsX(); i++) {
    int ibin = h_crit_68->GetXaxis()->FindBin(h_crit_ratio->GetXaxis()->GetBinCenter(i));
    //cout << i << endl;
    for(int j = 1; j <= h_crit_ratio->GetNbinsY(); j++) {
      h_crit_ratio->SetBinContent(i,j,0);
      int jbin = h_crit_68->GetYaxis()->FindBin(h_crit_ratio->GetYaxis()->GetBinCenter(j));
      double val = h_crit_68->GetBinContent(ibin, jbin);
      double all = h_crit_68_all->GetBinContent(i,j);
      double frac = 2*(val-all)/(all+val);
      if(val < 1E-8 && all < 1E-8)
	frac = 0;
      //cout << i << "\t" << j << "\t" << val << "\t" << all << "\t" << frac << endl;
      h_crit_ratio->SetBinContent(i,j,frac);
    }
  }
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.003);
  osc_param_space->Draw("");
  h_crit_ratio->SetContour(255);
  h_crit_ratio->Draw("colzsame");
  c->SaveAs("compare_fc_crit_28_2N.eps");

  h_crit_68_all->GetZaxis()->SetRangeUser(1,3);
  osc_param_space->Draw("");
  h_crit_68_all->Draw("colzsame");
  c->SaveAs("allcrit.eps");
  

  c->SetRightMargin(0.05);
  //compare T2K result with expected sensitivity (const dchisq)
  TH2D * h_sens_68 = 0;
  f_sens.GetObject("dchisq_for_cont1_68",h_sens_68);
  TH2D * h_sens_90 = 0;
  f_sens.GetObject("dchisq_for_cont1_90",h_sens_90);
  h_sens_68->SetLineStyle(2);
  h_sens_68->SetLineWidth(linewidth_val);
  h_sens_68->SetLineColor(kBlue-3);
  h_sens_90->SetLineColor(kBlue-3);
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();
  h_dc_28_fine_68->SetLineColor(kBlack);
  h_dc_28_fine_90->SetLineColor(kBlack);
  //h_fc_123_68->SetLineColor(kBlue);
  //h_fc_123_90->SetLineColor(kBlue);
  h_sens_68->SetLineColor(kRed);
  h_sens_90->SetLineColor(kRed);
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  //h_fc_123_68->Draw("cont3same");
  //h_fc_123_90->Draw("cont3same");
  h_sens_68->Draw("cont3same");
  h_sens_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_28_fine_68, "Run 1-3 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "Run 1-3 const #Delta#chi^{2} 90% CL","l");
  //leg->AddEntry(h_fc_123_68, "FC 68% CL","l");
  //leg->AddEntry(h_fc_123_90, "FC 90% CL","l");
  leg->AddEntry(h_sens_68, "Sensitivity const #Delta#chi^{2} 68% CL", "l");
  leg->AddEntry(h_sens_90, "Sensitivity const #Delta#chi^{2} 90% CL", "l");
  leg->Draw();
  //bestfit.Draw();
  c->SaveAs("compare_constchisq_sensitivity.eps");
  

  //Run 1-3 2+8 FC only
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->Draw();  
  h_fc_123_68->SetLineColor(kBlack);
  h_fc_123_90->SetLineColor(kBlack);
  h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_fc_123_68, "2+8 FC 68% CL","l");
  leg->AddEntry(h_fc_123_90, "2+8 FC 90% CL","l");
  leg->Draw();
  bestfit.SetMarkerStyle(29);
  bestfit.SetMarkerColor(kBlack);
  bestfit.Draw();
  c->SaveAs("28_FC.eps");

  //Run 1-3 2+N FC only
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->Draw();  
  h_fcall_123_68->SetLineColor(kBlack);
  h_fcall_123_90->SetLineColor(kBlack);
  h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_fcall_123_68, "2+N FC 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "2+N FC 90% CL","l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.Draw();
  c->SaveAs("2N_FC.eps");

  //Run 1-3 2+N FC only 2nd octant
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->Draw();  
  h_fcall_octant2only23_68->SetLineColor(kBlack);
  h_fcall_octant2only23_90->SetLineColor(kBlack);
  h_fcall_octant2only23_68->Draw("cont3same");
  h_fcall_octant2only23_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_fcall_octant2only23_68, "2+N FC 68% CL","l");
  leg->AddEntry(h_fcall_octant2only23_90, "2+N FC 90% CL","l");
  leg->Draw();
  bestfit_2all_octant2.SetMarkerStyle(29);
  bestfit_2all_octant2.SetMarkerColor(kBlack);
  bestfit_2all_octant2.Draw();
  c->SaveAs("2NO2_FC.eps");

  //runs 123,12,3 constdchisq comparisons (2+8)
  osc_param_space->GetXaxis()->SetRangeUser(0.75,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0035);
  osc_param_space->Draw();
  h_dc_123_68->SetLineColor(kBlack);
  h_dc_123_90->SetLineColor(kBlack);
  h_dc_12_68->SetLineColor(kBlue);
  h_dc_12_90->SetLineColor(kBlue);
  h_dc_3_68->SetLineColor(kRed);
  h_dc_3_90->SetLineColor(kRed);
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  h_dc_12_68->Draw("cont3same");
  h_dc_12_90->Draw("cont3same");
  h_dc_3_68->Draw("cont3same");
  h_dc_3_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_28_fine_68, "Run1+2+3 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "Run1+2+3 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dc_12_68, "Run1+2 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_12_90, "Run1+2 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dc_3_68, "Run3 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_3_90, "Run3 const #Delta#chi^{2} 90% CL","l");
  leg->Draw();
  bestfit123.SetMarkerStyle(29);
  bestfit123.SetMarkerColor(kBlack);
  bestfit123.Draw();
  bestfit12.SetMarkerStyle(29);
  bestfit12.SetMarkerColor(kBlue);
  bestfit12.Draw();
  bestfit3.SetMarkerStyle(29);
  bestfit3.SetMarkerColor(kRed);
  bestfit3.Draw();
  c->SaveAs("compare_runs.eps");

  //runs 123,12,3 constdchisq comparisons (2+N)
  osc_param_space->GetXaxis()->SetRangeUser(0.75,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0035);
  osc_param_space->Draw();
  h_dcall_123_68->SetLineColor(kBlack);
  h_dcall_123_90->SetLineColor(kBlack);
  h_dcall_12_68->SetLineColor(kBlue);
  h_dcall_12_90->SetLineColor(kBlue);
  h_dcall_3_68->SetLineColor(kRed);
  h_dcall_3_90->SetLineColor(kRed);
  bestfit_2all.SetMarkerColor(kRed);
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  h_dcall_12_68->Draw("cont3same");
  h_dcall_12_90->Draw("cont3same");
  h_dcall_3_68->Draw("cont3same");
  h_dcall_3_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dcall_123_68, "Run1+2+3 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "Run1+2+3 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dcall_12_68, "Run1+2 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_12_90, "Run1+2 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dcall_3_68, "Run3 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_3_90, "Run3 const #Delta#chi^{2} 90% CL","l");
  leg->Draw();
  bestfit_2all.SetMarkerStyle(29);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit12_2all.SetMarkerStyle(29);
  bestfit12_2all.SetMarkerColor(kBlue);
  bestfit3_2all.SetMarkerStyle(29);
  bestfit3_2all.SetMarkerColor(kRed);
  bestfit_2all.Draw();
  bestfit12_2all.Draw();
  bestfit3_2all.Draw();
  c->SaveAs("compare_runs_2N.eps");


  // different format / legends for PAC
  /*
  //Run 1-3 28, 2all const dchisq with 2+0
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0015,0.0033);
  osc_param_space->Draw();
  h_dc_20_68->Draw("cont3same");
  h_dc_20_90->Draw("cont3same");
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_20_68, "2+0 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_20_90, "2+0 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dc_28_fine_68, "2+8 const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "2+8 const #Delta#chi^{2} 90% CL","l");
  leg->AddEntry(h_dcall_123_68, "2+N const #Delta#chi^{2} 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "2+N const #Delta#chi^{2} 90% CL","l");
  leg->Draw();
  bestfit.Draw();
  TMarker bestfit20(1.000,0.002457,29);
  bestfit20.SetMarkerSize(1.8);
  bestfit20.SetMarkerColor(kBlue);
  bestfit20.Draw();
  c->SaveAs("compare_contours_constdchisq_2all_28_20.jpg");
  c->SaveAs("compare_contours_constdchisq_2all_28_20.png");
  c->SaveAs("compare_contours_constdchisq_2all_28_20.eps");
  c->SaveAs("compare_contours_constdchisq_2all_28_20.pdf");

  //change colors for kit comparisons
  ll_minos->SetLineColor(kGreen+2);
  g_skle->SetLineColor(kBlue);
  g_skzen->SetLineColor(kCyan+1);
  g_run12_90_statsyst->SetLineColor(kMagenta+1);
  h_dc_28_fine_68->SetLineColor(kBlack);
  h_dc_28_fine_90->SetLineColor(kBlack);
  h_dc_kit_68->SetLineColor(kRed);
  h_dc_kit_90->SetLineColor(kRed);
  h_fc_123_68->SetLineColor(kBlack);
  h_fc_123_90->SetLineColor(kBlack);
  h_fc_kit_68->SetLineColor(kRed);
  h_fc_kit_90->SetLineColor(kRed);
  bestfit.SetMarkerColor(kBlack);
  bestfit_kit.SetMarkerColor(kRed);

  //Run 1-3 const dchisq KIT vs VaLOR
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  h_dc_kit_68->Draw("cont3same");
  h_dc_kit_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_dc_28_fine_68, "A const #Delta#chi^{2} 3#nu 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "A const #Delta#chi^{2} 3#nu 90% CL","l");
  leg->AddEntry(h_dc_kit_68, "B const #Delta#chi^{2} 3#nu 68% CL","l");
  leg->AddEntry(h_dc_kit_90, "B const #Delta#chi^{2} 3#nu 90% CL","l");
  leg->Draw();
  bestfit.Draw();
  bestfit_kit.Draw();
  c->SaveAs("compare_contours_AB_constdchisq.png");
  c->SaveAs("compare_contours_AB_constdchisq.jpg");
  c->SaveAs("compare_contours_AB_constdchisq.eps");
  c->SaveAs("compare_contours_AB_constdchisq.pdf");
  //again,with externals
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0038);
  osc_param_space->Draw();  
  //h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  //h_dc_kit_68->Draw("cont3same");
  h_dc_kit_90->Draw("cont3same");
  leg->Clear();
  //leg->AddEntry(h_dc_28_fine_68, "T2K A const #Delta#chi^{2} 3#nu 68% CL","l");
  leg->AddEntry(h_dc_28_fine_90, "T2K A const #Delta#chi^{2} 3#nu 90% CL","l");
  //leg->AddEntry(h_dc_kit_68, "T2K B const #Delta#chi^{2} 3#nu 68% CL","l");
  leg->AddEntry(h_dc_kit_90, "T2K B const #Delta#chi^{2} 3#nu 90% CL","l");
  g_run12_90_statsyst->Draw("l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.55);
  leg->Draw();
  bestfit.Draw();
  bestfit_kit.Draw();
  c->SaveAs("compare_contours_AB_constdchisq_ext.png");
  c->SaveAs("compare_contours_AB_constdchisq_ext.jpg");
  c->SaveAs("compare_contours_AB_constdchisq_ext.eps");
  c->SaveAs("compare_contours_AB_constdchisq_ext.pdf");

  //Run 1-3 FC KIT vs VaLOR
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  h_fc_kit_68->Draw("cont3same");
  h_fc_kit_90->Draw("cont3same");
  leg->Clear();
  leg->AddEntry(h_fc_123_68, "A FC 3#nu 68% CL","l");
  leg->AddEntry(h_fc_123_90, "A FC 3#nu 90% CL","l");
  leg->AddEntry(h_fc_kit_68, "B FC 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "B FC 3#nu 90% CL","l");
  leg->SetX2NDC(0.43);
  leg->Draw();
  bestfit.Draw();
  bestfit_kit.Draw();
  c->SaveAs("compare_contours_AB_FC.png");
  c->SaveAs("compare_contours_AB_FC.jpg");
  c->SaveAs("compare_contours_AB_FC.eps");
  c->SaveAs("compare_contours_AB_FC.pdf");
  //again,with externals
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0038);
  osc_param_space->Draw(); 
  //h_fc_123_68->Draw("cont3same");
  h_fc_123_90->Draw("cont3same");
  //h_fc_kit_68->Draw("cont3same");
  h_fc_kit_90->Draw("cont3same");
  leg->Clear();
  //leg->AddEntry(h_fc_123_68, "T2K A FC 3#nu 68% CL","l");
  leg->AddEntry(h_fc_123_90, "T2K A FC 3#nu 90% CL","l");
  //leg->AddEntry(h_fc_kit_68, "T2K B FC 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "T2K B FC 3#nu 90% CL","l");
  g_run12_90_statsyst->Draw("l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  bestfit.Draw();
  bestfit_kit.Draw();
  c->SaveAs("compare_contours_AB_FC_ext.png");
  c->SaveAs("compare_contours_AB_FC_ext.jpg");
  c->SaveAs("compare_contours_AB_FC_ext.eps");
  c->SaveAs("compare_contours_AB_FC_ext.pdf");
  */

  /*
  //compare fine grained surface with not
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();
  h_dc_28_fine_68->Draw("cont3same");
  h_dc_28_fine_90->Draw("cont3same");
  h_dc_123_68->Draw("cont3same");
  h_dc_123_90->Draw("cont3same");
  c->SaveAs("fine_not.pdf");
  */


  //Run 1-3 FC KIT vs VaLOR
  osc_param_space->GetXaxis()->SetRangeUser(0.85,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  osc_param_space->Draw();  
  h_fcall_123_68->SetLineColor(kBlack);;
  h_fcall_123_90->SetLineColor(kBlack);;
  h_fc_kit_68->SetLineColor(kRed);;
  h_fc_kit_90->SetLineColor(kRed);;
  h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  h_fc_kit_68->Draw("cont3same");
  h_fc_kit_90->Draw("cont3same");
  bestfit_2all.Draw();
  bestfit_kit.Draw();
  leg->Clear();
  leg->AddEntry(h_fcall_123_68, "A FC 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "A FC 3#nu 90% CL","l");
  leg->AddEntry(h_fc_kit_68, "B FC 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "B FC 3#nu 90% CL","l");
  leg->SetX2NDC(0.43);
  leg->Draw();
  c->SaveAs("compare_contours_AB_FC.eps");
  leg->Clear();
  leg->AddEntry(h_fcall_123_68, "VaLOR FC 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "VaLOR FC 3#nu 90% CL","l");
  leg->AddEntry(h_fc_kit_68, "KIT FC 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "KIT FC 3#nu 90% CL","l");
  leg->SetX2NDC(0.43);
  leg->Draw();
  c->SaveAs("compare_contours_KITVAL_FC.eps");
  leg->Clear();
  leg->AddEntry(h_fcall_123_68, "Likelihood ratio FC 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "Likelihood ratio FC 3#nu 90% CL","l");
  leg->AddEntry(h_fc_kit_68, "Max likelihood FC 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "Max likelihood FC 3#nu 90% CL","l");
  leg->SetX2NDC(0.43);
  leg->Draw();
  c->SaveAs("compare_contours_ratiomax_FC.eps");
  leg->Clear();
  leg->AddEntry(h_fcall_123_68, "Likelihood ratio FC 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "Likelihood ratio FC 3#nu 90% CL","l");
  leg->AddEntry(h_fc_kit_68, "Max likelihood FC 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "Max likelihood FC 3#nu 90% CL","l");
  leg->SetX2NDC(0.43);
  leg->Draw();
  c->SaveAs("compare_contours_ratiomax_FC.eps");

  //again,with externals
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0038);
  osc_param_space->Draw(); 
  //h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  //h_fc_kit_68->Draw("cont3same");
  h_fc_kit_90->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(2);
  g_run12_90_statsyst->SetLineColor(kBlack);
  ll_minos->SetLineColor(kBlack);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  //g_run12_90_statsyst->Draw("l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  bestfit_2all.Draw();
  bestfit_kit.Draw();
  leg->Clear();
  //leg->AddEntry(h_fcall_123_68, "T2K A 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K A 3#nu 90% CL","l");
  //leg->AddEntry(h_fc_kit_68, "T2K B 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "T2K B 3#nu 90% CL","l");
  //leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("compare_contours_AB_FC_ext.eps");
  leg->Clear();
  //leg->AddEntry(h_fcall_123_68, "T2K VaLOR 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K VaLOR 3#nu 90% CL","l");
  //leg->AddEntry(h_fc_kit_68, "T2K KIT 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "T2K KIT 3#nu 90% CL","l");
  //leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("compare_contours_KITVAL_FC_ext.eps");
  leg->Clear();
  //leg->AddEntry(h_fcall_123_68, "T2K Likelihood ratio 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K Likelihood ratio 3#nu 90% CL","l");
  //leg->AddEntry(h_fc_kit_68, "T2K Max likelihood 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "T2K Max likelihood 3#nu 90% CL","l");
  //leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("compare_contours_ratiomax_FC_ext.eps");

  //again,with externals + old t2k
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0038);
  osc_param_space->Draw(); 
  //h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  //h_fc_kit_68->Draw("cont3same");
  h_fc_kit_90->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(2);
  g_run12_90_statsyst->SetLineColor(kBlack);
  ll_minos->SetLineColor(kBlack);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  g_run12_90_statsyst->Draw("l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  bestfit_2all.Draw();
  bestfit_kit.Draw();
  leg->Clear();
  //leg->AddEntry(h_fcall_123_68, "T2K A 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K A 3#nu 90% CL","l");
  //leg->AddEntry(h_fc_kit_68, "T2K B 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "T2K B 3#nu 90% CL","l");
  leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("compare_contours_AB_FC_extold.eps");
  leg->Clear();
  //leg->AddEntry(h_fcall_123_68, "T2K VaLOR 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K VaLOR 3#nu 90% CL","l");
  //leg->AddEntry(h_fc_kit_68, "T2K KIT 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "T2K KIT 3#nu 90% CL","l");
  leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("compare_contours_KITVAL_FC_extold.eps");
  leg->Clear();
  //leg->AddEntry(h_fcall_123_68, "T2K Likelihood ratio 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K Likelihood ratio 3#nu 90% CL","l");
  //leg->AddEntry(h_fc_kit_68, "T2K Max likelihood 3#nu 68% CL","l");
  leg->AddEntry(h_fc_kit_90, "T2K Max likelihood 3#nu 90% CL","l");
  leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("compare_contours_ratiomax_FC_extold.eps");




  //Run 1-3 FC only with externals
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0038);
  osc_param_space->Draw(); 
  //h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(2);
  g_run12_90_statsyst->SetLineColor(kBlack);
  ll_minos->SetLineColor(kRed);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  //g_run12_90_statsyst->Draw("l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  bestfit_2all.Draw();
  leg->Clear();
  //leg->AddEntry(h_fcall_123_68, "T2K 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K 3#nu 90% CL","l");
  //leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("contour_compare_FC_external.pdf");
  c->SaveAs("contour_compare_FC_external.eps");
  c->SaveAs("contour_compare_FC_external.png");

  //Run 1-3 FC only with externals + old t2k
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0038);
  osc_param_space->Draw(); 
  //h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(2);
  g_run12_90_statsyst->SetLineColor(kBlack);
  ll_minos->SetLineColor(kRed);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  g_run12_90_statsyst->Draw("l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  bestfit_2all.Draw();
  leg->Clear();
  //leg->AddEntry(h_fcall_123_68, "T2K 3#nu 68% CL","l");
  leg->AddEntry(h_fcall_123_90, "T2K 3#nu 90% CL","l");
  leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK L/E 2012 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK zenith 2012 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("contour_compare_FC_external_old.pdf");
  c->SaveAs("contour_compare_FC_external_old.png");
  c->SaveAs("contour_compare_FC_external_old.eps");
  //PAPER
  c->SaveAs("contour_compare_external.pdf");
  c->SaveAs("contour_compare_external.png");
  c->SaveAs("contour_compare_external.eps");

  //Run 1-3 constdchisq only with externals
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0038);
  osc_param_space->Draw(); 
  //h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(2);
  g_run12_90_statsyst->SetLineColor(kBlack);
  ll_minos->SetLineColor(kRed);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  //g_run12_90_statsyst->Draw("l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  bestfit_2all.Draw();
  leg->Clear();
  //leg->AddEntry(h_dcall_123_68, "T2K 3#nu 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "T2K 3#nu 90% CL","l");
  //leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("contour_compare_constdchisq_external.pdf");
  c->SaveAs("contour_compare_constdchisq_external.eps");
  c->SaveAs("contour_compare_constdchisq_external.png");

  //Run 1-3 constdchisq only with externals + old t2k
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0038);
  osc_param_space->Draw(); 
  //h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(2);
  g_run12_90_statsyst->SetLineColor(kBlack);
  ll_minos->SetLineColor(kRed);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  g_run12_90_statsyst->Draw("l");
  ll_minos->Draw("cont3same");
  g_skle->Draw("l");
  g_skzen->Draw("l");
  bestfit_2all.Draw();
  leg->Clear();
  //leg->AddEntry(h_dcall_123_68, "T2K 3#nu 68% CL","l");
  leg->AddEntry(h_dcall_123_90, "T2K 3#nu 90% CL","l");
  leg->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg->AddEntry("g_skle", "SuperK 2012 L/E 2#nu 90% CL", "l");
  leg->AddEntry("g_skzen", "SuperK 2012 zenith 3#nu 90% CL", "l");
  leg->SetX2NDC(0.6);
  leg->Draw();
  c->SaveAs("contour_compare_constdchisq_external_old.pdf");
  c->SaveAs("contour_compare_constdchisq_external_old.eps");
  c->SaveAs("contour_compare_constdchisq_external_old.png");




  //PAPER
  TH2D * osc_param_space2 = new TH2D("osc_param_space2", "",
				     mixing23_nbins, 0.92, 1,
				     dm23sq_nbins, 0.00205, 0.00305);
  osc_param_space2->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space2->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  //osc_param_space2->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space2->GetYaxis()->SetTitleOffset(1.4);
  osc_param_space2->Draw();  
  h_fcall_123_68->SetLineColor(kBlack);;
  h_fcall_123_90->SetLineColor(kBlack);;
  h_fc_kit_68->SetLineColor(kRed);;
  h_fc_kit_90->SetLineColor(kRed);;
  h_fc_kit_68->Draw("cont3same");
  h_fc_kit_90->Draw("cont3same");
  h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  bestfit_2all.Draw();
  bestfit_kit.SetMarkerStyle(24);
  bestfit_kit.SetMarkerColor(kRed);
  bestfit_kit.Draw();
  TLegend * leg2 = new TLegend(0.17,0.7,0.5,0.92);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(h_fcall_123_68, "T2K 68% C.L.","l");
  leg2->AddEntry(h_fcall_123_90, "T2K 90% C.L.","l");
  leg2->AddEntry(h_fc_kit_68, "T2K 68% C.L. alternative analysis","l");
  leg2->AddEntry(h_fc_kit_90, "T2K 90% C.L. alternative analysis","l");
  leg2->SetX2NDC(0.43);
  cout << "text size" << leg2->GetTextSize() << endl;
  leg2->SetTextSize(osc_param_space2->GetXaxis()->GetTitleSize());
  leg2->SetTextSize(osc_param_space2->GetXaxis()->GetTitleSize()*2);
  leg2->SetTextSize(0.05);
  cout << "text size" << leg2->GetTextSize() << endl;
  leg2->Draw();
  c->SaveAs("contour_compare_methods.eps");
  c->SaveAs("contour_compare_methods.pdf");
  c->SaveAs("contour_compare_methods.png");
  leg2->Clear();
  delete leg2;
  leg2 = new TLegend(0.17,0.7,0.65,0.92);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->AddEntry(h_fcall_123_68, "T2K 68% C.L.","l");
  leg2->AddEntry(h_fcall_123_90, "T2K 90% C.L.","l");
  leg2->AddEntry(h_fc_kit_68, "T2K 68% C.L. (alternative analysis)","l");
  leg2->AddEntry(h_fc_kit_90, "T2K 90% C.L. (alternative analysis)","l");
  leg2->SetX2NDC(0.43);
  leg2->SetTextSize(0.05);
  leg2->Draw();
  c->SaveAs("contour_compare_methods_brackets.eps");
  c->SaveAs("contour_compare_methods_brackets.pdf");
  c->SaveAs("contour_compare_methods_brackets.png");

  //PAPER
  //Run 1-3 FC only with externals + old t2k
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0042);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  //h_fcall_123_68->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(3);
  g_run12_90_statsyst->SetLineColor(kBlack);
  h_fcall_octant2only23_90->SetLineStyle(2);
  ll_minos->SetLineColor(kRed);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  g_skzen->Draw("l");
  g_skle->Draw("l");
  ll_minos->Draw("cont3same");
  g_run12_90_statsyst->Draw("l");
  h_fcall_123_90->Draw("cont3same");
  h_fcall_octant2only23_90->Draw("cont3same");
  bestfit_2all.Draw();
  bestfit_2all_octant2.SetMarkerStyle(24);
  bestfit_2all_octant2.Draw();
  leg2->Clear();
  delete leg2;
  leg2 = new TLegend(0.185,0.58,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  cout << leg2->GetTextFont() << endl;
  leg2->SetTextFont(132);
  cout << leg2->GetTextFont() << endl;
  //leg2->AddEntry(h_fcall_123_68, "T2K 3#nu 68% CL","l");
  leg2->AddEntry(h_fcall_123_90, "T2K 3#nu (#theta_{23} #leq #pi/4) 90% CL","l");
  leg2->AddEntry(h_fcall_octant2only23_90, "T2K 3#nu (#theta_{23} #geq #pi/4) 90% CL","l");
  leg2->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg2->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg2->AddEntry("g_skzen", "SK zenith 2012 3#nu 90% CL", "l");
  leg2->AddEntry("g_skle", "SK L/E 2012 2#nu 90% CL", "l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(0.06);
  leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.17, 0.45, 0.30);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bestfit_2all, "T2K 3#nu (#theta_{23} #leq #pi/4) best fit","p");
  legbf->AddEntry(&bestfit_2all_octant2, "T2K 3#nu (#theta_{23} #geq #pi/4) best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(0.06);
  legbf->Draw();
  cout << "text size " << leg2->GetTextSize() << endl;
  c->SaveAs("contour_compare_external.pdf");
  c->SaveAs("contour_compare_external.png");
  c->SaveAs("contour_compare_external.eps");
  cout << osc_param_space->GetXaxis()->GetLabelSize() << endl;
  cout << osc_param_space->GetXaxis()->GetTitleSize() << endl;
  double thick = 2, thick_t2k = 3;
  g_skzen->SetLineWidth(thick);
  g_skle->SetLineWidth(thick);
  ll_minos->SetLineWidth(thick);
  h_fcall_123_90->SetLineWidth(thick_t2k);
  h_fcall_123_68->SetLineWidth(thick_t2k);
  g_run12_90_statsyst->SetLineWidth(thick);
  c->SaveAs("contour_compare_external_thick.pdf");
  c->SaveAs("contour_compare_external_thick.png");
  c->SaveAs("contour_compare_external_thick.eps");


  //PAPER
  //Run 1-3 const dchisq only with externals + old t2k
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0042);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  //h_dcall_123_68->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(3);
  g_run12_90_statsyst->SetLineColor(kBlack);
  h_dcall_octant2only23_90->SetLineStyle(2);
  h_dcall_octant2only23_90->SetLineColor(kBlack);
  ll_minos->SetLineColor(kRed);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  g_skzen->Draw("l");
  g_skle->Draw("l");
  ll_minos->Draw("cont3same");
  g_run12_90_statsyst->Draw("l");
  h_dcall_123_90->Draw("cont3same");
  h_dcall_octant2only23_90->Draw("cont3same");
  bestfit_2all.SetMarkerStyle(23);
  bestfit_2all.Draw();
  bestfit_2all_octant2.SetMarkerStyle(24);
  bestfit_2all_octant2.Draw();
  leg2->Clear();
  delete leg2;
  leg2 = new TLegend(0.185,0.58,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  cout << leg2->GetTextFont() << endl;
  leg2->SetTextFont(132);
  cout << leg2->GetTextFont() << endl;
  //leg2->AddEntry(h_dcall_123_68, "T2K 3#nu 68% CL","l");
  leg2->AddEntry(h_dcall_123_90, "T2K 3#nu (#theta_{23} #leq #pi/4) 90% CL","l");
  leg2->AddEntry(h_dcall_octant2only23_90, "T2K 3#nu (#theta_{23} #geq #pi/4) 90% CL","l");
  leg2->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% CL", "l");
  leg2->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg2->AddEntry("g_skzen", "SK zenith 2012 3#nu 90% CL", "l");
  leg2->AddEntry("g_skle", "SK L/E 2012 2#nu 90% CL", "l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(0.06);
  leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.17, 0.45, 0.30);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bestfit_2all, "T2K 3#nu (#theta_{23} #leq #pi/4) best fit","p");
  legbf->AddEntry(&bestfit_2all_octant2, "T2K 3#nu (#theta_{23} #geq #pi/4) best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(0.06);
  legbf->Draw();
  cout << "text size " << leg2->GetTextSize() << endl;
  c->SaveAs("contour_compare_constdchisq_external.pdf");
  c->SaveAs("contour_compare_constdchisq_external.png");
  c->SaveAs("contour_compare_constdchisq_external.eps");
  cout << osc_param_space->GetXaxis()->GetLabelSize() << endl;
  cout << osc_param_space->GetXaxis()->GetTitleSize() << endl;
  double thick = 2, thick_t2k = 3;
  g_skzen->SetLineWidth(thick);
  g_skle->SetLineWidth(thick);
  ll_minos->SetLineWidth(thick);
  h_dcall_123_90->SetLineWidth(thick_t2k);
  h_dcall_123_68->SetLineWidth(thick_t2k);
  g_run12_90_statsyst->SetLineWidth(thick);
  c->SaveAs("contour_compare_constdchisq_external_thick.pdf");
  c->SaveAs("contour_compare_constdchisq_external_thick.png");
  c->SaveAs("contour_compare_constdchisq_external_thick.eps");




  //PAPER 2nd submission. Using sinsqtheta-marginalised contour PRL ARXIV
  TFile fsinsq("2012results/DATA/sintheta23/combine/combination/sinsq2theta_contour_conversion.root");
  TH2D * hsinsq68 = 0, * hsinsq90 = 0;
  fsinsq.GetObject("hdchisq_68",hsinsq68);
  fsinsq.GetObject("hdchisq_68",hsinsq90);
  hsinsq68->SetContourLevel(0,2.30);
  hsinsq90->SetContourLevel(0,4.61);
  hsinsq68->SetLineWidth(3);
  hsinsq90->SetLineWidth(3);

  //Run 1-3 const dchisq only with externals + old t2k
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0042);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  //h_dcall_123_68->Draw("cont3same");
  g_run12_90_statsyst->SetLineStyle(3);
  g_run12_90_statsyst->SetLineColor(kBlack);
  ll_minos->SetLineColor(kRed);
  g_skle->SetLineColor(kGreen+2);
  g_skzen->SetLineColor(kBlue);
  hsinsq90->SetLineStyle(1);
  hsinsq90->SetLineColor(kBlack);
  g_skzen->Draw("l");
  g_skle->Draw("l");
  ll_minos->Draw("cont3same");
  g_run12_90_statsyst->Draw("l");
  hsinsq90->Draw("cont3same");
  bestfit_2all_from_sinsqtheta.SetMarkerStyle(23);
  bestfit_2all_from_sinsqtheta.Draw();
  leg2->Clear();
  delete leg2;
  leg2 = new TLegend(0.185,0.58,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  cout << leg2->GetTextFont() << endl;
  leg2->SetTextFont(132);
  cout << leg2->GetTextFont() << endl;
  //leg2->AddEntry(h_dcall_123_68, "T2K 3#nu 68% CL","l");
  //leg2->AddEntry(h_dcall_123_90, "T2K 3#nu (#theta_{23} #leq #pi/4) 90% CL","l");
  //leg2->AddEntry(h_dcall_octant2only23_90, "T2K 3#nu (#theta_{23} #geq #pi/4) 90% CL","l");
  leg2->AddEntry(hsinsq90, "T2K 3#nu 90% CL","l");
  leg2->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% C.L.", "l");
  leg2->AddEntry(ll_minos, "MINOS 2013 2#nu 90% C.L.", "l");
  leg2->AddEntry("g_skzen", "SK zenith 2012 3#nu 90% C.L.", "l");
  leg2->AddEntry("g_skle", "SK L/E 2012 2#nu 90% C.L.", "l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(0.06);
  leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.17, 0.45, 0.30);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  //legbf->AddEntry(&bestfit_2all, "T2K 3#nu (#theta_{23} #leq #pi/4) best fit","p");
  //legbf->AddEntry(&bestfit_2all_octant2, "T2K 3#nu (#theta_{23} #geq #pi/4) best fit","p");
  legbf->AddEntry(&bestfit_2all_from_sinsqtheta, "T2K 3#nu best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(0.06);
  legbf->Draw();
  cout << "text size " << leg2->GetTextSize() << endl;
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta.pdf");
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta.png");
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta.eps");

  //add external best-fit points
  TMarker m_minos(0.950,0.00241,24);
  m_minos.SetMarkerColor(kRed);
  m_minos.SetMarkerSize(1.7);
  ll_minos->SetMarkerStyle(24);
  ll_minos->SetMarkerColor(kRed);
  ll_minos->SetMarkerSize(1.7);
  TMarker m_skle(1.000,0.0025,25);
  m_skle.SetMarkerColor(kGreen+2);
  m_skle.SetMarkerSize(1.7);
  g_skle->SetMarkerStyle(25);
  g_skle->SetMarkerColor(kGreen+2);
  g_skle->SetMarkerSize(1.7);
  TMarker m_skzen(0.9775,0.00266,26);
  m_skzen.SetMarkerColor(kBlue);
  m_skzen.SetMarkerSize(1.7);
  g_skzen->SetMarkerStyle(26);
  g_skzen->SetMarkerColor(kBlue);
  g_skzen->SetMarkerSize(1.7);
  TMarker m_t2k12(0.98,0.00265,27);
  m_t2k12.SetMarkerColor(kBlack);
  m_t2k12.SetMarkerSize(1.7);
  g_run12_90_statsyst->SetMarkerStyle(27);
  g_run12_90_statsyst->SetMarkerColor(kBlack);
  g_run12_90_statsyst->SetMarkerSize(1.7);
  hsinsq90->SetMarkerStyle(23);
  hsinsq90->SetMarkerColor(kBlack);
  hsinsq90->SetMarkerSize(1.8);
  //draw
  osc_param_space->Draw();
  m_minos.Draw();
  m_skle.Draw();
  m_skzen.Draw();
  m_t2k12.Draw();
  g_skzen->Draw("l");
  g_skle->Draw("l");
  ll_minos->Draw("cont3same");
  g_run12_90_statsyst->Draw("l");
  hsinsq90->Draw("cont3same");
  bestfit_2all_from_sinsqtheta.SetMarkerStyle(23);
  bestfit_2all_from_sinsqtheta.Draw();
  delete leg2;
  leg2 = new TLegend(0.185,0.58,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(0.06);
  leg2->AddEntry(hsinsq90, "T2K 3#nu 90% C.L.","lp");
  leg2->AddEntry("g_run12_90_statsyst", "T2K 2011 2#nu 90% C.L.", "lp");
  leg2->AddEntry(ll_minos, "MINOS 2013 2#nu 90% C.L.", "lp");
  leg2->AddEntry("g_skzen", "SK zenith 2012 3#nu 90% C.L.", "lp");
  leg2->AddEntry("g_skle", "SK L/E 2012 2#nu 90% C.L.", "lp");
  leg2->Draw();
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta_markers.pdf");
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta_markers.png");
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta_markers.eps");

  //also draw the sinsq2theta contours, to check the difference between the
  //margninalised & 2nd octant contour
  h_dcall_octant2only23_90->SetLineStyle(2);
  h_dcall_octant2only23_90->SetLineColor(kBlack);
  h_dcall_123_90->Draw("cont3same");
  h_dcall_octant2only23_90->Draw("cont3same");
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta_check.pdf");

  //return;

  //THE PRL PLOT, converting y-axis
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  // dm23sq, mixing23 range and step
  double xmixing23_min   = 0.6;
  double xmixing23_max   = 1.5;
  double xmixing23_step  = 0.005;
  double xdm23sq_min     = 1;
  double xdm23sq_max     = 6;
  double xdm23sq_step    = 0.1;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int xmixing23_nbins = 1 + TMath::CeilNint ((xmixing23_max - xmixing23_min) / xmixing23_step);
  int xdm23sq_nbins   = 1 + TMath::CeilNint ((xdm23sq_max   - xdm23sq_min  ) / xdm23sq_step  );
  TH2D * osc_param_spacex = new TH2D("osc_param_space", "",
       xmixing23_nbins, 
       xmixing23_min - 0.5*xmixing23_step, 
       xmixing23_max + 0.5*xmixing23_step,
       xdm23sq_nbins, 
       xdm23sq_min - 0.5*xdm23sq_step, 
       xdm23sq_max + 0.5*xdm23sq_step);
  osc_param_spacex->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_spacex->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (10^{-3} eV^{2}/c^{4})");
  //osc_param_spacex->GetXaxis()->SetTitleOffset(1.2);
  osc_param_spacex->GetYaxis()->SetTitleOffset(0.85);
  osc_param_spacex->GetXaxis()->SetRangeUser(0.82,1.0);
  osc_param_spacex->GetYaxis()->SetRangeUser(1.8,4.2);
  osc_param_spacex->GetXaxis()->SetTitleFont(132);
  osc_param_spacex->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (10^{-3} eV^{2}/c^{4})");
  osc_param_spacex->Draw();
  TMarker xm_minos(0.950,2.41,24);
  xm_minos.SetMarkerColor(kRed);
  xm_minos.SetMarkerSize(1.7);
  TH2D * xll_minos = new TH2D("","",
			      ll_minos->GetNbinsX(), ll_minos->GetXaxis()->GetBinLowEdge(1),        ll_minos->GetXaxis()->GetBinUpEdge(ll_minos->GetNbinsX()),
			      ll_minos->GetNbinsY(), ll_minos->GetYaxis()->GetBinLowEdge(1) * 1000, ll_minos->GetYaxis()->GetBinUpEdge(ll_minos->GetNbinsY()) * 1000);
  for(int ix = 1; ix <=  ll_minos->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= ll_minos->GetNbinsY(); iy++) {
      xll_minos->SetBinContent(ix,iy, ll_minos->GetBinContent(ix,iy));
      xll_minos->SetBinContent(ix,iy, ll_minos->GetBinContent(ix,iy));
    }//iy
  }//ix
  xll_minos->SetContour(1);
  xll_minos->SetContourLevel(0,4.61/2);
  xll_minos->SetLineStyle(1);
  xll_minos->SetMarkerStyle(24);
  xll_minos->SetMarkerColor(kRed);
  xll_minos->SetMarkerSize(1.7);
  xll_minos->SetLineColor(kRed);
  xll_minos->SetLineWidth(ll_minos->GetLineWidth());
  TMarker xm_skle(1.000,2.5,25);
  xm_skle.SetMarkerColor(kGreen+2);
  xm_skle.SetMarkerSize(1.7);
  TGraph * xg_skle = new TGraph(*g_skle);
  for(int i = 0; i < g_skle->GetN(); i++) {
    double xx, yy;
    g_skle ->GetPoint(i, xx, yy);
    xg_skle->SetPoint(i, xx, yy*1000);
  }
  xg_skle->SetMarkerStyle(25);
  xg_skle->SetMarkerColor(kGreen+2);
  xg_skle->SetMarkerSize(1.7);
  TMarker xm_skzen(0.9775,2.66,26);
  xm_skzen.SetMarkerColor(kBlue);
  xm_skzen.SetMarkerSize(1.7);
  TGraph * xg_skzen = new TGraph(*g_skzen);
  for(int i = 0; i < g_skzen->GetN(); i++) {
    double xx, yy;
    g_skzen ->GetPoint(i, xx, yy);
    xg_skzen->SetPoint(i, xx, yy*1000);
  }
  xg_skzen->SetMarkerStyle(26);
  xg_skzen->SetMarkerColor(kBlue);
  xg_skzen->SetMarkerSize(1.7);
  TMarker xm_t2k12(0.98,2.65,27);
  xm_t2k12.SetMarkerColor(kBlack);
  xm_t2k12.SetMarkerSize(1.7);
  TGraph * xg_run12_90_statsyst = new TGraph(*g_run12_90_statsyst);
  for(int i = 0; i < g_run12_90_statsyst->GetN(); i++) {
    double xx, yy;
    g_run12_90_statsyst ->GetPoint(i, xx, yy);
    xg_run12_90_statsyst->SetPoint(i, xx, yy*1000);
  }
  xg_run12_90_statsyst->SetMarkerStyle(27);
  xg_run12_90_statsyst->SetMarkerColor(kBlack);
  xg_run12_90_statsyst->SetMarkerSize(1.7);
#line 0
  Double_t xaxis[1000];//[hsinsq90->GetNbinsY()+1];
  for(int ix = 0; ix <= hsinsq90->GetNbinsX(); ix++) {
    cout << hsinsq90->GetXaxis()->GetBinUpEdge(ix) << endl;
    xaxis[ix] = hsinsq90->GetXaxis()->GetBinUpEdge(ix);
    cout << xaxis[ix] << endl;
  }
  //yaxis[hsinsq90->GetNbinsY()+1] = 1000;
  TH2D * xhsinsq90 = new TH2D("","",
			      hsinsq90->GetNbinsX(), xaxis,//hsinsq90->GetXaxis()->GetBinLowEdge(1),        hsinsq90->GetXaxis()->GetBinUpEdge(hsinsq90->GetNbinsX()),
			      hsinsq90->GetNbinsY(), hsinsq90->GetYaxis()->GetBinLowEdge(1) * 1000, hsinsq90->GetYaxis()->GetBinUpEdge(hsinsq90->GetNbinsY()) * 1000);
  for(int ix = 1; ix <=  hsinsq90->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= hsinsq90->GetNbinsY(); iy++) {
      xhsinsq90->SetBinContent(ix,iy, hsinsq90->GetBinContent(ix,iy));
      xhsinsq90->SetBinContent(ix,iy, hsinsq90->GetBinContent(ix,iy));
    }//iy
  }//ix
  xhsinsq90->SetMarkerStyle(23);
  xhsinsq90->SetMarkerColor(kBlack);
  xhsinsq90->SetMarkerSize(1.8);
  xhsinsq90->SetContour(1);
  xhsinsq90->SetContourLevel(0,4.61);
  xhsinsq90->SetLineWidth(3);
  //draw
  osc_param_spacex->Draw();
  xm_minos.Draw();
  xm_skle.Draw();
  xm_skzen.Draw();
  xm_t2k12.Draw();
  xg_skzen->Draw("l");
  xg_skle->Draw("l");
  xll_minos->Draw("cont3same");
  xg_run12_90_statsyst->Draw("l");
  xhsinsq90->Draw("cont3same");
  //xhsinsq90->Draw("cont3colz");
  xbestfit_2all_from_sinsqtheta.SetMarkerStyle(23);
  xbestfit_2all_from_sinsqtheta.Draw();
  delete leg2;
  leg2 = new TLegend(0.13,0.59,0.4,0.96);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(0.06);
  leg2->AddEntry(xhsinsq90, "T2K 3#nu 90% C.L.","lp");
  xg_run12_90_statsyst->SetName("xg_run12_90_statsyst");
  leg2->AddEntry("xg_run12_90_statsyst", "T2K 2011 2#nu 90% C.L.", "lp");
  leg2->AddEntry(xll_minos, "MINOS 2013 2#nu 90% C.L.", "lp");
  xg_skzen->SetName("xg_skzen");
  leg2->AddEntry("xg_skzen", "SK zenith 2012 3#nu 90% C.L.", "lp");
  xg_skle->SetName("xg_skle");
  leg2->AddEntry("xg_skle", "SK L/E 2012 2#nu 90% C.L.", "lp");
  leg2->Draw();
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta_markers.pdf");
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta_markers.png");
  c->SaveAs("contour_compare_constdchisq_external_sinsqtheta_markers.eps");



  //LEPTON PHOTON
  //Run 1-3 FC only with externals (including icecube)
  osc_param_space->GetXaxis()->SetRangeUser(0.72,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.00134,0.00534);
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  h_fcall_123_68->SetLineWidth(2);
  h_fcall_123_90->SetLineWidth(2);
  ll_minos->SetLineColor(kRed);
  ll_minos->SetLineWidth(2);
  g_skle->SetLineColor(kBlue);
  g_skle->SetLineWidth(2);
  g_skle->SetLineStyle(2);
  g_skzen->SetLineColor(kBlue);
  g_skzen->SetLineWidth(2);
  g_ice->SetLineColor(kGreen+2);
  g_ice->SetLineWidth(2);
  //draw
  //TGraphPainter paint;
  //paint.Smooth(g_ice2, g_ice->GetN(), g_ice->GetX(), g_ice->GetY(), 0); //DOES NOTHING
  g_ice->Draw("l");
  g_skzen->Draw("l");
  g_skle->Draw("l");
  ll_minos->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  //h_fcall_123_68->Draw("cont3same");
  //bestfit_2all.Draw();
  leg2->Clear();
  delete leg2;
  leg2 = new TLegend(0.185,0.63,0.45,0.92);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(132);
  //leg2->AddEntry(h_fcall_123_68, "T2K 2013 3#nu 68% CL","l");
  leg2->AddEntry(h_fcall_123_90, "T2K 2013 3#nu 90% CL","l");
  leg2->AddEntry(ll_minos, "MINOS 2013 2#nu 90% CL", "l");
  leg2->AddEntry("g_skzen", "SK zenith 2012 3#nu 90% CL", "l");
  leg2->AddEntry("g_skle", "SK L/E 2012 2#nu 90% CL", "l");
  leg2->AddEntry("g_ice", "IceCube-79 2013 2#nu 90% CL", "l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(0.06);
  leg2->Draw();
  legbf->Clear();
  delete legbf;
  legbf = new TLegend(0.185, 0.11, 0.45, 0.35);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bestfit_2all, "T2K 3#nu best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(0.06);
  //legbf->Draw();
  c->SaveAs("contour_compare_external_icecube.pdf");
  c->SaveAs("contour_compare_external_icecube.png");
  c->SaveAs("contour_compare_external_icecube.eps");
  //make lines thicker
  double thicker = 5;
  g_ice->SetLineWidth(thicker);
  g_skzen->SetLineWidth(thicker);
  g_skle->SetLineWidth(thicker);
  ll_minos->SetLineWidth(thicker);
  h_fcall_123_90->SetLineWidth(thicker);
  h_fcall_123_68->SetLineWidth(thicker);
  c->SaveAs("contour_compare_external_icecube_thick.pdf");
  c->SaveAs("contour_compare_external_icecube_thick.png");
  c->SaveAs("contour_compare_external_icecube_thick.eps");

  //compare octant 1, octant 2, 2 flavour
  TFile f2nu("2012results/DATA/2nu/combination/fitresults_3.root");
  TH2D * h_dcall_2nu = 0;
  f2nu.GetObject("hdchisq", h_dcall_2nu);
  TH2D * h_dcall_2nu_68 = new TH2D(*h_dcall_2nu);
  TH2D * h_dcall_2nu_90 = new TH2D(*h_dcall_2nu);
  h_dcall_2nu_68->SetContour(1);
  h_dcall_2nu_68->SetContourLevel(0,2.30);
  h_dcall_2nu_68->SetLineStyle(2);
  h_dcall_2nu_68->SetLineWidth(linewidth_val);
  h_dcall_2nu_90->SetContour(1);
  h_dcall_2nu_90->SetContourLevel(0,4.61);
  h_dcall_2nu_90->SetLineStyle(1);
  h_dcall_2nu_90->SetLineWidth(linewidth_val);
  //style1
  h_dcall_123_90->SetLineWidth(linewidth_val);
  h_dcall_123_90->SetLineStyle(1);
  h_dcall_123_90->SetLineColor(kBlack);
  h_dcall_123_68->SetLineWidth(linewidth_val);
  h_dcall_123_68->SetLineStyle(2);
  h_dcall_123_68->SetLineColor(kBlack);
  bestfit_2all.SetMarkerColor(kBlack);
  bestfit_2all.SetMarkerStyle(25);
  bestfit_2all.SetMarkerSize(1.8);
  //style2
  h_dcall_octant2only23_90->SetLineWidth(linewidth_val);
  h_dcall_octant2only23_90->SetLineStyle(1);
  h_dcall_octant2only23_90->SetLineColor(kRed);
  h_dcall_octant2only23_68->SetLineWidth(linewidth_val);
  h_dcall_octant2only23_68->SetLineStyle(2);
  h_dcall_octant2only23_68->SetLineColor(kRed);
  bestfit_2all_2nu.SetMarkerColor(kRed);
  bestfit_2all_2nu.SetMarkerStyle(24);
  bestfit_2all_2nu.SetMarkerSize(1.8);
  //style3
  h_dcall_2nu_90->SetLineWidth(linewidth_val);
  h_dcall_2nu_90->SetLineStyle(1);
  h_dcall_2nu_90->SetLineColor(kBlue);
  h_dcall_2nu_68->SetLineWidth(linewidth_val);
  h_dcall_2nu_68->SetLineStyle(2);
  h_dcall_2nu_68->SetLineColor(kBlue);
  bestfit_2all_octant2.SetMarkerColor(kBlue);
  bestfit_2all_octant2.SetMarkerStyle(26);
  bestfit_2all_octant2.SetMarkerSize(1.8);
  //draw
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.14);
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1.0);
  osc_param_space->GetYaxis()->SetRangeUser(0.0018,0.0032);
  cout <<   osc_param_space->GetXaxis()->GetTitleFont() << endl;
  osc_param_space->GetXaxis()->SetTitleFont(132);
  osc_param_space->Draw(); 
  h_dcall_123_68->Draw("cont3same");
  h_dcall_123_90->Draw("cont3same");
  h_dcall_octant2only23_68->Draw("cont3same");
  h_dcall_octant2only23_90->Draw("cont3same");
  h_dcall_2nu_68->Draw("cont3same");
  h_dcall_2nu_90->Draw("cont3same");
  bestfit_2all.Draw("P");
  bestfit_2all_2nu.Draw("P");
  bestfit_2all_octant2.Draw("P");
  leg2->Clear();
  delete leg2;
  leg2 = new TLegend(0.185,0.78,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  cout << leg2->GetTextFont() << endl;
  leg2->SetTextFont(132);
  cout << leg2->GetTextFont() << endl;
  //leg2->AddEntry(h_dcall_123_68, "T2K 3#nu 68% CL","l");
  leg2->AddEntry(h_dcall_123_90,           "T2K #theta_{23} #leq #pi/4 90% CL","l");
  leg2->AddEntry(h_dcall_octant2only23_90, "T2K #theta_{23} #geq #pi/4 90% CL","l");
  leg2->AddEntry(h_dcall_2nu_90,           "T2K 2-flavour 90% CL",             "l");
  leg2->SetX2NDC(0.6);
  leg2->SetTextSize(0.05);
  leg2->Draw();
  TLegend * legbf = new TLegend(0.185, 0.17, 0.45, 0.30);
  legbf->SetFillColor(0);
  legbf->SetFillStyle(0);
  legbf->SetBorderSize(0);
  legbf->SetTextFont(132);
  legbf->AddEntry(&bestfit_2all,         "T2K #theta_{23} #leq #pi/4 best fit","p");
  legbf->AddEntry(&bestfit_2all_octant2, "T2K #theta_{23} #geq #pi/4 best fit","p");
  legbf->AddEntry(&bestfit_2all_2nu,     "T2K 2-flavour best fit","p");
  legbf->SetX2NDC(0.6);
  legbf->SetTextSize(0.05);
  legbf->Draw();
  c->SaveAs("contour_compare_octants_2nu.pdf");
  c->SaveAs("contour_compare_octants_2nu.png");
  c->SaveAs("contour_compare_octants_2nu.eps");




  ///checking contour changes
  h_fcall_123_68->SetLineWidth(1);
  h_fcall_123_90->SetLineWidth(1);
  TH2D * h_fcall_123_copy_68 = new TH2D(*h_fcall_123_68);
  TH2D * h_fcall_123_copy_90 = new TH2D(*h_fcall_123_90);
  h_fcall_123_copy_68->SetContourLevel(0, h_fcall_123_copy_68->GetContourLevel(0) + 0.01);
  h_fcall_123_copy_90->SetContourLevel(0, h_fcall_123_copy_90->GetContourLevel(0) + 0.01);
  h_fcall_123_copy_68->SetLineColor(kRed);
  h_fcall_123_copy_90->SetLineColor(kRed);
  osc_param_space->GetXaxis()->SetRangeUser(0.9,1);
  osc_param_space->GetYaxis()->SetRangeUser(0.002,0.003);
  osc_param_space->Draw();
  h_fcall_123_copy_68->Draw("cont3same");
  h_fcall_123_copy_90->Draw("cont3same");
  h_fcall_123_68->Draw("cont3same");
  h_fcall_123_90->Draw("cont3same");
  c->SaveAs("contour_check.pdf");





}
