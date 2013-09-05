//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

/*
#include <iomanip.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include "Conventions.h"
#include "Utilities.h"
#include "SuperkExpectation.h"
#include "Opt.h"
#include "BANFF.h"
#include "FitParamList.h"

#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH1.h"
#include "TSystem.h"
#include "TString.h"

#ifndef __CINT__
void test_sk_ereco_spectrum_nominal_withsys
(
  const char * sk_sample            ,
  //const char * hypothesis           ,
  double       input_mixing12       ,
  double       input_mixing13       ,
  double       input_mixing23       ,
  double       input_deltacp        ,
  double       input_dm12sq         ,
  double       input_dm23sq         ,
  bool         disable_banff_tuning ,
  const char * config_file          ,
  bool         print_spectrum       ,
  //sk
  double       input_sk_escale ,
  double       input_sk_ccqeoth ,
  double       input_sk_ccqerc0 ,
  double       input_sk_ccqerc1 ,
  double       input_sk_ccqerc2 ,
  double       input_sk_ccnqe ,
  double       input_sk_nc ,
  double       input_sk_nue ,
  //banff flux
  double       input_banff_00 ,
  double       input_banff_01 ,
  double       input_banff_02 ,
  double       input_banff_03 ,
  double       input_banff_04 ,
  double       input_banff_05 ,
  double       input_banff_06 ,
  double       input_banff_07 ,
  double       input_banff_08 ,
  double       input_banff_09 ,
  double       input_banff_10 ,
  double       input_banff_11 ,
  double       input_banff_12 ,
  double       input_banff_13 ,
  double       input_banff_14 ,
  double       input_banff_15 ,
  //banff xsec
  double       input_banff_16 ,
  double       input_banff_17 ,
  double       input_banff_18 ,
  double       input_banff_19 ,
  double       input_banff_20 ,
  double       input_banff_21 ,
  double       input_banff_22 ,
  //xsec
  double       input_xsec_pf ,
  double       input_xsec_sf ,
  double       input_xsec_wshape ,
  double       input_xsec_ccoth ,
  double       input_xsec_piless ,
  double       input_xsec_binde ,
  double       input_xsec_cccoh ,
  double       input_xsec_nc1pic ,
  double       input_xsec_ncoth ,
  double       input_xsec_nuetonumu ,
  double       input_xsec_antinutonu ,
  //fsi
  double       input_fsi_0 ,
  double       input_fsi_1 ,
  double       input_fsi_2 ,
  double       input_fsi_3 ,
  double       input_fsi_4 ,
  double       input_fsi_5
 );
#endif
*/

double * create_syst_array(int mode = 0)
{
  double sys[48] = {0};

  //skdet
  sys[0] = 0;
  sys[1] = 0;
  sys[2] = 0;
  sys[3] = 0;
  sys[4] = 0;
  sys[5] = 0;
  sys[6] = 0;
  sys[7] = 0;
  //banff
  sys[8] =     0.994975356;
  sys[9] =     1.006686263;
  sys[10] =     1.001717758;
  sys[11] =     0.962504019;
  sys[12] =     0.900565463;
  sys[13] =     0.909521076;
  sys[14] =     0.936021362;
  sys[15] =     0.967047925;
  sys[16] =     0.988344347;
  sys[17] =     0.997112586;
  sys[18] =     1.003840461;
  sys[19] =     0.980269599;
  sys[20] =     0.957568592;
  sys[21] =     0.946720515;
  sys[22] =     0.953040190;
  sys[23] =     0.984581514;
  sys[24] =     1.048517110;
  sys[25] =     1.010808890;
  sys[26] =     0.950922419;
  sys[27] =     0.708408251;
  sys[28] =     1.346600558;
  sys[29] =     1.369718196;
  sys[30] =     1.015765518;
  //xsec
  sys[31] = 0;
  sys[32] = 0;
  sys[33] = 0;
  sys[34] = 0;
  sys[35] = 0;
  sys[36] = 0;
  sys[37] = 0;
  sys[38] = 0;
  sys[39] = 0;
  sys[40] = 0;
  sys[41] = 0;
  //fsi
  sys[42] = 1;
  sys[43] = 1;
  sys[44] = 1;
  sys[45] = 1;
  sys[46] = 1;
  sys[47] = 1;

  if(mode == 1) { //2+8 run1-3 best fit

    sys[0]  = -1.93489e-05;
    sys[5]  =  -0.044303;
    sys[6]  = 0.254787;
    sys[24] =  1.07436;
    sys[26] = 0.94079;
    sys[32] = 1.30807e-07;
    sys[36] = 0.064003;
    sys[35] =  -0.165958;

    double sigma_maqe  = (1.07436 - 1.048517110) / 0.02564;
    double sigma_ccqe0 = (0.94079 - 0.950922419) / 0.007343;

    sys[8] += (sigma_maqe * -0.004611) + (sigma_ccqe0 * -0.003242);
    sys[9] += (sigma_maqe * -0.005068) + (sigma_ccqe0 * -0.003668);
    sys[10] += (sigma_maqe * -0.004867) + (sigma_ccqe0 * -0.003606);
    sys[11] += (sigma_maqe * -0.004987) + (sigma_ccqe0 * -0.003691);
    sys[12] += (sigma_maqe * -0.005563) + (sigma_ccqe0 * -0.004084);
    sys[13] += (sigma_maqe * -0.004923) + (sigma_ccqe0 * -0.003193);
    sys[14] += (sigma_maqe * -0.004009) + (sigma_ccqe0 * -0.002345);
    sys[15] += (sigma_maqe * -0.003085) + (sigma_ccqe0 * -0.001678);
    sys[16] += (sigma_maqe * -0.002521) + (sigma_ccqe0 * -0.001363);
    sys[17] += (sigma_maqe * -0.002203) + (sigma_ccqe0 * -0.000759);
    sys[18] += (sigma_maqe * -0.001638) + (sigma_ccqe0 * -0.0001584);
    sys[19] += (sigma_maqe * -0.004763) + (sigma_ccqe0 * -0.003229);
    sys[20] += (sigma_maqe * -0.004421) + (sigma_ccqe0 * -0.00293);
    sys[21] += (sigma_maqe * -0.004349) + (sigma_ccqe0 * -0.002869);
    sys[22] += (sigma_maqe * -0.003582) + (sigma_ccqe0 * -0.002077);
    sys[23] += (sigma_maqe * -0.002025) + (sigma_ccqe0 * -0.0009498);

    sys[25] += (sigma_maqe * -0.001327) + (sigma_ccqe0 * 0.0009846);

    sys[27] += (sigma_maqe * -0.0158) + (sigma_ccqe0 * 0.006794);
    sys[28] += (sigma_maqe * 0.0006522) + (sigma_ccqe0 * 0.00107);
    sys[29] += (sigma_maqe * 0.003963) + (sigma_ccqe0 * 0.003639);
    sys[30] += (sigma_maqe * -0.002518) + (sigma_ccqe0 * 0.001492);
  }
  else if(mode == 2) { //2+all run1-3 best fit
    sys[0] =    -3.3381e-06     ;
    sys[2] =      0.0125947       ;
    sys[3] =      0.0043772       ;
    sys[4] =      0.121606        ;
    sys[1] =      -0.000131565    ;
    sys[5] =        -0.049803       ;
    sys[6] =   0.246003        ;
    sys[7] =        -0.00195244     ;
    sys[8] =      0.991546        ;
    sys[9] =      1.00095 ;
    sys[10] =      0.995547        ;
    sys[11] =      0.957774        ;
    sys[12] =      0.896189        ;
    sys[13] =      0.910326        ;
    sys[14] =      0.938223        ;
    sys[15] =      0.968709        ;
    sys[16] =      0.987996        ;
    sys[17] =      0.995195        ;
    sys[18] =     1.00104 ;
    sys[19] =     0.979558        ;
    sys[20] =     0.957982        ;
    sys[21] =     0.947926        ;
    sys[22] =     0.955453        ;
    sys[23] =     0.984811        ;
    sys[24] =     1.07357 ;
    sys[25] =     1.01054 ;
    sys[26] =     0.945412        ;
    sys[27] =     0.642925        ;
    sys[28] =     1.38322 ;
    sys[29] =     1.37186 ;
    sys[30] =     1.01332 ;
    sys[31] =     -0.0145648      ;
    sys[32] =     5.06505e-09     ;
    sys[36] =  0.0624528       ;
    sys[33] =         -0.0496601      ;
    sys[34] =     -0.011389       ;
    sys[35] =         -0.166358       ;
    sys[37] =     0.00470345      ;
    sys[38] =         0.0328214       ;
    sys[39] =     0.0507135       ;
    sys[40] =    0.00363545      ;
    sys[41] =       0.0337973       ;
    sys[42] =        1       ;
    sys[43] =        0.99999 ;
    sys[44] =        0.999951        ;
    sys[45] =        0.998348        ;
    sys[46] =        1.00002 ;
    sys[47] =        1.00171 ;
  }//mode == 2

  else if(mode == 3) {//fake data set 1
    sys[6] = 1;
  }
  else if(mode == 4) {//fake data set 2
    sys[31] = 1;
    sys[24] = 1;
  }
  else if(mode == 5) {//fake data set 3
    sys[0] = -1;
  }
  else if(mode == 6) {//fake data set 4
    sys[8] = sys[8]*1.05;
    sys[9] = sys[9]*1.05;
    sys[10] = sys[10]*1.05;
    sys[11] = sys[11]*1.05;
  }

  else if(mode == 7) { //2+all run1-3 best fit in second octant
    sys[0] =    -1.91287e-06     ; //sk escale
    sys[2] =      0.0132393       ; //sk ccqe rc ...
    sys[3] =      0.00527075       ;
    sys[4] =      0.12097        ;
    sys[1] =      0.000482819    ; //sk ccqe other
    sys[5] =        -0.0480843       ; //sk ccnqe
    sys[6] =   0.250507        ; //sk nc
    sys[7] =        -0.00188268     ; //sk nue cc
    sys[8] =      0.991577        ; //banff ...
    sys[9] =      1.001 ;
    sys[10] =      0.99561        ;
    sys[11] =      0.957817        ;
    sys[12] =      0.896185        ;
    sys[13] =      0.910331        ;
    sys[14] =      0.938242        ;
    sys[15] =      0.968715        ;
    sys[16] =      0.988003        ;
    sys[17] =      0.99524        ;
    sys[18] =     1.00111 ;
    sys[19] =     0.979593        ;
    sys[20] =     0.958027        ;
    sys[21] =     0.947962        ;
    sys[22] =     0.955518        ;
    sys[23] =     0.98489        ;
    sys[24] =     1.07366 ;
    sys[25] =     1.01059 ;
    sys[26] =     0.945387        ;
    sys[27] =     0.642753        ;
    sys[28] =     1.38331 ;
    sys[29] =     1.37194 ;
    sys[30] =     1.01316 ;
    sys[31] =     -0.0145107      ; //pf
    sys[32] =     1.76098e-09     ; //sf
    sys[36] =  0.0631316       ; //bindE
    sys[33] =         -0.0483071      ; //wshape
    sys[34] =     -0.0113211       ; //ccoth shape
    sys[35] =         -0.16334       ; //pdd
    sys[37] =     0.00472817      ; //cccoh
    sys[38] =         0.0338255       ; //nc1pic
    sys[39] =     0.0514926       ;//ncoth
    sys[40] =    0.000669648      ;//nuetonumu
    sys[41] =       0.0307822       ;//nutonubar
    sys[42] =        1       ; //fsi ...
    sys[43] =        0.99999 ;
    sys[44] =        0.999953        ;
    sys[45] =        0.9984        ;
    sys[46] =        1.00003 ;
    sys[47] =        1.00174 ;
  }//mode == 7

  else if(mode == 8) { //2+all run1-3 best fit for sinsqtheta
    sys[0] =    -1.23446e-05     ; //sk escale
    sys[2] =      0.0133159       ; //sk ccqe rc ...
    sys[3] =      0.00542347       ;
    sys[4] =      0.121166        ;
    sys[1] =      0.000613943    ; //sk ccqe other
    sys[5] =        -0.0479978       ; //sk ccnqe
    sys[6] =   0.250777        ; //sk nc
    sys[7] =        -0.00231493     ; //sk nue cc
    sys[8] =      0.991613        ; //banff ...
    sys[9] =      1.00103 ;
    sys[10] =      0.995631        ;
    sys[11] =      0.957844        ;
    sys[12] =      0.896238        ;
    sys[13] =      0.910379        ;
    sys[14] =      0.938265        ;
    sys[15] =      0.968748        ;
    sys[16] =      0.988027        ;
    sys[17] =      0.995229        ;
    sys[18] =     1.00106	 ;
    sys[19] =     0.97962       ;
    sys[20] =     0.958024        ;
    sys[21] =     0.947954        ;
    sys[22] =     0.955464        ;
    sys[23] =     0.984801        ;
    sys[24] =     1.07351 ;
    sys[25] =     1.01057 ;
    sys[26] =     0.945421        ;
    sys[27] =     0.642812        ;
    sys[28] =     1.38331 ;
    sys[29] =     1.37192 ;
    sys[30] =     1.01324 ;
    sys[31] =     -0.0143919      ; //pf
    sys[32] =     6.11723e-07     ; //sf
    sys[36] =  0.0632484       ; //bindE
    sys[33] =         -0.0482414      ; //wshape
    sys[34] =     -0.0112547       ; //ccoth shape
    sys[35] =         -0.163483       ; //pdd
    sys[37] =     0.00481022      ; //cccoh
    sys[38] =         0.0338681       ; //nc1pic
    sys[39] =     0.0515436       ;//ncoth
    sys[40] =    -6.9191e-05      ;//nuetonumu
    sys[41] =       0.0303501       ;//nutonubar
    sys[42] =        1       ; //fsi ...
    sys[43] =        0.99999 ;
    sys[44] =        0.999953        ;
    sys[45] =        0.998393        ;
    sys[46] =        1.00003 ;
    sys[47] =        1.00179 ;
  }//mode == 8

  return sys;
}

void test_sk_ereco_spectrum_nominal_withsys
(
 int systmode = 0,
 //  double *     inputs_syst          = 0, /* see order below */
  //const char * hypothesis           = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12       = -99,                            /* */
  double       input_mixing13       = -99,                            /* */
  double       input_mixing23       = -99,                            /* */
  double       input_deltacp        = -99,                            /* */
  double       input_dm12sq         = -99,                            /* */
  double       input_dm23sq         = -99,                            /* */
  const char * sk_sample            = "1Rmu",                         /* options: "1Rmu", "1Re" */
  //bool         disable_banff_tuning = false,                          /* */
  const char * config_file          = "analysis_configuration.xml",   /* */
  bool         print_spectrum       = false,                          /* */
  const char * mixparamchoice       = "sinsq2theta",
  const char * hierarchy            = "normal"
  /*
  //sk
  double       input_sk_escale = 0,
  double       input_sk_ccqeoth = 0,
  double       input_sk_ccqerc0 = 0,
  double       input_sk_ccqerc1 = 0,
  double       input_sk_ccqerc2 = 0,
  double       input_sk_ccnqe = 0,
  double       input_sk_nc = 0,
  double       input_sk_nue = 0,
  //banff flux
  double       input_banff_00 = 0,
  double       input_banff_01 = 0,
  double       input_banff_02 = 0,
  double       input_banff_03 = 0,
  double       input_banff_04 = 0,
  double       input_banff_05 = 0,
  double       input_banff_06 = 0,
  double       input_banff_07 = 0,
  double       input_banff_08 = 0,
  double       input_banff_09 = 0,
  double       input_banff_10 = 0,
  double       input_banff_11 = 0,
  double       input_banff_12 = 0,
  double       input_banff_13 = 0,
  double       input_banff_14 = 0,
  double       input_banff_15 = 0,
  //banff xsec
  double       input_banff_16 = 0,
  double       input_banff_17 = 0,
  double       input_banff_18 = 0,
  double       input_banff_19 = 0,
  double       input_banff_20 = 0,
  double       input_banff_21 = 0,
  double       input_banff_22 = 0,
  //xsec
  double       input_xsec_pf = 0,
  double       input_xsec_sf = 0,
  double       input_xsec_wshape = 0,
  double       input_xsec_ccoth = 0,
  double       input_xsec_piless = 0,
  double       input_xsec_binde = 0,
  double       input_xsec_cccoh = 0,
  double       input_xsec_nc1pic = 0,
  double       input_xsec_ncoth = 0,
  double       input_xsec_nuetonumu = 0,
  double       input_xsec_antinutonu = 0,
  //fsi
  double       input_fsi_0 = 0,
  double       input_fsi_1 = 0,
  double       input_fsi_2 = 0,
  double       input_fsi_3 = 0,
  double       input_fsi_4 = 0,
  double       input_fsi_5 = 0
  */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
  FitParamList * fparm = FitParamList::NamedParameterList("nuosc3f_numudisapp_fix4+nuisance");
  if(fparm->PhysicsHypothesis() == conv::kHypNuMuDisapp2F) {
    fparm->Update(conv::kPhypMixing23,  input_mixing23);
    fparm->Update(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(fparm->PhysicsHypothesis() == conv::kHypOsc3F) {
    fparm->Update(conv::kPhypMixing12,  input_mixing12);
    fparm->Update(conv::kPhypMixing13,  input_mixing13);
    fparm->Update(conv::kPhypMixing23,  input_mixing23);
    fparm->Update(conv::kPhypDeltaCP,   input_deltacp);
    fparm->Update(conv::kPhypDm12sq,    input_dm12sq);
    fparm->Update(conv::kPhypDm23sq,    input_dm23sq);
  }
  //set systematics to nominal
  fparm->Update(conv::kNuisXSecSpecFunc, 0);//spec func starts at 0.5!

  double * inputs_syst = create_syst_array(systmode);

  //sk
  fparm->Update(conv::kNuisSuperkEnergyScale,inputs_syst[0]);
  fparm->Update(conv::kNuisSuperkEff1RmuCCQEoth,inputs_syst[1]);
  fparm->Update(conv::kNuisSuperkEff1RmuCCQErc0u,inputs_syst[2]);
  fparm->Update(conv::kNuisSuperkEff1RmuCCQErc1u,inputs_syst[3]);
  fparm->Update(conv::kNuisSuperkEff1RmuCCQErc2u,inputs_syst[4]);
  fparm->Update(conv::kNuisSuperkEff1RmuCCnonQE,inputs_syst[5]);
  fparm->Update(conv::kNuisSuperkEff1RmuNC,inputs_syst[6]);
  fparm->Update(conv::kNuisSuperkEff1RmuNueCC,inputs_syst[7]);
  //banff flux
  fparm->Update(conv::kNuisBANFF00,inputs_syst[8]);
  fparm->Update(conv::kNuisBANFF01,inputs_syst[9]);
  fparm->Update(conv::kNuisBANFF02,inputs_syst[10]);
  fparm->Update(conv::kNuisBANFF03,inputs_syst[11]);
  fparm->Update(conv::kNuisBANFF04,inputs_syst[12]);
  fparm->Update(conv::kNuisBANFF05,inputs_syst[13]);
  fparm->Update(conv::kNuisBANFF06,inputs_syst[14]);
  fparm->Update(conv::kNuisBANFF07,inputs_syst[15]);
  fparm->Update(conv::kNuisBANFF08,inputs_syst[16]);
  fparm->Update(conv::kNuisBANFF09,inputs_syst[17]);
  fparm->Update(conv::kNuisBANFF10,inputs_syst[18]);
  fparm->Update(conv::kNuisBANFF11,inputs_syst[19]);
  fparm->Update(conv::kNuisBANFF12,inputs_syst[20]);
  fparm->Update(conv::kNuisBANFF13,inputs_syst[21]);
  fparm->Update(conv::kNuisBANFF14,inputs_syst[22]);
  fparm->Update(conv::kNuisBANFF15,inputs_syst[23]);
  //banff xsec
  fparm->Update(conv::kNuisBANFFMAQE,inputs_syst[24]);
  fparm->Update(conv::kNuisBANFFMARES,inputs_syst[25]);
  fparm->Update(conv::kNuisBANFFCCQEE0,inputs_syst[26]);
  fparm->Update(conv::kNuisBANFFCCQEE1,inputs_syst[27]);
  fparm->Update(conv::kNuisBANFFCCQEE2,inputs_syst[28]);
  fparm->Update(conv::kNuisBANFFCC1piE0,inputs_syst[29]);
  fparm->Update(conv::kNuisBANFFCC1piE1,inputs_syst[30]);
  //xsec
  fparm->Update(conv::kNuisXSecFermiMom,inputs_syst[31]);
  fparm->Update(conv::kNuisXSecSpecFunc,inputs_syst[32]);
  fparm->Update(conv::kNuisXSecWshape,inputs_syst[33]);
  fparm->Update(conv::kNuisXSecShapeCCoth,inputs_syst[34]);
  fparm->Update(conv::kNuisXSecPiLessDelta,inputs_syst[35]);
  fparm->Update(conv::kNuisXSecBindE,inputs_syst[36]);
  fparm->Update(conv::kNuisXSecNormCCcoh,inputs_syst[37]);
  fparm->Update(conv::kNuisXSecNormNC1piC,inputs_syst[38]);
  fparm->Update(conv::kNuisXSecNormNCoth,inputs_syst[39]);
  fparm->Update(conv::kNuisXSecNormNuEToNuMu,inputs_syst[40]);
  fparm->Update(conv::kNuisXSecNormAntiNuToNu,inputs_syst[41]);
  //fsi
  fparm->Update(conv::kNuisFSI000,inputs_syst[42]);
  fparm->Update(conv::kNuisFSI001,inputs_syst[43]);
  fparm->Update(conv::kNuisFSI002,inputs_syst[44]);
  fparm->Update(conv::kNuisFSI003,inputs_syst[45]);
  fparm->Update(conv::kNuisFSI004,inputs_syst[46]);
  fparm->Update(conv::kNuisFSI005,inputs_syst[47]);

  string oscparams = util::AsString(fparm->PhysicsHypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Plotting spectra for sample: " << util::AsString(s) << endl;
  cout << "- Fit param list: " << *fparm << endl;
  cout << *Opt::Instance() << endl;

  // Histogram booking
  const int nbins = util::NumErecoBins(s);
  const double * ereco_bin_edges = util::ErecoBinEdges(s);
  TH1D * hErecoTot = new TH1D("","",nbins,ereco_bin_edges);
  TH1D * hErecoMode[conv::knm]; //< contribution from each mode  
  const int nbins_coarse = util::NumVariableXsEtrueBins(s,conv::kSKmNuMuCCQE);
  const double * ereco_bin_edges_coarse = util::VariableXsEtrueBinEdges(s,conv::kSKmNuMuCCQE);
  TH1 * hErecoTot_coarse;
  TH1 * hErecoMode_coarse[conv::knm]; //< contribution from each mode

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // Generate prediction
  cout << "- Generating predictions..." << endl;
  skexp.Generate(*fparm,true);
  skexp.PrintSystematicVariations();

  // Get energy spectra
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) << endl;
    TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
    assert(h);
    hErecoMode[im] = new TH1D(*h);
    hErecoMode_coarse[im] = hErecoMode[im]->Rebin(nbins_coarse,Form("%s_coarse",hErecoMode[im]->GetName()),ereco_bin_edges_coarse);
  }
  hErecoTot = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  hErecoTot_coarse = hErecoTot->Rebin(nbins_coarse,Form("%s_coarse",hErecoTot->GetName()),ereco_bin_edges_coarse);

  //
  // Draw
  //

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TLegend * legend = new TLegend(0.6, 0.4, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
  if(input_mixing23 == 1.0) hErecoTot->GetYaxis()->SetRangeUser(0,5);
  hErecoTot->Draw();
  hErecoTot->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot->GetYaxis()->SetTitle("Number of events");
  THStack * stack = new THStack("stack","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack->Add(hErecoMode[im]);
  }
  stack->Draw("same");
  legend->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend->Draw();
  c->Update();
  c->SaveAs(Form("sk_ereco_spectrum-%f-%f.eps",input_mixing23,input_dm23sq));
  hErecoTot->GetXaxis()->UnZoom();

  //draw coarse binned histograms
  TCanvas * c_coarse = new TCanvas("c_coarse","",20,20,500,500);
  TLegend * legend_coarse = new TLegend(0.6, 0.4, 0.9, 0.9);
  legend_coarse->SetName("legend_coarse");
  legend_coarse->SetFillColor(0);
  legend_coarse->SetBorderSize(0);
  hErecoTot_coarse->GetXaxis()->SetRangeUser(0,7-conv::kASmallNum);
  hErecoTot_coarse->Draw();
  hErecoTot_coarse->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot_coarse->GetYaxis()->SetTitle("Number of events");
  THStack * stack_coarse = new THStack("stack_coarse","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack_coarse->Add(hErecoMode_coarse[im]);
  }
  stack_coarse->Draw("same");
  legend_coarse->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend_coarse->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend_coarse->Draw();
  c_coarse->Update();
  c_coarse->SaveAs(Form("sk_ereco_spectrum-%f-%f_coarse.eps",input_mixing23,input_dm23sq));


  //can save this to file &/or print it
  stringstream streamer;

  streamer << fixed;
  double ntot = hErecoTot->Integral();
  double ntot_coarse = hErecoTot->Integral();
  if(TMath::Abs(ntot - ntot_coarse) > conv::kASmallNum) {
    cerr << "Integrals are different between differently binned histograms "
	 << "50MeV: " << ntot
	 << "BANFF binning: " << ntot_coarse << endl;
    exit(1);
  }

  streamer << "Total number of events: " 
       << setw(10) << setfill(' ') << setprecision(6) << ntot << endl;
  const int nslices=6;
  double slice_max[nslices] = { 0.50, 0.75, 1.00, 1.25, 1.50, 2.00 };
  for(int islice=0; islice<nslices; islice++) {
    double esmax = slice_max[islice];
    streamer << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << esmax << " GeV): " 
         << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(esmax)) << endl;
  }
  for(int im = 0; im < conv::knm; im++) {
     conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
     if(m == conv::kSKmUndefined) continue;
     double nmode = hErecoMode[im]->Integral();
     streamer << "Contrib from mode: " 
          << setw(40) << setfill(' ') << util::AsStringTex(m) << " : " 
          << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
          << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  }

  //
  // Save histograms in a ROOT file
  //

  TFile f(Form("sk_ereco_spectrum-%s.root",oscparams.c_str()) ,"recreate");
  hErecoTot -> Write(util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str());
  hErecoTot_coarse->Write(Form("%s_coarse",util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str()));
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    hErecoMode[im] -> Write(util::AsStringCompact(s,m,"sk_","").c_str());
    hErecoMode_coarse[im] -> Write(Form("sk_%s_coarse",util::AsStringCompact(s,m).c_str()));
  }
  stack->Write();
  legend->Write();
  stack_coarse->Write();
  legend_coarse->Write();
  f.Close();

  //
  // Print-out numbers requested during the analysis review
  //

  if(print_spectrum) {  
    int nprec=5;
    streamer << fixed;

    streamer << endl;
    for(int ibin = 1; ibin <= hErecoTot_coarse->GetXaxis()->GetNbins(); ibin++) {
      streamer << setw(6) << setfill(' ') << ibin
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinLowEdge(ibin)
	       << " - "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinUpEdge(ibin)
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot_coarse->GetBinContent(ibin)
	       << endl;
    }//ibin

    streamer << endl;
    for(int ibin = 1; ibin <= hErecoTot->GetXaxis()->GetNbins(); ibin++) {
      streamer << setw(6) << setfill(' ') << ibin
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinLowEdge(ibin)
	       << " - "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinUpEdge(ibin)
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot->GetBinContent(ibin)
	       << endl;
    }//ibin
  }

  //actually print
  cout << streamer.str() << endl;
  
  // get size of streamer
  long size;
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);

  //write to file
  TString outname = Form("sk_ereco_spectrum-%s.txt",oscparams.c_str());
  ofstream outfile (outname);
  outfile.write (streamer.str().c_str(),size);
  outfile.close();

}
