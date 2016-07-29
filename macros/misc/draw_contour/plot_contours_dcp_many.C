
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

#include <vector>

using namespace std;
using namespace valor;
using namespace valor::conv;

const double CLfrac_1sigma = 0.6827;  // 1sigma% CL
const double CLfrac_90CL   = 0.9;  // 90% CL
const double CLfrac_95CL   = 0.95;  // 95% CL
const double CLfrac_2sigma = 0.9545;  // 2sigma% CL

void plot_contours_dcp_many
( 
 const bool SaveRootFile = true,
 const string tag = "dcp_run17c_react_Data",

 TString titleX = "#delta_{CP} (radians)",
 //TString titleX = "sin^{2}#theta_{23}",
 //TString titleX = "sin^{2}#theta_{13}",
 TString titleY = "-2lnL",

 const bool doBrazilian = false,
 
 const bool doGlobalMin  = true,
 const bool same_mh_color = true, // use same color for NH and IH
 const bool same_mh_style = true, // use same style for NH and IH
 
 const bool draw_preliminary = true,
 TString str_preliminary = "T2K Run1-7c preliminary",

 const bool do_fc = true,
 const bool draw_fc_err = false,
 const bool do_fc_conf_interv = false,

 const bool do_const = false,

 const bool draw_1sigmaCL = false,
 const bool draw_90CL     = true,
 const bool draw_2sigmaCL = true,
 const bool draw_3sigmaCL = false,
 const bool draw_4sigmaCL = false,
 const bool draw_5sigmaCL = false,

 const bool draw_reactors = false,
 
 const bool do_range_1sigma = false,
 const double critval = 1.,
 
 // Compare true chisq
 
 // Asimov A best-fit (global)
 const double truepar0 = kBadResult,
 const string truemh0 = "normal",

 // // run 1-6 best-fit (global)
 // const double truepar0 = -1.8535397,
 // const string truemh0 = "normal",

 // // Asimov A best-fit (global)
 // const double truepar0 = -1.601,
 // const string truemh0 = "normal",

 // // Asimov B best-fit (global)
 // const double truepar0 = -2.82743339, 
 // const string truemh0 = "inverted",

 // True point
 const double truepar1 = 0., //0., //-1.351,
 const string truemh1 = "normal", //"inverted",
 
 //double mybestfit = 3.14, 
 double mybestfit = kBadResult
  )
{

  
  //
  // Run 1-7b
  //

  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b_FLUXRUN16/Raj/fit_nh_dcpminuspi_normal_-1.571_Raj.root";
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b/fullchisq/config_run16fluxtun_margtoysrun17fluxtun/fit_fitmarg_fc_nh_10k_mpiover2_normal_-1.571.root";
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b/fullchisq/config_run16fluxtun_margtoysrun17fluxtun/fit_fitmarg_fc_nh_10k_m1p759_normal_-1.759.root";  
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b/fullchisq/config_run16fluxtun_margtoysrun17fluxtun/fit_fitmarg_fc_nh_10k_0_normal_0.000.root";
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b/fullchisq/config_run16fluxtun_margtoysrun17fluxtun/fit_fitmarg_fc_ih_10k_mpiover2_inverted_-1.571.root";
  
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b_FLUXRUN17/fullchisq/run17fluxtunig_all/fit_fitmarg_fc_nh_10k_fullchisq_run17fluxtun_0.000.root";  
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b_FLUXRUN17/fullchisq/run17fluxtunig_all/fit_fitmarg_fc_nh_10k_fullchisq_run17fluxtun_-1.571.root";  
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b_FLUXRUN17/fullchisq/run17fluxtunig_all/fit_fitmarg_fc_nh_10k_fullchisq_run17fluxtun_-1.759.root ";    
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17b_FLUXRUN17/fullchisq/run17fluxtunig_all/fit_fitmarg_fc_ih_10k_fullchisq_run17fluxtun_-1.571.root";  

  //
  // Run 1-7c
  //
  const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17c_FINAL/fullchisq/fit_fitmarg_fc_nh_10k_run17c_final_FullChisq_normal_-1.571.root";  
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17c_FINAL/fullchisq/fit_fitmarg_fc_nh_10k_run17c_final_FullChisq_normal_-1.791.root";  
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17c_FINAL/fullchisq/fit_fitmarg_fc_nh_10k_run17c_final_FullChisq_normal_0.000.root";  
  //const char *infilename_toys = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17c_FINAL/fullchisq/fit_fitmarg_fc_nh_10k_run17c_final_FullChisq_inverted_-1.571.root";  

  const int Nfiles = 1;
  const int Nmh = 2;
  char *infilename[Nfiles][Nmh] = {
    


    //////////////////////////////////////////////
    //                                          //
    //         Run 1-7c - Data - Final          // 
    //                                          //
    //////////////////////////////////////////////

    // deltaCP - Run 1-7c w/ reactors - Final
    {
      "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_dcp_react_nh_Data_Run17c_nexp_1_0.000_-999.000-job0-_0_0.root",
      "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_dcp_react_ih_Data_Run17c_nexp_1_0.000_-999.000-job0-_0_0.root"
    },

    // // dmsq - Run 1-7c w/ reactors - Final
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_dmsq_react_nh_Data_Run17c_nexp_1_0.003_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_dmsq_react_ih_Data_Run17c_nexp_1_0.003_-999.000-job0-_0_0.root"
    // },

    // // sinsqth23 - Run 1-7c w/ reactors - Final
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth23_react_nh_Data_Run17c_nexp_1_0.450_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth23_react_ih_Data_Run17c_nexp_1_0.450_-999.000-job0-_0_0.root"
    // },

    // // deltaCP - Run 1-7c T2K only - Final
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/datafit/dcp/data_fit_nh_a-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/datafit/dcp/data_fit_ih_a-_0_0.root"
    // },

    // // sinsqth13 - Run 1-7c T2K only - Final
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/datafit/ss13/data_fit_nh_a-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/datafit/ss13/data_fit_ih_a-_0_0.root"
    // },
    
    // // deltaCP - Run 1-7c w/ reactors - Ereco only - Final
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/"
    // },

    // // deltaCP - Run 1-7c w/ T2K only - Ereco only - Final
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17c/files/
    // },                           





    //////////////////////////////////////////////
    //                                          //
    //      Run 1-7c -Asimov - Run 1-7 Flux      // 
    //                                          //
    //////////////////////////////////////////////
    
    // // deltaCP - Run 1-7 Asimov A w/ reactors - final
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_dcp_react_nh_Asimov_A_Run17c_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_dcp_react_ih_Asimov_A_Run17c_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },

    // // deltaCP - Run 1-7 Asimov B w/ reactors - final
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_dcp_react_nh_Asimov_B_Run17c_nexp_1_0.000_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_dcp_react_ih_Asimov_B_Run17c_nexp_1_0.000_-999.000-job0-_0_0.root"
    // },

    // // deltaCP - Run 1-7 Asimov A T2K only - final
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/asimovs_run17c_Raj/dcp/dcp_nh_a.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/asimovs_run17c_Raj/dcp/dcp_ih_a.root"
    // },

    // // deltaCP - Run 1-7 Asimov B T2K only - final
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/asimovs_run17c_Raj/dcp/dcp_nh_b.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/asimovs_run17c_Raj/dcp/dcp_ih_b.root"
    // },
  
    
    // // deltaCP - FHC: 14.882 x 10^20 POT, RHC: 9.871 x 10^20 POT - Asimov A w/ reactors - final
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_dcp_react_nh_Asimov_A_FHC14p882_RHC9p871_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_dcp_react_ih_Asimov_A_FHC14p882_RHC9p871_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },
    // // deltaCP - FHC: 12.482 x 10^20 POT, RHC: 12.271 x 10^20 POT - Asimov A w/ reactors - final
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_dcp_react_nh_Asimov_A_FHC12p482_RHC12p271_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_dcp_react_ih_Asimov_A_FHC12p482_RHC12p271_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },
    // // deltaCP - FHC: 17.282 x 10^20 POT, RHC: 7.471 x 10^20 POT - Asimov A w/ reactors - final
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_dcp_react_nh_Asimov_A_FHC17p282_RHC7p471_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_dcp_react_ih_Asimov_A_FHC17p282_RHC7p471_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },







 



    //////////////////////////////////////////////
    //                                          //
    //      Run 1-7b - Data - Run 1-7 Flux       // 
    //                                          //
    //////////////////////////////////////////////

    // // deltaCP - Run 1-7 w/ reactors - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dcp_react_nh_Data_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dcp_react_ih_Data_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },

    // // dmsq - Run 1-7 w/ reactors - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dmsq_react_nh_Data_Run17_nexp_1_0.003_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dmsq_react_ih_Data_Run17_nexp_1_0.003_-999.000-job0-_0_0.root"
    // },

    // // sinsqth23 - Run 1-7 w/ reactors - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_sinsqth23_react_nh_Data_Run17_nexp_1_0.528_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_sinsqth23_react_ih_Data_Run17_nexp_1_0.528_-999.000-job0-_0_0.root"
    // },

    // // deltaCP - Run 1-7 T2K only - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/data_fit_dcp_t2konly_nh_run17.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/data_fit_dcp_t2konly_ih_run17.root"
    // },

    // // sinsqth13 - Run 1-7 T2K only - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/data_fit_sinsqth13_t2konly_nh_run17.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/data_fit_sinsqth13_t2konly_ih_run17.root"
    // },

    // // deltaCP - Run 1-7 T2K only - Flux Run 1-7 (CHECK)
    //  {
    //    "../TECH_NOTE_14_3_16/PROVA_FC/fit_fitmarg_joint_dcp_t2konly_nh_Data_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //    "../TECH_NOTE_14_3_16/PROVA_FC/fit_fitmarg_joint_dcp_t2konly_ih_Data_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root"
    //  },

    // // deltaCP - Run 1-7 w/ reactors - Ereco only - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dcp_react_nh_Data_Run17_ErecoOnly_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dcp_react_ih_Data_Run17_ErecoOnly_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },

    // // deltaCP - Run 1-7 w/ T2K only - Ereco only - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dcp_t2konly_nh_Data_Run17_ErecoOnly_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dcp_t2konly_ih_Data_Run17_ErecoOnly_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },                           
                                                                         
                                                                                                                                             


    //////////////////////////////////////////////
    //                                          //
    //      Run 1-7b -Asimov - Run 1-7 Flux      // 
    //                                          //
    //////////////////////////////////////////////

    // // deltaCP - Run 1-7 Asimov A T2K only - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/asimov_fit_nh_a_dcp_t2konly_run17.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/asimov_fit_ih_a_dcp_t2konly_run17.root"
    // },

    // // deltaCP - Run 1-7 Asimov B T2K only - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/asimov_fit_nh_b_dcp_t2konly_run17.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/asimov_fit_ih_b_dcp_t2konly_run17.root"
    // },
    
    // // deltaCP - Run 1-7 Asimov A w/ reactors - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_dcp_react_nh_Asimov_A_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_dcp_react_ih_Asimov_A_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },

    // // deltaCP - Run 1-7 Asimov B w/ reactors - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_dcp_react_nh_Asimov_B_Run17_nexp_1_0.000_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_dcp_react_ih_Asimov_B_Run17_nexp_1_0.000_-999.000-job0-_0_0.root"
    // },
  
    // // // deltaCP - Run 1-7 Asimov A w/ reactors - Flux Run 1-6
    //  {
    //    "toyfits_dcp_datanh-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
    //    "toyfits_dcp_dataih-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root"
    //  },
    // // {
    // //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun16/files/fit_fitmarg_joint_dcp_react_nh_Asimov_A_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root",
    // //   "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun16/files/fit_fitmarg_joint_dcp_react_ih_Asimov_A_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // // }
    
    
    //////////////////////////////
    //                          //
    //      Run 1-6 -Asimov     // 
    //                          //
    //////////////////////////////

    // // deltaCP - Run 1-6 Asimov A w/ reactors
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN16/files/valor/fit_fitmarg_joint_dcp_react_nh_Run16_Asimov_A_nexp_1_-1.601_-999.000-job0-_0_0_TH1D.root", // NH 
    //   "../TECH_NOTE_14_3_16/SENSI/RUN16/files/valor/fit_fitmarg_joint_dcp_react_ih_Run16_Asimov_A_nexp_1_-1.601_-999.000-job0-_0_0_TH1D.root"  // IH
    //   // "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_nh_A_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH 
    //   // "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_ih_A_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
    //   //"fit_PROVA_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // }, 
    
    // // deltaCP - FHC: 14.4e20, RHC: 7.47e20   - Asimov A w/ reactors
    // {
    //   "fit_fitmarg_joint_dcp_react_nh_fhc_14p4_rhc_7p47_Asimov_A_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "fit_fitmarg_joint_dcp_react_ih_fhc_14p4_rhc_7p47_Asimov_A_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },
    // // deltaCP - FHC: 10.07e20, RHC: 11.17e20 - Asimov A w/ reactors
    // {
    //   "fit_fitmarg_joint_dcp_react_nh_fhc_10p07_rhc_11p17_Asimov_A_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "fit_fitmarg_joint_dcp_react_ih_fhc_10p07_rhc_11p17_Asimov_A_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },
    // // deltaCP - FHC: 7.00e20, RHC: 14.87e20  - Asimov A w/ reactors
    // {
    //   "fit_fitmarg_joint_dcp_react_nh_fhc_7p00_rhc_14p87_Asimov_A_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "fit_fitmarg_joint_dcp_react_ih_fhc_7p00_rhc_14p87_Asimov_A_nexp_1_-1.601_-999.000-job0-_0_0.root",
    // },

    // // deltaCP - Run 1-6 Asimov B w/ reactors
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/Run16/files/valor/fit_fitmarg_joint_dcp_react_nh_Run16_Asimov_B_nexp_1_0.000_-999.000-job0-_0_0_TH1D.root", // NH 
    //   "../TECH_NOTE_14_3_16/SENSI/Run16/files/valor/fit_fitmarg_joint_dcp_react_ih_Run16_Asimov_B_nexp_1_0.000_-999.000-job0-_0_0_TH1D.root"  // IH
    // //   "../TECH_NOTE_14_3_16/SENSI/Run16/files/valor/fit_fitmarg_joint_dcp_react_nh_B_nexp_1_0.000_-999.000-job0-_0_0.root", // NH 
    // //  "../TECH_NOTE_14_3_16/SENSI/Run16/files/valor/fit_fitmarg_joint_dcp_react_ih_B_nexp_1_0.000_-999.000-job0-_0_0.root"  // IH
    // }
    // ,

    // // deltaCP - FHC: 14.4e20, RHC: 7.47e20   - Asimov A w/ reactors
    // {
    //   "fit_fitmarg_joint_dcp_react_nh_fhc_14p4_rhc_7p47_Asimov_B_nexp_1_0.000_-999.000-job0-_0_0.root",
    //   "fit_fitmarg_joint_dcp_react_ih_fhc_14p4_rhc_7p47_Asimov_B_nexp_1_0.000_-999.000-job0-_0_0.root"
    // },
    // // deltaCP - FHC: 10.07e20, RHC: 11.17e20 - Asimov A w/ reactors
    // {
    //   "fit_fitmarg_joint_dcp_react_nh_fhc_10p07_rhc_11p17_Asimov_B_nexp_1_0.000_-999.000-job0-_0_0.root",
    //   "fit_fitmarg_joint_dcp_react_ih_fhc_10p07_rhc_11p17_Asimov_B_nexp_1_0.000_-999.000-job0-_0_0.root"
    // },
    // // deltaCP - FHC: 7.00e20, RHC: 14.87e20  - Asimov A w/ reactors
    // {
    //   "fit_fitmarg_joint_dcp_react_nh_fhc_7p00_rhc_14p87_Asimov_B_nexp_1_0.000_-999.000-job0-_0_0.root",
    //   "fit_fitmarg_joint_dcp_react_ih_fhc_7p00_rhc_14p87_Asimov_B_nexp_1_0.000_-999.000-job0-_0_0.root",
    // },


    // // deltaCP - Run 1-6 Asimov A w/o reactors
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN16/files/valor/toyfits_dcp_anh-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/SENSI/RUN16/files/valor/toyfits_dcp_aih-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root"  // IH
    // },
    
    // // deltaCP - Run 1-6 Asimov B w/o reactors
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN16/files/valor/toyfits_dcp_cnh-sinsqtheta[12_0.304,13_0.022,23_0.450],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/SENSI/RUN16/files/valor/toyfits_dcp_cih-sinsqtheta[12_0.304,13_0.022,23_0.450],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root"  // IH
    // }
    // ,

    // // deltaCP - FHC: 12.5e21, RHC: 12.5e21 w/ reactors - Asimov A w/ reactors
    // {
    //   "fit_fitmarg_joint_dcp_react_nh_Asimov_A_Both12p5E21POT_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "fit_fitmarg_joint_dcp_react_ih_Asimov_A_Both12p5E21POT_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },


    
    //////////////////////////////
    //                          //
    //      Run 1-6 - Data      // 
    //                          //
    //////////////////////////////

    // // deltaCP - Run 1-6 Data w/ reactors
    // {
    //   "../TECH_NOTE_14_3_16/DATA/RUN16/files/fit_fitmarg_joint_dcp_react_nh_Data_Run16_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/DATA/RUN16/files/fit_fitmarg_joint_dcp_react_ih_Data_Run16_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
    //   //"fit_fitmarg_fc_50k_normal_nexp_1_-999.000_-999.000-job0-_0_0.root"
    //   //"fit_fit_marg_PROVA_dcp_nexp_1_0.000_-999.000-job0-_0_0.root",
    //   //"fit_fit_marg_PROVA_dcp_nexp_1_0.000_-999.000-job0-_0_0.root"
    // }
    // ,

    // // deltaCP - Run 1-6 Data w/ reactors - e-like 15 bins only for RHC
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_nh_Data_Run16_elike_RHConly_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_ih_Data_Run16_elike_RHConly_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
    // }
    // ,
  
    // // deltaCP - Run 1-6 Data T2K only
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/toyfits_dcp_datanh-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/toyfits_dcp_dataih-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root"  // IH
    // }
    // ,
  
    // // deltaCP - Run 1-6 Data w/ reactors - FHC only
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_nh_Data_Run16_FHC_nexp_1_-1.601_-999.000-job0-_0_0.root", // FHC - NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_ih_Data_Run16_FHC_nexp_1_-1.601_-999.000-job0-_0_0.root", // FHC - IH 
    // }
    // ,
    // // deltaCP - Run 1-6 Data w/ reactors - RHC only
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_nh_Data_Run16_RHC_nexp_1_-1.601_-999.000-job0-_0_0.root", // RHC - NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_ih_Data_Run16_RHC_nexp_1_-1.601_-999.000-job0-_0_0.root"  // RHC - IH
    // }

    // // sinsqth13 - Run 1-6 Data T2K only
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth13_t2konly_nh_Data_Run16_nexp_1_0.022_-999.000-job0-_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth13_t2konly_ih_Data_Run16_nexp_1_0.022_-999.000-job0-_0_0.root", // IH 
    // }
    // ,
    
    // // sinsqth23 - Run 1-6 Data w/ reactors
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23_react_nh_Data_Run16_nexp_1_0.528_-999.000-job0-_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23_react_ih_Data_Run16_nexp_1_0.528_-999.000-job0-_0_0.root", // IH 
    // }
    // ,
    
    // // dmsq - Run 1-6 Data w/ reactors
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dmsq_react_nh_Data_Run16_nexp_1_0.003_-999.000-job0-_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dmsq_react_ih_Data_Run16_nexp_1_0.003_-999.000-job0-_0_0.root", // IH 
    // }
    // ,
    
    // // deltaCP - Run 1-6 Data w/ reactors - rate-only
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_nh_Data_Run16_rateonly_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_ih_Data_Run16_rateonly_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
    // }
    // ,

  };
  
 
  // File names of Feldman-Cousins
  char *fcname = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17c_FLUXMERGE/FCgraphs_fluxmerge_dcp_3p142.root";
  //char *fcname = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17c_FLUXRUN17b/FCgraphs_Run17c_fluxRun17b.root";
  //char *fcname = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN17_FLUXRUN16/FCgraphs_10k_Run17_FluxRun16.root";
  //char *fcname = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/RUN16/FCgraphs_p8k.root";
  //char *fcname = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/FCgraphs_4k.root";

  //char *fcname = "FCgraphs.root";

  char *latexname[Nfiles][Nmh] = {
    {"Normal Hierarchy","Inverted Hierarchy"},
        
    //{"Run 1-7c - NH", "Run 1-7c - IH"},
    // {"Run 1-7b - NH", "Run 1-7b - IH"},

    // {"E_{reco} Vs #theta - NH", "E_{reco} Vs #theta - IH"},
    // {"E_{reco} only - NH", "E_{reco} only - IH"},
    
    // {"#nu: 7.48e20  +  #bar{#nu}: 7.47e20 (Run 1-7c)", "IH"},
    // {"#nu: 14.9e20  +  #bar{#nu}: 9.9e20 ", "IH"},
    // {"#nu: 12.5e20  +  #bar{#nu}: 12.3e20","IH"},
    // {"#nu: 17.3e20  +  #bar{#nu}: 7.5e20 ", "IH"}
    
    // {"Run 1-7c - NH"},
    // {"14.882e20 FHC + 9.871e20 RHC - NH"}, 
    // {"12.482e20 FHC + 12.271e20 RHC - NH"}, 
    // {"17.282e20 FHC + 7.471e20 RHC - NH"}, 
  };
  

  char name[200];
  
  TFile *filein[Nfiles][Nmh];
  TTree *treein[Nfiles][Nmh];

  TH1D *hchisq1D[Nfiles][Nmh];
  TH1D *hdeltachisq1D[Nfiles][Nmh];
  TH1D *hdeltachisq1D_GlobMin[Nfiles][Nmh];

  double bestfit_x = -999.;
  double bestfit_mh = -999.;

  TCanvas *cChisq; //[Nfiles][Nmh];

  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){      
    
      cout << endl;
      cout << "Get file " << infilename[ifile][imh] << endl;
      cout << endl;

      filein[ifile][imh] = new TFile(infilename[ifile][imh],"READ");
      if(!filein[ifile][imh]->IsOpen()){
	cerr << endl;
	cerr << "File " << infilename[ifile][imh] << " is not open!" << endl;
	cerr << endl;
	exit(1);
      }
      filein[ifile][imh]->cd();
      treein[ifile][imh] = (TTree*)filein[ifile][imh]->Get("fitnt");
      
      sprintf(name,"hchisq1D_file%d_mh%d",ifile,imh);
      hchisq1D[ifile][imh] = new TH1D(name,name,1,0,1);
      treein[ifile][imh]->SetBranchAddress("hchisq1D", &hchisq1D[ifile][imh] );
      
      sprintf(name,"hdeltachisq1D_file%d_mh%d",ifile,imh);
      hdeltachisq1D[ifile][imh] = new TH1D(name,name,1,0,1);
      treein[ifile][imh]->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D[ifile][imh] );

      // Get event
      treein[ifile][imh]->GetEntry(0);
      
      // Print chisq distributions
      cout << std::setprecision(9);
      cout << endl;
      cout << "Chisq file " << ifile << " - mh = " << imh << " :" << endl;
      cout << "# of entries = " << hchisq1D[ifile][imh]->GetEntries() << endl;
      for(int i=0;i<hchisq1D[ifile][imh]->GetNbinsX();i++){
	cout << hchisq1D[ifile][imh]->GetBinContent(i+1) << " ";
      }
      cout << endl;
      cout << "DeltaChisq: # of entries = " << hdeltachisq1D[ifile][imh]->GetEntries() << endl;
      for(int i=0;i<hdeltachisq1D[ifile][imh]->GetNbinsX();i++){
      	cout << hdeltachisq1D[ifile][imh]->GetBinContent(i+1) << " ";
      }
      cout << endl;  
    
      sprintf(name,"hdeltachisq1D_GlobMin_file%d_mh%d",ifile,imh);
      hdeltachisq1D_GlobMin[ifile][imh] = (TH1D*)hchisq1D[ifile][imh]->Clone(name);
      hdeltachisq1D_GlobMin[ifile][imh]->Reset();

    } // end loop Nmh

    
    // Shift the chisq wrt to the global minimum
    
    double ChisqMin_0 = hchisq1D[ifile][0]->GetBinContent(hchisq1D[ifile][0]->GetMinimumBin());

    double ChisqMin_1 = 999999;
    //if(hchisq1D[ifile][1]) 
    ChisqMin_1 = hchisq1D[ifile][1]->GetBinContent(hchisq1D[ifile][1]->GetMinimumBin());
    double ChisqMin_Glob = ChisqMin_0; 
    if( ChisqMin_1 < ChisqMin_Glob ) ChisqMin_Glob = ChisqMin_1;
    
    int XminBin_0 = hchisq1D[ifile][0]->GetMinimumBin();
    double Xmin_0 = hchisq1D[ifile][0]->GetBinCenter(XminBin_0);
    
    int XminBin_1 = -999;
    double Xmin_1 = -999;   
    if(hchisq1D[ifile][1]){
      XminBin_1 = hchisq1D[ifile][1]->GetMinimumBin();
      Xmin_1 = hchisq1D[ifile][1]->GetBinCenter(XminBin_1);
    }

    cout << endl;
    cout << "Minimum for " << latexname[ifile][0] << " :" << endl;
    cout << " - X = " << Xmin_0 << endl;
    cout << " - ChisqMin " << latexname[ifile][0] << " = " << ChisqMin_0 << endl;

    if(hchisq1D[ifile][1]){
      cout << "Minimum for " << latexname[ifile][1] << " :" << endl;
      cout << " - X = " << Xmin_1 << endl;
      cout << " - ChisqMin " << latexname[ifile][1] << " = " << ChisqMin_1 << endl;
      cout << endl;
    }
    if(ChisqMin_0<ChisqMin_1){
      cout << "Best-fit = " << Xmin_0 << endl;
    }
    else{
      cout << "Best-fit = " << Xmin_1 << endl;
    }
    cout << "ChisqMin Glob = " << ChisqMin_Glob << endl;
    cout << endl;
    
    
    for(int imh=0;imh<Nmh;imh++){      
      double ChisqMin_Curr = hchisq1D[ifile][imh]->GetBinContent(hchisq1D[ifile][imh]->GetMinimumBin());
      for(int i=0;i<hdeltachisq1D_GlobMin[ifile][imh]->GetNbinsX();i++){
	double cont = hchisq1D[ifile][imh]->GetBinContent(i+1); 
	if(doGlobalMin) hdeltachisq1D_GlobMin[ifile][imh]->SetBinContent(i+1,cont-ChisqMin_Glob);
	else            hdeltachisq1D_GlobMin[ifile][imh]->SetBinContent(i+1,cont-ChisqMin_Curr);
      }
  
      // for(int i=0;i<hdeltachisq1D_GlobMin[ifile][imh]->GetNbinsX();i++){
      // 	cout << hdeltachisq1D_GlobMin[ifile][imh]->GetBinContent(i+1) << " ";
      // }
      // cout << endl;

    }

  } // end loop Nfiles
  

  


  // Calculate 1 sigma ranges

  if(do_range_1sigma){
    cout << endl;
    cout << "Calculate range for critical value: " << critval << endl;
    cout << endl;
    
    for(int ifile=0;ifile<Nfiles;ifile++){
      for(int imh=0;imh<Nmh;imh++){      
	
	cout << endl;
	cout << infilename[ifile][imh] << endl;
	
	cout << "mybestfit = " << mybestfit << endl;
	double Status_range = CalcBestFitErrorsMINOS(hdeltachisq1D_GlobMin[ifile][imh],critval,mybestfit);
	
	if(Status_range==conv::kBadResult){
	  cout << "Range not found" << endl;
	}
	
      }
    }
  }


  // Marginalize deltaCP
  
  cout << endl;
  cout << "Get chisq(NH) - chisq(IH)" << endl;
  
  double chisqMH[Nfiles][Nmh] = {conv::kBadResult};

  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){      

      double like_ave = 0.;
      for(int i=0;i<hchisq1D[ifile][imh]->GetNbinsX();i++){
	double chisq = hchisq1D[ifile][imh]->GetBinContent(i+1);
	double like = TMath::Exp(-chisq / 2.);
	like_ave += like;
      }
      like_ave /= hchisq1D[ifile][imh]->GetNbinsX();
      
      chisqMH[ifile][imh] = - 2 * TMath::Log(like_ave);

      cout << " - " << latexname[ifile][imh] << " :" << endl;
      cout << "    chisq = " << chisqMH[ifile][imh] << endl;

    }
    cout << "Delta Chisq = " << chisqMH[ifile][0] - chisqMH[ifile][1] << endl;
  }
  cout << endl; 
  
  // Init histogram

  int NbinsX = hdeltachisq1D_GlobMin[0][0]->GetNbinsX(); 
  int NbinsY = hdeltachisq1D_GlobMin[0][0]->GetNbinsY();
  double lowedgeX = hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1);
  double upedgeX = hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(NbinsX);
  char *TitleX = hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetTitle();
  char *TitleY = hdeltachisq1D_GlobMin[0][0]->GetYaxis()->GetTitle();
  TH1D *hInit = new TH1D("hInit","hInit",NbinsX,lowedgeX,upedgeX);
  hInit->GetXaxis()->SetTitle(TitleX);
  hInit->GetYaxis()->SetTitle(TitleY);



  // Get file of Feldman-Cousins
  
  TGraph *gFC90CL_nh;
  TGraph *gFC90CL_ih;
  TGraph *gFC2sigma_nh;
  TGraph *gFC2sigma_ih;
  TFile *fin_fc_nh;
  TFile *fin_fc_ih;
  
  if(do_fc){

    fin_fc = new TFile(fcname,"READ");
    
    gFC90CL_nh      = (TGraph*) fin_fc->Get("gFC90CL_nh");
    gFC90CL_nh_p1   = (TGraph*) fin_fc->Get("gFC90CL_nh_p1");
    gFC90CL_nh_m1   = (TGraph*) fin_fc->Get("gFC90CL_nh_m1");
    gFC2sigma_nh    = (TGraph*) fin_fc->Get("gFC2sigma_nh");
    gFC2sigma_nh_p1 = (TGraph*) fin_fc->Get("gFC2sigma_nh_p1");
    gFC2sigma_nh_m1 = (TGraph*) fin_fc->Get("gFC2sigma_nh_m1");
    gFC90CL_ih      = (TGraph*) fin_fc->Get("gFC90CL_ih");
    gFC90CL_ih_p1   = (TGraph*) fin_fc->Get("gFC90CL_ih_p1");
    gFC90CL_ih_m1   = (TGraph*) fin_fc->Get("gFC90CL_ih_m1");
    gFC2sigma_ih    = (TGraph*) fin_fc->Get("gFC2sigma_ih");
    gFC2sigma_ih_p1 = (TGraph*) fin_fc->Get("gFC2sigma_ih_p1");
    gFC2sigma_ih_m1 = (TGraph*) fin_fc->Get("gFC2sigma_ih_m1");

    //fin_fc_nh = new TFile(fcname_nh,"READ");
    //fin_fc_ih = new TFile(fcname_ih,"READ");    
    //gFC90CL_nh = (TGraph*)fin_fc_nh->Get("gFC90CL_nh");
    //gFC2sigma_nh = (TGraph*)fin_fc_nh->Get("gFC2sigma_nh");
    //gFC90CL_ih = (TGraph*)fin_fc_ih->Get("gFC90CL_ih");
    //gFC2sigma_ih = (TGraph*)fin_fc_ih->Get("gFC2sigma_ih");

    
    if(do_fc_conf_interv){
      gFC90CL_nh   ->SetLineColor(kBlack);
      gFC90CL_nh_p1->SetLineColor(kBlack);
      gFC90CL_nh_m1->SetLineColor(kBlack);
      gFC90CL_ih   ->SetLineColor(kRed);
      gFC90CL_ih_p1->SetLineColor(kRed);
      gFC90CL_ih_m1->SetLineColor(kRed);
    }
    else{
      gFC90CL_nh   ->SetLineColor(kBlue);
      gFC90CL_nh_p1->SetLineColor(kBlue);
      gFC90CL_nh_m1->SetLineColor(kBlue);
      gFC90CL_ih   ->SetLineColor(kBlue);
      gFC90CL_ih_p1->SetLineColor(kBlue);
      gFC90CL_ih_m1->SetLineColor(kBlue);
    }
    
    gFC2sigma_nh   ->SetLineColor(kOrange-4);
    gFC2sigma_nh_p1->SetLineColor(kOrange-4);
    gFC2sigma_nh_m1->SetLineColor(kOrange-4);
    gFC2sigma_ih   ->SetLineColor(kOrange-4);
    gFC2sigma_ih_p1->SetLineColor(kOrange-4);
    gFC2sigma_ih_m1->SetLineColor(kOrange-4);

    gFC90CL_nh   ->SetMarkerSize(0.5);
    gFC90CL_nh_p1->SetMarkerSize(0.5);
    gFC90CL_nh_m1->SetMarkerSize(0.5);
    gFC90CL_ih   ->SetMarkerSize(0.5);
    gFC90CL_ih_p1->SetMarkerSize(0.5);
    gFC90CL_ih_m1->SetMarkerSize(0.5);
    gFC2sigma_nh   ->SetMarkerSize(0.5);
    gFC2sigma_nh_p1->SetMarkerSize(0.5);
    gFC2sigma_nh_m1->SetMarkerSize(0.5);
    gFC2sigma_ih   ->SetMarkerSize(0.5);
    gFC2sigma_ih_p1->SetMarkerSize(0.5);
    gFC2sigma_ih_m1->SetMarkerSize(0.5);
    
    gFC90CL_nh     ->SetLineWidth(3);
    gFC90CL_nh_p1  ->SetLineWidth(3);
    gFC90CL_nh_m1  ->SetLineWidth(3);
    gFC90CL_ih     ->SetLineWidth(3);
    gFC90CL_ih_p1  ->SetLineWidth(3);
    gFC90CL_ih_m1  ->SetLineWidth(3);
    gFC2sigma_nh   ->SetLineWidth(3);
    gFC2sigma_nh_p1->SetLineWidth(3);
    gFC2sigma_nh_m1->SetLineWidth(3);
    gFC2sigma_ih   ->SetLineWidth(3);
    gFC2sigma_ih_p1->SetLineWidth(3);
    gFC2sigma_ih_m1->SetLineWidth(3);

    if(do_fc_conf_interv){
      gFC90CL_nh     ->SetLineStyle(2);
      gFC90CL_nh_p1  ->SetLineStyle(2);
      gFC90CL_ih_m1  ->SetLineStyle(2);
      gFC2sigma_nh   ->SetLineStyle(2);
      gFC2sigma_nh_p1->SetLineStyle(2);
      gFC2sigma_nh_m1->SetLineStyle(2);
      gFC90CL_ih     ->SetLineStyle(2);
      gFC90CL_ih_p1  ->SetLineStyle(2);
      gFC90CL_ih_m1  ->SetLineStyle(2);
      gFC2sigma_ih   ->SetLineStyle(2);
      gFC2sigma_ih_p1->SetLineStyle(2);
      gFC2sigma_ih_m1->SetLineStyle(2);
    }
  }





  //
  // Compute CL intervals at 90%CL
  //
  
  TLine *min90CL_nh;
  TLine *max90CL_nh;
  TLine *min90CL_nh_2;
  TLine *max90CL_nh_2;
  TLine *min90CL_ih;
  TLine *max90CL_ih;
  TBox *hatch_min90CL_nh;
  TBox *hatch_max90CL_nh;
  TBox *hatch_min90CL_nh_2;
  TBox *hatch_max90CL_nh_2;
  TBox *hatch_min90CL_ih;
  TBox *hatch_max90CL_ih;

  if(do_fc && do_fc_conf_interv){
    
    cout << endl;
    cout << "/////////// Compute confidence intervals w/ Feldman-Cousins at 90%CL ///////////" << endl;
    
    TH1D *hdeltachisq1D_shifted_fc[Nfiles][Nmh]; // shifted wrt Feldman-Cousins critical values
    double lim_min_nh = 0;
    double lim_max_nh = 0;
    double lim_min_nh_2 = 0;
    double lim_max_nh_2 = 0;
    double lim_min_ih = 0;
    double lim_max_ih = 0;

    // Normal Hierarchy interval
    
    cout << endl;
    cout << "Normal hierarchy:" << endl;
    
    hdeltachisq1D_shifted_fc[0][0] = (TH1D*)hdeltachisq1D_GlobMin[0][0]->Clone(TString::Format("hdeltachisq1D_shifted_fc_file0_mh0"));
    
    for(int ibin=0;ibin<hdeltachisq1D_shifted_fc[0][0]->GetNbinsX();ibin++){
      double bincenter = hdeltachisq1D_shifted_fc[0][0]->GetBinCenter(ibin+1);
      double bincont = hdeltachisq1D_shifted_fc[0][0]->GetBinContent(ibin+1);
      double bincont_new = bincont - gFC90CL_nh->Eval(bincenter); // bincont=0 --> boundary of CL
      hdeltachisq1D_shifted_fc[0][0]->SetBinContent(ibin+1,bincont_new);
    }
    
    //int bin_bf = hdeltachisq1D_GlobMin[0][0]->GetMinimumBin();
    //double dcp_bestfit = hdeltachisq1D_GlobMin[0][0]->GetBinCenter(bin_bf);
    //double dcp_bestfit_cont = hdeltachisq1D_GlobMin[0][0]->GetBinContent(bin_bf);
    int bin_bf = hdeltachisq1D_shifted_fc[0][0]->GetMinimumBin();
    double dcp_bestfit = hdeltachisq1D_shifted_fc[0][0]->GetBinCenter(bin_bf);
    double dcp_bestfit_cont = hdeltachisq1D_shifted_fc[0][0]->GetBinContent(bin_bf);
    //cout << "dcp_bestfit = " << dcp_bestfit << endl;
    //cout << "dcp_bestfit_cont = " << dcp_bestfit_cont << endl;
    
    double Status_range = CalcBestFitErrorsMINOS_FC(hdeltachisq1D_shifted_fc[0][0],0.,dcp_bestfit,lim_min_nh,lim_max_nh);
    if(Status_range==conv::kBadResult){
      cout << "Range not found" << endl;
    }
    //cout << "lim_min_nh = " << lim_min_nh << endl;
    //cout << "lim_max_nh = " << lim_max_nh << endl;

 
    
    // // Normal Hierarchy interval (second part)
    // cout << endl;
    // cout << "Normal hierarchy (second part):" << endl;
    // bin_bf = hdeltachisq1D_shifted_fc[0][0]->FindBin(+TMath::Pi());
    // dcp_bestfit = hdeltachisq1D_shifted_fc[0][0]->GetBinCenter(bin_bf);
    // dcp_bestfit_cont = hdeltachisq1D_shifted_fc[0][0]->GetBinContent(bin_bf);
    // Status_range = CalcBestFitErrorsMINOS_FC(hdeltachisq1D_shifted_fc[0][0],0.,dcp_bestfit,lim_min_nh_2,lim_max_nh_2);
    // if(Status_range==conv::kBadResult){
    //   cout << "Range not found" << endl;
    // }
    // //cout << "lim_min_nh_2 = " << lim_min_nh_2 << endl;
    // //cout << "lim_max_nh_2 = " << lim_max_nh_2 << endl;


   

    // Inverted Hierarchy interval
    
    cout << endl;
    cout << "Inverted hierarchy:" << endl;
    cout << endl;
    
    hdeltachisq1D_shifted_fc[0][1] = (TH1D*)hdeltachisq1D_GlobMin[0][1]->Clone(TString::Format("hdeltachisq1D_shifted_fc_file0_mh1"));
    
    for(int ibin=0;ibin<hdeltachisq1D_shifted_fc[0][1]->GetNbinsX();ibin++){
      double bincenter = hdeltachisq1D_shifted_fc[0][1]->GetBinCenter(ibin+1);
      double bincont = hdeltachisq1D_shifted_fc[0][1]->GetBinContent(ibin+1);
      double bincont_new = bincont - gFC90CL_ih->Eval(bincenter); // bincont=0 --> boundary of CL
      hdeltachisq1D_shifted_fc[0][1]->SetBinContent(ibin+1,bincont_new);
    }
    
    //int bin_bf = hdeltachisq1D_GlobMin[0][1]->GetMinimumBin();
    //double dcp_bestfit = hdeltachisq1D_GlobMin[0][1]->GetBinCenter(bin_bf);
    //double dcp_bestfit_cont = hdeltachisq1D_GlobMin[0][1]->GetBinContent(bin_bf);
    int bin_bf = hdeltachisq1D_shifted_fc[0][1]->GetMinimumBin();
    double dcp_bestfit = hdeltachisq1D_shifted_fc[0][1]->GetBinCenter(bin_bf);
    double dcp_bestfit_cont = hdeltachisq1D_shifted_fc[0][1]->GetBinContent(bin_bf);
    //cout << "dcp_bestfit = " << dcp_bestfit << endl;
    //cout << "dcp_bestfit_cont = " << dcp_bestfit_cont << endl;
    
    double Status_range = CalcBestFitErrorsMINOS_FC(hdeltachisq1D_shifted_fc[0][1],0.,dcp_bestfit,lim_min_ih,lim_max_ih);
    if(Status_range==conv::kBadResult){
      cout << "Range not found" << endl;
    }
    //cout << "lim_min_ih = " << lim_min_ih << endl;
    //cout << "lim_max_ih = " << lim_max_ih << endl;


    // // 
    // // Useful to check whether the method is correct
    // //
    // TCanvas *cdeltachisq1D_shifted_fc_nh = new TCanvas("cdeltachisq1D_shifted_fc_nh","cdeltachisq1D_shifted_fc_nh");
    // hdeltachisq1D_shifted_fc[0][0]->Draw("histo text");
    // TLine *prova_interp = new TLine(-3.04734487,0.111296385,-3.01592895,-0.00225649251); // bins before and after y=0
    // TLine *prova_0 = new TLine(-3.14,0.,+3.14,0.); // line at y=0
    // TLine *prova_limit = new TLine(lim_min_nh,-10.,lim_min_nh,10.); // interval limit
    // prova_interp->Draw();
    // prova_0->Draw();
    // prova_limit->Draw();

    double max_line_y = 25; //10; //25;
    min90CL_nh = new TLine(lim_min_nh,0.,lim_min_nh,max_line_y);
    max90CL_nh = new TLine(lim_max_nh,0.,lim_max_nh,max_line_y);
    // min90CL_nh_2 = new TLine(lim_min_nh_2,0.,lim_min_nh_2,max_line_y);
    // //max90CL_nh_2 = new TLine(lim_max_nh_2,0.,lim_max_nh_2,max_line_y);
    min90CL_ih = new TLine(lim_min_ih,0.,lim_min_ih,max_line_y);
    max90CL_ih = new TLine(lim_max_ih,0.,lim_max_ih,max_line_y);

    double min_box = lim_min_nh+0.25;
    if(min_box < -TMath::Pi()) min_box = -TMath::Pi();
    hatch_min90CL_nh = new TBox(lim_min_nh,0.,min_box,max_line_y);
    double max_box = lim_max_nh-0.25;
    if(min_box > +TMath::Pi()) max_box = +TMath::Pi();
    hatch_max90CL_nh = new TBox(max_box,0.,lim_max_nh,max_line_y);
    
    min_box = lim_min_nh_2+0.25;
    if(min_box < -TMath::Pi()) min_box = -TMath::Pi();
    hatch_min90CL_nh_2 = new TBox(lim_min_nh_2,0.,min_box,max_line_y);
    // max_box = lim_max_nh_2-0.25;
    // if(min_box > +TMath::Pi()) max_box = +TMath::Pi();
    // hatch_max90CL_nh_2 = new TBox(max_box,0.,lim_max_nh_2,max_line_y);

    min_box = lim_min_ih+0.25;
    if(min_box < -TMath::Pi()) min_box = -TMath::Pi();
    hatch_min90CL_ih = new TBox(lim_min_ih,0.,min_box,max_line_y);
    max_box = lim_max_ih-0.25;
    if(min_box > +TMath::Pi()) max_box = +TMath::Pi();
    hatch_max90CL_ih = new TBox(max_box,0.,lim_max_ih,max_line_y);
    
  }
  







  /////////////////////////////////////////////
  //                                         //
  //            Draw histograms              //
  //                                         //
  /////////////////////////////////////////////


  // Chisq
  double Max_chisq = 0.;
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      if(hchisq1D[ifile][imh]->GetMaximum()>Max_chisq) 
	Max_chisq = hchisq1D[ifile][imh]->GetMaximum();
    }
  }
  sprintf(name,"cChisq1D");
  cChisq = new TCanvas(name,name);
  hInit->SetMaximum(Max_chisq*1.2);
  hInit->DrawCopy();
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      hchisq1D[ifile][imh]->SetMaximum(Max_chisq*1.2);
      if(same_mh_style) hchisq1D[ifile][imh]->SetLineStyle(1);
      else hchisq1D[ifile][imh]->SetLineStyle(imh+1);
      hchisq1D[ifile][imh]->SetLineWidth(2);
      if(same_mh_color) hchisq1D[ifile][imh]->SetLineColor(ifile+1);
      else hchisq1D[ifile][imh]->SetLineColor(ifile+imh+1);
      if((ifile+imh)==0) hchisq1D[ifile][imh]->Draw("same");
      else hchisq1D[ifile][imh]->Draw("same");
    }
  }
  TLegend *lchisq = new TLegend(0.6,0.45,0.8,0.65);
  lchisq->SetFillColor(0);
  lchisq->SetLineColor(0);
  lchisq->SetTextFont(132);
  lchisq->SetTextSize(0.07);
  for(int ifile=0;ifile<Nfiles;ifile++){
    //for(int imh=0;imh<1;imh++){
    for(int imh=0;imh<Nmh;imh++){
      lchisq->AddEntry(hchisq1D[ifile][imh],latexname[ifile][imh],"l");
    }
  }
  lchisq->Draw();


  // Delta Chisq 

  double Max_deltachisq = 0.;
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      if(hdeltachisq1D_GlobMin[ifile][imh]->GetMaximum()>Max_deltachisq) 
	Max_deltachisq = hdeltachisq1D_GlobMin[ifile][imh]->GetMaximum();
    }
  }

  sprintf(name,"cDeltaChisq1D");
  cDeltaChisq = new TCanvas(name,name);

  hInit->SetMaximum(Max_deltachisq*1.2);
  hInit->GetXaxis()->SetTitle(titleX);
  hInit->GetYaxis()->SetTitle(titleY);
  hInit->DrawCopy();

  // Draw reactors box
  const double sinsqth13_nom = 0.0219;
  const double sinsqth13_err = 0.0012;
  double X1 = sinsqth13_nom - sinsqth13_err;
  double X2 = sinsqth13_nom + sinsqth13_err;
  double Y1 = hInit->GetBinContent(hInit->GetMinimumBin());
  //double Y2 = hInit->GetMaximum();
  double Y2 = 45.;
  TBox *breactors = new TBox(X1,Y1,X2,Y2);
  breactors->SetFillColor(kYellow+1);
  if(draw_reactors) breactors->Draw("");


  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      hdeltachisq1D_GlobMin[ifile][imh]->SetMaximum(Max_deltachisq*1.2);
      if(same_mh_style) hdeltachisq1D_GlobMin[ifile][imh]->SetLineStyle(1);
      else hdeltachisq1D_GlobMin[ifile][imh]->SetLineStyle(imh+1);
      hdeltachisq1D_GlobMin[ifile][imh]->SetLineWidth(2);
      if(same_mh_color) hdeltachisq1D_GlobMin[ifile][imh]->SetLineColor(ifile+1);
      else hdeltachisq1D_GlobMin[ifile][imh]->SetLineColor(ifile+imh+1);
      if((ifile+imh)==0){ 
	hdeltachisq1D_GlobMin[ifile][imh]->DrawCopy("same");
	if(do_fc_conf_interv) hdeltachisq1D_GlobMin[ifile][imh]->DrawCopy("l same");
      }
      else{
	hdeltachisq1D_GlobMin[ifile][imh]->DrawCopy("same");
	if(do_fc_conf_interv) hdeltachisq1D_GlobMin[ifile][imh]->DrawCopy("l same");
      }
    }
  }
  
  TLegend *ldeltachisq = new TLegend(0.6,0.45,0.8,0.6);
  ldeltachisq->SetFillColor(0);
  ldeltachisq->SetLineColor(0);
  ldeltachisq->SetTextFont(132);
  ldeltachisq->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    //for(int imh=0;imh<1;imh++){
    for(int imh=0;imh<Nmh;imh++){
      ldeltachisq->AddEntry(hdeltachisq1D_GlobMin[ifile][imh],latexname[ifile][imh],"l");
    }
  }
  ldeltachisq->Draw();

  if(do_fc){
    if(draw_90CL){
      gFC90CL_nh   ->Draw("l same");
      gFC90CL_ih   ->Draw("l same");
      if(draw_fc_err){ 
	gFC90CL_nh_p1->Draw("l same");
	gFC90CL_nh_m1->Draw("l same");
	gFC90CL_ih_p1->Draw("l same");
	gFC90CL_ih_m1->Draw("l same");
      }

      if(do_fc_conf_interv){
	
	min90CL_nh->SetLineWidth(2);
	min90CL_nh->SetLineStyle(hdeltachisq1D_GlobMin[0][0]->GetLineStyle());
	min90CL_nh->SetLineColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	min90CL_nh->Draw(); 
	hatch_min90CL_nh->SetLineWidth(1);
	hatch_min90CL_nh->SetLineStyle(hdeltachisq1D_GlobMin[0][0]->GetLineStyle());
	hatch_min90CL_nh->SetFillStyle(3004); // 3004 3844
	hatch_min90CL_nh->SetLineColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	hatch_min90CL_nh->SetFillColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	hatch_min90CL_nh->Draw();
	
	max90CL_nh->SetLineWidth(2);
	max90CL_nh->SetLineStyle(hdeltachisq1D_GlobMin[0][0]->GetLineStyle());
	max90CL_nh->SetLineColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	max90CL_nh->Draw(); 
	hatch_max90CL_nh->SetLineWidth(1);
	hatch_max90CL_nh->SetLineStyle(hdeltachisq1D_GlobMin[0][0]->GetLineStyle());
	hatch_max90CL_nh->SetFillStyle(3004); // 3004 3844
	hatch_max90CL_nh->SetLineColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	hatch_max90CL_nh->SetFillColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	hatch_max90CL_nh->Draw();
	
	// min90CL_nh_2->SetLineWidth(2);
	// min90CL_nh_2->SetLineStyle(hdeltachisq1D_GlobMin[0][0]->GetLineStyle());
	// min90CL_nh_2->SetLineColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	// min90CL_nh_2->Draw(); 
	// hatch_min90CL_nh_2->SetLineWidth(1);
	// hatch_min90CL_nh_2->SetLineStyle(hdeltachisq1D_GlobMin[0][0]->GetLineStyle());
	// hatch_min90CL_nh_2->SetFillStyle(3004); // 3004 3844
	// hatch_min90CL_nh_2->SetLineColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	// hatch_min90CL_nh_2->SetFillColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	// hatch_min90CL_nh_2->Draw();
	
	// // max90CL_nh_2->SetLineWidth(2);
	// // max90CL_nh_2->SetLineStyle(hdeltachisq1D_GlobMin[0][0]->GetLineStyle());
	// // max90CL_nh_2->SetLineColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	// // max90CL_nh_2->Draw(); 
	// // hatch_max90CL_nh_2->SetLineWidth(1);
	// // hatch_max90CL_nh_2->SetLineStyle(hdeltachisq1D_GlobMin[0][0]->GetLineStyle());
	// // hatch_max90CL_nh_2->SetFillStyle(3004); // 3004 3844
	// // hatch_max90CL_nh_2->SetLineColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	// // hatch_max90CL_nh_2->SetFillColor(hdeltachisq1D_GlobMin[0][0]->GetLineColor());
	// // hatch_max90CL_nh_2->Draw();

	min90CL_ih->SetLineWidth(2);
	min90CL_ih->SetLineStyle(hdeltachisq1D_GlobMin[0][1]->GetLineStyle());
	min90CL_ih->SetLineColor(hdeltachisq1D_GlobMin[0][1]->GetLineColor());
	min90CL_ih->Draw(); 
	hatch_min90CL_ih->SetLineWidth(1);
	hatch_min90CL_ih->SetLineStyle(hdeltachisq1D_GlobMin[0][1]->GetLineStyle());
	hatch_min90CL_ih->SetFillStyle(3004); // 3005 3844
	hatch_min90CL_ih->SetLineColor(hdeltachisq1D_GlobMin[0][1]->GetLineColor());
	hatch_min90CL_ih->SetFillColor(hdeltachisq1D_GlobMin[0][1]->GetLineColor());
	hatch_min90CL_ih->Draw();
	
	max90CL_ih->SetLineWidth(2);
	max90CL_ih->SetLineStyle(hdeltachisq1D_GlobMin[0][1]->GetLineStyle());
	max90CL_ih->SetLineColor(hdeltachisq1D_GlobMin[0][1]->GetLineColor());
	max90CL_ih->Draw();      
	hatch_max90CL_ih->SetLineWidth(1);
	hatch_max90CL_ih->SetLineStyle(hdeltachisq1D_GlobMin[0][1]->GetLineStyle());
	hatch_max90CL_ih->SetFillStyle(3004); // 3005 3844
	hatch_max90CL_ih->SetLineColor(hdeltachisq1D_GlobMin[0][1]->GetLineColor());
	hatch_max90CL_ih->SetFillColor(hdeltachisq1D_GlobMin[0][1]->GetLineColor());
	hatch_max90CL_ih->Draw();
      }
    }
    if(draw_2sigmaCL){
      gFC2sigma_nh   ->Draw("l same");
      gFC2sigma_ih   ->Draw("l same");
      if(draw_fc_err){ 
	gFC2sigma_nh_p1->Draw("l same");
	gFC2sigma_nh_m1->Draw("l same");
	gFC2sigma_ih_p1->Draw("l same");
	gFC2sigma_ih_m1->Draw("l same");
      }

    }
  }

  TLegend *lcritval = new TLegend(0.6,0.45,0.8,0.75);
  lcritval->SetFillColor(0);
  lcritval->SetLineColor(0);
  lcritval->SetTextFont(132);
  lcritval->SetTextSize(0.05);
  if(do_fc){
    if(draw_90CL){
      //lcritval->AddEntry(gFC90CL_nh,"#Delta #chi^{2}_{crit} (90%CL) - NH","l");
      //lcritval->AddEntry(gFC90CL_ih,"#Delta #chi^{2}_{crit} (90%CL) - IH","l");
      lcritval->AddEntry(gFC90CL_nh,"-2lnL_{crit} (90% CL)","l");
    }
    if(draw_2sigmaCL){
      //lcritval->AddEntry(gFC2sigma_nh,"#Delta #chi^{2}_{crit} (2#sigma) - NH","l");
      //lcritval->AddEntry(gFC2sigma_ih,"#Delta #chi^{2}_{crit} (2#sigma) - IH","l");
      lcritval->AddEntry(gFC2sigma_nh,"-2lnL_{crit} (2#sigma CL)","l");
    }
    lcritval->Draw();
  }
  
  TLine *l1sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),1.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),1.);
  TLine *l90 = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),2.71,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),2.71);
  TLine *l2sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),4.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),4.);
  TLine *l3sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),9.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),9.);
  TLine *l4sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),16.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),16.);
  TLine *l5sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),25.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),25.);
  l1sigma->SetLineWidth(2);
  l90->SetLineWidth(2);
  l2sigma->SetLineWidth(2);
  l3sigma->SetLineWidth(2);
  l4sigma->SetLineWidth(2);
  l5sigma->SetLineWidth(2);
  l1sigma->SetLineColor(kBlue);
  l90->SetLineColor(kBlue);
  l2sigma->SetLineColor(kBlue);
  l3sigma->SetLineColor(kBlue);
  l4sigma->SetLineColor(kBlue);
  l5sigma->SetLineColor(kBlue);
  if(do_const){
    if(draw_1sigmaCL) l1sigma->Draw();
    if(draw_90CL)     l90->Draw();
    if(draw_2sigmaCL) l2sigma->Draw();
    if(draw_3sigmaCL) l3sigma->Draw();
    if(draw_4sigmaCL) l4sigma->Draw();
    if(draw_5sigmaCL) l5sigma->Draw();
  }
  
  TLatex *latLine1Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),1.1,"1#sigma");
  TLatex *latLine90CL = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),2.81,"90% CL");
  TLatex *latLine2Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),4.1,"2#sigma");
  TLatex *latLine3Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),9.1,"3#sigma");
  TLatex *latLine4Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),16.1,"4#sigma");
  TLatex *latLine5Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),25.1,"5#sigma");
  latLine1Sigma->SetTextColor(kBlue);
  latLine90CL->SetTextColor(kBlue);
  latLine2Sigma->SetTextColor(kBlue);
  latLine3Sigma->SetTextColor(kBlue);
  latLine4Sigma->SetTextColor(kBlue);
  latLine5Sigma->SetTextColor(kBlue);
  latLine1Sigma->SetTextSize(0.07);
  latLine90CL->SetTextSize(0.07);
  latLine2Sigma->SetTextSize(0.07);
  latLine3Sigma->SetTextSize(0.07);
  latLine4Sigma->SetTextSize(0.07);
  latLine5Sigma->SetTextSize(0.07);
  if(do_const){
    if(draw_1sigmaCL) latLine1Sigma->Draw("same");
    if(draw_90CL)     latLine90CL->Draw("same");
    if(draw_2sigmaCL) latLine2Sigma->Draw("same");
    if(draw_3sigmaCL) latLine3Sigma->Draw("same");
    if(draw_4sigmaCL) latLine4Sigma->Draw("same");
    if(draw_5sigmaCL) latLine5Sigma->Draw("same");
  }

  TLatex *latPreliminary = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),1.1,str_preliminary);
  latPreliminary->SetTextSize(0.05);
  if(draw_preliminary) latPreliminary->Draw("same");

  // Store in root file
  if(SaveRootFile){
    TFile* file = new TFile(Form("deltachisq_%s.root",tag.c_str()),"recreate");    
    for(int ifile=0;ifile<Nfiles;ifile++){
      for(int imh=0;imh<Nmh;imh++){
  	sprintf(name,"hdeltachisq_File_%i_MH_%i",ifile,imh);
  	hdeltachisq1D_GlobMin[ifile][imh]->Write();//(name);
      }
    }    
    file->Close();
  }
  


  //

  if(truepar0 != kBadResult && truepar1 != kBadResult){

    int bfbin = -999;
    int truebin0=-999;
    int truebin1=-999;
    string mh0 = "";
    string mh1 = "";
    int mhbin0 = -999;
    int mhbin1 = -999;
    
    if(truemh0=="normal"){        
      mhbin0 = 0;
      mh0 = "NH";
    }
    else if(truemh0=="inverted"){
      mhbin0 = 1;
      mh0 = "IH";
    }
    
    if(truemh1=="normal"){        
      mhbin1 = 0;
      mh1 = "NH";
    }
    else if(truemh1=="inverted"){
      mhbin1 = 1;
      mh1 = "IH";
    }
    
    truebin0 = hchisq1D[0][mhbin0]->FindBin(truepar0);
    truebin1 = hchisq1D[0][mhbin1]->FindBin(truepar1);
    
    double truechisq0 = hchisq1D[0][mhbin0]->GetBinContent(truebin0);
    double truechisq1 = hchisq1D[0][mhbin1]->GetBinContent(truebin1);
    
    cout << std::setprecision(10);
    cout << "chisq(" << truepar0 << "," << mh0 << ") - " 
	 << "chisq(" << truepar1 << "," << mh1 << ") = "
	 << truechisq0 - truechisq1 
	 << endl;
  }

























  ////////////////////////////////////////////////////
  //                                                //
  //             Brazilian flag plot                //
  //                                                //
  ////////////////////////////////////////////////////
  
  if(doBrazilian){
    
    if(do_fc==false){
      cout << endl;
      cout << "You cannot do brazilian flag plots w/o Feldman-Cousins critical values" << endl;
      cout << endl;
      return;
    }

    cout << endl;
    cout << "Make the Brazilian plot" << endl;
    cout << endl;
    
    if(Nfiles!=1){
      cerr << endl;
      cerr << "The number of files is " << Nfiles << endl;
      cerr << "You can make the Brazilian plot only with 1 file" << endl;
      cerr << endl;
    }
    if(Nmh!=2){
      cerr << endl;
      cerr << "The number of mass hierarchies is " << Nmh << endl;
      cerr << "You can make the Brazilian plot only with both mass hierarchies" << endl;
      cerr << endl;
    }

    TFile *filetoys = new TFile(infilename_toys,"READ");
    if(!filetoys->IsOpen()){
      cerr << endl;
      cerr << "File " << infilename << " is not open" << endl;
      cerr << endl;
    }
    cout << endl;
    cout << "Reading file " << infilename_toys << endl;
    cout << endl;
    TTree *treetoys = (TTree*)filetoys->Get("fitnt");
    if(!treetoys){
      cerr << endl;
      cerr << "Tree object is empty" << endl;
      cerr << endl;
    }

    TH1D *hdeltachisq1D_toy_nh; // take it from hchisq1D_toy_nh at each event
    TH1D *hdeltachisq1D_toy_ih; // take it from hchisq1D_toy_nh at each event

    TH1D *hNevts_1Rmu    = new TH1D("hNevts_1Rmu","hNevts_1Rmu",500,50,200);
    TH1D *hNevts_RHC1Rmu = new TH1D("hNevts_RHC1Rmu","hNevts_RHC1Rmu",240,20,100);
    TH1D *hNevts_1Re    = new TH1D("hNevts_1Re","hNevts_1Re",120,0,60);
    TH1D *hNevts_RHC1Re = new TH1D("hNevts_RHC1Re","hNevts_RHC1Re",40,0,20);

    TH1D *hchisq1D_toy_nh = new TH1D("hchisq1D_toy_nh","hchisq1D_toy_nh",1,0,1);
    TH1D *hchisq1D_toy_ih = new TH1D("hchisq1D_toy_ih","hchisq1D_toy_ih",1,0,1);
    TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
    TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
    TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
    TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);

    treetoys->SetBranchAddress("hchisq1D_nh",        &hchisq1D_toy_nh      );
    treetoys->SetBranchAddress("hchisq1D_ih",        &hchisq1D_toy_ih      );
    treetoys->SetBranchAddress("hErecoToy_1Rmu",     &hErecoToy_1Rmu     );
    treetoys->SetBranchAddress("hErecoToy_1Re",      &hErecoToy_1Re      );
    treetoys->SetBranchAddress("hErecoToy_RHC1Rmu",  &hErecoToy_RHC1Rmu  );
    treetoys->SetBranchAddress("hErecoToy_RHC1Re",   &hErecoToy_RHC1Re   );

    int Nentries = treetoys->GetEntries();
  
    treetoys->GetEntry(0);

    int NbinsChisq = hchisq1D_toy_nh->GetNbinsX();

    vector<double> *vec_deltachisq_bin_nh = new vector<double>[NbinsChisq];
    vector<double> *vec_deltachisq_bin_ih = new vector<double>[NbinsChisq];
    vector<double> *vec_like_penalty      = new vector<double>[NbinsChisq];

    // (delta)chisq histogram for each bin of dcp and mh --> calculate the +-1,2sigma ranges
    TH1D **hchisq1D_toy_nh_bin      = new TH1D[NbinsChisq];
    TH1D **hchisq1D_toy_ih_bin      = new TH1D[NbinsChisq];
    TH1D **hdeltachisq1D_toy_nh_bin = new TH1D[NbinsChisq];
    TH1D **hdeltachisq1D_toy_ih_bin = new TH1D[NbinsChisq];
    for(int ibin=0;ibin<NbinsChisq;ibin++){
      hchisq1D_toy_nh_bin[ibin]      = new TH1D(Form("hchisq1D_toy_nh_bin%d",ibin),     Form("hchisq1D_toy_nh_bin%d",ibin),500000,0,500);
      hchisq1D_toy_ih_bin[ibin]      = new TH1D(Form("hchisq1D_toy_ih_bin%d",ibin),     Form("hchisq1D_toy_ih_bin%d",ibin),500000,0,500);
      hdeltachisq1D_toy_nh_bin[ibin] = new TH1D(Form("hdeltachisq1D_toy_nh_bin%d",ibin),Form("hdeltachisq1D_toy_nh_bin%d",ibin),500000,0,30);
      hdeltachisq1D_toy_ih_bin[ibin] = new TH1D(Form("hdeltachisq1D_toy_ih_bin%d",ibin),Form("hdeltachisq1D_toy_ih_bin%d",ibin),500000,0,30);
    }

    // // store each (delta)chisq histogram of each toy --> just to plot
    // TH1D **hchisq1D_toy_nh_store      = new TH1D[Nentries];
    // TH1D **hchisq1D_toy_ih_store      = new TH1D[Nentries];
    // TH1D **hdeltachisq1D_toy_nh_store = new TH1D[Nentries];
    // TH1D **hdeltachisq1D_toy_ih_store = new TH1D[Nentries];

    // mean sensitivity
    TH1D *hchisq1D_nh_mean = hchisq1D_toy_nh->Clone("hchisq1D_nh_mean");
    hchisq1D_nh_mean->Reset();
    TH1D *hchisq1D_ih_mean = hchisq1D_toy_ih->Clone("hchisq1D_ih_mean");
    hchisq1D_ih_mean->Reset();

    
    int nbins = hchisq1D_ih_mean->GetXaxis()->GetNbins();
    double lowedge = hchisq1D_ih_mean->GetXaxis()->GetBinLowEdge(1);
    double upedge = hchisq1D_ih_mean->GetXaxis()->GetBinUpEdge(nbins);
    //cout << nbins << " " << lowedge << " " << upedge << endl;
    //TH2D *hchisq1D_Vs_deltachisq_nh = new TH2D("hchisq1D_Vs_deltachisq_nh","hchisq1D_Vs_deltachisq_nh",hchisq1D_ih_mean->GetNbinsX(),hchisq1D_ih_mean->GetXaxis()->GetXbins()->GetArray(),100,0,20);
    TH2D *hdcp_Vs_deltachisq_nh = new TH2D("hdcp_Vs_deltachisq_nh","hdcp_Vs_deltachisq_nh",nbins,lowedge,upedge,100,0,40);
    TH2D *hdcp_Vs_deltachisq_ih = new TH2D("hdcp_Vs_deltachisq_ih","hdcp_Vs_deltachisq_ih",nbins,lowedge,upedge,100,0,40);
    
    cout << "Enter the loop of " << Nentries << " entries" << endl;

    for(int ientry=0;ientry<Nentries;ientry++){
    
      treetoys->GetEntry(ientry);    
    
      if(ientry%1000==0) cout << "Entry " << ientry << endl;
    
      // hchisq1D_toy_nh_store[ientry]      = (TH1D*) hchisq1D_toy_nh->Clone(Form("hchisq1D_toy%d_nh_store",ientry));
      // hchisq1D_toy_ih_store[ientry]      = (TH1D*) hchisq1D_toy_ih->Clone(Form("hchisq1D_toy%d_ih_store",ientry));    
		    
      // Take global minimum
      double min_glob = 0.;
      double min_nh = hchisq1D_toy_nh->GetBinContent(hchisq1D_toy_nh->GetMinimumBin());
      double min_ih = hchisq1D_toy_ih->GetBinContent(hchisq1D_toy_ih->GetMinimumBin());
      if(min_nh < min_ih){
	min_glob = min_nh;
	//continue;
      }
      else{
	min_glob = min_ih;
	//continue;
      }

      // Build mean sensitivity
      hchisq1D_nh_mean->Add(hchisq1D_toy_nh,1.);
      hchisq1D_ih_mean->Add(hchisq1D_toy_ih,1.);
      
      //cout << "Toy " << ientry << " --> " << min_glob << endl;

      // Take global deltachisq and shift deltachisq of nh and ih
    
      hdeltachisq1D_toy_nh = (TH1D*) hchisq1D_toy_nh ->Clone("hchisq1D_toy_nh"); 
      hdeltachisq1D_toy_ih = (TH1D*) hchisq1D_toy_ih ->Clone("hchisq1D_toy_ih");

      for(int ibin=0;ibin<NbinsChisq;ibin++){
      
	double bincont = hdeltachisq1D_toy_nh->GetBinContent(ibin+1);
	hdeltachisq1D_toy_nh->SetBinContent(ibin+1,bincont-min_glob);
								       
	bincont = hdeltachisq1D_toy_ih->GetBinContent(ibin+1);
	hdeltachisq1D_toy_ih->SetBinContent(ibin+1,bincont-min_glob);
    
	// hdeltachisq1D_toy_nh_store[ientry] = (TH1D*) hdeltachisq1D_toy_nh->Clone(Form("hdeltachisq1D_toy%d_nh_store",ientry));
	// hdeltachisq1D_toy_ih_store[ientry] = (TH1D*) hdeltachisq1D_toy_ih->Clone(Form("hdeltachisq1D_toy%d_ih_store",ientry));


      }


      
      // Apply cut on deltachisq
      
      int bin_mpiover2 = hdeltachisq1D_toy_nh->FindBin(-TMath::Pi()/2.);
      double deltachisq_fluct_mpiover2_nh = hdeltachisq1D_toy_nh->GetBinContent(bin_mpiover2);
      double deltachisq_fluct_mpiover2_ih = hdeltachisq1D_toy_ih->GetBinContent(bin_mpiover2);
      
      int bin_ppiover2 = hdeltachisq1D_toy_nh->FindBin(+TMath::Pi()/2.);
      double deltachisq_fluct_ppiover2_nh = hdeltachisq1D_toy_nh->GetBinContent(bin_ppiover2);
      double deltachisq_fluct_ppiover2_ih = hdeltachisq1D_toy_ih->GetBinContent(bin_ppiover2);

      //if(deltachisq_fluct_mpiover2_nh < 4) continue; 
      //if(deltachisq_fluct_ppiover2_nh < 4) continue; 

      ///////




      // Fill deltaCP bin histograms

      for(int ibin=0;ibin<NbinsChisq;ibin++){
	double chisqval = hchisq1D_toy_nh->GetBinContent(ibin+1);
	hchisq1D_toy_nh_bin[ibin]->Fill(chisqval);
	double deltachisqval = hdeltachisq1D_toy_nh->GetBinContent(ibin+1);
	hdeltachisq1D_toy_nh_bin[ibin]->Fill(deltachisqval);
	vec_deltachisq_bin_nh[ibin].push_back(deltachisqval);
	
	hdcp_Vs_deltachisq_nh->Fill(hdeltachisq1D_toy_nh->GetBinCenter(ibin+1),deltachisqval);

	chisqval = hchisq1D_toy_ih->GetBinContent(ibin+1);
	hchisq1D_toy_ih_bin[ibin]->Fill(chisqval);
	deltachisqval = hdeltachisq1D_toy_ih->GetBinContent(ibin+1);
	hdeltachisq1D_toy_ih_bin[ibin]->Fill(deltachisqval);
	vec_deltachisq_bin_ih[ibin].push_back(deltachisqval);
	
	hdcp_Vs_deltachisq_ih->Fill(hdeltachisq1D_toy_ih->GetBinCenter(ibin+1),deltachisqval);
	
	vec_like_penalty[ibin].push_back(1.); // no posterior predictive
      }

      

      // Fill histograms of # of events
      
      double nevents_rhc1re = hErecoToy_RHC1Re->Integral();    
      double nevents_1re = hErecoToy_1Re->Integral();    
      double nevents_rhc1rmu = hErecoToy_RHC1Rmu->Integral();    
      double nevents_1rmu = hErecoToy_1Rmu->Integral();    

      // cout << "nevents_1rmu = " << nevents_1rmu << endl;
      // cout << "nevents_rhc1rmu = " << nevents_rhc1rmu << endl;
      // cout << "nevents_1re = " << nevents_1re << endl;
      // cout << "nevents_rhc1re = " << nevents_rhc1re << endl;
      // cout << endl;

      hNevts_1Rmu->Fill(nevents_1rmu);
      hNevts_1Re->Fill(nevents_1re);
      hNevts_RHC1Rmu->Fill(nevents_rhc1rmu);
      hNevts_RHC1Re->Fill(nevents_rhc1re);

    } // end loop over events

    // Scale mean chisq to # of toys and get mean deltachisq (global mean nh and ih)

    hchisq1D_nh_mean->Scale(1./Nentries);
    hchisq1D_ih_mean->Scale(1./Nentries);

    TH1D *hdeltachisq1D_nh_mean = (TH1D*) hchisq1D_nh_mean->Clone("hdeltachisq1D_nh_mean");
    TH1D *hdeltachisq1D_ih_mean = (TH1D*) hchisq1D_ih_mean->Clone("hdeltachisq1D_ih_mean");

    double glob_min_mean = -999;
    double min_mean_nh = hchisq1D_nh_mean->GetBinContent(hchisq1D_nh_mean->GetMinimumBin());
    double min_mean_ih = hchisq1D_ih_mean->GetBinContent(hchisq1D_ih_mean->GetMinimumBin());
    if(min_mean_nh < min_mean_ih) glob_min_mean = min_mean_nh;
    else                          glob_min_mean = min_mean_ih;

    for(int ibin=0;ibin<NbinsChisq;ibin++){
      double bincont = hdeltachisq1D_nh_mean->GetBinContent(ibin+1);
      hdeltachisq1D_nh_mean->SetBinContent(ibin+1,bincont-glob_min_mean);
      bincont = hdeltachisq1D_ih_mean->GetBinContent(ibin+1);
      hdeltachisq1D_ih_mean->SetBinContent(ibin+1,bincont-glob_min_mean);
    }
    
    // Build the +-1,2 sigma histograms

    TH1D *hchisq1D_nh_p1sigma = hchisq1D_toy_nh->Clone("hchisq1D_nh_p1sigma");
    hchisq1D_nh_p1sigma->Reset();
    TH1D *hchisq1D_ih_p1sigma = hchisq1D_toy_ih->Clone("hchisq1D_ih_p1sigma");
    hchisq1D_ih_p1sigma->Reset();

    TH1D *hchisq1D_nh_p2sigma = hchisq1D_toy_nh->Clone("hchisq1D_nh_p1sigma");
    hchisq1D_nh_p2sigma->Reset();
    TH1D *hchisq1D_ih_p2sigma = hchisq1D_toy_ih->Clone("hchisq1D_ih_p1sigma");
    hchisq1D_ih_p2sigma->Reset();

    TH1D *hdeltachisq1D_nh_p1sigma = hchisq1D_toy_nh->Clone("hdeltachisq1D_nh_p1sigma");
    hdeltachisq1D_nh_p1sigma->Reset();
    TH1D *hdeltachisq1D_ih_p1sigma = hchisq1D_toy_ih->Clone("hdeltachisq1D_ih_p1sigma");
    hdeltachisq1D_ih_p1sigma->Reset();

    TH1D *hdeltachisq1D_nh_p2sigma = hchisq1D_toy_nh->Clone("hdeltachisq1D_nh_p1sigma");
    hdeltachisq1D_nh_p2sigma->Reset();
    TH1D *hdeltachisq1D_ih_p2sigma = hchisq1D_toy_ih->Clone("hdeltachisq1D_ih_p1sigma");
    hdeltachisq1D_ih_p2sigma->Reset();
 
    // Build p-value histograms
    TH1D *hPvalueChisq_nh = hchisq1D_toy_nh->Clone("hPvalueChisq_nh");
    hPvalueChisq_nh->Reset();
    TH1D *hPvalueChisq_ih = hchisq1D_toy_ih->Clone("hPvalueChisq_ih");
    hPvalueChisq_ih->Reset();

    TH1D *hZscoreChisq_nh = hchisq1D_toy_nh->Clone("hZscoreChisq_nh");
    hZscoreChisq_nh->Reset();
    TH1D *hZscoreChisq_ih = hchisq1D_toy_ih->Clone("hZscoreChisq_ih");
    hZscoreChisq_ih->Reset();

    TH1D *hProbChisq_90CL_nh = hchisq1D_toy_nh->Clone("hProbChisq_90CL_nh");
    hProbChisq_90CL_nh->Reset();
    TH1D *hProbChisq_90CL_ih = hchisq1D_toy_nh->Clone("hProbChisq_90CL_ih");
    hProbChisq_90CL_ih->Reset();

    TH1D *hProbChisq_2sigma_nh = hchisq1D_toy_nh->Clone("hProbChisq_2sigma_nh");
    hProbChisq_2sigma_nh->Reset();
    TH1D *hProbChisq_2sigma_ih = hchisq1D_toy_nh->Clone("hProbChisq_2sigma_ih");
    hProbChisq_2sigma_ih->Reset();

    //

    for(int ibin=0;ibin<NbinsChisq;ibin++){

      // Get +1,2sigma critical values of data wrt sensitivity on 1 true point of dcp 
      double p1sigma = CalcDLcrit_posterior(vec_deltachisq_bin_nh[ibin],vec_like_penalty[ibin],CLfrac_1sigma*100.,0.);
      double p2sigma = CalcDLcrit_posterior(vec_deltachisq_bin_nh[ibin],vec_like_penalty[ibin],CLfrac_2sigma*100.,0.);
      hdeltachisq1D_nh_p1sigma->SetBinContent(ibin+1,p1sigma);
      hdeltachisq1D_nh_p2sigma->SetBinContent(ibin+1,p2sigma);  
    
      p1sigma = CalcDLcrit_posterior(vec_deltachisq_bin_ih[ibin],vec_like_penalty[ibin],CLfrac_1sigma*100.,0.);
      p2sigma = CalcDLcrit_posterior(vec_deltachisq_bin_ih[ibin],vec_like_penalty[ibin],CLfrac_2sigma*100.,0.);
      hdeltachisq1D_ih_p1sigma->SetBinContent(ibin+1,p1sigma);
      hdeltachisq1D_ih_p2sigma->SetBinContent(ibin+1,p2sigma);  

      // P-value of data wrt sensitivity on 1 true point of dcp 
      TH1D *h_nh = (TH1D*) hdeltachisq1D_GlobMin[0][0]; // 0 --> NH
      double chisq_data = h_nh->GetBinContent(ibin+1);
      double pvalue = GetPvalue_rightside(vec_deltachisq_bin_nh[ibin],chisq_data);
      hPvalueChisq_nh->SetBinContent(ibin+1,pvalue);
      TF1 fnh("fnh","ROOT::Math::gaussian_quantile(1.-(x/2.),1.) ",0.,1.);
      double zscore = fnh.Eval(pvalue); 
      hZscoreChisq_nh->SetBinContent(ibin+1,zscore);

      TH1D *h_ih = (TH1D*) hdeltachisq1D_GlobMin[0][1]; // 1 --> IH
      chisq_data = h_ih->GetBinContent(ibin+1);
      pvalue = GetPvalue_rightside(vec_deltachisq_bin_ih[ibin],chisq_data);
      hPvalueChisq_ih->SetBinContent(ibin+1,pvalue);
      TF1 fih("fih","ROOT::Math::gaussian_quantile(1.-(x/2.),1.) ",0.,1.);
      zscore = fih.Eval(pvalue); 
      hZscoreChisq_ih->SetBinContent(ibin+1,zscore);


      // Probability to exclude dcp at 90% or 2sigma CL from sensitivity
      double critval_pt = gFC90CL_nh->Eval(hProbChisq_90CL_nh->GetBinCenter(ibin+1));
      double prob = GetPvalue_rightside(vec_deltachisq_bin_nh[ibin],critval_pt);
      hProbChisq_90CL_nh->SetBinContent(ibin+1,prob);
      critval_pt = gFC2sigma_nh->Eval(hProbChisq_2sigma_nh->GetBinCenter(ibin+1));
      prob = GetPvalue_rightside(vec_deltachisq_bin_nh[ibin],critval_pt);
      hProbChisq_2sigma_nh->SetBinContent(ibin+1,prob);

      double critval_pt = gFC90CL_ih->Eval(hProbChisq_90CL_ih->GetBinCenter(ibin+1));
      double prob = GetPvalue_rightside(vec_deltachisq_bin_ih[ibin],critval_pt);
      hProbChisq_90CL_ih->SetBinContent(ibin+1,prob);
      critval_pt = gFC2sigma_ih->Eval(hProbChisq_2sigma_ih->GetBinCenter(ibin+1));
      prob = GetPvalue_rightside(vec_deltachisq_bin_ih[ibin],critval_pt);
      hProbChisq_2sigma_ih->SetBinContent(ibin+1,prob);
    }
    
    
    // Draw histograms

    // TCanvas *cchisqAll = new TCanvas("cchisqAll","cchisqAll");
    // hchisq1D_toy_nh_store[0]->DrawCopy();
    // for(int itoy=0;itoy<Nentries;itoy++){
    //   hchisq1D_toy_nh_store[itoy]->Draw("same");
    // }

    // TCanvas *cdeltachisqAll = new TCanvas("cdeltachisqAll","cdeltachisqAll");
    // hdeltachisq1D_toy_nh_store[0]->DrawCopy();
    // for(int itoy=0;itoy<Nentries;itoy++){
    //   hdeltachisq1D_toy_nh_store[itoy]->Draw("same");
    // }

  
    // Mean sensitivity

    TCanvas *cchisq1D_mean = new TCanvas("cchisq1D_mean","cchisq1D_mean");
    hchisq1D_nh_mean->SetLineColor(kBlack);
    hchisq1D_ih_mean->SetLineColor(kRed);
    hchisq1D_nh_mean->DrawCopy();
    hchisq1D_ih_mean->Draw("same");
  
    TLegend *lchisq1D_mean = new TLegend(0.6,0.1,0.8,0.65);
    lchisq1D_mean->SetFillColor(0);
    lchisq1D_mean->SetLineColor(0);
    lchisq1D_mean->SetTextFont(132);
    lchisq1D_mean->SetTextSize(0.06);
    lchisq1D_mean->AddEntry(hchisq1D_nh_mean,"Normal Hierarchy","l");
    lchisq1D_mean->AddEntry(hchisq1D_ih_mean,"Inverted Hierarchy","l");
    lchisq1D_mean->Draw();

    hchisq1D_nh_mean->SetLineColor(kRed);
    hchisq1D_ih_mean->SetLineColor(kRed);

    // Brazilian plot Normal Hierarchy
  
    TCanvas *cBrazilianPlot_nh = new TCanvas("cBrazilianPlot_nh","cBrazilianPlot_nh");
    hInit->GetYaxis()->SetRangeUser(0.,18.);
    hInit->DrawCopy();
    hdeltachisq1D_nh_p2sigma->SetLineColor(kCyan);
    hdeltachisq1D_nh_p2sigma->SetFillColor(kCyan);
    hdeltachisq1D_nh_p2sigma->Draw("same");
    hdeltachisq1D_nh_p1sigma->SetLineColor(kGreen);
    hdeltachisq1D_nh_p1sigma->SetFillColor(kGreen);
    hdeltachisq1D_nh_p1sigma->Draw("same");
    hdeltachisq1D_nh_mean->SetLineWidth(2);
    hdeltachisq1D_nh_mean->SetLineColor(kRed);
    hdeltachisq1D_nh_mean->SetLineStyle(1);
    hdeltachisq1D_nh_mean->Draw("same");
    //hdeltachisq1D_GlobMin[0][0]->SetLineStyle(1); // 0 --> Normal Hierarchy
    // hdeltachisq1D_GlobMin[0][0]->SetLineWidth(2);
    // hdeltachisq1D_GlobMin[0][0]->SetLineColor(1);
    hdeltachisq1D_GlobMin[0][0]->Draw("same");
    if(do_fc){
      if(draw_90CL)     gFC90CL_nh->Draw("l same");
      if(draw_2sigmaCL) gFC2sigma_nh->Draw("l same");
    }
    gPad->RedrawAxis();
    
    if(draw_preliminary) latPreliminary->Draw("same");

    TLegend *lBrazilianPlot_nh = new TLegend(0.6,0.1,0.8,0.65);
    lBrazilianPlot_nh->SetFillColor(0);
    lBrazilianPlot_nh->SetLineColor(0);
    lBrazilianPlot_nh->SetTextFont(132);
    lBrazilianPlot_nh->SetTextSize(0.06);
    //lBrazilianPlot_nh->AddEntry(hdeltachisq1D_nh_p1sigma,"Exp. #Delta #chi^{2} 1#sigma spread","f");
    //lBrazilianPlot_nh->AddEntry(hdeltachisq1D_nh_p2sigma,"Exp. #Delta #chi^{2} 2#sigma spread","f");
    lBrazilianPlot_nh->AddEntry(hdeltachisq1D_nh_p1sigma,"68.27\% of toys MC","f");
    lBrazilianPlot_nh->AddEntry(hdeltachisq1D_nh_p2sigma,"95.45\% of toys MC","f");
    lBrazilianPlot_nh->AddEntry(hdeltachisq1D_nh_mean,"Mean Exp. -2lnL","l");
    //lBrazilianPlot_nh->AddEntry(hdeltachisq1D_nh_mean,"Mean Exp. #Delta #chi^{2}","l");
    if(draw_90CL){
      lBrazilianPlot_nh->AddEntry(gFC90CL_nh,"-2lnL_{crit} (90% CL)","l");
      //lBrazilianPlot_nh->AddEntry(gFC90CL_nh,"#Delta #chi^{2}_{crit} (90% CL)","l");
    }
    if(draw_2sigmaCL){
      lBrazilianPlot_nh->AddEntry(gFC2sigma_nh,"-2lnL_{crit} (2#sigma CL)","l");
    }
    lBrazilianPlot_nh->AddEntry(hdeltachisq1D_GlobMin[0][0],"Data","l");
    lBrazilianPlot_nh->Draw();

    TLatex *latNormal = new TLatex(hdeltachisq1D_nh_mean->GetXaxis()->GetBinLowEdge(5),1.1,"Normal Hierarchy");
    latNormal->SetTextSize(0.06);
    latNormal->Draw("same");


    // Brazilian plot Inverted Hierarchy
  
    TCanvas *cBrazilianPlot_ih = new TCanvas("cBrazilianPlot_ih","cBrazilianPlot_ih");
    hInit->GetYaxis()->SetRangeUser(0.,18.);
    hInit->DrawCopy();
    hdeltachisq1D_ih_p2sigma->SetLineColor(kCyan);
    hdeltachisq1D_ih_p2sigma->SetFillColor(kCyan);
    hdeltachisq1D_ih_p2sigma->Draw("same");
    hdeltachisq1D_ih_p1sigma->SetLineColor(kGreen);
    hdeltachisq1D_ih_p1sigma->SetFillColor(kGreen);
    hdeltachisq1D_ih_p1sigma->Draw("same");
    hdeltachisq1D_ih_mean->SetLineWidth(2);
    hdeltachisq1D_ih_mean->SetLineColor(kRed);
    hdeltachisq1D_ih_mean->SetLineStyle(1);
    hdeltachisq1D_ih_mean->Draw("same");
    //hdeltachisq1D_GlobMin[0][1]->SetLineStyle(1); // 1 --> Inverted Hierarchy
    //hdeltachisq1D_GlobMin[0][1]->SetLineWidth(2);
    //hdeltachisq1D_GlobMin[0][1]->SetLineColor(1);
    hdeltachisq1D_GlobMin[0][1]->Draw("same");
    if(do_fc){
      if(draw_90CL)     gFC90CL_ih->Draw("l same");
      if(draw_2sigmaCL) gFC2sigma_ih->Draw("l same");
    }
    gPad->RedrawAxis();

    if(draw_preliminary) latPreliminary->Draw("same");

    TLegend *lBrazilianPlot_ih = new TLegend(0.6,0.1,0.8,0.65);
    lBrazilianPlot_ih->SetFillColor(0);
    lBrazilianPlot_ih->SetLineColor(0);
    lBrazilianPlot_ih->SetTextFont(132);
    lBrazilianPlot_ih->SetTextSize(0.06);
    //lBrazilianPlot_ih->AddEntry(hdeltachisq1D_ih_p1sigma,"Exp. #Delta #chi^{2} 1#sigma spread","f");
    //lBrazilianPlot_ih->AddEntry(hdeltachisq1D_ih_p2sigma,"Exp. #Delta #chi^{2} 2#sigma spread","f");
    lBrazilianPlot_ih->AddEntry(hdeltachisq1D_ih_p1sigma,"68.27\% of toys MC","f");
    lBrazilianPlot_ih->AddEntry(hdeltachisq1D_ih_p2sigma,"95.45\% of toys MC","f");
    //lBrazilianPlot_ih->AddEntry(hdeltachisq1D_ih_mean,"Mean Exp. #Delta #chi^{2}","l");
    lBrazilianPlot_ih->AddEntry(hdeltachisq1D_ih_mean,"Mean Exp. -2lnL","l");
    if(draw_90CL)
      lBrazilianPlot_ih->AddEntry(gFC90CL_ih,"-2lnL_{crit} (90% CL)","l");
    if(draw_2sigmaCL)
      lBrazilianPlot_ih->AddEntry(gFC2sigma_ih,"-2lnL_{crit} (2#sigma CL)","l");
    lBrazilianPlot_ih->AddEntry(hdeltachisq1D_GlobMin[0][1],"Data","l");
    lBrazilianPlot_ih->Draw();

    TLatex *latInverted = new TLatex(hdeltachisq1D_ih_mean->GetXaxis()->GetBinLowEdge(5),1.1,"Inverted Hierarchy");
    latInverted->SetTextSize(0.06);
    latInverted->Draw("same");

    //
  
    TCanvas *cPvalueChisq = new TCanvas("cPvalueChisq","cPvalueChisq");  
    cPvalueChisq->SetLogy();
    cPvalueChisq->SetGridx();
    cPvalueChisq->SetGridy();
    hInit->GetYaxis()->SetRangeUser(0.01,1.1);
    hInit->GetYaxis()->SetTitle("P-value");
    hInit->DrawCopy();
    hPvalueChisq_nh->SetLineWidth(2);
    hPvalueChisq_nh->SetLineColor(kBlack);
    hPvalueChisq_nh->Draw("same");
    hPvalueChisq_ih->SetLineWidth(2);
    hPvalueChisq_ih->SetLineColor(kRed);
    hPvalueChisq_ih->Draw("same");
    if(draw_preliminary) latPreliminary->Draw("same");

    TLegend *lPvalueChisq = new TLegend(0.6,0.4,0.8,0.65);
    lPvalueChisq->SetFillColor(0);
    lPvalueChisq->SetLineColor(0);
    lPvalueChisq->SetTextFont(132);
    lPvalueChisq->SetTextSize(0.06);
    lPvalueChisq->AddEntry(hPvalueChisq_nh,"Normal Hierarchy","l");
    lPvalueChisq->AddEntry(hPvalueChisq_ih,"Inverted Hierarchy","l");
    lPvalueChisq->Draw();

    //

    TCanvas *cZscoreChisq = new TCanvas("cZscoreChisq","cZscoreChisq");  
    cZscoreChisq->SetGridx();
    cZscoreChisq->SetGridy();
    hInit->GetYaxis()->SetRangeUser(0,2.5);
    hInit->GetYaxis()->SetTitle("# of Std. Deviations");
    hInit->DrawCopy();
    hZscoreChisq_nh->SetLineWidth(2);
    hZscoreChisq_nh->SetLineColor(kBlack);
    hZscoreChisq_nh->Draw("same");
    hZscoreChisq_ih->SetLineWidth(2);
    hZscoreChisq_ih->SetLineColor(kRed);
    hZscoreChisq_ih->Draw("same");
    if(draw_preliminary) latPreliminary->Draw("same");
  
    TLegend *lZscoreChisq = new TLegend(0.6,0.4,0.8,0.65);
    lZscoreChisq->SetFillColor(0);
    lZscoreChisq->SetLineColor(0);
    lZscoreChisq->SetTextFont(132);
    lZscoreChisq->SetTextSize(0.06);
    lZscoreChisq->AddEntry(hZscoreChisq_nh,"Normal Hierarchy","l");
    lZscoreChisq->AddEntry(hZscoreChisq_ih,"Inverted Hierarchy","l");
    lZscoreChisq->Draw();
    if(draw_preliminary) latPreliminary->Draw("same");
  
    //

    TCanvas *cProbChisq = new TCanvas("cProbChisq","cProbChisq");
    cProbChisq->SetLogy();
    cProbChisq->SetGridx();
    cProbChisq->SetGridy();
    hInit->GetYaxis()->SetRangeUser(0.01,1.1);
    hInit->GetYaxis()->SetTitle("Prob. to exclude #delta_{CP}");
    hInit->DrawCopy();
    hProbChisq_90CL_nh->SetLineWidth(2);
    hProbChisq_90CL_nh->SetLineColor(kBlack);
    hProbChisq_90CL_nh->SetLineStyle(1);
    hProbChisq_90CL_nh->Draw("same");
    hProbChisq_2sigma_nh->SetLineWidth(2);
    hProbChisq_2sigma_nh->SetLineColor(kRed);
    hProbChisq_2sigma_nh->SetLineStyle(1);
    hProbChisq_2sigma_nh->Draw("same");
    hProbChisq_90CL_ih->SetLineWidth(2);
    hProbChisq_90CL_ih->SetLineColor(kBlack);
    hProbChisq_90CL_ih->SetLineStyle(2);
    hProbChisq_90CL_ih->Draw("same");
    hProbChisq_2sigma_ih->SetLineWidth(2);
    hProbChisq_2sigma_ih->SetLineColor(kRed);
    hProbChisq_2sigma_ih->SetLineStyle(2);
    hProbChisq_2sigma_ih->Draw("same");
    if(draw_preliminary) latPreliminary->Draw("same");
  
    TLegend *lProbChisq = new TLegend(0.6,0.4,0.8,0.65);
    lProbChisq->SetFillColor(0);
    lProbChisq->SetLineColor(0);
    lProbChisq->SetTextFont(132);
    lProbChisq->SetTextSize(0.06);
    lProbChisq->AddEntry(hProbChisq_90CL_nh,"90% CL - Normal","l");
    lProbChisq->AddEntry(hProbChisq_90CL_ih,"90% CL - Inverted","l");
    lProbChisq->AddEntry(hProbChisq_2sigma_nh,"2#sigma CL - Normal","l");
    lProbChisq->AddEntry(hProbChisq_2sigma_ih,"2#sigma CL - Inverted","l");
    lProbChisq->Draw();

    //

    TCanvas *cdcp_Vs_deltachisq_nh = new TCanvas("cdcp_Vs_deltachisq_nh","cdcp_Vs_deltachisq_nh");
    hdcp_Vs_deltachisq_nh->GetXaxis()->SetTitle("#delta_{CP}");
    hdcp_Vs_deltachisq_nh->GetYaxis()->SetTitle("-2lnL");
    //hdcp_Vs_deltachisq_nh->GetYaxis()->SetTitle("#Delta #chi^{2}");
    hdcp_Vs_deltachisq_nh->Draw("colz");
    TLatex *latNormal_dcp_dchisq = new TLatex(hdeltachisq1D_nh_mean->GetXaxis()->GetBinLowEdge(5),1.1,"Normal Hierarchy");
    latNormal_dcp_dchisq->SetTextSize(0.06);
    latNormal_dcp_dchisq->Draw("same");
    
    TCanvas *cdcp_Vs_deltachisq_ih = new TCanvas("cdcp_Vs_deltachisq_ih","cdcp_Vs_deltachisq_ih");
    hdcp_Vs_deltachisq_ih->GetXaxis()->SetTitle("#delta_{CP}");
    hdcp_Vs_deltachisq_ih->GetYaxis()->SetTitle("-2lnL");
    //hdcp_Vs_deltachisq_ih->GetYaxis()->SetTitle("#Delta #chi^{2}");
    hdcp_Vs_deltachisq_ih->Draw("colz");
    TLatex *latInverted_dcp_dchisq = new TLatex(hdeltachisq1D_ih_mean->GetXaxis()->GetBinLowEdge(5),1.1,"Inverted Hierarchy");
    latInverted_dcp_dchisq->SetTextSize(0.06);
    latInverted_dcp_dchisq->Draw("same");
    
    //

    double dcp_test = 0.;
    int ibin = hdeltachisq1D_toy_nh_bin[0]->FindBin(0.);
  
    TCanvas *hdeltachisq1D_nh_bin = new TCanvas("hdeltachisq1D_nh_bin","hdeltachisq1D_nh_bin");
    hdeltachisq1D_toy_nh_bin[ibin]->DrawCopy();

    TCanvas *hchisq1D_nh_bin = new TCanvas("hchisq1D_nh_bin","hchisq1D_nh_bin");
    hchisq1D_toy_nh_bin[ibin]->DrawCopy();
    
    // 

    TCanvas *cNevts = new TCanvas("cNevts","cNevts");
    cNevts->Divide(2,2);
    cNevts->cd(1); hNevts_1Rmu->Draw();
    cNevts->cd(2); hNevts_1Re->Draw();
    cNevts->cd(3); hNevts_RHC1Rmu->Draw();
    cNevts->cd(4); hNevts_RHC1Re->Draw();


    // Print

    int bin_dcp0 = hProbChisq_2sigma_ih->FindBin(0.);
    double prob_excl_dcp0_90CL_nh = hProbChisq_90CL_nh->GetBinContent(bin_dcp0); 
    double prob_excl_dcp0_90CL_ih = hProbChisq_90CL_ih->GetBinContent(bin_dcp0); 
    double prob_excl_dcp0_2sigma_nh = hProbChisq_2sigma_nh->GetBinContent(bin_dcp0); 
    double prob_excl_dcp0_2sigma_ih = hProbChisq_2sigma_ih->GetBinContent(bin_dcp0); 
    int bin_dcpPi = hProbChisq_2sigma_ih->FindBin(TMath::Pi());
    double prob_excl_dcpPi_90CL_nh = hProbChisq_90CL_nh->GetBinContent(bin_dcpPi); 
    double prob_excl_dcpPi_90CL_ih = hProbChisq_90CL_ih->GetBinContent(bin_dcpPi); 
    double prob_excl_dcpPi_2sigma_nh = hProbChisq_2sigma_nh->GetBinContent(bin_dcpPi); 
    double prob_excl_dcpPi_2sigma_ih = hProbChisq_2sigma_ih->GetBinContent(bin_dcpPi); 

    cout << endl;
    cout << "Normal Hierarchy:" << endl;
    cout << " - Prob. to exclude dcp=0 at 90%CL: " << prob_excl_dcp0_90CL_nh << endl;
    cout << " - Prob. to exclude dcp=0 at 2sigma: " << prob_excl_dcp0_2sigma_nh << endl;
    cout << " - Prob. to exclude dcp=pi at 90%CL:" << prob_excl_dcpPi_90CL_nh << endl;
    cout << " - Prob. to exclude dcp=pi at 2sigma:" << prob_excl_dcpPi_2sigma_nh << endl;
    cout << "Inverted Hierarchy:" << endl;
    cout << " - Prob. to exclude dcp=0 at 90%CL: " << prob_excl_dcp0_90CL_ih << endl;
    cout << " - Prob. to exclude dcp=0 at 2sigma: " << prob_excl_dcp0_2sigma_ih << endl;
    cout << " - Prob. to exclude dcp=pi at 90CL:" << prob_excl_dcpPi_90CL_ih << endl;
    cout << " - Prob. to exclude dcp=pi at 2sigma:" << prob_excl_dcpPi_2sigma_ih << endl;
    cout << endl;
    cout << "# of events 1Rmu:" << endl;
    cout << " - Mean = " << hNevts_1Rmu->GetMean() << endl;
    cout << " - RMS = " << hNevts_1Rmu->GetRMS() << endl;
    cout << "# of events RHC1Rmu:" << endl;
    cout << " - Mean = " << hNevts_RHC1Rmu->GetMean() << endl;
    cout << " - RMS = " << hNevts_RHC1Rmu->GetRMS() << endl;
    cout << "# of events 1Re:" << endl;
    cout << " - Mean = " << hNevts_1Re->GetMean() << endl;
    cout << " - RMS = " << hNevts_1Re->GetRMS() << endl;
    cout << "# of events RHC1Re:" << endl;
    cout << " - Mean = " << hNevts_RHC1Re->GetMean() << endl;
    cout << " - RMS = " << hNevts_RHC1Re->GetRMS() << endl;
    cout << endl;
  }
  
  return;
}

double CalcBestFitErrorsMINOS(TH1D *histo, double CritVal, double val_bf){ 
  
  cout << endl;
  cout << "Compute MINOS errors" << endl;
  cout << endl;

  //double error = conv::kBadResult;
  const double CritChisq_1D = CritVal;
  
  // Find 68.27% CL contour (1sigma)
  
  double errneg = conv::kBadResult;
  double errpos = conv::kBadResult;

  int bin_bf = conv::kBadResult;        // best-fit bin
  double chisq_min = conv::kBadResult;  // best-fit chisq

  if(val_bf==conv::kBadResult){
    chisq_min = histo->GetBinContent(histo->GetMinimumBin());
    bin_bf = histo->GetMinimumBin();
    val_bf = histo->GetBinCenter(bin_bf);

    //cout << endl;
    cout << "Find the global minimum:" << endl;
    cout << "best-fit value = " << val_bf << endl;
    cout << "best-fit chisq = " << chisq_min << endl;
    //cout << endl;
  }
  else{    
    bin_bf = histo->FindBin(val_bf);             // best-fit bin
    chisq_min = histo->GetBinContent(bin_bf); // best-fit chisq

    //cout << endl;
    cout << "Find the global minimum:" << endl;
    cout << "best-fit value = " << val_bf << endl;
    cout << "best-fit chisq = " << chisq_min << endl;
    //cout << endl;
  }


  // cout << endl;
  // cout << "val_bf = " << val_bf << endl;
  // cout << "bin_bf = " << bin_bf << endl;
  // cout << "chisq_min = " << chisq_min << endl;
  // cout << endl;

  // Bins before (b) and after (a) the intersection with critical value
  int bin_negsigma_a = conv::kBadResult; 
  int bin_negsigma_b = conv::kBadResult;
  int bin_possigma_a = conv::kBadResult;
  int bin_possigma_b = conv::kBadResult;
  
  // Find param value at DeltaChisq=1 (negative error)
  
  for(int ibin=bin_bf;ibin>=0;ibin--){ // histogram convention     
    
    double bincont = histo->GetBinContent(ibin);
    // if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Bin content is negative!" << endl;
    //   cerr << "bincont = " << bincont << endl;
    //   cerr << endl;
    // }
    
    // cout << "ibin = " << ibin 
    // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
    // 	   << " - bincont = " << bincont << endl;
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1 (following x direction)
      bin_negsigma_b = ibin; 
      bin_negsigma_a = bin_negsigma_b+1; 
      break;
    }
    // else if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Negative DeltaChisq 1D:" << endl;
    //   cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
    //   cerr << endl;
    //   //exit(1);
    //   return conv::kBadResult;
    // }
        
  }
  
   
  // Positive error
  
  // cout << endl;
  // cout << "Positive error" << endl;
  // cout << endl;
  
  for(int ibin=bin_bf;ibin<=histo->GetNbinsX();ibin++){ // histogram convention 
    double bincont = histo->GetBinContent(ibin);
    // if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Negative DeltaChisq 1D:" << endl;
    //   cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
    //   cerr << endl;
    //   //exit(1);
    //   return conv::kBadResult;
    // }
    
    // cout << "ibin = " << ibin
    // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
    //  	   << " - bincont = " << bincont << endl;
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1
      bin_possigma_a = ibin; 
      bin_possigma_b = bin_possigma_a-1; 
      break;
    }
  }
  //cout << "bin_possigma_b = " << bin_possigma_b << endl;
  //cout << "bin_possigma_a = " << bin_possigma_a << endl;
  
  
  // If no bin is found it means the 1sigma bin is outside the physical region
  // and bin variables are still set to default value.
  // If so set left extreme to 1 (first bin of the histogram) and right extreme
  // to Nbins (last bin of the histogram)
  
  if(bin_negsigma_a == conv::kBadResult){
    bin_negsigma_a = 1;
    bin_negsigma_b = 1;
  }
  if(bin_possigma_b == conv::kBadResult){
    bin_possigma_a = histo->GetNbinsX();
    bin_possigma_b = histo->GetNbinsX();
  }
  
  // Check if the interval bins are correctly found
  if(bin_negsigma_a > bin_possigma_b){
    cerr << endl;
    cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
    cerr << "Left extreme bin is higher than right extreme bin" << endl;
    cerr << "The 1sigma error interval is messed up:" << endl;
    cerr << "bin_possigma_b = " << bin_possigma_b << endl;
    cerr << "bin_negsigma_a = " << bin_negsigma_a << endl;
    cerr << endl;      
    //exit(1);
    return conv::kBadResult;
  }
  
  double X1 = conv::kBadResult;
  double Y1 = conv::kBadResult;
  double X2 = conv::kBadResult;
  double Y2 = conv::kBadResult;
  double Xout = conv::kBadResult;
  
  // Extrapolation (left extreme). Do it only if 1sigma bin is in the physical region
  if(bin_negsigma_a>1){
    X1 = histo->GetXaxis()->GetBinCenter(bin_negsigma_b);
    Y1 = histo->GetBinContent(bin_negsigma_b);
    X2 = histo->GetXaxis()->GetBinCenter(bin_negsigma_a);
    Y2 = histo->GetBinContent(bin_negsigma_a);
    Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
  }
  else{
    Xout = histo->GetBinCenter(bin_negsigma_a);
  }

  errneg = fabs(val_bf - Xout); 
  //cout << "fparam_bf_errneg[0] = " << fparam_bf_errneg[0] << endl;
  
  // Extrapolation (right extreme). Do it only if 1sigma bin is in the physical region
  if(bin_possigma_b < histo->GetNbinsX()){
    X1 = histo->GetXaxis()->GetBinCenter(bin_possigma_b);
    Y1 = histo->GetBinContent(bin_possigma_b);
    X2 = histo->GetXaxis()->GetBinCenter(bin_possigma_a);
    Y2 = histo->GetBinContent(bin_possigma_a);
    Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
  }
  else{
    Xout = histo->GetBinCenter(bin_possigma_b);
  }
  errpos[0] = fabs(val_bf - Xout); 
  //cout << "fparam_bf_errpos[0] = " << fparam_bf_errpos[0] << endl;
  
  cout << endl;
  cout << "Best-fit = " << val_bf << " + " << errpos << " - " << errneg << endl;
  cout << endl;
  cout << "1 sigma range:"<< endl;
  cout << "[ " << val_bf - errneg << " ; " << val_bf + errpos << " ] " << endl;
  cout << endl;
}
 

double CalcBestFitErrorsMINOS_FC(TH1D *histo, double CritVal, double val_bf,double &lim_x_min, double &lim_x_max){ 
  
  cout << endl;
  cout << "Feldman-Cousins confidence intervals" << endl;
  cout << endl;

  //double error = conv::kBadResult;
  const double CritChisq_1D = CritVal;
  
  // Find 68.27% CL contour (1sigma)
  
  double errneg = conv::kBadResult;
  double errpos = conv::kBadResult;

  int bin_bf = conv::kBadResult;        // best-fit bin
  double chisq_min = conv::kBadResult;  // best-fit chisq

  if(val_bf==conv::kBadResult){
    chisq_min = histo->GetBinContent(histo->GetMinimumBin());
    bin_bf = histo->GetMinimumBin();
    val_bf = histo->GetBinCenter(bin_bf);

    //cout << endl;
    cout << "Find the global minimum:" << endl;
    cout << "best-fit value = " << val_bf << endl;
    cout << "best-fit chisq = " << chisq_min << endl;
    //cout << endl;
  }
  else{    
    bin_bf = histo->FindBin(val_bf);             // best-fit bin
    chisq_min = histo->GetBinContent(bin_bf); // best-fit chisq

    //cout << endl;
    cout << "Find the global minimum:" << endl;
    cout << "best-fit value = " << val_bf << endl;
    cout << "best-fit chisq = " << chisq_min << endl;
    //cout << endl;
  }


  // cout << endl;
  // cout << "val_bf = " << val_bf << endl;
  // cout << "bin_bf = " << bin_bf << endl;
  // cout << "chisq_min = " << chisq_min << endl;
  // cout << endl;

  // Bins before (b) and after (a) the intersection with critical value
  int bin_negsigma_a = conv::kBadResult; 
  int bin_negsigma_b = conv::kBadResult;
  int bin_possigma_a = conv::kBadResult;
  int bin_possigma_b = conv::kBadResult;
  
  // Find param value at DeltaChisq=1 (negative error)
  
  for(int ibin=bin_bf;ibin>=0;ibin--){ // histogram convention     
    
    double bincont = histo->GetBinContent(ibin);
    // if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Bin content is negative!" << endl;
    //   cerr << "bincont = " << bincont << endl;
    //   cerr << endl;
    // }
    
    //
    // Useful for checks
    //
    // cout << "ibin = " << ibin 
    // 	   << " - X = " << histo->GetXaxis()->GetBinCenter(ibin)
    // 	   << " - bincont = " << bincont << endl;
    
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1 (following x direction)
      bin_negsigma_b = ibin; 
      bin_negsigma_a = bin_negsigma_b+1; 
      break;
    }
    // else if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Negative DeltaChisq 1D:" << endl;
    //   cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
    //   cerr << endl;
    //   //exit(1);
    //   return conv::kBadResult;
    // }
        
  }
  
  
  //
  // Useful for checks
  //
  // double x1 = histo->GetBinCenter(bin_negsigma_b);
  // double y1 = histo->GetBinContent(bin_negsigma_b);
  // double x2 = histo->GetBinCenter(bin_negsigma_a);
  // double y2 = histo->GetBinContent(bin_negsigma_a);
  // cout << "x1 = " << x1 << ", y1 = " << y1 << endl;
  // cout << "x2 = " << x2 << ", y2 = " << y2 << endl;
  // cout << " - interpolation = " << GetInterpolationX(x1,y1,x2,y2,0.) << endl;



   
  // Positive error
  
  // cout << endl;
  // cout << "Positive error" << endl;
  // cout << endl;
  
  for(int ibin=bin_bf;ibin<=histo->GetNbinsX();ibin++){ // histogram convention 
    double bincont = histo->GetBinContent(ibin);
    // if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Negative DeltaChisq 1D:" << endl;
    //   cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
    //   cerr << endl;
    //   //exit(1);
    //   return conv::kBadResult;
    // }
    
    // cout << "ibin = " << ibin
    // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
    //  	   << " - bincont = " << bincont << endl;
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1
      bin_possigma_a = ibin; 
      bin_possigma_b = bin_possigma_a-1; 
      break;
    }
  }
  //cout << "bin_possigma_b = " << bin_possigma_b << endl;
  //cout << "bin_possigma_a = " << bin_possigma_a << endl;
  
  
  // If no bin is found it means the 1sigma bin is outside the physical region
  // and bin variables are still set to default value.
  // If so set left extreme to 1 (first bin of the histogram) and right extreme
  // to Nbins (last bin of the histogram)
  
  if(bin_negsigma_a == conv::kBadResult){
    bin_negsigma_a = 1;
    bin_negsigma_b = 1;
  }
  if(bin_possigma_b == conv::kBadResult){
    bin_possigma_a = histo->GetNbinsX();
    bin_possigma_b = histo->GetNbinsX();
  }
  
  // Check if the interval bins are correctly found
  if(bin_negsigma_a > bin_possigma_b){
    cerr << endl;
    cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
    cerr << "Left extreme bin is higher than right extreme bin" << endl;
    cerr << "The 1sigma error interval is messed up:" << endl;
    cerr << "bin_possigma_b = " << bin_possigma_b << endl;
    cerr << "bin_negsigma_a = " << bin_negsigma_a << endl;
    cerr << endl;      
    //exit(1);
    return conv::kBadResult;
  }
  
  double X1 = conv::kBadResult;
  double Y1 = conv::kBadResult;
  double X2 = conv::kBadResult;
  double Y2 = conv::kBadResult;
  double Xout = conv::kBadResult;
  
  // Extrapolation (left extreme). Do it only if 1sigma bin is in the physical region
  if(bin_negsigma_a>1){
    X1 = histo->GetXaxis()->GetBinCenter(bin_negsigma_b);
    Y1 = histo->GetBinContent(bin_negsigma_b);
    X2 = histo->GetXaxis()->GetBinCenter(bin_negsigma_a);
    Y2 = histo->GetBinContent(bin_negsigma_a);
    Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
  }
  else{
    Xout = histo->GetBinCenter(bin_negsigma_a);
  }

  errneg = fabs(val_bf - Xout); 
  //cout << "fparam_bf_errneg[0] = " << fparam_bf_errneg[0] << endl;
  
  // Extrapolation (right extreme). Do it only if 1sigma bin is in the physical region
  if(bin_possigma_b < histo->GetNbinsX()){
    X1 = histo->GetXaxis()->GetBinCenter(bin_possigma_b);
    Y1 = histo->GetBinContent(bin_possigma_b);
    X2 = histo->GetXaxis()->GetBinCenter(bin_possigma_a);
    Y2 = histo->GetBinContent(bin_possigma_a);
    Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
  }
  else{
    Xout = histo->GetBinCenter(bin_possigma_b);
  }
  errpos[0] = fabs(val_bf - Xout); 
  //cout << "fparam_bf_errpos[0] = " << fparam_bf_errpos[0] << endl;
  
  cout << endl;
  cout << "Best-fit = " << val_bf << " + " << errpos << " - " << errneg << endl;
  cout << endl;
  cout << "The confidence interval is:"<< endl;
  cout << "[ " << val_bf - errneg << " ; " << val_bf + errpos << " ] " << endl;
  cout << endl;

  lim_x_min = val_bf - errneg;
  lim_x_max = val_bf + errpos;
  
}
 


double GetInterpolationX(double X1,double Y1,double X2, double Y2,double Yin){
  
  double slope = (Y2-Y1)/(X2-X1);
  double Xout = (Yin-Y1)/slope + X1;
  return Xout;

  // cout << "Yin = " << Yin << endl;
  // cout << "X1 = " << X1 << " - X2 = " << X2 << endl;
  // cout << "Y1 = " << Y1 << " - Y2 = " << Y2 << endl;
  // cout << "Xout = " << Xout << endl;
  // cout << endl;
  // double Xin = Xout;
  // Yin = (Y2-Y1)/(X2-X1)*(Xin-X1) + Y1;
  // cout << "Yin = " << Yin << endl;
 
  // TGraph *gr = new TGraph();
  // gr->SetPoint(0,X1,Y1);
  // gr->SetPoint(1,X2,Y2);
  // cout << "Eval = " << gr->Eval(Xin) << endl;
  // return gr->Eval(Xin);
}


double CalcDLcrit_posterior(vector<double> dchisq, vector<double> weight, double conflev, double sigmacl){ 

  // Get error on coverage: sigmacl = sqrt(N*p*(1-p))/N

  int nelem = dchisq.size();
  
  double *new_dchisq = new double[nelem];
  double *new_weight = new double[nelem];
  int *index = new int[nelem];
  
  // copy vector into array
  for(int i=0;i<nelem;i++){
    new_dchisq[i] = dchisq[i];
    new_weight[i] = weight[i];
    index[i] = i;
  }
  
  sort(dchisq.begin(), dchisq.end());
  
  // store the sorted indices --> use for both dchisq and weights
  //index = TMath::Sort(nelem,new_dchisq,index,0);
  TMath::Sort(nelem,new_dchisq,index,0);

  double tot_weight = 0.;
  for(int i=0;i<nelem;i++) tot_weight += new_weight[index[i]]; 

  //cout << "conflev = " << conflev << endl;
  
  double sum = 0.;
  double critval = -999.;
  for(int i=0;i<nelem;i++){
    sum += new_weight[index[i]]; 
    
    //cout << std::setprecision(16);
    //cout << sum << " ";
    
    double epsilon = 0.000001; // for the double precision in the break condition

    if( sum >= ( (0.01 * conflev + sigmacl)*tot_weight - epsilon ) ){
      
      // cout << endl;
      // cout << "sum = " << sum << endl;
      // cout << "new_weight[index[i]] = " << new_weight[index[i]] << endl;
      // cout << "before: sum = " << (sum - new_weight[index[i]]) << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      // cout << "current: sum = " << sum << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      // cout << "after: sum = " << (sum + new_weight[index[i]]) << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      
      critval = new_dchisq[index[i]]; 
      break;
    }
  }
    
  return critval;
}


double CalcDLcrit(vector<double> dchisq, double conflev, double sigmacl){

  // Get error on critical value: sigmacl = sqrt(N*p*(1-p))

  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector   
  vector<double>::iterator dchisqIter;
  dchisqIter = dchisq.begin() + (int)(0.01*conflev*dchisq.size() + sigmacl);
  double critval = (*dchisqIter);
  return critval;
}

double CalcDLcrit(TH1D *h_delta_true,double FCfrac_CL){
  
  //Feldman-Cousins: Calculate chi2 value at 90%CL************************
  double ntoys_total = 0;
  double ntoys = 0;
  double ntoys_frac = 0;
  int chi2_bin = 0;
  double chi2_val = 0;
  double dlval = 0;
  
  ntoys_total = h_delta_true->Integral();

  for(int ibin=0;ibin<h_delta_true->GetNbinsX();ibin++){
    ntoys += h_delta_true->GetBinContent(ibin+1);
    ntoys_frac = ntoys/ntoys_total;
    
    if(ntoys_frac>=FCfrac_CL) {
      chi2_val = h_delta_true->GetBinCenter(ibin+1);
      break;
    }
  }
  return chi2_val;
}


double GetPvalue_rightside(vector<double> vec_toys,double data){

  double ntot_toys = vec_toys.size();

  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it >= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}

double GetPvalue_posterior_rightside(vector<double> vec_toys,vector<double> weight,double data){

  double ntot_toys = vec_toys.size();
  
  double tot = 0.; 
  double nfrac=0;
  for(int i=0;i<vec_toys.size();i++){
    tot += weight[i];
    double cont = vec_toys[i];
    if(cont >= data) nfrac += weight[i];
  }
  double pvalue =  nfrac / tot;

  return pvalue;
}
