//
// Use this macro to draw 60% and 90% confidence regions from an experiment,
// macro assumes 1 experiment per root file and dumps the output
// into a root file and PDF.
//
// Give it the location of the root file, witout the .root at the end eg. $dir/toyA
//
// If the input contains either "ss23" or "cp13" in the name, the bf point will be  added

void draw_contour
(
 //const char tag[200] = "fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC12p482_RHC12p271_react_ih",
 //const char tag[200] = "fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC14p882_RHC9p871_react_ih",
 //const char tag[200] = "fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC17p282_RHC7p471_react_ih",
 const char tag[200] = "fitmarg_joint_sinsqth23Vsdmsq_Run17c_t2konly_nh",
 

 //////////////////////////////////////////////
 //                                          //
 //         Run 1-7c - Data - Final          // 
 //                                          //
 //////////////////////////////////////////////

 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Data_Run17c_nexp_1_0.450_0.003-job0-_0_0.root", // NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Data_Run17c_nexp_1_0.450_0.003-job0-_0_0.root", // IH - reactors

 const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth23Vsdmsq_t2konly_nh_Data_Run17c_nexp_1_0.528_0.003-job0-_0_0.root", // NH - T2K only
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth23Vsdmsq_t2konly_ih_Data_Run17c_nexp_1_0.528_0.003-job0-_0_0.root", // IH - T2K only
 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_Data_Run17c_nexp_1_0.022_0.000-job0-_0_0.root", // NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_Data_Run17c_nexp_1_0.022_0.000-job0-_0_0.root", // IH - reactors 
 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_nh_Data_Run17c_nexp_1_0.450_0.000-job0-_0_0.root", // NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_ih_Data_Run17c_nexp_1_0.450_0.000-job0-_0_0.root", // IH - reactors
 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_dmsqVsdcp_react_nh_Data_Run17c_nexp_1_0.003_0.000-job0-_0_0.root", // NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/fit_fitmarg_joint_dmsqVsdcp_react_ih_Data_Run17c_nexp_1_0.003_0.000-job0-_0_0.root", // IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/datafit/cp13/data_fit_nh_a-_0_0.root", // NH - T2K only
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17c/files/datafit/cp13/data_fit_ih_a-_0_0.root", // IH - T2K only
 

 //////////////////////////////////////////////
 //                                          //
 //        Run 1-7c - Asimov - Final         // 
 //                                          //
 //////////////////////////////////////////////

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_Asimov_A_Run17c_nexp_1_0.022_-1.601-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_Asimov_A_Run17c_nexp_1_0.022_-1.601-job0-_0_0.root", // Asimov A - IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_Asimov_B_Run17c_nexp_1_0.022_0.000-job0-_0_0.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_Asimov_B_Run17c_nexp_1_0.022_0.000-job0-_0_0.root", // Asimov B - IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Asimov_A_Run17c_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Asimov_A_Run17c_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Asimov_B_Run17c_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Asimov_B_Run17c_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - IH - reactors
 
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Asimov_A_FHC12p482_RHC12p271_nexp_1_0.528_0.003-job0-_0_0.root", // 
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Asimov_A_FHC12p482_RHC12p271_nexp_1_0.528_0.003-job0-_0_0.root", // 
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Asimov_A_FHC14p882_RHC9p871_nexp_1_0.528_0.003-job0-_0_0.root", // 
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Asimov_A_FHC14p882_RHC9p871_nexp_1_0.528_0.003-job0-_0_0.root", // 
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Asimov_A_FHC17p282_RHC7p471_nexp_1_0.528_0.003-job0-_0_0.root", //  
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/HighStat/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Asimov_A_FHC17p282_RHC7p471_nexp_1_0.528_0.003-job0-_0_0.root", // 

 //////////////////////////////////////////////
 //                                          //
 //     Run 1-7b - Data - Run 1-7b Flux      // 
 //                                          //
 //////////////////////////////////////////////
 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_Data_Run17_nexp_1_0.022_-1.601-job0-_0_0.root", 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_Data_Run17_nexp_1_0.022_-1.601-job0-_0_0.root", 

 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/data_fit_sinsqth13Vsdcp_t2konly_nh_run17.root", 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/data_fit_sinsqth13Vsdcp_t2konly_ih_run17.root", 

 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dmsqVsdcp_react_nh_Data_Run17_nexp_1_0.003_-1.601-job0-_0_0.root", 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_dmsqVsdcp_react_ih_Data_Run17_nexp_1_0.003_-1.601-job0-_0_0.root", 

 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Data_Run17_nexp_1_0.528_0.003-job0-_0_0.root", 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Data_Run17_nexp_1_0.528_0.003-job0-_0_0.root", 

 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_nh_Data_Run17_nexp_1_0.528_-1.601-job0-_0_0.root", 
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN17b/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_ih_Data_Run17_nexp_1_0.528_-1.601-job0-_0_0.root", 


 //////////////////////////////////////////////
 //                                          //
 //     Run 1-7b -Asimov - Run 1-7b Flux     // 
 //                                          //
 //////////////////////////////////////////////

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_Asimov_A_Run17_nexp_1_0.022_-1.601-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_Asimov_A_Run17_nexp_1_0.022_-1.601-job0-_0_0.root", // Asimov A - IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_Asimov_B_Run17_nexp_1_0.022_0.000-job0-_0_0.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_Asimov_B_Run17_nexp_1_0.022_0.000-job0-_0_0.root", // Asimov B - IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_nh_Asimov_A_Run17_nexp_1_0.528_-1.601-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_ih_Asimov_A_Run17_nexp_1_0.528_-1.601-job0-_0_0.root", // Asimov A - IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Asimov_A_Run17_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Asimov_A_Run17_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Asimov_B_Run17_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Asimov_B_Run17_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/asimov_fit_nh_a_sinsqth13Vsdcp_t2konly_run17.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/asimov_fit_ih_a_sinsqth13Vsdcp_t2konly_run17.root", // Asimov A - IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/asimov_fit_nh_b_sinsqth13Vsdcp_t2konly_run17.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/files/asimov_fit_ih_b_sinsqth13Vsdcp_t2konly_run17.root", // Asimov B - IH - reactors


 /////////////////////////////////
 //                             //
 // Run 1-6 -Asimov NO FIT BIAS // 
 //                             //
 /////////////////////////////////

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/files/fit_fitmarg_joint_sinsqth13Vsdcp_nofitbias_Asimov_A_react_nh_nexp_1_0.022_-1.601-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/files/fit_fitmarg_joint_sinsqth13Vsdcp_nofitbias_Asimov_A_react_ih_nexp_1_0.022_-1.601-job0-_0_0.root", // Asimov A - IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/files/fit_fitmarg_joint_sinsqth13Vsdcp_nofitbias_Asimov_B_react_nh_nexp_1_0.022_0.000-job0-_0_0.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/files/fit_fitmarg_joint_sinsqth13Vsdcp_nofitbias_Asimov_B_react_ih_nexp_1_0.022_0.000-job0-_0_0.root", // Asimov A - IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/files/fit_fitmarg_joint_sinsqth23Vsdmsq_nofitbias_Asimov_A_react_nh_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/files/fit_fitmarg_joint_sinsqth23Vsdmsq_nofitbias_Asimov_A_react_ih_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/files/fit_fitmarg_joint_sinsqth23Vsdmsq_nofitbias_Asimov_B_react_nh_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/files/fit_fitmarg_joint_sinsqth23Vsdmsq_nofitbias_Asimov_B_react_ih_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov A - IH - reactors

 //////////////////////////////
 //                          //
 //      Run 1-6 -Asimov     // 
 //                          //
 //////////////////////////////

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_A_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_A_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_B_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_B_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - IH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth23Vsdmsq_statonly_react_nh_A_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - NH - reactors - statonly
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth23Vsdmsq_statonly_react_ih_A_nexp_1_0.528_0.003-job0-_0_0.root", // Asimov A - IH - reactors - statonly
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth23Vsdmsq_statonly_react_nh_B_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - NH - reactors - statonly
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth23Vsdmsq_statonly_react_ih_B_nexp_1_0.450_0.003-job0-_0_0.root", // Asimov B - IH - reactors - statonly

 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_A_nexp_1_0.022_-1.601-job0-_0_0.root", // Asimov A - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_A_nexp_1_0.022_-1.601-job0-_0_0.root", // Asimov A - IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_B_nexp_1_0.022_0.000-job0-_0_0.root", // Asimov B - NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_B_nexp_1_0.022_0.000-job0-_0_0.root", // Asimov B - IH - reactors
 
 //////////////////////////////
 //                          //
 //      Run 1-6 - Data      // 
 //                          //
 //////////////////////////////

 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_nh_Data_Run16_nexp_1_0.528_0.003-job0-_0_0.root", // NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23Vsdmsq_react_ih_Data_Run16_nexp_1_0.528_0.003-job0-_0_0.root", // IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_nh_Data_Run16_nexp_1_0.022_-1.601-job0-_0_0.root", // NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth13Vsdcp_react_ih_Data_Run16_nexp_1_0.022_-1.601-job0-_0_0.root", // IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_nh_Data_Run16_nexp_1_0.528_-1.601-job0-_0_0.root", // NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_ih_Data_Run16_nexp_1_0.528_-1.601-job0-_0_0.root", // IH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_nh_Data_Run16_like_nexp_1_0.528_-1.601-job0-_0_0.root", // NH - reactors
 //const char * filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_sinsqth23Vsdcp_react_nh_fullchisq_nexp_1_-999.000_0.000-job0-_0_0.root",
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dmsqVsdcp_react_nh_Data_Run16_nexp_1_0.003_-1.601-job0-_0_0.root", // NH - reactors
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dmsqVsdcp_react_ih_Data_Run16_nexp_1_0.003_-1.601-job0-_0_0.root", // NH - reactors

 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN16/run_16_datafit_Raj/cp13/datafit_dcp_nh-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // NH - T2K only
 //const char* filename = "../TECH_NOTE_14_3_16/DATA/RUN16/run_16_datafit_Raj/cp13/datafit_dcp_ih-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // IH - T2K only
 


 //const bool doRasterScan = false,
 const char* opt = "cp13", //"cp13","ss23"
 const bool SavePDF = false
 )
{
  //   char* filename = "toyA_joint";
  
  TCanvas* c = new TCanvas("c","c",0,0,900,600);
  TFile f(Form("%s",filename));
  TTree* t = f.Get("fitnt");
  
  TH2D* h1;
  TH2D* h2;
  TH2D* h3;
  TH2D* h4;
  TH2D *h1_rasterscan;
  TH2D *h2_rasterscan;

  valor::FitParamList* fp = 0;
  valor::PhysParamList* ip = 0;
  
  t->SetBranchAddress("hchisq2D",&h1);
  t->SetBranchAddress("best_fit_params",&fp);
  t->SetBranchAddress("input_params",&ip);    
  t->GetEntry(0);
  
  f.Close();
  
  
  string file2 = filename;
  string option = opt;
  double input_x,input_y;
  if(file2.find("ss23")!= string::npos || option.find("ss23")!= string::npos )
    {
      input_x = ip->Value(valor::conv::kPhypMixing23);
      input_y = ip->Value(valor::conv::kPhypDm23sq);
    }
  else if ( file2.find("cp13")!= string::npos || option.find("cp13")!= string::npos)
    {
      input_x = ip->Value(valor::conv::kPhypMixing13);
      input_y = ip->Value(valor::conv::kPhypCPV);
    }
  else{
    cout << "\n\n\tERROR: Could not find which parameters are being used\n\n" <<endl;
    input_x = 0;
    input_y = 0;
  }
  
  cout << "x: " << input_x << endl;
  cout << "y: " << input_y << endl;
  
  
  h3 = (TH2D*)h1->Clone();
  h2 = (TH2D*)h1->Clone();
  h4 = (TH2D*)h1->Clone();
  //h1_rasterscan = (TH2D*)h1->Clone("h1_rasterscan");
  double chisq_min = h1->GetBinContent(h1->GetMinimumBin());

  cout << "chisq_min = " << chisq_min << endl;

  h1->SetContour(1);
  h3->SetContour(1);
  
  h1->SetContourLevel(0,chisq_min+2.3);
  h3->SetContourLevel(0,chisq_min+4.61);
  h1->SetLineWidth(2);
  h3->SetLineWidth(2);
  h1->SetLineColor(kRed);
  h3->SetLineColor(kRed);
  h1->SetLineStyle(2);
  



  int minbin = h1->GetMinimumBin();
  int binx;
  int biny;
  h1->GetBinXYZ(minbin,binx,biny,0);
  
  double xmin = h1->GetXaxis()->GetBinCenter(binx);
  double ymin = h1->GetYaxis()->GetBinCenter(biny);
  
  cout << "best-fit x = " << xmin << endl;
  cout << "best-fit y = " << ymin << endl;

  TMarker * bf = new TMarker(xmin,ymin,kFullStar);
  bf->SetMarkerColor(kRed);
  bf->SetMarkerSize(1.5);
  
  TMarker * bf2 = new TMarker(input_x,input_y,kFullStar);
  bf2->SetMarkerColor(kViolet);
  bf2->SetMarkerSize(1.5);
  
  TLegend* leg = new TLegend(0.7,0.8,0.9,0.98);
  leg->AddEntry(bf,"Best fit point","p");
  //leg->AddEntry(bf2,"True asimov point","p");
  leg->AddEntry(h1,"68%","l");
  leg->AddEntry(h3,"90%","l");
  h2->Draw("colz");
  h1->Draw("cont3same");
  h3->Draw("cont3same");
  bf->Draw();
  bf2->Draw();
  
  leg->Draw();



  TCanvas *cCont = new TCanvas("cCont","cCont");
  h1->Draw("cont3");
  h3->Draw("cont3 same");
  bf->Draw();
  //bf2->Draw();
  
  
  // Produce raster scan  
  // TCanvas *cRasterScan = new TCanvas("cRasterScan","cRasterScan");
  // if(doRasterScan){
  //   TGraph *bf_raster = new TGraph();
  //   for(int iy=0;iy<h1->GetYaxis()->GetNbins();iy++){
  //     // Fine the minimum for each value of y
  //     double raster_min = 99999.;
  //     for(int ix=0;ix<h1->GetXaxis()->GetNbins();ix++){
  // 	double cont = h1_rasterscan->GetBinContent(ix+1,iy+1);
  // 	if(cont < raster_min) raster_min = cont;
  //     }
  //     cout << raster_min << " ";
  //     // Shift wrt the minimum for each value of y
  //     for(int ix=0;ix<h1->GetXaxis()->GetNbins();ix++){
  // 	double cont = h1_rasterscan->GetBinContent(ix+1,iy+1);
  // 	//cout << "cont = " << cont << " - raster_min = " << raster_min << endl;
  // 	if(cont == raster_min){
  // 	  double cont_x = h1_rasterscan->GetXaxis()->GetBinCenter(ix+1);
  // 	  double cont_y = h1_rasterscan->GetYaxis()->GetBinCenter(iy+1);
  // 	  bf_raster->SetPoint(iy,cont_x,cont_y);
  // 	}
  // 	cont = cont - raster_min;
  // 	h1_rasterscan->SetBinContent(ix+1,iy+1,cont);
  //     }
  //   }
  //   //h1_rasterscan->DrawCopy("colz");
  //   h2_rasterscan = (TH2D*)h1_rasterscan->Clone("h2_rasterscan");
  //   h1_rasterscan->SetContour(1);
  //   h2_rasterscan->SetContour(1);
  //   h1_rasterscan->SetContourLevel(0,1.0);
  //   h2_rasterscan->SetContourLevel(0,2.7);
  //   h1_rasterscan->SetLineStyle(2);
  //   h1_rasterscan->SetLineWidth(2);
  //   h1_rasterscan->SetLineColor(kRed);
  //   h2_rasterscan->SetLineWidth(2);
  //   h2_rasterscan->SetLineColor(kRed);
  //   h1_rasterscan->Draw("cont3");
  //   h2_rasterscan->Draw("cont3 same");
  //    bf_raster->Draw("p");
  //   bf_raster->SetLineColor(kRed);
  //   bf_raster->SetMarkerColor(kRed);
  //   bf_raster->SetMarkerStyle(kFullStar);
  // } // doRasterScan
  // // Save pdf files
  //if(SavePDF) c->SaveAs(Form("%s.pdf",tag));
  if(SavePDF) cCont->SaveAs(Form("contour_%s.pdf",tag));
  //if(SavePDF) cRasterScan->SaveAs(Form("rasterscan_%s.pdf",tag));
  

  // Save contours in root file

  TFile* f2 = new TFile(Form("%s_contour.root",tag),"recreate");
  
  save_contour(f2,h1,h3);
  h2->Write("chisq_map_2D");
  
  bf->Write("best_fit_point_marker");
  f2->Close();
  
}

void save_contour(TFile* f2, TH2D* h1=0, TH2D* h3=0)
{
    TCanvas* c2 = new TCanvas("c2","c",0,0,900,600);
    h3->Draw("CONT LIST");
    c2->Update();
    TObjArray *contours = 0;
    contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
    
    TList *list90 = (TList*)contours->At(0);
    TGraph *g;
    int it=0;
    TIter next(list90);
    while ((g = (TGraph*)next()))
    {
        char name[200];
        sprintf(name, "g90_valor_%d",it);
        g->SetLineWidth(2);
        g->SetLineStyle(1);
        g->SetLineColor(kRed);
        g->Write(name);
        it++;
    }
    
    c2->Clear();
    h1->Draw("CONT LIST");
    c2->Update();
    contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList *list68 = (TList*)contours->At(0);
    TIter next68(list68);
    it=0;
    while ((g = (TGraph*)next68()))
    {
        char name[200];
        sprintf(name, "g68_valor_%d",it);
        g->SetLineWidth(2);
        g->SetLineStyle(2);
        g->SetLineColor(kRed);
        g->Write(name);
        it++;
    }
    
    delete c2;
    delete contours;
    delete list68;
    delete list90;
    delete g;
    
}
