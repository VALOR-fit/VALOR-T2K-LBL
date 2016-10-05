
// Produce VALOR contours with macros/misc/draw_contours.C
// and read the output files

#define SwapXY // Swap x and y axis of contours // does this to the chisq2d, contour and the TMarker!

using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareContours()
{
  const bool scalegraphX = false; // scale axis of graphs (e.g. dmsq x 10+3)
  const bool scalegraphY = false; // scale axis of graphs (e.g. dmsq x 10+3)

  const bool setinit = true; // true --> set init histo by hand, false --> use already defined set up
  
  // // sinsqth23 Vs dmsq
  // const double x_low = 0.3;
  // const double x_up = 0.7;
  // const double y_low = 1.95;
  // const double y_up = 3.05;
  // const char *title_x = "sin^{2} #theta_{23}";
  // const char *title_y = "|#Delta m^{2}_{32}| (10^{-3} eV^{2}/c^{4})";
  
  // // sinsqth23 Vs dcp (reactors)
  // const double x_low = 0.3;
  // const double x_up = 0.7;
  // const double y_low = -TMath::Pi();
  // const double y_up = +TMath::Pi();
  // const char *title_x = "sin^{2} #theta_{23}";
  // const char *title_y = "#delta_{CP} (radians)";
  
  // // sinsqth13 Vs dcp (reactors)
  // const double x_low = 0.014999;
  // const double x_up = 0.0360001;
  // const double y_low = -TMath::Pi();
  // const double y_up = +TMath::Pi();
  // const char *title_x = "sin^{2} #theta_{13}";
  // const char *title_y = "#delta_{CP} (radians)";

  // // sinsqth13 Vs dmsq (reactors)
  // const double x_low = 1.95;
  // const double x_up = 3.05;
  // const double y_low = -TMath::Pi();
  // const double y_up = +TMath::Pi();
  // const char *title_x = "|#Delta m^{2}_{32}| (10^{-3} eV^{2}/c^{4})";
  // const char *title_y = "#delta_{CP} (radians)";

  // sinsqth13 Vs dcp (t2k only)
  const double x_low = 0.;
  const double x_up = 0.105;
  const double y_low = -TMath::Pi();
  const double y_up = +TMath::Pi();
  const char *title_x = "sin^{2} #theta_{13}";
  const char *title_y = "#delta_{CP} (radians)";

  const bool draw_preliminary = false;
  TString str_preliminary = "T2K Run1-7c preliminary";

  bool draw_reactor_band = false;

  //////////////////////////// NEW VALOR contours /////////////////////
  
  // Global scan
  
  const int Nfiles = 4;
  const char *inname[Nfiles] = {


    //////////////////////////////////////////////
    //                                          //
    //         Run 1-7c - Data - Final         // 
    //                                          //
    //////////////////////////////////////////////

    // // sinsqth13 Vs dcp w/ reactors - Data - Run 1-7c
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17c_Data_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17c_Data_react_ih_contour.root",

    // // sinsqth23 Vs dmsq w/ reactors - Data - Run 1-7c
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17c_Data_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17c_Data_react_ih_contour.root",
    
    // // sinsqth23 Vs dcp w/ reactors - Data - Run 1-7c
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth23Vsdcp_react/fitmarg_joint_sinsqth23Vsdcp_Run17c_Data_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth23Vsdcp_react/fitmarg_joint_sinsqth23Vsdcp_Run17c_Data_react_ih_contour.root",
    
    // // dmsq Vs dcp w/ reactors - Data - Run 1-7c
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_dmsqVsdcp_react/fitmarg_joint_dmsqVsdcp_Run17c_Data_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_dmsqVsdcp_react/fitmarg_joint_dmsqVsdcp_Run17c_Data_react_ih_contour.root",

    // // sinsqth13 Vs dcp T2K only  - Data - Run 1-7c
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth13Vsdcp_t2konly/fitmarg_joint_sinsqth13Vsdcp_Run17c_Data_t2konly_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth13Vsdcp_t2konly/fitmarg_joint_sinsqth13Vsdcp_Run17c_Data_t2konly_ih_contour.root",


    //////////////////////////////////////////////
    //                                          //
    //         Run 1-7c -Asimov - Final         // 
    //                                          //
    //////////////////////////////////////////////

    // // sinsqth13 Vs dcp w/ reactors - Asimov A - Run 1-7c
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17c_Asimov_A_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17c_Asimov_A_react_ih_contour.root",

    // // sinsqth13 Vs dcp w/ reactors - Asimov B - Run 1-7c
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17c_Asimov_B_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17c_Asimov_B_react_ih_contour.root",
    
    // // sinsqth23 Vs dmsq w/ reactors - Asimov A - Run 1-7c
    //"../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17c_Asimov_A_react_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17c_Asimov_A_react_ih_contour.root",
    
    // // sinsqth23 Vs dmsq w/ reactors - Asimov B - Run 1-7c
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17c_Asimov_B_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17c_Asimov_B_react_ih_contour.root",

    // // sinsqth13 Vs dcp T2K only  - Asimov A - Run 1-7c
    //"../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/asimovs_run17c_Raj/cp13/asimov_fit_nh_a-_0_0_contour.root",
    "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/asimovs_run17c_Raj/cp13/asimov_fit_ih_a-_0_0_contour.root",

    // // sinsqth13 Vs dcp T2K only  - Asimov B - Run 1-7c
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/asimovs_run17c_Raj/cp13/asimov_fit_nh_b-_0_0_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/files/asimovs_run17c_Raj/cp13/asimov_fit_ih_b-_0_0_contour.root",

    // High Stat (sensi)
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC12p482_RHC12p271_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC14p882_RHC9p871_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC17p282_RHC7p471_react_nh_contour.root",

    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC12p482_RHC12p271_react_ih_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC14p882_RHC9p871_react_ih_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_FHC17p282_RHC7p471_react_ih_contour.root",

    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/ss13Vsdcp_t2konly/asimov_fit_12_12_nh-_0_0_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/ss13Vsdcp_t2konly/asimov_fit_14_9_nh-_0_0_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/ss13Vsdcp_t2konly/asimov_fit_17_7_nh-_0_0_contour.root",

    "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/ss13Vsdcp_t2konly/asimov_fit_12_12_ih-_0_0_contour.root",
    "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/ss13Vsdcp_t2konly/asimov_fit_14_9_ih-_0_0_contour.root",
    "../TECH_NOTE_14_3_16/SENSI/RUN17c/final/contours_highstat/ss13Vsdcp_t2konly/asimov_fit_17_7_ih-_0_0_contour.root",

    //////////////////////////////////////////////
    //                                          //
    //     Run 1-7b - Data - Run 1-7b Flux      // 
    //                                          //
    //////////////////////////////////////////////

    // // sinsqth13 Vs dcp T2K only - Run 1-7b
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run17_Data_t2konly_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run17_Data_t2konly_ih_contour.root",

    // // dmsq Vs dcp reactors - Run 1-7b
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_dmsqVsdcp/fitmarg_joint_dmsqVsdcp_Run17_Data_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_dmsqVsdcp/fitmarg_joint_dmsqVsdcp_Run17_Data_react_ih_contour.root",

    // // sinsqth23 Vs dmsq reactors - Run 1-7b
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth23Vsdmsq/fitmarg_joint_sinsqth23Vsdmsq_Run17_Data_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth23Vsdmsq/fitmarg_joint_sinsqth23Vsdmsq_Run17_Data_react_ih_contour.root",

    // // sinsqth23 Vs dcp reactors - Run 1-7b
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth23Vsdcp/fitmarg_joint_sinsqth23Vsdcp_Run17_Data_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth23Vsdcp/fitmarg_joint_sinsqth23Vsdcp_Run17_Data_react_ih_contour.root",
    
    // // sinsqth13 Vs dcp T2K reactors - Run 1-7b
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run17_Data_react_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run17_Data_react_ih_contour.root",

    
  
    //////////////////////////////////////////////
    //                                          //
    //     Run 1-7b -Asimov - Run 1-7b Flux     // 
    //                                          //
    //////////////////////////////////////////////

    // // sinsqth13 Vs dcp w/ reactors - Asimov A - Run 1-7b
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_A_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_A_reactors_ih_contour.root",

    // // sinsqth13 Vs dcp w/ reactors - Asimov B - Run 1-7b
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_B_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_B_reactors_ih_contour.root",

    // // sinsqth23 Vs dmsq w/ reactors - Asimov A - Run 1-7b
    //"../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17_Asimov_A_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17_Asimov_A_reactors_ih_contour.root",

    // // sinsqth23 Vs dmsq w/ reactors - Asimov B - Run 1-7b
    //"../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17_Asimov_B_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17_Asimov_B_reactors_ih_contour.root",

    // // sinsqth13 Vs dcp T2K only  - Asimov A - Run 1-7b
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth13Vsdcp_t2konly/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_A_t2konly_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth13Vsdcp_t2konly/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_A_t2konly_ih_contour.root",

    // // sinsqth13 Vs dcp T2K only  - Asimov B - Run 1-7b
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth13Vsdcp_t2konly/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_B_t2konly_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17b/FluxRun17_Conference/contours_sinsqth13Vsdcp_t2konly/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_B_t2konly_ih_contour.root",
 

    //////////////////////////////
    //                          //
    //      Run 1-6 -Asimov     // 
    //                          //
    //////////////////////////////

    // // sinsqth13 Vs dcp w/ reactors - Asimov A
    // "../TECH_NOTE_14_3_16/SENSI/RUN16/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Asimov_A_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN16/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Asimov_A_reactors_ih_contour.root",

    // sinsqth13 Vs dcp w/ reactors - Asimov B
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Asimov_B_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Asimov_B_reactors_ih_contour.root",
    
    // // sinsqth23 Vs dmsq w/ reactors - Asimov A
    // "../TECH_NOTE_14_3_16/SENSI/RUN16/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN16/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_reactors_ih_contour.root",
    
    // // sinsqth23 Vs dmsq w/ reactors - Asimov B
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Asimov_B_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Asimov_B_reactors_ih_contour.root",

    // // sinsqth13 Vs dcp T2K only - Asimov A
    // "../TECH_NOTE_14_3_16/SENSI/RUN16/asimovs_16_Raj/cp13/toyfits_dcp_anh-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN16/asimovs_16_Raj/cp13/toyfits_dcp_aih-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_contour.root",

    // // sinsqth13 Vs dcp T2K only - Asimov B
    // "../TECH_NOTE_14_3_16/SENSI/RUN16/asimovs_16_Raj/cp13/toyfits_dcp_bnh-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN16/asimovs_16_Raj/cp13/toyfits_dcp_bih-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_contour.root",

    
    

    //////////////////////////////
    //                          //
    //      Run 1-6 - Data      // 
    //                          //
    //////////////////////////////

    // // sinsqth23 Vs dmsq w/ reactors - Data
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run16_Data_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run16_Data_reactors_ih_contour.root"

    // // sinsqth13 Vs dcp T2K only - Data
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_data_t2konly_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_data_t2konly_ih_contour.root"

    // // sinsqth13 Vs dcp w/ reactors - Data
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run16_Data_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run16_Data_reactors_ih_contour.root"

    // // dmsq Vs dcp w/ reactors - Data
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_dmsqVsdcp_react/fitmarg_joint_dmsqVsdcp_Run16_Data_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_dmsqVsdcp_react/fitmarg_joint_dmsqVsdcp_Run16_Data_reactors_ih_contour.root"

    // // sinsqth23 Vs dcp w/ reactors - Data
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_sinsqth23Vsdcp_react/fitmarg_joint_sinsqth23Vsdcp_Run16_Data_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/RUN16/contours_sinsqth23Vsdcp_react/fitmarg_joint_sinsqth23Vsdcp_Run16_Data_reactors_ih_contour.root"

    /////////////// NO FIT BIAS ////////////

    // sinsqth13 Vs dcp T2K reactors - Asimov A 
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth13Vsdcp_NoFitBias_Run16_Asimov_A_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth13Vsdcp_NoFitBias_Run16_Asimov_A_reactors_ih_contour.root"

    // sinsqth13 Vs dcp T2K reactors - Asimov B 
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth13Vsdcp_NoFitBias_Run16_Asimov_B_reactors_nh_contour.root"
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth13Vsdcp_NoFitBias_Run16_Asimov_B_reactors_ih_contour.root"

    // sinsqth23 Vs dmsq T2K reactors - Asimov A 
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth23Vsdmsq_NoFitBias_Run16_Asimov_A_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth23Vsdmsq_NoFitBias_Run16_Asimov_A_reactors_ih_contour.root"

    // sinsqth23 Vs dmsq T2K reactors - Asimov B
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth23Vsdmsq_NoFitBias_Run16_Asimov_B_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/RUN16/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth23Vsdmsq_NoFitBias_Run16_Asimov_B_reactors_ih_contour.root"
    
  };
  const int Ncont_68 = 2;
  const char *contour_68_name[Ncont_68] = { // different graph for each contour (+-)
    "g68_valor_0",
    "g68_valor_1"
  };
  const int Ncont_90 = 2;
  const char *contour_90_name[Ncont_90] = { // different graph for each contour (+-)
    "g90_valor_0",
    "g90_valor_1"
  };
  const char *latexname[Nfiles] = {
    //"Normal Hierarchy",
    //"Inverted Hierarchy",
    //"Run 1-7c",
    "#nu: 7.48e20 + #bar{#nu}: 7.47e20 (Run 1-7c)",
    "#nu: 12.5e20 + #bar{#nu}: 12.3e20",
    "#nu: 14.9e20 + #bar{#nu}: 9.9e20", 
    "#nu: 17.3e20 + #bar{#nu}: 7.5e20"
  };
  const char *chisq_name = "chisq_map_2D";
  
  const char *bestfitpt_name = "best_fit_point_marker";
    

  ////////////////////////////////////////////////////////////////////////
  
  
  // Read contours
  
  TFile *filein[Nfiles];
  TGraph *g68_valor[Ncont_68][Nfiles] = {0};
  TGraph *g90_valor[Ncont_90][Nfiles] = {0};
  TH2D *hchisq_valor[Nfiles] = {0};
  TH2D *hchisq_glob_valor[Nfiles] = {0};
  TH2D *hdeltachisq_glob_valor[Nfiles] = {0};
  TH2D *hdeltachisq_glob_valor_68[Nfiles] = {0};
  TH2D *hdeltachisq_glob_valor_90[Nfiles] = {0};
  TMarker *gbest_fit_global[Nfiles] = {0};

  for(int ifile=0;ifile<Nfiles;ifile++){

    filein[ifile] = new TFile(inname[ifile],"READ");
    
    for(int icont=0;icont<Ncont_68;icont++){
    
      g68_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_68_name[icont]);
#ifdef SwapXY
      if(g68_valor[icont][ifile] != NULL)
      SwapXY(g68_valor[icont][ifile]);
#endif
    }
    for(int icont=0;icont<Ncont_90;icont++){
      g90_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_90_name[icont]);
#ifdef SwapXY
      if(g90_valor[icont][ifile] != NULL)
      SwapXY(g90_valor[icont][ifile]);
#endif
    }
    hchisq_valor[ifile] = (TH2D*) filein[ifile]->Get(chisq_name);
    if(hchisq_valor[ifile] == NULL){
      std::cout << "Histogram not found! Exiting" << std::endl;
      exit(0);
    }
    
    hchisq_glob_valor[ifile] = (TH2D*) hchisq_valor[ifile]->Clone(Form("hchisq_glob_%d",ifile));

    hchisq_valor[ifile]->Reset(); // used to set the full range though tgraphs
#ifdef SwapXY
    hchisq_valor[ifile] = SwapXY(hchisq_valor[ifile]);
    hchisq_glob_valor[ifile] = SwapXY(hchisq_glob_valor[ifile]);
#endif
    gbest_fit_global[ifile] = (TMarker*)filein[ifile]->Get(bestfitpt_name);
#ifdef SwapXY
    SwapXY(gbest_fit_global[ifile]);
#endif
  }




  // Init histogram

  TH2D *hInit;
  if(setinit){
    hInit = new TH2D("hInit","hInit",201,x_low,x_up,201,y_low,y_up);
    hInit->GetXaxis()->SetTitle(title_x);
    hInit->GetYaxis()->SetTitle(title_y);
  }
  else{
    int NbinsX = hchisq_glob_valor[0]->GetNbinsX(); 
    double lowedgeX = hchisq_glob_valor[0]->GetXaxis()->GetBinLowEdge(1);
    double upedgeX = hchisq_glob_valor[0]->GetXaxis()->GetBinUpEdge(NbinsX);
    int NbinsY = hchisq_glob_valor[0]->GetNbinsY();
    double lowedgeY = hchisq_glob_valor[0]->GetYaxis()->GetBinLowEdge(1);
    double upedgeY = hchisq_glob_valor[0]->GetYaxis()->GetBinUpEdge(NbinsY);
    char *TitleX = hchisq_glob_valor[0]->GetXaxis()->GetTitle();
    char *TitleY = hchisq_glob_valor[0]->GetYaxis()->GetTitle();
    TH2D *hInit = new TH2D("hInit","hInit",NbinsX,lowedgeX,upedgeX,NbinsY,lowedgeY,upedgeY);
    hInit->GetXaxis()->SetTitle(TitleX);
    hInit->GetYaxis()->SetTitle(TitleY);

    cout << "lowedgeX = " << lowedgeX << endl;
    cout << "upedgeX = " << upedgeX << endl;
    cout << "NbinsX = " << NbinsX << endl;
    cout << "first X bin = " << hInit->GetXaxis()->GetBinCenter(1) << endl;
    cout << "last X bin = " << hInit->GetXaxis()->GetBinCenter(NbinsX) << endl;
    cout << "lowedgeY = " << lowedgeY << endl;
    cout << "upedgeY = " << upedgeY << endl;
    cout << "NbinsY = " << NbinsY << endl;
    cout << "first Y bin = " << hInit->GetYaxis()->GetBinCenter(1) << endl;
    cout << "last Y bin = " << hInit->GetYaxis()->GetBinCenter(NbinsY) << endl;
  }


  // Get global best-fit contour

  double chisq_min_glob = 99999999.;
  double chisq_min[Nfiles] = 0.;
  for(int ifile=0;ifile<Nfiles;ifile++){
    chisq_min[ifile] = hchisq_glob_valor[ifile]->GetBinContent(hchisq_glob_valor[ifile]->GetMinimumBin());
    if(chisq_min[ifile] < chisq_min_glob) chisq_min_glob = chisq_min[ifile];

    cout << "File " << inname[ifile] << " :" << endl;
    cout << "Min chisq = " << chisq_min[ifile] << endl;
    cout << endl;
  }
  cout << "Global min chisq = " << chisq_min_glob << endl;
  cout << endl;

  
  for(int ifile=0;ifile<Nfiles;ifile++){
    
    hdeltachisq_glob_valor_68[ifile] = (TH2D*) hchisq_glob_valor[ifile]->Clone(Form("hdeltachisq_glob_68_%d",ifile));
    hdeltachisq_glob_valor_90[ifile] = (TH2D*) hchisq_glob_valor[ifile]->Clone(Form("hdeltachisq_glob_90_%d",ifile));
    
    for(int i=0;i<hchisq_glob_valor[ifile]->GetNbinsX();i++){
      for(int j=0;j<hchisq_glob_valor[ifile]->GetNbinsY();j++){
	
	double cont = (double) hchisq_glob_valor[ifile]->GetBinContent(i+1,j+1);
	hdeltachisq_glob_valor_68[ifile]->SetBinContent(i+1,j+1,cont-chisq_min_glob);
	hdeltachisq_glob_valor_90[ifile]->SetBinContent(i+1,j+1,cont-chisq_min_glob);
	//hdeltachisq_glob_valor_68[ifile]->SetBinContent(i+1,j+1,cont-chisq_min[ifile]);
	//hdeltachisq_glob_valor_90[ifile]->SetBinContent(i+1,j+1,cont-chisq_min[ifile]);
	//cout << "cont = " << cont << endl;
      }
    }
  
    hdeltachisq_glob_valor[ifile] = (TH2D*) hdeltachisq_glob_valor_68[ifile]->Clone(Form("hdeltachisq_glob_%d",ifile));

    hdeltachisq_glob_valor_68[ifile]->SetContour(1);
    //hdeltachisq_glob_valor_68[ifile]->SetContourLevel(0,chisq_min+2.3);
    hdeltachisq_glob_valor_68[ifile]->SetContourLevel(0,2.3);
    hdeltachisq_glob_valor_68[ifile]->SetLineWidth(2);
    hdeltachisq_glob_valor_68[ifile]->SetLineColor(ifile+1);
    hdeltachisq_glob_valor_68[ifile]->SetLineStyle(2);
    
    hdeltachisq_glob_valor_90[ifile]->SetContour(1);
    //hdeltachisq_glob_valor_90[ifile]->SetContourLevel(0,chisq_min+4.61);
    hdeltachisq_glob_valor_90[ifile]->SetContourLevel(0,4.61);
    hdeltachisq_glob_valor_90[ifile]->SetLineWidth(2);
    hdeltachisq_glob_valor_90[ifile]->SetLineColor(ifile+1);
    hdeltachisq_glob_valor_90[ifile]->SetLineStyle(1);
  }
 
  
  // Draw Likelihood
  TCanvas *cLike[Nfiles];
  for(int ifile=0;ifile<Nfiles;ifile++){
    cLike[ifile] = new TCanvas(Form("cLike_file_%d",ifile),Form("cLike_file_%d",ifile));
    cLike[ifile]->cd();
    //hdeltachisq_glob_valor_90[ifile]->Draw("colz");
    hchisq_glob_valor[ifile]->Draw("colz");
    hdeltachisq_glob_valor[ifile]->Draw("colz");
  }
  

  // Contours independent hierarchy

  TCanvas *ccontours = new TCanvas("ccontours","ccontours_independent");
  hInit->Draw();

  double posX = hInit->GetXaxis()->GetBinCenter(2);
  double posY = hInit->GetYaxis()->GetBinCenter(2);

  TLatex *latPreliminary = new TLatex(posX,posY,str_preliminary);
  latPreliminary->SetTextSize(0.05);
  if(draw_preliminary) latPreliminary->Draw("same");



  // Draw reactors box
  
  //double X1 = kGlobFitSinsqTheta13 - kGlobFitSinsqTheta13*kGlobFitErrSinsqTheta13;
  //double X2 = kGlobFitSinsqTheta13 + kGlobFitSinsqTheta13*kGlobFitErrSinsqTheta13;
  //cout << kGlobFitSinsqTheta13 << " +- " << kGlobFitSinsqTheta13*kGlobFitErrSinsqTheta13 << endl;

  const double sinsqth13_nom = 0.0219;
  const double sinsqth13_err = 0.0012;
  double X1 = sinsqth13_nom + sinsqth13_err;
  double X2 = sinsqth13_nom - sinsqth13_err;
  
  double Y1 = hchisq_valor[0]->GetYaxis()->GetBinLowEdge(hchisq_valor[0]->GetYaxis()->GetFirst());
  //double Y2 = hInit->GetMaximum();
  double Y2 = hchisq_valor[0]->GetYaxis()->GetBinUpEdge(hchisq_valor[0]->GetYaxis()->GetLast());
 
  TBox *breactors = new TBox(X1,Y1,X2,Y2);
  breactors->SetFillColor(kYellow+1);
  breactors->SetLineColor(kYellow+1);
  if(draw_reactor_band) breactors->Draw("same");
  
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int icont=0;icont<Ncont_90;icont++){
      if(g90_valor[icont][ifile]){
	if(scalegraphX) ScaleGraphAxis(g90_valor[icont][ifile],1000.,1.);
	if(scalegraphY) ScaleGraphAxis(g90_valor[icont][ifile],1.,1000.);
        g90_valor[icont][ifile]->SetLineWidth(2);
        g90_valor[icont][ifile]->SetLineColor(ifile+1);
        g90_valor[icont][ifile]->SetLineStyle(1);
        g90_valor[icont][ifile]->DrawClone("l same");
      }
    }
    for(int icont=0;icont<Ncont_68;icont++){
      if(g68_valor[icont][ifile]){
	if(scalegraphX) ScaleGraphAxis(g68_valor[icont][ifile],1000.,1.);
	if(scalegraphY) ScaleGraphAxis(g68_valor[icont][ifile],1.,1000.);
        g68_valor[icont][ifile]->SetLineWidth(2);
        g68_valor[icont][ifile]->SetLineColor(ifile+1);
        g68_valor[icont][ifile]->SetLineStyle(2);
        g68_valor[icont][ifile]->DrawClone("l same");
      }
    }
    if(scalegraphX) ScaleMarkerAxis(gbest_fit_global[ifile],1000.,1.);
    if(scalegraphY) ScaleMarkerAxis(gbest_fit_global[ifile],1.,1000.);
    gbest_fit_global[ifile]->Draw("same");
    gbest_fit_global[ifile]->SetMarkerColor(ifile+1);
  }



  TLegend *lcontours2 = new TLegend(0.5,0.5,0.8,0.65);
  lcontours2->SetBorderSize(0);
  lcontours2->SetFillColor(0);
  lcontours2->SetLineColor(0);
  lcontours2->SetFillColor(kWhite);
  lcontours2->SetFillStyle(1);
  lcontours2->SetTextFont(132);
  lcontours2->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    lcontours2->AddEntry(g90_valor[0][ifile],latexname[ifile],"l");
  }
  lcontours2->Draw();

    
  TLegend *lcontours = new TLegend(0.5,0.45,0.8,0.65);
  lcontours->SetBorderSize(0);
  lcontours->SetFillColor(0);
  lcontours->SetLineColor(0);
  lcontours->SetFillColor(kWhite);
  lcontours->SetFillStyle(1);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.05);
  lcontours->AddEntry(g68_valor[0][0],"68.27%CL","l");
  lcontours->AddEntry(g90_valor[0][0],"90%CL","l");  
  //lcontours->AddEntry(g68_valor[0][0],"68.27%CL (-2ln L_{crit} = 2.3)","l");
  //lcontours->AddEntry(g90_valor[0][0],"90%CL (-2ln L_{crit} = 4.61)","l");  
  lcontours->AddEntry(gbest_fit_global[0],"Best-fit","p");
  if(draw_reactor_band) lcontours->AddEntry(breactors,"PDG 2015","f");
  lcontours->Draw();

  double posX = hInit->GetXaxis()->GetBinCenter(3);
  double posY = hInit->GetYaxis()->GetBinCenter(3);
  TLatex *latcontours = new TLatex(posX,posY,"Fixed Mass Hierarchy");
  latcontours->SetTextSize(0.06);
  latcontours->Draw("same");
  
  
  return;






  // // TCanvas *cDeltaChisq[Nfiles];
  // // for(int ifile=0;ifile<Nfiles;ifile++){
  // //   cDeltaChisq[ifile] = new TCanvas(Form("cDeltaChisq_%d",ifile),Form("cDeltaChisq_%d",ifile));
  // //   if(ifile==0) hdeltachisq_glob_valor_68[ifile]->Draw("colz");
  // //   else hdeltachisq_glob_valor_68[ifile]->Draw("colz same");
  // //   hdeltachisq_glob_valor_90[ifile]->Draw("colz same");
  // // }

  // TCanvas *cContGlobMin = new TCanvas("cContGlobMin","cContGlobMin_global");
  // hInit->Draw();
  // for(int ifile=0;ifile<Nfiles;ifile++){
  //   hdeltachisq_glob_valor_68[ifile]->Draw("cont3 same");
  //   hdeltachisq_glob_valor_90[ifile]->Draw("cont3 same");
  //   gbest_fit_global[ifile]->Draw("same");
  //   gbest_fit_global[ifile]->SetMarkerColor(ifile+1);
  // }
  // if(draw_reactor_band) breactors->Draw("");
  // lcontours->Draw();
  

  // // Print
  // cout << endl;
  // cout << "The global best-fit is" << endl;
  // for(int ifile=0;ifile<Nfiles;ifile++){
  //   cout << " - " << latexname[ifile] << endl;
  //   cout << "   bf_x = " << gbest_fit_global[ifile]->GetX() << endl;
  //   cout << "   bf_y = " << gbest_fit_global[ifile]->GetY() << endl;
  //   cout << "   chisq = " << chisq_min[ifile] << endl;
  // }
  // cout << endl;

  return;
}



void ConvertHistoFromSinsqthToSinsq2th(TGraph *graph){
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    double sinsq2th = util::FromSinsqthToSinsq2th(x);
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,sinsq2th,y);
  }
  return;
}

void SwapXY(TGraph *graph){
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,y,x);
  }
  return;
}

TH2D* SwapXY(TH2D *graph){
   
  int nbinsX = graph->GetNbinsX();
  int nbinsY = graph->GetNbinsY();
  double xmin = graph->GetXaxis()->GetXmin();
  double xmax = graph->GetXaxis()->GetXmax();
  double ymin = graph->GetYaxis()->GetXmin();
  double ymax = graph->GetYaxis()->GetXmax();
 
  TH2D* graph_out = new TH2D("h_temp","",nbinsY,ymin,ymax,nbinsX,xmin,xmax);
 
  //TArrayD *arrayx = graph->GetXaxis()->GetXbins();
  //TArrayD *arrayy = graph->GetYaxis()->GetXbins();
  //TH2D* graph_out = new TH2D("h_temp","",nbinsY,arrayy->GetArray(),nbinsX,arrayx->GetArray());

  graph_out->GetXaxis()->SetTitle(graph->GetYaxis()->GetTitle());
  graph_out->GetYaxis()->SetTitle(graph->GetXaxis()->GetTitle());
  
  for(int i=0;i<graph_out->GetNbinsX();i++){
    for(int j=0;j<graph_out->GetNbinsY();j++){  
      double bincont = graph->GetBinContent(j+1,i+1);
      graph_out->SetBinContent(i+1,j+1,bincont);
    }
  }  
  return graph_out;

}

void SwapXY(TMarker *graph){
  double x = graph->GetX();
  double y = graph->GetY();
  graph->SetX(y);
  graph->SetY(x);
  
  return;
}

void ScaleMarkerAxis(TMarker *marker,double scale_x,double scale_y){
  double x = marker->GetX();
  double y = marker->GetY();
  double x_new = x * scale_x;
  double y_new = y * scale_y;
  marker->SetX(x_new);
  marker->SetY(y_new);
  return;
}
void ScaleGraphAxis(TGraph *graph,double scale_x,double scale_y){
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    double x_new = x * scale_x;
    double y_new = y * scale_y;
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,x_new,y_new);
  }
  return;
}
