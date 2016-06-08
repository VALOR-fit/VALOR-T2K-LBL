
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;

void plot_contours_dcp_many
( 
 const bool SaveRootFile = true,
 const string tag = "dcp_run17_react_Asimov_B",

 const bool doGlobalMin  = true,
 const bool same_mh_color = true, // use same color for NH and IH
 const bool same_mh_style = false, // use same style for NH and IH
 const bool draw_preliminary = false,

 const bool do_fc = false,
 const bool draw_fc_err = false,

 const bool do_const = true,

 const bool draw_1sigmaCL = false,
 const bool draw_90CL     = true,
 const bool draw_2sigmaCL = true,
 const bool draw_3sigmaCL = false,
 const bool draw_4sigmaCL = false,
 const bool draw_5sigmaCL = false,

 const bool draw_reactors = false,

 const bool do_range_1sigma = true,
 const double critval = 1.,
 
 // Compare true chisq
 
 // // run 1-6 best-fit (global)
 // const double truepar0 = -1.8535397,
 // const string truemh0 = "normal",

 // Asimov A best-fit (global)
 const double truepar0 = -1.601,
 const string truemh0 = "normal",

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

  // NIWG2015c Summer 2016 statistics
  const int Nfiles = 1;
  const int Nmh = 2;
  char *infilename[Nfiles][Nmh] = {
    

    //////////////////////////////////////////////
    //                                          //
    //      Run 1-7 -Asimov - Run 1-7 Flux      // 
    //                                          //
    //////////////////////////////////////////////

    // // deltaCP - Run 1-7 Asimov A T2K only - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/files/asimov_fit_nh_a_dcp_t2konly_run17.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/files/asimov_fit_ih_a_dcp_t2konly_run17.root"
    // },

    // // deltaCP - Run 1-7 Asimov B T2K only - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/files/asimov_fit_nh_b_dcp_t2konly_run17.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/files/asimov_fit_ih_b_dcp_t2konly_run17.root"
    // },
    
    // // deltaCP - Run 1-7 Asimov A w/ reactors - Flux Run 1-7
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/files/fit_fitmarg_joint_dcp_react_nh_Asimov_A_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/files/fit_fitmarg_joint_dcp_react_ih_Asimov_A_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },

    // deltaCP - Run 1-7 Asimov B w/ reactors - Flux Run 1-7
    {
      "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/files/fit_fitmarg_joint_dcp_react_nh_Asimov_B_Run17_nexp_1_0.000_-999.000-job0-_0_0.root",
      "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/files/fit_fitmarg_joint_dcp_react_ih_Asimov_B_Run17_nexp_1_0.000_-999.000-job0-_0_0.root"
    },
  
    // // deltaCP - Run 1-7 Asimov A w/ reactors - Flux Run 1-6
    // {
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun16/files/fit_fitmarg_joint_dcp_react_nh_Asimov_A_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun16/files/fit_fitmarg_joint_dcp_react_ih_Asimov_A_Run17_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // }
    
    
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

    
    //////////////////////////////
    //                          //
    //      Run 1-6 - Data      // 
    //                          //
    //////////////////////////////

    // // deltaCP - Run 1-6 Data w/ reactors
    // {
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_nh_Data_Run16_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH 
    //   "../TECH_NOTE_14_3_16/DATA/files/fit_fitmarg_joint_dcp_react_ih_Data_Run16_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
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
  char *fcname = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/FCgraphs_p8k.root";
  //char *fcname = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/FCgraphs_4k.root";


  char *latexname[Nfiles][Nmh] = {
    {"Normal Hierarchy","Inverted Hierarchy"},
        
    //{"Run 1-7 - NH", "Run 1-7 - IH"},
    //{"Run 1-6 - NH", "Run 1-6 - IH"},
    
    // {"14.4e20 FHC + 7.47e20 RHC - NH",   "14.4e20 FHC + 7.47e20 RHC - IH"},
    // {"10.07e20 FHC + 11.17e20 RHC - NH", "10.07e20 FHC + 11.17e20 RHC - IH"},
    // {"7.00e20 FHC + 14.87e20 RHC - NH",  "7.00e20 FHC + 14.87e20 RHC - IH"}
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
    
    double ChisqMin_0 = hchisq1D[ifile][0]->GetMinimum();

    double ChisqMin_1 = 999999;
    if(hchisq1D[ifile][1]) ChisqMin_1 = hchisq1D[ifile][1]->GetMinimum();
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
      double ChisqMin_Curr = hchisq1D[ifile][imh]->GetMinimum();
      for(int i=0;i<hdeltachisq1D_GlobMin[ifile][imh]->GetNbinsX();i++){
	double cont = hchisq1D[ifile][imh]->GetBinContent(i+1); 
	if(doGlobalMin) hdeltachisq1D_GlobMin[ifile][imh]->SetBinContent(i+1,cont-ChisqMin_Glob);
	else            hdeltachisq1D_GlobMin[ifile][imh]->SetBinContent(i+1,cont-ChisqMin_Curr);
      }
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

    gFC90CL_nh   ->SetLineColor(kBlue);
    gFC90CL_nh_p1->SetLineColor(kBlue);
    gFC90CL_nh_m1->SetLineColor(kBlue);
    gFC90CL_ih   ->SetLineColor(kBlue);
    gFC90CL_ih_p1->SetLineColor(kBlue);
    gFC90CL_ih_m1->SetLineColor(kBlue);
    gFC2sigma_nh   ->SetLineColor(kGreen);
    gFC2sigma_nh_p1->SetLineColor(kGreen);
    gFC2sigma_nh_m1->SetLineColor(kGreen);
    gFC2sigma_ih   ->SetLineColor(kGreen);
    gFC2sigma_ih_p1->SetLineColor(kGreen);
    gFC2sigma_ih_m1->SetLineColor(kGreen);

    gFC90CL_nh   ->SetMarkerColor(kBlue);
    gFC90CL_nh_p1->SetMarkerColor(kBlue);
    gFC90CL_nh_m1->SetMarkerColor(kBlue);
    gFC90CL_ih   ->SetMarkerColor(kBlue);
    gFC90CL_ih_p1->SetMarkerColor(kBlue);
    gFC90CL_ih_m1->SetMarkerColor(kBlue);
    gFC2sigma_nh   ->SetMarkerColor(kGreen);
    gFC2sigma_nh_p1->SetMarkerColor(kGreen);
    gFC2sigma_nh_m1->SetMarkerColor(kGreen);
    gFC2sigma_ih   ->SetMarkerColor(kGreen);
    gFC2sigma_ih_p1->SetMarkerColor(kGreen);
    gFC2sigma_ih_m1->SetMarkerColor(kGreen);

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
    
    gFC90CL_nh   ->SetLineWidth(2);
    gFC90CL_nh_p1->SetLineWidth(2);
    gFC90CL_nh_m1->SetLineWidth(2);
    gFC90CL_ih    ->SetLineWidth(2);
    gFC90CL_ih_p1->SetLineWidth(2);
    gFC90CL_ih_m1->SetLineWidth(2);
    gFC2sigma_nh   ->SetLineWidth(2);
    gFC2sigma_nh_p1->SetLineWidth(2);
    gFC2sigma_nh_m1->SetLineWidth(2);
    gFC2sigma_ih   ->SetLineWidth(2);
    gFC2sigma_ih_p1->SetLineWidth(2);
    gFC2sigma_ih_m1->SetLineWidth(2);

    gFC90CL_ih    ->SetLineStyle(2);
    gFC90CL_ih_p1->SetLineStyle(2);
    gFC90CL_ih_m1->SetLineStyle(2);
    gFC2sigma_ih   ->SetLineStyle(2);
    gFC2sigma_ih_p1->SetLineStyle(2);
    gFC2sigma_ih_m1->SetLineStyle(2);
  }



  // Draw histograms
  
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
      if(hdeltachisq1D[ifile][imh]->GetMaximum()>Max_deltachisq) 
	Max_deltachisq = hdeltachisq1D[ifile][imh]->GetMaximum();
    }
  }

  sprintf(name,"cDeltaChisq1D");
  cDeltaChisq = new TCanvas(name,name);

  hInit->SetMaximum(Max_deltachisq*1.2);
  hInit->DrawCopy();

  // Draw reactors box
  const double sinsqth13_nom = 0.0219;
  const double sinsqth13_err = 0.0012;
  double X1 = sinsqth13_nom - sinsqth13_err;
  double X2 = sinsqth13_nom + sinsqth13_err;
  double Y1 = hInit->GetMinimum();
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
      if((ifile+imh)==0) hdeltachisq1D_GlobMin[ifile][imh]->Draw("same");
      else hdeltachisq1D_GlobMin[ifile][imh]->Draw("same");
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
      lcritval->AddEntry(gFC90CL_nh,"#Delta #chi^{2}_{crit} (90%CL) - NH","l");
      lcritval->AddEntry(gFC90CL_ih,"#Delta #chi^{2}_{crit} (90%CL) - IH","l");
    }
    if(draw_2sigmaCL){
      lcritval->AddEntry(gFC2sigma_nh,"#Delta #chi^{2}_{crit} (2#sigma) - NH","l");
      lcritval->AddEntry(gFC2sigma_ih,"#Delta #chi^{2}_{crit} (2#sigma) - IH","l");
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

  TLatex *latPreliminary = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),1.1,"PRELIMINARY");
  if(draw_preliminary) latPreliminary->Draw("same");

  // Store in root file
  if(SaveRootFile){
    TFile* f2 = new TFile(Form("deltachisq_%s.root",tag.c_str()),"recreate");    
    for(int ifile=0;ifile<Nfiles;ifile++){
      for(int imh=0;imh<Nmh;imh++){
  	sprintf(name,"hdeltachisq_File_%i_MH_%i",ifile,imh);
  	hdeltachisq1D_GlobMin[ifile][imh]->Write();//(name);
      }
    }    
    f2->Close();
  }
  
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
    chisq_min = histo->GetMinimum();
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
    if(bincont<0.){
      cerr << endl;
      cerr << "CalcBestFitErrorsMINOS" << endl;
      cerr << "Bin content is negative!" << endl;
      cerr << "bincont = " << bincont << endl;
      cerr << endl;
    }
    
    // cout << "ibin = " << ibin 
    // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
    // 	   << " - bincont = " << bincont << endl;
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1 (following x direction)
      bin_negsigma_b = ibin; 
      bin_negsigma_a = bin_negsigma_b+1; 
      break;
    }
    else if(bincont<0.){
      cerr << endl;
      cerr << "CalcBestFitErrorsMINOS" << endl;
      cerr << "Negative DeltaChisq 1D:" << endl;
      cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
      cerr << endl;
      //exit(1);
      return conv::kBadResult;
    }
        
  }
  
   
  // Positive error
  
  // cout << endl;
  // cout << "Positive error" << endl;
  // cout << endl;
  
  for(int ibin=bin_bf;ibin<=histo->GetNbinsX();ibin++){ // histogram convention 
    double bincont = histo->GetBinContent(ibin);
    if(bincont<0.){
      cerr << endl;
      cerr << "CalcBestFitErrorsMINOS" << endl;
      cerr << "Negative DeltaChisq 1D:" << endl;
      cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
      cerr << endl;
      //exit(1);
      return conv::kBadResult;
    }
    
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
  //TGraph *gr = new TGraph();
  //gr->SetPoint(0,X1,Y1);
  //gr->SetPoint(1,X2,Y2);
  //cout << "Eval = " << gr->Eval(Xin) << endl;
 }
