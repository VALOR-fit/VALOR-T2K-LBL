//
// Sensitivity study w/ 1D or 2D plost
//


void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style);

void draw_sensitivity_plots_toys_NEW(
				 const bool PrintPlots = false,
				 const char * output_file = "toys_sensitivity_plots.root",

				 const bool NuBarMode = true,

				 const bool doShiftMin = false,  // Don't shift wrt global minimum if Asimov dataset is used!!!

				 const double CLval = 4.61, //2.71, // 2.71, //4.61, // 90%CL
				 const bool reactor = false, 
				 const bool theta2theta = false,
				 const char * parname1_def = "mixing23", //"mixing23",
				 const char * parname2_def = "dmsq", //"dmsq",
				 const char * parname3_def = "none",
				 
				 const bool ChangeOrder_nh = false, // change order of params (2D plots)
				 const bool ChangeOrder_ih = false, // change order of params (2D plots) 
				 
				 const bool DrawLeg = true,
				 const char * LegendName_ih = "IH", 
				 const char * LegendName_nh = "NH",
				 
				 //const char *tag = "t13_Asimov_0p170e21_NuOscBar",
				 const char *tag = "t13_fakedataset_F_0p170e21_NuOscBar",
				 
				 //const char *tag = "t23Vsdmsq_Asimov_0p170e21_NuOscBar",
				 //const char *tag = "t23Vsdmsq_fakedataset_F_0p170e21_NuOscBar",



				 // th13 - Asimov - NuOscBar Fix NuOsc - Stat only

				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fitgrid_NuOscBar_0p170e21_Asimov_NH_mixing13_none_none_BinsTot_101_BinsTot_1_BinsTot_1_BinsPart_1_101_BinsPart_1_1_BinsPart_1_1",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fitgrid_NuOscBar_0p170e21_Asimov_IH_mixing13_none_none_BinsTot_101_BinsTot_1_BinsTot_1_BinsPart_1_101_BinsPart_1_1_BinsPart_1_1",				 

				 // th13 - fakedataset - NuOscBar Fix NuOsc - Stat only
				 
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fakedataset_0p170e21_201411/fitgrid_NuOscBar_0p170e21_fakedata_D_NH_mixing13_none_none_BinsTot_101_BinsTot_1_BinsTot_1_BinsPart_1_101_BinsPart_1_1_BinsPart_1_1",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fakedataset_0p170e21_201411/fitgrid_NuOscBar_0p170e21_fakedata_D_IH_mixing13_none_none_BinsTot_101_BinsTot_1_BinsTot_1_BinsPart_1_101_BinsPart_1_1_BinsPart_1_1",

				 // th23 Vs dmsq - Asimov - NuOscBar Fix NuOsc - Stat only

				 const char * name_nh = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fitgrid_NuOscBar_0p170e21_Asimov_NH_mixing23_dmsq_none_BinsTot_51_BinsTot_51_BinsTot_1_BinsPart_1_51_BinsPart_1_51_BinsPart_1_1",
				 const char * name_ih = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fitgrid_NuOscBar_0p170e21_Asimov_IH_mixing23_dmsq_none_BinsTot_51_BinsTot_51_BinsTot_1_BinsPart_1_51_BinsPart_1_51_BinsPart_1_1",

				 // th23 Vs dmsq - fakedataset - NuOscBar Fix NuOsc - Stat only
				 
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fakedataset_0p170e21_201411/fitgrid_NuOscBar_0p170e21_fakedata_D_NH_mixing23_dmsq_none_BinsTot_51_BinsTot_51_BinsTot_1_BinsPart_1_51_BinsPart_1_51_BinsPart_1_1",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fakedataset_0p170e21_201411/fitgrid_NuOscBar_0p170e21_fakedata_D_IH_mixing23_dmsq_none_BinsTot_51_BinsTot_51_BinsTot_1_BinsPart_1_51_BinsPart_1_51_BinsPart_1_1",

				 // const char * name_nh = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fakedataset_0p170e21_201411/fitgrid_NuOscBar_0p170e21_Asimov_NH_mixing23_dmsq_none_BinsTot_51_BinsTot_51_BinsTot_1_BinsPart_1_51_BinsPart_1_51_BinsPart_1_1_LargeDm2",
				 // const char * name_ih = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fakedataset_0p170e21_201411/fitgrid_NuOscBar_0p170e21_Asimov_NH_mixing23_dmsq_none_BinsTot_51_BinsTot_51_BinsTot_1_BinsPart_1_51_BinsPart_1_51_BinsPart_1_1_LargeDm2",






				 // Fit systematics
				 
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/sensi-Sensi_FHCpRHC_DataRun1to4_MCRun5to6_th13DB_NH_500_all",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/sensi-Sensi_FHCpRHC_DataRun1to4_MCRun5to6_th13DB_IH_500_all",

				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/sensi-Sensi_FHCpRHC_DataRun1to4_MCRun5to6_ResultsRun1to4_NH_500",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/sensi-Sensi_FHCpRHC_DataRun1to4_MCRun5to6_ResultsRun1to4_IH_500",


				 // 1 ring studies

				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/SystEffVsStat/fitgrid_Run1to4_1ring13a_Asimov_NH_dcp_Stat_deltacp_none_none",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/SystEffVsStat/fitgrid_Run1to4_1ring13a_Asimov_NH_dcp_StatSyst_deltacp_none_none",

				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/SystEffVsStat/fitgrid_Run1to4_1ring13a_Asimov_NH_th23Vsdmsq_Stat_mixing23_dmsq_none",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/SystEffVsStat/fitgrid_Run1to4_1ring13a_Asimov_NH_th23Vsdmsq_StatSyst_mixing23_dmsq_none",

				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/1ring13a_15_9_14/fitgrid_Run1to4_1ring13a_Asimov_NH_mixing230p4_dmsq_7p8x21_mixing23_dmsq_none",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/1ring13a_15_9_14/fitgrid_Run1to4_1ring13a_Asimov_IH_mixing230p4_dmsq_7p8x21_mixing23_dmsq_none",


				 // MultiRing studies w/ Asimov data set - deltacp only

				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_nr_Asimov_NH_Run14_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_nr_Asimov_NH_Run14_erecdelta_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_Asimov_NH_Run14_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_ElikeMue_Asimov_NH_Run14_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2rpi0_2rMulike_Asimov_NH_Run14_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_Asimov_NH_FullT2K_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_ElikeFitQun_2rMulike_Asimov_NH_FullT2K_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_ElikeMue_2rMulike_Asimov_NH_FullT2K_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_ElikeFitQun_2rMulike_Asimov_NH_FullT2K_deltacp_none_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/7_10_14/fitgrid_1r_2r_ElikeMue_fiTQun_Asimov_NH_FullT2K_deltacp_none_none", 
				 
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1rOnly_Asimov_NH_Run14_deltacp_none_none",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1rOnly_Asimov_NH_FullT2K_deltacp_none_none",
				 
				 
				 // MultiRing studies w/ Asimov data set - sinsqth23 Vs dmsq32
				 
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_nr_Asimov_NH_Run14_sin2th23Vsdmsq_mixing23_dmsq_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_nr_Asimov_NH_Run14_erecdelta_mixing23_dmsq_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_Asimov_NH_Run14_mixing23_dmsq_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_ElikeMue_Asimov_NH_Run14_mixing23_dmsq_none",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1r_2r_Asimov_NH_FullT2K_mixing23_dmsq_none",
				 
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1rOnly_Asimov_NH_Run14_mixing23_dmsq_none",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/MultiRing/2_10_14/fitgrid_1rOnly_Asimov_NH_FullT2K_mixing23_dmsq_none",


		 


				 // MultiRing studies w/ Asimov data set - sinsqth23 Vs dmsq32 (numu disapp only)
				 // OFF-AXIS OK
				 // //const char * name_nh = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p0_NH_dcp_51_deltacp_none_none_Asimov_FluxOK",
				 // const char * name_nh = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p25_NH_dcp_51_deltacp_none_none_Asimov_FluxOK",
				 // const char * name_ih = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p5_NH_dcp_51_deltacp_none_none_Asimov_FluxOK",
				 // const char * name_nh = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p0_NH_dcp_51_deltacp_none_none_Asimov_FluxOK_Steve",
				 // //const char * name_nh = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p25_NH_dcp_51_deltacp_none_none_Asimov_FluxOK_Steve",
				 // const char * name_ih = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p5_NH_dcp_51_deltacp_none_none_Asimov_FluxOK_Steve",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p0_IH_dcp_51_deltacp_none_none_Asimov_FluxOK",
				 // //const char * name_nh = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p25_IH_dcp_51_deltacp_none_none_Asimov_FluxOK",
				 //const char * name_ih = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p5_IH_dcp_51_deltacp_none_none_Asimov_FluxOK",
				 //const char * name_nh = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p0_IH_dcp_51_deltacp_none_none_Asimov_FluxOK_Steve",
				 // const char * name_nh = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p25_IH_dcp_51_deltacp_none_none_Asimov_FluxOK_Steve",
				 // const char * name_ih = "$T2KOSCANA/outputs/sensi/OffAxis/NEWTEST_OFFFLUX/fitgrid_Run1to6_joint_OA2p5_IH_dcp_51_deltacp_none_none_Asimov_FluxOK_Steve",
				 )
{
  using namespace std;
  using namespace std::string;
  using namespace valor;
  using namespace valor::conv;

  const char * filename_nh = Form("%s.root",name_nh);
  FitParamList * fparm_nh = 0;
  PhysParamList * iparm_nh = new PhysParamList;
  int nbins_param1_nh, nbins_param2_nh;
  int nbins_tot_param1_nh, nbins_tot_param2_nh;
  int nbins1_nh = 0;
  int nbins2_nh = 0;
  string *parname1_nh=new string();
  string *parname2_nh=new string();

  const char * filename_ih = Form("%s.root",name_ih);
  FitParamList * fparm_ih = 0;
  PhysParamList * iparm_ih = new PhysParamList;
  int nbins_param1_ih, nbins_param2_ih;
  int nbins_tot_param1_ih, nbins_tot_param2_ih;
  int nbins1_ih = 0;
  int nbins2_ih = 0;
  string *parname1_ih=new string();
  string *parname2_ih=new string();
  
  int max_bin_param1, min_bin_param1, max_bin_param2, min_bin_param2;
  double max_val_param1, min_val_param1, max_val_param2, min_val_param2;
  double step1, step2;

  // 1D: dChi2 vs param1  (Normal Hierarchy)

  cout << endl;
  cout << "Read file: " << filename_nh << endl;
  cout << endl;

  const char * filename_nh = Form("%s.root",name_nh);
  TChain ch_nh("fitnt");

  if(ChangeOrder_nh){
    ch_nh.Add(filename_nh);
    ch_nh.SetBranchAddress("best_fit_params", &fparm_nh); 
    ch_nh.SetBranchAddress("nbins_param1", &nbins_param2_nh); 
    ch_nh.SetBranchAddress("nbins_param2", &nbins_param1_nh); 
    ch_nh.SetBranchAddress("nbins_tot_param1", &nbins_tot_param2_nh); 
    ch_nh.SetBranchAddress("nbins_tot_param2", &nbins_tot_param1_nh); 
    ch_nh.SetBranchAddress("best_fit_params", &fparm_nh); 
    ch_nh.SetBranchAddress("parname1", &parname2_nh); // name of the grid parameter 2
    ch_nh.SetBranchAddress("parname2", &parname1_nh); // name of the grid parameter 1
    ch_nh.SetBranchAddress("max_val_param1",&max_val_param2);
    ch_nh.SetBranchAddress("min_val_param1",&min_val_param2);
    ch_nh.SetBranchAddress("max_val_param2",&max_val_param1);
    ch_nh.SetBranchAddress("min_val_param2",&min_val_param1);
    ch_nh.SetBranchAddress("step1",&step2);
    ch_nh.SetBranchAddress("step2",&step1);
  }
  else{
    ch_nh.Add(filename_nh);
    ch_nh.SetBranchAddress("best_fit_params", &fparm_nh); 
    ch_nh.SetBranchAddress("nbins_param1", &nbins_param1_nh); 
    ch_nh.SetBranchAddress("nbins_param2", &nbins_param2_nh); 
    ch_nh.SetBranchAddress("nbins_tot_param1", &nbins_tot_param1_nh); 
    ch_nh.SetBranchAddress("nbins_tot_param2", &nbins_tot_param2_nh); 
    ch_nh.SetBranchAddress("best_fit_params", &fparm_nh); 
    ch_nh.SetBranchAddress("parname1", &parname1_nh); // name of the grid parameter 1
    ch_nh.SetBranchAddress("parname2", &parname2_nh); // name of the grid parameter 2
    ch_nh.SetBranchAddress("max_val_param1",&max_val_param1);
    ch_nh.SetBranchAddress("min_val_param1",&min_val_param1);
    ch_nh.SetBranchAddress("max_val_param2",&max_val_param2);
    ch_nh.SetBranchAddress("min_val_param2",&min_val_param2);
    ch_nh.SetBranchAddress("step1",&step1);
    ch_nh.SetBranchAddress("step2",&step2);
  }
  
  if(*parname1_nh==*parname2_nh && *parname1_nh!=""){
   cout << "Name parameter 1 == Name parameter 2 = " << *parname1_nh << endl; 
   exit(0);
  }


  const int nentries_nh = ch_nh.GetEntries();
  cout << "nentries (NH) = " << nentries_nh << endl;
  
  ch_nh.GetEntry(0);
  int nbins1_nh = TMath::Max(nbins_tot_param1_nh,nbins_param1_nh); // if using make_toygrid_toys_joint.C nbins_tot_param* is used instead of nbins_param*
  int nbins2_nh = TMath::Max(nbins_tot_param2_nh,nbins_param2_nh); // otherwise nbins_tot_param* = 1 < nbins_param*
  
  int nDim=0;
  if(nbins1_nh>1 || nbins2_nh>1) nDim = 1;
  if(nbins1_nh>1 && nbins2_nh>1) nDim = 2;
  
   if(*parname1_nh!=""){
    cout << endl;
    cout << "- Parameter 1: " << *parname1_nh << endl;
    cout << "  # bins = " << nbins1_nh << endl;
    cout << "  min = " << min_val_param1 << endl;
    cout << "  max = " << max_val_param1 << endl;
    cout << "  step = " << step1 << endl;
    cout << endl;
  } 
  if(*parname2_nh!=""){
  cout << "- Parameter 2: " << *parname2_nh << endl;
  cout << "  # bins = " << nbins2_nh << endl;
  cout << "  min = " << min_val_param2 << endl;
  cout << "  max = " << max_val_param2 << endl;
  cout << "  step = " << step2 << endl;
  cout << endl;
  }

  if(*parname1_nh=="" && *parname2_nh==""){
    cout << "# bins parameter 1 = " << nbins1_nh << endl; 
    cout << "# bins parameter 2 = " << nbins2_nh << endl; 
    cout << endl;
  }

  const int ntoys_nh = nentries_nh / (nbins1_nh * nbins2_nh);
  cout << "ntoys (NH) = " << ntoys_nh << endl;
  cout << endl;
  
  //Set up bins and histograms
  double mixing13_min = 0;
  double mixing13_max = 0;
  double mixing13_step = 0;
  double mixing23_min = 0;
  double mixing23_max = 0;
  double mixing23_step = 0;
  double dmsq_min = 0;
  double dmsq_max = 0;
  double dmsq_step = 0;
  double deltacp_min = 0;
  double deltacp_max = 0;
  double deltacp_step = 0;

  if(!theta2theta){
    if(reactor){
      mixing13_min = 0.015;
      mixing13_max = 0.035;
      //mixing13_step = 0.001;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_nh-1);
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.1;
      //mixing13_step = 0.005;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_nh-1);
    }
  }
  else{
    if(reactor){
      mixing13_min = 0.06;
      mixing13_max = 0.14;
      //mixing13_step = 0.004;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_nh-1);
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.4;
      //mixing13_step = 0.02;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_nh-1);
    }
  }
  mixing23_min = 0.;
  mixing23_max = 1.;
  //mixing23_min = 0.3;
  //mixing23_max = 0.7;
  // //mixing23_step = 0.01;
  mixing23_step = (mixing23_max-mixing23_min)/(nbins1_nh-1);
  



  //dmsq_min = 2.0;//E-3;
  //dmsq_max = 3.0;//E-3;
  // //dmsq_step = 0.05;//E-3;

  dmsq_min = 3.0;//E-4;
  dmsq_max = 1.0;//E-4;
  //dmsq_step = 0.05;//E-2;
  dmsq_step = (dmsq_max-dmsq_min)/(nbins1_nh-1);




  deltacp_min = -conv::kpi;
  deltacp_max = conv::kpi;
  //deltacp_step = 2*conv::kpi/50;
  deltacp_step = (deltacp_max-deltacp_min)/(nbins1_nh-1);
  

  // Take steps and min-max values from the tree if they exist

  if(*parname1_nh==""){
    
    parname1_nh[0].append(parname1_def); // fill parname2 w/ default if not in the tree
    
    if(parname1_def=="deltacp"){
      min_val_param1 = deltacp_min;
      max_val_param1 = deltacp_max;
      step1 = deltacp_step;
    }
    else if(parname1_def=="mixing13"){
      min_val_param1 = mixing13_min;
      max_val_param1 = mixing13_max;
      step1 = mixing13_step;
    }
    else if(parname1_def=="mixing23"){
      min_val_param1 = mixing23_min;
      max_val_param1 = mixing23_max;
      step1 = mixing23_step;
    }
    else if(parname1_def=="dmsq"){
      min_val_param1 = dmsq_min;
      max_val_param1 = dmsq_max;
      step1 = dmsq_step;
    }
    else if(parname1_def=="none"){
      min_val_param1 = 0.;
      max_val_param1 = 0.;
      step1 = 0.;
    }
  }
  
  if(*parname2_nh==""){

    parname2_nh[0].append(parname2_def); // fill parname2 w/ default if not in the tree

    if(parname2_def=="deltacp"){
      min_val_param2 = deltacp_min;
      max_val_param2 = deltacp_max;
      step2 = deltacp_step;
    }
    else if(parname2_def=="mixing13"){
      min_val_param2 = mixing13_min;
      max_val_param2 = mixing13_max;
      step2 = mixing13_step;
    }
    else if(parname2_def=="mixing23"){
      min_val_param2 = mixing23_min;
      max_val_param2 = mixing23_max;
      step2 = mixing23_step;
    }
    else if(parname2_def=="dmsq"){
      min_val_param2 = dmsq_min;
      max_val_param2 = dmsq_max;
      step2 = dmsq_step;
    }
    else if(parname2_def=="none"){
    }
  }
   
  //
  char histoname[200];

  // 1D histograms

  Double_t** bestfit_points_nh_x = new Double_t*[ntoys_nh]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    bestfit_points_nh_x[itoy] = new Double_t[nbins1_nh];
  }
  Double_t** bestfit_points_nh_y = new Double_t*[ntoys_nh]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    bestfit_points_nh_y[itoy] = new Double_t[nbins1_nh];
  }

  sprintf(histoname,"h_chisq_%s_Ave_nh",parname1_nh[0].c_str());
  TH1D * h_chisq_par1_Ave_nh = new TH1D(histoname,histoname,nbins1_nh,min_val_param1-0.5*step1,max_val_param1+0.5*step1);

  TH1D ** h_chisq_par1_nh = new TH1D*[ntoys_nh];
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    sprintf(histoname,"h_%s_nh_Toy%d",parname1_nh[0].c_str(),itoy);
    h_chisq_par1_nh[itoy] = new TH1D(histoname,histoname,nbins1_nh,min_val_param1-0.5*step1,max_val_param1+0.5*step1);
  }  
  
  TH1D ** h_chisq_par1_nh_bin = new TH1D*[nbins1_nh];
  for(int ibin=0;ibin<nbins1_nh;ibin++){
    char histoname[200];
    sprintf(histoname,"h_chisq_%s_nh_bin%d",parname1_nh[0].c_str(),ibin);
    h_chisq_par1_nh_bin[ibin] = new TH1D(histoname,histoname,50,0,500);    
  }

  // 2D histograms

  sprintf(histoname,"h_chisq_%s_%s_Ave_nh",parname2_nh[0].c_str(),parname1_nh[0].c_str());
  TH2D * h_chisq_par2_par1_Ave_nh = new TH2D(histoname,histoname,nbins1_nh,min_val_param1-0.5*step1,max_val_param1+0.5*step1,nbins2_nh,min_val_param2-0.5*step2,max_val_param2+0.5*step2);

  TH2D ** h_chisq_par2_par1_nh= new TH2D*[ntoys_nh];
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    char histoname[200];
    sprintf(histoname,"h_chisq_%s_%s_nh_Toy%d",parname2_nh[0].c_str(),parname1_nh[0].c_str(),itoy);
    h_chisq_par2_par1_nh[itoy] = new TH2D(histoname,histoname,nbins1_nh,min_val_param1-0.5*step1,max_val_param1+0.5*step1,nbins2_nh,min_val_param2-0.5*step2,max_val_param2+0.5*step2);
  }  

  TH1D *** h_chisq_par2_par1_nh_bin = new TH1D**[nbins1_nh];
  for(int ibin1=0;ibin1<nbins1_nh;ibin1++){
    h_chisq_par2_par1_nh_bin[ibin1] = new TH1D*[nbins2_nh];
    for(int ibin2=0;ibin2<nbins2_nh;ibin2++){
      char histoname[200];
      sprintf(histoname,"h_chisq_%s_%s_nh_bin%d_bin%d",parname2_nh[0].c_str(),parname1_nh[0].c_str(),ibin2,ibin1);
      h_chisq_par2_par1_nh_bin[ibin1][ibin2] = new TH1D(histoname,histoname,50,0,500);    
    }
  }

  // Loop over toys in the tree

  int iEvtToy = 0;
  for(int i = 0; i < nentries_nh; i++){
    iEvtToy = i / (nbins1_nh*nbins2_nh);
    if( !(i%(nbins1_nh*nbins2_nh)) ) cout << "# Toy " << iEvtToy << endl;

    ch_nh.GetEntry(i);
    int migradminstat_nh = fparm_nh->GetStatusCode("migrad_min");
    if (migradminstat_nh == 4) cout << "Fit failed!" << endl;
    //cout << "migradminstat_nh = " << migradminstat_nh << endl;

    double mix13_nh = fparm_nh->BestFitValue(conv::kPhypMixing13); // Parameters are fixed in the fit
    double mix23_nh = fparm_nh->BestFitValue(conv::kPhypMixing23);
    double dmsq_nh = fparm_nh->BestFitValue(conv::kPhypDm23sq);
    double deltacp_nh = fparm_nh->BestFitValue(conv::kPhypCPV);
    double chisq_nh = fparm_nh->GetBestFitStatistic();
    if(NuBarMode){
      mix13_nh = fparm_nh->BestFitValue(conv::kPhypMixing13Bar); // Parameters are fixed in the fit
      mix23_nh = fparm_nh->BestFitValue(conv::kPhypMixing23Bar);
      dmsq_nh = fparm_nh->BestFitValue(conv::kPhypDm23sqBar);
      deltacp_nh = fparm_nh->BestFitValue(conv::kPhypCPVBar);
      chisq_nh = fparm_nh->GetBestFitStatistic();
    }
 
    double bf1_nh=0;
    if(*parname1_nh=="deltacp") bf1_nh = deltacp_nh;
    else if(*parname1_nh=="mixing13") bf1_nh = mix13_nh;
    else if(*parname1_nh=="mixing23") bf1_nh = mix23_nh;
    else if(*parname1_nh=="dmsq") bf1_nh = dmsq_nh;
    double bf2_nh=0;
    if(*parname2_nh=="deltacp") bf2_nh = deltacp_nh;
    else if(*parname2_nh=="mixing13") bf2_nh = mix13_nh;
    else if(*parname2_nh=="mixing23") bf2_nh = mix23_nh;
    else if(*parname2_nh=="dmsq") bf2_nh = dmsq_nh;

    int bin_idx=0;
    double bincont=0.;

    if(nbins1_nh>1 && nbins2_nh<=1){ // 1D param1
      bin_idx = h_chisq_par1_nh[iEvtToy]->FindBin(bf1_nh);
      bincont = h_chisq_par1_nh[iEvtToy]->GetBinContent(bin_idx);
    }
    else if(nbins1_nh<=1 && nbins2_nh>1){ // 1D param2 
      bin_idx = h_chisq_par1_nh[iEvtToy]->FindBin(bf1_nh);
      bincont = h_chisq_par1_nh[iEvtToy]->GetBinContent(bin_idx);
    }
    else if(nbins1_nh>1 && nbins2_nh>1){ // 2D param2 Vs param1
      bin_idx = h_chisq_par2_par1_nh[iEvtToy]->FindBin(bf1_nh,bf2_nh);
      bincont = h_chisq_par2_par1_nh[iEvtToy]->GetBinContent(bin_idx);
    }
  
    if(bincont!=0){
      if(bincont>chisq_nh){
	if(chisq_nh!=0){
	  if(nDim==1) h_chisq_par1_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
	  else if(nDim==2) h_chisq_par2_par1_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
	  //bestfit_points_nh_x[iEvtToy][bin_idx-1]  = bf1_nh;//sinsq2theta13;
	  //bestfit_points_nh_y[iEvtToy][bin_idx-1]  = bf2_nh;
	}
      }
    } 
    else{ 
      if(chisq_nh!=0){
	if(nDim==1) h_chisq_par1_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
	else if(nDim==2) h_chisq_par2_par1_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
	//bestfit_points_nh_x[iEvtToy][bin_idx-1] = bf1_nh;//sinsq2theta13;
	//bestfit_points_nh_y[iEvtToy][bin_idx-1] = bf2_nh;
      }
    }
  } // end loop over toys


  // Get the toy distribution for each param value
  if(nDim==1){
    cout << endl << endl;
    for(int ibin=0;ibin<nbins1_nh;ibin++){
      for(int itoy=0;itoy<ntoys_nh;itoy++){
	double bincont = h_chisq_par1_nh[itoy]->GetBinContent(ibin+1);
	h_chisq_par1_nh_bin[ibin]->Fill(bincont);
	cout << bincont << " ";
      }
      cout << endl;
    }
    cout << endl << endl;
  }
  if(nDim==2){
    for(int ibin1=0;ibin1<nbins1_nh;ibin1++){
      for(int ibin2=0;ibin2<nbins2_nh;ibin2++){
  	for(int itoy=0;itoy<ntoys_nh;itoy++){
	  int ibin = h_chisq_par2_par1_nh[itoy]->GetBin(ibin1+1,ibin2+1);
  	  double bincont = h_chisq_par2_par1_nh[itoy]->GetBinContent(ibin);
	  h_chisq_par2_par1_nh_bin[ibin1][ibin2]->Fill(bincont);

	  cout << bincont << " ";
  	}
      }
      cout << endl;
    }
  }
  
  // Compute the average sensitivity    
  for(int itoy=0;itoy<ntoys_nh;itoy++){
    h_chisq_par1_Ave_nh->Add(h_chisq_par1_nh[itoy]);
    h_chisq_par2_par1_Ave_nh->Add(h_chisq_par2_par1_nh[itoy]);
  }
  h_chisq_par1_Ave_nh->Scale(1./ntoys_nh);
  h_chisq_par2_par1_Ave_nh->Scale(1./ntoys_nh);
  
  // Find the minimum chi2 (global scan)
  
  cout << endl;
  cout << "Find the minimum chi2 (global scan)" << endl;
  cout << endl;

  double *min_chisq_nh = new double[ntoys_nh];
  double min_chisq_Ave_nh=999.;
    
  if(nDim==1){ 
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      min_chisq_nh[itoy] = 999.;
      //cout << "# Toy " << itoy << endl;    
      for(int i = 1; i <= nbins1_nh; i++){
  	double bin = h_chisq_par1_nh[itoy]->GetBinContent(i);
  	if (bin < min_chisq_nh[itoy]) min_chisq_nh[itoy] = bin;
  	//cout << bin << endl;
      }
      //cout << "For param1, min chisq  = " << min_chisq_nh[itoy] << endl;
      //cout << endl;
    }
    for(int i = 1; i <= nbins1_nh; i++){
      double bin = h_chisq_par1_Ave_nh->GetBinContent(i);
      if (bin < min_chisq_Ave_nh) min_chisq_Ave_nh = bin;
    }
    cout << "Global scanned chi2 minimum = " << min_chisq_Ave_nh << endl;     
  }
  else if(nDim==2){
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      min_chisq_nh[itoy] = 999.;
      //cout << "# Toy " << itoy << endl;    
      for(int i1 = 1; i1 <= nbins1_nh; i1++){
  	for(int i2 = 1; i2 <= nbins2_nh; i2++){
  	  double bin = h_chisq_par2_par1_nh[itoy]->GetBinContent(i1,i2);
  	  if (bin < min_chisq_nh[itoy]) min_chisq_nh[itoy] = bin;
  	  //cout << bin << endl;
  	}
      }
      //cout << "For param1 Vs param2, min chisq  = " << min_chisq_nh[itoy] << endl;
      //cout << endl;
    }

    for(int i1 = 1; i1 <= nbins1_nh; i1++){
      for(int i2 = 1; i2 <= nbins2_nh; i2++){
  	double bin = h_chisq_par2_par1_Ave_nh->GetBinContent(i1,i2);
  	if (bin < min_chisq_Ave_nh) min_chisq_Ave_nh = bin;
	//cout << bin << " ";
      }
      //cout << endl;
    }
    cout << "Global scanned chi2 minimum = " << min_chisq_Ave_nh << endl; 
  } // nDim==2
  



  // 1D: dChi2 vs param1  (Inverted Hierarchy)

  cout << endl;
  cout << "Read file: " << filename_ih << endl;
  cout << endl;

  const char * filename_ih = Form("%s.root",name_ih);
  TChain ch_ih("fitnt");

  if(ChangeOrder_ih){
    ch_ih.Add(filename_ih);
    ch_ih.SetBranchAddress("best_fit_params", &fparm_ih); 
    ch_ih.SetBranchAddress("nbins_param1", &nbins_param2_ih); 
    ch_ih.SetBranchAddress("nbins_param2", &nbins_param1_ih); 
    ch_ih.SetBranchAddress("nbins_tot_param1", &nbins_tot_param2_ih); 
    ch_ih.SetBranchAddress("nbins_tot_param2", &nbins_tot_param1_ih); 
    ch_ih.SetBranchAddress("best_fit_params", &fparm_ih); 
    ch_ih.SetBranchAddress("parname1", &parname2_ih); // name of the grid parameter 2
    ch_ih.SetBranchAddress("parname2", &parname1_ih); // name of the grid parameter 1
    ch_ih.SetBranchAddress("max_val_param1",&max_val_param2);
    ch_ih.SetBranchAddress("min_val_param1",&min_val_param2);
    ch_ih.SetBranchAddress("max_val_param2",&max_val_param1);
    ch_ih.SetBranchAddress("min_val_param2",&min_val_param1);
    ch_ih.SetBranchAddress("step1",&step2);
    ch_ih.SetBranchAddress("step2",&step1);
  }
  else{
    ch_ih.Add(filename_ih);
    ch_ih.SetBranchAddress("best_fit_params", &fparm_ih); 
    ch_ih.SetBranchAddress("nbins_param1", &nbins_param1_ih); 
    ch_ih.SetBranchAddress("nbins_param2", &nbins_param2_ih); 
    ch_ih.SetBranchAddress("nbins_tot_param1", &nbins_tot_param1_ih); 
    ch_ih.SetBranchAddress("nbins_tot_param2", &nbins_tot_param2_ih); 
    ch_ih.SetBranchAddress("best_fit_params", &fparm_ih); 
    ch_ih.SetBranchAddress("parname1", &parname1_ih); // name of the grid parameter 1
    ch_ih.SetBranchAddress("parname2", &parname2_ih); // name of the grid parameter 2
    ch_ih.SetBranchAddress("max_val_param1",&max_val_param1);
    ch_ih.SetBranchAddress("min_val_param1",&min_val_param1);
    ch_ih.SetBranchAddress("max_val_param2",&max_val_param2);
    ch_ih.SetBranchAddress("min_val_param2",&min_val_param2);
    ch_ih.SetBranchAddress("step1",&step1);
    ch_ih.SetBranchAddress("step2",&step2);
  }
  
  if(*parname1_ih==*parname2_ih && *parname1_ih!=""){
   cout << "Name parameter 1 == Name parameter 2 = " << *parname1_ih << endl; 
   exit(0);
  }

  const int nentries_ih = ch_ih.GetEntries();
  cout << "nentries (IH) = " << nentries_ih << endl;
  
  ch_ih.GetEntry(0);
  int nbins1_ih = TMath::Max(nbins_tot_param1_ih,nbins_param1_ih); // if using make_toygrid_toys_joint.C nbins_tot_param* is used instead of nbins_param*
  int nbins2_ih = TMath::Max(nbins_tot_param2_ih,nbins_param2_ih); // otherwise nbins_tot_param* = 1 < nbins_param*
  
  int nDim=0;
  if(nbins1_ih>1 || nbins2_ih>1) nDim = 1;
  if(nbins1_ih>1 && nbins2_ih>1) nDim = 2;
 
  if(*parname1_ih!=""){
    cout << endl;
    cout << "- Parameter 1: " << *parname1_ih << endl;
    cout << "  # bins = " << nbins1_ih << endl;
    cout << "  min = " << min_val_param1 << endl;
    cout << "  max = " << max_val_param1 << endl;
    cout << "  step = " << step1 << endl;
    cout << endl;
  } 
  if(*parname2_ih!=""){
  cout << "- Parameter 2: " << *parname2_ih << endl;
  cout << "  # bins = " << nbins2_ih << endl;
  cout << "  min = " << min_val_param2 << endl;
  cout << "  max = " << max_val_param2 << endl;
  cout << "  step = " << step2 << endl;
  cout << endl;
  }

  if(*parname1_ih=="" && *parname2_ih==""){
    cout << "# bins parameter 1 = " << nbins1_ih << endl; 
    cout << "# bins parameter 2 = " << nbins2_ih << endl; 
    cout << endl;
  }

  const int ntoys_ih = nentries_ih / (nbins1_ih * nbins2_ih);
  cout << "ntoys (IH) = " << ntoys_ih << endl;
  cout << endl;
  
  //Set up bins and histograms
  double mixing13_min = 0;
  double mixing13_max = 0;
  double mixing13_step = 0;
  double mixing23_min = 0;
  double mixing23_max = 0;
  double mixing23_step = 0;
  double dmsq_min = 0;
  double dmsq_max = 0;
  double dmsq_step = 0;
  double deltacp_min = 0;
  double deltacp_max = 0;
  double deltacp_step = 0;

  if(!theta2theta){
    if(reactor){
      mixing13_min = 0.015;
      mixing13_max = 0.035;
      //mixing13_step = 0.001;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_ih-1);
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.1;
      //mixing13_step = 0.005;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_ih-1);
    }
  }
  else{
    if(reactor){
      mixing13_min = 0.06;
      mixing13_max = 0.14;
      //mixing13_step = 0.004;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_ih-1);
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.4;
      //mixing13_step = 0.02;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_ih-1);
    }
  }
  mixing23_min = 0.3;
  mixing23_max = 0.7;
  //mixing23_step = 0.01;
  mixing23_step = (mixing23_max-mixing23_min)/(nbins1_ih-1);
  
  dmsq_min = 2.0;//E-3;
  dmsq_max = 3.0;//E-3;
  //dmsq_step = 0.05;//E-3;
  dmsq_step = (dmsq_max-dmsq_min)/(nbins1_ih-1);

  deltacp_min = -conv::kpi;
  deltacp_max = conv::kpi;
  //deltacp_step = 2*conv::kpi/50;
  deltacp_step = (deltacp_max-deltacp_min)/(nbins1_ih-1);
  

  // Take steps and min-max values from the tree if they exist

  if(*parname1_ih==""){
    
    parname1_ih[0].append(parname1_def); // fill parname2 w/ default if not in the tree
    
    if(parname1_def=="deltacp"){
      min_val_param1 = deltacp_min;
      max_val_param1 = deltacp_max;
      step1 = deltacp_step;
    }
    else if(parname1_def=="mixing13"){
      min_val_param1 = mixing13_min;
      max_val_param1 = mixing13_max;
      step1 = mixing13_step;
    }
    else if(parname1_def=="mixing23"){
      min_val_param1 = mixing23_min;
      max_val_param1 = mixing23_max;
      step1 = mixing23_step;
    }
    else if(parname1_def=="dmsq"){
      min_val_param1 = dmsq_min;
      max_val_param1 = dmsq_max;
      step1 = dmsq_step;
    }
    else if(parname1_def=="none"){
      min_val_param1 = 0.;
      max_val_param1 = 0.;
      step1 = 0.;
    }
  }
  
  if(*parname2_ih==""){

    parname2_ih[0].append(parname2_def); // fill parname2 w/ default if not in the tree

    if(parname2_def=="deltacp"){
      min_val_param2 = deltacp_min;
      max_val_param2 = deltacp_max;
      step2 = deltacp_step;
    }
    else if(parname2_def=="mixing13"){
      min_val_param2 = mixing13_min;
      max_val_param2 = mixing13_max;
      step2 = mixing13_step;
    }
    else if(parname2_def=="mixing23"){
      min_val_param2 = mixing23_min;
      max_val_param2 = mixing23_max;
      step2 = mixing23_step;
    }
    else if(parname2_def=="dmsq"){
      min_val_param2 = dmsq_min;
      max_val_param2 = dmsq_max;
      step2 = dmsq_step;
    }
    else if(parname2_def=="none"){
    }
  }
    
  //
  char histoname[200];

  // 1D histograms

  Double_t** bestfit_points_ih_x = new Double_t*[ntoys_ih]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_ih; itoy++){
    bestfit_points_ih_x[itoy] = new Double_t[nbins1_ih];
  }
  Double_t** bestfit_points_ih_y = new Double_t*[ntoys_ih]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_ih; itoy++){
    bestfit_points_ih_y[itoy] = new Double_t[nbins1_ih];
  }

  sprintf(histoname,"h_chisq_%s_Ave_ih",parname1_ih[0].c_str());
  TH1D * h_chisq_par1_Ave_ih = new TH1D(histoname,histoname,nbins1_ih,min_val_param1-0.5*step1,max_val_param1+0.5*step1);

  TH1D ** h_chisq_par1_ih = new TH1D*[ntoys_ih];
  for(int itoy=0; itoy<ntoys_ih; itoy++){
    sprintf(histoname,"h_%s_ih_Toy%d",parname1_ih[0].c_str(),itoy);
    h_chisq_par1_ih[itoy] = new TH1D(histoname,histoname,nbins1_ih,min_val_param1-0.5*step1,max_val_param1+0.5*step1);
  }  
  
  TH1D ** h_chisq_par1_ih_bin = new TH1D*[nbins1_ih];
  for(int ibin=0;ibin<nbins1_ih;ibin++){
    char histoname[200];
    sprintf(histoname,"h_chisq_%s_ih_bin%d",parname1_ih[0].c_str(),ibin);
    h_chisq_par1_ih_bin[ibin] = new TH1D(histoname,histoname,50,0,500);    
  }

  // 2D histograms

  sprintf(histoname,"h_chisq_%s_%s_Ave_ih",parname2_ih[0].c_str(),parname1_ih[0].c_str());
  TH2D * h_chisq_par2_par1_Ave_ih = new TH2D(histoname,histoname,nbins1_ih,min_val_param1-0.5*step1,max_val_param1+0.5*step1,nbins2_ih,min_val_param2-0.5*step2,max_val_param2+0.5*step2);

  TH2D ** h_chisq_par2_par1_ih= new TH2D*[ntoys_ih];
  for(int itoy=0; itoy<ntoys_ih; itoy++){
    char histoname[200];
    sprintf(histoname,"h_chisq_%s_%s_ih_Toy%d",parname2_ih[0].c_str(),parname1_ih[0].c_str(),itoy);
    h_chisq_par2_par1_ih[itoy] = new TH2D(histoname,histoname,nbins1_ih,min_val_param1-0.5*step1,max_val_param1+0.5*step1,nbins2_ih,min_val_param2-0.5*step2,max_val_param2+0.5*step2);
  }  

  TH1D *** h_chisq_par2_par1_ih_bin = new TH1D**[nbins1_ih];
  for(int ibin1=0;ibin1<nbins1_ih;ibin1++){
    h_chisq_par2_par1_ih_bin[ibin1] = new TH1D*[nbins2_ih];
    for(int ibin2=0;ibin2<nbins2_ih;ibin2++){
      char histoname[200];
      sprintf(histoname,"h_chisq_%s_%s_ih_bin%d_bin%d",parname2_ih[0].c_str(),parname1_ih[0].c_str(),ibin2,ibin1);
      h_chisq_par2_par1_ih_bin[ibin1][ibin2] = new TH1D(histoname,histoname,50,0,500);    
    }
  }

  // Loop over toys in the tree

  int iEvtToy = 0;
  for(int i = 0; i < nentries_ih; i++){
    iEvtToy = i / (nbins1_ih*nbins2_ih);
    if( !(i%(nbins1_ih*nbins2_ih)) ) cout << "# Toy " << iEvtToy << endl;
        
    ch_ih.GetEntry(i);
    int migradminstat_ih = fparm_ih->GetStatusCode("migrad_min");
    if (migradminstat_ih == 4) cout << "Fit failed!" << endl;
    //cout << "migradminstat_ih = " << migradminstat_ih << endl;

    double mix13_ih = fparm_ih->BestFitValue(conv::kPhypMixing13); // Parameters are fixed in the fit
    double mix23_ih = fparm_ih->BestFitValue(conv::kPhypMixing23);
    double dmsq_ih = fparm_ih->BestFitValue(conv::kPhypDm23sq);
    double deltacp_ih = fparm_ih->BestFitValue(conv::kPhypCPV);
    double chisq_ih = fparm_ih->GetBestFitStatistic();
    if(NuBarMode){
      mix13_ih = fparm_ih->BestFitValue(conv::kPhypMixing13Bar); // Parameters are fixed in the fit
      mix23_ih = fparm_ih->BestFitValue(conv::kPhypMixing23Bar);
      dmsq_ih = fparm_ih->BestFitValue(conv::kPhypDm23sqBar);
      deltacp_ih = fparm_ih->BestFitValue(conv::kPhypCPVBar);
      chisq_ih = fparm_ih->GetBestFitStatistic();
    }

    double bf1_ih=0;
    if(*parname1_ih=="deltacp") bf1_ih = deltacp_ih;
    else if(*parname1_ih=="mixing13") bf1_ih = mix13_ih;
    else if(*parname1_ih=="mixing23") bf1_ih = mix23_ih;
    else if(*parname1_ih=="dmsq") bf1_ih = dmsq_ih;
    double bf2_ih=0;
    if(*parname2_ih=="deltacp") bf2_ih = deltacp_ih;
    else if(*parname2_ih=="mixing13") bf2_ih = mix13_ih;
    else if(*parname2_ih=="mixing23") bf2_ih = mix23_ih;
    else if(*parname2_ih=="dmsq") bf2_ih = dmsq_ih;

    int bin_idx=0;
    double bincont=0.;

    if(nbins1_ih>1 && nbins2_ih<=1){ // 1D param1
      bin_idx = h_chisq_par1_ih[iEvtToy]->FindBin(bf1_ih);
      bincont = h_chisq_par1_ih[iEvtToy]->GetBinContent(bin_idx);
    }
    else if(nbins1_ih<=1 && nbins2_ih>1){ // 1D param2 
      bin_idx = h_chisq_par1_ih[iEvtToy]->FindBin(bf1_ih);
      bincont = h_chisq_par1_ih[iEvtToy]->GetBinContent(bin_idx);
    }
    else if(nbins1_ih>1 && nbins2_ih>1){ // 2D param2 Vs param1
      bin_idx = h_chisq_par2_par1_ih[iEvtToy]->FindBin(bf1_ih,bf2_ih);
      bincont = h_chisq_par2_par1_ih[iEvtToy]->GetBinContent(bin_idx);
    }
  
    if(bincont!=0){
      if(bincont>chisq_ih){
	if(chisq_ih!=0){
	  if(nDim==1) h_chisq_par1_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
	  else if(nDim==2) h_chisq_par2_par1_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
	  //bestfit_points_ih_x[iEvtToy][bin_idx-1]  = bf1_ih;//sinsq2theta13;
	  //bestfit_points_ih_y[iEvtToy][bin_idx-1]  = bf2_ih;
	}
      }
    } 
    else{ 
      if(chisq_ih!=0){
	if(nDim==1) h_chisq_par1_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
	else if(nDim==2) h_chisq_par2_par1_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
	//bestfit_points_ih_x[iEvtToy][bin_idx-1] = bf1_ih;//sinsq2theta13;
	//bestfit_points_ih_y[iEvtToy][bin_idx-1] = bf2_ih;
      }
    }
  } // end loop over toys
 
  // TGraph ** gbestfit_ih_Toy = new TGraph*[ntoys_ih]; 
  // for(int itoy=0;itoy<ntoys_ih;itoy++){
  //   char graphname[200];
  //   sprintf(graphname,"gbestfit_ih_Toy%d",itoy);
  //   gbestfit_ih_Toy[itoy] = new TGraph( nbins1_ih, bestfit_points_ih_x[itoy] ,bestfit_points_ih_y[itoy] ); // dcp Vs mix13  
  //   gbestfit_ih_Toy[itoy]->SetName(graphname);
  //   //cout << graphname << endl;
  // }  

  // Get the toy distribution for each param value
  if(nDim==1){
    cout << endl << endl;
    for(int ibin=0;ibin<nbins1_ih;ibin++){
      for(int itoy=0;itoy<ntoys_ih;itoy++){
	double bincont = h_chisq_par1_ih[itoy]->GetBinContent(ibin+1);
	h_chisq_par1_ih_bin[ibin]->Fill(bincont);
	cout << bincont << " ";
      }
    }
    cout << endl << endl;
  }
  if(nDim==2){
    for(int ibin1=0;ibin1<nbins1_ih;ibin1++){
      for(int ibin2=0;ibin2<nbins2_ih;ibin2++){
  	for(int itoy=0;itoy<ntoys_ih;itoy++){
	  int ibin = h_chisq_par2_par1_ih[itoy]->GetBin(ibin1+1,ibin2+1);
  	  double bincont = h_chisq_par2_par1_ih[itoy]->GetBinContent(ibin);
	  h_chisq_par2_par1_ih_bin[ibin1][ibin2]->Fill(bincont);
  	}
      }
    }
  }
  
  // Compute the average sensitivity    
  for(int itoy=0;itoy<ntoys_ih;itoy++){
    h_chisq_par1_Ave_ih->Add(h_chisq_par1_ih[itoy]);
    h_chisq_par2_par1_Ave_ih->Add(h_chisq_par2_par1_ih[itoy]);
  }
  h_chisq_par1_Ave_ih->Scale(1./ntoys_ih);
  h_chisq_par2_par1_Ave_ih->Scale(1./ntoys_ih);
  
  // Find the minimum chi2 (global scan)
  
  cout << endl;
  cout << "Find the minimum chi2 (global scan)" << endl;
  cout << endl;

  double *min_chisq_ih = new double[ntoys_ih];
  double min_chisq_Ave_ih=999.;
    
  if(nDim==1){ 
    for(int itoy=0;itoy<ntoys_ih;itoy++){
      min_chisq_ih[itoy] = 999.;
      //cout << "# Toy " << itoy << endl;    
      for(int i = 1; i <= nbins1_ih; i++){
  	double bin = h_chisq_par1_ih[itoy]->GetBinContent(i);
  	if (bin < min_chisq_ih[itoy]) min_chisq_ih[itoy] = bin;
  	//cout << bin << endl;
      }
      //cout << "For param1, min chisq  = " << min_chisq_ih[itoy] << endl;
      //cout << endl;
    }
    for(int i = 1; i <= nbins1_ih; i++){
      double bin = h_chisq_par1_Ave_ih->GetBinContent(i);
      if (bin < min_chisq_Ave_ih) min_chisq_Ave_ih = bin;
    }
    cout << "Global scanned chi2 minimum = " << min_chisq_Ave_ih << endl;     
  } // nDim==1
  else if(nDim==2){
    for(int itoy=0;itoy<ntoys_ih;itoy++){
      min_chisq_ih[itoy] = 999.;
      //cout << "# Toy " << itoy << endl;    
      for(int i1 = 1; i1 <= nbins1_ih; i1++){
  	for(int i2 = 1; i2 <= nbins2_ih; i2++){
  	  double bin = h_chisq_par2_par1_ih[itoy]->GetBinContent(i1,i2);
  	  if (bin < min_chisq_ih[itoy]) min_chisq_ih[itoy] = bin;
  	  //cout << bin << endl;
  	}
      }
      //cout << "For param1 Vs param2, min chisq  = " << min_chisq_ih[itoy] << endl;
      //cout << endl;
    }

    for(int i1 = 1; i1 <= nbins1_ih; i1++){
      for(int i2 = 1; i2 <= nbins2_ih; i2++){
  	double bin = h_chisq_par2_par1_Ave_ih->GetBinContent(i1,i2);
  	if (bin < min_chisq_Ave_ih) min_chisq_Ave_ih = bin;
	cout << bin << " ";
      }
      cout << endl;
    }
    cout << "Global scanned chi2 minimum = " << min_chisq_Ave_ih << endl; 
  } // nDim==2
  



  ///////////////////



  // Shift the chi2 wrt the minimum
  double ShiftMinGlobal = 99999.;
  if(doShiftMin){
    ShiftMinGlobal = TMath::Min(min_chisq_Ave_ih,min_chisq_Ave_nh);
    cout << "ShiftMinGlobal = " << ShiftMinGlobal << endl;
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      ShiftMinGlobal = TMath::Min(min_chisq_nh[itoy],min_chisq_ih[itoy]);
    }
  }
  else 
    ShiftMinGlobal = 0.;

  // NH
  
  if(nDim==1){
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      for(int i = 1; i <= nbins1_nh; i++){
	double bin = h_chisq_par1_nh[itoy]->GetBinContent(i);
	h_chisq_par1_nh[itoy]->SetBinContent(i,bin-ShiftMinGlobal); // scanned global minimum
	//h_chisq_par1_nh[itoy]->SetBinContent(i,bin-min_chisq_nh[itoy]); // scanned global minimum
	//h_chisq_par1_nh[itoy]->SetBinContent(i,bin-chisq_bf_nh); // fiited global minimum
      }
    }
    for(int i = 1; i <= nbins1_nh; i++){
      double bin = h_chisq_par1_Ave_nh->GetBinContent(i);
      h_chisq_par1_Ave_nh->SetBinContent(i,bin-ShiftMinGlobal); // scanned global minimum
      //h_chisq_par1_Ave_nh->SetBinContent(i,bin-min_chisq_Ave_nh); // scanned global minimum
    }
  }   
  else if(nDim==2){
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      for(int i1 = 1; i1 <= nbins1_nh; i1++){
	for(int i2 = 1; i2 <= nbins2_nh; i2++){
	  double bin = h_chisq_par2_par1_nh[itoy]->GetBinContent(i1,i2);
	  h_chisq_par2_par1_nh[itoy]->SetBinContent(i1,i2,bin-ShiftMinGlobal); // scanned global minimum
	  //h_chisq_par2_par1_nh[itoy]->SetBinContent(i1,i2,bin-min_chisq_nh[itoy]); // scanned global minimum
	  //h_chisq_par1_nh[itoy]->SetBinContent(i,bin-chisq_bf_nh); // fiited global minimum
	}
      }
    }
    for(int i1 = 1; i1 <= nbins1_nh; i1++){
      for(int i2 = 1; i2 <= nbins2_nh; i2++){
	double bin = h_chisq_par2_par1_Ave_nh->GetBinContent(i1,i2);
	h_chisq_par2_par1_Ave_nh->SetBinContent(i1,i2,bin-ShiftMinGlobal); // scanned global minimum
	//h_chisq_par2_par1_Ave_nh->SetBinContent(i1,i2,bin-min_chisq_Ave_nh); // scanned global minimum
      }
    }
  }
  
  // IH
  
  if(nDim==1){
    for(int itoy=0;itoy<ntoys_ih;itoy++){
      for(int i = 1; i <= nbins1_ih; i++){
	double bin = h_chisq_par1_ih[itoy]->GetBinContent(i);
	h_chisq_par1_ih[itoy]->SetBinContent(i,bin-ShiftMinGlobal); // scanned global minimum
	//h_chisq_par1_ih[itoy]->SetBinContent(i,bin-min_chisq_ih[itoy]); // scanned global minimum
	//h_chisq_par1_ih[itoy]->SetBinContent(i,bin-chisq_bf_ih); // fitted global minimum
      }
    }
    for(int i = 1; i <= nbins1_ih; i++){
      double bin = h_chisq_par1_Ave_ih->GetBinContent(i);
      h_chisq_par1_Ave_ih->SetBinContent(i,bin-ShiftMinGlobal); // scanned global minimum
      //h_chisq_par1_Ave_ih->SetBinContent(i,bin-min_chisq_Ave_ih); // scanned global minimum
    }
  }
  else if(nDim==2){
    for(int itoy=0;itoy<ntoys_ih;itoy++){
      for(int i1 = 1; i1 <= nbins1_ih; i1++){
	for(int i2 = 1; i2 <= nbins2_ih; i2++){
	  double bin = h_chisq_par2_par1_ih[itoy]->GetBinContent(i1,i2);
	  h_chisq_par2_par1_ih[itoy]->SetBinContent(i1,i2,bin-ShiftMinGlobal); // scanned global minimum
	  //h_chisq_par2_par1_ih[itoy]->SetBinContent(i1,i2,bin-min_chisq_ih[itoy]); // scanned global minimum
	  //h_chisq_par1_ih[itoy]->SetBinContent(i,bin-chisq_bf_ih); // fitted global minimum
	}
      }
    }
    for(int i1 = 1; i1 <= nbins1_ih; i1++){
      for(int i2 = 1; i2 <= nbins2_ih; i2++){
	double bin = h_chisq_par2_par1_Ave_ih->GetBinContent(i1,i2);
	h_chisq_par2_par1_Ave_ih->SetBinContent(i1,i2,bin-ShiftMinGlobal); // scanned global minimum
	//h_chisq_par2_par1_Ave_ih->SetBinContent(i1,i2,bin-min_chisq_Ave_ih); // scanned global minimum
      }
    }
  }

  




  // Draw plots

  // Normal Hierarchy

  // 2-dimensional plots
  if(nDim==2){

    if (!NuBarMode) {
      h_chisq_par2_par1_Ave_nh->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
      h_chisq_par2_par1_Ave_nh->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
      h_chisq_par2_par1_Ave_ih->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
      h_chisq_par2_par1_Ave_ih->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
    }
    else {
      h_chisq_par2_par1_Ave_nh->GetXaxis()->SetTitle("sin^{2}#bar{#theta_{23}}");
      h_chisq_par2_par1_Ave_nh->GetYaxis()->SetTitle("#bar{#Delta m^{2}_{32}}");
      h_chisq_par2_par1_Ave_ih->GetXaxis()->SetTitle("sin^{2}#bar{#theta_{23}}");
      h_chisq_par2_par1_Ave_ih->GetYaxis()->SetTitle("#bar{#Delta m^{2}_{32}}");
    }

    TCanvas *c_chisq_par2_par1_nh = new TCanvas("c_chisq_par2_par1_nh","c_chisq_par2_par1_nh");
    h_chisq_par2_par1_Ave_nh->Draw("colz");
    gPad->Update();
    TPaletteAxis *palette_nh = (TPaletteAxis*)h_chisq_par2_par1_Ave_nh->GetListOfFunctions()->FindObject("palette");
    palette_nh->SetX1NDC(0.9);
    palette_nh->SetX2NDC(0.95);
    //palette_nh->SetY1NDC(0.);
    //palette_nh->SetY2NDC(.8);
    gPad->Modified();
    gPad->Update();
    //h_chisq_par2_par1_Ave_nh->DrawClone("colz");

    
    TCanvas *c_chisq_par2_par1_ih = new TCanvas("c_chisq_par2_par1_ih","c_chisq_par2_par1_ih");
    h_chisq_par2_par1_Ave_ih->Draw("colz");
    gPad->Update();
    TPaletteAxis *palette_ih = (TPaletteAxis*)h_chisq_par2_par1_Ave_ih->GetListOfFunctions()->FindObject("palette");
    palette_ih->SetX1NDC(0.9);
    palette_ih->SetX2NDC(0.95);
    // //palette_ih->SetY1NDC(0.);
    // //palette_ih->SetY2NDC(.8);
    gPad->Modified();
    gPad->Update();
    //h_chisq_par2_par1_Ave_ih->DrawClone("colz");

    TCanvas *c_cont_par2_par1 = new TCanvas("c_cont_par2_par1","c_cont_par2_par1");
    c_cont_par2_par1->cd();
    c_cont_par2_par1->SetGridx();
    c_cont_par2_par1->SetGridy();

    h_chisq_par2_par1_Ave_nh_Cont = (TH2F*)h_chisq_par2_par1_Ave_nh->Clone("");
    h_chisq_par2_par1_Ave_nh_Cont->GetXaxis()->SetTitle("sin^{2}#bar{#theta}_{23}}");
    h_chisq_par2_par1_Ave_nh_Cont->GetYaxis()->SetTitle("#bar{#Delta m^{2}}_{32}");
    h_chisq_par2_par1_Ave_nh_Cont->SetContour(1);
    h_chisq_par2_par1_Ave_nh_Cont->SetContourLevel(0,CLval);
    h_chisq_par2_par1_Ave_nh_Cont->SetLineColor(kBlue);
    h_chisq_par2_par1_Ave_nh_Cont->Draw("cont3");
 
    h_chisq_par2_par1_Ave_ih_Cont = (TH2F*)h_chisq_par2_par1_Ave_ih->Clone("");
    h_chisq_par2_par1_Ave_ih_Cont->SetContour(1);
    h_chisq_par2_par1_Ave_ih_Cont->SetLineColor(kBlack);
    h_chisq_par2_par1_Ave_ih_Cont->SetContourLevel(0,CLval);
    h_chisq_par2_par1_Ave_ih_Cont->Draw("cont3 same");

    if(DrawLeg){
      TLegend *l_cont_par2_par1 = new TLegend(0.7,0.4,0.89,0.69);
      l_cont_par2_par1->SetFillColor(0);
      l_cont_par2_par1->SetLineColor(0);
      l_cont_par2_par1->AddEntry(h_chisq_par2_par1_Ave_nh_Cont,LegendName_nh,"L");
      l_cont_par2_par1->AddEntry(h_chisq_par2_par1_Ave_ih_Cont,LegendName_ih,"L");
      l_cont_par2_par1->Draw();
    }
  }
  // 1-dimensional plots
  else if(nDim==1){
    
    TCanvas *c_chisq_par1_nh = new TCanvas("c_chisq_par1_nh","c_chisq_par1_nh");
    double min = TMath::Min(0.,(double)h_chisq_par1_Ave_nh->GetMinimum());
    min = TMath::Max(0.,min);
    double max = h_chisq_par1_Ave_nh->GetMaximum() * 1.25;
    h_chisq_par1_Ave_nh->SetMinimum(max);
    h_chisq_par1_Ave_nh->SetMinimum(0);
    h_chisq_par1_Ave_nh->GetXaxis()->SetRangeUser(0,0.4);
    h_chisq_par1_Ave_nh->DrawClone("histo");

    TCanvas *c_chisq_par1_ih = new TCanvas("c_chisq_par1_ih","c_chisq_par1_ih");
    double min = TMath::Min(0.,(double)h_chisq_par1_Ave_ih->GetMinimum());
    min = TMath::Max(0.,min);
    double max = h_chisq_par1_Ave_ih->GetMaximum() * 1.25;
    h_chisq_par1_Ave_ih->SetMinimum(max);
    h_chisq_par1_Ave_ih->SetMinimum(0);
    h_chisq_par1_Ave_ih->GetXaxis()->SetRangeUser(0,0.4);
    h_chisq_par1_Ave_ih->DrawClone("histo");

    //

    TCanvas *c_chisq_par1_mhall = new TCanvas("c_chisq_par1_mhall","c_chisq_par1_mhall");
    h_chisq_par1_Ave_nh->SetMaximum(10);
    h_chisq_par1_Ave_nh->SetMinimum(0);
    h_chisq_par1_Ave_nh->SetLineColor(kBlack);
    h_chisq_par1_Ave_nh->SetLineWidth(2);
    if(*parname1_nh=="deltacp") h_chisq_par1_Ave_nh->GetXaxis()->SetTitle("#delta_{CP}");
    h_chisq_par1_Ave_nh->GetYaxis()->SetTitle("#Delta #chi^{2}");
    h_chisq_par1_Ave_nh->DrawClone("histo");

    h_chisq_par1_Ave_ih->SetLineColor(kRed);
    h_chisq_par1_Ave_ih->SetLineWidth(2);
    h_chisq_par1_Ave_ih->GetXaxis()->SetTitle(parname2_nh[0].c_str());
    h_chisq_par1_Ave_ih->DrawClone("histo same");
    
    TLine *lineConstCL = new TLine(0,CLval,0.4,CLval);
    lineConstCL->SetLineWidth(2);
    lineConstCL->Draw("same");
    
    if(DrawLeg){
      TLegend *l_cont_par1 = new TLegend(0.6,0.4,0.89,0.69);
      l_cont_par1->SetFillColor(0);
      l_cont_par1->SetLineColor(0);
      l_cont_par1->AddEntry(h_chisq_par1_Ave_nh,LegendName_nh,"L");
      l_cont_par1->AddEntry(h_chisq_par1_Ave_ih,LegendName_ih,"L");
      l_cont_par1->Draw();
    }
    
  
  } // NDim == 1


  // Print plots
  
  if(PrintPlots){
    
    char name[200];
    
    if(nDim==1){
      sprintf(name,"%s%s%s","Like1D_NH_",tag,".pdf");
      c_chisq_par1_nh->Print(name);
      sprintf(name,"%s%s%s","Like1D_IH_",tag,".pdf");
      c_chisq_par1_ih->Print(name);
      sprintf(name,"%s%s%s","Cont1D_",tag,".pdf");
      c_chisq_par1_mhall->Print(name);
     }
    else if(nDim==2){
      sprintf(name,"%s%s%s","Like2D_NH_",tag,".pdf");
      c_chisq_par2_par1_nh->Print(name);
      sprintf(name,"%s%s%s","Like2D_IH_",tag,".pdf");
      c_chisq_par2_par1_ih->Print(name);
      sprintf(name,"%s%s%s","Cont2D_",tag,".pdf");
      c_cont_par2_par1->Print(name);
    }
    
  }


  // Write histograms to root file
  if (output_file) {
    TFile hist_file(output_file,"RECREATE");
    hist_file.cd();
    if(nDim==2){
      h_chisq_par2_par1_Ave_nh->Write();
      h_chisq_par2_par1_Ave_ih->Write();
    }
    if(nDim==1) {
      h_chisq_par1_Ave_nh->Write();
      h_chisq_par1_Ave_ih->Write();
    }   
    cout<<"Histograms written to "<<output_file<<endl;
  }
 
  return;
}

  
void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style){
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
  contour->SetContour(1);//number of contour levels                                                                                                                                    
  contour->SetContourLevel(0,0.0);//SetContourLevel(int level, double value)                                                                                                           
  contour->SetStats(0);
} 
