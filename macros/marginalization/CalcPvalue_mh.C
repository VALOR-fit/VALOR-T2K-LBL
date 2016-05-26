
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_mh_marginalization.cxx
// Compute p-values using deltachisq01_glob
// Store final results in a output root file to be used by other macros
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

void CalcPvalue_mh
(
 // Fit of Run 1-6. Best-fit = Normal Hierarchy
 TString mytag = "Pvalue_IH_100kToys_Run16_Asimov_B",

 //const char * rootfile = "data_run1to6/Asimov_A_Run16.root", 
 const char * rootfile = "data_run1to6/Asimov_B_Run16.root", 
 //const char * rootfile = "data_run1to6/SKdata_Run16.root", 
 
 const bool use_posterior_predictive = true,
 const bool use_post_pred_fhc_mulike = true,
 const bool use_post_pred_fhc_elike  = true,
 const bool use_post_pred_rhc_mulike = true,
 const bool use_post_pred_rhc_elike  = true,

 //double deltachisq_nh_ih_data = -2.0816342, // chisq(NH)-chisq(IH) - Run 1-6 data set
 //double deltachisq_nh_ih_data = -1.0468992, // Asimov A (Normal Hierarchy)
 double deltachisq_nh_ih_data = 0.03464154, // Asimov B (Normal Hierarchy)
 //double deltachisq_nh_ih_data = -0.28990745, // Asimov A (dcp=0, Normal Hierarchy)
 
 //
 // Standard fake data sets
 //
 
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_nh_normal_nexp_10k_beta1.root" // 10k fake data sets, 100k marg toys, Normal Hierarchy
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_ih_inverted_nexp_10k_beta1.root" // 10k fake data sets, 100k marg toys, Inverted Hierarchy
 
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_nh_2_normal_nexp_100k_1.root" // 10k fake data sets, 100k marg toys, Normal Hierarchy
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_ih_2_inverted_nexp_100k_1.root" // 10k fake data sets, 100k marg toys, Inverted Hierarchy
 
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/50k/fit_fitmarg_masshier_20k_ih_all_inverted_nexp_100k_1.root" // 20k fake data sets, 100k marg toys, Inverted Hierarchy

 //
 // Many toys
 //
 
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/50k/fit_fitmarg_masshier_50k_nh_normal_1.000.root", // 50k toys
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/50k/fit_fitmarg_masshier_50k_ih_nexp_100_1.000.root", // 50k toys

 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/50k/fit_fitmarg_masshier_50k_nh_normal_1.000_other.root", // 50 toys
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/50k/fit_fitmarg_masshier_50k_ih_inverted_1.000_other.root", // 50k toys

 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/50k/fit_fitmarg_masshier_50k_normal_Tot100k.root",  // 100k toys
 char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/50k/fit_fitmarg_masshier_50k_inverted_Tot100k.root", // 100k toys
 
 //
 // Fix deltaCP in fake data sets generation
 //
 
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_nh_dcp_0_normal_nexp_10k_1.000.root" // 10k fake data sets, 100k marg toys, Normal Hierarchy, fix dcp=0
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_ih_dcp_0_inverted_nexp_10k_1.000.root" // 10k fake data sets, 100k marg toys, Inverted Hierarchy, fix dcp=0
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_nh_dcp_m_pi_over_2_normal_nexp_10k_1.000.root" // 10k fake data sets, 100k marg toys, Normal Hierarchy, fix dcp=-pi/2
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_ih_dcp_m_pi_over_2_inverted_nexp_10k_1.000.root" // 10k fake data sets, 100k marg toys, Inverted Hierarchy, fix dcp=-pi/2

 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_nh_dcp_p_pi_over_2_normal_nexp_10k_1.000.root" // 10k fake data sets, 100k marg toys, Normal Hierarchy, fix dcp=+pi/2
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_ih_dcp_p_pi_over_2_inverted_nexp_10k_1.000.root" // 10k fake data sets, 100k marg toys, Inverted Hierarchy, fix dcp=+pi/2

 //
 // Single fit
 //

 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_data_normal_nexp_1_1.000_-999.000-job0-_0_0.root" // Run 1-6 data
 )
{

  // Set data Feldman-Cousins deltachisq later from deltachisq_nh_ih 
  double deltachisq_fc_data = conv::kBadResult; // chisq(True)-chisq(BF) - Run 1-6 data set 
  
  // const double MaxHisto = 20.;
  // const double MinHisto = -20.;
  // const double NBinHisto = 200.;

  const double MaxHisto = 10.;
  const double MinHisto = -10.;
  const double NBinHisto = 200.;
  const double NBinHistoMany = 2000000.;
  const double newbinwidth = 10.;

  char name[200];

  TFile *filein = new TFile(infilename,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  
  vector<double> vec_chisq_nh_ih;
  vector<double> vec_chisq_fc;
  
  vector<double> vec_like_penalty;
  TH1D *hdeltachisq_nh_ih_posterior = new TH1D("hdeltachisq_nh_ih_posterior","hdeltachisq_nh_ih_posterior",NBinHistoMany,MinHisto,MaxHisto);
  //TH2D *hdeltachisq_nh_ih_Vs_chisqweight_posterior = new TH2D("hdeltachisq_nh_ih_Vs_chisqweight_posterior","hdeltachisq_nh_ih_Vs_chisqweight_posterior",NBinHisto,MinHisto,MaxHisto,NBinHisto,250,350);
  TH2D *hdeltachisq_nh_ih_Vs_chisqweight_posterior = new TH2D("hdeltachisq_nh_ih_Vs_chisqweight_posterior","hdeltachisq_nh_ih_Vs_chisqweight_posterior",NBinHisto,MinHisto,MaxHisto,NBinHisto,0,30);
  TH1D *hWeight_posterior = new TH1D("hWeight_posterior","hWeight_posterior",300,200,500);
  //TH2D *hbestfit_Vs_chisqweight_posterior = new TH2D("hbestfit_Vs_chisqweight_posterior","hbestfit_Vs_chisqweight_posterior",2,0,2,NBinHisto,250,350);
  TH2D *hbestfit_Vs_chisqweight_posterior = new TH2D("hbestfit_Vs_chisqweight_posterior","hbestfit_Vs_chisqweight_posterior",2,0,2,NBinHisto,0,30);
  TH1D *hbestfit_posterior_1D = new TH1D("hbestfit_posterior_1D","hbestfit_posterior_1D",2,0,2);

  // Tree variables. Set globals because used by both trees
  bool IsFit2D_nh=false;
  bool IsFit2D_ih=false;
  int run=-999;
  string *tag;
  int nexpts=-999;
  int ToyID=-999;
  double truepar0=-999;
  double truepar1=-999;    
  string *trueMH;
  string *trueparname0;
  string *trueparname1;
  string *MH_glob;
  double deltachisq_fc_glob=conv::kALargeNum;  
  double deltachisq_nh_ih_glob=conv::kALargeNum;  
  double chisqmin_glob=conv::kALargeNum;
  double chisqmin_nh=conv::kALargeNum;
  double chisqmin_ih=conv::kALargeNum;
  TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_sample[conv::kns] = {};
  TH2D *hErecoExp_1Rmu = new TH2D("hErecoExp2D_1Rmu","hErecoExp2D_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_1Re = new TH2D("hErecoExp2D_1Re","hErecoExp2D_1Re",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Rmu = new TH2D("hErecoExp2D_RHC1Rmu","hErecoExp2D_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Re = new TH2D("hErecoExp2D_RHC1Re","hErecoExp2D_RHC1Re",1,0,1,1,0,1); 
 
  treein->SetBranchAddress("ToyID",&ToyID);
  treein->SetBranchAddress("Job",&run);
  treein->SetBranchAddress("Tag",&tag);
  treein->SetBranchAddress("nexpts",&nexpts);
  treein->SetBranchAddress("truepar0", &truepar0);
  treein->SetBranchAddress("truepar1", &truepar1);
  treein->SetBranchAddress("trueMH",   &trueMH);
  treein->SetBranchAddress("trueparname0",&trueparname0);
  treein->SetBranchAddress("trueparname1",&trueparname1);
  treein->SetBranchAddress("deltachisq0BF_glob",&deltachisq_fc_glob);
  treein->SetBranchAddress("deltachisq01_glob",&deltachisq_nh_ih_glob);
  treein->SetBranchAddress("MH_glob",   &MH_glob);
  treein->SetBranchAddress("chisqmin_glob",&chisqmin_glob);
  treein->SetBranchAddress("chisqmin_nh",&chisqmin_nh);
  treein->SetBranchAddress("chisqmin_ih",&chisqmin_ih);
  treein->SetBranchAddress("hErecoToy_1Rmu",     &hErecoToy_1Rmu     );
  treein->SetBranchAddress("hErecoToy_1Re",      &hErecoToy_1Re      );
  treein->SetBranchAddress("hErecoToy_RHC1Rmu",  &hErecoToy_RHC1Rmu  );
  treein->SetBranchAddress("hErecoToy_RHC1Re",   &hErecoToy_RHC1Re   );
  treein->SetBranchAddress("hErecoExp_1Rmu",     &hErecoExp_1Rmu     );
  treein->SetBranchAddress("hErecoExp_1Re",      &hErecoExp_1Re      );
  treein->SetBranchAddress("hErecoExp_RHC1Rmu",  &hErecoExp_RHC1Rmu  );
  treein->SetBranchAddress("hErecoExp_RHC1Re",   &hErecoExp_RHC1Re   );

  // Data files

  TH2D *h_data_1rmu;
  TH2D *h_data_1re;
  TH2D *h_data_rhc1rmu;
  TH2D *h_data_rhc1re;
  
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));
  SuperkObservation * skobs_fakedata = new SuperkObservation;
  TFile *infile;
  if(use_posterior_predictive){
    cout << endl;
    cout << "PERFORM A POSTERIOR PREDICTIVE ANALYSIS:" << endl;
    cout << " - use root file" << rootfile2 << endl;
    cout << endl;
    infile = new TFile(rootfile2.Data(),"READ");
    if(!infile->IsOpen()){
      cout << endl;
      cout << "The data root file doesn't exist: don't draw data spectra" << endl;
      cout << endl;
    }
    infile->Close();
    

    conv::SuperkSample_t s_1rmu    = util::SKSampleFromString("1Rmu"); 
    conv::SuperkSample_t s_1re     = util::SKSampleFromString("1Re"); 
    conv::SuperkSample_t s_rhc1rmu = util::SKSampleFromString("RHC1Rmu"); 
    conv::SuperkSample_t s_rhc1re  = util::SKSampleFromString("RHC1Re"); 
    
    //int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
    bool using_xml_file = false;
    int iexperiment = 0;
    skobs_fakedata->Read(rootfile2.Data(),iexperiment,using_xml_file);
    h_data_1rmu    = (TH2D*) skobs_fakedata->TwoDSpectrum(s_1rmu)   ->Clone("h_data_1rmu");
    h_data_1re     = (TH2D*) skobs_fakedata->TwoDSpectrum(s_1re)    ->Clone("h_data_1re");
    h_data_rhc1rmu = (TH2D*) skobs_fakedata->TwoDSpectrum(s_rhc1rmu)->Clone("h_data_rhc1rmu");
    h_data_rhc1re  = (TH2D*) skobs_fakedata->TwoDSpectrum(s_rhc1re) ->Clone("h_data_rhc1re"); 
  }
    
  cout << "# of entries 1Rmu = " << h_data_1rmu->Integral() << endl;
  cout << "# of entries 1Re = " << h_data_1re->Integral() << endl;
  cout << "# of entries RHC1Rmu = " << h_data_rhc1rmu->Integral() << endl;
  cout << "# of entries RHC1Re = " << h_data_rhc1re->Integral() << endl;

  // Declare histograms
  TH1D *hnevents_rhc1re = new TH1D("hnevents_rhc1re","hnevents_rhc1re",NBinHisto,0.,MaxHisto);  
  TH1D *hnevents_1re = new TH1D("hnevents_1re","hnevents_1re",NBinHisto,0.,MaxHisto);  
  TH1D *hnevents_rhc1rmu = new TH1D("hnevents_rhc1rmu","hnevents_rhc1rmu",NBinHisto,0.,MaxHisto);  
  TH1D *hnevents_1rmu = new TH1D("hnevents_1rmu","hnevents_1rmu",NBinHisto,0.,MaxHisto);  
  
  TH1D *hdeltachisq_nh_ih = new TH1D("hdeltachisq_nh_ih","hdeltachisq_nh_ih",NBinHisto,MinHisto,MaxHisto);
  TH1D *hdeltachisq_fc = new TH1D("hdeltachisq_fc","hdeltachisq_fc",NBinHisto,MinHisto,MaxHisto);
  TH1D *hchisq_glob = new TH1D("hchisq_glob","hchisq_glob",400,0,400);
  TH1D *hchisq_nh = new TH1D("hchisq_nh","hchisq_nh",400,0,400);
  TH1D *hchisq_ih = new TH1D("hchisq_ih","hchisq_ih",400,0,400);
  TH1D *hMH_glob = new TH1D("hMH_glob","hMH_glob",2,0.,2.);
  
  // loop over the tree  
  
  int Nentries = treein->GetEntries();
  //Nentries = 10000;
  
  for(int ientry=0;ientry<Nentries;ientry++){
  
    treein->GetEntry(ientry);

    if(!(ientry%500)) cout << "Entry " << ientry << endl;

    //

    // double nevents_rhc1re = (double)hErecoToy_RHC1Re->Integral();    
    // double nevents_1re = (double)hErecoToy_1Re->Integral();    
    // double nevents_rhc1rmu = (double)hErecoToy_RHC1Rmu->Integral();    
    // double nevents_1rmu = (double)hErecoToy_1Rmu->Integral();    
    // if(
    //    (nevents_rhc1re<0) ||
    //    (nevents_1re<0) ||
    //    (nevents_rhc1rmu<0) ||
    //    (nevents_1rmu<0) 
    //    ){
    //   cerr << endl;
    //   cerr << "ToyID: " << ToyID << endl;
    //   cerr << "# of events RHC1Re  = " << nevents_rhc1re  << " < 0" << endl;
    //   cerr << "# of events 1Re     = " << nevents_1re     << " < 0" << endl;
    //   cerr << "# of events RHC1Rmu = " << nevents_rhc1rmu << " < 0" << endl;
    //   cerr << "# of events 1Rmu    = " << nevents_1rmu    << " < 0" << endl;
    //   cerr << endl;
    //   exit(1);
    // }
    // hnevents_rhc1re->Fill(nevents_rhc1re);    
    // hnevents_1re->Fill(nevents_1re);    
    // hnevents_rhc1rmu->Fill(nevents_rhc1rmu);    
    // hnevents_1rmu->Fill(nevents_1rmu);

    if(deltachisq_nh_ih_glob<MinHisto || deltachisq_nh_ih_glob>MaxHisto){
      cout << "deltachisq_nh_ih_glob = " << deltachisq_nh_ih_glob << endl;
      //exit(1);
    }
    if(deltachisq_fc_glob<0.){
      cout << "deltachisq_fc_glob = " << deltachisq_fc_glob << endl;
      //exit(1);
    }
    
    //
    
    // cout << "chisqmin_nh = " << chisqmin_nh
    //  	 << " | chisqmin_ih = " << chisqmin_ih
    //  	 << " | chisqmin_glob = " << chisqmin_glob
    // 	 << " | deltachisq_nh_ih_glob = " << deltachisq_nh_ih_glob 
    // 	 << " | deltachisq_fc_glob = " << deltachisq_fc_glob 
    //  	 << " |  bestfit MH = " << *MH_glob 
    // 	 << endl;
    
    double my_deltachisq_nh_ih = chisqmin_nh - chisqmin_ih;
    if(my_deltachisq_nh_ih != deltachisq_nh_ih_glob){
      cout << "my_deltachisq_nh_ih = " << my_deltachisq_nh_ih << endl;
      cout << "deltachisq_nh_ih_glob = " << deltachisq_nh_ih_glob << endl;
      exit(1);
    }
    hdeltachisq_nh_ih->Fill(deltachisq_nh_ih_glob);
    hdeltachisq_fc->Fill(deltachisq_fc_glob);
    
    // Minimized chisq distribution
    hchisq_glob->Fill(chisqmin_glob);
    hchisq_nh->Fill(chisqmin_nh);
    hchisq_ih->Fill(chisqmin_ih);
          
    if(*MH_glob=="normal"){
      hMH_glob->Fill(0.1);
    }
    else if(*MH_glob=="inverted"){
      hMH_glob->Fill(1.1);
    }
    else{
      cerr << endl;
      cerr << "MH_glob = " << *MH_glob << endl;
      cerr << "Best-fit Mass Hierarchy has an incorrect value!" << endl;
      cerr << endl;
      exit(1);
    }

    // Fill vectors for pvalue calculations
    vec_chisq_nh_ih.push_back(deltachisq_nh_ih_glob);
    vec_chisq_fc.push_back(deltachisq_fc_glob);

    
    // Get Likelihood weights --> penalty term
    double toy_weight = 1.;
    double toy_chisq = 0;
    if(use_posterior_predictive){
      if(use_post_pred_fhc_mulike) toy_chisq += calc_chisq(h_data_1rmu,   hErecoExp_1Rmu);    // FHC 1Rmu
      if(use_post_pred_fhc_elike)  toy_chisq += calc_chisq(h_data_1re,    hErecoExp_1Re);     // FHC 1Re
      if(use_post_pred_rhc_mulike) toy_chisq += calc_chisq(h_data_rhc1rmu,hErecoExp_RHC1Rmu); // RHC 1Rmu
      if(use_post_pred_rhc_elike)  toy_chisq += calc_chisq(h_data_rhc1re, hErecoExp_RHC1Re);  // RHC 1Re
      toy_weight = TMath::Exp(-1*toy_chisq/2);
    }
    vec_like_penalty.push_back(toy_weight);

    hdeltachisq_nh_ih_posterior->Fill(deltachisq_nh_ih_glob,toy_weight);    
    hdeltachisq_nh_ih_Vs_chisqweight_posterior->Fill(deltachisq_nh_ih_glob,toy_chisq);
    hWeight_posterior->Fill(toy_weight);

    if(*MH_glob=="normal"){
        hbestfit_posterior_1D->Fill(0.1,toy_weight);
        hbestfit_Vs_chisqweight_posterior->Fill(0.1,toy_chisq);
    }
    else if(*MH_glob=="inverted"){
      hbestfit_posterior_1D->Fill(1.1,toy_weight);
      hbestfit_Vs_chisqweight_posterior->Fill(1.1,toy_chisq);
    }
    
    //cout << "toy_chisq = " << toy_chisq << " - toy_weight = " << toy_weight << endl;
  } // end loop over toys
  
  
  //
  // Check the Data deltachisq value, depending on the True Mass Hierarchy in the toys
  // If trueMH = MH_glob --> set deltachisq_fc_data to 0
  // Otherwise --> set deltachisq_fc_data to |deltachisq_nh_ih_data|
  //
  
  string databfMH = ""; // best-fit Mass Hierarchy from fit - Run 1-6 data
  if(deltachisq_nh_ih_data < 0) databfMH = "normal";

  cout << endl;
  cout << "//////////// WARNING ////////////" << endl;
  if(*trueMH == databfMH){
    cout << endl;
    cout << "In Data set:" << endl;
    cout << "True Mass Hierarchy = Best-Fit Data Mass Hierarchy" << endl;
    cout << "Set deltachisq_fc_data (Feldman-Cousins) to 0" << endl;
    cout << endl;    
    deltachisq_fc_data = 0.;
  }
  else{
    cout << endl;
    cout << "In Data set:" << endl;
    cout << "True Mass Hierarchy != Best-Fit Data Mass Hierarchy" << endl;
    cout << "Set deltachisq_fc_data (Feldman-Cousins) to |chisq(NH)-chisq(IH)|" << endl;
    cout << endl;
    deltachisq_fc_data = fabs(deltachisq_nh_ih_data);
  }
  cout << "////////////////////////////////" << endl;
  cout << endl;
  
  
  
  
  // Compute p-values (NH-IH)
  
  double Ntoys = (double)Nentries;
  double pvalue_nh_ih_right = GetPvalue_rightside(vec_chisq_nh_ih,Ntoys,deltachisq_nh_ih_data);
  double pvalue_nh_ih_left  = GetPvalue_leftside(vec_chisq_nh_ih,Ntoys,deltachisq_nh_ih_data);

  double pvalue_nh_ih_glob  = -999;  
  if(*trueMH == "normal"){
    pvalue_nh_ih_glob = pvalue_nh_ih_right;
  }
  else if(*trueMH == "inverted"){    
    pvalue_nh_ih_glob = pvalue_nh_ih_left;    
  }
  else{
    cout << endl;
    cout << "The true Mass Hierarchy is not known: " << *trueMH << endl;
    cout << endl;
  }

  // Calculate the quantile (p-value = gauss area of both tails --> divide by 2)
  TF1 f2("f2","ROOT::Math::gaussian_quantile(1.-(x/2.),1.) ",0.,1.);
  double zscore_nh_ih_glob = f2.Eval(pvalue_nh_ih_glob); 


  // Compute p-values (a-la-Feldman-Cousins)

  double pvalue_fc_glob = GetPvalue_rightside(vec_chisq_fc,Ntoys,deltachisq_fc_data);
  
  double FCCritVal_1sigma_fc = CalcDLcrit(vec_chisq_fc,CLfrac_1sigma*100.,0.); 
  double FCCritVal_90CL_fc   = CalcDLcrit(vec_chisq_fc,CLfrac_90CL*100.,0.); 
  double FCCritVal_95CL_fc   = CalcDLcrit(vec_chisq_fc,CLfrac_95CL*100.,0.); 
  double FCCritVal_2sigma_fc = CalcDLcrit(vec_chisq_fc,CLfrac_2sigma*100.,0.); 

  double zscore_fc_glob = f2.Eval(pvalue_fc_glob); 
  

  // Compute p-value (use likelihood weights from data sets)
  
  //double pvalue_nh_ih_posterior_left = GetPvalue_Histo_leftside(hdeltachisq_nh_ih_posterior,deltachisq_nh_ih_data);
  //double pvalue_nh_ih_posterior_right = GetPvalue_Histo_rightside(hdeltachisq_nh_ih_posterior,deltachisq_nh_ih_data);
  double pvalue_nh_ih_posterior_left = GetPvalue_posterior_leftside(vec_chisq_nh_ih,vec_like_penalty,deltachisq_nh_ih_data);
  double pvalue_nh_ih_posterior_right = GetPvalue_posterior_rightside(vec_chisq_nh_ih,vec_like_penalty,deltachisq_nh_ih_data);
  double pvalue_nh_ih_posterior_glob = -999;
  if(*trueMH == "normal"){
    pvalue_nh_ih_posterior_glob = pvalue_nh_ih_posterior_right;
  }
  else if(*trueMH == "inverted"){    
    pvalue_nh_ih_posterior_glob = pvalue_nh_ih_posterior_left;    
  }
  else{
    cout << endl;
    cout << "The true Mass Hierarchy is not known: " << *trueMH << endl;
    cout << endl;
  }
  
  double zscore_nh_ih_posterior_glob = f2.Eval(pvalue_nh_ih_posterior_glob); 

  double pvalue_nh_ih_posterior_err = get_pvalue_error_vector(vec_chisq_nh_ih,vec_like_penalty,deltachisq_nh_ih_data,pvalue_nh_ih_posterior_glob); 
  double zscore_nh_ih_posterior_perr = f2.Eval(pvalue_nh_ih_posterior_glob - pvalue_nh_ih_posterior_err); // smaller p-value gives larger significance 
  zscore_nh_ih_posterior_perr = fabs(zscore_nh_ih_posterior_glob - zscore_nh_ih_posterior_perr); 
  double zscore_nh_ih_posterior_merr = f2.Eval(pvalue_nh_ih_posterior_glob + pvalue_nh_ih_posterior_err); // larger p-value gives lower significance 
  zscore_nh_ih_posterior_merr = fabs(zscore_nh_ih_posterior_glob - zscore_nh_ih_posterior_merr); 


  // Draw histograms
  
  // treein->GetEntry(493);
  // //treein->GetEntry(4739);
  // TCanvas *cErecoToy_GetEntry = new TCanvas("cErecoToy_GetEntry","cErecoToy_GetEntry");
  // hErecoToy_RHC1Re->DrawCopy();

  TCanvas *cMH = new TCanvas("cMH","cMH");
  hMH_glob->GetYaxis()->SetTitle("Number of entries");
  hMH_glob->GetXaxis()->SetTitle("Mass Hierarchy best-fit");
  hMH_glob->DrawCopy();
  hMH_glob->Draw("same");

  TCanvas *cchisq_nh = new TCanvas("cchisq_nh","cchisq_nh");
  hchisq_nh->GetYaxis()->SetTitle("Number of entries");
  hchisq_nh->GetXaxis()->SetTitle("#chi^{2} (Normal Hierarchy)");
  hchisq_nh->DrawCopy();
  hchisq_nh->Draw("same");
  //cchisq_nh->SetLogy();

  TCanvas *cchisq_ih = new TCanvas("cchisq_ih","cchisq_ih");
  hchisq_ih->GetYaxis()->SetTitle("Number of entries");
  hchisq_ih->GetXaxis()->SetTitle("#chi^{2} (Inverted Hierarchy)");
  hchisq_ih->DrawCopy();
  hchisq_ih->Draw("same");
  //cchisq_ih->SetLogy();

  TCanvas *cchisq_glob = new TCanvas("cchisq_glob","cchisq_glob");
  hchisq_glob->GetYaxis()->SetTitle("Number of entries");
  hchisq_glob->GetXaxis()->SetTitle("#chi^{2}_{glob}");
  hchisq_glob->DrawCopy();
  hchisq_glob->Draw("same");
  //cchisq_glob->SetLogy();

  TCanvas *cPvalue_NH_IH = new TCanvas("cPvalue_NH_IH","cPvalue_NH_IH");
  //cPvalue->SetLogy();
  hdeltachisq_nh_ih->GetYaxis()->SetTitle("Number of entries");
  hdeltachisq_nh_ih->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(NH) - #chi^{2}(IH)");
  hdeltachisq_nh_ih->DrawCopy();  
  TLine *ldeltachisq_nh_ih = new TLine(deltachisq_nh_ih_data,0.,deltachisq_nh_ih_data,hdeltachisq_nh_ih->GetMaximum());
  ldeltachisq_nh_ih->SetLineColor(kRed);
  ldeltachisq_nh_ih->SetLineWidth(2);
  ldeltachisq_nh_ih->Draw();
  TLatex *latPvalue_nh_ih = new TLatex(hdeltachisq_nh_ih->GetXaxis()->GetBinLowEdge(5),hdeltachisq_nh_ih->GetMaximum()*0.7, TString::Format("P-value = %1.3f",pvalue_nh_ih_glob).Data() );
  latPvalue_nh_ih->Draw("same");
  latPvalue_nh_ih->SetTextSize(0.06);
  TLine *lPvalue_nh_ih = new TLine(pvalue_nh_ih_glob,pvalue_nh_ih_glob,pvalue_nh_ih_glob,pvalue_nh_ih_glob);  
  TLatex *latZscore_nh_ih = new TLatex(hdeltachisq_nh_ih->GetXaxis()->GetBinLowEdge(5),hdeltachisq_nh_ih->GetMaximum()*0.5, TString::Format("Z = %1.2f#sigma",zscore_nh_ih_glob).Data() );
  latZscore_nh_ih->Draw("same");
  latZscore_nh_ih->SetTextSize(0.06);
  TLine *lZscore_nh_ih = new TLine(zscore_nh_ih_glob,zscore_nh_ih_glob,zscore_nh_ih_glob,zscore_nh_ih_glob);
  TLegend *legPvalue_nh_ih = new TLegend(0.5,0.6,0.8,0.8);
  legPvalue_nh_ih->SetFillColor(0);
  legPvalue_nh_ih->SetLineColor(0);
  legPvalue_nh_ih->SetTextFont(132);
  legPvalue_nh_ih->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_nh_ih_data); 
  legPvalue_nh_ih->AddEntry(ldeltachisq_nh_ih,name,"l");
  legPvalue_nh_ih->Draw();  

  TCanvas *cPvalue_FC = new TCanvas("cPvalue_FC","cPvalue_FC");
  //cPvalue->SetLogy();
  hdeltachisq_fc->GetYaxis()->SetTitle("Number of entries");
  hdeltachisq_fc->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}_{True} - #chi^{2}_{Best-Fit}");
  hdeltachisq_fc->DrawCopy(); 
  TLine *ldeltachisq_fc = new TLine(deltachisq_fc_data,0.,deltachisq_fc_data,hdeltachisq_fc->GetMaximum());
  ldeltachisq_fc->SetLineColor(kRed);
  ldeltachisq_fc->SetLineWidth(2);
  ldeltachisq_fc->Draw();
  TLatex *latPvalue_fc = new TLatex(hdeltachisq_fc->GetXaxis()->GetBinLowEdge(5),hdeltachisq_fc->GetMaximum()*0.7, TString::Format("P-value = %1.3f",pvalue_fc_glob).Data() );
  latPvalue_fc->Draw("same");
  latPvalue_fc->SetTextSize(0.06);
  TLine *lPvalue_fc = new TLine(pvalue_fc_glob,pvalue_fc_glob,pvalue_fc_glob,pvalue_fc_glob);  
  TLatex *latZscore_fc = new TLatex(hdeltachisq_fc->GetXaxis()->GetBinLowEdge(5),hdeltachisq_fc->GetMaximum()*0.5, TString::Format("Z = %1.2f#sigma",zscore_fc_glob).Data() );
  latZscore_fc->Draw("same");
  latZscore_fc->SetTextSize(0.06);
  TLine *lZscore_fc = new TLine(zscore_fc_glob,zscore_fc_glob,zscore_fc_glob,zscore_fc_glob);
  TLegend *legPvalue_fc = new TLegend(0.5,0.6,0.8,0.8);
  legPvalue_fc->SetFillColor(0);
  legPvalue_fc->SetLineColor(0);
  legPvalue_fc->SetTextFont(132);
  legPvalue_fc->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_fc_data); 
  legPvalue_fc->AddEntry(ldeltachisq_fc,name,"l");
  legPvalue_fc->Draw();    
  
  TLine *lFCCritVal_1sigma_fc = new TLine(FCCritVal_1sigma_fc,0.,FCCritVal_1sigma_fc,hdeltachisq_fc->GetMaximum());  
  TLine *lFCCritVal_90CL_fc = new TLine(FCCritVal_90CL_fc,0.,FCCritVal_90CL_fc,hdeltachisq_fc->GetMaximum());  
  TLine *lFCCritVal_2sigma_fc = new TLine(FCCritVal_2sigma_fc,0.,FCCritVal_2sigma_fc,hdeltachisq_fc->GetMaximum());  

  lFCCritVal_1sigma_fc->SetLineWidth(2);
  lFCCritVal_1sigma_fc->SetLineColor(kBlue);
  lFCCritVal_1sigma_fc->Draw();
  lFCCritVal_90CL_fc->SetLineWidth(2);
  lFCCritVal_90CL_fc->SetLineColor(kGreen);
  lFCCritVal_90CL_fc->Draw();
  lFCCritVal_2sigma_fc->SetLineWidth(2);
  lFCCritVal_2sigma_fc->SetLineColor(kRed+3);
  lFCCritVal_2sigma_fc->Draw();
  
  

  // Draw posterior predictive p-value

  hdeltachisq_nh_ih_posterior->Scale(1. / hdeltachisq_nh_ih_posterior->Integral());
  double mergebin = NBinHistoMany/(MaxHisto-MinHisto) / newbinwidth; 
  hdeltachisq_nh_ih_posterior->Rebin(mergebin);

  TCanvas *cPvalue_NH_IH_posterior = new TCanvas("cPvalue_NH_IH_posterior","cPvalue_NH_IH_posterior");
  //cPvalue->SetLogy();
  hdeltachisq_nh_ih_posterior->GetYaxis()->SetTitle("p.d.f.");
  hdeltachisq_nh_ih_posterior->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(NH) - #chi^{2}(IH)");
  hdeltachisq_nh_ih_posterior->DrawCopy();  
  ldeltachisq_nh_ih->Draw();
  TLatex *latPvalue_nh_ih_posterior = new TLatex(hdeltachisq_nh_ih_posterior->GetXaxis()->GetBinLowEdge(5),hdeltachisq_nh_ih_posterior->GetMaximum()*0.7, TString::Format("P-value = %1.3f",pvalue_nh_ih_posterior_glob).Data() );
  latPvalue_nh_ih->Draw("same");
  latPvalue_nh_ih->SetTextSize(0.06);
  TLine *lPvalue_nh_ih_posterior = new TLine(pvalue_nh_ih_posterior_glob,pvalue_nh_ih_posterior_glob,pvalue_nh_ih_posterior_glob,pvalue_nh_ih_posterior_glob);  
  TLatex *latZscore_nh_ih_posterior = new TLatex(hdeltachisq_nh_ih_posterior->GetXaxis()->GetBinLowEdge(5),hdeltachisq_nh_ih_posterior->GetMaximum()*0.5, TString::Format("Z = %1.2f#sigma",zscore_nh_ih_posterior_glob).Data() );
  latZscore_nh_ih_posterior->Draw("same");
  latZscore_nh_ih_posterior->SetTextSize(0.06);
  TLine *lZscore_nh_ih_posterior = new TLine(zscore_nh_ih_posterior_glob,zscore_nh_ih_posterior_glob,zscore_nh_ih_posterior_glob,zscore_nh_ih_posterior_glob);

  cout << "zscore_nh_ih_posterior_glob = " << zscore_nh_ih_posterior_glob << endl;
  
  TLegend *legPvalue_nh_ih_posterior = new TLegend(0.5,0.6,0.8,0.8);
  legPvalue_nh_ih_posterior->SetFillColor(0);
  legPvalue_nh_ih_posterior->SetLineColor(0);
  legPvalue_nh_ih_posterior->SetTextFont(132);
  legPvalue_nh_ih_posterior->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_nh_ih_data); 
  legPvalue_nh_ih_posterior->AddEntry(ldeltachisq_nh_ih,name,"l");
  legPvalue_nh_ih_posterior->Draw();  
  
  TCanvas *cWeight_posterior = new TCanvas("cWeight_posterior","cWeight_posterior");
  hWeight_posterior->GetXaxis()->SetTitle("L (data,fake exp)");
  hWeight_posterior->GetYaxis()->SetTitle("# of entries");
  hWeight_posterior->Draw();

  TCanvas *cdeltachisq_nh_ih_Vs_chisqweight_posterior = new TCanvas("cdeltachisq_nh_ih_Vs_chisqweight_posterior","cdeltachisq_nh_ih_Vs_chisqweight_posterior");
  hdeltachisq_nh_ih_Vs_chisqweight_posterior->GetXaxis()->SetTitle("test statistic");
  hdeltachisq_nh_ih_Vs_chisqweight_posterior->GetYaxis()->SetTitle("#chi^{2} (data,fake exp)");
  hdeltachisq_nh_ih_Vs_chisqweight_posterior->Draw("colz");

  TCanvas *cbestfit_posterior_1D = new TCanvas("cbestfit_posterior_1D","cbestfit_posterior_1D");
  hbestfit_posterior_1D->GetXaxis()->SetTitle("Mass Hierarchy posterior best-fit");
  hbestfit_posterior_1D->Draw();

  TCanvas *cbestfit_Vs_chisqweight_posterior = new TCanvas("cbestfit_Vs_chisqweight_posterior","cbestfit_Vs_chisqweight_posterior");
  hbestfit_Vs_chisqweight_posterior->GetXaxis()->SetTitle("Mass Hierarchy best-fit");
  hbestfit_Vs_chisqweight_posterior->GetYaxis()->SetTitle("#chi^{2} (data,fake exp)");
  hbestfit_Vs_chisqweight_posterior->Draw("colz");

  //


  TLine *lPvalue_nh_ih_posterior_err  = new TLine(pvalue_nh_ih_posterior_err,1.,1.,1.);  
  TLine *lZscore_nh_ih_posterior_perr = new TLine(zscore_nh_ih_posterior_perr,1.,1.,1.);
  TLine *lZscore_nh_ih_posterior_merr = new TLine(zscore_nh_ih_posterior_merr,1.,1.,1.);



  char outname[200];
  sprintf(outname,"%s.root",mytag.Data());
  TFile *outfile = new TFile(outname,"RECREATE");

  hdeltachisq_nh_ih_posterior ->Write("hdeltachisq_Toys_nh_ih_post");
  lPvalue_nh_ih_posterior     ->Write("Pvalue_nh_ih_posterior");
  lZscore_nh_ih_posterior     ->Write("Zscore_nh_ih_posterior");

  cout << "tline = " << lZscore_nh_ih_posterior->GetX1() << endl;

  lPvalue_nh_ih_posterior_err ->Write("lPvalue_nh_ih_posterior_err");
  lZscore_nh_ih_posterior_perr->Write("lzscore_nh_ih_posterior_perr");
  lZscore_nh_ih_posterior_merr->Write("lzscore_nh_ih_posterior_merr");

  hdeltachisq_nh_ih           ->Write("hdeltachisq_Toys_nh_ih");
  ldeltachisq_nh_ih           ->Write("ldeltachisq_Data_nh_ih");
  lPvalue_nh_ih               ->Write("Pvalue_nh_ih");
  lZscore_nh_ih               ->Write("Zscore_nh_ih");
  hdeltachisq_fc              ->Write("hdeltachisq_Toys_fc");
  ldeltachisq_fc              ->Write("ldeltachisq_Data_fc");
  lPvalue_fc                  ->Write("Pvalue_fc");
  lZscore_fc                  ->Write("Zscore_fc");
  outfile->Close();

  // Print results
  
  cout << endl;
  cout << std::setprecision(5);
  cout << "Fit of: " << *trueparname0 << " Vs " << *trueparname1 << endl;
  cout << "- Toy experiments:" << endl;
  cout << "   Toy experiments are generated without marginalizing the Mass Hierarchy" << endl; 
  cout << "   True Mass Hierarchy: " << *trueMH << endl;
  cout << "   Total # of toy experiments = " << Ntoys << endl;
  cout << "   truepar0 = " << truepar0 << endl;
  cout << "   truepar1 = " << truepar1 << endl;
  cout << endl;
  cout << "- Data: " << endl;
  cout << "   Best-fit Mass Hierarchy: " << databfMH << endl;
  cout << "   DeltaChisq = Chisq(NH) - Chisq(IH) = " << deltachisq_nh_ih_data << endl;
  cout << "   DeltaChisq = Chisq(True) - Chisq(Best-Fit) = " << deltachisq_fc_data << endl;
  cout << endl;
  cout << "- Calculate right sided p-values:" << endl;
  cout << "   p-value = " << pvalue_nh_ih_right << endl;
  cout << "- Calculate left sided p-values:" << endl;
  cout << "   p-value = " << pvalue_nh_ih_left << endl;
  cout << "- The global p-value is calculated as ";
  if     (*trueMH == "normal")   cout << "right-side" << endl;
  else if(*trueMH == "inverted") cout << "left-side" << endl;
  cout << "- DeltaChisq = Chisq(NH) - Chisq(IH):" << endl;
  cout << "   p-value = " << pvalue_nh_ih_glob << endl;
  cout << "   z-score = " << zscore_nh_ih_glob << " sigma" << endl;
  cout << "- DeltaChisq = Chisq(True) - Chisq(Best-Fit):" << endl;
  cout << "   p-value = " << pvalue_fc_glob << endl;
  cout << "   z-score = " << zscore_fc_glob << " sigma" << endl;
  cout << endl;
  cout << " - Calculate posterior predictive p-value:" << endl;
  cout << "   Use the data samples:" << endl;
  if(use_post_pred_fhc_mulike) cout << "     - FHC mu-like" << endl; 
  if(use_post_pred_fhc_elike) cout  << "     - FHC e-like"  << endl; 
  if(use_post_pred_fhc_mulike) cout << "     - RHC mu-like" << endl; 
  if(use_post_pred_rhc_elike) cout  << "     - RHC e-like"  << endl; 
  cout << "   DeltaChisq = Chisq(NH) - Chisq(IH):" << endl;
  cout << "   p-value = " << pvalue_nh_ih_posterior_glob 
       << " +- " << pvalue_nh_ih_posterior_err
       << endl;
  cout << "   z-score = " << zscore_nh_ih_posterior_glob
       << " + " << zscore_nh_ih_posterior_perr
       << " - " << zscore_nh_ih_posterior_merr
       << " sigma" << endl;  
  cout << endl;
  cout << "Critical values:" << endl;
  cout << " - 1sigma = " << FCCritVal_1sigma_fc << endl;
  cout << " - 90% CL = " << FCCritVal_90CL_fc << endl;
  cout << " - 2sigma = " << FCCritVal_2sigma_fc << endl;
  cout << endl;
  // cout << "Use my function:" << endl;
  // cout << CalcDLcrit_posterior(vec_chisq_nh_ih,vec_like_penalty,CLfrac_1sigma*100.,0.) << endl;
  // cout << CalcDLcrit_posterior(vec_chisq_nh_ih,vec_like_penalty,CLfrac_90CL*100.,0.) << endl;
  // cout << CalcDLcrit_posterior(vec_chisq_nh_ih,vec_like_penalty,CLfrac_2sigma*100.,0.) << endl;
  // cout << "Use histogram:" << endl;
  // cout << CalcDLcrit(hdeltachisq_nh_ih_posterior,CLfrac_1sigma) << endl;
  // cout << CalcDLcrit(hdeltachisq_nh_ih_posterior,CLfrac_90CL) << endl;
  // cout << CalcDLcrit(hdeltachisq_nh_ih_posterior,CLfrac_2sigma) << endl;
  // cout << endl;
  return;
}



double GetPvalue_rightside(vector<double> vec_toys,int ntoys,double data){
  double ntot_toys = ntoys;
  if(ntot_toys != vec_toys.size()){    
    cout << endl;
    cout << "Warning:" << endl;
    cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
    cout << "Input total number of toys = " << ntot_toys << endl; 
    ntot_toys = vec_toys.size();
    cout << "Set total number of toys to vector size: " << ntot_toys << endl;
    cout << endl;
  }
  if(ntot_toys<=0.){
    cerr << endl;
    cerr << "ntot_toys = " << ntot_toys << endl;
    cerr << endl;
    exit(1);
  }

  cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
  cout << "Input total number of toys = " << ntot_toys << endl; 
  ntot_toys = vec_toys.size();
  cout << "Set total number of toys to vector size: " << ntot_toys << endl;
  cout << endl;

  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it >= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}

double GetPvalue_leftside(vector<double> vec_toys,int ntoys,double data){
  double ntot_toys = ntoys;
  if(ntot_toys != vec_toys.size()){    
    cout << endl;
    cout << "Warning:" << endl;
    cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
    cout << "Input total number of toys = " << ntot_toys << endl; 
    ntot_toys = vec_toys.size();
    cout << "Set total number of toys to vector size: " << ntot_toys << endl;
    cout << endl;
  }
  if(ntot_toys<=0.){
    cerr << endl;
    cerr << "ntot_toys = " << ntot_toys << endl;
    cerr << endl;
    exit(1);
  }

  cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
  cout << "Input total number of toys = " << ntot_toys << endl; 
  ntot_toys = vec_toys.size();
  cout << "Set total number of toys to vector size: " << ntot_toys << endl;
  cout << endl;

  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it <= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}

double GetPvalue_Histo_rightside(TH1D *histo,double data){
  int bin_data = histo->FindBin(data);
  double pvalue = histo->Integral(bin_data,histo->GetNbinsX());
  pvalue /= histo->Integral();
  
  return pvalue;
}

double GetPvalue_Histo_leftside(TH1D *histo,double data){
  int bin_data = histo->FindBin(data);
  double pvalue = histo->Integral(1,bin_data);
  pvalue /= histo->Integral();
  
  return pvalue;
}

double CalcDLcrit(vector<double> dchisq, double conflev, double sigmacl){
  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector   
  vector<double>::iterator dchisqIter;
  dchisqIter = dchisq.begin() + (int)(0.01*conflev*dchisq.size() + sigmacl);
  //cout << 0.01*conflev*dchisq.size() << endl;
  double critval = (*dchisqIter);
  return critval;
}

double calc_chisq(TH2D* hobs, TH2D* hexp){
  int nbins_x = hobs->GetNbinsX();
  int nbins_y = hobs->GetNbinsY();
  double chisq = 0;
  
  for(int i=1;i<= nbins_x;i++)
    {
      for(int j=1;j<= nbins_y;j++)
	{
	  double obs = hobs->GetBinContent(i,j);
	  double exp = hexp->GetBinContent(i,j);
	  chisq += valor::util::ChisqPoisson(obs,exp);
	}
    } 
  return chisq;
}

double CalcDLcrit_posterior(vector<double> dchisq, vector<double> weight, double conflev, double sigmacl){

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
  index = TMath::Sort(nelem,new_dchisq,index,0);
  
  // get the sum
  double tot = 0.;
  for(int i=0;i<nelem;i++){
    tot += new_dchisq[index[i]];
  }
      
  double tot_weight = 0.;
  for(int i=0;i<nelem;i++) tot_weight += new_weight[index[i]]; 
  
  double sum = 0.;
  double critval = -999.;
  for(int i=0;i<nelem;i++){
    sum += new_weight[index[i]] / tot_weight; 
    
    critval = new_dchisq[index[i]];
    //cout << critval << " ";
    //cout << sum << " ";
    
    if(sum > 0.01 * conflev){
      critval = new_dchisq[index[i]]; 
      //cout << "sum = " << sum + new_weight[index[i]] / tot_weight << endl;
      //cout << "sum = " << sum << endl;
      break;
    }
  }
      
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
      
      // chi2_bin = ibin+1;  // histgram rule for the # of bins!!!
      // chi2_val = chi2_bin*h_delta_true->GetBinWidth(1) + MinDL;
      // cout << "chi2_val = " << chi2_val << endl;

      // chi2_bin = ibin;  // histogram rule for the # of bins!!!
      // chi2_val = chi2_bin*h_delta_true->GetBinWidth(1)  
      // 	+ h_delta_true->GetXaxis()->GetBinLowEdge(1);

      chi2_val = h_delta_true->GetBinCenter(ibin+1);
 
      break;
    }
  }
  return chi2_val;
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


double GetPvalue_posterior_leftside(vector<double> vec_toys,vector<double> weight,double data){
  double ntot_toys = vec_toys.size();  
  double tot = 0.; 
  double nfrac=0;
  for(int i=0;i<vec_toys.size();i++){
    tot += weight[i];
    double cont = vec_toys[i];
    if(cont <= data) nfrac += weight[i];
  }
  double pvalue =  nfrac / tot;

  return pvalue;
}

double get_pvalue_error_vector(vector<double> vec_dchisq,vector<double> vec_weight, double data_critval, double pvalue) // takes the p-value not CL!!!
{
  double error = 0;
  double total_right=0;
  double total_left=0;
  
  std::vector<double> _right_tail_weight;
  std::vector<double> _left_tail_weight;

  for(int i=0;i<vec_dchisq.size();i++){
    
    double weight = vec_weight[i]; 

    if(vec_dchisq[i]<data_critval){
      _left_tail_weight.push_back(weight);
    }
    else{
      _right_tail_weight.push_back(weight);
    }
  }
  
  for(int i=0;i<_right_tail_weight.size();i++){
    error+= TMath::Power(_right_tail_weight[i]*(1-pvalue),2);
    total_right+= _right_tail_weight[i];
  }
  
  for(int i=0;i<_left_tail_weight.size();i++){
    error+= TMath::Power(_left_tail_weight[i]*(pvalue),2);
    total_left+= _left_tail_weight[i];
  }
  
  error = error / (TMath::Power(total_left + total_right ,2));
  return TMath::Sqrt(error);  
}
