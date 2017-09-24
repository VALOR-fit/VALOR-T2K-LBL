
// ------------------------------------------------------------------------------------------------
//
// Read output from run_fit_marginalization.cxx (deltachisq distribution
// and run_fc_marginalization.cxx (critical values)
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------

#include "Build.h"
#include "AnalysisOptions.h"
#include <vector>
#include "TString.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

const double CLfrac_1sigma = 0.6827;  // 1sigma% CL
const double CLfrac_90CL   = 0.9;  // 90% CL
const double CLfrac_95CL   = 0.95;  // 95% CL
const double CLfrac_2sigma = 0.9545;  // 2sigma% CL
const double CLfrac_3sigma = 0.9973;  // 3sigma% CL

void fc_critical_values(){

  cout << "Making the F&C Critical Values" << endl;
  
  const int Nmh = 2;
  const int NFCpts = 8; // # of points of Feldman-Cousins
  bool fixed_hierarchy = 0;
  
  // RUN 1-8
  // Root path to toys location, for convenience specifying path in "infilename_fc"
  string path_fc_nh = "/home/phrdqd/VaLOR/analysis/tn/final/";
  string path_fc_ih = "/home/phrdqd/VaLOR/analysis/tn/final/";
    
  TH1D* h1[NFCpts][Nmh] = {0};
  TH1D* h1_bf[NFCpts][Nmh] = {0};
  
  char* hierarch[2] = {"nh","ih"};
  
  TH1D* h_draw_68 = new TH1D("h1","",1000,0,25);
  TH1D* h_draw_90 = new TH1D("h2","",1000,0,25);;
  TH1D* h_draw_2sig = new TH1D("h3","",1000,0,25);;
  TH1D* h_draw_3sig = new TH1D("h4","",1000,0,25);;
  
  TGraph* h_coverage_3sigma_nh = new TGraph();
  TGraph* h_coverage_2sigma_nh = new TGraph();
  TGraph* h_coverage_1sigma_nh = new TGraph();
  TGraph* h_coverage_90pc_nh = new TGraph();
  
  TGraph* h_coverage_3sigma_ih = new TGraph();
  TGraph* h_coverage_2sigma_ih = new TGraph();
  TGraph* h_coverage_1sigma_ih = new TGraph();
  TGraph* h_coverage_90pc_ih = new TGraph();
  
  h_draw_68->SetLineColor(kRed);
  h_draw_90->SetLineColor(kBlue);
  h_draw_2sig->SetLineColor(kGreen+3);
  h_draw_3sig->SetLineColor(kOrange-3);
  h_draw_68->SetLineWidth(2);
  h_draw_90->SetLineWidth(2);
  h_draw_2sig->SetLineWidth(2);
  h_draw_3sig->SetLineWidth(2);
  
  ofstream * outfilep = new ofstream("FC_output.txt",ios_base::trunc);
  ofstream & outfile = *outfilep;
  stringstream streamer;
  streamer << "\\begin{table}[t]"   << endl;
  streamer << "\\centering"   << endl;
  streamer << "\\begin{tabular}{| c | c | c|  c |}" << endl;
  streamer << "\\cline{1-4}" << endl;
  streamer << "$\\delta_{CP}$ & "
  << "CL & "
  << "Critical Value & "
  << "Coverage (Gaus) \\\\ "
  << endl;
  streamer << "\\cline{1-4}" << endl;
  
  cout << "Reading file names with F&C Fits" << endl;
  
  const char * infilename_fc[NFCpts][Nmh] = {
    
    ///////////////////////
    //                   //
    //      Run 1-8      //
    //                   //
    ///////////////////////
    

    {"fitmarg_fc_dcp_normal_final_10k_0.0.root","fitmarg_fc_dcp_inverted_final_10k_0.0.root"},
    {"fitmarg_fc_dcp_normal_final_10k_-0.785.root","fitmarg_fc_dcp_inverted_final_11k_-0.785.root"},
    {"fitmarg_fc_dcp_normal_final_10k_0.785.root","fitmarg_fc_dcp_inverted_final_10k_0.785.root"},
    {"fitmarg_fc_dcp_normal_final_15k_-1.571.root","fitmarg_fc_dcp_inverted_final_10k_-1.571.root"},
    {"fitmarg_fc_dcp_normal_final_31k_1.571.root","fitmarg_fc_dcp_inverted_final_14k_1.571.root"},
    {"fitmarg_fc_dcp_normal_final_10k_-2.356.root","fitmarg_fc_dcp_inverted_final_10k_-2.356.root"},
    {"fitmarg_fc_dcp_normal_final_10k_2.356.root","fitmarg_fc_dcp_inverted_final_10k_2.356.root"},    
    {"fitmarg_fc_dcp_normal_final_10k_3.142.root","fitmarg_fc_dcp_inverted_final_10k_3.142.root"}
  };


  
  const int nptcpy = 1;
  double copypt[nptcpy][2] = {
    { TMath::Pi(), -TMath::Pi() } // copy crit val of TMath::Pi() in -TMath::Pi()
  };
  
  TGraph *gFC90CL_nh = new TGraph();
  gFC90CL_nh->SetName("gFC90CL_nh");
  gFC90CL_nh->SetTitle("gFC90CL_nh");
  TGraph *gFC90CL_nh_p1 = new TGraph();
  gFC90CL_nh_p1->SetName("gFC90CL_nh_p1");
  gFC90CL_nh_p1->SetTitle("gFC90CL_nh_p1");
  TGraph *gFC90CL_nh_m1 = new TGraph();
  gFC90CL_nh_m1->SetName("gFC90CL_nh_m1");
  gFC90CL_nh_m1->SetTitle("gFC90CL_nh_m1");
  
  TGraph *gFC90CL_ih = new TGraph();
  gFC90CL_ih->SetName("gFC90CL_ih");
  gFC90CL_ih->SetTitle("gFC90CL_ih");
  TGraph *gFC90CL_ih_p1 = new TGraph();
  gFC90CL_ih_p1->SetName("gFC90CL_ih_p1");
  gFC90CL_ih_p1->SetTitle("gFC90CL_ih_p1");
  TGraph *gFC90CL_ih_m1 = new TGraph();
  gFC90CL_ih_m1->SetName("gFC90CL_ih_m1");
  gFC90CL_ih_m1->SetTitle("gFC90CL_ih_m1");

  TGraph *gFC1sigma_nh = new TGraph();
  gFC1sigma_nh->SetName("gFC1sigma_nh");
  gFC1sigma_nh->SetTitle("gFC1sigma_nh");
  TGraph *gFC1sigma_nh_p1 = new TGraph();
  gFC1sigma_nh_p1->SetName("gFC1sigma_nh_p1");
  gFC1sigma_nh_p1->SetTitle("gFC1sigma_nh_p1");
  TGraph *gFC1sigma_nh_m1 = new TGraph();
  gFC1sigma_nh_m1->SetName("gFC1sigma_nh_m1");
  gFC1sigma_nh_m1->SetTitle("gFC1sigma_nh_m1");
  
  TGraph *gFC1sigma_ih = new TGraph();
  gFC1sigma_ih->SetName("gFC1sigma_ih");
  gFC1sigma_ih->SetTitle("gFC1sigma_ih");
  TGraph *gFC1sigma_ih_p1 = new TGraph();
  gFC1sigma_ih_p1->SetName("gFC1sigma_ih_p1");
  gFC1sigma_ih_p1->SetTitle("gFC1sigma_ih_p1");
  TGraph *gFC1sigma_ih_m1 = new TGraph();
  gFC1sigma_ih_m1->SetName("gFC1sigma_ih_m1");
  gFC1sigma_ih_m1->SetTitle("gFC1sigma_ih_m1");
  
  TGraph *gFC2sigma_nh = new TGraph();
  gFC2sigma_nh->SetName("gFC2sigma_nh");
  gFC2sigma_nh->SetTitle("gFC2sigma_nh");
  TGraph *gFC2sigma_nh_p1 = new TGraph();
  gFC2sigma_nh_p1->SetName("gFC2sigma_nh_p1");
  gFC2sigma_nh_p1->SetTitle("gFC2sigma_nh_p1");
  TGraph *gFC2sigma_nh_m1 = new TGraph();
  gFC2sigma_nh_m1->SetName("gFC2sigma_nh_m1");
  gFC2sigma_nh_m1->SetTitle("gFC2sigma_nh_m1");
  
  TGraph *gFC2sigma_ih = new TGraph();
  gFC2sigma_ih->SetName("gFC2sigma_ih");
  gFC2sigma_ih->SetTitle("gFC2sigma_ih");
  TGraph *gFC2sigma_ih_p1 = new TGraph();
  gFC2sigma_ih_p1->SetName("gFC2sigma_ih_p1");
  gFC2sigma_ih_p1->SetTitle("gFC2sigma_ih_p1");
  TGraph *gFC2sigma_ih_m1 = new TGraph();
  gFC2sigma_ih_m1->SetName("gFC2sigma_ih_m1");
  gFC2sigma_ih_m1->SetTitle("gFC2sigma_ih_m1");
  
  TGraph *gFC3sigma_nh = new TGraph();
  gFC3sigma_nh->SetName("gFC3sigma_nh");
  gFC3sigma_nh->SetTitle("gFC3sigma_nh");
  TGraph *gFC3sigma_nh_p1 = new TGraph();
  gFC3sigma_nh_p1->SetName("gFC3sigma_nh_p1");
  gFC3sigma_nh_p1->SetTitle("gFC3sigma_nh_p1");
  TGraph *gFC3sigma_nh_m1 = new TGraph();
  gFC3sigma_nh_m1->SetName("gFC3sigma_nh_m1");
  gFC3sigma_nh_m1->SetTitle("gFC3sigma_nh_m1");
  
  TGraph *gFC3sigma_ih = new TGraph();
  gFC3sigma_ih->SetName("gFC3sigma_ih");
  gFC3sigma_ih->SetTitle("gFC3sigma_ih");
  TGraph *gFC3sigma_ih_p1 = new TGraph();
  gFC3sigma_ih_p1->SetName("gFC3sigma_ih_p1");
  gFC3sigma_ih_p1->SetTitle("gFC3sigma_ih_p1");
  TGraph *gFC3sigma_ih_m1 = new TGraph();
  gFC3sigma_ih_m1->SetName("gFC3sigma_ih_m1");
  gFC3sigma_ih_m1->SetTitle("gFC3sigma_ih_m1");
  
  TFile *filein_fc[NFCpts][Nmh];
  TTree *treein_fc[NFCpts][Nmh];
  
  // Tree variables
  double truepar0=-999;
  double truepar1=-999;
  string *trueMH;
  string *trueparname0;
  string *trueparname1;
  double chisqtrue0=conv::kALargeNum;
  double deltachisq0BF_glob=conv::kALargeNum;
  double chisqmin_nh=conv::kALargeNum;
  double chisqmin_ih=conv::kALargeNum;
  double bestfit0_nh=conv::kALargeNum;
  double bestfit0_ih=conv::kALargeNum;
  double bestfit0_glob=conv::kALargeNum;
  
  for(int ifc=0;ifc<NFCpts;ifc++){ // Loop over Feldman-Cousins files
    for(int imh=0;imh<Nmh;imh++){  // Loop over the Mass Hierarchy (0=NH, 1=IH)
      TString pathtot="";
      
      if(imh==0){
        //pathtot = TString::Format("%s%s",path_fc_nh.c_str(),infilename_fc[ifc][imh].c_str());
        pathtot = TString::Format("%s%s",path_fc_nh.c_str(),infilename_fc[ifc][imh]);
      }
      else if(imh==1){
        //pathtot = TString::Format("%s%s",path_fc_ih.c_str(),infilename_fc[ifc][imh].c_str());
        pathtot = TString::Format("%s%s",path_fc_ih.c_str(),infilename_fc[ifc][imh]);
      }
      else{
        cerr << endl;
        cerr << "You are using more than 2 mass hierarchy values!!!" << endl;
        cerr << endl;
        exit(1);
      }
      
      filein_fc[ifc][imh] = new TFile(pathtot.Data(),"READ");
      if(!filein_fc[ifc][imh]->IsOpen()){
        cerr << endl;
        cerr << "File " << pathtot.Data() << " is not open" << endl;
        cerr << endl;
      }
      treein_fc[ifc][imh] = (TTree*)filein_fc[ifc][imh]->Get("fitnt");
      if(!treein_fc[ifc][imh]){
        cerr << endl;
        cerr << "Tree object is empty" << endl;
        cerr << endl;
      }
      
      cout << endl;
      
      cout << "MH File: " << imh << endl;       
      cout << "FC File: " << ifc << endl;
      cout << "Reading file: " << pathtot.Data() << endl;
      
      treein_fc[ifc][imh]->SetBranchAddress("truepar0", &truepar0 );
      treein_fc[ifc][imh]->SetBranchAddress("truepar1", &truepar1 );
      treein_fc[ifc][imh]->SetBranchAddress("trueMH",   &trueMH   );
      treein_fc[ifc][imh]->SetBranchAddress("trueparname0",   &trueparname0   );
      treein_fc[ifc][imh]->SetBranchAddress("trueparname1",   &trueparname1   );
      treein_fc[ifc][imh]->SetBranchAddress("chisqtrue0",&chisqtrue0);
      treein_fc[ifc][imh]->SetBranchAddress("deltachisq0BF_glob",&deltachisq0BF_glob);
      treein_fc[ifc][imh]->SetBranchAddress("chisqmin_nh",&chisqmin_nh);
      treein_fc[ifc][imh]->SetBranchAddress("chisqmin_ih",&chisqmin_ih);
      treein_fc[ifc][imh]->SetBranchAddress("bestfit0_nh",&bestfit0_nh);
      treein_fc[ifc][imh]->SetBranchAddress("bestfit0_ih",&bestfit0_ih);
      treein_fc[ifc][imh]->SetBranchAddress("bestfit0_glob",&bestfit0_glob);
      
      
      vector<double> myvector;         // vector to compute the critical value for each file
      vector<double> vec_like_penalty; // vector of penalty weights
      
      int Nentries = treein_fc[ifc][imh]->GetEntries();
      cout << "Tot # of entries = " << Nentries << endl;
      for(int ientry=0;ientry<Nentries;ientry++){ // Loop over toy experiments
//        for(int ientry=0;ientry<500;ientry++){ // Loop over toy experiments
      
        
        if(ientry%1000==0){
          cout << "Entry " << ientry << endl;
        }
        treein_fc[ifc][imh]->GetEntry(ientry);
        
        if(ientry == 0){
          h1[ifc][imh] = new TH1D(Form("h_%d_%d",ifc,imh),Form("hist_%s_%1.3lf",hierarch[imh],truepar0),100,0,8);
          h1_bf[ifc][imh] = new TH1D(Form("h_%d_%d_bf",ifc,imh),Form("hist_%s_%1.3lf",hierarch[imh],truepar0),100,-TMath::Pi(),TMath::Pi());
          
        }
        double deltachisq0BF_used = deltachisq0BF_glob;
        double bestfit0_used = bestfit0_glob;
        
        if(fixed_hierarchy){
          if(imh>0){
            deltachisq0BF_used = chisqtrue0 - chisqmin_ih;
            bestfit0_used = bestfit0_ih;
          }
          else{
            deltachisq0BF_used = chisqtrue0 - chisqmin_nh;
            bestfit0_used = bestfit0_nh;
          }
        }
        
        myvector.push_back(deltachisq0BF_used);
        h1[ifc][imh]->Fill(deltachisq0BF_used);
        h1_bf[ifc][imh]->Fill(bestfit0_used);
        
        // if(chisqtrue0>500){
        //   cout << "chisqtrue0 = " << chisqtrue0 << endl;
        //   cout << "deltachisq0BF = " << deltachisq0BF_used << endl;
        //   exit(1);
        // }
        
        double weight = 1.;
        vec_like_penalty.push_back(weight);
        
        //cout << deltachisq0BF_glob << endl;
        
      } // loop over ientry
      
      // Get critical values
      
      //double err = sqrt(Nentries);
      double err_1sigma = sqrt(Nentries * CLfrac_1sigma * (1.-CLfrac_1sigma) ); //
      double err_90CL   = sqrt(Nentries * CLfrac_90CL   * (1.-CLfrac_90CL) );
      double err_95CL   = sqrt(Nentries * CLfrac_95CL   * (1.-CLfrac_95CL) );
      double err_2sigma = sqrt(Nentries * CLfrac_2sigma * (1.-CLfrac_2sigma) );
      double err_3sigma = sqrt(Nentries * CLfrac_3sigma * (1.-CLfrac_3sigma) );
      
      
      //CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_90CL*100.,0.);
      //exit(1);
      
      
      double FCCritVal_post_1sigma    = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_1sigma*100.,0.);
      double FCCritVal_post_1sigma_p1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_1sigma*100.,+err_1sigma / Nentries );
      double FCCritVal_post_1sigma_m1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_1sigma*100.,-err_1sigma / Nentries );
      
      double FCCritVal_post_90CL    = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_90CL*100.,0.);
      double FCCritVal_post_90CL_p1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_90CL*100.,+err_90CL / Nentries );
      double FCCritVal_post_90CL_m1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_90CL*100.,-err_90CL / Nentries );
      
      double FCCritVal_post_95CL    = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_95CL*100.,0.);
      double FCCritVal_post_95CL_p1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_95CL*100.,+err_95CL / Nentries );
      double FCCritVal_post_95CL_m1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_95CL*100.,-err_95CL / Nentries );
      
      double FCCritVal_post_2sigma    = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_2sigma*100.,0.);
      double FCCritVal_post_2sigma_p1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_2sigma*100.,+err_2sigma / Nentries );
      double FCCritVal_post_2sigma_m1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_2sigma*100.,-err_2sigma / Nentries );
      
      double FCCritVal_post_3sigma    = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_3sigma*100.,0.);
      double FCCritVal_post_3sigma_p1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_3sigma*100.,+err_3sigma / Nentries );
      double FCCritVal_post_3sigma_m1 = CalcDLcrit_posterior(myvector,vec_like_penalty,CLfrac_3sigma*100.,-err_3sigma / Nentries );
      
      double FCCritVal_1sigma    = CalcDLcrit(myvector,CLfrac_1sigma*100.,0.);
      double FCCritVal_1sigma_p1 = CalcDLcrit(myvector,CLfrac_1sigma*100.,+err_1sigma);
      double FCCritVal_1sigma_m1 = CalcDLcrit(myvector,CLfrac_1sigma*100.,-err_1sigma);
      
      double FCCritVal_90CL    = CalcDLcrit(myvector,CLfrac_90CL*100.,0.);
      double FCCritVal_90CL_p1 = CalcDLcrit(myvector,CLfrac_90CL*100.,+err_90CL);
      double FCCritVal_90CL_m1 = CalcDLcrit(myvector,CLfrac_90CL*100.,-err_90CL);
      
      double FCCritVal_95CL    = CalcDLcrit(myvector,CLfrac_95CL*100.,0.);
      double FCCritVal_95CL_p1 = CalcDLcrit(myvector,CLfrac_95CL*100.,+err_95CL);
      double FCCritVal_95CL_m1 = CalcDLcrit(myvector,CLfrac_95CL*100.,-err_95CL);
      
      double FCCritVal_2sigma    = CalcDLcrit(myvector,CLfrac_2sigma*100.,0.);
      double FCCritVal_2sigma_p1 = CalcDLcrit(myvector,CLfrac_2sigma*100.,+err_2sigma);
      double FCCritVal_2sigma_m1 = CalcDLcrit(myvector,CLfrac_2sigma*100.,-err_2sigma);
      
      double FCCritVal_3sigma    = CalcDLcrit(myvector,CLfrac_3sigma*100.,0.);
      double FCCritVal_3sigma_p1 = CalcDLcrit(myvector,CLfrac_3sigma*100.,+err_3sigma);
      double FCCritVal_3sigma_m1 = CalcDLcrit(myvector,CLfrac_3sigma*100.,-err_3sigma);
      
      
      streamer << "\\multirow{3}{*}{\\parbox{2cm}{\\centering"
      <<  std::setprecision(4) << truepar0 << "\\\\ (" << *trueMH << ")}} " << endl;
      streamer << "& 1$\\sigma$ & " << FCCritVal_1sigma << " & " << CalcCoverage(myvector,1) << " \\\\ \\cline{2-4}" << endl;
      streamer << "& 90\\% & " << FCCritVal_90CL << " & " << CalcCoverage(myvector,2.71) << " \\\\ \\cline{2-4}" << endl;
      streamer << "& 2$\\sigma$ & " << FCCritVal_2sigma << " & " << CalcCoverage(myvector,4) << " \\\\ \\cline{1-4}" << endl << endl;
      
      
      
      
      if(*trueMH=="normal"){
        gFC90CL_nh   ->SetPoint(gFC90CL_nh   ->GetN(),truepar0,FCCritVal_90CL);
        gFC90CL_nh_p1->SetPoint(gFC90CL_nh_p1->GetN(),truepar0,FCCritVal_90CL_p1);
        gFC90CL_nh_m1->SetPoint(gFC90CL_nh_m1->GetN(),truepar0,FCCritVal_90CL_m1);

        gFC1sigma_nh   ->SetPoint(gFC1sigma_nh   ->GetN(),truepar0,FCCritVal_1sigma);
        gFC1sigma_nh_p1->SetPoint(gFC1sigma_nh_p1->GetN(),truepar0,FCCritVal_1sigma_p1);
        gFC1sigma_nh_m1->SetPoint(gFC1sigma_nh_m1->GetN(),truepar0,FCCritVal_1sigma_m1);
        
        gFC2sigma_nh   ->SetPoint(gFC2sigma_nh   ->GetN(),truepar0,FCCritVal_2sigma);
        gFC2sigma_nh_p1->SetPoint(gFC2sigma_nh_p1->GetN(),truepar0,FCCritVal_2sigma_p1);
        gFC2sigma_nh_m1->SetPoint(gFC2sigma_nh_m1->GetN(),truepar0,FCCritVal_2sigma_m1);
        
        gFC3sigma_nh   ->SetPoint(gFC3sigma_nh   ->GetN(),truepar0,FCCritVal_3sigma);
        gFC3sigma_nh_p1->SetPoint(gFC3sigma_nh_p1->GetN(),truepar0,FCCritVal_3sigma_p1);
        gFC3sigma_nh_m1->SetPoint(gFC3sigma_nh_m1->GetN(),truepar0,FCCritVal_3sigma_m1);
        
        h_coverage_1sigma_nh->SetPoint(h_coverage_1sigma_nh->GetN(),truepar0,1 - CalcCoverage(myvector,1));
        h_coverage_2sigma_nh->SetPoint(h_coverage_2sigma_nh->GetN(),truepar0,1 - CalcCoverage(myvector,4));
        h_coverage_90pc_nh->SetPoint(h_coverage_90pc_nh->GetN(),truepar0, 1 - CalcCoverage(myvector,2.71));

      }
      else if(*trueMH=="inverted"){
        gFC90CL_ih   ->SetPoint(gFC90CL_ih   ->GetN(),truepar0,FCCritVal_90CL);
        gFC90CL_ih_p1->SetPoint(gFC90CL_ih_p1->GetN(),truepar0,FCCritVal_90CL_p1);
        gFC90CL_ih_m1->SetPoint(gFC90CL_ih_m1->GetN(),truepar0,FCCritVal_90CL_m1);

        gFC1sigma_ih   ->SetPoint(gFC1sigma_ih   ->GetN(),truepar0,FCCritVal_1sigma);
        gFC1sigma_ih_p1->SetPoint(gFC1sigma_ih_p1->GetN(),truepar0,FCCritVal_1sigma_p1);
        gFC1sigma_ih_m1->SetPoint(gFC1sigma_ih_m1->GetN(),truepar0,FCCritVal_1sigma_m1);
        
        gFC2sigma_ih   ->SetPoint(gFC2sigma_ih   ->GetN(),truepar0,FCCritVal_2sigma);
        gFC2sigma_ih_p1->SetPoint(gFC2sigma_ih_p1->GetN(),truepar0,FCCritVal_2sigma_p1);
        gFC2sigma_ih_m1->SetPoint(gFC2sigma_ih_m1->GetN(),truepar0,FCCritVal_2sigma_m1);
        
        gFC3sigma_ih   ->SetPoint(gFC3sigma_ih   ->GetN(),truepar0,FCCritVal_3sigma);
        gFC3sigma_ih_p1->SetPoint(gFC3sigma_ih_p1->GetN(),truepar0,FCCritVal_3sigma_p1);
        gFC3sigma_ih_m1->SetPoint(gFC3sigma_ih_m1->GetN(),truepar0,FCCritVal_3sigma_m1);
        
        h_coverage_1sigma_ih->SetPoint(h_coverage_1sigma_ih->GetN(),truepar0,1 - CalcCoverage(myvector,1));
        h_coverage_2sigma_ih->SetPoint(h_coverage_2sigma_ih->GetN(),truepar0,1 - CalcCoverage(myvector,4));
        h_coverage_90pc_ih->SetPoint(h_coverage_90pc_ih->GetN(),truepar0,1 - CalcCoverage(myvector,2.71));

      }
      else {
        cerr << endl;
        cerr << "The mass hierarchy is not defined: " << *trueMH << endl;
        cerr << endl;
      }
      
      double histogram_maximum = h1[ifc][imh]->GetMaximum();
      int bin = h_draw_68->FindBin(FCCritVal_1sigma);
      h_draw_68->Reset();
      h_draw_68->SetBinContent(bin,histogram_maximum*0.8);
      bin = h_draw_90->FindBin(FCCritVal_90CL);
      h_draw_90->Reset();
      h_draw_90->SetBinContent(bin,histogram_maximum*0.8);
      bin = h_draw_2sig->FindBin(FCCritVal_2sigma);
      h_draw_2sig->Reset();
      h_draw_2sig->SetBinContent(bin,histogram_maximum*0.8);
      h_draw_3sig->Reset();
      h_draw_3sig->SetBinContent(bin,histogram_maximum*0.8);

      TLegend *l = new TLegend(0.65,0.65,0.85,0.9);
      l->SetLineColor(0);
      l->AddEntry(h1[ifc][imh],"#Delta#chi^{2}","l");
      l->AddEntry(h_draw_68,Form("1#sigmaCL: %1.2lf",FCCritVal_1sigma),"l");
      l->AddEntry(h_draw_90,Form("90%%CL: %1.2lf",FCCritVal_90CL),"l");
      l->AddEntry(h_draw_2sig,Form("2#sigmaCL: %1.2lf",FCCritVal_2sigma),"l");
      l->AddEntry(h_draw_3sig,Form("3#sigmaCL: %1.2lf",FCCritVal_3sigma),"l");
      
      TCanvas* c = new TCanvas();
      h1[ifc][imh]->GetXaxis()->SetTitle("#Delta#chi^{2}");
      h1[ifc][imh]->GetYaxis()->SetTitle("N Toys");
      h1[ifc][imh]->SetLineWidth(2);
      h1_bf[ifc][imh]->GetXaxis()->SetTitle("#delta_{CP}");
      h1_bf[ifc][imh]->GetYaxis()->SetTitle("N Toys");
      h1_bf[ifc][imh]->SetLineWidth(2);
      h1[ifc][imh]->Draw();
      h_draw_68->Draw("same");
      h_draw_90->Draw("same");
      h_draw_2sig->Draw("same");
      h_draw_3sig->Draw("same");
      h1[ifc][imh]->Draw("axissame");
      l->Draw();
      c->SaveAs(Form("%s.pdf",h1[ifc][imh]->GetTitle()));
      
      TCanvas* c2 = new TCanvas();
      h1_bf[ifc][imh]->Draw();
      c2->SaveAs(Form("%s_bestfits.pdf",h1[ifc][imh]->GetTitle()));
      
      cout << endl;
      cout << "True Mass Hierarchy: " << *trueMH << endl;
      cout << *trueparname0 << " = " << truepar0 << endl;
      cout << *trueparname1 << " = " << truepar1 << endl;
      cout << "Critical values:"<< endl;
      cout << " - 1sigma = " << FCCritVal_1sigma
      << " + " << fabs(FCCritVal_1sigma-FCCritVal_1sigma_p1)
      << " - " << fabs(FCCritVal_1sigma-FCCritVal_1sigma_m1)
      << endl;
      cout << " - 90% CL = " << FCCritVal_90CL
      << " + " << fabs(FCCritVal_90CL-FCCritVal_90CL_p1)
      << " - " << fabs(FCCritVal_90CL-FCCritVal_90CL_m1)
      << endl;
      cout << " - 95% CL = " << FCCritVal_95CL
      << " + " << fabs(FCCritVal_95CL-FCCritVal_95CL_p1)
      << " - " << fabs(FCCritVal_95CL-FCCritVal_95CL_m1)
      << endl;
      cout << " - 2sigma = " << FCCritVal_2sigma
      << " + " << fabs(FCCritVal_2sigma-FCCritVal_2sigma_p1)
      << " - " << fabs(FCCritVal_2sigma-FCCritVal_2sigma_m1)
      << endl;
      cout << " - 3sigma = " << FCCritVal_3sigma
      << " + " << fabs(FCCritVal_3sigma-FCCritVal_3sigma_p1)
      << " - " << fabs(FCCritVal_3sigma-FCCritVal_3sigma_m1)
      << endl;
      
      cout << "Critical values posterior-like weight=1 (Used):"<< endl;
      cout << " - 1sigma = " << FCCritVal_post_1sigma
      << " + " << fabs(FCCritVal_post_1sigma-FCCritVal_post_1sigma_p1)
      << " - " << fabs(FCCritVal_post_1sigma-FCCritVal_post_1sigma_m1)
      << endl;
      cout << " - 90% CL = " << FCCritVal_post_90CL
      << " + " << fabs(FCCritVal_post_90CL-FCCritVal_post_90CL_p1)
      << " - " << fabs(FCCritVal_post_90CL-FCCritVal_post_90CL_m1)
      << endl;
      cout << " - 95% CL = " << FCCritVal_post_95CL
      << " + " << fabs(FCCritVal_post_95CL-FCCritVal_post_95CL_p1)
      << " - " << fabs(FCCritVal_post_95CL-FCCritVal_post_95CL_m1)
      << endl;
      cout << " - 2sigma = " << FCCritVal_post_2sigma
      << " + " << fabs(FCCritVal_post_2sigma-FCCritVal_post_2sigma_p1)
      << " - " << fabs(FCCritVal_post_2sigma-FCCritVal_post_2sigma_m1)
      << endl;
      cout << " - 3sigma = " << FCCritVal_post_3sigma
      << " + " << fabs(FCCritVal_post_3sigma-FCCritVal_post_3sigma_p1)
      << " - " << fabs(FCCritVal_post_3sigma-FCCritVal_post_3sigma_m1)
      << endl;
      
    } // loop over ifc
  } // loop over imh
  
  
  streamer << "\\end{tabular}"   << endl;
  streamer << "\\end{table}" << endl;


  // get size of file
  long size;
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);
  // write to outfile
  outfile.write (streamer.str().c_str(),size);
  outfile.close();
  
  
  // Copy +pi in -pi
  
  cout << endl;
  
  for(int ipt=0;ipt<gFC90CL_nh->GetN();ipt++){
    
    double x_90CL_nh,      y_90CL_nh=0;
    double x_90CL_nh_p1,   y_90CL_nh_p1=0;
    double x_90CL_nh_m1,   y_90CL_nh_m1=0;
    double x_90CL_ih,      y_90CL_ih=0;
    double x_90CL_ih_p1,   y_90CL_ih_p1=0;
    double x_90CL_ih_m1,   y_90CL_ih_m1=0;

    double x_1sigma_nh,    y_1sigma_nh=0;
    double x_1sigma_nh_p1, y_1sigma_nh_p1=0;
    double x_1sigma_nh_m1, y_1sigma_nh_m1=0;
    double x_1sigma_ih,    y_1sigma_ih=0;
    double x_1sigma_ih_p1, y_1sigma_ih_p1=0;
    double x_1sigma_ih_m1, y_1sigma_ih_m1=0;
    
    double x_2sigma_nh,    y_2sigma_nh=0;
    double x_2sigma_nh_p1, y_2sigma_nh_p1=0;
    double x_2sigma_nh_m1, y_2sigma_nh_m1=0;
    double x_2sigma_ih,    y_2sigma_ih=0;
    double x_2sigma_ih_p1, y_2sigma_ih_p1=0;
    double x_2sigma_ih_m1, y_2sigma_ih_m1=0;

    double x_3sigma_nh,    y_3sigma_nh=0;
    double x_3sigma_nh_p1, y_3sigma_nh_p1=0;
    double x_3sigma_nh_m1, y_3sigma_nh_m1=0;
    double x_3sigma_ih,    y_3sigma_ih=0;
    double x_3sigma_ih_p1, y_3sigma_ih_p1=0;
    double x_3sigma_ih_m1, y_3sigma_ih_m1=0;
    
    gFC90CL_nh->GetPoint(ipt, x_90CL_nh, y_90CL_nh);
    gFC90CL_nh_p1->GetPoint(ipt, x_90CL_nh_p1, y_90CL_nh_p1);
    gFC90CL_nh_m1->GetPoint(ipt, x_90CL_nh_m1, y_90CL_nh_m1);
    gFC90CL_ih->GetPoint(ipt, x_90CL_ih, y_90CL_ih);
    gFC90CL_ih_p1->GetPoint(ipt, x_90CL_ih_p1, y_90CL_ih_p1);
    gFC90CL_ih_m1->GetPoint(ipt, x_90CL_ih_m1, y_90CL_ih_m1);

    gFC1sigma_nh->GetPoint(ipt, x_1sigma_nh, y_1sigma_nh);
    gFC1sigma_nh_p1->GetPoint(ipt, x_1sigma_nh_p1, y_1sigma_nh_p1);
    gFC1sigma_nh_m1->GetPoint(ipt, x_1sigma_nh_m1, y_1sigma_nh_m1);
    gFC1sigma_ih->GetPoint(ipt, x_1sigma_ih, y_1sigma_ih);
    gFC1sigma_ih_p1->GetPoint(ipt, x_1sigma_ih_p1, y_1sigma_ih_p1);
    gFC1sigma_ih_m1->GetPoint(ipt, x_1sigma_ih_m1, y_1sigma_ih_m1);

    gFC2sigma_nh->GetPoint(ipt, x_2sigma_nh, y_2sigma_nh);
    gFC2sigma_nh_p1->GetPoint(ipt, x_2sigma_nh_p1, y_2sigma_nh_p1);
    gFC2sigma_nh_m1->GetPoint(ipt, x_2sigma_nh_m1, y_2sigma_nh_m1);
    gFC2sigma_ih->GetPoint(ipt, x_2sigma_ih, y_2sigma_ih);
    gFC2sigma_ih_p1->GetPoint(ipt, x_2sigma_ih_p1, y_2sigma_ih_p1);
    gFC2sigma_ih_m1->GetPoint(ipt, x_2sigma_ih_m1, y_2sigma_ih_m1);

    gFC3sigma_nh->GetPoint(ipt, x_3sigma_nh, y_3sigma_nh);
    gFC3sigma_nh_p1->GetPoint(ipt, x_3sigma_nh_p1, y_3sigma_nh_p1);
    gFC3sigma_nh_m1->GetPoint(ipt, x_3sigma_nh_m1, y_3sigma_nh_m1);
    gFC3sigma_ih->GetPoint(ipt, x_3sigma_ih, y_3sigma_ih);
    gFC3sigma_ih_p1->GetPoint(ipt, x_3sigma_ih_p1, y_3sigma_ih_p1);
    gFC3sigma_ih_m1->GetPoint(ipt, x_3sigma_ih_m1, y_3sigma_ih_m1);
    
    for(int iptcpy=0;iptcpy<nptcpy;iptcpy++){
      if( fabs( x_90CL_nh - copypt[iptcpy][0] ) < 0.0001 ){
        
        cout << "Copy critical value of point " << copypt[iptcpy][0]
        << " in " << copypt[iptcpy][1] << endl;
        cout << endl;
        
        gFC90CL_nh   ->SetPoint(gFC90CL_nh   ->GetN(),copypt[iptcpy][1],y_90CL_nh);
        gFC90CL_nh_p1->SetPoint(gFC90CL_nh_p1->GetN(),copypt[iptcpy][1],y_90CL_nh_p1);
        gFC90CL_nh_m1->SetPoint(gFC90CL_nh_m1->GetN(),copypt[iptcpy][1],y_90CL_nh_m1);
        gFC90CL_ih   ->SetPoint(gFC90CL_ih   ->GetN(),copypt[iptcpy][1],y_90CL_ih);
        gFC90CL_ih_p1->SetPoint(gFC90CL_ih_p1->GetN(),copypt[iptcpy][1],y_90CL_ih_p1);
        gFC90CL_ih_m1->SetPoint(gFC90CL_ih_m1->GetN(),copypt[iptcpy][1],y_90CL_ih_m1);

        gFC1sigma_nh   ->SetPoint(gFC1sigma_nh   ->GetN(),copypt[iptcpy][1],y_1sigma_nh);
        gFC1sigma_nh_p1->SetPoint(gFC1sigma_nh_p1->GetN(),copypt[iptcpy][1],y_1sigma_nh_p1);
        gFC1sigma_nh_m1->SetPoint(gFC1sigma_nh_m1->GetN(),copypt[iptcpy][1],y_1sigma_nh_m1);
        gFC1sigma_ih   ->SetPoint(gFC1sigma_ih   ->GetN(),copypt[iptcpy][1],y_1sigma_ih);
        gFC1sigma_ih_p1->SetPoint(gFC1sigma_ih_p1->GetN(),copypt[iptcpy][1],y_1sigma_ih_p1);
        gFC1sigma_ih_m1->SetPoint(gFC1sigma_ih_m1->GetN(),copypt[iptcpy][1],y_1sigma_ih_m1);

        gFC2sigma_nh   ->SetPoint(gFC2sigma_nh   ->GetN(),copypt[iptcpy][1],y_2sigma_nh);
        gFC2sigma_nh_p1->SetPoint(gFC2sigma_nh_p1->GetN(),copypt[iptcpy][1],y_2sigma_nh_p1);
        gFC2sigma_nh_m1->SetPoint(gFC2sigma_nh_m1->GetN(),copypt[iptcpy][1],y_2sigma_nh_m1);
        gFC2sigma_ih   ->SetPoint(gFC2sigma_ih   ->GetN(),copypt[iptcpy][1],y_2sigma_ih);
        gFC2sigma_ih_p1->SetPoint(gFC2sigma_ih_p1->GetN(),copypt[iptcpy][1],y_2sigma_ih_p1);
        gFC2sigma_ih_m1->SetPoint(gFC2sigma_ih_m1->GetN(),copypt[iptcpy][1],y_2sigma_ih_m1);

        gFC3sigma_nh   ->SetPoint(gFC3sigma_nh   ->GetN(),copypt[iptcpy][1],y_3sigma_nh);
        gFC3sigma_nh_p1->SetPoint(gFC3sigma_nh_p1->GetN(),copypt[iptcpy][1],y_3sigma_nh_p1);
        gFC3sigma_nh_m1->SetPoint(gFC3sigma_nh_m1->GetN(),copypt[iptcpy][1],y_3sigma_nh_m1);
        gFC3sigma_ih   ->SetPoint(gFC3sigma_ih   ->GetN(),copypt[iptcpy][1],y_3sigma_ih);
        gFC3sigma_ih_p1->SetPoint(gFC3sigma_ih_p1->GetN(),copypt[iptcpy][1],y_3sigma_ih_p1);
        gFC3sigma_ih_m1->SetPoint(gFC3sigma_ih_m1->GetN(),copypt[iptcpy][1],y_3sigma_ih_m1);
      }
    }
  }
  
  //
  
  h_coverage_1sigma_nh->Sort();
  h_coverage_2sigma_nh->Sort();
  h_coverage_90pc_nh->Sort();
  
  h_coverage_1sigma_ih->Sort();
  h_coverage_2sigma_ih->Sort();
  h_coverage_90pc_ih->Sort();
  
  gFC90CL_nh   ->Sort();
  gFC90CL_nh_p1->Sort();
  gFC90CL_nh_m1->Sort();
  gFC90CL_ih   ->Sort();
  gFC90CL_ih_p1->Sort();
  gFC90CL_ih_m1->Sort();

  gFC1sigma_nh   ->Sort();
  gFC1sigma_nh_p1->Sort();
  gFC1sigma_nh_m1->Sort();
  gFC1sigma_ih   ->Sort();
  gFC1sigma_ih_p1->Sort();
  gFC1sigma_ih_m1->Sort();
  
  gFC2sigma_nh   ->Sort();
  gFC2sigma_nh_p1->Sort();
  gFC2sigma_nh_m1->Sort();
  gFC2sigma_ih   ->Sort();
  gFC2sigma_ih_p1->Sort();
  gFC2sigma_ih_m1->Sort();

  gFC3sigma_nh   ->Sort();
  gFC3sigma_nh_p1->Sort();
  gFC3sigma_nh_m1->Sort();
  gFC3sigma_ih   ->Sort();
  gFC3sigma_ih_p1->Sort();
  gFC3sigma_ih_m1->Sort();
  
  gFC90CL_nh   ->SetLineColor(kRed);
  gFC90CL_nh_p1->SetLineColor(kRed);
  gFC90CL_nh_m1->SetLineColor(kRed);
  gFC90CL_ih   ->SetLineColor(kRed);
  gFC90CL_ih_p1->SetLineColor(kRed);
  gFC90CL_ih_m1->SetLineColor(kRed);

  gFC1sigma_nh   ->SetLineColor(kOrange-3);
  gFC1sigma_nh_p1->SetLineColor(kOrange-3);
  gFC1sigma_nh_m1->SetLineColor(kOrange-3);
  gFC1sigma_ih   ->SetLineColor(kOrange-3);
  gFC1sigma_ih_p1->SetLineColor(kOrange-3);
  gFC1sigma_ih_m1->SetLineColor(kOrange-3);

  gFC2sigma_nh   ->SetLineColor(kBlue);
  gFC2sigma_nh_p1->SetLineColor(kBlue);
  gFC2sigma_nh_m1->SetLineColor(kBlue);
  gFC2sigma_ih   ->SetLineColor(kBlue);
  gFC2sigma_ih_p1->SetLineColor(kBlue);
  gFC2sigma_ih_m1->SetLineColor(kBlue);

  gFC3sigma_nh   ->SetLineColor(kGreen+3);
  gFC3sigma_nh_p1->SetLineColor(kGreen+3);
  gFC3sigma_nh_m1->SetLineColor(kGreen+3);
  gFC3sigma_ih   ->SetLineColor(kGreen+3);
  gFC3sigma_ih_p1->SetLineColor(kGreen+3);
  gFC3sigma_ih_m1->SetLineColor(kGreen+3);
  
  h_coverage_1sigma_nh->SetLineColor(kRed);
  h_coverage_2sigma_nh->SetLineColor(kBlue);
  h_coverage_90pc_nh->SetLineColor(kGreen+3);
  
  h_coverage_1sigma_ih->SetLineColor(kRed);
  h_coverage_2sigma_ih->SetLineColor(kBlue);
  h_coverage_90pc_ih->SetLineColor(kGreen+3);
  
  h_coverage_1sigma_ih->SetLineStyle(2);
  h_coverage_2sigma_ih->SetLineStyle(2);
  h_coverage_90pc_ih->SetLineStyle(2);
  
  gFC90CL_ih   ->SetLineStyle(2);
  gFC90CL_ih_p1->SetLineStyle(2);
  gFC90CL_ih_m1->SetLineStyle(2);
  gFC1sigma_ih   ->SetLineStyle(2);
  gFC1sigma_ih_p1->SetLineStyle(2);
  gFC1sigma_ih_m1->SetLineStyle(2);
  gFC2sigma_ih   ->SetLineStyle(2);
  gFC2sigma_ih_p1->SetLineStyle(2);
  gFC2sigma_ih_m1->SetLineStyle(2);
  gFC3sigma_ih   ->SetLineStyle(2);
  gFC3sigma_ih_p1->SetLineStyle(2);
  gFC3sigma_ih_m1->SetLineStyle(2);
  
  TCanvas *cDeltaChisq = new TCanvas("cDeltaChisq","cDeltaChisq");
  gFC90CL_nh->GetXaxis()->SetTitle("#delta_{CP}");
  gFC90CL_nh->GetYaxis()->SetTitle("#Delta #chi^{2}_{crit}");
//  gFC90CL_nh->GetYaxis()->SetRangeUser(0.8,5);
  gFC90CL_nh->GetYaxis()->SetRangeUser(1.5,10);

  gFC90CL_nh   ->Draw("apl");
  gFC90CL_nh_p1->Draw("pl same");
  gFC90CL_nh_m1->Draw("pl same");
  gFC2sigma_nh   ->Draw("pl same");
  gFC2sigma_nh_p1->Draw("pl same");
  gFC2sigma_nh_m1->Draw("pl same");
//  gFC3sigma_nh   ->Draw("pl same");
//  gFC3sigma_nh_p1->Draw("pl same");
//  gFC3sigma_nh_m1->Draw("pl same");
  gFC90CL_ih   ->Draw("pl same");
  gFC90CL_ih_p1->Draw("pl same");
  gFC90CL_ih_m1->Draw("pl same");
  gFC2sigma_ih   ->Draw("pl same");
  gFC2sigma_ih_p1->Draw("pl same");
  gFC2sigma_ih_m1->Draw("pl same");
//  gFC3sigma_ih   ->Draw("pl same");
//  gFC3sigma_ih_p1->Draw("pl same");
//  gFC3sigma_ih_m1->Draw("pl same");
  
  TLegend *lDeltaChisq = new TLegend(0.75,0.2,0.9,0.35);
  lDeltaChisq->SetLineColor(0);
  lDeltaChisq->SetFillStyle(0);
  lDeltaChisq->AddEntry(gFC90CL_nh,"90%","l");
  lDeltaChisq->AddEntry(gFC2sigma_nh,"2#sigma","l");
  lDeltaChisq->AddEntry(gFC3sigma_nh,"3#sigma","l");
  lDeltaChisq->AddEntry(gFC90CL_nh,"NH","l");

  lDeltaChisq->AddEntry(gFC90CL_ih,"90%","l");
  lDeltaChisq->AddEntry(gFC2sigma_ih,"2#sigma","l");
  lDeltaChisq->AddEntry(gFC3sigma_ih,"3#sigma","l");
  lDeltaChisq->AddEntry(gFC90CL_ih,"IH","l");

  lDeltaChisq->Draw();
//  cDeltaChisq->SaveAs("FC_critical_values.pdf");


  TCanvas *cDeltaChisq_p1sigma = new TCanvas("cDeltaChisq_p1sigma","cDeltaChisq_p1sigma");

  gFC90CL_nh->GetXaxis()->SetTitle("#delta_{CP}");
  gFC90CL_nh->GetYaxis()->SetTitle("#Delta #chi^{2}_{crit}");
//  gFC90CL_nh->GetYaxis()->SetRangeUser(0.8,5);
  gFC90CL_nh->GetYaxis()->SetRangeUser(0,5);

  gFC90CL_nh   ->Draw("apl");
  gFC90CL_nh_p1->Draw("pl same");
  gFC90CL_nh_m1->Draw("pl same");
  gFC1sigma_nh   ->Draw("pl same");
  gFC1sigma_nh_p1->Draw("pl same");
  gFC1sigma_nh_m1->Draw("pl same");
  gFC2sigma_nh   ->Draw("pl same");
  gFC2sigma_nh_p1->Draw("pl same");
  gFC2sigma_nh_m1->Draw("pl same");
//  gFC3sigma_nh   ->Draw("pl same");
//  gFC3sigma_nh_p1->Draw("pl same");
//  gFC3sigma_nh_m1->Draw("pl same");

  gFC90CL_ih   ->Draw("pl same");
  gFC90CL_ih_p1->Draw("pl same");
  gFC90CL_ih_m1->Draw("pl same");
  gFC1sigma_ih   ->Draw("pl same");
  gFC1sigma_ih_p1->Draw("pl same");
  gFC1sigma_ih_m1->Draw("pl same");
  gFC2sigma_ih   ->Draw("pl same");
  gFC2sigma_ih_p1->Draw("pl same");
  gFC2sigma_ih_m1->Draw("pl same");
//  gFC3sigma_ih   ->Draw("pl same");
//  gFC3sigma_ih_p1->Draw("pl same");
//  gFC3sigma_ih_m1->Draw("pl same");
  
  TLegend *lDeltaChisq_p1sigma = new TLegend(0.75,0.2,1.0,0.35);
  lDeltaChisq_p1sigma->SetLineColor(0);
  lDeltaChisq_p1sigma->SetFillStyle(0);
  lDeltaChisq_p1sigma->AddEntry(gFC1sigma_nh,"1#sigma","l");  
  lDeltaChisq_p1sigma->AddEntry(gFC90CL_nh,"90%","l");
  lDeltaChisq_p1sigma->AddEntry(gFC2sigma_nh,"2#sigma","l");
//  lDeltaChisq_p1sigma->AddEntry(gFC3sigma_nh,"3#sigma","l");
  lDeltaChisq_p1sigma->AddEntry(gFC2sigma_nh,"NH","l");

  lDeltaChisq_p1sigma->AddEntry(gFC2sigma_ih,"IH","l");

  lDeltaChisq_p1sigma->Draw();

  TLatex* title = new TLatex(0.0, 1.01 * 5.0, "T2K Run1-8 Preliminary");
  title->SetTextSize(0.05);
  title->Draw();
  TLatex* caveat = new TLatex(-2.5, 1.01 * 5.0, "Final systematics pending");
  caveat->SetTextSize(0.025);
  caveat->Draw();

  cDeltaChisq_p1sigma->SaveAs("FC_critical_values.pdf");
  cDeltaChisq_p1sigma->SaveAs("FC_critical_values.png");
  cDeltaChisq_p1sigma->SaveAs("FC_critical_values.eps");

  TFile *fileout_fc = new TFile("FC_critvals.root","RECREATE");

  cDeltaChisq->Write();
  cDeltaChisq_p1sigma->Write();  
  fileout_fc->Close();
  
  TCanvas *CoverageCanvas = new TCanvas("CoverageCanvas","CoverageCanvas");
  h_coverage_1sigma_nh->GetXaxis()->SetTitle("#delta_{CP}");
  h_coverage_1sigma_nh->GetYaxis()->SetTitle("1 - Coverage");
  CoverageCanvas->SetLogy();
  TF1 *fa1 = new TF1("fa1","1-0.68",-4,4);
  TF1 *fa2 = new TF1("fa2","1-0.90",-4,4);
  TF1 *fa3 = new TF1("fa3","1-0.95",-4,4);

  h_coverage_1sigma_nh->GetYaxis()->SetRangeUser(0.001,1);
//  h_coverage_1sigma_nh->GetXaxis()->SetRangeUser(-3.1412,3.1412);
  h_coverage_1sigma_nh->Draw("apl");
  fa1->SetLineColorAlpha(17, 0.35);
  fa2->SetLineColorAlpha(17, 0.35);
  fa3->SetLineColorAlpha(17, 0.35);
  fa1->SetLineWidth(1);
  fa2->SetLineWidth(1);
  fa3->SetLineWidth(1);

  fa1->Draw("lsame");
  fa2->Draw("lsame");
  fa3->Draw("lsame");
  h_coverage_1sigma_nh->Draw("pl");
  h_coverage_2sigma_nh->Draw("pl");
  h_coverage_90pc_nh->Draw("pl");
  h_coverage_1sigma_ih->Draw("pl");
  h_coverage_2sigma_ih->Draw("pl");
  h_coverage_90pc_ih->Draw("pl");
  
  TLegend *leg = new TLegend(0.6,0.2,0.88,0.4);
  leg->SetLineColor(0);
  leg->AddEntry(h_coverage_1sigma_nh,"1#sigma (68%)","l");
  leg->AddEntry(h_coverage_90pc_nh,"90%","l");
  leg->AddEntry(h_coverage_2sigma_nh,"2#sigma (95%)","l");
  leg->AddEntry(h_coverage_1sigma_ih,"IH","l");
  leg->Draw();
  CoverageCanvas->SaveAs("Coverage.pdf");
  CoverageCanvas->SaveAs("Coverage.root");


  cout << "Writing TGraph File" << endl;
  
  TFile *fileout_fc = new TFile("FCgraphs.root","RECREATE");
  gFC90CL_nh    ->Write();
  gFC90CL_nh_p1->Write();
  gFC90CL_nh_m1->Write();
  gFC90CL_ih   ->Write();
  gFC90CL_ih_p1->Write();
  gFC90CL_ih_m1->Write();
  
  gFC1sigma_nh   ->Write();
  gFC1sigma_nh_p1->Write();
  gFC1sigma_nh_m1->Write();
  gFC1sigma_ih   ->Write();
  gFC1sigma_ih_p1->Write();
  gFC1sigma_ih_m1->Write();

  gFC2sigma_nh   ->Write();
  gFC2sigma_nh_p1->Write();
  gFC2sigma_nh_m1->Write();
  gFC2sigma_ih   ->Write();
  gFC2sigma_ih_p1->Write();
  gFC2sigma_ih_m1->Write();

  gFC3sigma_nh   ->Write();
  gFC3sigma_nh_p1->Write();
  gFC3sigma_nh_m1->Write();
  gFC3sigma_ih   ->Write();
  gFC3sigma_ih_p1->Write();
  gFC3sigma_ih_m1->Write();
  
  h_coverage_1sigma_nh->Write("h_coverage_1sigma_nh");
  h_coverage_2sigma_nh->Write("h_coverage_2sigma_nh");
  h_coverage_90pc_nh->Write("h_coverage_90pc_nh");
  h_coverage_1sigma_ih->Write("h_coverage_1sigma_ih");
  h_coverage_2sigma_ih->Write("h_coverage_2sigma_ih");
  h_coverage_90pc_ih->Write("h_coverage_90pc_ih");
  fileout_fc->Close();

  cout << "Finished processing" << endl;
  return;
}


double CalcDLcrit(vector<double> dchisq, double conflev, double sigmacl){
  
  // Get error on critical value: sigmacl = sqrt(N*p*(1-p))
  
  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector
  vector<double>::iterator dchisqIter;
  dchisqIter = dchisq.begin() + (int)(0.01*conflev*dchisq.size() + sigmacl);  // (int)(...) gets one toy more after the X%CL edge (if even number)
  double critval = (*dchisqIter);
  return critval;
}

double CalcCoverage(vector<double> dchisq, double crit_dchisq){
  
  // Get error on critical value: sigmacl = sqrt(N*p*(1-p))
  
  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector
  for (unsigned i=0; i < dchisq.size(); i++) {
    if( crit_dchisq < dchisq[i])
      return double(i)/dchisq.size();
  }

  return -999;
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
      //if( sum >= (0.01 * conflev + sigmacl)*tot_weight ){
      
      //cout << endl;
      //cout << "sum = " << sum << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      
      critval = new_dchisq[index[i]];
      break;
    }
  }
  
  return critval;
}
