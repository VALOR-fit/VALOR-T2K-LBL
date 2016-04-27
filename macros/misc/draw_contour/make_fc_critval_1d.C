
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

void make_fc_critval_1d(){

  const int Nmh = 2;
  const int NFCpts = 8; // # of points of Feldman-Cousins

  string path_fc_nh = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/NH/";
  string path_fc_ih = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/IH/";

  const char * infilename_fc[NFCpts][Nmh] = {
    
    // // Inverted Hierarchy - 5k Toys
    // "total_tree_-0.785.root",
    // "total_tree_-2.356.root",
    // "total_tree_0.785.root",
    // "total_tree_2.356.root",
    // "total_tree_halfpi.root",
    // "total_tree_mhalfpi.root",
    // "total_tree_pi.root",
    // "total_tree_zero.root"

    // // Normal Hierarchy - 5k Toys
    // "fit_fitmarg_fc_dcp_normal_nexp_20_-2.356.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_20_-1.571.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_20_-0.785.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_20_0.000.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_20_0.785.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_20_1.571.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_20_2.356.root",
    // //"fit_fitmarg_fc_dcp_normal_nexp_20_3.141.root" // +-pi
    // "fit_fitmarg_fc_dcp_normal_nexp_20_3.142.root", // +-pi 

    // // Normal Hierarchy - 10k Toys 
    // "fit_fitmarg_fc_dcp_normal_nexp_10k_-0.785.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_10k_-1.571.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_10k_-2.356.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_10k_0.000.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_10k_0.785.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_10k_1.571.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_10k_2.356.root",
    // "fit_fitmarg_fc_dcp_normal_nexp_10k_3.142.root" // +-pi

    // Normal Hierarhcy 10k Toys - Inverted Hierarchy 5k Toys 
    {"fit_fitmarg_fc_dcp_normal_nexp_10k_-0.785.root", "total_tree_-0.785.root"},
    {"fit_fitmarg_fc_dcp_normal_nexp_10k_-2.356.root", "total_tree_-2.356.root"},
    {"fit_fitmarg_fc_dcp_normal_nexp_10k_-1.571.root", "total_tree_mhalfpi.root"},
    {"fit_fitmarg_fc_dcp_normal_nexp_10k_-0.785.root", "total_tree_-0.785.root"},
    {"fit_fitmarg_fc_dcp_normal_nexp_10k_0.000.root" , "total_tree_zero.root"},
    {"fit_fitmarg_fc_dcp_normal_nexp_10k_1.571.root" , "total_tree_halfpi.root"},
    {"fit_fitmarg_fc_dcp_normal_nexp_10k_2.356.root" , "total_tree_-2.356.root"},
    {"fit_fitmarg_fc_dcp_normal_nexp_10k_3.142.root" , "total_tree_pi.root"} // +-pi

  };

  const int nptcpy = 1;
  double copypt[nptcpy][2] = {
    //{ TMath::Pi(), -TMath::Pi() } // copy crit val of TMath::Pi() in -TMath::Pi()
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
      cout << "Reading file: " << pathtot.Data() << endl;

      treein_fc[ifc][imh]->SetBranchAddress("truepar0", &truepar0 );
      treein_fc[ifc][imh]->SetBranchAddress("truepar1", &truepar1 );
      treein_fc[ifc][imh]->SetBranchAddress("trueMH",   &trueMH   );
      treein_fc[ifc][imh]->SetBranchAddress("trueparname0",   &trueparname0   );
      treein_fc[ifc][imh]->SetBranchAddress("trueparname1",   &trueparname1   );
      treein_fc[ifc][imh]->SetBranchAddress("chisqtrue0",&chisqtrue0);
      treein_fc[ifc][imh]->SetBranchAddress("deltachisq0BF_glob",&deltachisq0BF_glob);
    

      vector<double> myvector; // vector to compute the critical value for each file

      int Nentries = treein_fc[ifc][imh]->GetEntries();

      for(int ientry=0;ientry<Nentries;ientry++){ // Loop over toy experiments
	
	if(ientry%1000==0){
	  cout << "Entry " << ientry << endl;
	}
	treein_fc[ifc][imh]->GetEntry(ientry);
	double deltachisq0BF_used = deltachisq0BF_glob;
	myvector.push_back(deltachisq0BF_used);
	
      } // loop over ientry
            
      // Get critical values 
      
      //double err = sqrt(Nentries);
      double err_1sigma = sqrt(Nentries * CLfrac_1sigma * (1.-CLfrac_1sigma) );
      double err_90CL   = sqrt(Nentries * CLfrac_90CL   * (1.-CLfrac_90CL) );
      double err_95CL   = sqrt(Nentries * CLfrac_95CL   * (1.-CLfrac_95CL) );
      double err_2sigma = sqrt(Nentries * CLfrac_2sigma * (1.-CLfrac_2sigma) );
      
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

      if(*trueMH=="normal"){
      	gFC90CL_nh   ->SetPoint(gFC90CL_nh   ->GetN(),truepar0,FCCritVal_90CL);
      	gFC90CL_nh_p1->SetPoint(gFC90CL_nh_p1->GetN(),truepar0,FCCritVal_90CL_p1);
      	gFC90CL_nh_m1->SetPoint(gFC90CL_nh_m1->GetN(),truepar0,FCCritVal_90CL_m1);

      	gFC2sigma_nh   ->SetPoint(gFC2sigma_nh   ->GetN(),truepar0,FCCritVal_2sigma);
      	gFC2sigma_nh_p1->SetPoint(gFC2sigma_nh_p1->GetN(),truepar0,FCCritVal_2sigma_p1);
      	gFC2sigma_nh_m1->SetPoint(gFC2sigma_nh_m1->GetN(),truepar0,FCCritVal_2sigma_m1);
      }
      else if(*trueMH=="inverted"){
      	gFC90CL_ih   ->SetPoint(gFC90CL_ih   ->GetN(),truepar0,FCCritVal_90CL);
      	gFC90CL_ih_p1->SetPoint(gFC90CL_ih_p1->GetN(),truepar0,FCCritVal_90CL_p1);
      	gFC90CL_ih_m1->SetPoint(gFC90CL_ih_m1->GetN(),truepar0,FCCritVal_90CL_m1);

      	gFC2sigma_ih   ->SetPoint(gFC2sigma_ih   ->GetN(),truepar0,FCCritVal_2sigma);
      	gFC2sigma_ih_p1->SetPoint(gFC2sigma_ih_p1->GetN(),truepar0,FCCritVal_2sigma_p1);
      	gFC2sigma_ih_m1->SetPoint(gFC2sigma_ih_m1->GetN(),truepar0,FCCritVal_2sigma_m1);
      }
      else {
	cerr << endl;
	cerr << "The mass hierarchy is not defined: " << *trueMH << endl;
	cerr << endl;
      }

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

    } // loop over ifc
  } // loop over imh

  
  // Copy +pi in -pi
 
  cout << endl;
  
  for(int ipt=0;ipt<gFC90CL_nh->GetN();ipt++){
    
    double x_90CL_nh,      y_90CL_nh=0;
    double x_90CL_nh_p1,   y_90CL_nh_p1=0;
    double x_90CL_nh_m1,   y_90CL_nh_m1=0;
    double x_90CL_ih,      y_90CL_ih=0;
    double x_90CL_ih_p1,   y_90CL_ih_p1=0;
    double x_90CL_ih_m1,   y_90CL_ih_m1=0;
    double x_2sigma_nh,    y_2sigma_nh=0;
    double x_2sigma_nh_p1, y_2sigma_nh_p1=0;
    double x_2sigma_nh_m1, y_2sigma_nh_m1=0;
    double x_2sigma_ih,    y_2sigma_ih=0;
    double x_2sigma_ih_p1, y_2sigma_ih_p1=0;
    double x_2sigma_ih_m1, y_2sigma_ih_m1=0;
    
    gFC90CL_nh->GetPoint(ipt, x_90CL_nh, y_90CL_nh);
    gFC90CL_nh_p1->GetPoint(ipt, x_90CL_nh_p1, y_90CL_nh_p1);
    gFC90CL_nh_m1->GetPoint(ipt, x_90CL_nh_m1, y_90CL_nh_m1);
    gFC90CL_ih->GetPoint(ipt, x_90CL_ih, y_90CL_ih);
    gFC90CL_ih_p1->GetPoint(ipt, x_90CL_ih_p1, y_90CL_ih_p1);
    gFC90CL_ih_m1->GetPoint(ipt, x_90CL_ih_m1, y_90CL_ih_m1);
    gFC2sigma_nh->GetPoint(ipt, x_2sigma_nh, y_2sigma_nh);
    gFC2sigma_nh_p1->GetPoint(ipt, x_2sigma_nh_p1, y_2sigma_nh_p1);
    gFC2sigma_nh_m1->GetPoint(ipt, x_2sigma_nh_m1, y_2sigma_nh_m1);
    gFC2sigma_ih->GetPoint(ipt, x_2sigma_ih, y_2sigma_ih);
    gFC2sigma_ih_p1->GetPoint(ipt, x_2sigma_ih_p1, y_2sigma_ih_p1);
    gFC2sigma_ih_m1->GetPoint(ipt, x_2sigma_ih_m1, y_2sigma_ih_m1);
    
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
 	gFC2sigma_nh   ->SetPoint(gFC2sigma_nh   ->GetN(),copypt[iptcpy][1],y_2sigma_nh);
 	gFC2sigma_nh_p1->SetPoint(gFC2sigma_nh_p1->GetN(),copypt[iptcpy][1],y_2sigma_nh_p1);
 	gFC2sigma_nh_m1->SetPoint(gFC2sigma_nh_m1->GetN(),copypt[iptcpy][1],y_2sigma_nh_m1);
	gFC2sigma_ih   ->SetPoint(gFC2sigma_ih   ->GetN(),copypt[iptcpy][1],y_2sigma_ih);
	gFC2sigma_ih_p1->SetPoint(gFC2sigma_ih_p1->GetN(),copypt[iptcpy][1],y_2sigma_ih_p1);
	gFC2sigma_ih_m1->SetPoint(gFC2sigma_ih_m1->GetN(),copypt[iptcpy][1],y_2sigma_ih_m1);
      }
    }
  }

  //

  gFC90CL_nh   ->Sort();
  gFC90CL_nh_p1->Sort();
  gFC90CL_nh_m1->Sort();
  gFC90CL_ih   ->Sort();
  gFC90CL_ih_p1->Sort();
  gFC90CL_ih_m1->Sort();
  gFC2sigma_nh   ->Sort();
  gFC2sigma_nh_p1->Sort();
  gFC2sigma_nh_m1->Sort();
  gFC2sigma_ih   ->Sort();
  gFC2sigma_ih_p1->Sort();
  gFC2sigma_ih_m1->Sort();

  gFC90CL_nh   ->SetLineColor(kRed);
  gFC90CL_nh_p1->SetLineColor(kRed);
  gFC90CL_nh_m1->SetLineColor(kRed);
  gFC90CL_ih   ->SetLineColor(kRed);
  gFC90CL_ih_p1->SetLineColor(kRed);
  gFC90CL_ih_m1->SetLineColor(kRed);
  gFC2sigma_nh   ->SetLineColor(kBlue);
  gFC2sigma_nh_p1->SetLineColor(kBlue);
  gFC2sigma_nh_m1->SetLineColor(kBlue);
  gFC2sigma_ih   ->SetLineColor(kBlue);
  gFC2sigma_ih_p1->SetLineColor(kBlue);
  gFC2sigma_ih_m1->SetLineColor(kBlue);

  gFC90CL_ih   ->SetLineStyle(2);
  gFC90CL_ih_p1->SetLineStyle(2);
  gFC90CL_ih_m1->SetLineStyle(2);
  gFC2sigma_ih   ->SetLineStyle(2);
  gFC2sigma_ih_p1->SetLineStyle(2);
  gFC2sigma_ih_m1->SetLineStyle(2);
  
  TCanvas *cDeltaChisq = new TCanvas("cDeltaChisq","cDeltaChisq");
  gFC90CL_nh   ->Draw("apl");
  gFC90CL_nh_p1->Draw("pl same");
  gFC90CL_nh_m1->Draw("pl same");
  gFC2sigma_nh   ->Draw("pl same");
  gFC2sigma_nh_p1->Draw("pl same");
  gFC2sigma_nh_m1->Draw("pl same");
  gFC90CL_ih   ->Draw("pl same");
  gFC90CL_ih_p1->Draw("pl same");
  gFC90CL_ih_m1->Draw("pl same");
  gFC2sigma_ih   ->Draw("pl same");
  gFC2sigma_ih_p1->Draw("pl same");
  gFC2sigma_ih_m1->Draw("pl same");

  TFile *fileout_fc = new TFile("FCgraphs.root","RECREATE");
  gFC90CL_nh    ->Write();
  gFC90CL_nh_p1->Write();
  gFC90CL_nh_m1->Write();  
  gFC90CL_ih   ->Write();
  gFC90CL_ih_p1->Write();
  gFC90CL_ih_m1->Write();  
  gFC2sigma_nh   ->Write();
  gFC2sigma_nh_p1->Write();
  gFC2sigma_nh_m1->Write();  
  gFC2sigma_ih   ->Write();
  gFC2sigma_ih_p1->Write();
  gFC2sigma_ih_m1->Write();  
  fileout_fc->Close();
  
  return;
}


double CalcDLcrit(vector<double> dchisq, double conflev, double sigmacl){
  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector   
  vector<double>::iterator dchisqIter;
  dchisqIter = dchisq.begin() + (int)(0.01*conflev*dchisq.size() + sigmacl);
  double critval = (*dchisqIter);
  return critval;
}
