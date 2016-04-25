
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

  const int Nmh = 1;
  const int NFCpts = 8; // # of points of Feldman-Cousins

  string path_fc_nh = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/";
  string path_fc_ih = "../TECH_NOTE_14_3_16/FELDMAN_COUSINS/1D_DCP_VS_MH/";
  
  string infilename_fc[NFCpts][Nmh] = {
    "fit_fitmarg_fc_dcp_normal_nexp_20_-2.356.root",
    "fit_fitmarg_fc_dcp_normal_nexp_20_-1.571.root",
    "fit_fitmarg_fc_dcp_normal_nexp_20_-0.785.root",
    "fit_fitmarg_fc_dcp_normal_nexp_20_0.000.root",
    "fit_fitmarg_fc_dcp_normal_nexp_20_0.785.root",
    "fit_fitmarg_fc_dcp_normal_nexp_20_1.571.root",
    "fit_fitmarg_fc_dcp_normal_nexp_20_2.356.root",
    "fit_fitmarg_fc_dcp_normal_nexp_20_3.141.root" // +-pi
  };

  TGraph *gFC90CL_nh = new TGraph();
  gFC90CL_nh->SetName("gFC90CL_nh");
  gFC90CL_nh->SetTitle("gFC90CL_nh");
  TGraph *gFC90CL_ih = new TGraph();
  gFC90CL_ih->SetName("gFC90CL_ih");
  gFC90CL_ih->SetTitle("gFC90CL_ih");

  TGraph *gFC2sigma_nh = new TGraph();
  gFC2sigma_nh->SetName("gFC2sigma_nh");
  gFC2sigma_nh->SetTitle("gFC2sigma_nh");
  TGraph *gFC2sigma_ih = new TGraph();
  gFC2sigma_ih->SetName("gFC2sigma_ih");
  gFC2sigma_ih->SetTitle("gFC2sigma_ih");
  
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
	pathtot = TString::Format("%s%s",path_fc_nh.c_str(),infilename_fc[ifc][imh].c_str());
      }
      else if(imh==1){
	pathtot = TString::Format("%s%s",path_fc_ih.c_str(),infilename_fc[ifc][imh].c_str());
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

      cout << "Reading file: " << pathtot.Data() << endl;

      treein_fc[ifc][imh]->SetBranchAddress("truepar0", &truepar0 );
      treein_fc[ifc][imh]->SetBranchAddress("truepar1", &truepar1 );
      treein_fc[ifc][imh]->SetBranchAddress("trueMH",   &trueMH   );
      treein_fc[ifc][imh]->SetBranchAddress("trueparname0",   &trueparname0   );
      treein_fc[ifc][imh]->SetBranchAddress("trueparname1",   &trueparname1   );
      treein_fc[ifc][imh]->SetBranchAddress("chisqtrue0",&chisqtrue0);
      treein_fc[ifc][imh]->SetBranchAddress("deltachisq0BF_glob",&deltachisq0BF_glob);
    

      vector<double> myvector; // vector to compute the critical value for each file

      //int Nentries = treein_fc[ifc][imh]->GetEntries();
      int Nentries = 100;
 
      for(int ientry=0;ientry<Nentries;ientry++){ // Loop over toy experiments
	
	treein_fc[ifc][imh]->GetEntry(ientry);
	double deltachisq0BF_used = deltachisq0BF_glob;
	myvector.push_back(deltachisq0BF_used);
	
      } // loop over ientry
            
      // Get critical values 
      double FCCritVal_1sigma = CalcDLcrit(myvector,CLfrac_1sigma*100.,0.); 
      double FCCritVal_90CL   = CalcDLcrit(myvector,CLfrac_90CL*100.,0.); 
      double FCCritVal_95CL   = CalcDLcrit(myvector,CLfrac_95CL*100.,0.); 
      double FCCritVal_2sigma = CalcDLcrit(myvector,CLfrac_2sigma*100.,0.); 

      if(*trueMH=="normal"){
      	gFC90CL_nh  ->SetPoint(gFC90CL_nh->GetN(),  truepar0,FCCritVal_90CL);
	gFC2sigma_nh->SetPoint(gFC2sigma_nh->GetN(),truepar0,FCCritVal_2sigma);
      }
      else if(*trueMH=="inverted"){
	gFC90CL_ih  ->SetPoint(gFC90CL_ih->GetN(),  truepar0,FCCritVal_90CL);
	gFC2sigma_ih->SetPoint(gFC2sigma_ih->GetN(),truepar0,FCCritVal_2sigma);
      }
      else {
	cerr << endl;
	cerr << "The mass hierarchy is not defined: " << *trueMH << endl;
	cerr << endl;
      }

    } // loop over ifc
  } // loop over imh

  gFC90CL_nh->Sort();
  gFC90CL_ih->Sort();
  gFC2sigma_nh->Sort();
  gFC2sigma_ih->Sort();

  gFC90CL_nh->SetLineColor(kRed);
  gFC90CL_ih->SetLineColor(kRed);
  gFC2sigma_nh->SetLineColor(kBlue);
  gFC2sigma_ih->SetLineColor(kBlue);

  gFC90CL_ih->SetLineStyle(2);
  gFC2sigma_ih->SetLineStyle(2);
  
  TCanvas *cDeltaChisq = new TCanvas("cDeltaChisq","cDeltaChisq");
  gFC90CL_nh->Draw("apl");
  gFC2sigma_nh->Draw("pl same");

  TFile *fileout_fc = new TFile("FCgraphs.root","RECREATE");
  gFC90CL_nh->Write();
  gFC90CL_ih->Write();
  gFC2sigma_nh->Write();
  gFC2sigma_ih->Write();
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
