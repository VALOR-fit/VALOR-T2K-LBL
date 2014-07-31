#include<iostream>
#include<fstream>
#include<stdio.h>
#include<iomanip>
#include<math.h>
#include<TH1.h>
#include<TH1F.h>
#include<TH1D.h>
#include<TH2.h>
#include<TH2F.h>
#include<TH2D.h>
#include<TF1.h>
#include<TSystem.h>
#include<TFile.h>
#include<TCanvas.h>
#include<TMinuit.h>
#include<TMath.h>
#include<TString.h>
#include<TGraph.h>
#include<TLegend.h>
#include<TLegendEntry.h>
#include<TRandom3.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TChain.h>
#include<TRandom3.h>
#include<TFile.h>
#include<TFitResult.h>
#include<TLine.h>
#include<TStyle.h>
#include<TVector.h>
#include<TVectorD.h>

using namespace std;

//const char infilename[200] = "selected/neut/13a/sk13a_nRnueCC.root";
const char infilename[200] = "selected/neut/13a/sk13a_nRnumuCC.root";

TH1F *herec = new TH1F("herec","herec",100,0,10);
TH1F *hmode = new TH1F("hmode","hmode",200,-100,100);
TH1F *henu = new TH1F("henu","henu",100,0,10);
TH1F *hwgt = new TH1F("hwgt","hwgt",100,0,3);
TH1F *hplep = new TH1F("hplep","hplep",100,0,10);
TH1F *hcostheta = new TH1F("hcostheta","hcostheta",100,-1,1);
TH1F *hQ2 = new TH1F("hQ2","hQ2",100,0,10);
TH1F *hevis = new TH1F("hevis","hevis",100,0,10);
TH1I *hnring = new TH1I("hnring","hnring",10,0,10);

void DrawTreeOutput(){

  TFile *infile = new TFile(infilename,"Read");

  // Read the Tree
  
  TChain *cnue = new TChain("tnue"); 
  cnue->AddFile(infilename);
  TTree *tnue = (TTree*)cnue;

  TChain *cnuebar = new TChain("tnuebar"); 
  cnuebar->AddFile(infilename);
  TTree *tnuebar = (TTree*)cnuebar;

  TChain *cnuesig = new TChain("tnuesig"); 
  cnuesig->AddFile(infilename);
  TTree *tnuesig = (TTree*)cnuesig;

  TChain *cnumu = new TChain("tnumu"); 
  cnumu->AddFile(infilename);
  TTree *tnumu = (TTree*)cnumu;

  TChain *cnumubar = new TChain("tnumubar"); 
  cnumubar->AddFile(infilename);
  TTree *tnumubar = (TTree*)cnumubar;

  //

  float erec_nue;
  float erec_nuebar;
  float erec_nuesig;
  float erec_numu;
  float erec_numubar;
  int mode_nue;
  int mode_nuebar;
  int mode_nuesig;
  int mode_numu;
  int mode_numubar;
  float enu_nue;
  float enu_nuebar;
  float enu_nuesig;
  float enu_numu;
  float enu_numubar;
  float wgt_nue;
  float wgt_nuebar;
  float wgt_nuesig;
  float wgt_numu;
  float wgt_numubar;
  float plep_nue;
  float plep_nuebar;
  float plep_nuesig;
  float plep_numu;
  float plep_numubar;
  float costheta_nue;
  float costheta_nuebar;
  float costheta_nuesig;
  float costheta_numu;
  float costheta_numubar;
  float Q2_nue;
  float Q2_nuebar;
  float Q2_nuesig;
  float Q2_numu;
  float Q2_numubar;
  float evis_nue;
  float evis_nuebar;
  float evis_nuesig;
  float evis_numu;
  float evis_numubar;
  int nring_nue;
  int nring_nuebar;
  int nring_nuesig;
  int nring_numu;
  int nring_numubar;

  TBranch *b_erec_nue;   
  TBranch *b_erec_nuebar;   
  TBranch *b_erec_nuesig;   
  TBranch *b_erec_numu;   
  TBranch *b_erec_numubar;   
  TBranch *b_mode_nue;   
  TBranch *b_mode_nuebar;   
  TBranch *b_mode_nuesig;   
  TBranch *b_mode_numu;   
  TBranch *b_mode_numubar;   
  TBranch *b_enu_nue;   
  TBranch *b_enu_nuebar;   
  TBranch *b_enu_nuesig;   
  TBranch *b_enu_numu;   
  TBranch *b_enu_numubar;   
  TBranch *b_wgt_nue;   
  TBranch *b_wgt_nuebar;   
  TBranch *b_wgt_nuesig;   
  TBranch *b_wgt_numu;   
  TBranch *b_wgt_numubar;   
  TBranch *b_plep_nue;   
  TBranch *b_plep_nuebar;   
  TBranch *b_plep_nuesig;   
  TBranch *b_plep_numu;   
  TBranch *b_plep_numubar;   
  TBranch *b_costheta_nue;   
  TBranch *b_costheta_nuebar;   
  TBranch *b_costheta_nuesig;   
  TBranch *b_costheta_numu;   
  TBranch *b_costheta_numubar;   
  TBranch *b_Q2_nue;   
  TBranch *b_Q2_nuebar;   
  TBranch *b_Q2_nuesig;   
  TBranch *b_Q2_numu;   
  TBranch *b_Q2_numubar;   
  TBranch *b_evis_nue;   
  TBranch *b_evis_nuebar;   
  TBranch *b_evis_nuesig;   
  TBranch *b_evis_numu;   
  TBranch *b_evis_numubar;   
  TBranch *b_nring_nue;   
  TBranch *b_nring_nuebar;   
  TBranch *b_nring_nuesig;   
  TBranch *b_nring_numu;   
  TBranch *b_nring_numubar;   

  // tnue 
  cout << endl;
  cout << "Read tree: tnue" << endl;
  cout << endl;
  tnue->SetBranchStatus("*", 0);	
  tnue->SetBranchAddress("erec", &erec_nue, &b_erec_nue);
  tnue->SetBranchAddress("mode", &mode_nue, &b_mode_nue);
  tnue->SetBranchAddress("enu", &enu_nue, &b_enu_nue);
  tnue->SetBranchAddress("wgt", &wgt_nue, &b_wgt_nue);
  tnue->SetBranchAddress("plep", &plep_nue, &b_plep_nue);
  tnue->SetBranchAddress("costheta", &costheta_nue, &b_costheta_nue);
  tnue->SetBranchAddress("Q2", &Q2_nue, &b_Q2_nue);
  tnue->SetBranchAddress("evis", &evis_nue, &b_evis_nue);
  tnue->SetBranchAddress("nring", &nring_nue, &b_nring_nue);
  Int_t nevent = (Int_t)tnue->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnue->GetEntry(i);	 // Get the Event			
    herec->Fill(erec_nue/1000.);
    hmode->Fill(mode_nue);
    henu->Fill(enu_nue);
    hwgt->Fill(wgt_nue);
    hplep->Fill(plep_nue);
    hcostheta->Fill(costheta_nue);
    hQ2->Fill(Q2_nue);
    hevis->Fill(evis_nue);
    hnring->Fill(nring_nue);
  }

  // tnuebar 
  cout << endl;
  cout << "Read tree: tnuebar" << endl;
  cout << endl;
  tnuebar->SetBranchStatus("*", 0);	
  tnuebar->SetBranchAddress("erec", &erec_nuebar, &b_erec_nuebar);
  tnuebar->SetBranchAddress("mode", &mode_nuebar, &b_mode_nuebar);
  tnuebar->SetBranchAddress("enu", &enu_nuebar, &b_enu_nuebar);
  tnuebar->SetBranchAddress("wgt", &wgt_nuebar, &b_wgt_nuebar);
  tnuebar->SetBranchAddress("plep", &plep_nuebar, &b_plep_nuebar);
  tnuebar->SetBranchAddress("costheta", &costheta_nuebar, &b_costheta_nuebar);
  tnuebar->SetBranchAddress("Q2", &Q2_nuebar, &b_Q2_nuebar);
  tnuebar->SetBranchAddress("evis", &evis_nuebar, &b_evis_nuebar);
  tnuebar->SetBranchAddress("nring", &nring_nuebar, &b_nring_nuebar);
  Int_t nevent = (Int_t)tnuebar->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnuebar->GetEntry(i);	 // Get the Event			
    herec->Fill(erec_nuebar/1000.);
    hmode->Fill(mode_nuebar);
    henu->Fill(enu_nuebar);
    hwgt->Fill(wgt_nuebar);
    hplep->Fill(plep_nuebar);
    hcostheta->Fill(costheta_nuebar);
    hQ2->Fill(Q2_nuebar);
    hevis->Fill(evis_nuebar);
    hnring->Fill(nring_nuebar);
  }

  // tnuesig 
  cout << endl;
  cout << "Read tree: tnuesig" << endl;
  cout << endl;
  tnuesig->SetBranchStatus("*", 0);	
  tnuesig->SetBranchAddress("erec", &erec_nuesig, &b_erec_nuesig);
  tnuesig->SetBranchAddress("mode", &mode_nuesig, &b_mode_nuesig);
  tnuesig->SetBranchAddress("enu", &enu_nuesig, &b_enu_nuesig);
  tnuesig->SetBranchAddress("wgt", &wgt_nuesig, &b_wgt_nuesig);
  tnuesig->SetBranchAddress("plep", &plep_nuesig, &b_plep_nuesig);
  tnuesig->SetBranchAddress("costheta", &costheta_nuesig, &b_costheta_nuesig);
  tnuesig->SetBranchAddress("Q2", &Q2_nuesig, &b_Q2_nuesig);
  tnuesig->SetBranchAddress("evis", &evis_nuesig, &b_evis_nuesig);
  tnuesig->SetBranchAddress("nring", &nring_nuesig, &b_nring_nuesig);
  Int_t nevent = (Int_t)tnuesig->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnuesig->GetEntry(i);	 // Get the Event			
    herec->Fill(erec_nuesig/1000.);
    hmode->Fill(mode_nuesig);
    henu->Fill(enu_nuesig);
    hwgt->Fill(wgt_nuesig);
    hplep->Fill(plep_nuesig);
    hcostheta->Fill(costheta_nuesig);
    hQ2->Fill(Q2_nuesig);
    hevis->Fill(evis_nuesig);
    hnring->Fill(nring_nuesig);
  }

  // tnumu 
  cout << endl;
  cout << "Read tree: tnumu" << endl;
  cout << endl;
  tnumu->SetBranchStatus("*", 0);	
  tnumu->SetBranchAddress("erec", &erec_numu, &b_erec_numu);
  tnumu->SetBranchAddress("mode", &mode_numu, &b_mode_numu);
  tnumu->SetBranchAddress("enu", &enu_numu, &b_enu_numu);
  tnumu->SetBranchAddress("wgt", &wgt_numu, &b_wgt_numu);
  tnumu->SetBranchAddress("plep", &plep_numu, &b_plep_numu);
  tnumu->SetBranchAddress("costheta", &costheta_numu, &b_costheta_numu);
  tnumu->SetBranchAddress("Q2", &Q2_numu, &b_Q2_numu);
  tnumu->SetBranchAddress("evis", &evis_numu, &b_evis_numu);
  tnumu->SetBranchAddress("nring", &nring_numu, &b_nring_numu);
  Int_t nevent = (Int_t)tnumu->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnumu->GetEntry(i);	 // Get the Event			
    herec->Fill(erec_numu/1000.);
    hmode->Fill(mode_numu);
    henu->Fill(enu_numu);
    hwgt->Fill(wgt_numu);
    hplep->Fill(plep_numu);
    hcostheta->Fill(costheta_numu);
    hQ2->Fill(Q2_numu);
    hevis->Fill(evis_numu);
    hnring->Fill(nring_numu);
  }

  // tnumubar 
  cout << endl;
  cout << "Read tree: tnumubar" << endl;
  cout << endl;
  tnumubar->SetBranchStatus("*", 0);	
  tnumubar->SetBranchAddress("erec", &erec_numubar, &b_erec_numubar);
  tnumubar->SetBranchAddress("mode", &mode_numubar, &b_mode_numubar);
  tnumubar->SetBranchAddress("enu", &enu_numubar, &b_enu_numubar);
  tnumubar->SetBranchAddress("wgt", &wgt_numubar, &b_wgt_numubar);
  tnumubar->SetBranchAddress("plep", &plep_numubar, &b_plep_numubar);
  tnumubar->SetBranchAddress("costheta", &costheta_numubar, &b_costheta_numubar);
  tnumubar->SetBranchAddress("Q2", &Q2_numubar, &b_Q2_numubar);
  tnumubar->SetBranchAddress("evis", &evis_numubar, &b_evis_numubar);
  tnumubar->SetBranchAddress("nring", &nring_numubar, &b_nring_numubar);
  Int_t nevent = (Int_t)tnumubar->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnumubar->GetEntry(i);	 // Get the Event			
    herec->Fill(erec_numubar/1000.);
    hmode->Fill(mode_numubar);
    henu->Fill(enu_numubar);
    hwgt->Fill(wgt_numubar);
    hplep->Fill(plep_numubar);
    hcostheta->Fill(costheta_numubar);
    hQ2->Fill(Q2_numubar);
    hevis->Fill(evis_numubar);
    hnring->Fill(nring_numubar);
  }


  TCanvas *cerec = new TCanvas("cerec","cerec");
  herec->Draw();

  TCanvas *cmode = new TCanvas("cmode","cmode");
  hmode->Draw();

  TCanvas *cenu = new TCanvas("cenu","cenu");
  henu->Draw();

  TCanvas *cwgt = new TCanvas("cwgt","cwgt");
  hwgt->Draw();

  TCanvas *cplep = new TCanvas("cplep","cplep");
  hplep->Draw();

  TCanvas *ccostheta = new TCanvas("ccostheta","ccostheta");
  hcostheta->Draw();

  TCanvas *cQ2 = new TCanvas("cQ2","cQ2");
  hQ2->Draw();

  TCanvas *cevis = new TCanvas("cevis","cevis");
  hevis->Draw();

  TCanvas *cnring = new TCanvas("cnring","cnring");
  hnring->Draw();

  return;
}
