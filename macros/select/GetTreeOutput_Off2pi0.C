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

const char infilename[200] = "$T2KOSCANA/inputs/official/2013/superk/mc/neut/13a/smallnt_nue2rpi0_ap13a_fqv3r0_flx11bv32run1234_fqpi0cut1250_v20130423.root";
const char outfilename[200] = "$T2KOSCANA/inputs/test/sk13a_2rpi0_fitTQun_Nishimura.root";

TH1F *herec = new TH1F("herec","herec",100,0,10);
TH1F *hmode = new TH1F("hmode","hmode",200,-100,100);
TH1F *henu = new TH1F("henu","henu",100,0,10);
TH1F *hwgt = new TH1F("hwgt","hwgt",100,0,3);

void GetTreeOutput_Off2pi0(){

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
  
  // output trees

  TFile *outfile = new TFile(outfilename,"RECREATE"); // file containing the tree

  TTree *tnue_out = new TTree("tnue","recreate");
  TTree *tnuebar_out = new TTree("tnuebar","recreate");
  TTree *tnuesig_out = new TTree("tnuesig","recreate");
  TTree *tnumu_out = new TTree("tnumu","recreate");
  TTree *tnumubar_out = new TTree("tnumubar","recreate");

  //

  int evflag;
  float erec;
  int mode;
  float enu;
  float wgt;

  TBranch *b_evflag;
  TBranch *b_erec;   
  TBranch *b_mode;   
  TBranch *b_enu;   
  TBranch *b_wgt;   

  int evflag_out;
  float erec_out;
  int mode_out;
  float enu_out;
  double wgt_out;
  
  // tnue 

  tnue_out->Branch("evflag",&evflag_out,"evflag_out/I");
  tnue_out->Branch("erec",&erec_out,"erec/F");
  tnue_out->Branch("mode",&mode_out,"mode/I");
  tnue_out->Branch("enu",&enu_out,"enu/F");
  tnue_out->Branch("wgt",&wgt_out,"wgt/D");

  cout << endl;
  cout << "Read tree: tnue" << endl;
  cout << endl;
  tnue->SetBranchStatus("*", 0);	
  tnue->SetBranchAddress("evflag", &evflag, &b_evflag);
  tnue->SetBranchAddress("erec", &erec, &b_erec);
  tnue->SetBranchAddress("mode", &mode, &b_mode);
  tnue->SetBranchAddress("enu", &enu, &b_enu);
  tnue->SetBranchAddress("wgt", &wgt, &b_wgt);
  Int_t nevent = (Int_t)tnue->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnue->GetEntry(i);	 // Get the Event			

    if( !(evflag & 0x10) ) continue; // select 2pi0 events

    // fill output tree
    evflag_out = (int) evflag;
    erec_out =  (float) erec;
    mode_out =  (int) mode;
    enu_out =  (float) enu;
    wgt_out =  (double) wgt;
    tnue_out->Fill();

    herec->Fill(erec_out/1000.);
    hmode->Fill(mode_out);
    henu->Fill(enu_out);
    hwgt->Fill(wgt_out);
  }

  // tnuebar 

  tnuebar_out->Branch("evflag",&evflag_out,"evflag_out/I");
  tnuebar_out->Branch("erec",&erec_out,"erec/F");
  tnuebar_out->Branch("mode",&mode_out,"mode/I");
  tnuebar_out->Branch("enu",&enu_out,"enu/F");
  tnuebar_out->Branch("wgt",&wgt_out,"wgt/D");

  cout << endl;
  cout << "Read tree: tnuebar" << endl;
  cout << endl;
  tnuebar->SetBranchStatus("*", 0);
  tnuebar->SetBranchAddress("evflag", &evflag, &b_evflag);	
  tnuebar->SetBranchAddress("erec", &erec, &b_erec);
  tnuebar->SetBranchAddress("mode", &mode, &b_mode);
  tnuebar->SetBranchAddress("enu", &enu, &b_enu);
  tnuebar->SetBranchAddress("wgt", &wgt, &b_wgt);
  Int_t nevent = (Int_t)tnuebar->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnuebar->GetEntry(i);	 // Get the Event			

    if( !(evflag & 0x10) ) continue; // select 2pi0 events

    // fill output tree
    evflag_out = (int) evflag;
    erec_out =  (float) erec;
    mode_out =  (int) mode;
    enu_out =  (float) enu;
    wgt_out =  (double) wgt;
    tnuebar_out->Fill();

    herec->Fill(erec_out/1000.);
    hmode->Fill(mode_out);
    henu->Fill(enu_out);
    hwgt->Fill(wgt_out);
  }

  // tnuesig 

  tnuesig_out->Branch("evflag",&evflag_out,"evflag_out/I");
  tnuesig_out->Branch("erec",&erec_out,"erec/F");
  tnuesig_out->Branch("mode",&mode_out,"mode/I");
  tnuesig_out->Branch("enu",&enu_out,"enu/F");
  tnuesig_out->Branch("wgt",&wgt_out,"wgt/D");

  cout << endl;
  cout << "Read tree: tnuesig" << endl;
  cout << endl;
  tnuesig->SetBranchStatus("*", 0);
  tnuesig->SetBranchAddress("evflag", &evflag, &b_evflag);	
  tnuesig->SetBranchAddress("erec", &erec, &b_erec);
  tnuesig->SetBranchAddress("mode", &mode, &b_mode);
  tnuesig->SetBranchAddress("enu", &enu, &b_enu);
  tnuesig->SetBranchAddress("wgt", &wgt, &b_wgt);
  Int_t nevent = (Int_t)tnuesig->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnuesig->GetEntry(i);	 // Get the Event

    if( !(evflag & 0x10) ) continue; // select 2pi0 events

    // fill output tree
    evflag_out = (int) evflag;
    erec_out =  (float) erec;
    mode_out =  (int) mode;
    enu_out =  (float) enu;
    wgt_out =  (double) wgt;
    tnuesig_out->Fill();

    herec->Fill(erec_out/1000.);
    hmode->Fill(mode_out);
    henu->Fill(enu_out);
    hwgt->Fill(wgt_out);
  }

  // tnumu 

  tnumu_out->Branch("evflag",&evflag_out,"evflag_out/I");
  tnumu_out->Branch("erec",&erec_out,"erec/F");
  tnumu_out->Branch("mode",&mode_out,"mode/I");
  tnumu_out->Branch("enu",&enu_out,"enu/F");
  tnumu_out->Branch("wgt",&wgt_out,"wgt/D");

  cout << endl;
  cout << "Read tree: tnumu" << endl;
  cout << endl;
  tnumu->SetBranchStatus("*", 0);
  tnumu->SetBranchAddress("evflag", &evflag, &b_evflag);	
  tnumu->SetBranchAddress("erec", &erec, &b_erec);
  tnumu->SetBranchAddress("mode", &mode, &b_mode);
  tnumu->SetBranchAddress("enu", &enu, &b_enu);
  tnumu->SetBranchAddress("wgt", &wgt, &b_wgt);
  Int_t nevent = (Int_t)tnumu->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnumu->GetEntry(i);	 // Get the Event			

    if( !(evflag & 0x10) ) continue; // select 2pi0 events

    // fill output tree
    evflag_out = (int) evflag;
    erec_out =  (float) erec;
    mode_out =  (int) mode;
    enu_out =  (float) enu;
    wgt_out =  (double) wgt;
    tnumu_out->Fill();

    herec->Fill(erec_out/1000.);
    hmode->Fill(mode_out);
    henu->Fill(enu_out);
    hwgt->Fill(wgt_out);
  }

  // tnumubar 

  tnumubar_out->Branch("evflag",&evflag_out,"evflag_out/I");
  tnumubar_out->Branch("erec",&erec_out,"erec/F");
  tnumubar_out->Branch("mode",&mode_out,"mode/I");
  tnumubar_out->Branch("enu",&enu_out,"enu/F");
  tnumubar_out->Branch("wgt",&wgt_out,"wgt/D");

  cout << endl;
  cout << "Read tree: tnumubar" << endl;
  cout << endl;
  tnumubar->SetBranchStatus("*", 0);
  tnumubar->SetBranchAddress("evflag", &evflag, &b_evflag);
  tnumubar->SetBranchAddress("erec", &erec, &b_erec);
  tnumubar->SetBranchAddress("mode", &mode, &b_mode);
  tnumubar->SetBranchAddress("enu", &enu, &b_enu);
  tnumubar->SetBranchAddress("wgt", &wgt, &b_wgt);
  Int_t nevent = (Int_t)tnumubar->GetEntries();					
  for(int i=0;i<nevent;i++){           
    if( (i%1000)==0 ) cout << "Event " << i << endl;
    tnumubar->GetEntry(i);	 // Get the Event			

    if( !(evflag & 0x10) ) continue; // select 2pi0 events

    // fill output tree
    evflag_out = (int) evflag;
    erec_out =  (float) erec;
    mode_out =  (int) mode;
    enu_out =  (float) enu;
    wgt_out =  (double) wgt;
    tnumubar_out->Fill();

    herec->Fill(erec_out/1000.);
    hmode->Fill(mode_out);
    henu->Fill(enu_out);
    hwgt->Fill(wgt_out);
  }

  // Write trees in file

  tnue_out->Write();
  tnuebar_out->Write();
  tnuesig_out->Write();
  tnumu_out->Write();
  tnumubar_out->Write();
  outfile->Close();
  
  // Draw tree content
  
  TCanvas *cerec = new TCanvas("cerec","cerec");
  herec->Draw();
  TCanvas *cmode = new TCanvas("cmode","cmode");
  hmode->Draw();
  TCanvas *cenu = new TCanvas("cenu","cenu");
  henu->Draw();
  TCanvas *cwgt = new TCanvas("cwgt","cwgt");
  hwgt->Draw();

  return;
}
