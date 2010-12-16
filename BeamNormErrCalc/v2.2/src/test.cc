#define INPUT_DIR "/home/kubo1/beam_svn/analysis_code/BeamNormErrCalc/trunk/share"

#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TTimeStamp.h>
#include <iostream>

#include "HadErrCalc.h"

#define PRINT 0

using namespace std;

int main(int argc, char * argv[]){
  
  TTimeStamp t1;

  HadErrCalc *HEC = new HadErrCalc(INPUT_DIR);
  
  /* read sample historgams */
  TFile * f1 = TFile::Open("../share/fc1rmu.root");
  TFile * f2 = TFile::Open("../share/fc1re.root");
  TFile * f3 = TFile::Open("../share/tmp_beam.root");
  TString head = "h_enu_int_";
  TString flav[4] = {"Numu_","Numb_","Nue_","SigNue_"};
  TString mode[7] = {"CCQE", "CC1PI", "CCCOH", "CCOTHER", "NC1PI0", "NCCOH", "NCOTHER"};
  
  TH1D *hinp_1rmu[3];
  
  hinp_1rmu[0] = new TH1D("h1rnum","h1rnum",200,0,10);
  hinp_1rmu[1] = new TH1D("h1rnmb","h1rnmb",200,0,10);
  hinp_1rmu[2] = new TH1D("h1rnue","h1rnue",200,0,10);
  
  TH1D *hndnum = (TH1D*)f3->Get("h_nd_beam_0");
  TH1D *hndnmb = (TH1D*)f3->Get("h_nd_beam_1");
  TH1D *hndnue = (TH1D*)f3->Get("h_nd_beam_2");
  TH1D *hbgnum = (TH1D*)f3->Get("h_bkg_beam_0");
  TH1D *hbgnmb = (TH1D*)f3->Get("h_bkg_beam_1");
  TH1D *hbgnue = (TH1D*)f3->Get("h_bkg_beam_2");
  TH1D *hsig = (TH1D*)f3->Get("h_sig_beam");
  
  TString hname;
  for(int i=0;i<3;i++){
    for(int j=0;j<7;j++){
      hname = head + flav[i] + mode[j];
      TH1D *h1 = (TH1D*)f1->Get(hname);
      hinp_1rmu[i]->Add(h1);
    }
  }
  
  TTimeStamp t2;
  cout << t2.GetSec() - t1.GetSec() << "sec for init." << endl;

  /* test N times with same histograms */
  double n = 10;

  for(int i=0;i<n;i++){
    /* set input histograms */
    HEC->SetNDHist(hndnum, hndnmb, hndnue);       /* required */
    HEC->Set1RmuHist(hinp_1rmu[0], hinp_1rmu[1], hinp_1rmu[2]);     /* optional (for numu analysis) */
    HEC->SetNueSigHist(hsig);                  /* optional (for nue analysis)  */
    HEC->SetNueBGHist(hbgnum, hbgnmb, hbgnue); /* optional (for nue analysis)  */
    /* Exec. calculation */
    HEC->Calc();
    /* now, HadErrCalc has calculated results of the loop*/
    if(i%100 == 99)
      cout << "loop : " << i+1 << "   1RmuError = " << HEC->GetSumErr(0)*100 << "%%"<< endl;
  }
  
  TTimeStamp t3;
  cout << n <<" calls / " << ( t3.GetSec() - t2.GetSec() ) <<" sec." << endl;

  // print out the last results //
  cout <<  endl;
  cout << "NSK/NND errors -----------" << endl;
  cout << "1Rmu/ND error \t\t" << HEC->GetSumErr(0)*100 << " %%" << endl;
  cout << "NueSig/ND error\t\t" << HEC->GetSumErr(1)*100 << " %%"<< endl;
  cout << "NueBG/ND error \t\t" << HEC->GetSumErr(2)*100 << " %%" << endl;
  cout << "NueSig+BG/ND error\t" << HEC->GetSumErr(3)*100 << " %%" << endl;
  cout <<  endl;

  cout << "NSK,NND single errors --------" << endl;
  cout << "NND error \t" << HEC->GetNDErr()*100 << " %%" << endl;
  cout << "1Rmu error \t" << HEC->GetSKErr(0)*100 << " %%" << endl;
  cout << "NueSig error\t" << HEC->GetSKErr(1)*100 << " %%"<< endl;
  cout << "NueBG error \t" << HEC->GetSKErr(2)*100 << " %%" << endl;
  cout << "NueSig+BG error\t" << HEC->GetSKErr(3)*100 << " %%" << endl;

  cout <<  endl;
  // flavor by flavor //
  cout << "flavor by flavor --------" << endl;
  cout << "1Rmu/ND error numu\t" << HEC->GetErr(0,0)*100 << " %%"<< endl;
  cout << "1Rmu/ND error numub\t" << HEC->GetErr(0,1)*100 << " %%"<< endl;
  cout << "1Rmu/ND error nue\t" << HEC->GetErr(0,2)*100 << " %%"<< endl;
  cout << "NueBG/ND error numu\t" << HEC->GetErr(2,0)*100 << " %%"<< endl;
  cout << "NueBG/ND error numub\t" << HEC->GetErr(2,1)*100 << " %%"<< endl;
  cout << "NueBG/ND error nue\t" << HEC->GetErr(2,2)*100 << " %%"<< endl;
  cout <<  endl;
  
  // errors from each components //
  HEC->ErrorComponents();
  // verbose outputs //
  HEC->ShowDetail();

  if(PRINT){
    
    // print contribution histograms //
    TApplication app("app",0,0);
    double s[6] = {0.0, 0.2, 0.55, 0.6, 0.8, 1.0  };
    double r[6] = {0.3, 0.0, 1.0,  1.0, 1.0, 1.0  };
    double g[6] = {0.0, 0.0, 1.0,  1.0, 0.5, 0.0  };
    double b[6] = {1.0, 1.0, 1.0,  0.5, 0.0, 0.0  };
    
    double zmax;
    
    gStyle->SetPalette(1);
    //  int cc = TColor::CreateGradientColorTable(6, s, r, g, b, 50);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetStatColor(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.2);
    gStyle->SetTitleOffset(1.3,"Y");
    
    TH2D *hnd_pi = HEC->GetContHist(1,0);
    TH2D *hsk_pi_1rmu = HEC->GetContHist(0,0,0);
    TH2D *hsknd_pi_1rmu = HEC->GetContHist(2,0,0);
    zmax = max(hsk_pi_1rmu->GetMaximum(),hnd_pi->GetMaximum());
    hsk_pi_1rmu->SetMaximum(zmax);
    hnd_pi->SetMaximum(zmax);
    hsknd_pi_1rmu->SetMaximum(zmax);
    hsk_pi_1rmu->SetMinimum(-zmax);
    hnd_pi->SetMinimum(-zmax);
    hsknd_pi_1rmu->SetMinimum(-zmax);
    /*    
    cout << hsk_pi_1rmu->GetBinContent(221) << endl;
    cout << hnd_pi->GetBinContent(221) << endl;
    cout << hsknd_pi_1rmu->GetBinContent(221) << endl;
    */
    TCanvas *c11 = new TCanvas("c11","1rmu_pi",1000,700);
    c11->Divide(2,2);
    c11->cd(1);
    hsk_pi_1rmu->Draw("colz");
    c11->cd(2);
    hnd_pi->DrawClone("colz");
    c11->cd(3);
    hsknd_pi_1rmu->Draw("colz");
    
    TH2D *hnd_k = HEC->GetContHist(1,1);
    TH2D *hsk_k_1rmu = HEC->GetContHist(0,1,0);
    TH2D *hsknd_k_1rmu = HEC->GetContHist(2,1,0);
    zmax = max(hsk_k_1rmu->GetMaximum(),hnd_k->GetMaximum());
    hsk_k_1rmu->SetMaximum(zmax);
    hnd_k->SetMaximum(zmax);
    hsknd_k_1rmu->SetMaximum(zmax);
    hsk_k_1rmu->SetMinimum(-zmax);
    hnd_k->SetMinimum(-zmax);
    hsknd_k_1rmu->SetMinimum(-zmax);
    
    TCanvas *c12 = new TCanvas("c12","1rmu_k",1000,700);
    c12->Divide(2,2);
    c12->cd(1);
    hsk_k_1rmu->Draw("colz");
    c12->cd(2);
    hnd_k->DrawClone("colz");
    c12->cd(3);
    hsknd_k_1rmu->Draw("colz");
    
    TH2D *hsk_pi_nuetot = HEC->GetContHist(0,0,3);
    TH2D *hsknd_pi_nuetot = HEC->GetContHist(2,0,3);
    zmax = max(hsk_pi_nuetot->GetMaximum(),hnd_pi->GetMaximum());
    hsk_pi_nuetot->SetMaximum(zmax);
    hnd_pi->SetMaximum(zmax);
    hsknd_pi_nuetot->SetMaximum(zmax);
    hsk_pi_nuetot->SetMinimum(-zmax);
    hnd_pi->SetMinimum(-zmax);
    hsknd_pi_nuetot->SetMinimum(-zmax);
    
    TCanvas *c21 = new TCanvas("c21","nuetot_pi",1000,700);
    c21->Divide(2,2);
    c21->cd(1);
    hsk_pi_nuetot->Draw("colz");
    c21->cd(2);
    hnd_pi->DrawClone("colz");
    c21->cd(3);
    hsknd_pi_nuetot->Draw("colz");
    
    TH2D *hsk_k_nuetot = HEC->GetContHist(0,1,3);
    TH2D *hsknd_k_nuetot = HEC->GetContHist(2,1,3);
    zmax = max(hsk_k_nuetot->GetMaximum(),hnd_k->GetMaximum());
    hsk_k_nuetot->SetMaximum(zmax);
    hnd_k->SetMaximum(zmax);
    hsknd_k_nuetot->SetMaximum(zmax);
    hsk_k_nuetot->SetMinimum(-zmax);
    hnd_k->SetMinimum(-zmax);
    hsknd_k_nuetot->SetMinimum(-zmax);
    
    TCanvas *c22 = new TCanvas("c22","nuetot_k",1000,700);
    c22->Divide(2,2);
    c22->cd(1);
    hsk_k_nuetot->Draw("colz");
    c22->cd(2);
    hnd_k->DrawClone("colz");
    c22->cd(3);
    hsknd_k_nuetot->Draw("colz");

    app.Run();
  }
};
