//#define INPUT_DIR "/home/kubo1/beam_svn/analysis_code/BeamNormErrCalc/trunk/share"
#define INPUT_DIR "../share"

#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTimeStamp.h>
#include <iostream>

#include "HadErrCalc.h"

using namespace std;

int main(int argc, char * argv[]){
  
  TTimeStamp t1;

  HadErrCalc *HEC = new HadErrCalc(INPUT_DIR);
  
  /* read sample historgams */
  
  TFile * f = TFile::Open("../share/test_beaminput.root");
  TH1D *h1rnum = (TH1D*)f->Get("h_1rmu_num_sk");
  TH1D *h1rnmb = (TH1D*)f->Get("h_1rmu_nmb_sk");
  TH1D *h1rnue = (TH1D*)f->Get("h_1rmu_nue_sk");
  
  TH1D *hsignue = (TH1D*)f->Get("h_sig_nue_sk");
  
  TH1D *hbkgnum = (TH1D*)f->Get("h_bkg_num_sk");
  TH1D *hbkgnmb = (TH1D*)f->Get("h_bkg_nmb_sk");
  TH1D *hbkgnue = (TH1D*)f->Get("h_bkg_nue_sk");

  TH1D *hndnum = (TH1D*)f->Get("h_1trk_num_nd");
  TH1D *hndnmb = (TH1D*)f->Get("h_1trk_nmb_nd");
  TH1D *hndnue = (TH1D*)f->Get("h_1trk_nue_nd");

  TTimeStamp t2;
  cout << t2.GetSec() - t1.GetSec() << "sec for init." << endl;

  /* test N times with same histograms */
  double n = 10;

  for(int i=0;i<n;i++){
    /* set input histograms */
    HEC->SetNDHist(hndnum, hndnmb, hndnue);       /* required */
    HEC->Set1RmuHist(h1rnum, h1rnmb, h1rnue);     /* optional (for numu analysis) */
    HEC->SetNueSigHist(hsignue);                  /* optional (for nue analysis)  */
    HEC->SetNueBGHist(hbkgnum, hbkgnmb, hbkgnue); /* optional (for nue analysis)  */
    /* Exec. calculation */
    HEC->Calc();
    /* now, HadErrCalc has calculated results of the loop*/
    if(i%100 == 99)
    cout << "loop : " << i+1 << "   1RmuError = " << HEC->GetSumErr(0)*100 << "\%"<< endl;
  }
  
  TTimeStamp t3;
  cout << n <<" calls / " << ( t3.GetSec() - t2.GetSec() ) <<" sec." << endl;

  // print out the last results //
  cout <<  endl;
  cout << "NSK/NND errors -----------" << endl;
  cout << "1Rmu/ND error \t\t" << HEC->GetSumErr(0)*100 << " \%" << endl;
  cout << "NueSig/ND error\t\t" << HEC->GetSumErr(1)*100 << " \%"<< endl;
  cout << "NueBG/ND error \t\t" << HEC->GetSumErr(2)*100 << " \%" << endl;
  cout << "NueSig+BG/ND error\t" << HEC->GetSumErr(3)*100 << " \%" << endl;
  cout <<  endl;

  cout << "NSK,NND single errors --------" << endl;
  cout << "NND error \t" << HEC->GetNDErr()*100 << " \%" << endl;
  cout << "1Rmu error \t" << HEC->GetSKErr(0)*100 << " \%" << endl;
  cout << "NueSig error\t" << HEC->GetSKErr(1)*100 << " \%"<< endl;
  cout << "NueBG error \t" << HEC->GetSKErr(2)*100 << " \%" << endl;
  cout << "NueSig+BG error\t" << HEC->GetSKErr(3)*100 << " \%" << endl;

  cout <<  endl;
  // flavor by flavor //
  cout << "flavor by flavor --------" << endl;
  cout << "1Rmu/ND error numu\t" << HEC->GetErr(0,0)*100 << " \%"<< endl;
  cout << "1Rmu/ND error numub\t" << HEC->GetErr(0,1)*100 << " \%"<< endl;
  cout << "1Rmu/ND error nue\t" << HEC->GetErr(0,2)*100 << " \%"<< endl;
  cout << "NueBG/ND error numu\t" << HEC->GetErr(2,0)*100 << " \%"<< endl;
  cout << "NueBG/ND error numub\t" << HEC->GetErr(2,1)*100 << " \%"<< endl;
  cout << "NueBG/ND error nue\t" << HEC->GetErr(2,2)*100 << " \%"<< endl;
  cout <<  endl;
  
  // errors from each components //
  HEC->ErrorComponents();
  // verbose outputs //
  HEC->ShowDetail();

};
