#include <iostream>
#include <cmath>
#include <cstdlib>
#include "HadErrCalc.h"

HadErrCalc::HadErrCalc(const char * input_dir, bool ingrid){
  
  // read  p-theta error histograms
  TFile *f;
  
  TString dirname = input_dir;
  if(!ingrid){
    f = TFile::Open(dirname + "/"+DEF_ND5_INP);
  }
  else{
    f = TFile::Open(dirname + "/"+DEF_ING_INP);
  }

  
  if(!f){
    std::cout << "error database not found" << std::endl;
    std::cout << "please set correct input directory !! " << std::endl;
    exit(1);
  }
  
  TString pi_name[2] = {"pip","pim"};
  TString k_name[3] = {"kp","km","k0"};
  TString f_name[3] = {"numu","numub","nue"};
  TString hname;
  
  std::cout << "reading hadpro p-theta error table..." << std::endl;
  int PBIN, KBIN;
  
  for(int c=0;c<2;c++){
    if(!c){
      std::cout << "pi plus table " << std::endl;
      PBIN = PIPBIN;
    }
    else{
      std::cout << "pi minus table " << std::endl;
      PBIN = PIMBIN;
    }
    for(int m=0;m<3;m++){
      for(int i=0; i<PBIN; i++){
	TH1D *hsk;
	TH1D *hnd;
	hname.Form("hsk_%s_%s_%d",pi_name[c].Data(),f_name[m].Data(),i);
	hsk = (TH1D*)f->Get(hname);
	hname.Form("hnd_%s_%s_%d",pi_name[c].Data(),f_name[m].Data(),i);
	hnd = (TH1D*)f->Get(hname);
	for(int b=0;b<20;b++){
	  if(!c){
	    pip_econt_sk[i].data[m][b] = hsk->GetBinContent(b+1);
	    pip_econt_nd[i].data[m][b] = hnd->GetBinContent(b+1);
	  }
	  else{
	    pim_econt_sk[i].data[m][b] = hsk->GetBinContent(b+1);
	    pim_econt_nd[i].data[m][b] = hnd->GetBinContent(b+1);
	  }
	}
	delete hsk;
	delete hnd;
      }
    }
  }
  
  for(int c=0;c<3;c++){
    if(!c){
      std::cout << "K plus table " << std::endl;
      KBIN = KPBIN;
    }
    else if(c==1){
      std::cout << "K minus table " << std::endl;
      KBIN = KMBIN;
    }
    else{
      std::cout << "K0L table " << std::endl;
      KBIN = K0BIN;
    }
    for(int m=0;m<3;m++){
      for(int i=0; i<KBIN; i++){
	TH1D *hsk;
	TH1D *hnd;
	hname.Form("hsk_%s_%s_%d",k_name[c].Data(),f_name[m].Data(),i);
	hsk = (TH1D*)f->Get(hname);
	hname.Form("hnd_%s_%s_%d",k_name[c].Data(),f_name[m].Data(),i);
	hnd = (TH1D*)f->Get(hname);
	for(int b=0;b<20;b++){
	  if(!c){
	    kp_econt_sk[i].data[m][b] = hsk->GetBinContent(b+1);
	    kp_econt_nd[i].data[m][b] = hnd->GetBinContent(b+1);
	  }
	  else if(c==1){
	    km_econt_sk[i].data[m][b] = hsk->GetBinContent(b+1);
	    km_econt_nd[i].data[m][b] = hnd->GetBinContent(b+1);
	  }
	  else{
	    k0_econt_sk[i].data[m][b] = hsk->GetBinContent(b+1);
	    k0_econt_nd[i].data[m][b] = hnd->GetBinContent(b+1);
	  }
	}
	delete hsk;
	delete hnd;
      }
    }
  }
  
  std::cout << "reading other systematic sources..." << std::endl;
  
  if(!ingrid){
    nsys = 0;
    TFile *f_sec = TFile::Open(dirname+"/flux_xsec_1sigma.root");
    NewSys("had_xsec",f_sec,"xsec_1sigma");
    TFile *f_hcur = TFile::Open(dirname+"/horn_5kA_unc.root");
    NewSys("horn_current",f_hcur,"horn_5ka_unc");
    
    TFile *f_hali = TFile::Open(dirname+"/horn_misalign_unc.root");
    TFile *f_tali = TFile::Open(dirname+"/target_misalign_unc.root");
    NewSys("horn_align",f_hali,"horn_misalign_unc");
    NewSys("target_align",f_tali,"target_misalign_unc");
    
    // 2010/12/11, changed OA angle error size
    //  TFile *f_oa = TFile::Open(dirname+"/beam_1mrad_unc.root");
    //  NewSys("off-axis_angle",f_oa,"beam_1mrad_unc");  
    // 2011/5/23, changed OA angle error size
    TFile *f_oa = TFile::Open(dirname+"/beam_0.44mrad_unc.root");
    NewSys("off-axis_angle",f_oa,"beam_0.44mrad_unc");  
    //TFile *f_oa = TFile::Open(dirname+"/beam_0.4mrad_unc.root");
    //NewSys("off-axis_angle",f_oa,"beam_0.4mrad_unc");      

    // 2011/5/25, scaled pbeam error by 0.62/0.58
    //    TFile *f_pb = TFile::Open(dirname+"/pbeam_unc.root");
    TFile *f_pb = TFile::Open(dirname+"/pbeam_unc_run2.root");
    NewSys("proton beam",f_pb,"pbeam_unc");  

    TFile *f_pi_tert = TFile::Open(dirname+"/pion_tert_unc.root");
    NewSys("Pi x Scaling",f_pi_tert,"pitert_unc");

    //TFile *f_pi_bin = TFile::Open(dirname+"/pion_bin_unc.root");
    //NewSys("Pi Binning",f_pi_bin,"pibin_unc");

    TFile *f_k_norm = TFile::Open(dirname+"/k_norm_unc.root");
    NewSys("K norm 15%",f_k_norm,"k_norm_unc");

    TFile *f_sec_hixf2 = TFile::Open(dirname+"/sec_hixfp_2p0.root");
    NewSys("Hi xf proton",f_sec_hixf2,"sec_hixfp_unc");

    TFile *f_sec_nucl = TFile::Open(dirname+"/sec_nucl_prod_unc.root");
    NewSys("Sec. p,n prod.",f_sec_nucl,"sec_nucl_unc");

  }
  
  std::cout << "initializing..." << std::endl;
  Init();

  pibin = new NA61PiBin(dirname + "/na61_final.root");
  TFile *fk = TFile::Open(dirname + "/eichtenweight-fluka_rebin.root");
  hphase_k = (TH2D*)fk->Get("h_kp_xpt_weight");
  
  f->Close();
  std::cout << " done. " << std::endl;
};


void HadErrCalc::NewSys(const char * name, TFile *fsys, const char * basename){
  if(nsys>=MAX_NSYS){
    std::cerr << "too many systematic error sources: change MAX_NSYS !"  << std::endl;
    return;
  }
  TH1D *hsk_numu, *hsk_numub, *hsk_nue;
  TH1D *hnd_numu, *hnd_numub, *hnd_nue;
  TString base = basename;

  hsk_numu = (TH1D*)fsys->Get(base+"_sk_numu");
  hsk_numub = (TH1D*)fsys->Get(base+"_sk_numub");
  hsk_nue = (TH1D*)fsys->Get(base+"_sk_nue");
  hnd_numu = (TH1D*)fsys->Get(base+"_nd5_numu");
  hnd_numub = (TH1D*)fsys->Get(base+"_nd5_numub");
  hnd_nue = (TH1D*)fsys->Get(base+"_nd5_nue");
  
  Sys[nsys].SetName(name);
  Sys[nsys].SetSKHist(hsk_numu, hsk_numub, hsk_nue);
  Sys[nsys].SetNDHist(hnd_numu, hnd_numub, hnd_nue);
  nsys ++;
};

void HadErrCalc::CalcSKSum(){
  for(int i=0;i<4;i++){
    sumhaderr_pi_nsknd[i] = 0;
    sumhaderr_k_nsknd[i] = 0;
    sumhaderr_pi_nsk[i] = 0;
    sumhaderr_k_nsk[i] = 0;
  }
  int nflav = 3;
  double s = 0;
  double ss = 0;
  
  // pion +
  for(int j=0;j<PIPBIN;j++){ // pion-bin loop
    ss = 0;
    for(int mode=0;mode<3;mode++){
      if(mode==1) nflav = 1;
      else nflav = 3;
      s = 0;
      for(int m=0;m<nflav;m++){ // flavor loop
	s += err_pip_nsk[j][mode][m] * int_nsk[mode][m];
      }
      err_pip_nsk_all[j][mode] = s*inv_nsk[mode][3]; 
      err_pip_nsknd_all[j][mode] = fabs(err_pip_nsk_all[j][mode]-err_pip_nnd[j]);
      sumhaderr_pi_nsk[mode] += err_pip_nsk_all[j][mode];
      sumhaderr_pi_nsknd[mode] += err_pip_nsknd_all[j][mode];
      if(mode) ss += s;
      if(mode ==2){
	err_pip_nsk_all[j][3] = ss/(int_nsk[1][3]+int_nsk[2][3]);
	err_pip_nsknd_all[j][3] = fabs(err_pip_nsk_all[j][3]-err_pip_nnd[j]);
	sumhaderr_pi_nsk[3] += err_pip_nsk_all[j][3];
	sumhaderr_pi_nsknd[3] += err_pip_nsknd_all[j][3];
      }
    }
  }
  // pion -
  for(int j=0;j<PIMBIN;j++){ // pion-bin loop
    ss = 0;
    for(int mode=0;mode<3;mode++){
      if(mode==1) nflav = 1;
      else nflav = 3;
      s = 0;
      for(int m=0;m<nflav;m++){ // flavor loop
	s += err_pim_nsk[j][mode][m] * int_nsk[mode][m];
      }
      err_pim_nsk_all[j][mode] = s*inv_nsk[mode][3]; 
      err_pim_nsknd_all[j][mode] = fabs(err_pim_nsk_all[j][mode]-err_pim_nnd[j]);
      sumhaderr_pi_nsk[mode] += err_pim_nsk_all[j][mode];
      sumhaderr_pi_nsknd[mode] += err_pim_nsknd_all[j][mode];
      if(mode) ss += s;
      if(mode ==2){
	err_pim_nsk_all[j][3] = ss/(int_nsk[1][3]+int_nsk[2][3]);
	err_pim_nsknd_all[j][3] = fabs(err_pim_nsk_all[j][3]-err_pim_nnd[j]);
	sumhaderr_pi_nsk[3] += err_pim_nsk_all[j][3];
	sumhaderr_pi_nsknd[3] += err_pim_nsknd_all[j][3];
      }
    }
  }

  // kaon +
  for(int j=0;j<KPBIN;j++){ // kaon-bin loop
    ss = 0;
    for(int mode=0;mode<3;mode++){
      if(mode==1) nflav = 1;
      else nflav = 3;
      s = 0;
      for(int m=0;m<nflav;m++){ // flavor loop
	s += err_kp_nsk[j][mode][m] * int_nsk[mode][m];
      }
      err_kp_nsk_all[j][mode] = s*inv_nsk[mode][3]; 
      err_kp_nsknd_all[j][mode] = fabs(err_kp_nsk_all[j][mode]-err_kp_nnd[j]);
      sumhaderr_k_nsk[mode] += err_kp_nsk_all[j][mode];
      sumhaderr_k_nsknd[mode] += err_kp_nsknd_all[j][mode];
      if(mode) ss += s;
      if(mode ==2){
	err_kp_nsk_all[j][3] = ss/(int_nsk[1][3]+int_nsk[2][3]);
	err_kp_nsknd_all[j][3] = fabs(err_kp_nsk_all[j][3]-err_kp_nnd[j]);
	sumhaderr_k_nsk[3] += err_kp_nsk_all[j][3];
	sumhaderr_k_nsknd[3] += err_kp_nsknd_all[j][3];
      }
    }
  }

  // kaon -
  for(int j=0;j<KMBIN;j++){ // kaon-bin loop
    ss = 0;
    for(int mode=0;mode<3;mode++){
      if(mode==1) nflav = 1;
      else nflav = 3;
      s = 0;
      for(int m=0;m<nflav;m++){ // flavor loop
	s += err_km_nsk[j][mode][m] * int_nsk[mode][m];
      }
      err_km_nsk_all[j][mode] = s*inv_nsk[mode][3]; 
      err_km_nsknd_all[j][mode] = fabs(err_km_nsk_all[j][mode]-err_km_nnd[j]);
      sumhaderr_k_nsk[mode] += err_km_nsk_all[j][mode];
      sumhaderr_k_nsknd[mode] += err_km_nsknd_all[j][mode];
      if(mode) ss += s;
      if(mode ==2){
	err_km_nsk_all[j][3] = ss/(int_nsk[1][3]+int_nsk[2][3]);
	err_km_nsknd_all[j][3] = fabs(err_km_nsk_all[j][3]-err_km_nnd[j]);
	sumhaderr_k_nsk[3] += err_km_nsk_all[j][3];
	sumhaderr_k_nsknd[3] += err_km_nsknd_all[j][3];
      }
    }
  }

  // kaon 0
  for(int j=0;j<K0BIN;j++){ // kaon-bin loop
    ss = 0;
    for(int mode=0;mode<3;mode++){
      if(mode==1) nflav = 1;
      else nflav = 3;
      s = 0;
      for(int m=0;m<nflav;m++){ // flavor loop
	s += err_k0_nsk[j][mode][m] * int_nsk[mode][m];
      }
      err_k0_nsk_all[j][mode] = s*inv_nsk[mode][3]; 
      err_k0_nsknd_all[j][mode] = fabs(err_k0_nsk_all[j][mode]-err_k0_nnd[j]);
      sumhaderr_k_nsk[mode] += err_k0_nsk_all[j][mode];
      sumhaderr_k_nsknd[mode] += err_k0_nsknd_all[j][mode];
      if(mode) ss += s;
      if(mode ==2){
	err_k0_nsk_all[j][3] = ss/(int_nsk[1][3]+int_nsk[2][3]);
	err_k0_nsknd_all[j][3] = fabs(err_k0_nsk_all[j][3]-err_k0_nnd[j]);
	sumhaderr_k_nsk[3] += err_k0_nsk_all[j][3];
	sumhaderr_k_nsknd[3] += err_k0_nsknd_all[j][3];
      }
    }
  }
  
  for(int m=0;m<4;m++){ // mode
    AddNeutronErr(sumhaderr_pi_nsk[m]);
    AddNeutronErr(sumhaderr_pi_nsknd[m]);
    
    sumhaderr_nsknd[m] = 
      sqrt(sumhaderr_pi_nsknd[m]*sumhaderr_pi_nsknd[m]
	   +sumhaderr_k_nsknd[m]*sumhaderr_k_nsknd[m]);
    sumhaderr_nsk[m] = 
      sqrt(sumhaderr_pi_nsk[m]*sumhaderr_pi_nsk[m]
	   +sumhaderr_k_nsk[m]*sumhaderr_k_nsk[m]);
  }
};

void HadErrCalc::CalcSK(int mode){

  int nflav = 3;
  if(mode==1) nflav = 1;
  
  for(int i=0;i<3;i++){
    haderr_pi_nsk[mode][i] = 0;
    haderr_k_nsk[mode][i] = 0;
    haderr_pi_nsknd[mode][i] = 0;
    haderr_k_nsknd[mode][i] = 0;
  }

  // pion +
  for(int j=0;j<PIPBIN;j++){ // pion-bin loop
    for(int m=0;m<nflav;m++){ // flavor loop
      err_pip_nsk[j][mode][m] =0 ;
      err_pip_nsknd[j][mode][m] =0 ;
      for(int n=0;n<20;n++){ // enu loop
	err_pip_nsk[j][mode][m] += pip_econt_sk[j].data[m][n]*input_nsk[mode][m][n];
      }
      err_pip_nsk[j][mode][m] *= inv_nsk[mode][m];
      err_pip_nsknd[j][mode][m] = fabs(err_pip_nsk[j][mode][m]-err_pip_nnd[j]);
      haderr_pi_nsk[mode][m] += err_pip_nsk[j][mode][m];
      haderr_pi_nsknd[mode][m] += err_pip_nsknd[j][mode][m];
    }
  }
  
  // pion -
  for(int j=0;j<PIMBIN;j++){ // pion-bin loop
    for(int m=0;m<nflav;m++){ // flavor loop
      err_pim_nsk[j][mode][m] =0 ;
      err_pim_nsknd[j][mode][m] =0 ;
      for(int n=0;n<20;n++){ // enu loop
	err_pim_nsk[j][mode][m] += pim_econt_sk[j].data[m][n]*input_nsk[mode][m][n];
      }
      err_pim_nsk[j][mode][m] *= inv_nsk[mode][m];
      err_pim_nsknd[j][mode][m] = fabs(err_pim_nsk[j][mode][m]-err_pim_nnd[j]);
      haderr_pi_nsk[mode][m] += err_pim_nsk[j][mode][m];
      haderr_pi_nsknd[mode][m] += err_pim_nsknd[j][mode][m];
    }
  }
  
  // kaon +  
  for(int j=0;j<KPBIN;j++){ // kaon-bin loop
    for(int m=0;m<nflav;m++){ // flavor loop
      err_kp_nsk[j][mode][m] = 0;
      err_kp_nsknd[j][mode][m] =0 ;
      for(int n=0;n<20;n++){ // enu loop
	err_kp_nsk[j][mode][m] += kp_econt_sk[j].data[m][n]*input_nsk[mode][m][n];
      }
      err_kp_nsk[j][mode][m] *= inv_nsk[mode][m];
      err_kp_nsknd[j][mode][m] = fabs(err_kp_nsk[j][mode][m]-err_kp_nnd[j]);
      haderr_k_nsk[mode][m] += err_kp_nsk[j][mode][m];
      haderr_k_nsknd[mode][m] += err_kp_nsknd[j][mode][m];
    }
  }
  // kaon -
  for(int j=0;j<KMBIN;j++){ // kaon-bin loop
    for(int m=0;m<nflav;m++){ // flavor loop
      err_km_nsk[j][mode][m] = 0;
      err_km_nsknd[j][mode][m] =0 ;
      for(int n=0;n<20;n++){ // enu loop
	err_km_nsk[j][mode][m] += km_econt_sk[j].data[m][n]*input_nsk[mode][m][n];
      }
      err_km_nsk[j][mode][m] *= inv_nsk[mode][m];
      err_km_nsknd[j][mode][m] = fabs(err_km_nsk[j][mode][m]-err_km_nnd[j]);
      haderr_k_nsk[mode][m] += err_km_nsk[j][mode][m];
      haderr_k_nsknd[mode][m] += err_km_nsknd[j][mode][m];
    }
  }
  // kaon 0  
  for(int j=0;j<K0BIN;j++){ // kaon-bin loop
    for(int m=0;m<nflav;m++){ // flavor loop
      err_k0_nsk[j][mode][m] = 0;
      err_k0_nsknd[j][mode][m] =0 ;
      for(int n=0;n<20;n++){ // enu loop
	err_k0_nsk[j][mode][m] += k0_econt_sk[j].data[m][n]*input_nsk[mode][m][n];
      }
      err_k0_nsk[j][mode][m] *= inv_nsk[mode][m];
      err_k0_nsknd[j][mode][m] = fabs(err_k0_nsk[j][mode][m]-err_k0_nnd[j]);
      haderr_k_nsk[mode][m] += err_k0_nsk[j][mode][m];
      haderr_k_nsknd[mode][m] += err_k0_nsknd[j][mode][m];
    }
  }
  
  for(int m=0;m<nflav;m++){ // flavor loop
    AddNeutronErr(haderr_pi_nsk[mode][m]);
    AddNeutronErr(haderr_pi_nsknd[mode][m]);
    
    haderr_nsknd[mode][m] = 
      sqrt(haderr_pi_nsknd[mode][m]*haderr_pi_nsknd[mode][m]
	   +haderr_k_nsknd[mode][m]*haderr_k_nsknd[mode][m]);
  }
};

void HadErrCalc::CalcND(){

  haderr_pi_nnd = 0;
  haderr_k_nnd = 0;

  // pion+
  for(int j=0;j<PIPBIN;j++){ // pion-bin loop
    err_pip_nnd[j] =0 ;
    for(int m=0;m<3;m++){ // flavor loop
      for(int n=0;n<20;n++){ // enu loop
	err_pip_nnd[j] += pip_econt_nd[j].data[m][n]*input_nnd[m][n];
      }
    }
    err_pip_nnd[j] *= inv_nnd[3];
    haderr_pi_nnd += err_pip_nnd[j];
  }
  // pion-
  for(int j=0;j<PIMBIN;j++){ // pion-bin loop
    err_pim_nnd[j] =0 ;
    for(int m=0;m<3;m++){ // flavor loop
      for(int n=0;n<20;n++){ // enu loop
	err_pim_nnd[j] += pim_econt_nd[j].data[m][n]*input_nnd[m][n];
      }
    }
    err_pim_nnd[j] *= inv_nnd[3];
    haderr_pi_nnd += err_pim_nnd[j];
  }
  
  // kaon+
  for(int j=0;j<KPBIN;j++){ // kaon-bin loop
    err_kp_nnd[j] = 0;
    for(int m=0;m<3;m++){ // flavor loop
      for(int n=0;n<20;n++){ // enu loop
	err_kp_nnd[j] += kp_econt_nd[j].data[m][n]*input_nnd[m][n];
      }
    }
    err_kp_nnd[j] *= inv_nnd[3];
    haderr_k_nnd += err_kp_nnd[j];
  }
  // kaon-
  for(int j=0;j<KMBIN;j++){ // kaon-bin loop
    err_km_nnd[j] = 0;
    for(int m=0;m<3;m++){ // flavor loop
      for(int n=0;n<20;n++){ // enu loop
	err_km_nnd[j] += km_econt_nd[j].data[m][n]*input_nnd[m][n];
      }
    }
    err_km_nnd[j] *= inv_nnd[3];
    haderr_k_nnd += err_km_nnd[j];
  }
  // kaon+
  for(int j=0;j<K0BIN;j++){ // kaon-bin loop
    err_k0_nnd[j] = 0;
    for(int m=0;m<3;m++){ // flavor loop
      for(int n=0;n<20;n++){ // enu loop
	err_k0_nnd[j] += k0_econt_nd[j].data[m][n]*input_nnd[m][n];
      }
    }
    err_k0_nnd[j] *= inv_nnd[3];
    haderr_k_nnd += err_k0_nnd[j];
  }
  AddNeutronErr(haderr_pi_nnd);
};

void HadErrCalc::CalcSys(int id){
  for(int i=0;i<4;i++){
    Sys[id].sumsyserr_nsknd[i] = 0;
    Sys[id].sumsyserr_nsk[i] = 0;
  }
  Sys[id].syserr_nnd = 0;
  for(int m=0;m<3;m++){ // flavor loop
    for(int n=0;n<20;n++){ // enu loop
      Sys[id].syserr_nnd += input_nnd[m][n]*Sys[id].herrnd[m]->GetBinContent(n+1);
    }
  }
  Sys[id].syserr_nnd *= inv_nnd[3];
  
  for(int mode=0;mode<3;mode++){ // mode loop
    int nflav = 3;
    if(mode==1) nflav = 1;
    for(int m=0;m<nflav;m++){ // flavor loop
      Sys[id].syserr_nsk[mode][m] = 0;
      Sys[id].syserr_nsknd[mode][m] = 0;
      for(int n=0;n<20;n++){ // enu loop
	Sys[id].syserr_nsk[mode][m] += input_nsk[mode][m][n]*Sys[id].herrsk[m]->GetBinContent(n+1);
      }
      Sys[id].sumsyserr_nsk[mode] += Sys[id].syserr_nsk[mode][m];
      Sys[id].syserr_nsk[mode][m] *= inv_nsk[mode][m];
      Sys[id].syserr_nsknd[mode][m] = (1+Sys[id].syserr_nsk[mode][m])/(1+Sys[id].syserr_nnd)-1;
    }
    if(mode) Sys[id].sumsyserr_nsk[3] += Sys[id].sumsyserr_nsk[mode]; 
    Sys[id].sumsyserr_nsk[mode] *= inv_nsk[mode][3];
    Sys[id].sumsyserr_nsknd[mode] = (1+Sys[id].sumsyserr_nsk[mode])/(1+Sys[id].syserr_nnd)-1;
  }
  Sys[id].sumsyserr_nsk[3] /= (int_nsk[1][3]+int_nsk[2][3]);
  Sys[id].sumsyserr_nsknd[3] = (1+Sys[id].sumsyserr_nsk[3])/(1+Sys[id].syserr_nnd)-1;
};


void HadErrCalc::Calc(){
  if(!fFlag[0]){
    std::cout <<" ND histograms are not given in this cycle !" <<std::endl ;
    return;
  }

  for(int j=0;j<4;j++){
    if(int_nnd[j]>0) inv_nnd[j] = 1./int_nnd[j];
    else inv_nnd[j] = 0;
    for(int i=0;i<3;i++){
      if(int_nsk[i][j]>0) inv_nsk[i][j] = 1./int_nsk[i][j];
      else inv_nsk[i][j] = 0;
    }
  }
  
  CalcND();
  if(fFlag[1]) CalcSK(0);
  if(fFlag[2]) CalcSK(1);
  if(fFlag[3]) CalcSK(2);
  
  CalcSKSum();

  for(int i=0;i<nsys;i++){
    CalcSys(i);
  }
  
  // ND grand total error
  ggtoterr_nnd = haderr_pi_nnd * haderr_pi_nnd + haderr_k_nnd * haderr_k_nnd;
  for(int i=0;i<nsys;i++){
    ggtoterr_nnd += Sys[i].syserr_nnd * Sys[i].syserr_nnd;
  }
  ggtoterr_nnd = sqrt(ggtoterr_nnd);
  
  // SK/ND  total error
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      gtoterr_nsknd[i][j] = 
	haderr_nsknd[i][j] * haderr_nsknd[i][j];
      for(int k=0;k<nsys;k++){
	gtoterr_nsknd[i][j] +=
	  Sys[k].syserr_nsknd[i][j]*Sys[k].syserr_nsknd[i][j];
      }
      gtoterr_nsknd[i][j] = sqrt(gtoterr_nsknd[i][j]);
    }
  }
  // SK/ND grand total error
  for(int i=0;i<4;i++){
    ggtoterr_nsknd[i] = 
      sumhaderr_nsknd[i] * sumhaderr_nsknd[i];
    ggtoterr_nsk[i] = 
      sumhaderr_nsk[i] * sumhaderr_nsk[i];
    for(int k=0;k<nsys;k++){
      ggtoterr_nsknd[i] +=
	Sys[k].sumsyserr_nsknd[i]*Sys[k].sumsyserr_nsknd[i];
      ggtoterr_nsk[i] +=
	Sys[k].sumsyserr_nsk[i]*Sys[k].sumsyserr_nsk[i];
    }
    ggtoterr_nsknd[i] = sqrt(ggtoterr_nsknd[i]);
    ggtoterr_nsk[i] = sqrt(ggtoterr_nsk[i]);
  }
  
  ResetFlag();
};

void HadErrCalc::ResetFlag(){
  for(int i=0;i<4;i++){
    fFlag[i] = 0;
  }
};

void HadErrCalc::Clear(){
};

void HadErrCalc::Init(){
  for(int i=0;i<3;i++){
    int_nnd[i] = 0;
    for(int j=0;j<3;j++){
      int_nsk[i][j] = 0;
    }
  }
  for(int b=0;b<20;b++){
    for(int i=0;i<3;i++){
      input_nnd[i][b] = 0;
      for(int j=0;j<3;j++){
	input_nsk[i][j][b] = 0;
      }
    }
  }
};

int HadErrCalc::ChkHist(TH1D *h){
  if(h->GetNbinsX() != 200 ||
     h->GetBinLowEdge(1) != 0 ||
     h->GetBinLowEdge(h->GetNbinsX()+1) != 10){
    std::cout << "input histogram should be [0-10] GeV, 200bins !" << std::endl;
    std::cout << "your input was " << h->GetBinLowEdge(1) << " - " <<
      h->GetBinLowEdge(h->GetNbinsX()+1) << " GeV, " <<
      h->GetNbinsX() << "bins" << std::endl;
    return -1;
  }
  return 0;
};

void HadErrCalc::SetHist(int mode, int flav, TH1D *h){
  double *input;
  if(mode == 3) input = input_nnd[flav];
  else input = input_nsk[mode][flav];

  for(int i=0;i<20;i++){
    input[i] = 0;
  }
  
  if(h){
    for(int i=1;i<202;i++){
      input[GetEnuBin(i)-1] += h->GetBinContent(i);
    }
  }
};

void HadErrCalc::SetNueSigHist(TH1D *h_nue_sig){
  if(ChkHist(h_nue_sig)) return;
  
  if(h_nue_sig){
    SetHist(1,0,h_nue_sig);
    int_nsk[1][0] = h_nue_sig->Integral(1,201);
    int_nsk[1][3] = h_nue_sig->Integral(1,201);
  }
  SetHist(1,1);
  SetHist(1,2);
  if(int_nsk[1][3]>0) fFlag[2]=1;
};

void HadErrCalc::SetNueBGHist(TH1D *h_nue_bg_numu, TH1D *h_nue_bg_numub, TH1D *h_nue_bg_nue){
  if(ChkHist(h_nue_bg_numu)) return;
  if(ChkHist(h_nue_bg_numub)) return;
  if(ChkHist(h_nue_bg_nue)) return;

  int_nsk[2][3] = 0;
  if(h_nue_bg_numu){
    SetHist(2,0,h_nue_bg_numu);
    int_nsk[2][0] = h_nue_bg_numu->Integral(1,201);
    int_nsk[2][3] += h_nue_bg_numu->Integral(1,201);
  }
  if(h_nue_bg_numub){
    SetHist(2,1,h_nue_bg_numub);
    int_nsk[2][1] = h_nue_bg_numub->Integral(1,201);
    int_nsk[2][3] += h_nue_bg_numub->Integral(1,201);
  }
  if(h_nue_bg_nue){
    SetHist(2,2,h_nue_bg_nue);
    int_nsk[2][2] = h_nue_bg_nue->Integral(1,201);
    int_nsk[2][3] += h_nue_bg_nue->Integral(1,201);
  }

  if(int_nsk[2][3]>0) fFlag[3]=1;
};

void HadErrCalc::Set1RmuHist(TH1D *h_1rmu_numu, TH1D *h_1rmu_numub, TH1D *h_1rmu_nue){
  if(ChkHist(h_1rmu_numu)) return;
  if(ChkHist(h_1rmu_numub)) return;
  if(ChkHist(h_1rmu_nue)) return;
  int_nsk[0][3] = 0;
  
  if(h_1rmu_numu){
    SetHist(0,0,h_1rmu_numu);
    int_nsk[0][0] = h_1rmu_numu->Integral(1,201);
    int_nsk[0][3] += h_1rmu_numu->Integral(1,201);
  }
  if(h_1rmu_numub){
    SetHist(0,1,h_1rmu_numub);
    int_nsk[0][1] = h_1rmu_numub->Integral(1,201);
    int_nsk[0][3] += h_1rmu_numub->Integral(1,201);
  }
  if(h_1rmu_nue){
    SetHist(0,2,h_1rmu_nue);
    int_nsk[0][2] = h_1rmu_nue->Integral(1,201);
    int_nsk[0][3] += h_1rmu_nue->Integral(1,201);
  }
  if(int_nsk[0][3]>0) fFlag[1]=1;

};


void HadErrCalc::SetNDHist(TH1D *h_nd_numu, TH1D *h_nd_numub, TH1D *h_nd_nue){
  int_nnd[3] = 0;  
  if(ChkHist(h_nd_numu)) return;
  if(ChkHist(h_nd_numub)) return;
  if(ChkHist(h_nd_nue)) return;
  
  if(h_nd_numu){
    SetHist(3,0,h_nd_numu);
    int_nnd[0] = h_nd_numu->Integral(1,201);
    int_nnd[3] += h_nd_numu->Integral(1,201);
  }
  if(h_nd_numub){
    SetHist(3,1,h_nd_numub);
    int_nnd[1] = h_nd_numub->Integral(1,201);
    int_nnd[3] += h_nd_numub->Integral(1,201);
  }
  if(h_nd_nue){
    SetHist(3,2,h_nd_nue);
    int_nnd[2] = h_nd_nue->Integral(1,201);
    int_nnd[3] += h_nd_nue->Integral(1,201);
  }
  if(int_nnd[3]>0) fFlag[0]=1;
};

