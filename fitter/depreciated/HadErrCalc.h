#define DEF_ND5_INP "delta_flux_nd5_v3.root"
#define DEF_ING_INP "delta_flux_ing_v3.root"
#define MAX_NSYS 10 // maximum # of sys. error terms

#define PIPBIN 309   // Pion+ phase space # of bins
#define PIMBIN 361   // Pion- phase space # of bins
#define KPBIN 700    // Kaon+ phase space # of bins
#define KMBIN 700    // Kaon- phase space # of bins
#define K0BIN 700    // Kaon0 phase space # of bins

//#define NEUTRON_ERR 0.0015 // neutron production error size
//#define ISOSPIN_ERR 0.0045 // uncertainty on isospin symmetry
#define NEUTRON_ERR 0. // neutron production error size
#define ISOSPIN_ERR 0. // uncertainty on isospin symmetry

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TString.h>
#include "na61_pibin.h"

class BinErr{
 public:
  double data[3][20]; // [flavor],[enu]
};

class SysErr{
  // class for individual sys. error sources treated by 1sigma covariances
 private:
  TString fName;
 public:
  TH1D *herrsk[3];
  TH1D *herrnd[3];
  
  double syserr_nnd;
  double syserr_nsk[3][3]; // [mode][flavor]
  double syserr_nsknd[3][3]; // [mode][flavor]
  
  double sumsyserr_nsk[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  double sumsyserr_nsknd[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  
  SysErr(){  };
  
  void SetName(const char * name){
    fName = name;
  };
  const char * GetName(){return fName.Data();};
  
  void SetSKHist(TH1D *hsk_numu, TH1D *hsk_numub, TH1D *hsk_nue){
    herrsk[0] = hsk_numu;
    herrsk[1] = hsk_numub;
    herrsk[2] = hsk_nue;
  };
  void SetNDHist(TH1D *hnd_numu, TH1D *hnd_numub, TH1D *hnd_nue){
    herrnd[0] = hnd_numu;
    herrnd[1] = hnd_numub;
    herrnd[2] = hnd_nue;
  };

  void Show(){
    std::cout <<"\n-------------"<< fName << " error ---------------" << std::endl;
    std::cout << "[before cancellation]" << std::endl;
    std::cout << "\t\ttotal\tnumu\tnumub\tnue" << std::endl;
    printf("ND: \t\t%.2f%%\n",syserr_nnd * 100);
    printf("SK 1Rmu: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumsyserr_nsk[0]*100,syserr_nsk[0][0]*100,syserr_nsk[0][1]*100,syserr_nsk[0][2]*100);
    printf("SK NueSig: \t%.2f%%\n",syserr_nsk[1][0]*100);
    printf("SK NueBG: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumsyserr_nsk[2]*100,syserr_nsk[2][0]*100,syserr_nsk[2][1]*100,syserr_nsk[2][2]*100);
    std::cout << "[after cancellation]" << std::endl;
    std::cout << "\t\ttotal\tnumu\tnumub\tnue" << std::endl;
    printf("1Rmu/ND: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumsyserr_nsknd[0]*100,syserr_nsknd[0][0]*100,syserr_nsknd[0][1]*100,syserr_nsknd[0][2]*100);
    printf("NueSig/ND: \t%.2f%%\n",syserr_nsknd[1][0]*100);
    printf("NueBG/ND: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumsyserr_nsknd[2]*100,syserr_nsknd[2][0]*100,syserr_nsknd[2][1]*100,syserr_nsknd[2][2]*100);
    printf("NueTot/ND: \t%.2f%%\n",sumsyserr_nsknd[3]*100);
  };
};

class HadErrCalc{
 private:
  int nsys;
  SysErr Sys[MAX_NSYS]; // sys. errors other than  pi/K
  
  // error contribution inputs
  BinErr pip_econt_sk[PIPBIN]; // [pion+]
  BinErr pim_econt_sk[PIMBIN]; // [pion-]
  BinErr kp_econt_sk[KPBIN];  // [kaon+]
  BinErr km_econt_sk[KMBIN];  // [kaon-]
  BinErr k0_econt_sk[K0BIN];  // [kaon0]
  BinErr pip_econt_nd[PIPBIN]; // [pion+]
  BinErr pim_econt_nd[PIMBIN]; // [pion-]
  BinErr kp_econt_nd[KPBIN];  // [kaon+]
  BinErr km_econt_nd[KMBIN];  // [kaon-]
  BinErr k0_econt_nd[K0BIN];  // [kaon0]
  
  bool fFlag[4]; // [nd, 1rmu, nue_sig, nue_bkg]
  
  // original NSK/NND (integral of input histograms)
  double int_nsk[3][4]; // [1rmu, nue_sig, nue_bkg][numu,numub,nue,all]
  double int_nnd[4];    // [numu, numub, nue, all]
  double inv_nsk[3][4]; // inverse
  double inv_nnd[4];    // inverse
  
  // input histogram contents
  double input_nsk[3][3][20]; // [1rmu, nue_sig, nue_bkg][numu,numub,nue][enu]
  double input_nnd[3][20];    // [numu, numub, nue][enu]

  // error contribution outputs
  double err_pip_nsk[PIPBIN][3][3]; // [pion+],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_pim_nsk[PIMBIN][3][3]; // [pion-],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_kp_nsk[KPBIN][3][3]; // [kaon+],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_km_nsk[KMBIN][3][3]; // [kaon-],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_k0_nsk[K0BIN][3][3]; // [kaon0],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  
  double err_pip_nnd[PIPBIN]; // [pion+],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_pim_nnd[PIMBIN]; // [pion-],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_kp_nnd[KPBIN]; // [kaon+],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_km_nnd[KMBIN]; // [kaon-],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_k0_nnd[K0BIN]; // [kaon0],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]

  // Nsk/Nnd error outputs
  double err_pip_nsknd[PIPBIN][3][3]; // [pion+],[1rmu, nue_sig, nue_bkg][numu,numub,nue]
  double err_pim_nsknd[PIMBIN][3][3]; // [pion+],[1rmu, nue_sig, nue_bkg][numu,numub,nue]
  double err_kp_nsknd[KPBIN][3][3]; // [kaon],[1rmu, nue_sig, nue_bkg][numu,numub,nue]
  double err_km_nsknd[KMBIN][3][3]; // [kaon],[1rmu, nue_sig, nue_bkg][numu,numub,nue]
  double err_k0_nsknd[K0BIN][3][3]; // [kaon],[1rmu, nue_sig, nue_bkg][numu,numub,nue]
  
  
  // error contribution histograms //
  TH2D * hphase_k;
  NA61PiBin *pibin;
  
  TH2D * hcontnd[3]; //{pip, pim, k}
  TH2D * hcontsk[3][4]; // {pip,pim.k}{1rmu, nue_sig, nue_bg, nue_all}
  TH2D * hcontsknd[3][4]; // {pip,pim,k}{1rmu, nue_sig, nue_bg, nue_all}
  
  // sum of all flavors
  double err_pip_nsk_all[PIPBIN][4];
  double err_pip_nsknd_all[PIPBIN][4];
  double err_pim_nsk_all[PIMBIN][4];
  double err_pim_nsknd_all[PIMBIN][4];
  double err_kp_nsk_all[KPBIN][4];
  double err_kp_nsknd_all[KPBIN][4];
  double err_km_nsk_all[KMBIN][4];
  double err_km_nsknd_all[KMBIN][4];
  double err_k0_nsk_all[K0BIN][4];
  double err_k0_nsknd_all[K0BIN][4];

  // hadpro total error (mode by mode, flavor by flavor)
  double haderr_pi_nnd;
  double haderr_k_nnd;
  double haderr_pi_nsk[3][3];
  double haderr_k_nsk[3][3];
  double haderr_pi_nsknd[3][3];
  double haderr_k_nsknd[3][3];
  double haderr_nsknd[3][3];
  
  // hadpro total error (mode by mode)
  double sumhaderr_pi_nsk[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  double sumhaderr_k_nsk[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  double sumhaderr_nsk[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  double sumhaderr_pi_nsknd[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  double sumhaderr_k_nsknd[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  double sumhaderr_nsknd[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  
  // beam grand total error
  double gtoterr_nsknd[3][3]; // [mode][flavor]
  double ggtoterr_nsknd[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  double ggtoterr_nsk[4]; // [1rmu, nue_sig_ nue_bkg, nue_sig+bkg]
  double ggtoterr_nnd;

  void CalcSys(int id);
  void CalcSK(int mode);
  void CalcSKSum();
  void CalcND();
  void Init();
  void NewSys(const char * name, TFile *fsys, const char * basename);
  void Clear();
  void ResetFlag();
  int ChkHist(TH1D *h);
  void SetHist(int mode, int flav, TH1D *h = 0);

  void AddNeutronErr(double &val){
    //    val += NEUTRON_ERR;
    //    val = TMath::Sqrt(val*val + ISOSPIN_ERR*ISOSPIN_ERR);
    val *= 1.05;
  }
  
 public:
  
  // prototypes
  HadErrCalc(const char *input_dir = "../input/", bool ingrid = 0);
  void Calc();
  void SetNueSigHist(TH1D *h_nue_sig);
  void SetNueBGHist(TH1D *h_nue_bg_numu, TH1D *h_nue_bg_numub, TH1D *h_nue_bg_nue);
  void Set1RmuHist(TH1D *h_1rmu_numu, TH1D *h_1rmu_numub, TH1D *h_1rmu_nue);
  void SetNDHist(TH1D *h_nd_numu, TH1D *h_nd_numub, TH1D *h_nd_nue);

  // inlines, getters
  ~HadErrCalc(){};
  
  double GetErr(int mode, int flav){ // mode = {1rmu, nuesig, nueBG} // flav = {numu, numub, nue}
    return gtoterr_nsknd[mode][flav];
  };

  double GetSumErr(int mode){ // mode = {1rmu, nuesig, nueBG, nueTot}
    return ggtoterr_nsknd[mode];
  };

  double GetSKErr(int mode){
    return ggtoterr_nsk[mode];
  };

  double GetNDErr(){
    return ggtoterr_nnd;
  };

  void ShowPIPErr(int mode){
    double sk =0, nd=0, sknd=0;
    for(int i=0;i<PIPBIN;i++){
      sk += err_pip_nsk_all[i][mode];
      sknd += err_pip_nsknd_all[i][mode];
      nd += err_pip_nnd[i];
    }
    std::cout << sk << " " << nd <<  " " << sknd << std::endl;
  };
  void ShowPIMErr(int mode){
    double sk =0, nd=0, sknd=0;
    for(int i=0;i<PIMBIN;i++){
      sk += err_pim_nsk_all[i][mode];
      sknd += err_pim_nsknd_all[i][mode];
      nd += err_pim_nnd[i];
    }
    std::cout << sk << " " << nd <<  " " << sknd << std::endl;
  };

  TH2D * GetContHist(int det, int pik, int mode = 0){
    TString s_pik[3] = {"pip","pim","k"};
    TString s_mode[4] = {"1Rmu", "NueSig", "NueBG", "NueAll"};
    TString s_x[2] = {"P_{#pi} [GeV/c]", "X_{F}"};
    TString s_y[2] = {"#theta [mrad]", "P_{T} [GeV/c]"};
    int nbins[3] = {PIPBIN, PIMBIN, KPBIN};
    TH2D *h;
    if(det == 1){//ND
      h = hcontnd[pik];
      if(!h){
	TString name = "h_nd_" + s_pik[pik];
	if(pik ==2)
	  h = (TH2D*)hphase_k->Clone(name);
	else
	  h = new TH2D(name,name,201,0,20.1,22,0,440);
      }
      
      if(!pik){ // pip
	for(int x =0;x<201;x++){
	  for(int y=0;y<22;y++){
	    double cx = h->GetXaxis()->GetBinCenter(x+1);
	    double cy = h->GetYaxis()->GetBinCenter(y+1);
	    h->SetBinContent(x+1, y+1, err_pip_nnd[pibin->get_binID_pip(cy,cx)]);
	  }
	}
      }
      else if(pik==1){ // pim
	for(int x =0;x<201;x++){
	  for(int y=0;y<22;y++){
	    double cx = h->GetXaxis()->GetBinCenter(x+1);
	    double cy = h->GetYaxis()->GetBinCenter(y+1);
	    h->SetBinContent(x+1, y+1, err_pim_nnd[pibin->get_binID_pim(cy,cx)]);
	  }
	}
      }
      else
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_kp_nnd[i] + err_km_nnd[i] + err_k0_nnd[i]);
    }
    
    if(det == 0){//SK
      h = hcontsk[pik][mode];
      if(!h){
	TString name = "h_sk_" + s_mode[mode]+ "_" + s_pik[pik]; 
	if(pik ==2)
	  h = (TH2D*)hphase_k->Clone(name);
	else
	  h = new TH2D(name,name,201,0,20.1,22,0,440);
      }
      if(!pik){ // pip
	for(int x =0;x<201;x++){
	  for(int y=0;y<22;y++){
	    double cx = h->GetXaxis()->GetBinCenter(x+1);
	    double cy = h->GetYaxis()->GetBinCenter(y+1);
	    h->SetBinContent(x+1, y+1, err_pip_nsk_all[pibin->get_binID_pip(cy,cx)][mode]);
	  }
	}
      }
      else if(pik==1){ // pim
	for(int x =0;x<201;x++){
	  for(int y=0;y<22;y++){
	    double cx = h->GetXaxis()->GetBinCenter(x+1);
	    double cy = h->GetYaxis()->GetBinCenter(y+1);
	    h->SetBinContent(x+1, y+1, err_pim_nsk_all[pibin->get_binID_pim(cy,cx)][mode]);
	  }
	}
      }
      else
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_kp_nsk_all[i][mode] + 
						       err_km_nsk_all[i][mode] + err_k0_nsk_all[i][mode]);
    }
    
    if(det == 2){//SK/ND
      h = hcontsknd[pik][mode];
      if(!h){
	TString name = "h_sk_" + s_mode[mode]+ "_" + s_pik[pik]+"/nd"; 
	if(pik ==2)
	  h = (TH2D*)hphase_k->Clone(name);
	else
	  h = new TH2D(name,name,201,0,20.1,22,0,440);
      }
      if(!pik){ // pip
	for(int x =0;x<201;x++){
	  for(int y=0;y<22;y++){
	    double cx = h->GetXaxis()->GetBinCenter(x+1);
	    double cy = h->GetYaxis()->GetBinCenter(y+1);
	    h->SetBinContent(x+1, y+1, err_pip_nsk_all[pibin->get_binID_pip(cy,cx)][mode]
			     - err_pip_nnd[pibin->get_binID_pip(cy,cx)]);
	  }
	}
      }
      else if(pik==1){ // pim
	for(int x =0;x<201;x++){
	  for(int y=0;y<22;y++){
	    double cx = h->GetXaxis()->GetBinCenter(x+1);
	    double cy = h->GetYaxis()->GetBinCenter(y+1);
	    h->SetBinContent(x+1, y+1, err_pim_nsk_all[pibin->get_binID_pim(cy,cx)][mode]
			     - err_pim_nnd[pibin->get_binID_pim(cy,cx)]);
	  }
	}
      }
      else
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_kp_nsk_all[i][mode] - err_kp_nnd[i]
						       +err_km_nsk_all[i][mode] - err_km_nnd[i]
						       +err_k0_nsk_all[i][mode] - err_k0_nnd[i]);
    }
    
    h->SetLabelSize(0.06,"XYZ");
    h->SetNdivisions(505,"XYZ");
    h->GetXaxis()->SetTitle(s_x[pik]);
    h->GetYaxis()->SetTitle(s_y[pik]);
    h->SetTitleSize(0.06 ,"XYZ");
    return h;
  };
  
  void ShowSys(){
    for(int i=0;i<nsys;i++){
      Sys[i].Show();
    }
  };

  void ShowHad(){
    std::cout <<"\n---------- Pion production error ----------" << std::endl;
    std::cout << "[befor cancellation]" << std::endl;
    std::cout << "\t\ttotal\tnumu\tnumub\tnue" << std::endl;
    printf("ND: \t\t%.2f%% \n",haderr_pi_nnd * 100);
    printf("SK 1Rmu: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumhaderr_pi_nsk[0]*100,haderr_pi_nsk[0][0]*100,haderr_pi_nsk[0][1]*100,haderr_pi_nsk[0][2]*100);
    printf("SK NueSig: \t%.2f%%\n",haderr_pi_nsk[1][0]*100);
    printf("SK NueBG: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumhaderr_pi_nsk[2]*100,haderr_pi_nsk[2][0]*100,haderr_pi_nsk[2][1]*100,haderr_pi_nsk[2][2]*100);
    std::cout << "[after cancellation]" << std::endl;
    std::cout << "\t\ttotal\tnumu\tnumub\tnue" << std::endl;
    printf("1Rmu/ND: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumhaderr_pi_nsknd[0]*100,haderr_pi_nsknd[0][0]*100,haderr_pi_nsknd[0][1]*100,haderr_pi_nsknd[0][2]*100);
    printf("NueSig/ND: \t%.2f%%\n",haderr_pi_nsknd[1][0]*100);
    printf("NueBG/ND: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumhaderr_pi_nsknd[2]*100,haderr_pi_nsknd[2][0]*100,haderr_pi_nsknd[2][1]*100,haderr_pi_nsknd[2][2]*100);
    printf("NueTot/ND: \t%.2f%%\n",sumhaderr_pi_nsknd[3]*100);
    
    std::cout <<"\n---------- Kaon production error ----------" << std::endl;
    std::cout << "[befor cancellation]" << std::endl;
    std::cout << "\t\ttotal\tnumu\tnumub\tnue" << std::endl;
    printf("ND: \t\t%.2f%% \n",haderr_k_nnd * 100);
    printf("SK 1Rmu: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumhaderr_k_nsk[0]*100,haderr_k_nsk[0][0]*100,haderr_k_nsk[0][1]*100,haderr_k_nsk[0][2]*100);
    printf("SK NueSig: \t%.2f%%\n",haderr_k_nsk[1][0]*100);
    printf("SK NueBG: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumhaderr_k_nsk[2]*100,haderr_k_nsk[2][0]*100,haderr_k_nsk[2][1]*100,haderr_k_nsk[2][2]*100);
    std::cout << "[after cancellation]" << std::endl;
    std::cout << "\t\ttotal\tnumu\tnumub\tnue" << std::endl;
    printf("1Rmu/ND: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumhaderr_k_nsknd[0]*100,haderr_k_nsknd[0][0]*100,haderr_k_nsknd[0][1]*100,haderr_k_nsknd[0][2]*100);
    printf("NueSig/ND: \t%.2f%%\n",haderr_k_nsknd[1][0]*100);
    printf("NueBG/ND: \t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",sumhaderr_k_nsknd[2]*100,haderr_k_nsknd[2][0]*100,haderr_k_nsknd[2][1]*100,haderr_k_nsknd[2][2]*100);
    printf("NueTot/ND: \t%.2f%%\n",sumhaderr_k_nsknd[3]*100);
  };

  void ErrorComponents(){
    std::cout << "\n------------------- Error components --------------------------------------" << std::endl;
    std::cout << "\t\t1Rmu/ND  \tNueSig/ND\tNueBG/ND\tNueTot/ND" << std::endl;
    printf("Total \t\t %.2f%%\t\t %.2f%%\t\t %.2f%%\t\t %.2f%%\n",ggtoterr_nsknd[0]*100, ggtoterr_nsknd[1]*100, ggtoterr_nsknd[2]*100, ggtoterr_nsknd[3]*100);
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    printf("Pion Prod.\t %.2f%%\t\t %.2f%%\t\t %.2f%%\t\t %.2f%%\n",sumhaderr_pi_nsknd[0]*100, sumhaderr_pi_nsknd[1]*100, sumhaderr_pi_nsknd[2]*100, sumhaderr_pi_nsknd[3]*100);
    printf("Kaon Prod.\t %.2f%%\t\t %.2f%%\t\t %.2f%%\t\t %.2f%%\n",sumhaderr_k_nsknd[0]*100, sumhaderr_k_nsknd[1]*100, sumhaderr_k_nsknd[2]*100, sumhaderr_k_nsknd[3]*100);
    for(int i=0;i<nsys;i++){
      printf("%s\t %.2f%%\t\t %.2f%%\t\t %.2f%%\t\t %.2f%%\n",Sys[i].GetName(),TMath::Abs(Sys[i].sumsyserr_nsknd[0]*100), TMath::Abs(Sys[i].sumsyserr_nsknd[1]*100), TMath::Abs(Sys[i].sumsyserr_nsknd[2]*100), TMath::Abs(Sys[i].sumsyserr_nsknd[3]*100));
    }
  };
  
  void ShowDetail(){
    std::cout << "\n-------- detailed outputs -------------------- \n" << std::endl;
    std::cout << "-------- number of events (input) ------------" << std::endl;
    std::cout << "\t\ttotal\tnumu\tnumub\tnue\t " << std::endl;
    printf("ND: \t\t%.2f \n",int_nnd[3]);
    printf("SK 1Rmu: \t%.2f\t%.2f\t%.2f\t%.2f\n",int_nsk[0][3],int_nsk[0][0],int_nsk[0][1],int_nsk[0][2]);
    printf("SK NueSig: \t%.2f\n",int_nsk[1][3]);
    printf("SK NueBG: \t%.2f\t%.2f\t%.2f\t%.2f\n",int_nsk[2][3],int_nsk[2][0],int_nsk[2][1],int_nsk[2][2]);
    
    ShowHad();
    ShowSys();
  }

  int GetEnuBin(int bin_in){
    if(bin_in < 1) return 0;
    else if(bin_in <=16 ) return (bin_in+1)/2;
    else if(bin_in <=24 ) return (bin_in+3)/4 + 4;
    else if(bin_in <=30 ) return 11;
    else if(bin_in <=80 ) return (bin_in+9)/10 + 8;
    else if(bin_in <=100) return 17;
    else if(bin_in <=140) return 18;
    else if(bin_in <=200) return 19;
    else return 20;
  }
};
