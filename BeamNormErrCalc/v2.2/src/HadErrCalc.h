#define DEF_NUMU_INP "delta_flux_nd5_numu.root"
#define DEF_NUMUB_INP "delta_flux_nd5_numub.root"
#define DEF_NUE_INP "delta_flux_nd5_nue.root"

#define DEF_NUMU_INP_ING "delta_flux_ing13_numu.root"
#define DEF_NUMUB_INP_ING "delta_flux_ing13_numub.root"
#define DEF_NUE_INP_ING "delta_flux_ing13_nue.root"

#define MAX_NSYS 10 // maximum # of sys. error terms
#define PBIN 816   // Pion phase space # of bins
#define KBIN 700   // Kaon phase space # of bins

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

class BinErr{
 public:
  double data[3][20]; // [flavor],[enu]
};

class SysErr{
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
    std::cout << "[befor cancellation]" << std::endl;
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
  SysErr Sys[MAX_NSYS]; // sys. errors other than primary pi/K
  
  // error contribution inputs
  BinErr pi_econt_sk[2][PBIN]; // [+/-],[pion]
  BinErr k_econt_sk[3][KBIN];  // [+/-/0],[kaon]
  BinErr pi_econt_nd[2][PBIN]; // [+/-],[pion]
  BinErr k_econt_nd[3][KBIN];  // [+/-/0],[kaon]
  
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
  double err_pi_nsk[PBIN][3][3]; // [pion],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_pi_nnd[PBIN]; // [pion]
  double err_k_nsk[KBIN][3][3]; // [kaon],[1rmu, nue_sig, nue_bkg],[numu,numub,nue]
  double err_k_nnd[KBIN]; // [kaon]
  
  // Nsk/Nnd error outputs
  double err_pi_nsknd[PBIN][3][3]; // [pion],[1rmu, nue_sig, nue_bkg][numu,numub,nue]
  double err_k_nsknd[KBIN][3][3]; // [kaon],[1rmu, nue_sig, nue_bkg][numu,numub,nue]
  
  // error contribution histograms //
  TH2D * hphase[2]; // {pi,K}
  TH2D * hcontnd[2]; // {pi,K}
  TH2D * hcontsk[2][4]; // {pi,K}{1rmu, nue_sig, nue_bg, nue_all}
  TH2D * hcontsknd[2][4];  // {pi,K}{1rmu, nue_sig, nue_bg, nue_all}
  
  double err_pi_nsk_all[PBIN][4];
  double err_pi_nsknd_all[PBIN][4];
  double err_k_nsk_all[KBIN][4];
  double err_k_nsknd_all[KBIN][4];
  
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


  TH2D * GetContHist(int det, int pik, int mode = 0){
    TString s_pik[2] = {"pi","k"};
    TString s_mode[4] = {"1Rmu", "NueSig", "NueBG", "NueAll"};
    TString s_x[2] = {"P_{#pi} [GeV/c]", "X_{F}"};
    TString s_y[2] = {"#theta [mrad]", "P_{T} [GeV/c]"};
    int nbins[2] = {816, 700};
    TH2D *h;
    if(det == 1){//ND
      h = hcontnd[pik];
      if(!h){
	TString name = "h_nd_" + s_pik[pik]; 
	h = (TH2D*)hphase[pik]->Clone(name);
	h->SetTitle(name);
      }
      if(!pik)
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_pi_nnd[i]);
      else
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_k_nnd[i]);
    }
    
    if(det == 0){//SK
      h = hcontsk[pik][mode];
      if(!h){
	TString name = "h_sk_" + s_mode[mode]+ "_" + s_pik[pik]; 
	h = (TH2D*)hphase[pik]->Clone(name);
	h->SetTitle(name);
      }
      if(!pik)
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_pi_nsk_all[i][mode]);
      else
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_k_nsk_all[i][mode]);
    }
    
    if(det == 2){//SK/ND
      h = hcontsknd[pik][mode];
      if(!h){
	TString name = "h_sk_" + s_mode[mode]+ "_" + s_pik[pik]+"/nd"; 
	h = (TH2D*)hphase[pik]->Clone(name);
	h->SetTitle(name);
      }
      if(!pik)
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_pi_nsk_all[i][mode] - err_pi_nnd[i]);
      else
	for(int i=0;i<nbins[pik];i++) h->SetBinContent(i,err_k_nsk_all[i][mode] - err_k_nnd[i]);
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
