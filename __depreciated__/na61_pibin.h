//search NA61 pion binID from p and theta //

#ifndef NA61PIBIN_INCLUDED
#define NA61PIBIN_INCLUDED

#define THETA_NBIN 10
const int nbins_pip[THETA_NBIN]=
  {27, 45, 52, 47, 40, 33, 29, 15, 10, 10};
const int nbins_pim[THETA_NBIN]=
  {31, 48, 51, 48, 42, 37, 34, 27, 22, 20};

const double theta_bin_le[THETA_NBIN]= 
  {0, 20, 40, 60, 100, 140, 180, 240, 300, 360};
const double theta_bin[THETA_NBIN]= 
  {20, 40, 60, 100, 140, 180, 240, 300, 360, 420};

class NA61PiBin{
 private:
  int start_bin_pip[THETA_NBIN+1];
  int start_bin_pim[THETA_NBIN+1];
  TFile *file;
  TH1F * h_pip[THETA_NBIN];
  TH1F * h_pim[THETA_NBIN];
  TH1F * h_pip_e[THETA_NBIN];
  TH1F * h_pim_e[THETA_NBIN];
  
 public:
  
  NA61PiBin(const char * fname = "normerr_input/input/na61_final.root"){
    std::cout << "init : NA61 bin & error infos" << std::endl;

    int sbp=0, sbm=0;
    file = TFile::Open(fname);
    TString hname;
    
    for(int i=0;i<THETA_NBIN;i++){
      hname.Form("h_piPlus_%d_%d",int(theta_bin_le[i]),int(theta_bin[i]));
      h_pip[i] = (TH1F*)file->Get(hname);
      hname.Form("h_piMinus_%d_%d",int(theta_bin_le[i]),int(theta_bin[i]));
      h_pim[i] = (TH1F*)file->Get(hname);
      hname.Form("h_piPlus_%d_%d_sys_tot",int(theta_bin_le[i]),int(theta_bin[i]));
      h_pip_e[i] = (TH1F*)file->Get(hname);
      hname.Form("h_piMinus_%d_%d_sys_tot",int(theta_bin_le[i]),int(theta_bin[i]));
      h_pim_e[i] = (TH1F*)file->Get(hname);
      start_bin_pip[i] = sbp;
      start_bin_pim[i] = sbm;
      sbp += nbins_pip[i];
      sbm += nbins_pim[i];

      h_pip_e[i]->Divide(h_pip[i]);
      h_pim_e[i]->Divide(h_pim[i]);
    }
    start_bin_pip[THETA_NBIN] = sbp;
    start_bin_pim[THETA_NBIN] = sbm;
    //    check_startbins();
  };

  void check_startbins(){
    for(int i=0;i<THETA_NBIN+1;i++){
      std::cout <<start_bin_pip[i]<<" "<< start_bin_pim[i]<<std::endl;
    }
  };
  
  int get_histID(double mrad){
    if(mrad<0) return -1;
    int id = THETA_NBIN;
    for(int i=0;i<THETA_NBIN;i++){
      if(mrad<theta_bin[i]){
	id = i;
	break;
      }
    }
    //    std::cout << id <<std::endl;
    return id;
  };

  int get_binID(int q, int histID, double mom){
    if(q<0 || q>1 || histID<0) return -1;
    else if(!q){//pip
      if(histID>=THETA_NBIN) return start_bin_pip[THETA_NBIN];
      return h_pip[histID]->FindBin(mom) + start_bin_pip[histID];
    }
    else{//pim
      if(histID>=THETA_NBIN) return start_bin_pim[THETA_NBIN];
      return h_pim[histID]->FindBin(mom) + start_bin_pim[histID];
    }
  };
  
  int get_binID_pip(double mrad, double mom){
    return get_binID(0, get_histID(mrad), mom);
  };
  int get_binID_pim(double mrad, double mom){
    return get_binID(1, get_histID(mrad), mom);
  };

  double get_binErr_pip(double mrad, double mom){
    // 0 or minus value means out of range //
    // will be converted to 50% in NormErrInput.cc //
    int hid = get_histID(mrad);
    if(hid<0 || hid >= THETA_NBIN) return -1; // angle out of range

    TH1F *h = h_pip_e[get_histID(mrad)];
    int lbid = h->FindBin(mom);
    if(lbid < 1 || lbid > h->GetNbinsX()) return -2; // momentum out of range
    else return h->GetBinContent(h->FindBin(mom));
  };
  
  double get_binErr_pim(double mrad, double mom){
    // 0 or minus value means out of range //
    // will be converted to 50% in NormErrInput.cc //
    int hid = get_histID(mrad);
    if(hid<0 || hid >= THETA_NBIN) return -1; // angle out of range

    TH1F *h = h_pim_e[get_histID(mrad)];
    int lbid = h->FindBin(mom);
    if(lbid < 1 || lbid > h->GetNbinsX()) return -2; // momentum out of range
    else return h->GetBinContent(h->FindBin(mom));
  };

};

#endif
