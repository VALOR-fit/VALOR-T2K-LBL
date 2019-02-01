/*
  Take flux from a fluka tree, and format it in the usually 50 MeV binning and normalisation /1E21 POT
 */

#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

TTree * GetTree(const char * fname, const char * tname = "h2000")
{
  TFile * f = new TFile(fname);
  TTree * t = 0;
  const char * treename = tname;
  f->GetObject(treename,t);
  if(t)
    return t;
  cerr << "Could not get tree " << treename
       << " from file " << fname << endl;
  exit(-1);
}

TH1D * GetHisto(TTree * t, const char * cut, const char * hname)
{
  TH1D * hnew = new TH1D(hname, hname, 600, 0, 30);
  cout << hname << endl;
  cout << t->Draw(TString::Format("Enu>>%s", hname),cut) << " events drawn by TTree::Draw()" << endl;
  cout << "histogram has " << hnew->GetEntries() << " entries" << endl;
  return hnew;
}

void get_flux_from_fluka(const char * infhc="fluka_flux_b_2_hktoch_pos_all.root", const char * inrhc="fluka_flux_b_2_hktoch_neg_all.root", const char * fnameout = "hk_flux_tochi.root")
{
  TTree * tfhc = GetTree(infhc);
  TTree * trhc = GetTree(inrhc);

  TTree * tfhc_version = GetTree(infhc, "h1000");
  TTree * trhc_version = GetTree(inrhc, "h1000");
  int nfiles[2];
  nfiles[0] = tfhc_version->GetEntries();
  nfiles[1] = trhc_version->GetEntries();

  TFile fout(fnameout, "RECREATE");
  for(int ib = 0; ib < 2; ib++) {
    TTree * t;
    string hname_beam;
    switch(ib){
    case 0:
      t = tfhc;
      hname_beam = "nu_mode";
      break;
    case 1:
      t = trhc;
      hname_beam = "nubar_mode";
      break;
    }//switch(ib)
    for(int im = 0; im < 4; im++) {
      string hname_flav;
      string cut;
      switch(im) {
      case 0:
	cut = TString::Format("norm/%d*(mode >= 11 && mode < 16)", nfiles[ib]);
	hname_flav = "numu";
	break;
      case 1:
	cut = TString::Format("norm/%d*(mode >= 21 && mode < 26)", nfiles[ib]);
	hname_flav = "numubar";
	break;
      case 2:
	cut = TString::Format("norm/%d*(mode >= 31 && mode < 36)", nfiles[ib]);
	hname_flav = "nue";
	break;
      case 3:
	cut = TString::Format("norm/%d*(mode >= 41 && mode < 46)", nfiles[ib]);
	hname_flav = "nuebar";
	break;
      }//switch(im)
      TString hname = TString::Format("flux_%s_%s", hname_beam.c_str(), hname_flav.c_str()).Data();
      cout << hname << " made of " << nfiles[ib] << " fluka files" << endl;
      TH1D * h = GetHisto(t, cut.c_str(), hname.Data());
      h->Write();
    }//im neutrino flavour
  }//ib beam mode
}
