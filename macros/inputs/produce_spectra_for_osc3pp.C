#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH3D.h"
#include "TH2D.h"
#include "Binning.h"

using namespace std;
using namespace valor;
using namespace valor::conv;
using namespace valor::util;

#if 0
const int nbins_etrue = 98;
const double edges_etrue[nbins_etrue + 1] = {
  0.0, 0.05, 0.1, 0.15,
  0.2, 0.225, 0.25, 0.275,
  0.3, 0.325, 0.35, 0.375,
  0.4, 0.425, 0.45, 0.475,
  0.5, 0.525, 0.55, 0.575,
  0.6, 0.625, 0.65, 0.675,
  0.7, 0.725, 0.75, 0.775,
  0.8, 0.825, 0.85, 0.875,
  0.9, 0.925, 0.95, 0.975,
  1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 
  1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75,
  1.8, 1.85, 1.9, 1.95,
  2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 
  2.4, 2.45, 2.5,
  2.6, 2.7, 2.8, 2.9,
  3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
  4.0, 4.2, 4.4, 4.6, 4.8,
  5.0, 5.2, 5.4, 5.6, 5.8,
  6.0, 6.5, 7.0, 7.5, 
  8, 9, 10, 30
};
const int nbins_ereco = 57;
const double edges_ereco[nbins_ereco + 1] = {
  0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
  1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
  2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
  3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
  4, 4.2, 4.4, 4.6, 4.8,
  5, 5.2, 5.4, 5.6, 5.8,
  6, 6.5, 7, 7.5,
  8, 9, 10, 30
};
#else
const int nbins_etrue = kNEtrueBins1Rmu;
const double * edges_etrue = kEtrueBinEdges1Rmu;
const int nbins_ereco = kNErecoBins1Rmu;
const double * edges_ereco = kErecoBinEdges1Rmu;
#endif

double produce_a_spectrum(const char * intag,
			  TFile * outfile,
			  vector<int> modes,
			  const char * selection,
			  const char * valorselection,
			  int iselection,
			  const char * ints)
{
  double tot = 0;
  TFile f(Form("%s%s.root", intag, valorselection));
  TH2D * hnew = new TH2D(Form("enu_erec_%s_%s", selection, ints),
		      Form("enu_erec_%s_%s;E_{true};E_{reco}", selection,ints),
		      nbins_etrue, edges_etrue, nbins_ereco, edges_ereco);
  for(size_t im = 0; im < modes.size(); im++) {
    TH3D * h = 0;
    f.GetObject(Form("nsk_2d_%d_%d", iselection, modes[im]), h);
    if(!h) {
      f.ls();
      cout << "Looking for: " << Form("nsk_2d_%d_%d", iselection, modes[im]) << endl;
      cout << "For " << hnew->GetTitle() << " in " << outfile->GetName() << endl;
    }
    cout << "Adding events from " << h->GetTitle() << " in " << selection << endl
	 << "\tto " << hnew->GetTitle() << " in " << outfile->GetName() << endl;
    for(int ix = 1; ix <= h->GetNbinsX(); ix++) {
      for(int iy = 1; iy <= h->GetNbinsY(); iy++) {
	for(int iz = 1; iz <= h->GetNbinsZ(); iz++) {
	  int global_bin = h->GetBin(ix, iy, iz);
	  double ereco = h->GetYaxis()->GetBinCenter(iy);
	  double etrue = h->GetZaxis()->GetBinCenter(iz);
	  double value = h->GetBinContent(global_bin);
	  hnew->Fill(etrue, ereco, value);
	  tot += value;
	}//iz
      }//iy
    }//ix
  }//im
  outfile->cd();
  hnew->Write();
  return tot;
}

void produce_spectra_for_osc3pp(const char * intag = "Template2D_mcnorm_fluxtuned_banfftuned_syst_multi")
{
  const int nsamples = 2;
  string samples[nsamples] = {"1Rmu", "1Re"};
  const int nints = 3;
  string ints[nints] = {"CCQE", "CCnQE", "NC"};
  const int nfluxes = 12;
  string fluxes[nfluxes] = {"numu_x_numu", "numubar_x_numubar", "nue_x_nue", "numu_x_nue", "nuebar_x_nuebar", "numubar_x_nuebar",
			    "numu_crs_numu_anu", "numubar_crs_numubar_anu", "nue_crs_nue_anu", "numu_crs_nue_anu", "nuebar_crs_nuebar_anu", "numubar_crs_nuebar_anu"};
  double tot[4] = {0,0,0,0};
  for(int ifl = 0; ifl < nfluxes; ifl++) {
    TFile * outfile = new TFile(Form("enu_erec_%s.root", fluxes[ifl].c_str()), "RECREATE");
    for(int is = 0; is < nsamples; is++) {
      int isample = (ifl >= 6) ? is + 2 : is;
      const char * valorsample = (ifl >= 6) ? Form("RHC%s", samples[is].c_str()) : samples[is].c_str();
      for(int ii = 0; ii < nints; ii++) {
	int istart = -99;
	int nmodes = -99;
	//NC
	if(ii == 2) {
	  if(ifl % 6 == 3 || ifl % 6 == 5) {
	    continue; //skip NC for oscillated modes
	  }
	  istart = 5 + ((ifl % 6) * 10); //5, 15, 25
	  if(ifl % 6 == 4)
	    istart = 40; //nuebar
	  nmodes = 5;
	}//i == 2 // NC
	//CCQE
	else if(ii == 0) {
	  istart = ((ifl % 6) * 10); //0, 10, 20, 30
	  if(ifl % 6 == 4)
	    istart = 35; //nuebar
	  else if(ifl % 6 == 5)
	    istart = 45; //osc nuebar
	  nmodes = 1;
	}//i == 0 // CCQE
	//CCnQE
	else if(ii == 1) {
	  istart = 1 + ((ifl % 6) *10); //1, 11, 21, 31
	  if(ifl % 6 == 4)
	    istart = 36; //nuebar
	  else if(ifl % 6 == 5)
	    istart = 46; //osc nuebar
	  nmodes = 4;
	}//i==1 // CCnQE
	vector<int> modes;
	for(int imode = 0; imode < nmodes; imode++) {
	  modes.push_back(istart + imode);
	}//imode
	double thistot = produce_a_spectrum(intag, outfile, modes, samples[is].c_str(), valorsample, isample, ints[ii].c_str());
	//if(true)
	if(!(ifl % 6 == 3 || ifl % 6 == 5)) //don't add oscillated nue(bar) for no oscillations
	  tot[isample] += thistot;
      }//ii
    }//is
    outfile->Close();
  }//ifl
  for(int isample = 0; isample < 4; isample++) {
    cout << "Total for " << isample << (isample < 2 ? " FHC " : " RHC ") << samples[isample % 2] << " \t = " << tot[isample] << endl;
  }//isample
}

