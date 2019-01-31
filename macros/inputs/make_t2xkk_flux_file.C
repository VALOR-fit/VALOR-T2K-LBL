#include "TFile.h"
#include "TH1D.h"
#include <iostream>
#include "TMath.h"

/*
  Takes a flux file from HK that just has the absolute flux (for both FHC & RHC)
  Takes RHC/FHC flux files from T2K, that has the absolute flux the SKMC was made with
  Spits out a VALOR formatted flux file that has both the absolute weighted flux & the flux weights
 */

void make_t2xkk_flux_file(
			  bool verbose_binning = false,
			  bool verbose_weights = false,
			  const char * outfilename = "inputs/hk/flux/flux_weights_and_flux_2p5_295km.root",
			  const char * in_korea_filename = "inputs/hk/flux/flux_295km_2p5degree.root",
			  const char * in_fhc_filename = "$T2KOSCANA/inputs/official/2018/flux/tuned13av3/run1-9a_v1/sk_tuned13av3_13anom_run1-9a_numode_fine.root",
			  const char * in_rhc_filename = "$T2KOSCANA/inputs/official/2018/flux/tuned13av3/run5c-9d_v1/sk_tuned13av3_13anom_run5c-9d_antinumode_fine.root"
			  )
{
  TFile * ffhc = new TFile(in_fhc_filename);
  TFile * frhc = new TFile(in_rhc_filename);
  TFile * fkorea = new TFile(in_korea_filename);
  fkorea->ls();

  const int nspecies = 4;
  const string species[nspecies] = { "numu", "numub", "nue", "nueb" };
  const string species_korea[nspecies] = { "numu", "numubar", "nue", "nuebar" };
  const int nmodes = 2;
  const string modes[nmodes] = { "numode", "anumode" };
  const string modes_korea[nmodes] = { "nu_mode", "nubar_mode" };

  TFile fout(outfilename, "RECREATE");

  for(int im = 0; im < nmodes; im++) {
    cout << modes[im] << endl;
    TFile * f = (im == 1 ? frhc : ffhc);
    //f->ls();
    for(int is = 0; is < nspecies; is++) {
      cout << species[is] << endl;
      
      //these are finely binned histograms
      TH1D * tune13av11 = 0, * nom = 0, * ratio = 0;
      f->GetObject(TString::Format("enu_sk_tuned13a_%s", species[is].c_str()), tune13av11);
      f->GetObject(TString::Format("enu_sk_13a_nom_%s", species[is].c_str()), nom);
      f->GetObject(TString::Format("enu_sk_tuned13a_%s_ratio", species[is].c_str()), ratio);
      cout << tune13av11 << "\t" << nom << "\t" << ratio << endl;
      assert(tune13av11);
      assert(nom);
      assert(ratio);
      
      //this is (potentially) a coarsely binned histogram
      TH1D * tunekorea = 0;
      fkorea->GetObject(TString::Format("flux_%s_%s", modes_korea[im].c_str(), species_korea[is].c_str()), tunekorea);
      cout << tunekorea << endl;
      assert(tunekorea);

      if(verbose_binning) {
	for(int ib = 0; ib < tune13av11->GetNbinsX(); ib++) {
	  double centre = nom->GetXaxis()->GetBinCenter(ib);
	  cout << "Listing upper bin edges. Fine bin " << ib
	       << " nom " << nom->GetXaxis()->GetBinUpEdge(nom->GetXaxis()->FindBin(centre))
	       << " ratio " << ratio->GetXaxis()->GetBinUpEdge(ratio->GetXaxis()->FindBin(centre))
	       << " tune13av11 " << tune13av11->GetXaxis()->GetBinUpEdge(tune13av11->GetXaxis()->FindBin(centre))
	       << " tunekorea " << tunekorea->GetXaxis()->GetBinUpEdge(tunekorea->GetXaxis()->FindBin(centre))
	       << endl;
	}
      }
      /*
      assert(tune13av11->GetNbinsX() == nom->GetNbinsX());
      assert(tune13av11->GetNbinsX() == ratio->GetNbinsX());
      assert(tune13av11->GetNbinsX() == tunekorea->GetNbinsX());
      */

      //create new histograms using the FINE binning
      TH1D * newnom = new TH1D(*tune13av11);
      TH1D * newtune = new TH1D(*tune13av11);
      TH1D * newtuneratio = new TH1D(*tune13av11);

      /*    
      if(tunekorea->GetNbinsX() != tune13av11->GetNbinsX())
	cout << "NBins korea " << tunekorea->GetNbinsX() << " SK " << tune13av11->GetNbinsX() << endl;
      for(int ib = 0; ib <= tunekorea->GetNbinsX(); ib++) {
	double centerk = tunekorea->GetBinCenter(ib);
	double centers = tune13av11->GetBinCenter(ib);
	if(TMath::Abs(centerk - centers) > 1E-6) {
	  cout << "Bin center " << ib << " korea " << centerk << " SK " << centers << endl;
	}
      }
      */

      for(int ib = 0; ib <= newtune->GetNbinsX(); ib++) {
	double center = newtune->GetBinCenter(ib);
	double c_tune13av11   = tune13av11->GetBinContent(tune13av11->FindBin(center));
	double c_nom    = nom->GetBinContent(nom->FindBin(center));
	double c_ratio  = ratio->GetBinContent(ratio->FindBin(center));
	double c_nomkorea = tunekorea->GetBinContent(tunekorea->FindBin(center));
	double c_tunekorea = c_nomkorea * c_ratio;

	//tune = ratio * nom
	double c_newratio = c_tunekorea / c_nom;
	
	if(verbose_weights || (TMath::Abs(c_tune13av11 - c_nom * c_ratio) > 1E-6)) {
	  cout << "Bin centre" << center << endl
	       << "\tT2K: Tune = nom * ratio = nom x ratio \t" << c_tune13av11 << "\t=\t"
	       << c_nom * c_ratio << "\t=\t"
	       << c_nom << "\tx\t"
	       << c_ratio << endl
	       << "\tKOREA: Tune, Ratio \t" << c_tunekorea << "\t" << c_newratio << endl;
	}

	if(TMath::IsNaN(c_newratio)) c_newratio = 0;

	newnom->SetBinContent(ib, c_nomkorea);
	newtune->SetBinContent(ib, c_tunekorea);
	newtuneratio->SetBinContent(ib, c_newratio);
      }//ib

      fout.cd();
      newnom->Write(TString::Format("enu_sk_nomkorea_%s_%s", modes[im].c_str(), species[is].c_str()));
      newtune->Write(TString::Format("enu_sk_tunedkorea_%s_%s", modes[im].c_str(), species[is].c_str()));
      newtuneratio->Write(TString::Format("%s_ratio", tunekorea->GetName()));
    }//is
  }//im

  fout.ls();
}
