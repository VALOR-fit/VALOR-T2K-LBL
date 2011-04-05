//-----------------------------------------------------------------------------------------------------
// Show number of fits generated so far (for calculation of Feldman-Cousins \Delta\chi^2 surfaces)
// as function of \Deltam^2 and sin^2(2\theta).
// Plots are used to inform choices for further MC processing.
//-----------------------------------------------------------------------------------------------------

#include <vector>

void test_feldman_cousins_nfits(const char * directory)
{
  gROOT->Reset();
  gStyle->SetOptStat(11);

  //
  // Define grid
  //

  // dm23sq, sinsq2theta23 range and step
  double sinsq2theta23_min   = 0.0;
  double sinsq2theta23_max   = 1.0;
  double sinsq2theta23_step  = 0.01;
  double dm23sq_min          = 1E-3;
  double dm23sq_max          = 6E-3;
  double dm23sq_step         = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
  int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );

  cout << "sinsq2theta23_nbins = " << sinsq2theta23_nbins << ", dm23sq_nbins = " << dm23sq_nbins<<endl;

  TH2D* hNFitsVsOscPar = new TH2D("hNFitsVsOscPar", "",
          sinsq2theta23_nbins, 
          sinsq2theta23_min - 0.5*sinsq2theta23_step,
          sinsq2theta23_max + 0.5*sinsq2theta23_step,
          dm23sq_nbins,
          dm23sq_min - 0.5*dm23sq_step,
          dm23sq_max + 0.5*dm23sq_step);

  //          
  // Load results of toy-experiment fits
  //          

  TChain chain("fcnt");
  chain.Add(Form("%s/*root",directory));

  chain.Draw("input_dm23sq:input_sinsq2theta23>>hNFitsVsOscPar","","goff");

  TCanvas * c = new TCanvas("c", "", 20, 20, 600, 600);
  hNFitsVsOscPar->Draw("colz");  
}
