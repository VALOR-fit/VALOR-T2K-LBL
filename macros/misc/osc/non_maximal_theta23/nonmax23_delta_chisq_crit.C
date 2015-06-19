/*
This macro is run before runnning NonMax23Discovery.C. 

Run 5000 experiments at true points in oscillation parameter space for sinsq(2*theta23) = 1.0 and each value 
of |dm32sq| in steps of 1E-4 eV2 using feldman_cousins_fitting_job.C. Include stat, syst and non23 variations.
This macro runs over the output files from those jobs, and writes the critical delta chisq values for each 
point at 68%, 90%, 95% and 99% CL to TH2Ds in a root file. These critical values are read by NonMax23Discovery.C.

Please see NonMax23Discovery.C for more details of how to run these macros.

Nick Grant
*/

#include <vector>
#include <iostream>

#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TChain.h>

using namespace std;

void nonmax23_delta_chisq_crit(const char *input_directory = "production/NH_3E20POT")
{
  const int nConfidenceLevel = 4;
 
  int confidenceLevel[nConfidenceLevel] = {68, 90, 95, 99};
  int percentile[nConfidenceLevel];

  const int mixing23_nbins = 1;

  double dm23sq_min = 2.2E-3;
  double dm23sq_max = 2.8E-3;
  double dm23sq_step = 1E-4;

  int nbins = (int)(1 + 1.001*(dm23sq_max - dm23sq_min) / dm23sq_step);
  const int dm23sq_nbins = nbins;

  TH2D * hDeltaChisqSurface[nConfidenceLevel];
  for(int i = 0; i < nConfidenceLevel; i++)
    {
    const char * histname = Form("hDeltaChisqSurface%02d", confidenceLevel[i]);
    hDeltaChisqSurface[i] = new TH2D(histname, "", mixing23_nbins, 0.995, 1.005,
	                  dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
    } 

  vector<double> deltaChisq;
  vector<double>::iterator deltaChisqIter;

  double bestfit_mixing23;
  double bestfit_dm23sq;
  double delta_chisq;
  double chisq_for_input_params;
  double bestfit_chisq;

  int mixing23_tag, dm23sq_tag;

  for(int mixing23_bin=0; mixing23_bin<mixing23_nbins; mixing23_bin++)
    {
    mixing23_tag = 100;
    for(int dm23sq_bin=0; dm23sq_bin<dm23sq_nbins; dm23sq_bin++)
      {
      dm23sq_tag = (int)((dm23sq_min + dm23sq_bin * dm23sq_step) * 1.001E5);

      TChain *chain = new TChain("fcnt");
      chain->Add(Form("$T2KOSCANA/%s/fc_%03d_%03d/*.root", input_directory, mixing23_tag, dm23sq_tag));

      chain->SetBranchAddress("bestfit_mixing23",&bestfit_mixing23);
      chain->SetBranchAddress("bestfit_dm23sq",&bestfit_dm23sq);
      chain->SetBranchAddress("bestfit_chisq",&bestfit_chisq);
      chain->SetBranchAddress("chisq_for_input_params",&chisq_for_input_params);
      chain->SetBranchAddress("delta_chisq",&delta_chisq);

      for(Long64_t j = 0; j < chain->GetEntries(); j++)
        {
        chain->GetEntry(j);

        deltaChisq.push_back(delta_chisq);
        }

      sort(deltaChisq.begin(), deltaChisq.end());

    for(int i=0; i<nConfidenceLevel; i++) 
      {
      percentile[i] = (int)(0.01*confidenceLevel[i]*deltaChisq.size());
      deltaChisqIter = deltaChisq.begin() + percentile[i];
      hDeltaChisqSurface[i]->SetBinContent(mixing23_bin+1, dm23sq_bin+1, (*deltaChisqIter));
      }

    deltaChisq.clear();
    delete chain;

    }//end of for(int dm23sq_bin=0; dm23sq_bin<dm23sq_nbins; dm23sq_bin++)
  }//end of for(int mixing23_bin=0; mixing23_bin<mixing23_nbins; mixing23_bin++)  

  TFile *outfile = new TFile("DeltaChisqCrit.root", "RECREATE");
  outfile->cd();
  for(int i=0; i<nConfidenceLevel; i++)
    hDeltaChisqSurface[i]->Write();
  outfile->Close();

}//end of NonMax23DeltaChisqCrit()


