/*
This macro makes the plots of non-maximal theta23 discovery potential that were shown at the T2K collaboration meeting on 
29th Sept 2012. These plots are on slide 21 in

http://www.t2k.org/meet/collab/201209/talks/satam/Joint%20Oscillation%20Analysis%20Report%20-%20VaLOR/view

There is more information on the method on the backup slides 45 and 46 in that presentation.

Run 5000 experiments at true points in oscillation parameter space for sinsq(2*theta23) = 1.0 and each value of |dm32sq| 
in steps of 1E-4 eV2 using feldman_cousins_fitting_job.C. Include stat, syst and non23 variations. Run NonMax23DeltaChisqCrit.C 
over the output of those jobs - this writes the critical values of delta chisq for 68%, 90%, 95% and 99% CL for each point 
to TH2Ds in a root file.  

Now run 5000 fake experiments for each true oscillation parameter point using test_n_mock_data_fits.C. Use steps of 0.01 
in sinsq(2*theta23) and 1E-4 in |dm32sq|, and again include stat, syst and non23 variations. Run this macro over the output 
of those jobs - it finds for each point the fraction of experiments for which the entire confidence region is contained within 
sinsq(2*theta23) < 1.0. This is done experiment by experiment by finding the delta chisq for that experiment between the 
best-fit point and the threshold of 1.0. It is necessary to add these lines to test_n_mock_data_fits.C:                                                                                  

fparm->Update(kPhypMixing23, threshold);
double threshold_chisq = fit_inputs->Chisquare();                                                                                                                                  

where threshold = 1.0 and threshold_chisq is written to the output ntuple. Then threshold_chisq - bestfit_chisq is compared 
with the delta chisq critical value for the value of |dm32sq| from the feldman_cousins_fitting_job.C jobs that is nearest 
to the best-fit value of |dm32sq| for that experiment. 

This is repeated this for each of 68%, 90%, 95% and 99% CL.

Percentage contours for this fraction are found by linear interpolation between the oscillation parameter points (I tried fitting 
a TSpline3 but found that the fit was very poor and in many cases gave an incorrect interpolation). Contours for 10%, 30%, 50%, 
70% and 90% are drawn, but a contour is drawn only if the interpolation can be made over the entire range in |dm32sq|. 
Some of these contours might not be drawn if that probability does not fall in the |dm32sq| range of the test_n_mock_data_fits.C 
jobs.

This macro assumes that the |dm32sq| range of the feldman_cousins_fitting_job.C jobs extends one grid point both above and below 
the range in |dm32sq| of the test_n_mock_data_fits.C jobs. This is to allow for best-fit points in the test_n_mock_data_fits.C 
datasets that are in those regions. 

If plot_bestfit_prob is true, it also plots percentage contours of best-fit points that are < 1.0. 

Nick Grant
*/

#include <vector>
#include <iostream>

#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TChain.h>

using namespace std;

void FormatCanvas(TCanvas *canvas, double mixing23_min, double mixing23_max, double dm23sq_min, double dm23sq_max);

void FormatGraph(TGraph *graph, int colour, int style, int width);

void nonmax23_discovery(
const char *input_directory = "production/NH_3E20POT",
bool plot_bestfit_prob = false,
bool plot_68CL_prob = true,
bool plot_90CL_prob = true,
bool plot_95CL_prob = true,
bool plot_99CL_prob = true
)
{
  const int nCL = 4;

  const int nProbValues = 5;
  double probValues[nProbValues] = {0.1, 0.3, 0.5, 0.7, 0.9};

  double threshold = 0.99999;

  double mixing23_min = 0.85;
  double mixing23_max = 0.95;
  double mixing23_step = 0.01;

  int nbins = (int)(1 + 1.001*(mixing23_max - mixing23_min) / mixing23_step);
  const int mixing23_nbins = nbins;

  double dm23sq_min = 2.3E-3;
  double dm23sq_max = 2.7E-3;
  double dm23sq_step = 1E-4;

  nbins = (int)(1 + 1.001*(dm23sq_max - dm23sq_min) / dm23sq_step);
  const int dm23sq_nbins = nbins;

  double deltaChisqCrit[nCL][dm23sq_nbins+2];

  TFile *DChisq = TFile::Open("DeltaChisqCrit.root");
  DChisq->cd();

  for(int dm23sq_bin=0; dm23sq_bin<dm23sq_nbins+2; dm23sq_bin++)
    {
    deltaChisqCrit[0][dm23sq_bin] = hDeltaChisqSurface68->GetBinContent(1, dm23sq_bin+1); 
    deltaChisqCrit[1][dm23sq_bin] = hDeltaChisqSurface90->GetBinContent(1, dm23sq_bin+1);
    deltaChisqCrit[2][dm23sq_bin] = hDeltaChisqSurface95->GetBinContent(1, dm23sq_bin+1);
    deltaChisqCrit[3][dm23sq_bin] = hDeltaChisqSurface99->GetBinContent(1, dm23sq_bin+1);
    }

  double bestfit_mixing23, bestfit_dm23sq, bestfit_chisq, threshold_chisq; 

  double bestfitFraction;
  double CLFraction[nCL];

  bool bestfitFracFound;
  bool CLFracFound[nCL];

  vector<double> bestfitMix23;
  vector<double>::iterator bestfitMix23Iter;

  vector<double> deltaChisq[nCL];
  vector<double>::iterator deltaChisqIter[nCL];

  double mix23Points[mixing23_nbins];
  for(int p=0; p<mixing23_nbins; p++)
    mix23Points[p] = mixing23_min + mixing23_step * p;

  double dm23sqPoints[dm23sq_nbins];
  for(int p=0; p<dm23sq_nbins; p++)
    dm23sqPoints[p] = dm23sq_min + dm23sq_step * p;

  double CLFracPoints[nCL][mixing23_nbins];
  double bestfitFracPoints[mixing23_nbins];

  double interpolatedBestfitProb[nProbValues][dm23sq_nbins];
  int interpolatedBestfitPointsFound[nProbValues] = {0};
  double interpolatedCLProb[nCL][nProbValues][dm23sq_nbins];
  int interpolatedCLPointsFound[nCL][nProbValues] = {0};

  int mixing23_tag, dm23sq_tag;

  int bestfit_dm23sq_bin;

  for(int dm23sq_bin=0; dm23sq_bin<dm23sq_nbins; dm23sq_bin++)
    {
    //dm23sq_tag might need to be changed to fit the names of the directories
    dm23sq_tag = (int)((dm23sq_min + dm23sq_bin * dm23sq_step) * 1.001E5);
      
    for(int mixing23_bin=0; mixing23_bin<mixing23_nbins; mixing23_bin++)
      {
      //mixing_tag might need to be changed to fit the names of the directories
      mixing23_tag = (int)((mixing23_min + mixing23_bin * mixing23_step) * 100.1);  

      TChain *chain = new TChain("nfitnt");
      chain->Add(Form("$T2KOSCANA/%s/toyexptfits_%03d_%03d/*.root", input_directory, mixing23_tag, dm23sq_tag));

      chain->SetBranchAddress("bestfit_mixing23",&bestfit_mixing23);  
      chain->SetBranchAddress("bestfit_dm23sq",&bestfit_dm23sq);
      chain->SetBranchAddress("bestfit_chisq",&bestfit_chisq);
      chain->SetBranchAddress("threshold_chisq",&threshold_chisq);

      for(Long64_t j = 0; j < chain->GetEntries(); j++)
        {
        chain->GetEntry(j);

        bestfitMix23.push_back(bestfit_mixing23);

        if(bestfit_dm23sq < (dm23sq_min - 0.5 * dm23sq_step))
	  bestfit_dm23sq_bin = 0;
        else if(bestfit_dm23sq > (dm23sq_max + 0.5 * dm23sq_step))
	  bestfit_dm23sq_bin = dm23sq_nbins + 1;
        else
	  bestfit_dm23sq_bin = (int)(1 + ((bestfit_dm23sq + 0.5 * dm23sq_step - dm23sq_min) / dm23sq_step));

        for(int k=0; k<nCL; k++)
          deltaChisq[k].push_back(threshold_chisq - bestfit_chisq - deltaChisqCrit[k][bestfit_dm23sq_bin]);
        }

      sort(bestfitMix23.begin(), bestfitMix23.end());
      for(int k=0; k<nCL; k++)
        sort(deltaChisq[k].begin(), deltaChisq[k].end());

      assert(bestfitMix23.size() == deltaChisq[0].size());

      const int nexpts = bestfitMix23.size();

      bestfitFracFound = false;
      for(int k=0; k<nCL; k++)
        {
        CLFracFound[k] = false;
        }

      //find fraction of experiments whose best-fit points are at sinsq(2*theta23) < 1
      //to be precise, it is the fraction whose best-fit points are at sinsq(2*theta23) <= threshold
      for(int i=0; i<nexpts; i++)
        {
        bestfitMix23Iter = bestfitMix23.begin() + i;
        if((*bestfitMix23Iter) > threshold && !bestfitFracFound)
	  {
          bestfitFraction = (float)i / nexpts;
          bestfitFracFound = true;
	  }
       
      //find fraction of experiments whose confidence regions are completely contained within sinsq(2*theta23) <= 1
      for(int k=0; k<nCL; k++)
        {
        deltaChisqIter[k] = deltaChisq[k].begin() + i;
        if((*deltaChisqIter[k]) > 0.0 && !CLFracFound[k])
	  {
	  CLFraction[k] = (float)(nexpts - i) / nexpts;
	  CLFracFound[k] = true;
	  }
        }
      }//end of for(int i=0; i<nexpts; i++)

      if(!bestfitFracFound)
        bestfitFraction = 1.0;

      for(int k=0; k<nCL; k++)
        {
        if(!CLFracFound[k])
          CLFraction[k] = 1.0;
        }
  
      bestfitFracPoints[mixing23_bin] = bestfitFraction;

      for(int k=0; k<nCL; k++)
        {
        CLFracPoints[k][mixing23_bin] = CLFraction[k];
        }

      bestfitMix23.clear();
      for(int k=0; k<nCL; k++)
        deltaChisq[k].clear();
 
      delete chain;

      }//end of for(int mixing23_bin=0; mixing23_bin<mixing23_nbins; mixing23_bin++)

    //interpolate probabilities of best-fit points between grid points in sinsq(2*theta23)
    //only loop to mixing23_nbins-1 so that interpolation can be made between the penultimate and last grid points
    for(int p=0; p<(mixing23_nbins-1); p++)
      {
      for(int prob=0; prob<nProbValues; prob++)
        {
        if(bestfitFracPoints[p] > probValues[prob] && bestfitFracPoints[p+1] < probValues[prob])
	  {
	  interpolatedBestfitProb[prob][dm23sq_bin] = mix23Points[p] + (mix23Points[p+1] - mix23Points[p]) * (bestfitFracPoints[p] - probValues[prob]) / (bestfitFracPoints[p] - bestfitFracPoints[p+1]);
	  interpolatedBestfitPointsFound[prob] += 1;
	  }
        }//end of for(int prob=0; prob<nProbValues; prob++)
      }//end of for(int p=0; p<(mixing23_nbins-1); p++)  

    //interpolate probabilities of confidence regions between grid points in sinsq(2*theta23)  
    for(int k=0; k<nCL; k++)
      {
      for(int p=0; p<(mixing23_nbins-1); p++)
        {
        for(int prob=0; prob<nProbValues; prob++)
          {
          if(CLFracPoints[k][p] > probValues[prob] && CLFracPoints[k][p+1] < probValues[prob])
            {
            interpolatedCLProb[k][prob][dm23sq_bin] = mix23Points[p] + (mix23Points[p+1] - mix23Points[p]) * (CLFracPoints[k][p] - probValues[prob]) / (CLFracPoints[k][p] - CLFracPoints[k][p+1]);
            interpolatedCLPointsFound[k][prob] += 1;
            }
          }
        }//end of for(int p=0; p<(mixing23_nbins-1); p++)
      }//end of for(int k=0; k<nCL; k++)
      

    }//end of for(int dm23sq_bin=0; dm23sq_bin<dm23sq_nbins; dm23sq_bin++)

  
  if(plot_bestfit_prob)
    {
      TCanvas *ProbBestfit = new TCanvas("ProbBestfit", "ProbBestfit", 800, 600);
      FormatCanvas(ProbBestfit, mixing23_min, mixing23_max, dm23sq_min, dm23sq_max);
      TLegend *lProbBestfit = new TLegend(0.2, 0.7, 0.45, 0.85);
      //only plot this graph if interpolation was possible at all values of dm23sq in the range of jobs
      if(interpolatedBestfitPointsFound[0] == dm23sq_nbins)
        {
          TGraph *Prob_Bestfit_10 = new TGraph(dm23sq_nbins, interpolatedBestfitProb[0], dm23sqPoints);
          FormatGraph(Prob_Bestfit_10, 6, 1, 2);
          Prob_Bestfit_10->Draw("L");
          lProbBestfit->AddEntry(Prob_Bestfit_10, "10%", "L");
	}
      if(interpolatedBestfitPointsFound[1] == dm23sq_nbins)
        {
          TGraph *Prob_Bestfit_30 = new TGraph(dm23sq_nbins, interpolatedBestfitProb[1], dm23sqPoints);
          FormatGraph(Prob_Bestfit_30, 1, 1, 2);
          Prob_Bestfit_30->Draw("L");
          lProbBestfit->AddEntry(Prob_Bestfit_30, "30%", "L");
        }
      if(interpolatedBestfitPointsFound[2] == dm23sq_nbins)
        {
          TGraph *Prob_Bestfit_50 = new TGraph(dm23sq_nbins, interpolatedBestfitProb[2], dm23sqPoints);
          FormatGraph(Prob_Bestfit_50, 8, 1, 2);
          Prob_Bestfit_50->Draw("L");
          lProbBestfit->AddEntry(Prob_Bestfit_50, "50%", "L");
	}
      if(interpolatedBestfitPointsFound[3] == dm23sq_nbins)
        {
          TGraph *Prob_Bestfit_70 = new TGraph(dm23sq_nbins, interpolatedBestfitProb[3], dm23sqPoints);
          FormatGraph(Prob_Bestfit_70, 4, 1, 2);
          Prob_Bestfit_70->Draw("L");
          lProbBestfit->AddEntry(Prob_Bestfit_70, "70%", "L");
        }
      if(interpolatedBestfitPointsFound[4] == dm23sq_nbins)
        {
	  TGraph *Prob_Bestfit_90 = new TGraph(dm23sq_nbins, interpolatedBestfitProb[4], dm23sqPoints);
	  FormatGraph(Prob_Bestfit_90, 2, 1, 2);
	  Prob_Bestfit_90->Draw("L");
	  lProbBestfit->AddEntry(Prob_Bestfit_90, "90%", "L");
        }

      lProbBestfit->SetBorderSize(0);
      lProbBestfit->Draw();

      //ProbBestfit->SaveAs("ProbBestfit.png");

    }//end of if(plot_bestfit_prob) 
  

  if(plot_68CL_prob)
    {
      TCanvas *Prob68CL = new TCanvas("Prob68CL", "Prob68CL", 800, 600);
      FormatCanvas(Prob68CL, mixing23_min, mixing23_max, dm23sq_min, dm23sq_max);
      TLegend *lProb68CL = new TLegend(0.7, 0.68, 0.85, 0.85);
      if(interpolatedCLPointsFound[0][0] == dm23sq_nbins)
	{
	  TGraph *Prob_68_10 = new TGraph(dm23sq_nbins, interpolatedCLProb[0][0], dm23sqPoints);
	  FormatGraph(Prob_68_10, 6, 1, 2);
	  Prob_68_10->Draw("L");
	  lProb68CL->AddEntry(Prob_68_10, "10%", "L");
	}
      if(interpolatedCLPointsFound[0][1] == dm23sq_nbins)
	{
	  TGraph *Prob_68_30 = new TGraph(dm23sq_nbins, interpolatedCLProb[0][1], dm23sqPoints);
	  FormatGraph(Prob_68_30, 1, 1, 2);
	  Prob_68_30->Draw("L");
	  lProb68CL->AddEntry(Prob_68_30, "30%", "L");
	}
      if(interpolatedCLPointsFound[0][2] == dm23sq_nbins)
	{
	  TGraph *Prob_68_50 = new TGraph(dm23sq_nbins, interpolatedCLProb[0][2], dm23sqPoints);
	  FormatGraph(Prob_68_50, 8, 1, 2);
	  Prob_68_50->Draw("L");
	  lProb68CL->AddEntry(Prob_68_50, "50%", "L");
	}
      if(interpolatedCLPointsFound[0][3] == dm23sq_nbins)
	{
	  TGraph *Prob_68_70 = new TGraph(dm23sq_nbins, interpolatedCLProb[0][3], dm23sqPoints);
	  FormatGraph(Prob_68_70, 4, 1, 2);
	  Prob_68_70->Draw("L");
	  lProb68CL->AddEntry(Prob_68_70, "70%", "L");
	}
      if(interpolatedCLPointsFound[0][4] == dm23sq_nbins)
        {
        TGraph *Prob_68_90 = new TGraph(dm23sq_nbins, interpolatedCLProb[0][4], dm23sqPoints);
        FormatGraph(Prob_68_90, 2, 1, 2);
        Prob_68_90->Draw("L");
        lProb68CL->AddEntry(Prob_68_90, "90%", "L");
        }

     lProb68CL->SetBorderSize(0);
     lProb68CL->Draw();

     Prob68CL->SaveAs("Prob68CL.png");

  }//end of if(plot_68CL_prob)

  if(plot_90CL_prob)
    {
    TCanvas *Prob90CL = new TCanvas("Prob90CL", "Prob90CL", 800, 600);
    FormatCanvas(Prob90CL, mixing23_min, mixing23_max, dm23sq_min, dm23sq_max);
    TLegend *lProb90CL = new TLegend(0.65, 0.63, 0.87, 0.87);
    if(interpolatedCLPointsFound[1][0] == dm23sq_nbins)
      {
      TGraph *Prob_90_10 = new TGraph(dm23sq_nbins, interpolatedCLProb[1][0], dm23sqPoints);
      FormatGraph(Prob_90_10, 6, 1, 2);
      Prob_90_10->Draw("L");
      lProb90CL->AddEntry(Prob_90_10, "10%", "L");
      }
    if(interpolatedCLPointsFound[1][1] == dm23sq_nbins)
      {
	TGraph *Prob_90_30 = new TGraph(dm23sq_nbins, interpolatedCLProb[1][1], dm23sqPoints);
	FormatGraph(Prob_90_30, 1, 1, 2);
	Prob_90_30->Draw("L");
	lProb90CL->AddEntry(Prob_90_30, "30%", "L");
      }
    if(interpolatedCLPointsFound[1][2] == dm23sq_nbins)
      {
	TGraph *Prob_90_50 = new TGraph(dm23sq_nbins, interpolatedCLProb[1][2], dm23sqPoints);
	FormatGraph(Prob_90_50, 8, 1, 2);
	Prob_90_50->Draw("L");
	lProb90CL->AddEntry(Prob_90_50, "50%", "L");
      }
    if(interpolatedCLPointsFound[1][3] == dm23sq_nbins)
      {
	TGraph *Prob_90_70 = new TGraph(dm23sq_nbins, interpolatedCLProb[1][3], dm23sqPoints);
	FormatGraph(Prob_90_70, 4, 1, 2);
	Prob_90_70->Draw("L");
	lProb90CL->AddEntry(Prob_90_70, "70%", "L");
      }
    if(interpolatedCLPointsFound[1][4] == dm23sq_nbins)
      {
	TGraph *Prob_90_90 = new TGraph(dm23sq_nbins, interpolatedCLProb[1][4], dm23sqPoints);
	FormatGraph(Prob_90_90, 2, 1, 2);
	Prob_90_90->Draw("L");
	lProb90CL->AddEntry(Prob_90_90, "90%", "L");
      }    
 
    lProb90CL->SetBorderSize(0);
    lProb90CL->Draw();

    Prob90CL->SaveAs("Prob90CL.png");

    }//end of if(plot_90CL_prob)

  if(plot_95CL_prob)
    {
      TCanvas *Prob95CL = new TCanvas("Prob95CL", "Prob95CL", 800, 600);
      FormatCanvas(Prob95CL, mixing23_min, mixing23_max, dm23sq_min, dm23sq_max);
      TLegend *lProb95CL = new TLegend(0.7, 0.68, 0.85, 0.85);
      if(interpolatedCLPointsFound[2][0] == dm23sq_nbins)
	{
	  TGraph *Prob_95_10 = new TGraph(dm23sq_nbins, interpolatedCLProb[2][0], dm23sqPoints);
	  FormatGraph(Prob_95_10, 6, 1, 2);
	  Prob_95_10->Draw("L");
	  lProb95CL->AddEntry(Prob_95_10, "10%", "L");
	}
      if(interpolatedCLPointsFound[2][1] == dm23sq_nbins)
	{
	  TGraph *Prob_95_30 = new TGraph(dm23sq_nbins, interpolatedCLProb[2][1], dm23sqPoints);
	  FormatGraph(Prob_95_30, 1, 1, 2);
	  Prob_95_30->Draw("L");
	  lProb95CL->AddEntry(Prob_95_30, "30%", "L");
	}
      if(interpolatedCLPointsFound[2][2] == dm23sq_nbins)
	{
	  TGraph *Prob_95_50 = new TGraph(dm23sq_nbins, interpolatedCLProb[2][2], dm23sqPoints);
	  FormatGraph(Prob_95_50, 8, 1, 2);
	  Prob_95_50->Draw("L");
	  lProb95CL->AddEntry(Prob_95_50, "50%", "L");
	}
      if(interpolatedCLPointsFound[2][3] == dm23sq_nbins)
	{
	  TGraph *Prob_95_70 = new TGraph(dm23sq_nbins, interpolatedCLProb[2][3], dm23sqPoints);
	  FormatGraph(Prob_95_70, 4, 1, 2);
	  Prob_95_70->Draw("L");
	  lProb95CL->AddEntry(Prob_95_70, "70%", "L");
	}
      if(interpolatedCLPointsFound[2][4] == dm23sq_nbins)
	{
	  TGraph *Prob_95_90 = new TGraph(dm23sq_nbins, interpolatedCLProb[2][4], dm23sqPoints);
	  FormatGraph(Prob_95_90, 2, 1, 2);
	  Prob_95_90->Draw("L");
	  lProb95CL->AddEntry(Prob_95_90, "90%", "L");
	}

      lProb95CL->SetBorderSize(0);
      lProb95CL->Draw();

      Prob95CL->SaveAs("Prob95CL.png");

    }//end of if(plot_95CL_prob)

      if(plot_99CL_prob)
	{
	  TCanvas *Prob99CL = new TCanvas("Prob99CL", "Prob99CL", 800, 600);
	  FormatCanvas(Prob99CL, mixing23_min, mixing23_max, dm23sq_min, dm23sq_max);
	  TLegend *lProb99CL = new TLegend(0.7, 0.68, 0.85, 0.85);
	  if(interpolatedCLPointsFound[3][0] == dm23sq_nbins)
	    {
	      TGraph *Prob_99_10 = new TGraph(dm23sq_nbins, interpolatedCLProb[3][0], dm23sqPoints);
	      FormatGraph(Prob_99_10, 6, 1, 2);
	      Prob_99_10->Draw("L");
	      lProb99CL->AddEntry(Prob_99_10, "10%", "L");
	    }
	  if(interpolatedCLPointsFound[3][1] == dm23sq_nbins)
	    {
	      TGraph *Prob_99_30 = new TGraph(dm23sq_nbins, interpolatedCLProb[3][1], dm23sqPoints);
	      FormatGraph(Prob_99_30, 1, 1, 2);
	      Prob_99_30->Draw("L");
	      lProb99CL->AddEntry(Prob_99_30, "30%", "L");
	    }
	  if(interpolatedCLPointsFound[3][2] == dm23sq_nbins)
	    {
	      TGraph *Prob_99_50 = new TGraph(dm23sq_nbins, interpolatedCLProb[3][2], dm23sqPoints);
	      FormatGraph(Prob_99_50, 8, 1, 2);
	      Prob_99_50->Draw("L");
	      lProb99CL->AddEntry(Prob_99_50, "50%", "L");
	    }
	  if(interpolatedCLPointsFound[3][3] == dm23sq_nbins)
	    {
	      TGraph *Prob_99_70 = new TGraph(dm23sq_nbins, interpolatedCLProb[3][3], dm23sqPoints);
	      FormatGraph(Prob_99_70, 4, 1, 2);
	      Prob_99_70->Draw("L");
	      lProb99CL->AddEntry(Prob_99_70, "70%", "L");
	    }
          if(interpolatedCLPointsFound[3][4] == dm23sq_nbins)
	    {
	      TGraph *Prob_99_90 = new TGraph(dm23sq_nbins, interpolatedCLProb[3][4], dm23sqPoints);
	      FormatGraph(Prob_99_90, 2, 1, 2);
	      Prob_99_90->Draw("L");
	      lProb99CL->AddEntry(Prob_99_90, "90%", "L");
	    }

	  lProb99CL->SetBorderSize(0);
	  lProb99CL->Draw();

          Prob99CL->SaveAs("Prob99CL.png");

	}//end of if(plot_99CL_prob)        

}//end of NonMax23Discovery()

void FormatCanvas(TCanvas *canvas, double mixing23_min, double mixing23_max, double dm23sq_min, double dm23sq_max)
{
  TGaxis::SetMaxDigits(2);

  canvas->Divide(1);
  canvas->cd(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.12);

  //dummyGraph is made to draw the axes                                                                                                                                                          
  //this is done since it is not known which of the probability graphs will be drawn
  double dummyX[2] = {30.0, 31.0};
  double dummyY[2] = {30.0, 31.0};
  TGraph *dummyGraph = new TGraph(2, dummyX, dummyY);
  dummyGraph->SetTitle("");
  dummyGraph->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  dummyGraph->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2})");
  dummyGraph->GetXaxis()->SetLimits(0.85, 1.01);
  dummyGraph->GetYaxis()->SetRangeUser(0.95*dm23sq_min, 1.05*dm23sq_max);
  dummyGraph->Draw("AL");
  dummyGraph->GetXaxis()->SetTitleOffset(1.15);
  dummyGraph->GetYaxis()->SetTitleOffset(1.4);
  dummyGraph->GetXaxis()->CenterTitle();
  dummyGraph->GetYaxis()->CenterTitle();
}

void FormatGraph(TGraph *graph, int colour, int style, int width)
{
  graph->SetTitle("");
  graph->SetLineColor(colour);
  graph->SetLineStyle(style);
  graph->SetLineWidth(width);
}
