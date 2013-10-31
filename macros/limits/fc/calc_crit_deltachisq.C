#include <vector>

void calc_crit_deltachisq
(
  const char * oscparameter             = "mixing23",
  int bin                               = 26
)
{
  const int nConfidenceLevel = 2;

  double confidenceLevel[nConfidenceLevel] = {68.0, 90.0};

  cout<<setprecision(3);

  int iexp, irun, input_mixing23_bin, input_dm23sq_bin;
  int migradmin_mixing23_1N, migradmin_mixing23_1N_alt, migradmin_dm23sq_1N, migraderr_mixing23_1N, migraderr_mixing23_1N_alt, migraderr_dm23sq_1N;
  int migradmin_2N, migradmin_2N_alt, migraderr_2N, migraderr_2N_alt;
  int fittries_mixing23_1N, fittries_mixing23_1N_alt, fittries_dm23sq_1N, fittries_2N, fittries_2N_alt;
  double fittime_mixing23_1N, fittime_mixing23_1N_alt, fittime_dm23sq_1N, fittime_2N, fittime_2N_alt;
  double input_mixing23, input_dm23sq, bestfit_mixing23_1N, bestfit_mixing23_1N_alt, bestfit_dm23sq_1N;
  double bestfit_chisq_mixing23_1N, bestfit_chisq_mixing23_1N_alt, bestfit_chisq_dm23sq_1N;
  double bestfit_mixing23_2N, bestfit_mixing23_2N_alt, bestfit_dm23sq_2N, bestfit_dm23sq_2N_alt, bestfit_chisq_2N, bestfit_chisq_2N_alt;

  if(oscparameter != "mixing23" && oscparameter != "dm23sq")
    {
    cerr<<"oscparameter must be either mixing23 or dm23sq - exiting"<<endl;
    exit(1);
    }  

  TChain *fitChain = new TChain("fitnt");
  fitChain->Add(Form("/opt/ppd/t2k/users/nick/t2kana1DFCRun1to4/production/%sbin%02d/fc-run_*.root", oscparameter, bin));

  TString filename = Form("/opt/ppd/t2k/users/nick/t2kana1DFCRun1to4/production/%sbin%02d/fc-run_*.root", oscparameter, bin);
  cout<<filename<<endl;

  fitChain->SetBranchAddress("iexp",&iexp);
  fitChain->SetBranchAddress("irun",&irun);
  fitChain->SetBranchAddress("input_mixing23_bin",&input_mixing23_bin);
  fitChain->SetBranchAddress("input_dm23sq_bin",&input_dm23sq_bin);
  fitChain->SetBranchAddress("input_mixing23",&input_mixing23);
  fitChain->SetBranchAddress("input_dm23sq",&input_dm23sq);
  fitChain->SetBranchAddress("bestfit_mixing23_1N",&bestfit_mixing23_1N);
  fitChain->SetBranchAddress("bestfit_mixing23_1N_alt",&bestfit_mixing23_1N_alt);
  fitChain->SetBranchAddress("bestfit_dm23sq_1N",&bestfit_dm23sq_1N);
  fitChain->SetBranchAddress("bestfit_chisq_mixing23_1N",&bestfit_chisq_mixing23_1N);
  fitChain->SetBranchAddress("bestfit_chisq_mixing23_1N_alt",&bestfit_chisq_mixing23_1N_alt);
  fitChain->SetBranchAddress("bestfit_chisq_dm23sq_1N",&bestfit_chisq_dm23sq_1N);
  fitChain->SetBranchAddress("bestfit_mixing23_2N",&bestfit_mixing23_2N);
  fitChain->SetBranchAddress("bestfit_mixing23_2N_alt",&bestfit_mixing23_2N_alt);
  fitChain->SetBranchAddress("bestfit_dm23sq_2N",&bestfit_dm23sq_2N);
  fitChain->SetBranchAddress("bestfit_dm23sq_2N_alt",&bestfit_dm23sq_2N_alt);
  fitChain->SetBranchAddress("bestfit_chisq_2N",&bestfit_chisq_2N);
  fitChain->SetBranchAddress("bestfit_chisq_2N_alt",&bestfit_chisq_2N_alt);
  fitChain->SetBranchAddress("fittries_mixing23_1N",&fittries_mixing23_1N);
  fitChain->SetBranchAddress("fittries_mixing23_1N_alt",&fittries_mixing23_1N_alt);
  fitChain->SetBranchAddress("fittries_dm23sq_1N",&fittries_dm23sq_1N);
  fitChain->SetBranchAddress("fittries_2N",&fittries_2N);
  fitChain->SetBranchAddress("fittries_2N_alt",&fittries_2N_alt);
  fitChain->SetBranchAddress("fittime_mixing23_1N",&fittime_mixing23_1N);
  fitChain->SetBranchAddress("fittime_mixing23_1N_alt",&fittime_mixing23_1N_alt);
  fitChain->SetBranchAddress("fittime_dm23sq_1N",&fittime_dm23sq_1N);
  fitChain->SetBranchAddress("fittime_2N",&fittime_2N);
  fitChain->SetBranchAddress("fittime_2N_alt",&fittime_2N_alt);
  fitChain->SetBranchAddress("migradmin_mixing23_1N",&migradmin_mixing23_1N);
  fitChain->SetBranchAddress("migradmin_mixing23_1N_alt",&migradmin_mixing23_1N_alt);
  fitChain->SetBranchAddress("migradmin_dm23sq_1N",&migradmin_dm23sq_1N);
  fitChain->SetBranchAddress("migraderr_mixing23_1N",&migraderr_mixing23_1N);
  fitChain->SetBranchAddress("migraderr_mixing23_1N_alt",&migraderr_mixing23_1N_alt);
  fitChain->SetBranchAddress("migraderr_dm23sq_1N",&migraderr_dm23sq_1N);
  fitChain->SetBranchAddress("migradmin_2N",&migradmin_2N);
  fitChain->SetBranchAddress("migradmin_2N_alt",&migradmin_2N_alt);
  fitChain->SetBranchAddress("migraderr_2N",&migraderr_2N);
  fitChain->SetBranchAddress("migraderr_2N_alt",&migraderr_2N_alt);
  
  vector<double> deltaChisq;

  TH1D *hInputMixing23 = new TH1D("hInputMixing23", "", 101, 0.25, 0.75); 
  TH1D *hInputDm23sq = new TH1D("hInputDm23sq", "", 41, 0.0015, 0.0035);
  TH2D *hInputBestFit1NMixing23 = new TH2D("hInputBestFit1NMixing23", "", 101, 0.25, 0.75, 101, 0.25, 0.75);
  TH2D *hInputBestFit2NMixing23 = new TH2D("hInputBestFit2NMixing23", "", 101, 0.25, 0.75, 101, 0.25, 0.75);
  TH2D *hInputBestFit1NDm23sq = new TH2D("hInputBestFit1NDm23sq", "", 41, 0.0015, 0.0035, 41, 0.0015, 0.0035);
  TH2D *hInputBestFit2NDm23sq = new TH2D("hInputBestFit2NDm23sq", "", 41, 0.0015, 0.0035, 41, 0.0015, 0.0035);

  for(Long64_t i = 0; i < fitChain->GetEntries(); i++)
    {       
    fitChain->GetEntry(i);

    if(oscparameter == "mixing23")
      {
      if(migradmin_dm23sq_1N == 0 && migraderr_dm23sq_1N == 3 && migradmin_2N == 0 && migraderr_2N == 3)
	{
	deltaChisq.push_back(bestfit_chisq_dm23sq_1N - bestfit_chisq_2N);
	hInputMixing23->Fill(input_mixing23);
        hInputDm23sq->Fill(input_dm23sq);
        hInputBestFit1NMixing23->Fill(input_mixing23, bestfit_mixing23_1N);
	hInputBestFit2NMixing23->Fill(input_mixing23, bestfit_mixing23_2N);
        hInputBestFit1NDm23sq->Fill(input_dm23sq, bestfit_dm23sq_1N);
        hInputBestFit2NDm23sq->Fill(input_dm23sq, bestfit_dm23sq_2N);
	}
      }
    else if(oscparameter == "dm23sq")
      {
      if(migradmin_mixing23_1N == 0 && migraderr_mixing23_1N == 3 && migradmin_2N == 0 && migraderr_2N == 3)
	{
	deltaChisq.push_back(bestfit_chisq_mixing23_1N - bestfit_chisq_2N);
	hInputMixing23->Fill(input_mixing23);
	hInputDm23sq->Fill(input_dm23sq);
        hInputBestFit1NMixing23->Fill(input_mixing23, bestfit_mixing23_1N);
        hInputBestFit2NMixing23->Fill(input_mixing23, bestfit_mixing23_2N);
        hInputBestFit1NDm23sq->Fill(input_dm23sq, bestfit_dm23sq_1N);
        hInputBestFit2NDm23sq->Fill(input_dm23sq, bestfit_dm23sq_2N);
        }
      }

    }//end of loop over entries

  cout<<endl;
  cout<<"deltaChisq size = "<<deltaChisq.size()<<endl;
  cout<<endl;

  sort(deltaChisq.begin(), deltaChisq.end());
   
  vector<double>::iterator deltaChisqIter;

  /*
  deltaChisqIter = deltaChisq.begin();
  while(deltaChisqIter != deltaChisq.end())
    {
      cout<<"    "<<(*deltaChisqIter);
      ++deltaChisqIter;
    }

  cout<<endl;
  cout<<endl;
  */

  for(int j=0; j<nConfidenceLevel; j++)
    {
    deltaChisqIter = deltaChisq.begin() + (int)(0.01*confidenceLevel[j]*deltaChisq.size());
    cout<<endl;
    cout<<"Critical value of delta chisq at "<<confidenceLevel[j]<<"% confidence = "<<(*deltaChisqIter)<<endl;
    cout<<endl;  
    }

  TCanvas *InputMixing23 = new TCanvas("InputMixing23", "InputMixing23", 800, 600);
  hInputMixing23->SetStats(0);
  hInputMixing23->GetXaxis()->SetTitle("Input value of sin^{2}#theta_{23}");
  hInputMixing23->GetXaxis()->SetTitleOffset(1.15);
  hInputMixing23->GetXaxis()->CenterTitle();
  hInputMixing23->Draw();

  TCanvas *InputDm23sq = new TCanvas("InputDm23sq", "InputDm23sq", 800, 600);
  hInputDm23sq->SetStats(0);
  hInputDm23sq->GetXaxis()->SetTitle("Input value of |#Deltam_{32}^{2}|");
  hInputDm23sq->GetXaxis()->SetTitleOffset(1.15);
  hInputDm23sq->GetXaxis()->CenterTitle();
  hInputDm23sq->Draw();

  TCanvas *InputBestFit1NMixing23 = new TCanvas("InputBestFit1NMixing23", "InputBestFit1NMixing23", 800, 600);
  hInputBestFit1NMixing23->SetStats(0);
  hInputBestFit1NMixing23->GetXaxis()->SetTitle("Input sin^{2}#theta_{23}");
  hInputBestFit1NMixing23->GetYaxis()->SetTitle("Best-fit sin^{2}#theta_{23} (1+N)");   
  hInputBestFit1NMixing23->GetXaxis()->CenterTitle();
  hInputBestFit1NMixing23->GetYaxis()->CenterTitle();
  hInputBestFit1NMixing23->Draw("colz");  

  TCanvas *InputBestFit2NMixing23 = new TCanvas("InputBestFit2NMixing23", "InputBestFit2NMixing23", 800, 600);
  hInputBestFit2NMixing23->SetStats(0);
  hInputBestFit2NMixing23->GetXaxis()->SetTitle("Input sin^{2}#theta_{23}");
  hInputBestFit2NMixing23->GetYaxis()->SetTitle("Best-fit sin^{2}#theta_{23} (2+N)");
  hInputBestFit2NMixing23->GetXaxis()->CenterTitle();
  hInputBestFit2NMixing23->GetYaxis()->CenterTitle();
  hInputBestFit2NMixing23->Draw("colz");

  TCanvas *InputBestFit1NDm23sq = new TCanvas("InputBestFit1NDm23sq", "InputBestFit1NDm23sq", 800, 600);
  hInputBestFit1NDm23sq->SetStats(0);
  hInputBestFit1NDm23sq->GetXaxis()->SetTitle("Input #Deltam_{32}^{2}");
  hInputBestFit1NDm23sq->GetXaxis()->SetTitle("Best-fit #Deltam_{32}^{2} (1+N)");
  hInputBestFit1NDm23sq->GetYaxis()->CenterTitle();
  hInputBestFit1NDm23sq->GetYaxis()->CenterTitle();
  hInputBestFit1NDm23sq->Draw("colz");

  TCanvas *InputBestFit2NDm23sq = new TCanvas("InputBestFit2NDm23sq", "InputBestFit2NDm23sq", 800, 600);
  hInputBestFit2NDm23sq->SetStats(0);
  hInputBestFit2NDm23sq->GetXaxis()->SetTitle("Input #Deltam_{32}^{2}");
  hInputBestFit2NDm23sq->GetYaxis()->SetTitle("Best-fit #Deltam_{32}^{2} (2+N)");
  hInputBestFit2NDm23sq->GetXaxis()->CenterTitle();
  hInputBestFit2NDm23sq->GetYaxis()->CenterTitle();
  hInputBestFit2NDm23sq->Draw("colz");

  TFile *InputValues = new TFile("InputValues.root", "UPDATE");
  InputValues->cd();
  hInputDm23sq->Write();
  InputValues->Close();

  deltaChisq.clear();

}//end of calc_crit_dchisq_mixing23_test
