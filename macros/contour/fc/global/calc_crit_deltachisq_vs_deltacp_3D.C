#include <vector>

typedef enum EConfLevel {
  kCLUndefined = -1, 
  kCL68, 
  kCL90,
  kCL95, 
  kCL99, 
};

//const int nConfidenceLevel = 12;
//double confidenceLevel[nConfidenceLevel] = {68.0, 90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0};
const int nConfidenceLevel = 3;
double confidenceLevel[nConfidenceLevel] = {
  68.0, 
  90.0,
  95.0,
  //99.0,
};


void calc_crit_deltachisq_vs_deltacp_3D
(
 char FCfilename[200] = "prova_idcp_%d.root",
 bool SingleBin = true,
 int lastbin = 6,
 bool CheckConvergence = true,
 int StepConvergence = 100,
 const char tag[200] = "critval",
 )
{
  cout<<setprecision(3);
  
  //

  int startbin = 0;
  if(SingleBin){
    startbin = lastbin;
    cout << endl;
    cout << "Study 1 single delta_cp bin --> " << lastbin << endl;
  }  
  const int nbins = lastbin - startbin + 1;

  // Declare global variables

  vector<double> deltaChisq[nbins];
  int nbadfits=0;
  int negdchisq=0;
  ofstream outfile;
  char outfilename[200];
  

  // Declare histograms, graphs, canvases

  TH1D *hDeltaChisq[nbins];
  for(int ibin=0;ibin<nbins;ibin++){
    TString name = Form("hDeltaChisq_%d",ibin);
    hDeltaChisq[ibin] = new TH1D(name,name, 150, -5., 10.);
  }

  TGraph *gDLcritConv[nConfidenceLevel][nbins]; 
  for(int icl=0;icl<nConfidenceLevel;icl++){
    for(int ibin=0;ibin<nbins;ibin++){
      gDLcritConv[icl][ibin] = new TGraph();
      TString name = Form("gDLcritConv_%d",ibin);
      gDLcritConv[icl][ibin]->SetName(name);
    }
  }
 
  TCanvas *cDeltaChisq[nbins];
  for(int ibin=0;ibin<nbins;ibin++){
    TString name = Form("cDeltaChisq_%d",ibin);
    cDeltaChisq[ibin] = new TCanvas(name,name);
    hDeltaChisq[ibin]->Draw("histo");
  }
  
  TCanvas *cDLConv_68[nbins];
  TCanvas *cDLConv_90[nbins];
  TCanvas *cDLConv_95[nbins];
  for(int ibin=0;ibin<nbins;ibin++){
    TString name = Form("cDLConv_68_%d",ibin);
    cDLConv_68[ibin] = new TCanvas(name,name);
    TString name = Form("cDLConv_90_%d",ibin);
    cDLConv_90[ibin] = new TCanvas(name,name);
    TString name = Form("cDLConv_95_%d",ibin);
    cDLConv_95[ibin] = new TCanvas(name,name);
  }
  

  //TH1D *hInputMixing23 = new TH1D("hInputMixing23", "", 101, 0.25, 0.75); 
  //TH1D *hInputDm23sq = new TH1D("hInputDm23sq", "", 41, 0.0015, 0.0035);
  //TH2D *hInputBestFit1NMixing23 = new TH2D("hInputBestFit1NMixing23", "", 101, 0.25, 0.75, 101, 0.25, 0.75);
  //TH2D *hInputBestFit2NMixing23 = new TH2D("hInputBestFit2NMixing23", "", 101, 0.25, 0.75, 101, 0.25, 0.75);
  //TH2D *hInputBestFit1NDm23sq = new TH2D("hInputBestFit1NDm23sq", "", 41, 0.0015, 0.0035, 41, 0.0015, 0.0035);
  //TH2D *hInputBestFit2NDm23sq = new TH2D("hInputBestFit2NDm23sq", "", 41, 0.0015, 0.0035, 41, 0.0015, 0.0035);

    
  // Output logfile

  TString outname="";
  outname.Append(tag);
  outname.Append(".log");
  sprintf(outfilename,"%s",(char*)outname.Data());
  outfile.open(outname.Data());
  if (!outfile.is_open()) {    
    std::cout << "Error opening file" << endl;
    std::cout << outfilename << std::endl; 
    exit(0);
  }

  //

  for(int ifile = startbin; ifile<=lastbin; ifile++){ // loop over delta_cp files    

    int ibin = ifile - startbin; // delta_cp bins for arrays and histograms
    

    cout << ibin << " " << ifile << " " << startbin << " " << lastbin << endl;


    // Check if the file exists and the tree is not corrupted. If yes skip it

    TString infilename = Form(FCfilename,ifile); 
    bool DLcritFlag = true;
    TFile *MyFile = new TFile(infilename,"READ");      
    if(!MyFile || MyFile->IsOpen()==false){
      cout << "        File " << infilename << " not available!!! --> SKIP IT!!!" << endl << endl;	
      outfile << "        File " << infilename << " not available!!! --> SKIP IT!!!" << endl << endl;	
      DLcritFlag = false;

      MyFile->Close();
      MyFile->Delete();
    } 
    if(!DLcritFlag) continue;

    TTree *fitChain = (TTree*)MyFile->Get("fitnt");
    if(!fitChain){
      cout << "        Problem in reading the tree --> SKIP IT!!!" << endl << endl;
      outfile << "        Problem in reading the tree --> SKIP IT!!!" << endl << endl;
      DLcritFlag = false;
      MyFile->Close();
      MyFile->Delete();
    }
    if(!DLcritFlag) continue;

    //

    int nsk_1Rmu,nsk_1Re;
    int input_mixing23_bin,input_mixing13_bin, input_dmsq_bin;
    double input_mixing23,input_mixing13,input_dmsq,input_deltacp;
    int fittries_true, fittries_true_alt, fittries_min, fittries_min_oct, fittries_min_mh,fittries_min_mh_oct;
    double fittime_true, fittime_true_alt, fittime_min, fittime_min_oct, fittime_min_mh,fittime_min_mh_oct;
    int migradmin_true, migradmin_true_alt, migradmin_min, migradmin_min_oct, migradmin_min_mh,migradmin_min_mh_oct;
    int migraderr_true, migraderr_true_alt, migraderr_min, migraderr_min_oct, migraderr_min_mh,migraderr_min_mh_oct;
    int migradmin1, migradmin2, migradmin3, migradmin4, migradmin5;
    int migraderr1, migraderr2, migraderr3, migraderr4, migraderr5; 
    double bestfit_chisq_seed1, bestfit_chisq_seed2, bestfit_chisq_seed3,bestfit_chisq_seed4;
    double sinsqtheta23_bestfit, sinsqtheta13_bestfit, centre, reseed;
    double bestfit_chisq_true, bestfit_chisq_true_alt, bestfit_chisq_min, bestfit_chisq_min_oct, bestfit_chisq_min_mh,bestfit_chisq_min_mh_oct;
    double bestfit_mixing23_true, bestfit_mixing23_true_alt, bestfit_mixing23_min, bestfit_mixing23_min_oct, bestfit_mixing23_min_mh,bestfit_mixing23_min_mh_oct;
    double bestfit_mixing13_true, bestfit_mixing13_true_alt, bestfit_mixing13_min, bestfit_mixing13_min_oct, bestfit_mixing13_min_mh,bestfit_mixing13_min_mh_oct;
    double bestfit_dmsq_true, bestfit_dmsq_true_alt, bestfit_dmsq_min, bestfit_dmsq_min_oct, bestfit_dmsq_min_mh,bestfit_dmsq_min_mh_oct;
    double bestfit_deltacp_true, bestfit_deltacp_true_alt, bestfit_deltacp_min, bestfit_deltacp_min_oct, bestfit_deltacp_min_mh,bestfit_deltacp_min_mh_oct;

    fitChain->SetBranchAddress("nsk_1Rmu",                  &nsk_1Rmu);
    fitChain->SetBranchAddress("nsk_1Re",                   &nsk_1Re);
    fitChain->SetBranchAddress("input_mixing23_bin",        &input_mixing23_bin);
    fitChain->SetBranchAddress("input_mixing13_bin",        &input_mixing13_bin);
    fitChain->SetBranchAddress("input_dmsq_bin",            &input_dmsq_bin);
    fitChain->SetBranchAddress("input_mixing23",            &input_mixing23);
    fitChain->SetBranchAddress("input_mixing13",            &input_mixing13);
    fitChain->SetBranchAddress("input_dmsq",                &input_dmsq);
    fitChain->SetBranchAddress("input_deltacp",             &input_deltacp);
    fitChain->SetBranchAddress("fittries_true",             &fittries_true);
    fitChain->SetBranchAddress("fittries_true_alt",         &fittries_true_alt);
    fitChain->SetBranchAddress("fittries_min",              &fittries_min);
    fitChain->SetBranchAddress("fittries_min_oct",          &fittries_min_oct);
    fitChain->SetBranchAddress("fittries_min_mh",           &fittries_min_mh);
    fitChain->SetBranchAddress("fittries_min_mh_oct",       &fittries_min_mh_oct);
    fitChain->SetBranchAddress("fittime_true",              &fittime_true);
    fitChain->SetBranchAddress("fittime_true_alt",          &fittime_true_alt);
    fitChain->SetBranchAddress("fittime_min",               &fittime_min);
    fitChain->SetBranchAddress("fittime_min_oct",           &fittime_min_oct);
    fitChain->SetBranchAddress("fittime_min_mh",            &fittime_min_mh);
    fitChain->SetBranchAddress("fittime_min_mh_oct",        &fittime_min_mh_oct);
    fitChain->SetBranchAddress("migradmin_true",            &migradmin_true);
    fitChain->SetBranchAddress("migradmin_true_alt",        &migradmin_true_alt);
    fitChain->SetBranchAddress("migradmin_min",             &migradmin_min);
    fitChain->SetBranchAddress("migradmin_min_oct",         &migradmin_min_oct);
    fitChain->SetBranchAddress("migradmin_min_mh",          &migradmin_min_mh);
    fitChain->SetBranchAddress("migradmin_min_mh_oct",      &migradmin_min_mh_oct);
    fitChain->SetBranchAddress("migraderr_true",            &migraderr_true);
    fitChain->SetBranchAddress("migraderr_true_alt",        &migraderr_true_alt);
    fitChain->SetBranchAddress("migraderr_min",             &migraderr_min);
    fitChain->SetBranchAddress("migraderr_min_oct",         &migraderr_min_oct);
    fitChain->SetBranchAddress("migraderr_min_mh",          &migraderr_min_mh);
    fitChain->SetBranchAddress("migraderr_min_mh_oct",      &migraderr_min_mh_oct);
    fitChain->SetBranchAddress("bestfit_chisq_true",        &bestfit_chisq_true);
    fitChain->SetBranchAddress("bestfit_chisq_true_alt",    &bestfit_chisq_true_alt);
    fitChain->SetBranchAddress("bestfit_chisq_min",         &bestfit_chisq_min);
    fitChain->SetBranchAddress("bestfit_chisq_min_oct",     &bestfit_chisq_min_oct);
    fitChain->SetBranchAddress("bestfit_chisq_min_mh",      &bestfit_chisq_min_mh);
    fitChain->SetBranchAddress("bestfit_chisq_min_mh_oct",  &bestfit_chisq_min_mh_oct);
    fitChain->SetBranchAddress("bestfit_mixing23_true",     &bestfit_mixing23_true);
    fitChain->SetBranchAddress("bestfit_mixing23_true_alt", &bestfit_mixing23_true_alt);
    fitChain->SetBranchAddress("bestfit_mixing23_min",      &bestfit_mixing23_min);
    fitChain->SetBranchAddress("bestfit_mixing23_min_oct",  &bestfit_mixing23_min_oct);
    fitChain->SetBranchAddress("bestfit_mixing23_min_mh",   &bestfit_mixing23_min_mh);
    fitChain->SetBranchAddress("bestfit_mixing23_min_mh_oct",  &bestfit_mixing23_min_mh_oct);
    fitChain->SetBranchAddress("bestfit_mixing13_true",     &bestfit_mixing13_true);
    fitChain->SetBranchAddress("bestfit_mixing13_true_alt", &bestfit_mixing13_true_alt);
    fitChain->SetBranchAddress("bestfit_mixing13_min",      &bestfit_mixing13_min);
    fitChain->SetBranchAddress("bestfit_mixing13_min_oct",  &bestfit_mixing13_min_oct);
    fitChain->SetBranchAddress("bestfit_mixing13_min_mh",   &bestfit_mixing13_min_mh);
    fitChain->SetBranchAddress("bestfit_mixing13_min_mh_oct",  &bestfit_mixing13_min_mh_oct);
    fitChain->SetBranchAddress("bestfit_dmsq_true",         &bestfit_dmsq_true);
    fitChain->SetBranchAddress("bestfit_dmsq_true_alt",     &bestfit_dmsq_true_alt);
    fitChain->SetBranchAddress("bestfit_dmsq_min",          &bestfit_dmsq_min);
    fitChain->SetBranchAddress("bestfit_dmsq_min_oct",      &bestfit_dmsq_min_oct);
    fitChain->SetBranchAddress("bestfit_dmsq_min_mh",       &bestfit_dmsq_min_mh);
    fitChain->SetBranchAddress("bestfit_dmsq_min_mh_oct",   &bestfit_dmsq_min_mh_oct);
    fitChain->SetBranchAddress("bestfit_deltacp_true",      &bestfit_deltacp_true);
    fitChain->SetBranchAddress("bestfit_deltacp_true_alt",  &bestfit_deltacp_true_alt);
    fitChain->SetBranchAddress("bestfit_deltacp_min",       &bestfit_deltacp_min);
    fitChain->SetBranchAddress("bestfit_deltacp_min_oct",   &bestfit_deltacp_min_oct);
    fitChain->SetBranchAddress("bestfit_deltacp_min_mh",    &bestfit_deltacp_min_mh);
    fitChain->SetBranchAddress("bestfit_deltacp_min_mh_oct",  &bestfit_deltacp_min_mh_oct);

    //

    fitChain->GetEntry(0);

    cout << endl;  
    outfile << endl;  
    cout << "///////////////////    delta_cp = " << input_deltacp << "   /////////////////////" << endl;
    outfile << "///////////////////    delta_cp = " << input_deltacp << "   /////////////////////" << endl;

    cout << endl;
    outfile << endl;
    cout << "delta_cp point " << ifile << " --> File: " << infilename << endl;
    outfile << "delta_cp point " << ifile << " --> File: " << infilename << endl;

    const int nentries = fitChain->GetEntries();
    cout << "# of toys = " << nentries << endl;
    outfile << "# of toys = " << nentries << endl;
    cout << endl;
    outfile << endl;
    
    // 

    for(Long64_t i = 0; i < nentries; i++) {       

      fitChain->GetEntry(i);
      
      if((migradmin_true == 0 && migraderr_true == 3 && migradmin_min == 0 && migraderr_min == 3) || (migradmin_true == 0 && migraderr_true == 3 && migradmin_min == -1 && migraderr_min == 0)) {
	double dchisq = bestfit_chisq_true - bestfit_chisq_min;
	deltaChisq[ibin].push_back(dchisq);
	hDeltaChisq[ibin]->Fill(bestfit_chisq_true - bestfit_chisq_min);
	
	if(dchisq<0.) negdchisq++;
	
	if(CheckConvergence){
	  if(i%StepConvergence==0){
	    for(int icl=0;icl<nConfidenceLevel;icl++){    
	      double CurrCritVal = CalcCritVal(deltaChisq[ibin],confidenceLevel[icl],0.);
	      gDLcritConv[icl][ibin]->SetPoint(gDLcritConv[icl][ibin]->GetN(),gDLcritConv[icl][ibin]->GetN()*StepConvergence,CurrCritVal);
	    }
	  }	
	}
      }
      else {
	cout << endl;
	cout << "Fit ended w/o convergence or fit problems:" << endl;
	cout << " migradmin_true = " << migradmin_true << " migraderr_true = " << migraderr_true << " migradmin_min = " << migradmin_min << " migraderr_min = " << migraderr_min << endl;
	outfile << endl;
	outfile << "Fit ended w/o convergence or fit problems:" << endl;
	outfile << " migradmin_true = " << migradmin_true << " migraderr_true = " << migraderr_true << " migradmin_min = " << migradmin_min << " migraderr_min = " << migraderr_min << endl;
	nbadfits++;
      }
   
    }//end of loop over entries

    // Print output

    cout << "Print outputs:" << endl;
    cout << endl;
    cout << "# of bad fits = " << nbadfits << endl;
    cout << "Frac negative dchisq = " << negdchisq / deltaChisq[ibin].size() << endl;
    cout << "deltaChisq size = " << deltaChisq[ibin].size() << endl;
    cout << endl;
    outfile << "Print outputs:" << endl;
    outfile << endl;
    outfile << "# of bad fits = " << nbadfits << endl;
    outfile << "Frac negative dchisq = " << negdchisq / deltaChisq[ibin].size() << endl;
    outfile << "deltaChisq size = " << deltaChisq[ibin].size() << endl;
    outfile << endl;
  
    // Get binomial error on critical values based on toy statistics
    
    int sigmaCL[nConfidenceLevel];
    double prob=0.;
    for(int icl=0;icl<nConfidenceLevel;icl++){    
      prob = confidenceLevel[icl] / 100.;
      sigmaCL[icl] = (int) (sqrt(deltaChisq[ibin].size()*prob*(1-prob))) + 0.5;
      cout << icl << " --> sigmaCL ( " << confidenceLevel[icl] << "CL ) = " << sigmaCL[icl] << endl;
      outfile << icl << " --> sigmaCL ( " << confidenceLevel[icl] << "CL ) = " << sigmaCL[icl] << endl;
    }

    // Get critical values

    cout << endl;
    outfile << endl;
    for(int icl=0; icl<nConfidenceLevel; icl++) {
   
      double value = CalcCritVal(deltaChisq[ibin],confidenceLevel[icl],0.);
      double valuep = CalcCritVal(deltaChisq[ibin],confidenceLevel[icl],sigmaCL[icl]);
      double errorpos = valuep - value;
      double valuem = CalcCritVal(deltaChisq[ibin],confidenceLevel[icl],-sigmaCL[icl]);
      double errorneg = value - valuem;
     
      cout << "Critical value of delta chisq at "<<confidenceLevel[icl]<<"% confidence = "<< value
	   << " + " << errorpos << "  -  " << errorneg     
	   << endl << endl; 
      outfile << "Critical value of delta chisq at "<<confidenceLevel[icl]<<"% confidence = "<< value
	      << " + " << errorpos << "  -  " << errorneg     
	      << endl << endl; 
    }
    
    // Draw histograms
    
    hDeltaChisq[ibin]->Scale(1/hDeltaChisq[ibin]->Integral());
    hDeltaChisq[ibin]->GetXaxis()->SetTitle("#Delta#chi^{2}");
    hDeltaChisq[ibin]->GetYaxis()->SetTitle("Fraction of events");
    hDeltaChisq[ibin]->Draw("");
  
    // Check the convergence of the critical values
    if(CheckConvergence){
      cDLConv_68[ibin]->cd();
      gDLcritConv[kCL68][ibin]->Draw("APL");
      cDLConv_90[ibin]->cd();
      gDLcritConv[kCL90][ibin]->Draw("APL");
      cDLConv_95[ibin]->cd();
      gDLcritConv[kCL95][ibin]->Draw("APL");
    }

    /*
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

      deltaChisq[ibin].clear();
    */

    //
    
    MyFile->Close();
    MyFile->Delete();

  } // end loop deltacp bins


  // Write output file
  
  TString outname="";
  outname.Append(tag);
  outname.Append(".root");
  sprintf(outfilename,"%s",(char*)outname.Data());
  TFile *InputValues = new TFile(outfilename, "RECREATE");
  for(int ibin=0;ibin<nbins;ibin++){
    hDeltaChisq[ibin]->Write();
  }

  return;

}//end of calc_crit_dchisq_mixing23_test



double CalcCritVal(vector<double> dchisq, double conflev, double sigmacl){
  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector   
  vector<double>::iterator dchisqIter;
  dchisqIter = dchisq.begin() + (int)(0.01*conflev*dchisq.size() + sigmacl);
  double critval = (*dchisqIter);
  return critval;
}
 
