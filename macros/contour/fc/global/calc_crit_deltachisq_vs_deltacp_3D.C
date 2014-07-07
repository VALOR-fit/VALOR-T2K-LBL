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
 char FCfilename[200]        = "../fc_dcp_allsyst_locmin_3D_ntoys_1000_idcp_%d.root",
 //char FCfilename[200]        = "../fc_dcp_allsyst_3D_ntoys_1000_idcp_%d.root",
 bool SingleBin              = true,
 int lastbin                 = 8,
 bool CheckConvergence       = false,
 int StepConvergence         = 100,
 const char tag[200]         = "critval",
 bool GetAsimov              = false,
 char Asimovfilename[200]    = "Asimov_allsyst_NH_ntoys_1_idcp_%d_run_0.root",
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  cout<<setprecision(3);
  
  double Ngod_Min = 0.;
  double Ngod_Max = 10.;
  const int Ngod_Bins = 75;
  double SetPar_ndf = 0.4;

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
  double asimov_mu, asimov_sigma;
  double mu[nbins], sigma[nbins];

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
      TString name = Form("gDLcritConv_bin%d_cl%2.0f",ibin,confidenceLevel[icl]);
      gDLcritConv[icl][ibin]->SetName(name);
    }
  }

  TGraph *gCritVal[nConfidenceLevel];
  for(int icl=0;icl<nConfidenceLevel;icl++){
    gCritVal[icl] = new TGraph();
    TString name = Form("gCritVal_%2.0f",confidenceLevel[icl]);
    gCritVal[icl]->SetName(name);
  }
  
  TCanvas *cDeltaChisq[nbins];
  for(int ibin=0;ibin<nbins;ibin++){
    TString name = Form("cDeltaChisq_%d",ibin);
    cDeltaChisq[ibin] = new TCanvas(name,name);
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
  
  TCanvas *cCritVal[nConfidenceLevel];
  for(int icl=0;icl<nConfidenceLevel;icl++){
    TString name = Form("cCritVal_%2.0fCL",confidenceLevel[icl]);
    cCritVal[icl] = new TCanvas(name,name);
  }
    
  // Output logfile

  TString outlogname="";
  outlogname.Append(tag);
  outlogname.Append(".log");
  sprintf(outfilename,"%s",(char*)outlogname.Data());
  outfile.open(outlogname.Data());
  if (!outfile.is_open()) {    
    std::cout << "Error opening file" << endl;
    std::cout << outfilename << std::endl; 
    exit(0);
  }

  //

  for(int ifile = startbin; ifile<=lastbin; ifile++){ // loop over delta_cp files    

    int ibin = ifile - startbin; // delta_cp bins for arrays and histograms
   
    // Get Asimov data set
    if(GetAsimov){
      sigma[ibin]=0.;
      TString infilename = Form(Asimovfilename,ifile); 
      bool AsimovFlag = true;
      TFile *MyFileAsimov = new TFile(infilename,"READ");      
      if(!MyFileAsimov || MyFileAsimov->IsOpen()==false){
    	cout << "        File " << infilename << " not available!!! --> SKIP IT!!!" << endl << endl;	
    	outfile << "        File " << infilename << " not available!!! --> SKIP IT!!!" << endl << endl;	
    	AsimovFlag = false;
    	MyFileAsimov->Close();
    	MyFileAsimov->Delete();
      } 
      TTree *fitChainAsimov;
      if(AsimovFlag){
	fitChainAsimov = (TTree*)MyFileAsimov->Get("fitnt");
	if(!fitChainAsimov){
	  cout << "        Problem in reading the tree --> SKIP IT!!!" << endl << endl;
	  outfile << "        Problem in reading the tree --> SKIP IT!!!" << endl << endl;
	  AsimovFlag = false;
	  MyFileAsimov->Close();
	  MyFileAsimov->Delete();
	}
	fitChainAsimov->SetBranchAddress("asimov_sigma",               &asimov_sigma);
	fitChainAsimov->GetEntry(0);
	sigma[ifile-startbin] = asimov_sigma;
	
	cout << endl;
	cout << "Asimov data set from: " << infilename << endl;
	cout << "sigma = " << sigma[ifile-startbin] << endl;
	cout << endl;
      }
    }

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
    
    mu[ifile-startbin] = input_deltacp; // true value for non-central chi2  

    // 

    for(Long64_t i = 0; i < nentries; i++) {       

      fitChain->GetEntry(i);
      
      //if((migradmin_true == 0 && migraderr_true == 3 && migradmin_min == 0 && migraderr_min == 3) || (migradmin_true == 0 && migraderr_true == 3 && migradmin_min == -1 && migraderr_min == 0)) {
      if((migradmin_true == 0 && migradmin_min == 0) || (migradmin_true == 0 && migradmin_min == -1)) {
  	double dchisq = bestfit_chisq_true - bestfit_chisq_min;

  	if(dchisq<0.){
  	  if(dchisq>-0.0001) dchisq = 0.;
  	  else negdchisq++;
  	  //cout << "dchisq = " << dchisq << endl;
  	}

  	deltaChisq[ibin].push_back(dchisq);
  	hDeltaChisq[ibin]->Fill(dchisq);

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
    cout << "Frac negative dchisq = " << (double)negdchisq / (double)deltaChisq[ibin].size() << endl;
    cout << "deltaChisq size = " << deltaChisq[ibin].size() << endl;
    cout << endl;
    outfile << "Print outputs:" << endl;
    outfile << endl;
    outfile << "# of bad fits = " << nbadfits << endl;
    outfile << "Frac negative dchisq = " << (double)negdchisq / (double)deltaChisq[ibin].size() << endl;
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

      gCritVal[icl]->SetPoint(gCritVal[icl]->GetN(),input_deltacp,value);
    }
    
    MyFile->Close();
    MyFile->Delete();
    
  } // end loop deltacp bins


  // Draw histograms
  
  for(int icl=0;icl<nConfidenceLevel;icl++){
    cCritVal[icl]->cd();
    gCritVal[icl]->Draw("APL");
  }
  
  for(int ibin=0;ibin<nbins;ibin++) {

    if(!hDeltaChisq[ibin]->Integral()) continue;
    
    cDeltaChisq[ibin]->cd();
    //hDeltaChisq[ibin]->Scale(1/hDeltaChisq[ibin]->Integral());
    hDeltaChisq[ibin]->Scale(1./hDeltaChisq[ibin]->Integral()*Ngod_Bins/(Ngod_Max-Ngod_Min));
    hDeltaChisq[ibin]->GetXaxis()->SetTitle("#Delta#chi^{2}");
    hDeltaChisq[ibin]->GetYaxis()->SetTitle("Fraction of events");
    hDeltaChisq[ibin]->Draw("");
    

    // Fit w/ dchisq asymptotic formula (arXiv:1007.1727v2)
  
    if(sigma[ibin]<conv::kASmallNum){
      cout << endl;
      cout << "/////////////////////////////////////" << endl;
      cout << endl;
      cout << "Can't use Asimov informations for file # " << ibin << endl;
      cout << "sigma = 0" << endl;
      cout << endl;
      cout << "/////////////////////////////////////" << endl;
      cout << endl;
    }
    else{
      cout << endl;
      cout << "Get non-central chi2 w/ 1dof:" << endl; 
      cout << "mu = " << mu[ibin] << endl;
      cout << "sigma = " << sigma[ibin] << endl;
      cout << endl;
    }
    
    TF1 *f_chi2_1dof = new TF1("f_chi2_1dof",chi2_1dof,Ngod_Min,Ngod_Max,1);
    f_chi2_1dof->SetLineColor(kRed);
    f_chi2_1dof->SetParameter(0,SetPar_ndf);
    f_chi2_1dof->SetParameter(0,0.5);
    f_chi2_1dof->SetParNames("dof");
    f_chi2_1dof->Draw("same");    
      
    TF1 *f_chi2_1dof_noncentr = new TF1("f_chi2_1dof_noncentr",chi2_1dof_noncentr,Ngod_Min,Ngod_Max,1);
    f_chi2_1dof_noncentr->SetLineColor(kBlue);
    //f_chi2_1dof_noncentr->SetParameter(0,mu[ibin]); // true mu (true deltacp)
    //f_chi2_1dof_noncentr->SetParameter(1,sigma[ibin]); // sigma   
    f_chi2_1dof_noncentr->SetParNames("dof");
    f_chi2_1dof_noncentr->Draw("same");
    hDeltaChisq[ibin]->GetYaxis()->SetRangeUser(0.,1.);
    
    if(CheckConvergence){
      cDLConv_68[ibin]->cd();
      gDLcritConv[kCL68][ibin]->Draw("APL");
      cDLConv_90[ibin]->cd();
      gDLcritConv[kCL90][ibin]->Draw("APL");
      cDLConv_95[ibin]->cd();
      gDLcritConv[kCL95][ibin]->Draw("APL");
    }
      
  }

  // Write output file
  
  TString outname="";
  outname.Append(tag);
  outname.Append(".root");
  sprintf(outfilename,"%s",(char*)outname.Data());
  TFile *InputValues = new TFile(outfilename, "RECREATE");
  for(int ibin=0;ibin<nbins;ibin++){
    hDeltaChisq[ibin]->Write();
  }  

  for(int icl=0;icl<nConfidenceLevel;icl++){
    gCritVal[icl]->Write();
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


Double_t chi2_1dof(Double_t *x,Double_t *par) {
  Double_t fitval = (pow(x[0],par[0]/2.-1.)*(1./exp(x[0]/2.))) / ((pow(2.,par[0]/2.)*TMath::Gamma(par[0]/2.)));
  return fitval;
}

Double_t chi2_1dof_noncentr(Double_t *x,Double_t *par) { // (arXiv:1007.1727v2 eq.19)
  double fact1,fact2,term,term1,term2,chisq=0.;
  double tmu = x[0]; // test statistic 
  double mu = par[0]; // true value
  double sigma = par[1]; // sigma
  fact1 = 1. / sqrt( tmu*2*TMath::Pi() ); 
  term1 = exp(-0.5 * tmu);
  chisq = fact1*term1;
  return chisq;
}

Double_t chi2_1dof_noncentr_bound(Double_t *x,Double_t *par) { // (arXiv:1007.1727v2 eq.19)
  double fact1,fact2,term,term1,term2,chisq=0.;
  double tmu = x[0]; // test statistic 
  double mu = par[0]; // true value
  double sigma = par[1]; // sigma
  // With physical border
  if(tmu<= mu*mu/sigma*sigma){
    fact1 = 1. / sqrt( tmu*2*TMath::Pi() ); 
    term1 = exp(-0.5 * tmu);
    chisq = fact1*term1;
  }
  else{
    fact1 = 0.5 / sqrt( tmu*2*TMath::Pi() ); 
    term1 = exp(-0.5 * tmu);    
    fact2 = 1./(sqrt(2*TMath::Pi()) * 2*mu/sigma);
    term2 = exp(-0.5 * (tmu + mu*mu/(sigma*sigma))*(tmu + mu*mu/(sigma*sigma))/(2*mu/(sigma*sigma)) );
    chisq = fact1*term1 + fact2*term2;
  }
  return chisq;
}
