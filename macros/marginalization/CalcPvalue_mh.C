
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
// Compute p-values using deltaChisq_01, deltaChisq_0BF and # of events
// Store final results in a output root file to be used by other macros
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

#include <vector>

using namespace std;
using namespace valor;
using namespace valor::conv;

void CalcPvalue_mh
(
 TString tag = "prova",
 double deltachisq_beta01_data = 1.43461,
 char infilename[300] = "fit_fitmarg_mh_normal_nexp_2_1.000_-999.000-job0-_0_1.root",
 )
{

  // const double MaxHisto = 20.;
  // const double MinHisto = -20.;
  // const double NBinHisto = 200.;

  const double MaxHisto = 20.;
  const double MinHisto = -20.;
  const double NBinHisto = 200.;
  
  char name[200];

  TFile *filein = new TFile(infilename,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  
  vector<double> vec_chisq_01;
  vector<double> vec_chisq_0BF;
  vector<double> vec_evt_toys;
  vector<double> vec_like_toys;

  PhysParamList *iparm;
  FitParamList *fparm; 
  vector<double> *chisqmin_subrange;
  vector<int> *niterations_subrange;
  vector<int> *nfits_subrange;
  vector<int> *fitstatus_subrange;
  vector<double> *bestfit0_subrange;
  vector<double> *bestfit1_subrange;  
  double truepar0=-999;
  double truepar1=-999;    
  int nfits=-999;
  int nchisqcalls=-999;
  int niterations=-999;
  int fitstatus=-999;
  double chisqmin=conv::kALargeNum;
  double chisqtrue0=conv::kALargeNum;
  double chisqtrue1=conv::kALargeNum;
  double deltachisqtrue=conv::kALargeNum;
  TH2D *hchisq2D_out;
  TH1D *hchisq1D_out;
  TH2D *hstep2D_out;
  TH1D *hstep1D_out;
  double bestfit0=-999;
  double bestfit1=-999;
  double bestfiterrNeg0=-999;
  double bestfiterrNeg1=-999;
  double bestfiterrPos0=-999;
  double bestfiterrPos1=-999;
  double bestfiterrSym0=-999;
  double bestfiterrSym1=-999;  
  double bestfitFullChisq0=-999;
  double bestfitFullChisq1=-999;
  double chisqminFullChisq=conv::kALargeNum;
  int ToyID=-999;
  TH2D *hErecoToy_RHC1Re;

  TBranch        *b_chisqmin_subrange;   //!

  treein->SetBranchAddress("best_fit_params",&fparm);
  treein->SetBranchAddress("input_params",&iparm);
  treein->SetBranchAddress("truepar0", &truepar0 );
  treein->SetBranchAddress("truepar1", &truepar1 );
  treein->SetBranchAddress("bestfit0", &bestfit0 );
  treein->SetBranchAddress("bestfit1", &bestfit1 );
  treein->SetBranchAddress("hchisq2D", &hchisq2D_out );
  treein->SetBranchAddress("hchisq1D", &hchisq1D_out );
  treein->SetBranchAddress("hstep2D", &hstep2D_out );
  treein->SetBranchAddress("hstep1D", &hstep1D_out );
  treein->SetBranchAddress("fitstatus",&fitstatus);
  treein->SetBranchAddress("chisqmin",&chisqmin);
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("deltachisqtrue",&deltachisqtrue);
  treein->SetBranchAddress("nfits",&nfits);
  treein->SetBranchAddress("nchisqcalls",&nchisqcalls);
  treein->SetBranchAddress("niterations",&niterations);
  treein->SetBranchAddress("chisqmin_subrange",&chisqmin_subrange,&b_chisqmin_subrange);
  treein->SetBranchAddress("niterations_subrange",&niterations_subrange);
  treein->SetBranchAddress("nfits_subrange",&nfits_subrange);
  treein->SetBranchAddress("fitstatus_subrange",&fitstatus_subrange);
  treein->SetBranchAddress("bestfit0_subrange",&bestfit0_subrange);
  treein->SetBranchAddress("bestfit1_subrange",&bestfit1_subrange);
  treein->SetBranchAddress("bestfiterrNeg0",&bestfiterrNeg0);
  treein->SetBranchAddress("bestfiterrNeg1",&bestfiterrNeg1);
  treein->SetBranchAddress("bestfiterrPos0",&bestfiterrPos0);
  treein->SetBranchAddress("bestfiterrPos1",&bestfiterrPos1);
  treein->SetBranchAddress("bestfiterrSym0",&bestfiterrSym0);
  treein->SetBranchAddress("bestfiterrSym1",&bestfiterrSym1);
  treein->SetBranchAddress("bestfitFullChisq0", &bestfitFullChisq0 );
  treein->SetBranchAddress("bestfitFullChisq1", &bestfitFullChisq1 );
  treein->SetBranchAddress("chisqminFullChisq",&chisqminFullChisq);
  treein->SetBranchAddress("ToyID",&ToyID);
  treein->SetBranchAddress("hErecoToy_RHC1Re",&hErecoToy_RHC1Re);


  // Declare histograms
  
  TH1D *hnevents = new TH1D("hnevents","hnevents",NBinHisto,0.,MaxHisto);  
  TH1D *hdeltachisqDiff01 = new TH1D("hdeltachisqDiff01","hdeltachisqDiff01",NBinHisto,MinHisto,MaxHisto);
  TH1D *hchisqmin = new TH1D("hchisqmin","hchisqmin",400,0,100);
  TH1D *hchisq0 = new TH1D("hchisq0","hchisq0",400,0,100);
  TH1D *hchisq1 = new TH1D("hchisq1","hchisq1",400,0,100);

  // Check if 1D or 2D fit
  bool IsFit2D;
  treein->GetEntry(0);
  if(bestfit1==-999) IsFit2D = false;
  else               IsFit2D = true;

  // loop over the tree  
  
  int Nentries = treein->GetEntries();
  for(int ientry=0;ientry<Nentries;ientry++){
  
    treein->GetEntry(ientry);

    if(!(ientry%500)) cout << "Entry " << ientry << endl;

    double nevents = (double)hErecoToy_RHC1Re->Integral();    
    if(nevents<0){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events = " << nevents << " < 0" << endl;
      cerr << endl;
      exit(1);
    }

    hnevents->Fill(nevents);    
    cout << "chisqtrue0 = " << chisqtrue0 
	 << " _ chisqtrue1 = " << chisqtrue1
	 << " _ chisqmin = " << chisqmin
	 << endl;
    double deltachisq_01 = chisqtrue0 - chisqtrue1;
    hdeltachisqDiff01->Fill(deltachisq_01);

    // Minimized chisq distribution
    hchisqmin->Fill(chisqmin);
    hchisq0->Fill(chisqtrue0);
    hchisq1->Fill(chisqtrue1);
          
    // Fill vectors for pvalue calculations
    vec_evt_toys.push_back(nevents);
    vec_chisq_01.push_back(deltachisq_01);
       
  } // end loop over toys


  if(deltachisq_01<-5 || deltachisq_01>20){
    cout << "deltachisq_01 = " << deltachisq_01 << endl;
    exit(1);
  }


  
  // Compute p-values

  double Ntoys = (double)Nentries;
  
  double pvalue_01_right = GetPvalue_rightside(vec_chisq_01,Ntoys,deltachisq_beta01_data);
  double pvalue_01_left = GetPvalue_leftside(vec_chisq_01,Ntoys,deltachisq_beta01_data);

  // Draw histograms
  
  // treein->GetEntry(493);
  // //treein->GetEntry(4739);
  // TCanvas *cchisq_out_GetEntry = new TCanvas("cchisq_out_GetEntry","cchisq_out_GetEntry");
  // if(IsFit2D) hchisq2D_out->DrawCopy("colz");
  // else        hchisq1D_out->DrawCopy();
  // TCanvas *cErecoToy_GetEntry = new TCanvas("cErecoToy_GetEntry","cErecoToy_GetEntry");
  // hErecoToy_RHC1Re->DrawCopy();

  TCanvas *cchisqmin = new TCanvas("cchisqmin","cchisqmin");
  hchisqmin->GetYaxis()->SetTitle("Number of entries");
  hchisqmin->GetXaxis()->SetTitle("#chi^{2} (#bar{#beta}_{BF})");
  hchisqmin->DrawCopy();
  hchisqmin->Draw("same");
  hchisqmin->SetLogy();

  TCanvas *cchisq0 = new TCanvas("cchisq0","cchisq0");
  hchisq0->GetYaxis()->SetTitle("Number of entries");
  hchisq0->GetXaxis()->SetTitle("#chi^{2} (#bar{#beta=0})");
  hchisq0->DrawCopy();
  hchisq0->Draw("same");
  hchisq0->SetLogy();

  TCanvas *cchisq1 = new TCanvas("cchisq1","cchisq1");
  hchisq1->GetYaxis()->SetTitle("Number of entries");
  hchisq1->GetXaxis()->SetTitle("#chi^{2} (#bar{#beta=1})");
  hchisq1->DrawCopy();
  hchisq1->Draw("same");
  hchisq1->SetLogy();

  TCanvas *cPvalue_01 = new TCanvas("cPvalue_01","cPvalue_01");
  cPvalue_01->SetLogy();
  hdeltachisqDiff01->GetYaxis()->SetTitle("Number of entries");
  hdeltachisqDiff01->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)");
  hdeltachisqDiff01->DrawCopy();
  TLine *ldeltachisqDiff01 = new TLine(deltachisq_beta01_data,0.,deltachisq_beta01_data,hdeltachisqDiff01->GetMaximum());
  ldeltachisqDiff01->SetLineColor(kRed);
  ldeltachisqDiff01->SetLineWidth(2);
  ldeltachisqDiff01->Draw();
  TLegend *lPvalue_01 = new TLegend(0.5,0.6,0.9,0.9);
  lPvalue_01->SetFillColor(0);
  lPvalue_01->SetLineColor(0);
  lPvalue_01->SetTextFont(132);
  lPvalue_01->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_beta01_data); 
  lPvalue_01->AddEntry(ldeltachisqDiff01,name,"l");
  lPvalue_01->Draw();


  char outname[200];
  sprintf(outname,"%s.root",tag.Data());
  TFile *outfile = new TFile(outname,"RECREATE");
  hNdata->Write();
  hdeltachisqDiff01->Write("hdeltachisqDiff01");
  outfile->Close();

  // Print results
  
  cout << endl;
  cout << std::setprecision(10);
  cout << "truepar0 = " << truepar0 << endl;
  if(IsFit2D) cout << "truepar1 = " << truepar1 << endl;
  cout << endl;
  cout << "Calculate right sided p-values:" << endl;
  cout << "   p-value 01 (rate+shape) = " << pvalue_01_right << endl;
  cout << endl;
  cout << "Calculate left sided p-values:" << endl;
  cout << "   p-value 01 (rate+shape) = " << pvalue_01_left << endl;
  cout << endl;

  return;
}



double GetPvalue_rightside(vector<double> vec_toys,int ntoys,double data){
  double ntot_toys = ntoys;
  if(ntot_toys != vec_toys.size()){    
    cout << endl;
    cout << "Warning:" << endl;
    cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
    cout << "Input total number of toys = " << ntot_toys << endl; 
    ntot_toys = vec_toys.size();
    cout << "Set total number of toys to vector size: " << ntot_toys << endl;
    cout << endl;
  }
  if(ntot_toys<=0.){
    cerr << endl;
    cerr << "ntot_toys = " << ntot_toys << endl;
    cerr << endl;
    exit(1);
  }

  cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
  cout << "Input total number of toys = " << ntot_toys << endl; 
  ntot_toys = vec_toys.size();
  cout << "Set total number of toys to vector size: " << ntot_toys << endl;
  cout << endl;

  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it >= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}

double GetPvalue_leftside(vector<double> vec_toys,int ntoys,double data){
  double ntot_toys = ntoys;
  if(ntot_toys != vec_toys.size()){    
    cout << endl;
    cout << "Warning:" << endl;
    cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
    cout << "Input total number of toys = " << ntot_toys << endl; 
    ntot_toys = vec_toys.size();
    cout << "Set total number of toys to vector size: " << ntot_toys << endl;
    cout << endl;
  }
  if(ntot_toys<=0.){
    cerr << endl;
    cerr << "ntot_toys = " << ntot_toys << endl;
    cerr << endl;
    exit(1);
  }

  cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
  cout << "Input total number of toys = " << ntot_toys << endl; 
  ntot_toys = vec_toys.size();
  cout << "Set total number of toys to vector size: " << ntot_toys << endl;
  cout << endl;

  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it <= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}

double GetPvalue_Discrete(TH1D *histo,double data){

  int bin_data = histo->FindBin(data);
  double pvalue = histo->Integral(bin_data,histo->GetNbinsX());
  pvalue /= histo->Integral();
  
  return pvalue;
}
