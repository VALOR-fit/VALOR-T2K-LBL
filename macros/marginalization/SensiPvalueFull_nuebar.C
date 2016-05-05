
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

using namespace std;
using namespace valor;
using namespace valor::conv;

void SensiPvalueFull_nuebar
(
 // Run 16 - Asimov A
 const char * Npot = "4.011x10^20 POT",
 double expbkg = 1.718,
 double data = 3.218, // # of data events
 //double deltachisq_beta01_data = 1.41638, // data deltachisq 01
 double deltachisq_beta01_data = 1.4346108, // data deltachisq 01 
 double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 char infilename_beta0[200] = "/mnt/lnec/davides/VALOR/results/fit/BANFF_FITBIAS_14_3_16/NUEBAR/100k/fit_fitmarg_joint_betaBar_0_100kMargToys_react_100kToys.root",
 char infilename_beta1[200] = "/mnt/lnec/davides/VALOR/results/fit/BANFF_FITBIAS_14_3_16/NUEBAR/100k/fit_fitmarg_joint_betaBar_1_100kMargToys_react_100kToys.root",
  )
{
  char name[200];

  double Ntoys_beta0;
  double Ntoys_beta1;

  vector<double> vec_pvalue_rateonly;
  vector<double> vec_pvalue_rateshape;
  
  vector<double> vec_chisq_01_beta0;
  vector<double> vec_evt_toys_beta0;

  vector<double> vec_chisq_01_beta1;
  vector<double> vec_evt_toys_beta1;

  double pvalue_rateonly_mean_all=0.;
  double pvalue_rateshape_mean_all=0.;
 
  double chisqtrue0;
  double chisqtrue1;
  TH2D *hErecoToy_RHC1Re;
  int ToyID;

  double Frac_Pvalue_rateonly_Bkg = 0.;
  double Frac_Pvalue_rateshape_Bkg = 0.;
  double Frac_Pvalue_rateonly_0p05 = 0.;
  double Frac_Pvalue_rateshape_0p05 = 0.;
  double Frac_Pvalue_rateonly_0p10 = 0.;
  double Frac_Pvalue_rateshape_0p10 = 0.;

  double Frac_Pvalue_rateonly_higher_rateshape = 0.;
  double Frac_Pvalue_rateonly_higher_rateshape_Bkg = 0.;

  // Declare histograms

  TH1D *hPvalue_rateonly = new TH1D("hPvalue_rateonly","hPvalue_rateonly",110,0.,1.1); 
  TH1D *hPvalue_rateshape = new TH1D("hPvalue_rateshape","hPvalue_rateshape",110,0.,1.1); 

  TH1D *hPvalue_rateonly_Bkg = new TH1D("hPvalue_rateonly_Bkg","hPvalue_rateonly_Bkg",110,0.,1.1); 
  TH1D *hPvalue_rateshape_Bkg = new TH1D("hPvalue_rateshape_Bkg","hPvalue_rateshape_Bkg",110,0.,1.1); 

  TH1D *hPvalue_rateonly_3evt = new TH1D("hPvalue_rateonly_3evt","hPvalue_rateonly_3evt",110,0.,1.1); 
  TH1D *hPvalue_rateshape_3evt = new TH1D("hPvalue_rateshape_3evt","hPvalue_rateshape_3evt",110,0.,1.1); 

  TH1D *hPvalue_rateonly_4evt = new TH1D("hPvalue_rateonly_4evt","hPvalue_rateonly_4evt",110,0.,1.1); 
  TH1D *hPvalue_rateshape_4evt = new TH1D("hPvalue_rateshape_4evt","hPvalue_rateshape_4evt",110,0.,1.1); 

  TH1D *hPvalue_rateonly_5evt = new TH1D("hPvalue_rateonly_5evt","hPvalue_rateonly_5evt",110,0.,1.1); 
  TH1D *hPvalue_rateshape_5evt = new TH1D("hPvalue_rateshape_5evt","hPvalue_rateshape_5evt",110,0.,1.1); 

  TH1D *hPvalue_rateonly_6evt = new TH1D("hPvalue_rateonly_6evt","hPvalue_rateonly_6evt",110,0.,1.1); 
  TH1D *hPvalue_rateshape_6evt = new TH1D("hPvalue_rateshape_6evt","hPvalue_rateshape_6evt",110,0.,1.1); 

  TH2D *hPvalueVsNevt = new TH2D("hPvalueVsNevt","hPvalueVsNevt",110,0.,1.1,12,0,12);

  TH1D *hPvalue_rateonly_0p05 = new TH1D("hPvalue_rateonly_0p05","hPvalue_rateonly_0p05",110,0.,0.06); 
  TH1D *hPvalue_rateonly_0p10 = new TH1D("hPvalue_rateonly_0p10","hPvalue_rateonly_0p10",110,0.,0.11);   
  TH1D *hPvalue_rateshape_0p05 = new TH1D("hPvalue_rateshape_0p05","hPvalue_rateshape_0p05",110,0.,0.06); 
  TH1D *hPvalue_rateshape_0p10 = new TH1D("hPvalue_rateshape_0p10","hPvalue_rateshape_0p10",110,0.,0.11); 

  
  // Get first file (beta=0)

  TFile *filein = new TFile(infilename_beta0,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("hErecoToy_RHC1Re",&hErecoToy_RHC1Re);
  treein->SetBranchAddress("ToyID",&ToyID);

  int Nentries = treein->GetEntries();

  Ntoys_beta0 = (double) Nentries; // Get # of entries for p-value calculation

  for(int ientry=0;ientry<Nentries;ientry++){
    treein->GetEntry(ientry);
    double nevents = (double)hErecoToy_RHC1Re->Integral();    
    if(nevents<0){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events = " << nevents << " < 0" << endl;
      cerr << endl;
      exit(1);
    }
    
    //if(!(ientry%1000)) cout << "Entry " << ientry << endl;

    double deltachisq_01 = chisqtrue0 - chisqtrue1;

    vec_evt_toys_beta0.push_back(nevents);
    vec_chisq_01_beta0.push_back(deltachisq_01);
  }


  // Get first file (beta=1)

  filein = new TFile(infilename_beta1,"READ");
  treein = (TTree*)filein->Get("fitnt");
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("hErecoToy_RHC1Re",&hErecoToy_RHC1Re);
  treein->SetBranchAddress("ToyID",&ToyID);

  Nentries = treein->GetEntries();
  
  Ntoys_beta1 = (double) Nentries; // Get # of entries for p-value calculation

  for(int ientry=0;ientry<Nentries;ientry++){
    treein->GetEntry(ientry);
    if( !(ientry % 1000) ) cout << "Compute pvalues for toy " << ientry << endl;
    
    double nevents = (double)hErecoToy_RHC1Re->Integral();    
    if(nevents<0.){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events = " << nevents << " < 0" << endl;
      cerr << endl;
      exit(1);
    }

    double deltachisq_01 = chisqtrue0 - chisqtrue1;
    vec_evt_toys_beta1.push_back(nevents);
    vec_chisq_01_beta1.push_back(deltachisq_01);
    
    double pvalue_rateonly_itoy = GetPvalue(vec_evt_toys_beta0,Ntoys_beta0,nevents);
    if(pvalue_rateonly_itoy<0. || pvalue_rateonly_itoy>1.){
      cerr << endl;
      cerr << "ToyID = " << ToyID << endl;
      cerr << "pvalue_rateonly = " << pvalue_rateonly_itoy << endl;
      cerr << endl;
      exit(1);
    }

    hPvalue_rateonly->Fill(pvalue_rateonly_itoy);
    if(nevents==3) hPvalue_rateonly_3evt->Fill(pvalue_rateonly_itoy);
    else if(nevents==4) hPvalue_rateonly_4evt->Fill(pvalue_rateonly_itoy);
    else if(nevents==5) hPvalue_rateonly_5evt->Fill(pvalue_rateonly_itoy);
    else if(nevents==6) hPvalue_rateonly_6evt->Fill(pvalue_rateonly_itoy);   

    double pvalue_rateshape_itoy = GetPvalue(vec_chisq_01_beta0,Ntoys_beta0,deltachisq_01);
    if(pvalue_rateshape_itoy<0. || pvalue_rateshape_itoy>1.){
      cerr << endl;
      cerr << "ToyID = " << ToyID << endl;
      cerr << "pvalue_rateshape = " << pvalue_rateshape_itoy << endl;
      cerr << endl;
      exit(1);
    }    
    
    // if(pvalue_rateshape_itoy<0.0001){
    //   cout << "Toy " << ientry << endl;
    //   cout << "nevents = " << nevents << endl;
    //   cout << "deltachisq_01 = " << deltachisq_01 << endl;
    //   cout << "pvalue_rateonly = " << pvalue_rateonly_itoy << endl;
    //   cout << "pvalue_rateshape = " << pvalue_rateshape_itoy << endl;    
    //   cout << endl;
    // }
    // else continue;
     

    hPvalue_rateshape->Fill(pvalue_rateshape_itoy);
    
    if(pvalue_rateshape_itoy>0.02)
      hPvalueVsNevt->Fill(pvalue_rateshape_itoy,nevents);
    else 
      hPvalueVsNevt->Fill(0.02,nevents); // for logarithmic scale

    // Isolate pvalue distributions with different number of events

    if(nevents==3) hPvalue_rateshape_3evt->Fill(pvalue_rateshape_itoy);
    else if(nevents==4) hPvalue_rateshape_4evt->Fill(pvalue_rateshape_itoy);
    else if(nevents==5) hPvalue_rateshape_5evt->Fill(pvalue_rateshape_itoy);
    else if(nevents==6) hPvalue_rateshape_6evt->Fill(pvalue_rateshape_itoy);
   
    //

    pvalue_rateonly_mean_all += pvalue_rateonly_itoy / (double)Nentries;
    pvalue_rateshape_mean_all += pvalue_rateshape_itoy / (double)Nentries;

    // Isolate pvalue distributions with pvalue < 0.05 (0.10)
    if(pvalue_rateonly_itoy<=0.05){
      hPvalue_rateonly_0p05->Fill(pvalue_rateonly_itoy);
    }
    if(pvalue_rateonly_itoy<=0.10){
      hPvalue_rateonly_0p10->Fill(pvalue_rateonly_itoy);
      //if(pvalue_rateonly_itoy>0.05) cout << pvalue_rateonly_itoy << endl;
    }
    if(pvalue_rateshape_itoy<=0.05){
      hPvalue_rateshape_0p05->Fill(pvalue_rateshape_itoy);
    }
    if(pvalue_rateshape_itoy<=0.10){
      hPvalue_rateshape_0p10->Fill(pvalue_rateshape_itoy);
    }

    // Isolate p-value for # obs events > # expected bkg
    if(nevents > expbkg){
      hPvalue_rateonly_Bkg->Fill(pvalue_rateonly_itoy);
      hPvalue_rateshape_Bkg->Fill(pvalue_rateshape_itoy);
    }
    // else{
    //   cout << "nevents = " << nevents << endl << endl;
    //   cout << "expbkg = " << expbkg << endl;
    //   cout << endl;
    // }

    // Check # of toys with rate-only pvalue > rate+shape pvalue
    
    if(pvalue_rateonly_itoy > pvalue_rateshape_itoy){
      Frac_Pvalue_rateonly_higher_rateshape += 1. / (double)Nentries;
      if(nevents > expbkg) Frac_Pvalue_rateonly_higher_rateshape_Bkg += 1. / Nentries;
    }

  } // end loop


  // Compute p-values
  
  double pvalue_evt_beta0 = GetPvalue(vec_evt_toys_beta0,Ntoys_beta0,data);
  double pvalue_01_beta0 = GetPvalue(vec_chisq_01_beta0,Ntoys_beta0,deltachisq_beta01_data);

  double pvalue_evt_beta1 = GetPvalue(vec_evt_toys_beta1,Ntoys_beta1,data);
  double pvalue_01_beta1 = GetPvalue(vec_chisq_01_beta1,Ntoys_beta1,deltachisq_beta01_data);
  
  // Draw histograms
  
  TCanvas *cPvalue_rateonly = new TCanvas("cPvalue_rateonly","cPvalue_rateonly");
  cPvalue_rateonly->SetLogy();
  hPvalue_rateonly_3evt->SetFillColor(kBlue);
  hPvalue_rateonly_3evt->SetFillStyle(3305);
  hPvalue_rateonly_4evt->SetFillColor(kRed);
  hPvalue_rateonly_4evt->SetFillStyle(3345);
  hPvalue_rateonly_5evt->SetFillColor(kGreen);
  hPvalue_rateonly_5evt->SetFillStyle(3354);
  hPvalue_rateonly_6evt->SetFillColor(kCyan);
  hPvalue_rateonly_6evt->SetFillStyle(3395);
  hPvalue_rateonly->GetXaxis()->SetTitle("p-value");
  hPvalue_rateonly->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateonly->DrawCopy();
  // hPvalue_rateonly_3evt->DrawCopy("same");
  // hPvalue_rateonly_4evt->DrawCopy("same");
  // hPvalue_rateonly_5evt->DrawCopy("same");
  // hPvalue_rateonly_6evt->DrawCopy("same");
  TLine *linePvalue_rateonly_Asimov = new TLine(pvalue_evt_beta0,0.,pvalue_evt_beta0,hPvalue_rateonly->GetMaximum());
  linePvalue_rateonly_Asimov->SetLineWidth(2);
  linePvalue_rateonly_Asimov->SetLineColor(kRed);
  linePvalue_rateonly_Asimov->DrawClone();
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateonly->GetMean(),hPvalue_rateonly->GetRMS()); 
  TLatex *lPvalue_rateonly = new TLatex(0.70,hPvalue_rateonly->GetMaximum()*0.8,name);
  lPvalue_rateonly->SetTextSize(0.07);
  lPvalue_rateonly->Draw();
  
  //
  
  TCanvas *cPvalue_rateshape = new TCanvas("cPvalue_rateshape","cPvalue_rateshape");
  cPvalue_rateshape->SetLogy();
  hPvalue_rateshape_3evt->SetFillColor(kBlue);
  hPvalue_rateshape_3evt->SetFillStyle(3305);
  hPvalue_rateshape_4evt->SetFillColor(kRed);
  hPvalue_rateshape_4evt->SetFillStyle(3345);
  hPvalue_rateshape_5evt->SetFillColor(kGreen);
  hPvalue_rateshape_5evt->SetFillStyle(3354);
  hPvalue_rateshape_6evt->SetFillColor(kCyan);
  hPvalue_rateshape_6evt->SetFillStyle(3395);
  hPvalue_rateshape->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape->DrawCopy();
  // hPvalue_rateshape_3evt->DrawCopy("same");
  // hPvalue_rateshape_4evt->DrawCopy("same");
  // hPvalue_rateshape_5evt->DrawCopy("same");
  // hPvalue_rateshape_6evt->DrawCopy("same");
  TLine *linePvalue_rateshape_Asimov = new TLine(pvalue_01_beta0,0.,pvalue_01_beta0,hPvalue_rateshape->GetMaximum());
  linePvalue_rateshape_Asimov->SetLineWidth(2);
  linePvalue_rateshape_Asimov->SetLineColor(kRed);
  linePvalue_rateshape_Asimov->DrawClone();
  TLegend *legPvalue_rateshape = new TLegend(0.7,0.6,1.0,0.9);
  legPvalue_rateshape->SetFillColor(0);
  legPvalue_rateshape->SetLineColor(0);
  legPvalue_rateshape->SetTextFont(132);
  legPvalue_rateshape->SetTextSize(0.06);
  legPvalue_rateshape->AddEntry(hPvalue_rateshape,"Toy experiments","l");
  legPvalue_rateshape->AddEntry(linePvalue_rateshape_Asimov,"Asimov dataset","l");
  legPvalue_rateshape->Draw();
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape->GetMean(),hPvalue_rateshape->GetRMS()); 
  TLatex *lPvalue_rateshape = new TLatex(0.70,hPvalue_rateshape->GetMaximum()*0.8,name);
  lPvalue_rateshape->SetTextSize(0.06);
  lPvalue_rateshape->Draw();

  //

  TCanvas *cPvalue_rateshape_3evt = new TCanvas("cPvalue_rateshape_3evt","cPvalue_rateshape_3evt");
  cPvalue_rateshape_3evt->SetLogy();
  hPvalue_rateshape_3evt->SetFillStyle(0);
  hPvalue_rateshape_3evt->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape_3evt->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_3evt->DrawCopy();  
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_3evt->GetMean(),hPvalue_rateshape_3evt->GetRMS()); 
  TLatex *lPvalue_rateshape_3evt = new TLatex(0.70,hPvalue_rateshape_3evt->GetMaximum()*0.8,name);
  lPvalue_rateshape_3evt->SetTextSize(0.07);
  lPvalue_rateshape_3evt->Draw();

  //

  TCanvas *cPvalue_rateshape_4evt = new TCanvas("cPvalue_rateshape_4evt","cPvalue_rateshape_4evt");
  cPvalue_rateshape_4evt->SetLogy();
  hPvalue_rateshape_4evt->SetFillStyle(0);
  hPvalue_rateshape_4evt->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape_4evt->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_4evt->DrawCopy();  
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_4evt->GetMean(),hPvalue_rateshape_4evt->GetRMS()); 
  TLatex *lPvalue_rateshape_4evt = new TLatex(0.70,hPvalue_rateshape_4evt->GetMaximum()*0.8,name);
  lPvalue_rateshape_4evt->SetTextSize(0.07);
  lPvalue_rateshape_4evt->Draw();

  //

  TCanvas *cPvalue_rateshape_5evt = new TCanvas("cPvalue_rateshape_5evt","cPvalue_rateshape_5evt");
  cPvalue_rateshape_5evt->SetLogy();
  hPvalue_rateshape_5evt->SetFillStyle(0);
  hPvalue_rateshape_5evt->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape_5evt->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_5evt->DrawCopy();  
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_5evt->GetMean(),hPvalue_rateshape_5evt->GetRMS()); 
  TLatex *lPvalue_rateshape_5evt = new TLatex(0.70,hPvalue_rateshape_5evt->GetMaximum()*0.8,name);
  lPvalue_rateshape_5evt->SetTextSize(0.07);
  lPvalue_rateshape_5evt->Draw();

  //

  TCanvas *cPvalue_rateshape_6evt = new TCanvas("cPvalue_rateshape_6evt","cPvalue_rateshape_6evt");
  cPvalue_rateshape_6evt->SetLogy();
  hPvalue_rateshape_6evt->SetFillStyle(0);
  hPvalue_rateshape_6evt->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape_6evt->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_6evt->DrawCopy();  
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_6evt->GetMean(),hPvalue_rateshape_6evt->GetRMS()); 
  TLatex *lPvalue_rateshape_6evt = new TLatex(0.70,hPvalue_rateshape_6evt->GetMaximum()*0.8,name);
  lPvalue_rateshape_6evt->SetTextSize(0.07);
  lPvalue_rateshape_6evt->Draw();

  //

  TCanvas *cPvalueVsNevt = new TCanvas("hPvalueVsNevt","hPvalueVsNevt");
  hPvalueVsNevt->GetXaxis()->SetTitle("p-value");
  hPvalueVsNevt->GetYaxis()->SetTitle("# of events");
  hPvalueVsNevt->Draw("colz");

  // 

  TCanvas *cPvalue_rateonly_0p10 = new TCanvas("cPvalue_rateonly_0p10","cPvalue_rateonly_0p10");
  cPvalue_rateonly_0p10->SetLogy();
  hPvalue_rateonly_0p10->GetXaxis()->SetTitle("p-value (#leq 0.10)");
  hPvalue_rateonly_0p10->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateonly_0p10->Draw();
  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateonly_0p10->GetMean(),hPvalue_rateonly_0p10->GetRMS()); 
  TLatex *lPvalue_rateonly_0p10 = new TLatex(0.07,hPvalue_rateonly_0p10->GetMaximum()*0.8,name);
  lPvalue_rateonly_0p10->SetTextSize(0.07);
  lPvalue_rateonly_0p10->Draw();

  TCanvas *cPvalue_rateshape_0p10 = new TCanvas("cPvalue_rateshape_0p10","cPvalue_rateshape_0p10");
  cPvalue_rateshape_0p10->SetLogy();
  hPvalue_rateshape_0p10->GetXaxis()->SetTitle("p-value (#leq 0.10)");
  hPvalue_rateshape_0p10->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_0p10->Draw();
  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateshape_0p10->GetMean(),hPvalue_rateshape_0p10->GetRMS()); 
  TLatex *lPvalue_rateshape_0p10 = new TLatex(0.07,hPvalue_rateshape_0p10->GetMaximum()*0.8,name);
  lPvalue_rateshape_0p10->SetTextSize(0.07);
  lPvalue_rateshape_0p10->Draw();

  TCanvas *cPvalue_rateonly_0p05 = new TCanvas("cPvalue_rateonly_0p05","cPvalue_rateonly_0p05");
  cPvalue_rateonly_0p05->SetLogy();
  hPvalue_rateonly_0p05->GetXaxis()->SetTitle("p-value (#leq 0.05)");
  hPvalue_rateonly_0p05->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateonly_0p05->Draw();
  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateonly_0p05->GetMean(),hPvalue_rateonly_0p05->GetRMS()); 
  TLatex *lPvalue_rateonly_0p05 = new TLatex(0.03,hPvalue_rateonly_0p05->GetMaximum()*0.8,name);
  lPvalue_rateonly_0p05->SetTextSize(0.07);
  lPvalue_rateonly_0p05->Draw();

  TCanvas *cPvalue_rateshape_0p05 = new TCanvas("cPvalue_rateshape_0p05","cPvalue_rateshape_0p05");
  cPvalue_rateshape_0p05->SetLogy();
  hPvalue_rateshape_0p05->GetXaxis()->SetTitle("p-value (#leq 0.05)");
  hPvalue_rateshape_0p05->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_0p05->Draw();
  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateshape_0p05->GetMean(),hPvalue_rateshape_0p05->GetRMS()); 
  TLatex *lPvalue_rateshape_0p05 = new TLatex(0.03,hPvalue_rateshape_0p05->GetMaximum()*0.8,name);
  lPvalue_rateshape_0p05->SetTextSize(0.07);
  lPvalue_rateshape_0p05->Draw();





  // p-value for # obs events > # exp bkg
  
  TCanvas *cPvalue_rateonly_Bkg = new TCanvas("cPvalue_rateonly_Bkg","cPvalue_rateonly_Bkg");
  cPvalue_rateonly_Bkg->SetLogy();
  hPvalue_rateonly_Bkg->GetXaxis()->SetTitle("p-value (# obs. events > # exp bkg)");
  hPvalue_rateonly_Bkg->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateonly_Bkg->Draw();
  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateonly_Bkg->GetMean(),hPvalue_rateonly_Bkg->GetRMS()); 
  hPvalue_rateonly_Bkg->Draw();
  TLatex *lPvalue_rateonly_Bkg = new TLatex(0.7,hPvalue_rateonly_Bkg->GetMaximum()*0.8,name);
  lPvalue_rateonly_Bkg->SetTextSize(0.07);
  lPvalue_rateonly_Bkg->Draw();
  sprintf(name,"Exp Bkg = %1.3f",expbkg); 
  TLatex *lExpBkg = new TLatex(0.4,hPvalue_rateonly_Bkg->GetMaximum()*0.8,name);
  lExpBkg->SetTextSize(0.07);
  lExpBkg->Draw();


  TCanvas *cPvalue_rateshape_Bkg = new TCanvas("cPvalue_rateshape_Bkg","cPvalue_rateshape_Bkg");
  cPvalue_rateshape_Bkg->SetLogy();
  hPvalue_rateshape_Bkg->GetXaxis()->SetTitle("p-value (# obs. events > # exp bkg)");
  hPvalue_rateshape_Bkg->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_Bkg->Draw();
  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateshape_Bkg->GetMean(),hPvalue_rateshape_Bkg->GetRMS()); 
  hPvalue_rateshape_Bkg->Draw();
  TLatex *lPvalue_rateshape_Bkg = new TLatex(0.7,hPvalue_rateshape_Bkg->GetMaximum()*0.8,name);
  lPvalue_rateshape_Bkg->SetTextSize(0.07);
  lPvalue_rateshape_Bkg->Draw();
  sprintf(name,"Exp Bkg = %1.3f",expbkg); 
  TLatex *lExpBkg = new TLatex(0.4,hPvalue_rateshape_Bkg->GetMaximum()*0.8,name);
  lExpBkg->SetTextSize(0.07);
  lExpBkg->Draw();
  

  // Print results
  
  cout << std::setprecision(10);
  cout << endl;
  cout << "Total # of POT = " << Npot << endl;
  cout << endl;
  cout << "//////////////// Asimov dataset //////////////" << endl;
  cout << "Calculate p-values (beta=0):" << endl;
  cout << "p-value rate-only = " << pvalue_evt_beta0 << endl;
  cout << "p-value rate+shape = " << pvalue_01_beta0 << endl;
  cout << "Calculate p-values (beta=1):" << endl;
  cout << "p-value rate-only = " << pvalue_evt_beta1 << endl;
  cout << "p-value rate+shape = " << pvalue_01_beta1 << endl;
  cout << endl;
  cout << "//////////////// Toys p-value (all toys) //////////////" << endl;
  cout << endl;
  cout << "p-value rate-only:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateonly->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateonly->GetRMS() << endl;
  cout << "Mean = " << pvalue_rateonly_mean_all << endl;
  cout << "p-value rate+shape:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateshape->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateshape->GetRMS() << endl;
  cout << "Mean = " << pvalue_rateshape_mean_all << endl;
  cout << endl;
  cout << "Fraction of toys with rate-only p-value > rate+shape p-value = " << Frac_Pvalue_rateonly_higher_rateshape << endl;
  cout << endl;
  cout << "//////////////// Toys p-value (<0.10) //////////////" << endl;
  cout << endl;
  cout << "Fraction of toys (rate-only) = " 
       << hPvalue_rateonly_0p10->Integral() / hPvalue_rateonly->Integral() << endl;
  cout << "Fraction of toys (rate+shape) = " 
       << hPvalue_rateshape_0p10->Integral() / hPvalue_rateshape->Integral() << endl;
  cout << "p-value rate-only:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateonly_0p10->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateonly_0p10->GetRMS() << endl;
  cout << "p-value rate+shape:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateshape_0p10->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateshape_0p10->GetRMS() << endl;
  cout << endl;
  cout << "//////////////// Toys p-value (<0.05) //////////////" << endl;
  cout << endl;
  cout << "Fraction of toys (rate-only) = " 
       << hPvalue_rateonly_0p05->Integral() / hPvalue_rateonly->Integral() << endl;
  cout << "Fraction of toys (rate+shape) = " 
       << hPvalue_rateshape_0p05->Integral() / hPvalue_rateshape->Integral() << endl;
  cout << "p-value rate-only:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateonly_0p05->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateonly_0p05->GetRMS() << endl;
  cout << "p-value rate+shape:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateshape_0p05->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateshape_0p05->GetRMS() << endl;
  cout << endl;
  cout << "//////////////// Toys p-value (# obs evts > # exp bkg) //////////////" << endl;
  cout << endl;
  cout << "Fraction of toys (rate-only) = " 
       << hPvalue_rateonly_Bkg->Integral() / hPvalue_rateonly->Integral() << endl;
  cout << "Fraction of toys (rate+shape) = " 
       << hPvalue_rateshape_Bkg->Integral() / hPvalue_rateshape->Integral() << endl;
  cout << "p-value rate-only:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateonly_Bkg->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateonly_Bkg->GetRMS() << endl;
  cout << "p-value rate+shape:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateshape_Bkg->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateshape_Bkg->GetRMS() << endl;
  cout << endl;
  cout << "Fraction of toys with rate-only p-value > rate+shape p-value = " << Frac_Pvalue_rateonly_higher_rateshape_Bkg << endl;
  cout << endl;
  
  // Print plots

  cPvalue_rateonly->Print("cPvalue_rateonly.pdf");
  cPvalue_rateshape->Print("cPvalue_rateshape.pdf");
  cPvalue_rateshape_3evt->Print("cPvalue_rateshape_3evt.pdf");
  cPvalue_rateshape_4evt->Print("cPvalue_rateshape_4evt.pdf");
  cPvalue_rateshape_5evt->Print("cPvalue_rateshape_5evt.pdf");
  cPvalue_rateshape_6evt->Print("cPvalue_rateshape_6evt.pdf");
  cPvalueVsNevt->Print("cPvalueVsNevt.pdf");
  cPvalue_rateonly_0p10->Print("cPvalue_rateonly_0p10.pdf");
  cPvalue_rateshape_0p10->Print("cPvalue_rateshape_0p10.pdf");
  cPvalue_rateonly_0p05->Print("cPvalue_rateonly_0p05.pdf");
  cPvalue_rateshape_0p05->Print("cPvalue_rateshape_0p05.pdf");
  cPvalue_rateonly_Bkg->Print("cPvalue_rateonly_Bkg.pdf");
  cPvalue_rateshape_Bkg->Print("cPvalue_rateshape_Bkg.pdf");

  return;
}



double GetPvalue(vector<double> vec_toys,int ntoys,double data){
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
  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it >= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}
