
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
 double expbkg = 1.718,
 double data = 3.218, // # of data events
 char infilename_beta0[200] = "/mnt/lnec/davides/VALOR/results/fit/BANFF_FITBIAS_14_3_16/NUEBAR/100k/fit_fitmarg_joint_betaBar_0_100kMargToys_react_100kToys.root",
 char infilename_beta1[200] = "/mnt/lnec/davides/VALOR/results/fit/BANFF_FITBIAS_14_3_16/NUEBAR/100k/fit_fitmarg_joint_betaBar_1_100kMargToys_react_100kToys.root"
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
  TH2D *hErecoExp_RHC1Re;
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

  TH1D *hPValueBeta0 = new TH1D("hPValueBeta0","hPValueBeta0",200,0,1.5);
  TH1D *hPValueBeta0_toyevt0 = new TH1D("hPValueBeta0_toyevt0","hPValueBeta0_toyevt0",200,0,1.5);
  TH1D *hPValueBeta0_toyevt1 = new TH1D("hPValueBeta0_toyevt1","hPValueBeta0_toyevt1",200,0,1.5);
  TH1D *hPValueBeta0_toyevt2 = new TH1D("hPValueBeta0_toyevt2","hPValueBeta0_toyevt2",200,0,1.5);
  TH1D *hPValueBeta0_toyevt3 = new TH1D("hPValueBeta0_toyevt3","hPValueBeta0_toyevt3",200,0,1.5);
  TH1D *hPValueBeta0_toyevt4 = new TH1D("hPValueBeta0_toyevt4","hPValueBeta0_toyevt4",200,0,1.5);
  TH1D *hPValueBeta0_toyevt5 = new TH1D("hPValueBeta0_toyevt5","hPValueBeta0_toyevt5",200,0,1.5);
  TH1D *hPValueBeta0_toyevt6 = new TH1D("hPValueBeta0_toyevt6","hPValueBeta0_toyevt6",200,0,1.5);
  TH1F *hPValueBeta0_toyevt7 = new TH1F("hPValueBeta0_toyevt7","hPValueBeta0_toyevt7",200,0,1.5);
  TH1F *hPValueBeta0_toyevt8 = new TH1F("hPValueBeta0_toyevt8","hPValueBeta0_toyevt8",200,0,1.5);
  TH1F *hPValueBeta0_toyevt9 = new TH1F("hPValueBeta0_toyevt9","hPValueBeta0_toyevt9",200,0,1.5);

  TH1D *hPValueBeta0_rateshape = new TH1D("hPValueBeta0_rateshape","hPValueBeta0_rateshape",200,0,1.5);
  TH1D *hPValueBeta0_rateshape_toyevt0 = new TH1D("hPValueBeta0_rateshape_toyevt0","hPValueBeta0_rateshape_toyevt0",200,0,1.5);
  TH1D *hPValueBeta0_rateshape_toyevt1 = new TH1D("hPValueBeta0_rateshape_toyevt1","hPValueBeta0_rateshape_toyevt1",200,0,1.5);
  TH1D *hPValueBeta0_rateshape_toyevt2 = new TH1D("hPValueBeta0_rateshape_toyevt2","hPValueBeta0_rateshape_toyevt2",200,0,1.5);
  TH1D *hPValueBeta0_rateshape_toyevt3 = new TH1D("hPValueBeta0_rateshape_toyevt3","hPValueBeta0_rateshape_toyevt3",200,0,1.5);
  TH1D *hPValueBeta0_rateshape_toyevt4 = new TH1D("hPValueBeta0_rateshape_toyevt4","hPValueBeta0_rateshape_toyevt4",200,0,1.5);
  TH1D *hPValueBeta0_rateshape_toyevt5 = new TH1D("hPValueBeta0_rateshape_toyevt5","hPValueBeta0_rateshape_toyevt5",200,0,1.5);
  TH1D *hPValueBeta0_rateshape_toyevt6 = new TH1D("hPValueBeta0_rateshape_toyevt6","hPValueBeta0_rateshape_toyevt6",200,0,1.5);
  TH1F *hPValueBeta0_rateshape_toyevt7 = new TH1F("hPValueBeta0_rateshape_toyevt7","hPValueBeta0_rateshape_toyevt7",200,0,1.5);
  TH1F *hPValueBeta0_rateshape_toyevt8 = new TH1F("hPValueBeta0_rateshape_toyevt8","hPValueBeta0_rateshape_toyevt8",200,0,1.5);
  TH1F *hPValueBeta0_rateshape_toyevt9 = new TH1F("hPValueBeta0_rateshape_toyevt9","hPValueBeta0_rateshape_toyevt9",200,0,1.5);
  
  THStack *sPValueBeta0_rateonly = new THStack("sPValueBeta0_rateonly","sPValueBeta0_rateonly");
  THStack *sPValueBeta0_rateshape = new THStack("sPValueBeta0_rateshape","sPValueBeta0_rateshape");

  
  TH1D *hPvalue_rateshape = new TH1D("hPvalue_rateshape","hPvalue_rateshape",110,0.,1.5);
  TH1D *hPvalue_rateonly  = new TH1D("hPvalue_rateonly","hPvalue_rateonly",110,0.,1.5);

  TH1D *hPvalue_rateonly_Bkg = new TH1D("hPvalue_rateonly_Bkg","hPvalue_rateonly_Bkg",110,0.,1.5);
  TH1D *hPvalue_rateshape_Bkg = new TH1D("hPvalue_rateshape_Bkg","hPvalue_rateshape_Bkg",110,0.,1.5);
  
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
  treein->SetBranchAddress("hErecoExp_RHC1Re",&hErecoExp_RHC1Re);

  treein->SetBranchAddress("ToyID",&ToyID);

  int Nentries = treein->GetEntries();

  Ntoys_beta0 = (double) Nentries; // Get # of entries for p-value calculation

  //  for(int ientry=0;ientry<Nentries;ientry++){
  for(int ientry=0;ientry<10;ientry++){

    treein->GetEntry(ientry);
    double nevents = (double)hErecoToy_RHC1Re->Integral();
    double expevents = (double)hErecoExp_RHC1Re->Integral();
    if(nevents<0){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events = " << nevents << " < 0" << endl;
      cerr << endl;
      exit(1);
    }
    
    //if(!(ientry%1000)) cout << "Entry " << ientry << endl;

    double deltachisq_01 = chisqtrue0 - chisqtrue1;
    for(int stat_throw = 0; stat_throw<10000; stat_throw++){
      double fluct = Rnd::Instance()->RndToyMCstat().Poisson(expevents);
      vec_evt_toys_beta0.push_back(fluct);
    }
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

  //  for(int ientry=0;ientry<Nentries;ientry++){
  for(int ientry=0;ientry<10;ientry++){
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
    if(nevents==0)      hPValueBeta0_toyevt0->Fill(pvalue_rateonly_itoy);
    else if(nevents==1) hPValueBeta0_toyevt1->Fill(pvalue_rateonly_itoy);
    else if(nevents==2) hPValueBeta0_toyevt2->Fill(pvalue_rateonly_itoy);
    else if(nevents==3) hPValueBeta0_toyevt3->Fill(pvalue_rateonly_itoy);
    else if(nevents==4) hPValueBeta0_toyevt4->Fill(pvalue_rateonly_itoy);
    else if(nevents==5) hPValueBeta0_toyevt5->Fill(pvalue_rateonly_itoy);
    else if(nevents==6) hPValueBeta0_toyevt6->Fill(pvalue_rateonly_itoy);
    else if(nevents==7) hPValueBeta0_toyevt7->Fill(pvalue_rateonly_itoy);
    else if(nevents==8) hPValueBeta0_toyevt8->Fill(pvalue_rateonly_itoy);
    else if(nevents>8)  hPValueBeta0_toyevt9->Fill(pvalue_rateonly_itoy);

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
     

    
    if(pvalue_rateshape_itoy>0.02)
      hPvalueVsNevt->Fill(pvalue_rateshape_itoy,nevents);
    else 
      hPvalueVsNevt->Fill(0.02,nevents); // for logarithmic scale

    // Isolate pvalue distributions with different number of events

    hPvalue_rateshape->Fill(pvalue_rateshape_itoy);
    if(nevents==0)      hPValueBeta0_rateshape_toyevt0->Fill(pvalue_rateshape_itoy);
    else if(nevents==1) hPValueBeta0_rateshape_toyevt1->Fill(pvalue_rateshape_itoy);
    else if(nevents==2) hPValueBeta0_rateshape_toyevt2->Fill(pvalue_rateshape_itoy);
    else if(nevents==3) hPValueBeta0_rateshape_toyevt3->Fill(pvalue_rateshape_itoy);
    else if(nevents==4) hPValueBeta0_rateshape_toyevt4->Fill(pvalue_rateshape_itoy);
    else if(nevents==5) hPValueBeta0_rateshape_toyevt5->Fill(pvalue_rateshape_itoy);
    else if(nevents==6) hPValueBeta0_rateshape_toyevt6->Fill(pvalue_rateshape_itoy);
    else if(nevents==7) hPValueBeta0_rateshape_toyevt7->Fill(pvalue_rateshape_itoy);
    else if(nevents==8) hPValueBeta0_rateshape_toyevt8->Fill(pvalue_rateshape_itoy);
    else if(nevents>8)  hPValueBeta0_rateshape_toyevt9->Fill(pvalue_rateshape_itoy);
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
  
//  double pvalue_evt_beta0 = GetPvalue(vec_evt_toys_beta0,Ntoys_beta0,data);
//  double pvalue_01_beta0 = GetPvalue(vec_chisq_01_beta0,Ntoys_beta0,deltachisq_beta01_data);
//
//  double pvalue_evt_beta1 = GetPvalue(vec_evt_toys_beta1,Ntoys_beta1,data);
//  double pvalue_01_beta1 = GetPvalue(vec_chisq_01_beta1,Ntoys_beta1,deltachisq_beta01_data);
  
  // Draw histograms
  
  TCanvas *cPvalue_rateonly = new TCanvas("cPvalue_rateonly","cPvalue_rateonly");
  cPvalue_rateonly->SetLogy();
  hPValueBeta0_toyevt0->SetFillColor(kBlack);
  hPValueBeta0_toyevt0->SetLineColor(kBlack);
  hPValueBeta0_toyevt1->SetFillColor(kBlue);
  hPValueBeta0_toyevt1->SetLineColor(kBlue);
  hPValueBeta0_toyevt2->SetFillColor(kGreen+2);
  hPValueBeta0_toyevt2->SetLineColor(kGreen+2);
  hPValueBeta0_toyevt3->SetFillColor(kBlue+2);
  hPValueBeta0_toyevt3->SetLineColor(kBlue+2);
  hPValueBeta0_toyevt4->SetFillStyle(3001);
  hPValueBeta0_toyevt4->SetFillColor(kBlue);
  hPValueBeta0_toyevt4->SetLineColor(kBlue);
  hPValueBeta0_toyevt5->SetFillStyle(3001);
  hPValueBeta0_toyevt5->SetFillColor(kGreen+2);
  hPValueBeta0_toyevt5->SetLineColor(kGreen+2);
  hPValueBeta0_toyevt6->SetFillStyle(3001);
  hPValueBeta0_toyevt6->SetFillColor(kBlue+2);
  hPValueBeta0_toyevt6->SetLineColor(kBlue+2);
  
  hPValueBeta0_toyevt7->SetFillStyle(3013);
  hPValueBeta0_toyevt7->SetFillColor(kBlue);
  hPValueBeta0_toyevt7->SetLineColor(kBlue);
  hPValueBeta0_toyevt8->SetFillStyle(3013);
  hPValueBeta0_toyevt8->SetFillColor(kGreen+2);
  hPValueBeta0_toyevt8->SetLineColor(kGreen+2);
  hPValueBeta0_toyevt9->SetFillStyle(3013);
  hPValueBeta0_toyevt9->SetFillColor(kBlue+2);
  hPValueBeta0_toyevt9->SetLineColor(kBlue+2);
  
  hPvalue_rateonly->GetXaxis()->SetTitle("p-value");
  hPvalue_rateonly->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateonly->DrawCopy();

  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt0);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt1);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt2);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt3);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt4);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt5);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt6);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt7);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt8);
  sPValueBeta0_rateonly->Add(hPValueBeta0_toyevt9);
  
  sPValueBeta0_rateonly->Draw("same");
  hPvalue_rateonly->DrawCopy("same");

  TLine *linePvalue_rateonly_Asimov = new TLine(pvalue_evt_beta0,0.,pvalue_evt_beta0,hPvalue_rateonly->GetMaximum());
  linePvalue_rateonly_Asimov->SetLineWidth(2);
  linePvalue_rateonly_Asimov->SetLineColor(kRed);
  linePvalue_rateonly_Asimov->DrawClone();
  
  std::stringstream mean0_rateonly;
  std::stringstream rms0_rateonly;
  mean0_rateonly << fixed << setprecision( 3 ) << hPvalue_rateonly->GetMean();
  rms0_rateonly << fixed << setprecision( 3 ) << hPvalue_rateonly->GetRMS();

  
  TLegend *leg = new TLegend(0.66,0.40,0.89,0.94);
  leg->SetBorderSize(1);
  leg->SetLineWidth(2);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);
  leg->AddEntry(hPvalue_rateonly,"Total","L");
  leg->AddEntry(hPvalue_rateonly,Form("Mean: %s",mean0_rateonly.str().c_str()),"L");
  leg->AddEntry(hPvalue_rateonly,Form("RMS: %s",rms0_rateonly.str().c_str()),"L");
  leg->AddEntry(hPValueBeta0_toyevt0,"# events = 0","F");
  leg->AddEntry(hPValueBeta0_toyevt1,"# events = 1","F");
  leg->AddEntry(hPValueBeta0_toyevt2,"# events = 2","F");
  leg->AddEntry(hPValueBeta0_toyevt3,"# events = 3","F");
  leg->AddEntry(hPValueBeta0_toyevt4,"# events = 4","F");
  leg->AddEntry(hPValueBeta0_toyevt5,"# events = 5","F");
  leg->AddEntry(hPValueBeta0_toyevt6,"# events = 6","F");
  leg->AddEntry(hPValueBeta0_toyevt7,"# events = 7","F");
  leg->AddEntry(hPValueBeta0_toyevt8,"# events = 8","F");
  leg->AddEntry(hPValueBeta0_toyevt9,"# events > 8","F");
  leg->Draw();
  
  
  //
  
  TCanvas *cPvalue_rateshape = new TCanvas("cPvalue_rateshape","cPvalue_rateshape");
  cPvalue_rateshape->SetLogy();
  hPValueBeta0_rateshape_toyevt0->SetFillColor(kBlack);
  hPValueBeta0_rateshape_toyevt0->SetLineColor(kBlack);
  hPValueBeta0_rateshape_toyevt1->SetFillColor(kBlue);
  hPValueBeta0_rateshape_toyevt1->SetLineColor(kBlue);
  hPValueBeta0_rateshape_toyevt2->SetFillColor(kGreen+2);
  hPValueBeta0_rateshape_toyevt2->SetLineColor(kGreen+2);
  hPValueBeta0_rateshape_toyevt3->SetFillColor(kBlue+2);
  hPValueBeta0_rateshape_toyevt3->SetLineColor(kBlue+2);
  hPValueBeta0_rateshape_toyevt4->SetFillStyle(3001);
  hPValueBeta0_rateshape_toyevt4->SetFillColor(kBlue);
  hPValueBeta0_rateshape_toyevt4->SetLineColor(kBlue);
  hPValueBeta0_rateshape_toyevt5->SetFillStyle(3001);
  hPValueBeta0_rateshape_toyevt5->SetFillColor(kGreen+2);
  hPValueBeta0_rateshape_toyevt5->SetLineColor(kGreen+2);
  hPValueBeta0_rateshape_toyevt6->SetFillStyle(3001);
  hPValueBeta0_rateshape_toyevt6->SetFillColor(kBlue+2);
  hPValueBeta0_rateshape_toyevt6->SetLineColor(kBlue+2);
  
  hPValueBeta0_rateshape_toyevt7->SetFillStyle(3013);
  hPValueBeta0_rateshape_toyevt7->SetFillColor(kBlue);
  hPValueBeta0_rateshape_toyevt7->SetLineColor(kBlue);
  hPValueBeta0_rateshape_toyevt8->SetFillStyle(3013);
  hPValueBeta0_rateshape_toyevt8->SetFillColor(kGreen+2);
  hPValueBeta0_rateshape_toyevt8->SetLineColor(kGreen+2);
  hPValueBeta0_rateshape_toyevt9->SetFillStyle(3013);
  hPValueBeta0_rateshape_toyevt9->SetFillColor(kBlue+2);
  hPValueBeta0_rateshape_toyevt9->SetLineColor(kBlue+2);
  
  
  hPvalue_rateshape->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape->DrawCopy();

  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt0);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt1);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt2);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt3);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt4);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt5);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt6);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt7);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt8);
  sPValueBeta0_rateshape->Add(hPValueBeta0_rateshape_toyevt9);
  
  sPValueBeta0_rateshape->Draw("same");
  hPvalue_rateshape->DrawCopy("same");
  
  std::stringstream mean0_rateshape;
  std::stringstream rms0_rateshape;
  mean0_rateshape << fixed << setprecision( 3 ) << hPvalue_rateshape->GetMean();
  rms0_rateshape << fixed << setprecision( 3 ) << hPvalue_rateshape->GetRMS();

  
  
  TLegend *leg_rateshape = new TLegend(0.66,0.40,0.89,0.94);
  leg_rateshape->SetBorderSize(1);
  leg_rateshape->SetLineWidth(2);
  leg_rateshape->SetTextFont(132);
  leg_rateshape->SetTextSize(0.05);
  leg_rateshape->AddEntry(hPvalue_rateshape,"Total","L");
  leg_rateshape->AddEntry(hPvalue_rateshape,Form("Mean: %s",mean0_rateshape.str().c_str()),"L");
  leg_rateshape->AddEntry(hPvalue_rateshape,Form("RMS: %s",rms0_rateshape.str().c_str()),"L");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt0,"# events = 0","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt1,"# events = 1","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt2,"# events = 2","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt3,"# events = 3","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt4,"# events = 4","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt5,"# events = 5","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt6,"# events = 6","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt7,"# events = 7","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt8,"# events = 8","F");
  leg_rateshape->AddEntry(hPValueBeta0_rateshape_toyevt9,"# events > 8","F");
  leg_rateshape->Draw();
  
  // hPvalue_rateshape_3evt->DrawCopy("same");
  // hPvalue_rateshape_4evt->DrawCopy("same");
  // hPvalue_rateshape_5evt->DrawCopy("same");
  // hPvalue_rateshape_6evt->DrawCopy("same");
//  TLine *linePvalue_rateshape_Asimov = new TLine(pvalue_01_beta0,0.,pvalue_01_beta0,hPvalue_rateshape->GetMaximum());
//  linePvalue_rateshape_Asimov->SetLineWidth(2);
//  linePvalue_rateshape_Asimov->SetLineColor(kRed);
//  linePvalue_rateshape_Asimov->DrawClone();
//  TLegend *legPvalue_rateshape = new TLegend(0.7,0.6,1.0,0.9);
//  legPvalue_rateshape->SetFillColor(0);
//  legPvalue_rateshape->SetLineColor(0);
//  legPvalue_rateshape->SetTextFont(132);
//  legPvalue_rateshape->SetTextSize(0.06);
//  legPvalue_rateshape->AddEntry(hPvalue_rateshape,"Toy experiments","l");
//  legPvalue_rateshape->AddEntry(linePvalue_rateshape_Asimov,"Asimov dataset","l");
//  legPvalue_rateshape->Draw();
//  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape->GetMean(),hPvalue_rateshape->GetRMS()); 
//  TLatex *lPvalue_rateshape = new TLatex(0.70,hPvalue_rateshape->GetMaximum()*0.8,name);
//  lPvalue_rateshape->SetTextSize(0.06);
//  lPvalue_rateshape->Draw();
//
//  //
//
//  TCanvas *cPvalue_rateshape_3evt = new TCanvas("cPvalue_rateshape_3evt","cPvalue_rateshape_3evt");
//  cPvalue_rateshape_3evt->SetLogy();
//  hPvalue_rateshape_3evt->SetFillStyle(0);
//  hPvalue_rateshape_3evt->GetXaxis()->SetTitle("p-value");
//  hPvalue_rateshape_3evt->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_3evt->DrawCopy();  
//  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_3evt->GetMean(),hPvalue_rateshape_3evt->GetRMS()); 
//  TLatex *lPvalue_rateshape_3evt = new TLatex(0.70,hPvalue_rateshape_3evt->GetMaximum()*0.8,name);
//  lPvalue_rateshape_3evt->SetTextSize(0.07);
//  lPvalue_rateshape_3evt->Draw();
//
//  //
//
//  TCanvas *cPvalue_rateshape_4evt = new TCanvas("cPvalue_rateshape_4evt","cPvalue_rateshape_4evt");
//  cPvalue_rateshape_4evt->SetLogy();
//  hPvalue_rateshape_4evt->SetFillStyle(0);
//  hPvalue_rateshape_4evt->GetXaxis()->SetTitle("p-value");
//  hPvalue_rateshape_4evt->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_4evt->DrawCopy();  
//  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_4evt->GetMean(),hPvalue_rateshape_4evt->GetRMS()); 
//  TLatex *lPvalue_rateshape_4evt = new TLatex(0.70,hPvalue_rateshape_4evt->GetMaximum()*0.8,name);
//  lPvalue_rateshape_4evt->SetTextSize(0.07);
//  lPvalue_rateshape_4evt->Draw();
//
//  //
//
//  TCanvas *cPvalue_rateshape_5evt = new TCanvas("cPvalue_rateshape_5evt","cPvalue_rateshape_5evt");
//  cPvalue_rateshape_5evt->SetLogy();
//  hPvalue_rateshape_5evt->SetFillStyle(0);
//  hPvalue_rateshape_5evt->GetXaxis()->SetTitle("p-value");
//  hPvalue_rateshape_5evt->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_5evt->DrawCopy();  
//  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_5evt->GetMean(),hPvalue_rateshape_5evt->GetRMS()); 
//  TLatex *lPvalue_rateshape_5evt = new TLatex(0.70,hPvalue_rateshape_5evt->GetMaximum()*0.8,name);
//  lPvalue_rateshape_5evt->SetTextSize(0.07);
//  lPvalue_rateshape_5evt->Draw();
//
//  //
//
//  TCanvas *cPvalue_rateshape_6evt = new TCanvas("cPvalue_rateshape_6evt","cPvalue_rateshape_6evt");
//  cPvalue_rateshape_6evt->SetLogy();
//  hPvalue_rateshape_6evt->SetFillStyle(0);
//  hPvalue_rateshape_6evt->GetXaxis()->SetTitle("p-value");
//  hPvalue_rateshape_6evt->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_6evt->DrawCopy();  
//  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_6evt->GetMean(),hPvalue_rateshape_6evt->GetRMS()); 
//  TLatex *lPvalue_rateshape_6evt = new TLatex(0.70,hPvalue_rateshape_6evt->GetMaximum()*0.8,name);
//  lPvalue_rateshape_6evt->SetTextSize(0.07);
//  lPvalue_rateshape_6evt->Draw();
//
//  //
//  
//  TCanvas *cPvalue_rateshape_7evt = new TCanvas("cPvalue_rateshape_7evt","cPvalue_rateshape_7evt");
//  cPvalue_rateshape_7evt->SetLogy();
//  hPvalue_rateshape_7evt->SetFillStyle(0);
//  hPvalue_rateshape_7evt->GetXaxis()->SetTitle("p-value");
//  hPvalue_rateshape_7evt->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_7evt->DrawCopy();
//  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_7evt->GetMean(),hPvalue_rateshape_7evt->GetRMS());
//  TLatex *lPvalue_rateshape_7evt = new TLatex(0.70,hPvalue_rateshape_7evt->GetMaximum()*0.8,name);
//  lPvalue_rateshape_7evt->SetTextSize(0.07);
//  lPvalue_rateshape_7evt->Draw();
//
//  TCanvas *cPvalue_rateshape_8evt = new TCanvas("cPvalue_rateshape_8evt","cPvalue_rateshape_8evt");
//  cPvalue_rateshape_8evt->SetLogy();
//  hPvalue_rateshape_8evt->SetFillStyle(0);
//  hPvalue_rateshape_8evt->GetXaxis()->SetTitle("p-value");
//  hPvalue_rateshape_8evt->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_8evt->DrawCopy();
//  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_8evt->GetMean(),hPvalue_rateshape_8evt->GetRMS());
//  TLatex *lPvalue_rateshape_8evt = new TLatex(0.70,hPvalue_rateshape_8evt->GetMaximum()*0.8,name);
//  lPvalue_rateshape_8evt->SetTextSize(0.07);
//  lPvalue_rateshape_8evt->Draw();
//
//  TCanvas *cPvalue_rateshape_9evt = new TCanvas("cPvalue_rateshape_9evt","cPvalue_rateshape_9evt");
//  cPvalue_rateshape_9evt->SetLogy();
//  hPvalue_rateshape_9evt->SetFillStyle(0);
//  hPvalue_rateshape_9evt->GetXaxis()->SetTitle("p-value");
//  hPvalue_rateshape_9evt->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_9evt->DrawCopy();
//  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_9evt->GetMean(),hPvalue_rateshape_9evt->GetRMS());
//  TLatex *lPvalue_rateshape_9evt = new TLatex(0.70,hPvalue_rateshape_9evt->GetMaximum()*0.8,name);
//  lPvalue_rateshape_9evt->SetTextSize(0.07);
//  lPvalue_rateshape_9evt->Draw();
//
//
//  TCanvas *cPvalueVsNevt = new TCanvas("hPvalueVsNevt","hPvalueVsNevt");
//  hPvalueVsNevt->GetXaxis()->SetTitle("p-value");
//  hPvalueVsNevt->GetYaxis()->SetTitle("# of events");
//  hPvalueVsNevt->Draw("colz");

  // 

//  TCanvas *cPvalue_rateonly_0p10 = new TCanvas("cPvalue_rateonly_0p10","cPvalue_rateonly_0p10");
//  cPvalue_rateonly_0p10->SetLogy();
//  hPvalue_rateonly_0p10->GetXaxis()->SetTitle("p-value (#leq 0.10)");
//  hPvalue_rateonly_0p10->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateonly_0p10->Draw();
//  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateonly_0p10->GetMean(),hPvalue_rateonly_0p10->GetRMS()); 
//  TLatex *lPvalue_rateonly_0p10 = new TLatex(0.07,hPvalue_rateonly_0p10->GetMaximum()*0.8,name);
//  lPvalue_rateonly_0p10->SetTextSize(0.07);
//  lPvalue_rateonly_0p10->Draw();
//
//  TCanvas *cPvalue_rateshape_0p10 = new TCanvas("cPvalue_rateshape_0p10","cPvalue_rateshape_0p10");
//  cPvalue_rateshape_0p10->SetLogy();
//  hPvalue_rateshape_0p10->GetXaxis()->SetTitle("p-value (#leq 0.10)");
//  hPvalue_rateshape_0p10->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_0p10->Draw();
//  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateshape_0p10->GetMean(),hPvalue_rateshape_0p10->GetRMS()); 
//  TLatex *lPvalue_rateshape_0p10 = new TLatex(0.07,hPvalue_rateshape_0p10->GetMaximum()*0.8,name);
//  lPvalue_rateshape_0p10->SetTextSize(0.07);
//  lPvalue_rateshape_0p10->Draw();
//
//  TCanvas *cPvalue_rateonly_0p05 = new TCanvas("cPvalue_rateonly_0p05","cPvalue_rateonly_0p05");
//  cPvalue_rateonly_0p05->SetLogy();
//  hPvalue_rateonly_0p05->GetXaxis()->SetTitle("p-value (#leq 0.05)");
//  hPvalue_rateonly_0p05->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateonly_0p05->Draw();
//  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateonly_0p05->GetMean(),hPvalue_rateonly_0p05->GetRMS()); 
//  TLatex *lPvalue_rateonly_0p05 = new TLatex(0.03,hPvalue_rateonly_0p05->GetMaximum()*0.8,name);
//  lPvalue_rateonly_0p05->SetTextSize(0.07);
//  lPvalue_rateonly_0p05->Draw();
//
//  TCanvas *cPvalue_rateshape_0p05 = new TCanvas("cPvalue_rateshape_0p05","cPvalue_rateshape_0p05");
//  cPvalue_rateshape_0p05->SetLogy();
//  hPvalue_rateshape_0p05->GetXaxis()->SetTitle("p-value (#leq 0.05)");
//  hPvalue_rateshape_0p05->GetYaxis()->SetTitle("Number of entries");
//  hPvalue_rateshape_0p05->Draw();
//  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateshape_0p05->GetMean(),hPvalue_rateshape_0p05->GetRMS()); 
//  TLatex *lPvalue_rateshape_0p05 = new TLatex(0.03,hPvalue_rateshape_0p05->GetMaximum()*0.8,name);
//  lPvalue_rateshape_0p05->SetTextSize(0.07);
//  lPvalue_rateshape_0p05->Draw();





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
//  cout << "//////////////// Asimov dataset //////////////" << endl;
//  cout << "Calculate p-values (beta=0):" << endl;
//  cout << "p-value rate-only = " << pvalue_evt_beta0 << endl;
//  cout << "p-value rate+shape = " << pvalue_01_beta0 << endl;
//  cout << "Calculate p-values (beta=1):" << endl;
//  cout << "p-value rate-only = " << pvalue_evt_beta1 << endl;
//  cout << "p-value rate+shape = " << pvalue_01_beta1 << endl;
//  cout << endl;
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
  
  cout << "\n\n Average p-values! \n" <<endl;
  cout << "rate only" << endl;
  cout << " 0  events: " <<  hPValueBeta0_toyevt0->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt0->Integral() <<  endl;
  cout << " 1  event : " <<  hPValueBeta0_toyevt1->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt1->Integral() << endl;
  cout << " 2  events: " <<  hPValueBeta0_toyevt2->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt2->Integral() <<  endl;
  cout << " 3  events: " <<  hPValueBeta0_toyevt3->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt3->Integral() <<  endl;
  cout << " 4  events: " <<  hPValueBeta0_toyevt4->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt4->Integral() <<  endl;
  cout << " 5  events: " <<  hPValueBeta0_toyevt5->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt5->Integral() <<  endl;
  cout << " 6  events: " <<  hPValueBeta0_toyevt6->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt6->Integral() <<  endl;
  cout << " 7  events: " <<  hPValueBeta0_toyevt7->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt7->Integral() <<  endl;
  cout << " 8  events: " <<  hPValueBeta0_toyevt8->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt8->Integral() <<  endl;
  cout << ">8  events: " <<  hPValueBeta0_toyevt9->GetMean() << "\t\t Toys: " << hPValueBeta0_toyevt9->Integral() <<  endl;

  cout << "rate+shape" << endl;
  cout << " 0  events: " <<  hPValueBeta0_rateshape_toyevt0->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt0->Integral() <<  endl;
  cout << " 1  event : " <<  hPValueBeta0_rateshape_toyevt1->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt1->Integral() << endl;
  cout << " 2  events: " <<  hPValueBeta0_rateshape_toyevt2->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt2->Integral() <<  endl;
  cout << " 3  events: " <<  hPValueBeta0_rateshape_toyevt3->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt3->Integral() <<  endl;
  cout << " 4  events: " <<  hPValueBeta0_rateshape_toyevt4->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt4->Integral() <<  endl;
  cout << " 5  events: " <<  hPValueBeta0_rateshape_toyevt5->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt5->Integral() <<  endl;
  cout << " 6  events: " <<  hPValueBeta0_rateshape_toyevt6->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt6->Integral() <<  endl;
  cout << " 7  events: " <<  hPValueBeta0_rateshape_toyevt7->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt7->Integral() <<  endl;
  cout << " 8  events: " <<  hPValueBeta0_rateshape_toyevt8->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt8->Integral() <<  endl;
  cout << ">8  events: " <<  hPValueBeta0_rateshape_toyevt9->GetMean() << "\t\t Toys: " << hPValueBeta0_rateshape_toyevt9->Integral() <<  endl;
  


  
  // Print plots

  cPvalue_rateonly->Print("cPvalue_rateonly.pdf");
  cPvalue_rateshape->Print("cPvalue_rateshape.pdf");
//  cPvalue_rateshape_3evt->Print("cPvalue_rateshape_3evt.pdf");
//  cPvalue_rateshape_4evt->Print("cPvalue_rateshape_4evt.pdf");
//  cPvalue_rateshape_5evt->Print("cPvalue_rateshape_5evt.pdf");
//  cPvalue_rateshape_6evt->Print("cPvalue_rateshape_6evt.pdf");
//  cPvalue_rateshape_7evt->Print("cPvalue_rateshape_7evt.pdf");
//  cPvalue_rateshape_8evt->Print("cPvalue_rateshape_8evt.pdf");
//  cPvalue_rateshape_9evt->Print("cPvalue_rateshape_9evt.pdf");

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
  double ntot_toys = vec_toys.size();
  /*
  if(ntot_toys != vec_toys.size()){    
    cout << endl;
    cout << "Warning:" << endl;
    cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
    cout << "Input total number of toys = " << ntot_toys << endl; 
    ntot_toys = vec_toys.size();
    cout << "Set total number of toys to vector size: " << ntot_toys << endl;
    cout << endl;
  }
   */
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
