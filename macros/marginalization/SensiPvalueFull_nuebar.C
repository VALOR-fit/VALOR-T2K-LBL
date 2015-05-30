
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
 // Mach3 - 4.5 x 10^20 POT
 double data = 4.18289, // # of data events
 double deltachisq_beta01_data = 1.93242, // data deltachisq 01
 double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 char infilename_beta0[200] = "../../outputs/marginalization/nuebarapp/fitmarg_betabar_0_4p5E20POT_DEV-3f+normbar_10kToys.root",
 char infilename_beta1[200] = "../../outputs/marginalization/nuebarapp/fitmarg_betabar_1_4p5E20POT_DEV-3f+normbar_10kToys.root",

 // // Mach3 - 4 x 10^20 POT
 // double data = 3.71812, // # of data events
 // double deltachisq_beta01_data = 1.72185, // data deltachisq 01
 // double deltachisq_beta0BF_data = 99999, // data deltachisq 0BF
 // char infilename_beta0[200] = "../../outputs/marginalization/nuebarapp/fitmarg_betabar_0_4E20POT_DEV-3f+normbar_10kToys.root",
 // char infilename_beta1[200] = "../../outputs/marginalization/nuebarapp/fitmarg_betabar_1_4E20POT_DEV-3f+normbar_10kToys.root",
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

  double chisqtrue0;
  double chisqtrue1;
  TH1D *hErecoToy;
  int ToyID;

  // Declare histograms

  TH1D *hPvalue_rateonly = new TH1D("hPvalue_rateonly","hPvalue_rateonly",110,0.,1.1); 
  TH1D *hPvalue_rateshape = new TH1D("hPvalue_rateshape","hPvalue_rateshape",110,0.,1.1); 

  TH1D *hPvalue_rateonly_3evt = new TH1D("hPvalue_rateonly_3evt","hPvalue_rateonly_3evt",110,0.,1.1); 
  TH1D *hPvalue_rateshape_3evt = new TH1D("hPvalue_rateshape_3evt","hPvalue_rateshape_3evt",110,0.,1.1); 

  TH1D *hPvalue_rateonly_4evt = new TH1D("hPvalue_rateonly_4evt","hPvalue_rateonly_4evt",110,.,1.1); 
  TH1D *hPvalue_rateshape_4evt = new TH1D("hPvalue_rateshape_4evt","hPvalue_rateshape_4evt",110,0.,1.1); 

  TH1D *hPvalue_rateonly_5evt = new TH1D("hPvalue_rateonly_5evt","hPvalue_rateonly_5evt",110,0.,1.1); 
  TH1D *hPvalue_rateshape_5evt = new TH1D("hPvalue_rateshape_5evt","hPvalue_rateshape_5evt",110,0.,1.1); 

  TH1D *hPvalue_rateonly_6evt = new TH1D("hPvalue_rateonly_6evt","hPvalue_rateonly_6evt",110,0.,1.1); 
  TH1D *hPvalue_rateshape_6evt = new TH1D("hPvalue_rateshape_6evt","hPvalue_rateshape_6evt",110,0.,1.1); 

  TH2D *hPvalueVsNevt = new TH2D("hPvalueVsNevt","hPvalueVsNevt",110,0.,1.1,12,0,12);

  // Get first file (beta=0)

  TFile *filein = new TFile(infilename_beta0,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("hErecoToy",&hErecoToy);
  treein->SetBranchAddress("ToyID",&ToyID);

  int Nentries = treein->GetEntries();

  Ntoys_beta0 = (double) Nentries; // Get # of entries for p-value calculation

  for(int ientry=0;ientry<Nentries;ientry++){
    treein->GetEntry(ientry);
    double nevents = (double)hErecoToy->Integral();    
    if(nevents<0){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events = " << nevents << " < 0" << endl;
      cerr << endl;
      exit(1);
    }

    double deltachisq_01 = chisqtrue0 - chisqtrue1;

    vec_evt_toys_beta0.push_back(nevents);
    vec_chisq_01_beta0.push_back(deltachisq_01);
  }


  // Get first file (beta=1)

  filein = new TFile(infilename_beta1,"READ");
  treein = (TTree*)filein->Get("fitnt");
  treein->SetBranchAddress("chisqtrue0",&chisqtrue0);
  treein->SetBranchAddress("chisqtrue1",&chisqtrue1);
  treein->SetBranchAddress("hErecoToy",&hErecoToy);
  treein->SetBranchAddress("ToyID",&ToyID);

  int Nentries = treein->GetEntries();

  Ntoys_beta1 = (double) Nentries; // Get # of entries for p-value calculation

  for(int ientry=0;ientry<Nentries;ientry++){
    //for(int ientry=0;ientry<100;ientry++){
    treein->GetEntry(ientry);
    if( !(ientry % 100) ) cout << "Compute pvalues for toy " << ientry << endl;
    
    double nevents = (double)hErecoToy->Integral();    
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

    if(nevents==3) hPvalue_rateshape_3evt->Fill(pvalue_rateshape_itoy);
    else if(nevents==4) hPvalue_rateshape_4evt->Fill(pvalue_rateshape_itoy);
    else if(nevents==5) hPvalue_rateshape_5evt->Fill(pvalue_rateshape_itoy);
    else if(nevents==6) hPvalue_rateshape_6evt->Fill(pvalue_rateshape_itoy);
   
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
  hPvalue_rateonly->GetYaxis()->SetTitle("# of entries");
  hPvalue_rateonly->DrawCopy();
  // hPvalue_rateonly_3evt->DrawCopy("same");
  // hPvalue_rateonly_4evt->DrawCopy("same");
  // hPvalue_rateonly_5evt->DrawCopy("same");
  // hPvalue_rateonly_6evt->DrawCopy("same");
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
  hPvalue_rateshape->GetYaxis()->SetTitle("# of entries");
  hPvalue_rateshape->DrawCopy();
  // hPvalue_rateshape_3evt->DrawCopy("same");
  // hPvalue_rateshape_4evt->DrawCopy("same");
  // hPvalue_rateshape_5evt->DrawCopy("same");
  // hPvalue_rateshape_6evt->DrawCopy("same");
  TLine *linePvalue_rateshape_Asimov = new TLine(pvalue_01_beta0,0.,pvalue_01_beta0,hPvalue_rateshape->GetMaximum());
  linePvalue_rateshape_Asimov->SetLineWidth(2);
  linePvalue_rateshape_Asimov->SetLineColor(kRed);
  linePvalue_rateshape_Asimov->Draw();
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
  hPvalue_rateshape_3evt->GetYaxis()->SetTitle("# of entries");
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
  hPvalue_rateshape_4evt->GetYaxis()->SetTitle("# of entries");
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
  hPvalue_rateshape_5evt->GetYaxis()->SetTitle("# of entries");
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
  hPvalue_rateshape_6evt->GetYaxis()->SetTitle("# of entries");
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


  
  // Print results
  
  cout << std::setprecision(10);
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
  cout << "Mean = " << hPvalue_rateonly->GetMean() << endl;
  cout << "RMS = " << hPvalue_rateonly->GetRMS() << endl;
  cout << "p-value rate+shape:" << endl;
  cout << "Mean = " << hPvalue_rateshape->GetMean() << endl;
  cout << "RMS = " << hPvalue_rateshape->GetRMS() << endl;
  cout << endl;
    
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
