 

using namespace std;
using namespace valor;
using namespace valor::conv;

void test_pvalue_toys
(
 const int Ntoys = 100000,
 const double ExpSig = 2,
 const double ExpBkg = 1,
 const double data = 3,
 const double dataMax = 10
 ){  
  char name[200];

  vector<double> vec_chisq_01;
  vector<double> vec_chisq_0BF;
  vector<double> vec_evt_poiss;

  TRandom3 *rndm = new TRandom3(0);
  
  TH1D *hPoiss = new TH1D("hPoiss","hPoiss",100,0,10);
  
  TH1D *hChisq_beta0 = new TH1D("hChisq_beta0","hChisq_beta0",100,0,20); 
  TH1D *hChisq_beta0_toyevt0 = new TH1D("hChisq_beta0_toyevt0","hChisq_beta0_toyevt0",100,0,20); 
  TH1D *hChisq_beta0_toyevt1 = new TH1D("hChisq_beta0_toyevt1","hChisq_beta0_toyevt1",100,0,20); 
  THStack *sChisq_beta0 = new THStack("sChisq_beta0","sChisq_beta0"); 

  TH1D *hChisq_beta1 = new TH1D("hChisq_beta1","hChisq_beta1",100,0,20); 
  TH1D *hChisq_beta1_toyevt0 = new TH1D("hChisq_beta1_toyevt0","hChisq_beta1_toyevt0",100,0,20); 
  TH1D *hChisq_beta1_toyevt1 = new TH1D("hChisq_beta1_toyevt1","hChisq_beta1_toyevt1",100,0,20); 
  THStack *sChisq_beta1 = new THStack("sChisq_beta1","sChisq_beta1"); 

  TH1D *hChisq_betaBF = new TH1D("hChisq_betaBF","hChisq_betaBF",100,0,20); 
  TH1D *hChisq_betaBF_toyevt0 = new TH1D("hChisq_betaBF_toyevt0","hChisq_betaBF_toyevt0",100,0,20); 
  TH1D *hChisq_betaBF_toyevt1 = new TH1D("hChisq_betaBF_toyevt1","hChisq_betaBF_toyevt1",100,0,20); 
  THStack *sChisq_betaBF = new THStack("sChisq_betaBF","sChisq_betaBF"); 

  TH1D *hDeltaChisq_beta01 = new TH1D("hDeltaChisq_beta01","hDeltaChisq_beta01",100,-20,30); 
  TH1D *hDeltaChisq_beta01_toyevt0 = new TH1D("hDeltaChisq_beta01_toyevt0","hDeltaChisq_beta01_toyevt0",100,-20,30); 
  TH1D *hDeltaChisq_beta01_toyevt1 = new TH1D("hDeltaChisq_beta01_toyevt1","hDeltaChisq_beta01_toyevt1",100,-20,30); 
  THStack *sDeltaChisq_beta01 = new THStack("sDeltaChisq_beta01","sDeltaChisq_beta01"); 

  TH1D *hDeltaChisq_beta0BF = new TH1D("hDeltaChisq_beta0BF","hDeltaChisq_beta0BF",100,-20,30); 
  TH1D *hDeltaChisq_beta0BF_toyevt0 = new TH1D("hDeltaChisq_beta0BF_toyevt0","hDeltaChisq_beta0BF_toyevt0",100,-20,30); 
  TH1D *hDeltaChisq_beta0BF_toyevt1 = new TH1D("hDeltaChisq_beta0BF_toyevt1","hDeltaChisq_beta0BF_toyevt1",100,-20,30); 
  THStack *sDeltaChisq_beta0BF = new THStack("sDeltaChisq_beta0BF","sDeltaChisq_beta0BF"); 
 
  TGraph *gPvalue_01 = new TGraph();
  TGraph *gPvalue_0BF = new TGraph();
  TGraph *gPvalue_evt = new TGraph();

  double pvalue_01 = 0.;
  double pvalue_0BF = 0.;

  // Get data chisq
  
  double beta0_data = 0.;
  double nexp_beta0_data = GetNexpTot(beta0_data,ExpSig,ExpBkg);

  double beta1_data = 1.;
  double nexp_beta1_data = GetNexpTot(beta1_data,ExpSig,ExpBkg);

  double betaBF_data = GetBetaBF(data,ExpSig,ExpBkg);
  double nexp_betaBF_data = GetNexpTot(betaBF_data,ExpSig,ExpBkg);

  double chisq_beta0_data = valor::util::ChisqPoisson(data,nexp_beta0_data);
  double chisq_beta1_data = valor::util::ChisqPoisson(data,nexp_beta1_data);
  double chisq_betaBF_data = valor::util::ChisqPoisson(data,nexp_betaBF_data);

  double deltachisq_beta0BF_data = chisq_beta0_data - chisq_betaBF_data;
  double deltachisq_beta01_data = chisq_beta0_data - chisq_beta1_data;

  //if(0){
  if(1){
    cout << endl;
    cout << "# of expected signal = " << ExpSig << endl;
    cout << "# of expected background = " << ExpBkg << endl;
    cout << endl;
    cout << "Data:" << endl;
    cout << "data = " << data << endl;
    cout << "nexp_beta0_data = " << nexp_beta0_data << endl; 
    cout << "nexp_beta1_data = " << nexp_beta1_data << endl;
    cout << "nexp_betaBF_data = " << nexp_betaBF_data << endl;
    
    cout << "chisq_beta0_data = " << chisq_beta0_data << endl;
    cout << "chisq_beta1_data = " << chisq_beta1_data << endl;
    cout << "chisq_betaBF_data = " << chisq_betaBF_data << endl;
    
    cout << "deltachisq_beta0BF_data = " << deltachisq_beta0BF_data << endl;
    cout << "deltachisq_beta01_data = " << deltachisq_beta01_data << endl;
    cout << endl;
  }

  // Toys

  for(int itoy=0;itoy<Ntoys;itoy++){
    
    double ntoy_bkg = (double) rndm->Poisson(ExpBkg);
    double ntoy_tot = ntoy_bkg;
    hPoiss->Fill(ntoy_tot);
    
    double beta0 = 0.;
    double nexp_beta0 = GetNexpTot(beta0,ExpSig,ExpBkg);
  
    double beta1 = 1.;
    double nexp_beta1 = GetNexpTot(beta1,ExpSig,ExpBkg);

    double betaBF = GetBetaBF(ntoy_tot,ExpSig,ExpBkg);
    double nexp_betaBF = GetNexpTot(betaBF,ExpSig,ExpBkg);
 
    double chisq_beta0 = valor::util::ChisqPoisson(ntoy_tot,nexp_beta0);
    double chisq_beta1 = valor::util::ChisqPoisson(ntoy_tot,nexp_beta1);
    double chisq_betaBF = valor::util::ChisqPoisson(ntoy_tot,nexp_betaBF);

    double deltachisq_beta0BF = chisq_beta0 - chisq_betaBF;
    double deltachisq_beta01 = chisq_beta0 - chisq_beta1;
    

    if(0){
    //if(chisq_betaBF>1.){
    //if(chisq_beta0==2. || chisq_beta0==0.){
    //if(ntoy_tot<2.){
      //if(chisq_beta0>0.5 && chisq_beta0<1.){
      //if(deltachisq_beta01 > 0. && deltachisq_beta01<0.5){
      cout << endl;
      cout << "Toy " << itoy << endl;
      //cout << "ntoy_sig = " << ntoy_sig << endl;
      cout << "ntoy_bkg = " << ntoy_bkg << endl;
      cout << "ntoy_tot = " << ntoy_tot << endl;
      cout << "betaBF = " << betaBF << endl;
      
      cout << "ntoy_tot = " << ntoy_tot << endl;
      cout << "nexp_beta0 = " << nexp_beta0 << endl; 
      cout << "nexp_beta1 = " << nexp_beta1 << endl; 
      cout << "nexp_betaBF = " << nexp_betaBF << endl; 
      
      cout << "chisq_beta0 = " << chisq_beta0 << endl;
      cout << "chisq_beta1 = " << chisq_beta1 << endl;
      cout << "chisq_betaBF = " << chisq_betaBF << endl;

      cout << "chisq_beta0_data = " << chisq_beta0_data << endl;
      cout << "chisq_beta1_data = " << chisq_beta1_data << endl;
      cout << "chisq_betaBF_data = " << chisq_betaBF_data << endl;
 
      cout << "deltachisq_beta0BF = " << deltachisq_beta0BF << endl;
      cout << "deltachisq_beta01 = " << deltachisq_beta01 << endl;
     
      cout << "deltachisq_beta0BF_data = " << deltachisq_beta0BF_data << endl;
      cout << "deltachisq_beta01_data = " << deltachisq_beta01_data << endl;
      cout << endl;
    }
    // else{
    //   continue;
    // }

    hChisq_beta0->Fill(chisq_beta0);
    hChisq_beta1->Fill(chisq_beta1);
    hChisq_betaBF->Fill(chisq_betaBF);
    
    hDeltaChisq_beta0BF->Fill(deltachisq_beta0BF);
    hDeltaChisq_beta01->Fill(deltachisq_beta01);
    
    if(ntoy_tot==0){
      hChisq_beta0_toyevt0->Fill(chisq_beta0);
      hChisq_beta1_toyevt0->Fill(chisq_beta1);
      hChisq_betaBF_toyevt0->Fill(chisq_betaBF);
      hDeltaChisq_beta0BF_toyevt0->Fill(deltachisq_beta0BF);
      hDeltaChisq_beta01_toyevt0->Fill(deltachisq_beta01);      
    }
    else if(ntoy_tot==1){
      hChisq_beta0_toyevt1->Fill(chisq_beta0);
      hChisq_beta1_toyevt1->Fill(chisq_beta1);
      hChisq_betaBF_toyevt1->Fill(chisq_betaBF);
      hDeltaChisq_beta0BF_toyevt1->Fill(deltachisq_beta0BF);
      hDeltaChisq_beta01_toyevt1->Fill(deltachisq_beta01);      
    }

    vec_chisq_0BF.push_back(deltachisq_beta0BF);
    vec_chisq_01.push_back(deltachisq_beta01);     
    vec_evt_poiss.push_back(ntoy_tot);

  } // loop toys


  // Calculate the pvalue

  double pvalue_01 = GetPvalue(vec_chisq_01,Ntoys,deltachisq_beta01_data);
  double pvalue_0BF = GetPvalue(vec_chisq_0BF,Ntoys,deltachisq_beta0BF_data);
  double pvalue_evt = GetPvalue(vec_evt_poiss,Ntoys,data);
  
  // Print results

  cout << endl;
  cout << "p-value 01  = " << pvalue_01 << endl;
  cout << "p-value BF  = " << pvalue_0BF << endl;
  cout << "p-value evt = " << pvalue_evt << endl;
  cout << endl;


  // Get p-value for different # of data
  
  for(int idata=0;idata<=dataMax;idata++){
    
    // Get data chisq
    
    double betaBF_dataMax = GetBetaBF((double)idata,ExpSig,ExpBkg);
    double nexp_betaBF_dataMax = GetNexpTot(betaBF_dataMax,ExpSig,ExpBkg);
    
    double chisq_beta0_dataMax = valor::util::ChisqPoisson((double)idata,nexp_beta0);
    double chisq_beta1_dataMax = valor::util::ChisqPoisson((double)idata,nexp_beta1);
    double chisq_betaBF_dataMax = valor::util::ChisqPoisson((double)idata,nexp_betaBF_dataMax);
    
    double deltachisq_beta0BF_dataMax = chisq_beta0_dataMax - chisq_betaBF_dataMax;
    double deltachisq_beta01_dataMax = chisq_beta0_dataMax - chisq_beta1_dataMax;

    pvalue_01 = GetPvalue(vec_chisq_01,Ntoys,deltachisq_beta01_dataMax);
    pvalue_0BF = GetPvalue(vec_chisq_0BF,Ntoys,deltachisq_beta0BF_dataMax);
    pvalue_evt = GetPvalue(vec_evt_poiss,Ntoys,(double)idata);

    gPvalue_01->SetPoint(idata,(double)idata,pvalue_01);
    gPvalue_0BF->SetPoint(idata,(double)idata,pvalue_0BF);
    gPvalue_evt->SetPoint(idata,(double)idata,pvalue_evt);

    cout << endl;
    cout << "# of data = " << idata << endl;
    cout << "# Exp bkg = " << ExpBkg << endl;
    cout << "# Exp sig = " << ExpSig << endl;
    cout << "p-value 01 = " << pvalue_01 << endl;
    cout << "p-value BF = " << pvalue_0BF << endl;
    cout << "p-value evt = " << pvalue_evt << endl;
    cout << endl;
  }


  // Draw histograms
  
  TLine *lPoiss = new TLine(data,0,data,hPoiss->GetMaximum());
  lPoiss->SetLineColor(kRed);
  lPoiss->SetLineWidth(2);

  TLine *lDeltaChisq_beta01 = new TLine(deltachisq_beta01_data,0,deltachisq_beta01_data,hDeltaChisq_beta01->GetMaximum());
  lDeltaChisq_beta01->SetLineColor(kRed);
  lDeltaChisq_beta01->SetLineWidth(2);

  TLine *lDeltaChisq_beta0BF = new TLine(deltachisq_beta0BF_data,0,deltachisq_beta0BF_data,hDeltaChisq_beta0BF->GetMaximum());
  lDeltaChisq_beta0BF->SetLineColor(kRed);
  lDeltaChisq_beta0BF->SetLineWidth(2);
  
  TCanvas *cPoiss = new TCanvas("cPoiss","cPoiss");
  hPoiss->GetXaxis()->SetTitle("# of events");
  hPoiss->GetYaxis()->SetTitle("# of entries");
  hPoiss->DrawCopy();
  lPoiss->Draw();

  TCanvas *cChisq_beta0 = new TCanvas("cChisq_beta0","cChisq_beta0");
  hChisq_beta0->GetXaxis()->SetTitle("#chi^{2}(#bar{#beta}=0)");
  hChisq_beta0->GetYaxis()->SetTitle("# of entries");
  hChisq_beta0->DrawCopy();
  hChisq_beta0_toyevt0->SetFillColor(kGreen);
  //hChisq_beta0_toyevt0->DrawCopy("same");
  hChisq_beta0_toyevt1->SetFillColor(kBlue);
  //hChisq_beta0_toyevt1->DrawCopy("same");
  sChisq_beta0->Add(hChisq_beta0_toyevt0);
  sChisq_beta0->Add(hChisq_beta0_toyevt1);
  sChisq_beta0->Draw("same");

  TCanvas *cChisq_beta1 = new TCanvas("cChisq_beta1","cChisq_beta1");
  hChisq_beta1->GetXaxis()->SetTitle("#chi^{2}(#bar{#beta}=1)");
  hChisq_beta1->GetYaxis()->SetTitle("# of entries");
  hChisq_beta1->DrawCopy();
  hChisq_beta1_toyevt0->SetFillColor(kGreen);
  //hChisq_beta1_toyevt0->DrawCopy("same");
  hChisq_beta1_toyevt1->SetFillColor(kBlue);
  //hChisq_beta1_toyevt1->DrawCopy("same");
  sChisq_beta1->Add(hChisq_beta1_toyevt0);
  sChisq_beta1->Add(hChisq_beta1_toyevt1);
  sChisq_beta1->Draw("same");

  TCanvas *cChisq_betaBF = new TCanvas("cChisq_betaBF","cChisq_betaBF");
  hChisq_betaBF->GetXaxis()->SetTitle("#chi^{2}(#bar{#beta}_{BF})");
  hChisq_betaBF->GetYaxis()->SetTitle("# of entries");
  hChisq_betaBF->DrawCopy();
  hChisq_betaBF_toyevt0->SetFillColor(kGreen);
  //hChisq_betaBF_toyevt0->DrawCopy("same");
  hChisq_betaBF_toyevt1->SetFillColor(kBlue);
  //hChisq_betaBF_toyevt1->DrawCopy("same");
  sChisq_betaBF->Add(hChisq_betaBF_toyevt0);
  sChisq_betaBF->Add(hChisq_betaBF_toyevt1);
  sChisq_betaBF->Draw("same");

  TCanvas *cDeltaChisq_beta0BF = new TCanvas("cDeltaChisq_beta0BF","cDeltaChisq_beta0BF");
  hDeltaChisq_beta0BF->GetXaxis()->SetTitle("#chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}_{BF})");
  hDeltaChisq_beta0BF->GetYaxis()->SetTitle("# of entries");
  hDeltaChisq_beta0BF->DrawCopy();
  hDeltaChisq_beta0BF_toyevt0->SetFillColor(kGreen);
  //hDeltaChisq_beta0BF_toyevt0->DrawCopy("same");
  hDeltaChisq_beta0BF_toyevt1->SetFillColor(kBlue);
  //hDeltaChisq_beta0BF_toyevt1->DrawCopy("same");
  sDeltaChisq_beta0BF->Add(hDeltaChisq_beta0BF_toyevt0);
  sDeltaChisq_beta0BF->Add(hDeltaChisq_beta0BF_toyevt1);
  sDeltaChisq_beta0BF->Draw("same");
  lDeltaChisq_beta0BF->Draw();

  TCanvas *cDeltaChisq_beta01 = new TCanvas("cDeltaChisq_beta01","cDeltaChisq_beta01");
  hDeltaChisq_beta01->GetXaxis()->SetTitle("#chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)");
  hDeltaChisq_beta01->GetYaxis()->SetTitle("# of entries");
  hDeltaChisq_beta01->DrawCopy();
  hDeltaChisq_beta01_toyevt0->SetFillColor(kGreen);
  hDeltaChisq_beta01_toyevt0->DrawCopy("same");
  hDeltaChisq_beta01_toyevt1->SetFillColor(kBlue);
  hDeltaChisq_beta01_toyevt1->DrawCopy("same");
  sDeltaChisq_beta01->Add(hDeltaChisq_beta01_toyevt0);
  sDeltaChisq_beta01->Add(hDeltaChisq_beta01_toyevt1);
  sDeltaChisq_beta01->Draw("same");
  lDeltaChisq_beta01->Draw();

  TCanvas *cPvalues = new TCanvas("cPvalues","cPvalues");
  gPvalue_01->SetLineColor(kBlue);
  gPvalue_01->SetMarkerSize(2.);
  gPvalue_01->SetMarkerStyle(23);
  gPvalue_01->SetLineWidth(2);
  gPvalue_01->GetXaxis()->SetTitle("# of observed events");
  gPvalue_01->GetXaxis()->SetLabelOffset(0.02);
  gPvalue_01->GetYaxis()->SetTitle("p-value");
  gPvalue_01->GetYaxis()->SetRangeUser(0.,1.19);
  gPvalue_01->GetYaxis()->SetLabelOffset(0.02);
  gPvalue_0BF->SetLineColor(kRed);
  gPvalue_0BF->SetMarkerSize(2.);
  gPvalue_0BF->SetMarkerStyle(26);
  gPvalue_0BF->SetLineWidth(2);
  gPvalue_evt->SetLineColor(kViolet);
  gPvalue_evt->SetMarkerSize(2.);
  gPvalue_evt->SetMarkerStyle(24);
  gPvalue_evt->SetLineWidth(2);
  gPvalue_01->Draw("alp");
  gPvalue_0BF->Draw("lp same");
  gPvalue_evt->Draw("lp same");
  TLine *lexpbkg = new TLine(ExpBkg,0.,ExpBkg,1.19);
  lexpbkg->SetLineColor(kBlack);
  lexpbkg->SetLineWidth(2);
  lexpbkg->Draw();
  TLegend *lPvalues = new TLegend(0.5,0.6,0.9,0.9);
  lPvalues->SetFillColor(0);
  lPvalues->SetLineColor(0);
  lPvalues->SetTextFont(132);
  lPvalues->SetTextSize(0.06);
  lPvalues->AddEntry(gPvalue_01,"Q = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)","lp");
  lPvalues->AddEntry(gPvalue_0BF,"Q = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}_{BF})","lp");
  lPvalues->AddEntry(gPvalue_evt,"Q = # of events","lp");
  sprintf(name,"# of exp bkg = %d",ExpBkg); 
  lPvalues->AddEntry(lexpbkg,name,"l");
  lPvalues->Draw();
  
  return;
}

double GetNexpTot(double beta,double nexpsig,double nexpbkg){
    
  if(beta<0.){
    cerr << "Beta < 0" << endl;
    exit(1);
  }
  if(nexpsig<0.){
    cerr << "nexpsig = " << nexpsig << " < 0" << endl;
    exit(1);
  }
  if(nexpbkg<0.){
    cerr << "nexpbkg = " << nexpbkg << " < 0" << endl;
    exit(1);
  }

  double nexptot = nexpsig*beta + nexpbkg;
  if(nexptot < 0.){
    cerr << "nexptot = " << nexptot << " < 0 " << endl;
    exit(1);
  }

  return nexptot;
}

double GetBetaBF(double ndata,double nexpsig,double nexpbkg){
  
  if(ndata<0.){
    cerr << "ndata = " << ndata << " < 0" << endl;
    exit(1);
  }
  if(nexpsig<0.){
    cerr << "nexpsig = " << nexpsig << " < 0" << endl;
    exit(1);
  }
  if(nexpbkg<0.){
    cerr << "nexpbkg = " << nexpbkg << " < 0" << endl;
    exit(1);
  }

  double betabf = (ndata - nexpbkg) / nexpsig;
  
  if(betabf<0.) betabf = 0.;

  return betabf;
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
