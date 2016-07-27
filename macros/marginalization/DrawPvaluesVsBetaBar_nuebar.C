
void DrawPvaluesVsBetaBar_nuebar
(
 //const char *path = "../TECH_NOTE_14_3_16/BETABAR/PVALUES",
 const char *path = "../TECH_NOTE_14_3_16/BETABAR/100k",
 
 const double Asimov = 3.2183084
)
{
  char name[300];

  const int Nfiles = 5; //7;
  TString file[Nfiles] = {
    
    // 10k marg toys
    // "pvalue_Asimov_A_betabar_0p5.root",
    // "pvalue_Asimov_A_betabar_1p0.root",
    // "pvalue_Asimov_A_betabar_1p5.root",
    // "pvalue_Asimov_A_betabar_2p0.root",
    // "pvalue_Asimov_A_betabar_2p5.root"   

    // 100k marg toys
    "pvalue_Asimov_A_Marg100k_Toys100k_betabar_0p5.root", 
    "pvalue_Asimov_A_Marg100k_Toys100k_betabar_1p0.root",
    "pvalue_Asimov_A_Marg100k_Toys100k_betabar_1p5.root", 
    "pvalue_Asimov_A_Marg100k_Toys100k_betabar_2p0.root", 
    "pvalue_Asimov_A_Marg100k_Toys100k_betabar_2p5.root" 
  };
  
  double betabar[Nfiles] = {0.5, 1., 1.5, 2., 2.5}; //{0., 0.5, 1., 1.5, 2., 2.5}; //, 3.};

  TGraph *gPvalue01_Vs_betabar = new TGraph();
  TGraph *gPvalueEvt_Vs_betabar = new TGraph();
  TGraph *gPvalue01_Vs_nevents = new TGraph();
  TGraph *gPvalueEvt_Vs_nevents = new TGraph();

  for(int ifile=0;ifile<Nfiles;ifile++){

    char namein[200];
    sprintf(namein,"%s/%s",path,file[ifile].Data());
	        
    TFile *filein = new TFile(namein,"READ");
    
    TH1D *hNdata = (TH1D*) filein->Get("hNdata");
    TH1D *hdeltachisq_beta01 = (TH1D*) filein->Get("hdeltachisq_beta01");
    TH1D *hpvalue_rateonly = (TH1D*) filein->Get("hpvalue_rateonly_right");
    TH1D *hpvalue_rateshape01 = (TH1D*) filein->Get("hpvalue_rateshape01_right");
    
    double ndata = hNdata->GetBinContent(1);
    double deltachisq_beta01 = hdeltachisq_beta01->GetBinContent(1);
    double pvalueEvt_rateonly = hpvalue_rateonly->GetBinContent(1);
    double pvalue01_rateshape = hpvalue_rateshape01->GetBinContent(1);
 
    gPvalue01_Vs_betabar->SetPoint(ifile,betabar[ifile],pvalue01_rateshape);
    //gPvalueEvt_Vs_betabar->SetPoint(ifile,betabar[ifile],pvalueEvt_rateonly);

    gPvalue01_Vs_nevents->SetPoint(ifile,ndata,pvalue01_rateshape);
    //gPvalueEvt_Vs_nevents->SetPoint(ifile,ndata,pvalueEvt_rateonly);
    
    cout << endl;
    cout << "Read file " << ifile+1 << " : " << namein << endl;
    cout << "ndata = " << ndata << endl;
    cout << "deltachisq_beta01 = " << deltachisq_beta01 << endl;
    cout << "pvalue_rateonly = " << pvalueEvt_rateonly << endl;
    cout << "pvalue_rateshape01 = " << pvalue01_rateshape << endl;
    
   }


  TCanvas *cPvalues_Vs_betabar = new TCanvas("cPvalues_Vs_betabar","cPvalues_Vs_betabar");
  //cPvalues_Vs_betabar->SetLogy();
  gPvalue01_Vs_betabar->SetLineColor(kBlue);
  gPvalue01_Vs_betabar->SetMarkerSize(2.);
  gPvalue01_Vs_betabar->SetMarkerStyle(23);
  gPvalue01_Vs_betabar->SetLineWidth(2);
  gPvalue01_Vs_betabar->GetXaxis()->SetTitle("#bar{#beta}");
  gPvalue01_Vs_betabar->GetXaxis()->SetLabelOffset(0.02);
  gPvalue01_Vs_betabar->GetYaxis()->SetTitle("p-value");
  gPvalue01_Vs_betabar->GetYaxis()->SetRangeUser(0.,0.7);
  gPvalue01_Vs_betabar->GetYaxis()->SetLabelOffset(0.02);
  //gPvalueEvt_Vs_betabar->SetLineColor(kViolet);
  //gPvalueEvt_Vs_betabar->SetMarkerSize(2.);
  //gPvalueEvt_Vs_betabar->SetMarkerStyle(24);
  //gPvalueEvt_Vs_betabar->SetLineWidth(2);
  gPvalue01_Vs_betabar->Draw("alp");
  //gPvalueEvt_Vs_betabar->Draw("lp same");  
  TLine *lAsimov_betabar = new TLine(1,0.,1,0.25);
  lAsimov_betabar->SetLineColor(kRed);
  lAsimov_betabar->SetLineWidth(2);
  lAsimov_betabar->Draw();
  TLegend *lPvalues_betabar = new TLegend(0.5,0.6,0.9,0.9);
  lPvalues_betabar->SetFillColor(0);
  lPvalues_betabar->SetLineColor(0);
  lPvalues_betabar->SetTextFont(132);
  lPvalues_betabar->SetTextSize(0.06);
  lPvalues_betabar->AddEntry(gPvalue01_Vs_betabar,"Q = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)","lp");
  //lPvalues_betabar->AddEntry(gPvalueEvt_Vs_betabar,"Q = # of events","lp");
  sprintf(name,"%s","PMNS expectation"); 
  lPvalues_betabar->AddEntry(lAsimov_betabar,name,"l");
  lPvalues_betabar->Draw();




  TCanvas *cPvalues_Vs_nevents = new TCanvas("cPvalues_Vs_nevents","cPvalues_Vs_nevents");
  //cPvalues_Vs_nevents->SetLogy();
  gPvalue01_Vs_nevents->SetLineColor(kBlue);
  gPvalue01_Vs_nevents->SetMarkerSize(2.);
  gPvalue01_Vs_nevents->SetMarkerStyle(23);
  gPvalue01_Vs_nevents->SetLineWidth(2);
  gPvalue01_Vs_nevents->GetXaxis()->SetTitle("# of observed events");
  gPvalue01_Vs_nevents->GetXaxis()->SetLabelOffset(0.02);
  gPvalue01_Vs_nevents->GetYaxis()->SetTitle("p-value");
  gPvalue01_Vs_nevents->GetYaxis()->SetRangeUser(0.,0.7);
  gPvalue01_Vs_nevents->GetYaxis()->SetLabelOffset(0.02);
  //gPvalueEvt_Vs_nevents->SetLineColor(kViolet);
  //gPvalueEvt_Vs_nevents->SetMarkerSize(2.);
  //gPvalueEvt_Vs_nevents->SetMarkerStyle(24);
  //gPvalueEvt_Vs_nevents->SetLineWidth(2);
  gPvalue01_Vs_nevents->Draw("alp");
  //gPvalueEvt_Vs_nevents->Draw("lp same");
  TLine *lAsimov_nevents = new TLine(Asimov,0.,Asimov,0.25);
  lAsimov_nevents->SetLineColor(kRed);
  lAsimov_nevents->SetLineWidth(2);
  lAsimov_nevents->Draw();  
  TLegend *lPvalues_nevents = new TLegend(0.5,0.6,0.9,0.9);
  lPvalues_nevents->SetFillColor(0);
  lPvalues_nevents->SetLineColor(0);
  lPvalues_nevents->SetTextFont(132);
  lPvalues_nevents->SetTextSize(0.06);
  lPvalues_nevents->AddEntry(gPvalue01_Vs_nevents,"Q = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)","lp");
  //lPvalues_nevents->AddEntry(gPvalueEvt_Vs_nevents,"Q = # of events","lp");
  sprintf(name,"%s","PMNS expectation"); 
  lPvalues_nevents->AddEntry(lAsimov_nevents,name,"l");
  lPvalues_nevents->Draw();


}
