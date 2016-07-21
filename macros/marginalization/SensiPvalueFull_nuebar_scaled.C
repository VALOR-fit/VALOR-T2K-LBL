
//
// ------------------------------------------------------------------------------------------------
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
// Compute p-values using deltaChisq_01
// Store final results in a output root file to be used by other macros
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;

void SensiPvalueFull_nuebar_scaled
(
 // Run 16 - Asimov A
 // double expbkg = 1.718,
 // double data = 3.218, // # of data events
 //double deltachisq_beta01_data = 1.41638, // data deltachisq 01
 // double pvalue_01_beta0 = 0.01, // data p-value 01
 char infilename_beta0[200] = "Toys.root",
 char infilename_beta1[200] = "Data.root",
 bool draw_beta_0_distribution = 1,
 bool use_toy_weights = 0,
 bool left_side_p_value = 0
 )
{
  char name[200];
  double pvalue_rateshape_mean_all = 0;
  
  double Ntoys_beta0;
  double Ntoys_beta1;
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
  
  TH1D *hPvalue_rateshape_7evt = new TH1D("hPvalue_rateshape_7evt","hPvalue_rateshape_7evt",110,0.,1.1);
  TH1D *hPvalue_rateshape_8evt = new TH1D("hPvalue_rateshape_8evt","hPvalue_rateshape_8evt",110,0.,1.1);
  TH1D *hPvalue_rateshape_9evt = new TH1D("hPvalue_rateshape_9evt","hPvalue_rateshape_9evt",110,0.,1.1);
  
  TH2D *hPvalueVsNevt = new TH2D("hPvalueVsNevt","hPvalueVsNevt",110,0.,1.1,12,0,12);
  
  TH1D *hPvalue_rateonly_0p05 = new TH1D("hPvalue_rateonly_0p05","hPvalue_rateonly_0p05",110,0.,0.06);
  TH1D *hPvalue_rateonly_0p10 = new TH1D("hPvalue_rateonly_0p10","hPvalue_rateonly_0p10",110,0.,0.11);
  TH1D *hPvalue_rateshape_0p05 = new TH1D("hPvalue_rateshape_0p05","hPvalue_rateshape_0p05",110,0.,0.06);
  TH1D *hPvalue_rateshape_0p10 = new TH1D("hPvalue_rateshape_0p10","hPvalue_rateshape_0p10",110,0.,0.11);
  
  TH1D* hStatisticBeta0         = new TH1D("hStatisticBeta0","hStatisticBeta0",1000,-50,50);
  TH1D* hStatisticBeta0Unscaled = new TH1D("hStatisticBeta0Unscaled","hStatisticBeta0Unscaled",1000,-50,50);
  
  TH1D* hStatisticBeta0rateonly         = new TH1D("hStatisticBeta0rateonly","hStatisticBeta0rateonly",1000,-50,50);
  TH1D* hStatisticBeta0Unscaledrateonly = new TH1D("hStatisticBeta0Unscaledrateonly","hStatisticBeta0Unscaledrateonly",1000,-50,50);
  
  
  
  // Get first file (beta=0)
  double  statistic = 0;
  double  statistic_rateonly = 0;
  TH2D*   hErecoToy_RHC1Rmu;
  TH2D*   hErecoToy_1Re;
  TH2D*   hErecoToy_1Rmu;
  TH2D*   hErecoToy_RHC1Re;
  double  hchisqtrue0;
  double  hchisqtrue1;
  int     hToyID;
  std::vector<double> right_tail_rateonly_weight;
  std::vector<double> right_tail_weight;
  std::vector<double> left_tail_rateonly_weight;
  std::vector<double> left_tail_weight;
  
  TH2D*   dErecoToy_1Re;
  TH2D*   dErecoToy_1Rmu;
  TH2D*   dErecoToy_RHC1Rmu;
  TH2D*   dErecoToy_RHC1Re;
  TH2D*   hErecoExp_RHC1Re;
  double  dchisqtrue0;
  double  dchisqtrue1;
  int     dToyID;
  
  // Get first file (beta=1)
  
  dfilein = new TFile(infilename_beta1,"READ");
  dtreein = (TTree*)dfilein->Get("fitnt");
  dtreein->SetBranchAddress("chisqtrue0",&dchisqtrue0);
  dtreein->SetBranchAddress("chisqtrue1",&dchisqtrue1);
  dtreein->SetBranchAddress("hErecoToy_1Re",&dErecoToy_1Re);
  dtreein->SetBranchAddress("hErecoToy_RHC1Re",&dErecoToy_RHC1Re);
  dtreein->SetBranchAddress("hErecoToy_RHC1Rmu",&dErecoToy_RHC1Rmu);
  dtreein->SetBranchAddress("hErecoToy_1Rmu",&dErecoToy_1Rmu);
  
  
  hfilein = new TFile(infilename_beta0,"READ");
  htreein = (TTree*)hfilein->Get("fitnt");
  htreein->SetBranchAddress("chisqtrue0",&hchisqtrue0);
  htreein->SetBranchAddress("chisqtrue1",&hchisqtrue1);
  htreein->SetBranchAddress("hErecoExp_1Re",&hErecoToy_1Re);
  htreein->SetBranchAddress("hErecoToy_RHC1Re",&hErecoToy_RHC1Re);
  htreein->SetBranchAddress("hErecoExp_RHC1Re",&hErecoExp_RHC1Re);
  htreein->SetBranchAddress("hErecoExp_RHC1Rmu",&hErecoToy_RHC1Rmu);
  htreein->SetBranchAddress("hErecoExp_1Rmu",&hErecoToy_1Rmu);
  
  
  int dNentries = dtreein->GetEntries();
  int hNentries = htreein->GetEntries();
  
  double pvalue_rateshape_itoy         = 0;
  double pvalue_error                  = 0;
  double pvalue_error_vector           = 0;
  
  double pvalue_rateonly_itoy          = 0;
  double pvalue_rateonly_error         = 0;
  double pvalue_rateonly_error_vector  = 0;
  double ddeltachisq_01                = 0;
  
  for(int ientry=0;ientry<dNentries;ientry++){
    hStatisticBeta0->Reset();
    hStatisticBeta0Unscaled->Reset();
    hStatisticBeta0rateonly->Reset();
    hStatisticBeta0Unscaledrateonly->Reset();
    dtreein->GetEntry(ientry);
    ddeltachisq_01 = dchisqtrue0 - dchisqtrue1;
    
    if( !(ientry % 1000) ) cout << "Compute pvalues for toy " << ientry << endl;
    double dnevents = (double)dErecoToy_RHC1Re->Integral();
    if(dnevents<0.){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events = " << dnevents << " < 0" << endl;
      cerr << endl;
      exit(1);
    }
    
    for(int hientry=0;hientry<hNentries;hientry++){
//         for(int hientry=0;hientry<10;hientry++){
      htreein->GetEntry(hientry);
      double hnevents = (double)hErecoToy_RHC1Re->Integral();
      if(hnevents<0){
        cerr << endl;
        cerr << "ToyID: " << hToyID << endl;
        cerr << "# of events = " << hnevents << " < 0" << endl;
        cerr << endl;
        exit(1);
      }
      double hdeltachisq_01 = hchisqtrue0 - hchisqtrue1;
      //      cout << "nevetns: " << hnevents << "  statistic: " << hdeltachisq_01 << endl;
      double toy_weight = 1;
      if(use_toy_weights){
        double toy_chisq = 0;
        toy_chisq += calc_chisq(dErecoToy_RHC1Rmu,hErecoToy_RHC1Rmu);
        toy_chisq += calc_chisq(dErecoToy_1Rmu,hErecoToy_1Rmu);
        toy_chisq += calc_chisq(dErecoToy_1Re,hErecoToy_1Re);
        //toy_chisq += calc_chisq(dErecoToy_RHC1Re,hErecoExp_RHC1Re);
        toy_weight = TMath::Exp(-1*toy_chisq/2);
        
      }
      hStatisticBeta0->Fill(hdeltachisq_01,toy_weight);
      hStatisticBeta0Unscaled->Fill(hdeltachisq_01);
      
      if(hdeltachisq_01<ddeltachisq_01)
      {
        left_tail_weight.push_back(toy_weight);
      }
      else{
        right_tail_weight.push_back(toy_weight);
      }
      
      double number_nuebar_events=hErecoExp_RHC1Re->Integral();
      for(int n_rateonly_throws=0; n_rateonly_throws<10000;n_rateonly_throws++){
        double fluct = Rnd::Instance()->RndToyMCstat().Poisson(number_nuebar_events);
        hStatisticBeta0rateonly->Fill(fluct,toy_weight);
        hStatisticBeta0Unscaledrateonly->Fill(fluct);
        //        hStatisticBeta0rateonly->Fill(hnevents,toy_weight);
        //        hStatisticBeta0Unscaledrateonly->Fill(hnevents);
        
        if(left_side_p_value)
        {
          //          if(hnevents <= dnevents){
          if(fluct <= dnevents){
            left_tail_rateonly_weight.push_back(toy_weight);
          }
          else{
            right_tail_rateonly_weight.push_back(toy_weight);
          }
          
        }
        else{
          //        if(hnevents < dnevents){
          if(fluct < dnevents){
            left_tail_rateonly_weight.push_back(toy_weight);
          }
          else{
            right_tail_rateonly_weight.push_back(toy_weight);
          }
        }
      }
    }
    
    
    
    cout << "statistic: " << ddeltachisq_01 << endl;
    
    if(draw_beta_0_distribution){
      statistic = ddeltachisq_01;
      statistic_rateonly = dnevents;
      cout << "statistic: " << statistic << endl;
    }
    
    if(!left_side_p_value){ //calculate p values and erros in the regular way
      pvalue_rateshape_itoy = GetPvalue(ddeltachisq_01,right_tail_weight,left_tail_weight);
      pvalue_error_vector   = get_error_vector(right_tail_weight,left_tail_weight,pvalue_rateshape_itoy);
      pvalue_error          = get_error(ddeltachisq_01,hStatisticBeta0,hStatisticBeta0Unscaled);
      pvalue_rateonly_itoy  = GetPvalue(statistic_rateonly,right_tail_rateonly_weight,left_tail_rateonly_weight);
      pvalue_rateonly_error_vector = get_error_vector(right_tail_rateonly_weight,left_tail_rateonly_weight,pvalue_rateonly_itoy);
      pvalue_rateonly_error = get_error(ddeltachisq_01,hStatisticBeta0rateonly,hStatisticBeta0Unscaledrateonly);
    }
    else{ //swap left and right tails to get p value
      pvalue_rateshape_itoy = GetPvalue(ddeltachisq_01,left_tail_weight,right_tail_weight);
      pvalue_error_vector   = get_error_vector(left_tail_weight,right_tail_weight,pvalue_rateshape_itoy);
      pvalue_error = get_left_error(ddeltachisq_01,hStatisticBeta0,hStatisticBeta0Unscaled);
      pvalue_rateonly_itoy  = GetPvalue(statistic_rateonly,left_tail_rateonly_weight,right_tail_rateonly_weight);
      pvalue_rateonly_error_vector = get_error_vector(left_tail_rateonly_weight,right_tail_rateonly_weight,pvalue_rateonly_itoy);
      pvalue_rateonly_error = get_left_error(statistic_rateonly,hStatisticBeta0rateonly,hStatisticBeta0Unscaledrateonly);
    }
    
    if(pvalue_rateshape_itoy<0. || pvalue_rateshape_itoy>1.){
      cerr << endl;
      cerr << "ToyID = " << ToyID << endl;
      cerr << "pvalue_rateshape = " << pvalue_rateshape_itoy << endl;
      cerr << endl;
      exit(1);
    }
    
    
    hPvalue_rateshape->Fill(pvalue_rateshape_itoy);
    
    if(pvalue_rateshape_itoy>0.02)
    hPvalueVsNevt->Fill(pvalue_rateshape_itoy,dnevents);
    else
    hPvalueVsNevt->Fill(0.02,dnevents); // for logarithmic scale
    
    // Isolate pvalue distributions with different number of events
    
    if(dnevents==3) hPvalue_rateshape_3evt->Fill(pvalue_rateshape_itoy);
    else if(dnevents==4) hPvalue_rateshape_4evt->Fill(pvalue_rateshape_itoy);
    else if(dnevents==5) hPvalue_rateshape_5evt->Fill(pvalue_rateshape_itoy);
    else if(dnevents==6) hPvalue_rateshape_6evt->Fill(pvalue_rateshape_itoy);
    else if(dnevents==7) hPvalue_rateshape_7evt->Fill(pvalue_rateshape_itoy);
    else if(dnevents==8) hPvalue_rateshape_8evt->Fill(pvalue_rateshape_itoy);
    else if(dnevents==9) hPvalue_rateshape_9evt->Fill(pvalue_rateshape_itoy);
    
    //
    
    pvalue_rateshape_mean_all += pvalue_rateshape_itoy / (double)dNentries;
    
    if(pvalue_rateshape_itoy<=0.05){
      hPvalue_rateshape_0p05->Fill(pvalue_rateshape_itoy);
    }
    if(pvalue_rateshape_itoy<=0.10){
      hPvalue_rateshape_0p10->Fill(pvalue_rateshape_itoy);
    }
    
    
    
    // Isolate p-value for # obs events > # expected bkg
    //    if(dnevents > expbkg){
    //      hPvalue_rateshape_Bkg->Fill(pvalue_rateshape_itoy);
    //    }
    
  } // end loop
  
  
  
  // Draw histograms
  
  
  //
  
  TCanvas *cPvalue_rateshape = new TCanvas("cPvalue_rateshape","cPvalue_rateshape");
  cPvalue_rateshape->SetLogy();
  hPvalue_rateshape_3evt->SetFillColor(kBlue+2);
  hPvalue_rateshape_4evt->SetFillColor(kBlue);
  hPvalue_rateshape_4evt->SetFillStyle(3001);
  hPvalue_rateshape_5evt->SetFillColor(kGreen+2);
  hPvalue_rateshape_5evt->SetFillStyle(3001);
  hPvalue_rateshape_6evt->SetFillColor(kBlue+2);
  hPvalue_rateshape_6evt->SetFillStyle(3001);
  
  hPvalue_rateshape_7evt->SetFillColor(kBlue);
  hPvalue_rateshape_7evt->SetFillStyle(3013);
  hPvalue_rateshape_8evt->SetFillColor(kGreen+2);
  hPvalue_rateshape_8evt->SetFillStyle(3013);
  hPvalue_rateshape_9evt->SetFillColor(kBlue+2);
  hPvalue_rateshape_9evt->SetFillStyle(3013);
  hPvalue_rateshape->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape->DrawCopy();
  // hPvalue_rateshape_3evt->DrawCopy("same");
  // hPvalue_rateshape_4evt->DrawCopy("same");
  // hPvalue_rateshape_5evt->DrawCopy("same");
  // hPvalue_rateshape_6evt->DrawCopy("same");
  //  TLine *linePvalue_rateshape_Asimov = new TLine(pvalue_01_beta0,0.,pvalue_01_beta0,hPvalue_rateshape->GetMaximum());
  //  linePvalue_rateshape_Asimov->SetLineWidth(2);
  //  linePvalue_rateshape_Asimov->SetLineColor(kRed);
  //  linePvalue_rateshape_Asimov->DrawClone();
  TLegend *legPvalue_rateshape = new TLegend(0.7,0.6,1.0,0.9);
  legPvalue_rateshape->SetFillColor(0);
  legPvalue_rateshape->SetLineColor(0);
  legPvalue_rateshape->SetTextFont(132);
  legPvalue_rateshape->SetTextSize(0.06);
  legPvalue_rateshape->AddEntry(hPvalue_rateshape,"Toy experiments","l");
  //  legPvalue_rateshape->AddEntry(linePvalue_rateshape_Asimov,"Asimov dataset","l");
  legPvalue_rateshape->Draw();
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape->GetMean(),hPvalue_rateshape->GetRMS());
  TLatex *lPvalue_rateshape = new TLatex(0.70,hPvalue_rateshape->GetMaximum()*0.8,name);
  lPvalue_rateshape->SetTextSize(0.06);
  lPvalue_rateshape->Draw();
  
  //
  TCanvas *cPvalue_rateshape_fit;
  if(draw_beta_0_distribution){
    cPvalue_rateshape_fit = new TCanvas("cPvalue_rateshape_fit","cPvalue_rateshape_fit");
    cPvalue_rateshape_fit->SetLogy();
    //    hStatisticBeta0->Rebin(50000);
    hStatisticBeta0->Scale(1./hStatisticBeta0->Integral());
    hStatisticBeta0->SetFillStyle(0);
    hStatisticBeta0->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(#beta=0) - #chi^{2}(#beta=1)");
    hStatisticBeta0->GetYaxis()->SetTitle("Number of entries");
    hStatisticBeta0->GetYaxis()->SetRangeUser(0.0000001,hStatisticBeta0->GetMaximum()*1.2);
    hStatisticBeta0->GetXaxis()->SetRangeUser(-5,26);
    hStatisticBeta0->DrawCopy();
    TLine *linePvalue_rateshape_statistic = new TLine(statistic,0.,statistic,hStatisticBeta0->GetMaximum()*1.2);
    linePvalue_rateshape_statistic->SetLineWidth(3);
    linePvalue_rateshape_statistic->SetLineColor(kRed);
    linePvalue_rateshape_statistic->DrawClone();
    sprintf(name,"#Delta #chi^{2} = %1.3f",statistic);
    TLatex *fit_pvalue = new TLatex(4.5,hStatisticBeta0->GetMaximum()*0.1,name);
    fit_pvalue->SetTextSize(0.07);
    fit_pvalue->Draw();
    sprintf(name,"p value = %1.3f #pm %1.4f",pvalue_rateshape_itoy,pvalue_error_vector);
    TLatex *fit_pvalue_2 = new TLatex(4.5,hStatisticBeta0->GetMaximum()*0.01,name);
    fit_pvalue_2->SetTextSize(0.07);
    fit_pvalue_2->Draw();
    
    
    
    //    cPvalue_rateonly_fit = new TCanvas("cPvalue_rateonly_fit","cPvalue_rateonly_fit");
    //    cPvalue_rateonly_fit->SetLogy();
    //    //    hStatisticBeta0rateonly->Rebin(50000);
    //    hStatisticBeta0rateonly->SetFillStyle(0);
    //    hStatisticBeta0rateonly->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(#beta=0) - #chi^{2}(#beta=1)");
    //    hStatisticBeta0rateonly->GetYaxis()->SetTitle("Number of entries");
    //    hStatisticBeta0rateonly->GetYaxis()->SetRangeUser(0.01,hStatisticBeta0rateonly->GetMaximum()*1.2);
    //    hStatisticBeta0rateonly->GetXaxis()->SetRangeUser(-1,24);
    //    hStatisticBeta0rateonly->DrawCopy();
    //
    //
    TCanvas *cBeta0Distribution = new TCanvas("cBeta0Distribution","cBeta0Distribution");
    hStatisticBeta0rateonly->GetYaxis()->SetTitle("Number of entries");
    hStatisticBeta0rateonly->GetXaxis()->SetTitle("Number of events");
    hStatisticBeta0rateonly->SetLineColor(kBlue);
    hStatisticBeta0rateonly->SetLineWidth(5);
    hStatisticBeta0rateonly->Scale(1./hStatisticBeta0rateonly->Integral());
    hStatisticBeta0rateonly->GetYaxis()->SetRangeUser(0.0000001,hStatisticBeta0rateonly->GetMaximum()*1.2);
    hStatisticBeta0rateonly->GetYaxis()->SetTitleSize(0.06);
    hStatisticBeta0rateonly->GetXaxis()->SetTitleSize(0.06);
    hStatisticBeta0rateonly->GetYaxis()->SetTitleFont(132);
    hStatisticBeta0rateonly->GetXaxis()->SetTitleFont(132);
    hStatisticBeta0rateonly->GetYaxis()->SetLabelFont(132);
    hStatisticBeta0rateonly->GetXaxis()->SetLabelFont(132);
    hStatisticBeta0rateonly->GetYaxis()->SetLabelSize(0.05);
    hStatisticBeta0rateonly->GetXaxis()->SetLabelSize(0.05);
    hStatisticBeta0rateonly->Scale(1./hStatisticBeta0rateonly->Integral());
    hStatisticBeta0rateonly->DrawCopy();
    cBeta0Distribution->SetLogy();
    
    TLine *linePvalue_rateonly_statistic = new TLine(statistic_rateonly,0.,statistic_rateonly,hStatisticBeta0rateonly->GetMaximum()*1.2);
    linePvalue_rateonly_statistic->SetLineWidth(3);
    linePvalue_rateonly_statistic->SetLineColor(kRed);
    linePvalue_rateonly_statistic->DrawClone();
    sprintf(name,"Data Events = %1.3f",statistic_rateonly);
    TLatex *fit_rateonly_pvalue = new TLatex(4.5,hStatisticBeta0rateonly->GetMaximum()*0.5,name);
    fit_rateonly_pvalue->SetTextSize(0.07);
    fit_rateonly_pvalue->Draw();
    sprintf(name,"p value = %1.4f #pm %1.5f",pvalue_rateonly_itoy,pvalue_rateonly_error_vector);
    TLatex *fit_rateonly_pvalue_2 = new TLatex(4.5,hStatisticBeta0rateonly->GetMaximum()*0.1,name);
    fit_rateonly_pvalue_2->SetTextSize(0.07);
    fit_rateonly_pvalue_2->Draw();
    
  }
  
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
  
  TCanvas *cPvalue_rateshape_7evt = new TCanvas("cPvalue_rateshape_7evt","cPvalue_rateshape_7evt");
  cPvalue_rateshape_7evt->SetLogy();
  hPvalue_rateshape_7evt->SetFillStyle(0);
  hPvalue_rateshape_7evt->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape_7evt->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_7evt->DrawCopy();
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_7evt->GetMean(),hPvalue_rateshape_7evt->GetRMS());
  TLatex *lPvalue_rateshape_7evt = new TLatex(0.70,hPvalue_rateshape_7evt->GetMaximum()*0.8,name);
  lPvalue_rateshape_7evt->SetTextSize(0.07);
  lPvalue_rateshape_7evt->Draw();
  
  TCanvas *cPvalue_rateshape_8evt = new TCanvas("cPvalue_rateshape_8evt","cPvalue_rateshape_8evt");
  cPvalue_rateshape_8evt->SetLogy();
  hPvalue_rateshape_8evt->SetFillStyle(0);
  hPvalue_rateshape_8evt->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape_8evt->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_8evt->DrawCopy();
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_8evt->GetMean(),hPvalue_rateshape_8evt->GetRMS());
  TLatex *lPvalue_rateshape_8evt = new TLatex(0.70,hPvalue_rateshape_8evt->GetMaximum()*0.8,name);
  lPvalue_rateshape_8evt->SetTextSize(0.07);
  lPvalue_rateshape_8evt->Draw();
  
  TCanvas *cPvalue_rateshape_9evt = new TCanvas("cPvalue_rateshape_9evt","cPvalue_rateshape_9evt");
  cPvalue_rateshape_9evt->SetLogy();
  hPvalue_rateshape_9evt->SetFillStyle(0);
  hPvalue_rateshape_9evt->GetXaxis()->SetTitle("p-value");
  hPvalue_rateshape_9evt->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_9evt->DrawCopy();
  sprintf(name,"#splitline{Mean=%1.2f}{RMS=%1.2f}",hPvalue_rateshape_9evt->GetMean(),hPvalue_rateshape_9evt->GetRMS());
  TLatex *lPvalue_rateshape_9evt = new TLatex(0.70,hPvalue_rateshape_9evt->GetMaximum()*0.8,name);
  lPvalue_rateshape_9evt->SetTextSize(0.07);
  lPvalue_rateshape_9evt->Draw();
  
  
  TCanvas *cPvalueVsNevt = new TCanvas("hPvalueVsNevt","hPvalueVsNevt");
  hPvalueVsNevt->GetXaxis()->SetTitle("p-value");
  hPvalueVsNevt->GetYaxis()->SetTitle("# of events");
  hPvalueVsNevt->Draw("colz");
  
  //
  
  
  TCanvas *cPvalue_rateshape_0p10 = new TCanvas("cPvalue_rateshape_0p10","cPvalue_rateshape_0p10");
  cPvalue_rateshape_0p10->SetLogy();
  hPvalue_rateshape_0p10->GetXaxis()->SetTitle("p-value (#leq 0.10)");
  hPvalue_rateshape_0p10->GetYaxis()->SetTitle("Number of entries");
  hPvalue_rateshape_0p10->Draw();
  sprintf(name,"#splitline{Mean=%1.4f}{RMS=%1.4f}",hPvalue_rateshape_0p10->GetMean(),hPvalue_rateshape_0p10->GetRMS());
  TLatex *lPvalue_rateshape_0p10 = new TLatex(0.07,hPvalue_rateshape_0p10->GetMaximum()*0.8,name);
  lPvalue_rateshape_0p10->SetTextSize(0.07);
  lPvalue_rateshape_0p10->Draw();
  
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
  //  sprintf(name,"Exp Bkg = %1.3f",expbkg);
  //  TLatex *lExpBkg = new TLatex(0.4,hPvalue_rateshape_Bkg->GetMaximum()*0.8,name);
  //  lExpBkg->SetTextSize(0.07);
  //  lExpBkg->Draw();
  
  
  
  
  
  // Print results
  
  cout << std::setprecision(10);
  cout << endl;
  //  cout << "//////////////// Asimov dataset //////////////" << endl;
  //  cout << "Calculate p-values (beta=0):" << endl;
  //  cout << "p-value rate+shape = " << pvalue_01_beta0 << endl;
  //  cout << endl;
  cout << "//////////////// Toys p-value (all toys) //////////////" << endl;
  cout << endl;
  cout << "p-value rate+shape:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateshape->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateshape->GetRMS() << endl;
  cout << "Mean = " << pvalue_rateshape_mean_all << endl;
  cout << endl;
  cout << "//////////////// Toys p-value (<0.10) //////////////" << endl;
  cout << endl;
  cout << "Fraction of toys (rate+shape) = "
  << hPvalue_rateshape_0p10->Integral() / hPvalue_rateshape->Integral() << endl;
  cout << "p-value rate-only:" << pvalue_rateonly_itoy  <<  "  Error1: " << pvalue_rateonly_error <<  "  Error2: " << pvalue_rateonly_error_vector <<  endl;
  cout << "p-value rate+shape:" << pvalue_rateshape_itoy   <<"  Error1: " <<  pvalue_error << "  Error2: " << pvalue_error_vector  <<  endl;
  cout << "Mean (histo) = " << hPvalue_rateshape_0p10->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateshape_0p10->GetRMS() << endl;
  cout << endl;
  cout << "//////////////// Toys p-value (<0.05) //////////////" << endl;
  cout << endl;
  cout << "Fraction of toys (rate+shape) = "
  << hPvalue_rateshape_0p05->Integral() / hPvalue_rateshape->Integral() << endl;
  cout << "p-value rate+shape:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateshape_0p05->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateshape_0p05->GetRMS() << endl;
  cout << endl;
  cout << "//////////////// Toys p-value (# obs evts > # exp bkg) //////////////" << endl;
  cout << endl;
  cout << "Fraction of toys (rate+shape) = "
  << hPvalue_rateshape_Bkg->Integral() / hPvalue_rateshape->Integral() << endl;
  cout << "p-value rate+shape:" << endl;
  cout << "Mean (histo) = " << hPvalue_rateshape_Bkg->GetMean() << endl;
  cout << "RMS (histo) = " << hPvalue_rateshape_Bkg->GetRMS() << endl;
  cout << endl;
  
  // Print plots
  if(draw_beta_0_distribution){
    cPvalue_rateshape_fit->Print("cPvalue_rateshape_fit.pdf");
    cBeta0Distribution->Print("cBeta0Distribution.pdf");
  }
  cPvalue_rateshape->Print("cPvalue_rateshape.pdf");
  cPvalue_rateshape_3evt->Print("cPvalue_rateshape_3evt.pdf");
  cPvalue_rateshape_4evt->Print("cPvalue_rateshape_4evt.pdf");
  cPvalue_rateshape_5evt->Print("cPvalue_rateshape_5evt.pdf");
  cPvalue_rateshape_6evt->Print("cPvalue_rateshape_6evt.pdf");
  cPvalue_rateshape_7evt->Print("cPvalue_rateshape_7evt.pdf");
  cPvalue_rateshape_8evt->Print("cPvalue_rateshape_8evt.pdf");
  cPvalue_rateshape_9evt->Print("cPvalue_rateshape_9evt.pdf");
  
  cPvalueVsNevt->Print("cPvalueVsNevt.pdf");
  cPvalue_rateshape_0p10->Print("cPvalue_rateshape_0p10.pdf");
  cPvalue_rateshape_0p05->Print("cPvalue_rateshape_0p05.pdf");
  cPvalue_rateshape_Bkg->Print("cPvalue_rateshape_Bkg.pdf");
  
  if(draw_beta_0_distribution){
    TFile f_out("SensiPValueFullNueBarScaled.root","recreate");
    hStatisticBeta0->Write("hStatisticBeta0");
    hStatisticBeta0Unscaledrateonly->Write("hStatisticBeta0Unscaledrateonly");
    hStatisticBeta0Unscaled->Write("hStatisticBeta0Unscaled");
    hStatisticBeta0rateonly->Write("hStatisticBeta0rateonly");
    f_out.Close();
  }
  
  return;
}





double calc_chisq(TH2D* hobs, TH2D* hexp)
{
  int nbins_x = hobs->GetNbinsX();
  int nbins_y = hobs->GetNbinsY();
  
  double chisq = 0;
  
  for(int i=1;i<= nbins_x;i++)
  {
    for(int j=1;j<= nbins_y;j++)
    {
      double obs = hobs->GetBinContent(i,j);
      double exp = hexp->GetBinContent(i,j);
      
      chisq += valor::util::ChisqPoisson(obs,exp);
    }
  }
  
  return chisq;
}

double GetPvalue(double _statistic, TH1D* hstatistic)
{
  int bin = hstatistic->FindBin(_statistic);
  double p_value = hstatistic->Integral(bin,hstatistic->GetNbinsX()) / hstatistic->Integral();
  return p_value;
}

double GetPvalue(double _statistic, std::vector<double> _right_tail_weight, std::vector<double> _left_tail_weight)
{
  double p_value = 1;
  double right_tail_sum=0;
  double left_tail_sum = 0;
  
  for(int i =0;i<_right_tail_weight.size();i++){
    right_tail_sum += _right_tail_weight[i];
  }
  
  for(int i =0;i<_left_tail_weight.size();i++){
    left_tail_sum += _left_tail_weight[i];
  }
  
  p_value = right_tail_sum / (right_tail_sum+ left_tail_sum);
  return p_value;
}

double get_error(double obs, TH1D* _hbeta, TH1D* _hbetaUnscaled)
{
  int obs_bin = _hbeta->FindBin(obs);
  double pvalue = _hbeta->Integral(obs_bin,_hbeta->GetXaxis()->GetNbins()) / _hbeta->Integral();
  double err_sqrd = 0;
  //  for(int i=1;i<obs_bin;i++)
  //  {
  //    double pd = _hbeta->Integral();
  //    double wl = _hbeta->GetBinContent(i)/_hbetaUnscaled->GetBinContent(i);
  //    double bl = _hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wl* (pvalue) / pd, 2) * bl;
  //    err_sqrd = err_sqrd + err;
  //  }
  
  //  for(int i=obs_bin;i<=_hbeta->GetXaxis()->GetNbins();i++)
  //  {
  //    double pd = _hbeta->Integral();
  //    double wr = _hbeta->GetBinContent(i)/_hbetaUnscaled->GetBinContent(i);
  //    double br = _hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wr* (1 - pvalue) / pd, 2) * br;
  //    err_sqrd = err_sqrd + err;
  //  }
  
  
  //Left error
  double pd = _hbeta->Integral();
  double wl = _hbeta->Integral(1,obs_bin-1)/_hbetaUnscaled->Integral(1,obs_bin-1);
  double bl = _hbetaUnscaled->Integral(1,obs_bin-1);
  double err = TMath::Power(wl* (pvalue) / (pd), 2) * bl;
  err_sqrd = err_sqrd + err;
  
  //Right error
  double pd = _hbeta->Integral();
  double wr = _hbeta->Integral(obs_bin,_hbeta->GetXaxis()->GetNbins())/_hbetaUnscaled->Integral(obs_bin,_hbeta->GetXaxis()->GetNbins());
  double br = _hbetaUnscaled->Integral(obs_bin,_hbeta->GetXaxis()->GetNbins());
  double err = TMath::Power(wr* (1 - pvalue) / (pd), 2) * br;
  err_sqrd = err_sqrd + err;
  
  double final_error = TMath::Sqrt(err_sqrd);
  
  return final_error;
  
}

double get_error_vector( std::vector<double> _right_tail_weight, std::vector<double> _left_tail_weight, double _pvalue)
{
  double error = 0;
  double total_right=0;
  double total_left=0;
  for(int i=0;i<_right_tail_weight.size();i++)
  {
    error+= TMath::Power(_right_tail_weight[i]*(1-_pvalue),2);
    total_right+= _right_tail_weight[i];
  }
  
  for(int i=0;i<_left_tail_weight.size();i++)
  {
    error+= TMath::Power(_left_tail_weight[i]*(_pvalue),2);
    total_left+= _left_tail_weight[i];
  }
  
  error = error / (TMath::Power(total_left + total_right ,2));
  return TMath::Sqrt(error);
}

double get_left_error(double obs, TH1D* hbeta, TH1D* hbetaUnscaled)
{
  int obs_bin = hbeta->FindBin(obs);
  double pvalue = hbeta->Integral(1,obs_bin) / hbeta->Integral() ;
  double err_sqrd = 0;
  //  for(int i=1;i<obs_bin;i++)
  //  {
  //    double pd = hbeta->Integral();
  //    double wl = hbeta->GetBinContent(i)/hbetaUnscaled->GetBinContent(i);
  //    double bl = hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wl* (pvalue) / pd, 2) * bl;
  //    err_sqrd = err_sqrd + err;
  //  }
  
  //Right error
  double pd = hbeta->Integral();
  double wr = hbeta->Integral(obs_bin+1,hbeta->GetXaxis()->GetNbins())/hbetaUnscaled->Integral(obs_bin+1,hbeta->GetXaxis()->GetNbins());
  double br = hbetaUnscaled->Integral(obs_bin +1,hbeta->GetXaxis()->GetNbins());
  double err = TMath::Power(wr* (pvalue) / pd, 2) * br;
  err_sqrd = err_sqrd + err;
  
  //Left error
  double pd = hbeta->Integral();
  double wl = hbeta->Integral(1,obs_bin)/hbetaUnscaled->Integral(1,obs_bin);
  double bl = hbetaUnscaled->Integral(1,obs_bin);
  double err = TMath::Power(wl* (1 - pvalue) / pd, 2) * bl;
  err_sqrd = err_sqrd + err;
  
  //  for(int i=obs_bin;i<=hbeta->GetXaxis()->GetNbins();i++)
  //  {
  //    double pd = hbeta->Integral();
  //    double wr = hbeta->GetBinContent(i)/hbetaUnscaled->GetBinContent(i);
  //    double br = hbetaUnscaled->GetBinContent(i);
  //    double err = TMath::Power(wr* (1 - pvalue) / pd, 2) * br;
  //    err_sqrd = err_sqrd + err;
  //  }
  double final_error = TMath::Sqrt(err_sqrd);
  return final_error;
}
