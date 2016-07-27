
//
// ------------------------------------------------------------------------------------------------           
//
// Read and plot the output from CalcPvalue_mh.C
// Compute plot the deltachisq distributions and value from data
// of Normal Hierarchy and Inverted Hierarchy in the same plot
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

void DrawPvalues_mh_NHandIH
(
 //////////////
 //          //
 // Standard //
 //          //
 //////////////
 
 // Run 1-6
 // char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/Pvalue_NH_10kToys_Run16.root",
 // char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/Pvalue_IH_10kToys_Run16.root"

 // Asimov A - dcp=-1.601
 // char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/Pvalue_NH_10kToys_Asimov_A.root",
 // char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/Pvalue_IH_10kToys_Asimov_A.root"

 // // Asimov A - dcp=0
 // char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/Pvalue_NH_10kToys_Asimov_A_dcp0.root",
 // char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/Pvalue_IH_10kToys_Asimov_A_dcp0.root"
 

 //////////////////////////
 //                     //
 // Standard - 50k toys //
 //                     //
 /////////////////////////

 // Run 1-6
 // char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/50k/Pvalue_NH_100kToys_Run16_Total.root",
 // char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/50k/Pvalue_IH_100kToys_Run16_Total.root"

 // Asimov A
 //char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/50k/Pvalue_NH_100kToys_Asimov_A.root",
 //char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/50k/Pvalue_IH_100kToys_Asimov_A.root"

 // // Asimov B
 // char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/50k/Pvalue_NH_100kToys_Asimov_B.root",
 // char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/50k/Pvalue_IH_100kToys_Asimov_B.root"

 

 /////////////////
 //             //
 // Fix deltaCP //
 //             //
 /////////////////
 
 // Run 1-6 - fix dcp=-pi/2
 //char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_NH_10kToys_Run16_fix_dcp_m_pi_over_2.root",
 //char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_IH_10kToys_Run16_fix_dcp_m_pi_over_2.root"
 
 // Run 1-6 - fix dcp=0
 //char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_NH_10kToys_Run16_fix_dcp_0.root",
 //char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_IH_10kToys_Run16_fix_dcp_0.root"

 // Run 1-6 - fix dcp=+pi/2
 //char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_NH_10kToys_Run16_fix_dcp_p_pi_over_2.root",
 //char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_IH_10kToys_Run16_fix_dcp_p_pi_over_2.root"


 // Asimov A - fix dcp=-pi/2
 //char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_NH_10kToys_Asimov_A_fix_dcp_m_pi_over_2.root",
 //char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_IH_10kToys_Asimov_A_fix_dcp_m_pi_over_2.root"
 
 // Asimov A - fix dcp=0
 //char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_NH_10kToys_Asimov_A_fix_dcp_0.root",
 //char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_IH_10kToys_Asimov_A_fix_dcp_0.root"

 // Asimov A - fix dcp=+pi/2
 //char *infilename_nh = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_NH_10kToys_Asimov_A_fix_dcp_p_pi_over_2.root",
 //char *infilename_ih = "../TECH_NOTE_14_3_16/MASSHIERARCHY/PVALUE/fixdcp/Pvalue_IH_10kToys_Asimov_A_fix_dcp_p_pi_over_2.root"


 )
{  
  char name[300];

  // Read file of Normal Hierarchy 

  TFile *filein_NH = new TFile(infilename_nh,"READ");
  
  TH1D  *hdeltachisq_Toys_nh_ih_post_NH = (TH1D*)  filein_NH->Get("hdeltachisq_Toys_nh_ih_post");
  TLine *lPvalue_nh_ih_post_NH          = (TLine*) filein_NH->Get("Pvalue_nh_ih_posterior"); 
  TLine *lPvalue_nh_ih_post_err_NH      = (TLine*) filein_NH->Get("lPvalue_nh_ih_posterior_err");
  TLine *lZscore_nh_ih_post_NH          = (TLine*) filein_NH->Get("Zscore_nh_ih_posterior"); 
  TLine *lZscore_nh_ih_post_perr_NH     = (TLine*) filein_NH->Get("lzscore_nh_ih_posterior_perr");
  TLine *lZscore_nh_ih_post_merr_NH     = (TLine*) filein_NH->Get("lzscore_nh_ih_posterior_merr");

  double Pvalue_nh_ih_post_NH           = lPvalue_nh_ih_post_NH->GetX1();
  double Pvalue_nh_ih_post_err_NH       = lPvalue_nh_ih_post_err_NH->GetX1();
  double Zscore_nh_ih_post_NH           = lZscore_nh_ih_post_NH->GetX1();
  double Zscore_nh_ih_post_perr_NH      = lZscore_nh_ih_post_perr_NH->GetX1();
  double Zscore_nh_ih_post_merr_NH      = lZscore_nh_ih_post_merr_NH->GetX1();

  cout << "Pvalue posterior (NH) = " << Pvalue_nh_ih_post_NH << " +- " << Pvalue_nh_ih_post_err_NH << endl;
  cout << "Z-score posterior (NH) = " << Zscore_nh_ih_post_NH << " + " << Zscore_nh_ih_post_perr_NH << " - " << Zscore_nh_ih_post_merr_NH << endl;

  TH1D  *hdeltachisq_Toys_nh_ih_NH = (TH1D*)  filein_NH->Get("hdeltachisq_Toys_nh_ih");
  TLine *ldeltachisq_Data_nh_ih_NH = (TLine*) filein_NH->Get("ldeltachisq_Data_nh_ih");
  TLine *lPvalue_nh_ih_NH          = (TLine*) filein_NH->Get("Pvalue_nh_ih"); 
  TLine *lZscore_nh_ih_NH          = (TLine*) filein_NH->Get("Zscore_nh_ih"); 
  TH1D  *hdeltachisq_Toys_fc_NH = (TH1D*)  filein_NH->Get("hdeltachisq_Toys_fc");
  TLine *ldeltachisq_Data_fc_NH = (TLine*) filein_NH->Get("ldeltachisq_Data_fc");
  TLine *lPvalue_fc_NH          = (TLine*) filein_NH->Get("Pvalue_fc"); 
  TLine *lZscore_fc_NH          = (TLine*) filein_NH->Get("Zscore_fc"); 

  double deltachisq_nh_ih_data_NH = ldeltachisq_Data_nh_ih_NH->GetX1();
  double Pvalue_nh_ih_NH          = lPvalue_nh_ih_NH->GetX1();
  double Zscore_nh_ih_NH          = lZscore_nh_ih_NH->GetX1();
  double deltachisq_fc_data_NH = ldeltachisq_Data_fc_NH->GetX1();
  double Pvalue_fc_NH          = lPvalue_fc_NH->GetX1();
  double Zscore_fc_NH          = lZscore_fc_NH->GetX1();


  // Read file of Inverted Hierarchy 

  TFile *filein_IH = new TFile(infilename_ih,"READ");
  
  TH1D  *hdeltachisq_Toys_nh_ih_post_IH = (TH1D*)  filein_IH->Get("hdeltachisq_Toys_nh_ih_post");
  TLine *lPvalue_nh_ih_post_IH          = (TLine*) filein_IH->Get("Pvalue_nh_ih_posterior"); 
  TLine *lPvalue_nh_ih_post_err_IH      = (TLine*) filein_IH->Get("lPvalue_nh_ih_posterior_err");
  TLine *lZscore_nh_ih_post_IH          = (TLine*) filein_IH->Get("Zscore_nh_ih_posterior"); 
  TLine *lZscore_nh_ih_post_perr_IH     = (TLine*) filein_IH->Get("lzscore_nh_ih_posterior_perr");
  TLine *lZscore_nh_ih_post_merr_IH     = (TLine*) filein_IH->Get("lzscore_nh_ih_posterior_merr");

  double Pvalue_nh_ih_post_IH           = lPvalue_nh_ih_post_IH->GetX1();
  double Pvalue_nh_ih_post_err_IH       = lPvalue_nh_ih_post_err_IH->GetX1();
  double Zscore_nh_ih_post_IH           = lZscore_nh_ih_post_IH->GetX1();
  double Zscore_nh_ih_post_perr_IH      = lZscore_nh_ih_post_perr_IH->GetX1();
  double Zscore_nh_ih_post_merr_IH      = lZscore_nh_ih_post_merr_IH->GetX1();

  cout << "Pvalue posterior (IH) = " << Pvalue_nh_ih_post_IH << " +- " << Pvalue_nh_ih_post_err_IH << endl;
  cout << "Z-score posterior (IH) = " << Zscore_nh_ih_post_IH << " + " << Zscore_nh_ih_post_perr_IH << " - " << Zscore_nh_ih_post_merr_IH << endl;

  TH1D  *hdeltachisq_Toys_nh_ih_IH = (TH1D*)  filein_IH->Get("hdeltachisq_Toys_nh_ih");
  TLine *ldeltachisq_Data_nh_ih_IH = (TLine*) filein_IH->Get("ldeltachisq_Data_nh_ih");
  TLine *lPvalue_nh_ih_IH          = (TLine*) filein_IH->Get("Pvalue_nh_ih"); 
  TLine *lZscore_nh_ih_IH          = (TLine*) filein_IH->Get("Zscore_nh_ih"); 
  TH1D  *hdeltachisq_Toys_fc_IH = (TH1D*)  filein_IH->Get("hdeltachisq_Toys_fc");
  TLine *ldeltachisq_Data_fc_IH = (TLine*) filein_IH->Get("ldeltachisq_Data_fc");
  TLine *lPvalue_fc_IH          = (TLine*) filein_IH->Get("Pvalue_fc"); 
  TLine *lZscore_fc_IH          = (TLine*) filein_IH->Get("Zscore_fc"); 

  double deltachisq_nh_ih_data_IH = ldeltachisq_Data_nh_ih_IH->GetX1();
  double Pvalue_nh_ih_IH          = lPvalue_nh_ih_IH->GetX1();
  double Zscore_nh_ih_IH          = lZscore_nh_ih_IH->GetX1();
  double deltachisq_fc_data_IH = ldeltachisq_Data_fc_IH->GetX1();
  double Pvalue_fc_IH          = lPvalue_fc_IH->GetX1();
  double Zscore_fc_IH          = lZscore_fc_IH->GetX1();

  if(deltachisq_nh_ih_data_NH != deltachisq_nh_ih_data_IH){
    cerr << endl;
    cerr << "Two different data fit were performed!!!" << endl;
    cerr << "deltachisq_nh_ih_data_NH = " << deltachisq_nh_ih_data_NH << endl;
    cerr << "deltachisq_nh_ih_data_IH = " << deltachisq_nh_ih_data_IH << endl;
    cerr << endl;
    exit(1);
  }


  // Draw p-value for chisq(NH)-chisq(IH)

  hdeltachisq_Toys_nh_ih_NH->SetLineWidth(2);
  hdeltachisq_Toys_nh_ih_NH->SetLineColor(kRed);
  hdeltachisq_Toys_nh_ih_IH->SetLineWidth(2);
  hdeltachisq_Toys_nh_ih_IH->SetLineColor(kBlack);
  
  ldeltachisq_Data_nh_ih_NH->SetLineWidth(2);
  ldeltachisq_Data_nh_ih_NH->SetLineColor(kGreen);
  ldeltachisq_Data_nh_ih_IH->SetLineWidth(2);
  ldeltachisq_Data_nh_ih_IH->SetLineColor(kGreen);

  TCanvas *PValue_NH_IH = new TCanvas("PValue_NH_IH","PValue_NH_IH");
  hdeltachisq_Toys_nh_ih_NH->Draw();
  hdeltachisq_Toys_nh_ih_NH->SetMaximum( hdeltachisq_Toys_nh_ih_NH->GetMaximum()*1.3 );
  hdeltachisq_Toys_nh_ih_NH->GetXaxis()->SetRangeUser(-5,8);
  hdeltachisq_Toys_nh_ih_IH->Draw("same");

  ldeltachisq_Data_nh_ih_IH->Draw();

  TLegend *legPvalue_nh_ih = new TLegend(0.5,0.6,0.8,0.9);
  legPvalue_nh_ih->SetFillColor(0);
  legPvalue_nh_ih->SetLineColor(0);
  legPvalue_nh_ih->SetTextFont(132);
  legPvalue_nh_ih->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_nh_ih_data_NH); 
  legPvalue_nh_ih->AddEntry(ldeltachisq_Data_nh_ih_NH,name,"l");
  sprintf(name,"Normal #Delta #chi^{2}_{MC}"); //,deltachisq_nh_ih_data_NH); 
  legPvalue_nh_ih->AddEntry(hdeltachisq_Toys_nh_ih_NH,name,"l");
  sprintf(name,"Inverted #Delta #chi^{2}_{MC}"); //,deltachisq_nh_ih_data_IH); 
  legPvalue_nh_ih->AddEntry(hdeltachisq_Toys_nh_ih_IH,name,"l");
  legPvalue_nh_ih->Draw();  

  int binlat = hdeltachisq_Toys_nh_ih_NH->GetXaxis()->FindBin(0.);
  TLatex *lat_nh_ih_NH = new TLatex(hdeltachisq_Toys_nh_ih_NH->GetXaxis()->GetBinLowEdge(binlat),hdeltachisq_Toys_nh_ih_NH->GetMaximum()*0.5, 
				    TString::Format("#splitline{Normal Hierarchy}{#splitline{P-value = %1.3f}{Z-score = %1.3f#sigma}}",Pvalue_nh_ih_NH,Zscore_nh_ih_NH).Data() );
  TLatex *lat_nh_ih_IH = new TLatex(hdeltachisq_Toys_nh_ih_IH->GetXaxis()->GetBinLowEdge(binlat),hdeltachisq_Toys_nh_ih_IH->GetMaximum()*0.5, 
				    TString::Format("#splitline{Inverted Hierarchy}{#splitline{P-value = %1.3f}{Z-score = %1.3f#sigma}}",Pvalue_nh_ih_IH,Zscore_nh_ih_IH).Data() );
  
  lat_nh_ih_NH->Draw("same");
  lat_nh_ih_NH->SetTextSize(0.06);
  lat_nh_ih_IH->Draw("same");
  lat_nh_ih_IH->SetTextSize(0.06);
  
  gPad->RedrawAxis();






  // Draw p-value for chisq(NH)-chisq(IH) (Posterior)

  hdeltachisq_Toys_nh_ih_post_NH->SetLineWidth(2);
  hdeltachisq_Toys_nh_ih_post_NH->SetLineColor(kRed);
  hdeltachisq_Toys_nh_ih_post_IH->SetLineWidth(2);
  hdeltachisq_Toys_nh_ih_post_IH->SetLineColor(kBlack);
  
  TCanvas *PValue_NH_IH_post = new TCanvas("PValue_NH_IH_post","PValue_NH_IH_post");
  hdeltachisq_Toys_nh_ih_post_NH->Draw();
  hdeltachisq_Toys_nh_ih_post_NH->SetMaximum( hdeltachisq_Toys_nh_ih_post_NH->GetMaximum()*1.3 );
  hdeltachisq_Toys_nh_ih_post_NH->GetXaxis()->SetRangeUser(-5,8);
  hdeltachisq_Toys_nh_ih_post_IH->Draw("same");
  
  ldeltachisq_Data_nh_ih_IH->Draw();

  TLegend *legPvalue_nh_ih_post = new TLegend(0.5,0.6,0.8,0.9);
  legPvalue_nh_ih_post->SetFillColor(0);
  legPvalue_nh_ih_post->SetLineColor(0);
  legPvalue_nh_ih_post->SetTextFont(132);
  legPvalue_nh_ih_post->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_nh_ih_data_NH); 
  legPvalue_nh_ih_post->AddEntry(ldeltachisq_Data_nh_ih_NH,name,"l");
  sprintf(name,"Normal #Delta #chi^{2}_{MC}"); //,deltachisq_nh_ih_data_post_NH); 
  legPvalue_nh_ih_post->AddEntry(hdeltachisq_Toys_nh_ih_post_NH,name,"l");
  sprintf(name,"Inverted #Delta #chi^{2}_{MC}"); //,deltachisq_nh_ih_data_post_IH); 
  legPvalue_nh_ih_post->AddEntry(hdeltachisq_Toys_nh_ih_post_IH,name,"l");
  legPvalue_nh_ih_post->Draw();  

  int binlat = hdeltachisq_Toys_nh_ih_post_NH->GetXaxis()->FindBin(0.);
  TLatex *lat_nh_ih_post_NH = new TLatex(hdeltachisq_Toys_nh_ih_post_NH->GetXaxis()->GetBinLowEdge(binlat),hdeltachisq_Toys_nh_ih_post_NH->GetMaximum()*0.5, 
					 TString::Format("#splitline{Normal Hierarchy}{#splitline{P-value = %1.3f #pm %1.3f}{Z-score = %1.3f +%1.3f -%1.3f #sigma}}",
							 Pvalue_nh_ih_post_NH,Pvalue_nh_ih_post_err_NH,Zscore_nh_ih_post_NH,Zscore_nh_ih_post_perr_NH,Zscore_nh_ih_post_merr_NH).Data());

  TLatex *lat_nh_ih_post_IH = new TLatex(hdeltachisq_Toys_nh_ih_post_IH->GetXaxis()->GetBinLowEdge(binlat),hdeltachisq_Toys_nh_ih_post_IH->GetMaximum()*0.5, 
					 TString::Format("#splitline{Inverted Hierarchy}{#splitline{P-value = %1.3f #pm %1.3f}{Z-score = %1.3f +%1.3f - %1.3f #sigma}}",
							 Pvalue_nh_ih_post_IH,Pvalue_nh_ih_post_err_IH,Zscore_nh_ih_post_IH,Zscore_nh_ih_post_perr_IH,Zscore_nh_ih_post_merr_IH).Data());
    
  lat_nh_ih_post_NH->Draw("same");
  lat_nh_ih_post_NH->SetTextSize(0.06);
  lat_nh_ih_post_IH->Draw("same");
  lat_nh_ih_post_IH->SetTextSize(0.06);
  
  gPad->RedrawAxis();

  
  // Draw p-value for Feldman-Cousins - Normal Hierarchy

  hdeltachisq_Toys_fc_NH->SetLineWidth(2);
  hdeltachisq_Toys_fc_NH->SetLineColor(kBlack);
  
  ldeltachisq_Data_fc_NH->SetLineWidth(2);
  ldeltachisq_Data_fc_NH->SetLineColor(kRed);

  TCanvas *PValue_FC_NH = new TCanvas("PValue_FC_NH","PValue_FC_NH");
  hdeltachisq_Toys_fc_NH->Draw();
  hdeltachisq_Toys_fc_NH->SetMaximum( hdeltachisq_Toys_fc_NH->GetMaximum()*1.3 );
  hdeltachisq_Toys_fc_NH->GetXaxis()->SetRangeUser(-1,5);

  TLegend *legPvalue_fc_NH = new TLegend(0.5,0.6,0.8,0.9);
  legPvalue_fc_NH->SetFillColor(0);
  legPvalue_fc_NH->SetLineColor(0);
  legPvalue_fc_NH->SetTextFont(132);
  legPvalue_fc_NH->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_fc_data_NH); 
  legPvalue_fc_NH->AddEntry(ldeltachisq_Data_fc_NH,name,"l");
  sprintf(name,"Normal #Delta #chi^{2}_{MC}",deltachisq_fc_data_NH); 
  legPvalue_fc_NH->AddEntry(hdeltachisq_Toys_fc_NH,name,"l");
  legPvalue_fc_NH->Draw();  

  ldeltachisq_Data_fc_NH->Draw();
  
  int binlat = hdeltachisq_Toys_fc_NH->GetXaxis()->FindBin(0.);
  TLatex *lat_fc_NH = new TLatex(hdeltachisq_Toys_fc_NH->GetXaxis()->GetBinLowEdge(binlat),hdeltachisq_Toys_fc_NH->GetMaximum()*0.5, 
				    TString::Format("#splitline{Normal Hierarchy}{#splitline{P-value = %1.3f}{Z-score = %1.3f#sigma}}",Pvalue_fc_NH,Zscore_fc_NH).Data() );
   
  lat_fc_NH->Draw("same");
  lat_fc_NH->SetTextSize(0.06);
 
  

  // Draw p-value for Feldman-Cousins - Inverted Hierarchy

  hdeltachisq_Toys_fc_IH->SetLineWidth(2);
  hdeltachisq_Toys_fc_IH->SetLineColor(kBlack);
  
  ldeltachisq_Data_fc_IH->SetLineWidth(2);
  ldeltachisq_Data_fc_IH->SetLineColor(kRed);

  TCanvas *PValue_FC_IH = new TCanvas("PValue_FC_IH","PValue_FC_IH");
  hdeltachisq_Toys_fc_IH->Draw();
  hdeltachisq_Toys_fc_IH->SetMaximum( hdeltachisq_Toys_fc_IH->GetMaximum()*1.3 );
  hdeltachisq_Toys_fc_IH->GetXaxis()->SetRangeUser(-1,5);

  TLegend *legPvalue_fc_IH = new TLegend(0.5,0.6,0.8,0.9);
  legPvalue_fc_IH->SetFillColor(0);
  legPvalue_fc_IH->SetLineColor(0);
  legPvalue_fc_IH->SetTextFont(132);
  legPvalue_fc_IH->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_fc_data_IH); 
  legPvalue_fc_IH->AddEntry(ldeltachisq_Data_fc_IH,name,"l");
  sprintf(name,"Inverted #Delta #chi^{2}_{MC}",deltachisq_fc_data_IH); 
  legPvalue_fc_IH->AddEntry(hdeltachisq_Toys_fc_IH,name,"l");
  legPvalue_fc_IH->Draw();  

  ldeltachisq_Data_fc_IH->Draw();
  
  int binlat = hdeltachisq_Toys_fc_IH->GetXaxis()->FindBin(0.);
  TLatex *lat_fc_IH = new TLatex(hdeltachisq_Toys_fc_IH->GetXaxis()->GetBinLowEdge(binlat),hdeltachisq_Toys_fc_IH->GetMaximum()*0.5, 
				    TString::Format("#splitline{Inverted Hierarchy}{#splitline{P-value = %1.3f}{Z-score = %1.3f#sigma}}",Pvalue_fc_IH,Zscore_fc_IH).Data() );
  
  lat_fc_IH->Draw("same");
  lat_fc_IH->SetTextSize(0.06);

  return;









}
