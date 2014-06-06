void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style);

void draw_sensitivity_plots_toys(
			    bool reactor = true, 
			    bool theta2theta = false,
			    const char * name_nh = "../inputs/toyexptsensi_PROVA_3toys",
			    const char * name_ih = "../inputs/toyexptsensi_PROVA"
			    )
{
 
  using namespace std;
  using namespace std::string;
  using namespace valor;
  using namespace valor::conv;
   
  const char * filename = Form("%s.root",name_nh);
  FitParamList *  fparm = 0;
  FitParamList *  fparm_deltacp = 0;
  PhysParamList * iparm = new PhysParamList;
 
  //Set up bins and histograms
  double mixing13_min = 0;
  double mixing13_max = 0;
  double mixing13_step = 0;

  if(!theta2theta){
    if(reactor){
      mixing13_min = 0.015;
      mixing13_max = 0.035;
      mixing13_step = 0.001;
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.1;
      mixing13_step = 0.005;
    }
  }
  else{
    if(reactor){
      mixing13_min = 0.06;
      mixing13_max = 0.14;
      mixing13_step = 0.004;
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.4;
      mixing13_step = 0.02;
    }
  }

  double mixing23_min = 0.3;
  double mixing23_max = 0.7;
  double mixing23_step = 0.01;
  
  double dmsq_min = 2.0;//E-3;
  double dmsq_max = 3.0;//E-3;
  double dmsq_step = 0.05;//E-3;

  double deltacp_min = -conv::kpi;
  double deltacp_max = conv::kpi;
  double deltacp_step = (2*conv::kpi)/50.;

  ///////////////////////////////// TO CHANGE ///////////////////////
  
  int nbins_mixing13 = 21;
  int nbins_mixing23 = 41;
  int nbins_dmsq = 21;
  const int nbins_deltacp = 51;

  //////////////////////////////////////////////////////////////////
  
  // 1D vs deltacp

  const char * filename_deltacp = Form("%s.root",name_nh);
  TChain ch_deltacp("fitnt");
  ch_deltacp.Add(filename_deltacp);
  //"best_fit_params" is already the minimum after the two fits one of each starting in a different octant
  ch_deltacp.SetBranchAddress("best_fit_params", &fparm_deltacp); 

  const int nentries_deltacp = ch_deltacp.GetEntries();
  cout << "nentries deltacp = " << nentries_deltacp << endl;
  
  const int ntoys = nentries_deltacp / nbins_deltacp;
  cout << "ntoys = " << ntoys << endl;

  //Double_t* bestfit_points_x = new Double_t[nbins_deltacp];
  //Double_t* bestfit_points_y = new Double_t[nbins_deltacp];

  Double_t** bestfit_points_x = new Double_t*[ntoys]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys; itoy++){
    bestfit_points_x[itoy] = new Double_t[nbins_deltacp];
  }
  Double_t** bestfit_points_y = new Double_t*[ntoys]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys; itoy++){
    bestfit_points_y[itoy] = new Double_t[nbins_deltacp];
  }

  //  TH1D * h_chisq_deltacp = new TH1D("h_deltacp","h_deltacp",nbins_deltacp,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  //  TH1D * h_mix13_deltacp = new TH1D("h_mix13_deltacp","h_mix13_deltacp",nbins_deltacp,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  TH1D ** h_chisq_deltacp = new TH1D*[ntoys];
  for(int itoy=0; itoy<ntoys; itoy++){
    char histoname[200];
    sprintf(histoname,"h_deltacp_Toy%d",itoy);
    h_chisq_deltacp[itoy] = new TH1D(histoname,histoname,nbins_deltacp,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  }  
  TH1D ** h_mix13_deltacp= new TH1D*[ntoys];
  for(int itoy=0; itoy<ntoys; itoy++){
    char histoname[200];
    sprintf(histoname,"h_mix13_deltacp_Toy%d",itoy);
    h_mix13_deltacp[itoy] = new TH1D(histoname,histoname,nbins_deltacp,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  }  

  // for(int i = 0; i < nentries_deltacp; i++){ 

  int iEvtToy = 0;
  for(int i = 0; i < nentries_deltacp; i++){
   
    iEvtToy = i / nbins_deltacp;
    //if( !(i%nbins_deltacp) ) cout << "# Toy " << iEvtToy << endl;
    
    ch_deltacp.GetEntry(i);
    int migradminstat = fparm_deltacp->GetStatusCode("migrad_min");
    if (migradminstat == 4) cout << "Fit failed!" << endl;
    double mix13 = fparm_deltacp->BestFitValue(conv::kPhypMixing13);
    double mix23 = fparm_deltacp->BestFitValue(conv::kPhypMixing23);
    double dmsq = fparm_deltacp->BestFitValue(conv::kPhypDm23sq);
    double deltacp = fparm_deltacp->BestFitValue(conv::kPhypCPV);
    double chisq = fparm_deltacp->GetBestFitStatistic();
    
    int bin_idx = h_chisq_deltacp[iEvtToy]->FindBin(deltacp);
    //double bincont = h_chisq_deltacp[iEvtToy]->GetBinContent(bin_idx) / (nentries_deltacp / nbins_deltacp);
    double bincont = h_chisq_deltacp[iEvtToy]->GetBinContent(bin_idx);
    if(bincont!=0){
      if(bincont>chisq){
	if(chisq!=0){
	  h_chisq_deltacp[iEvtToy]->SetBinContent(bin_idx,chisq+bincont);
	  bestfit_points_x[iEvtToy][bin_idx-1]  = mix13;//sinsq2theta13;
	  bestfit_points_y[iEvtToy][bin_idx-1]  = deltacp;
	}
      }
    } 
    else{ 
      if(chisq!=0){
	h_chisq_deltacp[iEvtToy]->SetBinContent(bin_idx,chisq+bincont);
	bestfit_points_x[iEvtToy][bin_idx-1] = mix13;//sinsq2theta13;
	bestfit_points_y[iEvtToy][bin_idx-1] = deltacp;
      }
    }
  }

  //TGraph* gbestfit = new TGraph( nbins_deltacp, bestfit_points_x ,bestfit_points_y ); // dcp Vs mix13  
  TGraph ** gbestfit_Toy = new TGraph*[ntoys]; 

  for(int itoy=0;itoy<ntoys;itoy++){
    char graphname[200];
    sprintf(graphname,"gbestfit_Toy%d",itoy);
    gbestfit_Toy[itoy] = new TGraph( nbins_deltacp, bestfit_points_x[itoy] ,bestfit_points_y[itoy] ); // dcp Vs mix13  
    gbestfit_Toy[itoy]->SetName(graphname);
    //cout << graphname << endl;
  }
  
  cout << endl;
  
  // Find the minimum chi2
  double min_chisq[ntoys];
  for(int itoy=0;itoy<ntoys;itoy++){
    min_chisq[itoy] = 999.;
    cout << "# Toy " << itoy << endl;    
    for(int i = 1; i <= nbins_deltacp; i++){
      double bin = h_chisq_deltacp[itoy]->GetBinContent(i);
      if (bin < min_chisq[itoy]) min_chisq[itoy] = bin;
      //cout << bin << ", ";
    }
    cout << "For deltaCP, min chisq  = " << min_chisq[itoy] << endl;
    cout << endl;
  }

  

  // Shift the chi2 wrt the minimum
  for(int itoy=0;itoy<ntoys;itoy++){
    for(int i = 1; i <= nbins_deltacp; i++){
      double bin = h_chisq_deltacp[itoy]->GetBinContent(i);
      h_chisq_deltacp[itoy]->SetBinContent(i,bin-min_chisq[itoy]); // scanned global minimum
      //h_chisq_deltacp[itoy]->SetBinContent(i,bin-chisq_bf); // fiited global minimum
    }
  }

  TCanvas *cdeltacp = new TCanvas("cdeltacp","cdeltacp");
  for(int itoy=0;itoy<ntoys;itoy++){
    h_chisq_deltacp[itoy]->GetXaxis()->SetTitle("#delta_{CP}");
    //h_chisq_deltacp[itoy]->GetYaxis()->SetTitle("sin^{2}#theta_{23}");
    h_chisq_deltacp[itoy]->GetYaxis()->SetTitle("#Delta#chi^{2}");
    //h_chisq_deltacp[itoy]->GetYaxis()->SetTitle("#chi^{2}");
    //h_chisq_deltacp[itoy]->GetYaxis()->SetRangeUser(55,135);
    h_chisq_deltacp[itoy]->GetYaxis()->SetRangeUser(0,4);
    h_chisq_deltacp[itoy]->SetLineWidth(2);
    if(itoy==0) h_chisq_deltacp[itoy]->Draw("c");  
    else h_chisq_deltacp[itoy]->Draw("c same");  
  }

  TCanvas *cbestfit = new TCanvas("cbestfit","cbestfit");
  for(int itoy=0;itoy<ntoys;itoy++){  
    gbestfit_Toy[itoy]->GetXaxis()->SetTitle("#delta_{CP}");
    gbestfit_Toy[itoy]->GetYaxis()->SetTitle("#Delta#chi^{2}");
    gbestfit_Toy[itoy]->GetYaxis()->SetRangeUser(-3.2,3.2);
    gbestfit_Toy[itoy]->SetLineWidth(2);
    gbestfit_Toy[itoy]->SetLineColor(kBlue);
    if(itoy==0) gbestfit_Toy[itoy]->Draw("ac");
    else gbestfit_Toy[itoy]->Draw("c same");
  }












 //  const char * filename_sensitivity = "data_sensitivity_joint_final_ih_chisq_vs_deltacp.root";
 //  TChain ch_sensitivity("fitnt");
 //  ch_sensitivity.Add(filename_sensitivity);

 //  FitParamList *  fparm_sensitivity = 0;
 //  ch_sensitivity.SetBranchAddress("best_fit_params", &fparm_sensitivity); 
 //  const int nentries_sensitivity = ch_sensitivity.GetEntries();
 //  const int ntoys = (const int) nentries_sensitivity/nbins_deltacp;
 //  cout << "nentries_sensitivity = " << nentries_sensitivity << " and ntoys = " << ntoys << endl;
  

 //  //  double chisq_sensitivity[nbins_deltacp][ntoys];
 //  //double x_sensitivity[nbins_deltacp][ntoys];
 //  double chisq_sensitivity[51][487];
 //  double x_sensitivity[51][487];
 //  int j = 0;
 //  int k = 0;

 //  for(int i = 0; i < nentries_sensitivity; i++) 
 //    { 
 //      ch_sensitivity.GetEntry(i);
 //      double deltacp = fparm_sensitivity->BestFitValue(conv::kPhypCPV);
 //      double chisq = fparm_sensitivity->GetBestFitStatistic();
 //      //      cout << " deltacp = " << deltacp << " and chisq = " << chisq << endl;
 //      chisq_sensitivity[j][k] = chisq;
 //      x_sensitivity[j][k] = deltacp;
 //      //      cout << "chisq_sensitivity["<<j<<"]["<<k<<"] = "<< chisq_sensitivity[j][k]<<endl;
 //      j++;
 //      if(j == nbins_deltacp){
 // 	j = 0;
 // 	k ++; 
 //      }
 //    }

 //  TGraph* gchisq_toy[nbins_deltacp];
  
  
 //  for(int i = 0; i < nbins_deltacp; i++) 
 //    {
 //      gchisq_toy[i] = new TGraph( ntoys, x_sensitivity[i], chisq_sensitivity[i]);
 //      gchisq_toy[i]->SetMarkerStyle(1);
 //      gchisq_toy[i]->Draw("P");
 //    }

 //  double chisq_average[51];
  
 //  for(int j = 0; j < nbins_deltacp; j++)
 //    { 
 //      double average = 0;
 //      for(int i = 0; i < ntoys; i++) 
 // 	{ 
 // 	  average + = chisq_sensitivity[j][i];
 // 	}
 //      average /= ntoys;
 //      chisq_average[j] = average;
 //    }
  


 //  double average_min = 9999;
 //  for(int j = 0; j < nbins_deltacp; j++)
 //    { 
 //      if(chisq_average[j]<average_min) 
 //  	average_min = chisq_average[j];
 //    }
 //  for(int j = 0; j < nbins_deltacp; j++)
 //    { 
 //      chisq_average[j] -= average_min;
      
 //      cout << chisq_average[j] << endl;
 //    }
  





  

 //  TGraph* gchisq_average = new TGraph( 51, bestfit_points_y, chisq_average);
 //  gchisq_average->SetLineColor(kBlue);
 //  gchisq_average->SetLineWidth(2);
 //  gchisq_average->Draw("L");

 //  TLegend * legend = new TLegend(0.7, 0.6, 0.9, 0.9);
 //  legend->SetFillColor(0); 
 //  legend->SetBorderSize(0);
 //  legend->AddEntry(h_chisq_deltacp,"T2K Run1+2+3+4 dataset fit (IH)","L");
 //  legend->AddEntry(gchisq_average,"T2K Run1+2+3+4 sensitivity (IH)","L");
 //  // legend->Draw("same");

 //  //to calculate the 68% CL
  
 //  double chisq_above[51][487];
 //  double n_above[51];
 //  double chisq_below[51][487];
 //  double n_below[51];

 // for(int j = 0; j < nbins_deltacp; j++)
 //    { 
 //      int i_above = 0;
 //      int i_below = 0;
 //      for(int i = 0; i < ntoys; i++) 
 // 	{ 
 // 	  if(chisq_sensitivity[j][i]>chisq_average[j]){
 // 	    chisq_above[j][i_above] = chisq_sensitivity[j][i];
 // 	    i_above++;
 // 	  }
 // 	  else{
 // 	    chisq_below[j][i_below] = chisq_sensitivity[j][i];
 // 	    i_below++;
 // 	  }
 // 	}
 //      n_above[j] = i_above;
 //      n_below[j] = i_below;
 //      //      cout << " n_above = " << n_above[j] << " and n_below = " << n_below[j] << endl;
 //    }

 // double sigma_above[51];
 // double sigma_below[51];

 // for(int j = 0; j < nbins_deltacp; j++)
 //   { 
 //     Int_t idx[999999];
 //     TMath::Sort(n_below[j], chisq_below[j], idx, true); //the one above, sorted from largest to smallest
 //     int cut = (int)n_below[j]*0.68;
 //     sigma_below[j] = chisq_below[j][idx[cut-1]];

 //     Int_t idx2[999999];
 //     TMath::Sort(n_above[j], chisq_above[j], idx2, false); //the one below, sorted from smallest to largest
 //     int cut2 = (int)n_above[j]*0.68;
 //     sigma_above[j] = chisq_above[j][idx2[cut2-1]];
 //     //    cout << "cut = " << cut << " and = " << chisq_below[0][idx[cut-1]];
 //   }

 // TGraph* gchisq_above = new TGraph( 51, bestfit_points_y, sigma_above);
 // gchisq_above->SetLineColor(kBlue);
 // gchisq_above->SetLineWidth(2);
 // gchisq_above->Draw("L");

 // TGraph* gchisq_below = new TGraph( 51, bestfit_points_y, sigma_below);
 // gchisq_below->SetLineColor(kBlue);
 // gchisq_below->SetLineWidth(2);
 // gchisq_below->Draw("L");

 // TGraph *grshade = new TGraph(2*51);
 //  for(int j = 0; j < nbins_deltacp; j++)
 //   { 
 //     grshade->SetPoint(j,bestfit_points_y[j],sigma_above[j]);
 //     grshade->SetPoint(51+j,bestfit_points_y[50-j],sigma_below[50-j]);
 //   }
 //  grshade->SetFillStyle(3013);
 //  grshade->SetFillColor(kAzure+10);
 //  grshade->Draw("F");
  

 // /*   
 // for(int i = 0; i < ntoys; i++) 
 //   { 
 //     cout << chisq_below[0][idx[i]] << "  " ;
 //   }
 // */

   
 //  const char * filename_out = Form("%s_sensitivity.root",name);
 //  TFile * f_out = new TFile(filename_out, "recreate");
 //  h_chisq_deltacp->Write();
 
}


void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style)
{
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
  contour->SetContour(1);//number of contour levels                                                                                                                                    
  contour->SetContourLevel(0,0.0);//SetContourLevel(int level, double value)                                                                                                           
  contour->SetStats(0);

} 
