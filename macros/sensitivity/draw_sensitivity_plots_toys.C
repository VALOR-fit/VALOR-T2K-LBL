//
// Sensitivity study w/ 1D or 2D plost
//

#define doShiftMin false
#define doNH true
#define doIH false
#define doData false

void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style);

void draw_sensitivity_plots_toys(
			    bool reactor = true, 
			    bool theta2theta = false,
			    //const char * name_nh = "../inputs/toyexptsensi_PROVA",
			    //const char * name_nh = "../inputs/toyexptsensi_PROVA_3toys_New",
			    //const char * name_nh = "../inputs/toyexptsensi_Run4_FHConly_StatOnly",			  
			    const char * name_nh = "../inputs/data_statonly",
			    //const char * name_nh = "../inputs/toyexptsensi_Run1-6_FHCpRHC_Syst_NH_dcp_0",
			    const char * name_ih = "../inputs/toyexptsensi_Run1-4_FHC_Syst_NH_dcp_0"
			    )
{
  using namespace std;
  using namespace std::string;
  using namespace valor;
  using namespace valor::conv;

  const char * filename_nh = Form("%s.root",name_nh);
  FitParamList * fparm_nh = 0;
  FitParamList * fparm_deltacp_nh = 0;
  PhysParamList * iparm_nh = new PhysParamList;
  int nbins_param1_nh, nbins_param2_nh;
  int nbins1_nh = 0;
  int nbins2_nh = 0;

  const char * filename_ih = Form("%s.root",name_ih);
  FitParamList * fparm_ih = 0;
  FitParamList * fparm_deltacp_ih = 0;
  PhysParamList * iparm_ih = new PhysParamList;
  int nbins_param1_ih, nbins_param2_ih;
  int nbins1_ih = 0;
  int nbins2_ih = 0;

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

  // 1D vs deltacp (Normal Hierarchy)

  const char * filename_deltacp_nh = Form("%s.root",name_nh);
  TChain ch_deltacp_nh("fitnt");
  ch_deltacp_nh.Add(filename_deltacp_nh);
  ch_deltacp_nh.SetBranchAddress("best_fit_params", &fparm_deltacp_nh); 
  ch_deltacp_nh.SetBranchAddress("nbins_param1", &nbins_param1_nh); 
  ch_deltacp_nh.SetBranchAddress("nbins_param2", &nbins_param2_nh); 

  const int nentries_deltacp_nh = ch_deltacp_nh.GetEntries();
  cout << "nentries deltacp (NH) = " << nentries_deltacp_nh << endl;
  ch_deltacp_nh.GetEntry(0);
  const int nbins1_nh = nbins_param1_nh;
 
  const int ntoys_nh = nentries_deltacp_nh / nbins1_nh;
  cout << "ntoys (NH) = " << ntoys_nh << endl;


  double mixing23_min = 0.3;
  double mixing23_max = 0.7;
  double mixing23_step = 0.01;
  
  double dmsq_min = 2.0;//E-3;
  double dmsq_max = 3.0;//E-3;
  double dmsq_step = 0.05;//E-3;

  double deltacp_min = -conv::kpi;
  double deltacp_max = conv::kpi;
  double deltacp_step = (2*conv::kpi)/(nentries_deltacp_nh-1);


  Double_t** bestfit_points_nh_x = new Double_t*[ntoys_nh]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    bestfit_points_nh_x[itoy] = new Double_t[nbins1_nh];
  }
  Double_t** bestfit_points_nh_y = new Double_t*[ntoys_nh]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    bestfit_points_nh_y[itoy] = new Double_t[nbins1_nh];
  }

  TH1D * h_chisq_deltacp_Ave_nh = new TH1D("h_chisq_deltacp_Ave_nh","h_chisq_deltacp_Ave_nh",nbins1_nh,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  TH1D ** h_chisq_deltacp_nh = new TH1D*[ntoys_nh];
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    char histoname[200];
    sprintf(histoname,"h_deltacp_nh_Toy%d",itoy);
    h_chisq_deltacp_nh[itoy] = new TH1D(histoname,histoname,nbins1_nh,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  }  
  TH1D ** h_mix13_deltacp_nh= new TH1D*[ntoys_nh];
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    char histoname[200];
    sprintf(histoname,"h_mix13_deltacp_nh_Toy%d",itoy);
    h_mix13_deltacp_nh[itoy] = new TH1D(histoname,histoname,nbins1_nh,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  }  
  TH1D ** h_chisq_deltacp_nh_bin = new TH1D*[nbins1_nh];
  for(int ibin=0;ibin<nbins1_nh;ibin++){
    char histoname[200];
    sprintf(histoname,"h_chisq_deltacp_nh_bin%d",ibin);
    h_chisq_deltacp_nh_bin[ibin] = new TH1D(histoname,histoname,50,0,500);
  }

  // Loop over toys in the tree

  int iEvtToy = 0;
  for(int i = 0; i < nentries_deltacp_nh; i++){
    iEvtToy = i / nbins1_nh;
    if( !(i%nbins1_nh) ) cout << "# Toy " << iEvtToy << endl;
    
    ch_deltacp_nh.GetEntry(i);
    int migradminstat_nh = fparm_deltacp_nh->GetStatusCode("migrad_min");
    if (migradminstat_nh == 4) cout << "Fit failed!" << endl;
    double mix13_nh = fparm_deltacp_nh->BestFitValue(conv::kPhypMixing13);
    double mix23_nh = fparm_deltacp_nh->BestFitValue(conv::kPhypMixing23);
    double dmsq_nh = fparm_deltacp_nh->BestFitValue(conv::kPhypDm23sq);
    double deltacp_nh = fparm_deltacp_nh->BestFitValue(conv::kPhypCPV);
    double chisq_nh = fparm_deltacp_nh->GetBestFitStatistic();
   
    int bin_idx = h_chisq_deltacp_nh[iEvtToy]->FindBin(deltacp_nh);
    //double bincont = h_chisq_deltacp_nh[iEvtToy]->GetBinContent(bin_idx) / (nentries_deltacp_nh / nbins1_nh);
    double bincont = h_chisq_deltacp_nh[iEvtToy]->GetBinContent(bin_idx);
    if(bincont!=0){
      if(bincont>chisq_nh){
	if(chisq_nh!=0){
	  h_chisq_deltacp_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
	  bestfit_points_nh_x[iEvtToy][bin_idx-1]  = mix13_nh;//sinsq2theta13;
	  bestfit_points_nh_y[iEvtToy][bin_idx-1]  = deltacp_nh;
	}
      }
    } 
    else{ 
      if(chisq_nh!=0){
	h_chisq_deltacp_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
	bestfit_points_nh_x[iEvtToy][bin_idx-1] = mix13_nh;//sinsq2theta13;
	bestfit_points_nh_y[iEvtToy][bin_idx-1] = deltacp_nh;
      }
    }
  }

  TGraph ** gbestfit_nh_Toy = new TGraph*[ntoys_nh]; 

  for(int itoy=0;itoy<ntoys_nh;itoy++){
    char graphname[200];
    sprintf(graphname,"gbestfit_nh_Toy%d",itoy);
    gbestfit_nh_Toy[itoy] = new TGraph( nbins1_nh, bestfit_points_nh_x[itoy] ,bestfit_points_nh_y[itoy] ); // dcp Vs mix13  
    gbestfit_nh_Toy[itoy]->SetName(graphname);
    //cout << graphname << endl;
  }
  
  cout << endl;

  // Get the toy distribution for each dcp value
  for(int ibin=0;ibin<nbins1_nh;ibin++){
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      double bincont = h_chisq_deltacp_nh[itoy]->GetBinContent(ibin+1);
      h_chisq_deltacp_nh_bin[ibin]->Fill(bincont);
    }
  }
    
  // Compute the average sensitivity    
  for(int itoy=0;itoy<ntoys_nh;itoy++){
    h_chisq_deltacp_Ave_nh->Add(h_chisq_deltacp_nh[itoy]);
  }
  h_chisq_deltacp_Ave_nh->Scale(1./ntoys_nh);

  // Find the minimum chi2
  double min_chisq_nh[ntoys_nh];
  for(int itoy=0;itoy<ntoys_nh;itoy++){
    min_chisq_nh[itoy] = 999.;
    cout << "# Toy " << itoy << endl;    
    for(int i = 1; i <= nbins1_nh; i++){
      double bin = h_chisq_deltacp_nh[itoy]->GetBinContent(i);
      if (bin < min_chisq_nh[itoy]) min_chisq_nh[itoy] = bin;
      cout << bin << endl;
    }
    cout << "For deltaCP, min chisq  = " << min_chisq_nh[itoy] << endl;
    cout << endl;
  }

  double min_chisq_Ave_nh=999.;
  for(int i = 1; i <= nbins1_nh; i++){
    double bin = h_chisq_deltacp_Ave_nh->GetBinContent(i);
    if (bin < min_chisq_Ave_nh) min_chisq_Ave_nh = bin;
  }

  // Take RMS of dcp in each true bin
  // Chi2 distribution is gaussian (many bins)

  TCanvas **c_chisq_deltacp_nh_bin = new TCanvas*[nbins1_nh];
  double *RMS_nh = new double[nbins1_nh];
  double *RelErr_nh = new double[nbins1_nh];      

  for(int ibin=0;ibin<nbins1_nh;ibin++){
    char canvname[200];
    sprintf(canvname,"c_chisq_deltacp_nh_bin%d",ibin);
    //c_chisq_deltacp_nh_bin = new TCanvas(canvname,canvname);
    //h_chisq_deltacp_nh_bin[ibin]->Draw();    
    RMS_nh[ibin] = h_chisq_deltacp_nh_bin[ibin]->GetRMS();
    RelErr_nh[ibin] = h_chisq_deltacp_nh_bin[ibin]->GetRMS() / h_chisq_deltacp_nh_bin[ibin]->GetMean();
  }
  
  // Sensi with error bars

  TGraphErrors *g_chisq_deltacp_nh_Err = new TGraphErrors();
  TGraph *g_chisq_deltacp_nh = new TGraph();
  for(int ipt=0;ipt<nbins1_nh;ipt++){
    double contx = h_chisq_deltacp_Ave_nh->GetBinCenter(ipt+1);
    double conty = h_chisq_deltacp_Ave_nh->GetBinContent(ipt+1);
    if(doShiftMin) conty -= min_chisq_Ave_nh; // shift the chi2 wrt the minimum
    double errx = h_chisq_deltacp_Ave_nh->GetBinWidth(ipt+1)/2;
    double erry = RMS_nh[ipt];

    cout << contx << " +- " << errx << " - " << conty << " +- " << erry << endl;

    if(erry<=0.) erry = 0.000001;
    g_chisq_deltacp_nh->SetPoint(ipt,contx,conty);
    g_chisq_deltacp_nh_Err->SetPoint(ipt,contx,conty);
    g_chisq_deltacp_nh_Err->SetPointError(ipt,errx,erry);	
    //cout << contx << " +- " << errx << endl;
  }

  // Shift the chi2 wrt the minimum
  if(doShiftMin){
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      for(int i = 1; i <= nbins1_nh; i++){
	double bin = h_chisq_deltacp_nh[itoy]->GetBinContent(i);
	h_chisq_deltacp_nh[itoy]->SetBinContent(i,bin-min_chisq_nh[itoy]); // scanned global minimum
	//h_chisq_deltacp_nh[itoy]->SetBinContent(i,bin-chisq_bf_nh); // fiited global minimum
      }
    }
    for(int i = 1; i <= nbins1_nh; i++){
      double bin = h_chisq_deltacp_Ave_nh->GetBinContent(i);
      h_chisq_deltacp_Ave_nh->SetBinContent(i,bin-min_chisq_Ave_nh); // scanned global minimum
    }
  }
















  // // 1D vs deltacp (Inverted Hierarchy)

  // const char * filename_deltacp_ih = Form("%s.root",name_ih);
  // TChain ch_deltacp_ih("fitnt");
  // ch_deltacp_ih.Add(filename_deltacp_ih);
  // ch_deltacp_ih.SetBranchAddress("best_fit_params", &fparm_deltacp_ih); 
  // ch_deltacp_ih.SetBranchAddress("nbins_param1", &nbins_param1_ih); 
  // ch_deltacp_ih.SetBranchAddress("nbins_param2", &nbins_param2_ih); 

  // const int nentries_deltacp_ih = ch_deltacp_ih.GetEntries();
  // cout << "nentries deltacp (IH) = " << nentries_deltacp_ih << endl;
  // ch_deltacp_ih.GetEntry(0);
  // const int nbins1_ih = nbins_param1_ih;
 
  // const int ntoys_ih = nentries_deltacp_ih / nbins1_ih;
  // cout << "ntoys (IH) = " << ntoys_ih << endl;

  // mixing23_min = 0.3;
  // mixing23_max = 0.7;
  // double mixing23_step = 0.01;
  
  // dmsq_min = 2.0;//E-3;
  // dmsq_max = 3.0;//E-3;
  // dmsq_step = 0.05;//E-3;

  // double deltacp_min = -conv::kpi;
  // double deltacp_max = conv::kpi;
  // double deltacp_step = (2*conv::kpi)/(nentries_deltacp_nh-1);

  // Double_t** bestfit_points_ih_x = new Double_t*[ntoys_ih]; // bestfit points for each toy
  // for(int itoy=0; itoy<ntoys_ih; itoy++){
  //   bestfit_points_ih_x[itoy] = new Double_t[nbins1_ih];
  // }
  // Double_t** bestfit_points_ih_y = new Double_t*[ntoys_ih]; // bestfit points for each toy
  // for(int itoy=0; itoy<ntoys_ih; itoy++){
  //   bestfit_points_ih_y[itoy] = new Double_t[nbins1_ih];
  // }

  // TH1D * h_chisq_deltacp_Ave_ih = new TH1D("h_chisq_deltacp_Ave_ih","h_chisq_deltacp_Ave_ih",nbins1_ih,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  // TH1D ** h_chisq_deltacp_ih = new TH1D*[ntoys_ih];
  // for(int itoy=0; itoy<ntoys_ih; itoy++){
  //   char histoname[200];
  //   sprintf(histoname,"h_deltacp_ih_Toy%d",itoy);
  //   h_chisq_deltacp_ih[itoy] = new TH1D(histoname,histoname,nbins1_ih,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  // }  
  // TH1D ** h_mix13_deltacp_ih= new TH1D*[ntoys_ih];
  // for(int itoy=0; itoy<ntoys_ih; itoy++){
  //   char histoname[200];
  //   sprintf(histoname,"h_mix13_deltacp_ih_Toy%d",itoy);
  //   h_mix13_deltacp_ih[itoy] = new TH1D(histoname,histoname,nbins1_ih,deltacp_min-0.5*deltacp_step,deltacp_max+0.5*deltacp_step);
  // }  
  // TH1D ** h_chisq_deltacp_ih_bin = new TH1D*[nbins1_ih];
  // for(int ibin=0;ibin<nbins1_ih;ibin++){
  //   char histoname[200];
  //   sprintf(histoname,"h_chisq_deltacp_ih_bin%d",ibin);
  //   h_chisq_deltacp_ih_bin[ibin] = new TH1D(histoname,histoname,50,0,500);
  // }

  // // Loop over toys in the tree

  // int iEvtToy = 0;
  // for(int i = 0; i < nentries_deltacp_ih; i++){
  //   iEvtToy = i / nbins1_ih;
  //   if( !(i%nbins1_ih) ) cout << "# Toy " << iEvtToy << endl;
    
  //   ch_deltacp_ih.GetEntry(i);
  //   int migradminstat_ih = fparm_deltacp_ih->GetStatusCode("migrad_min");
  //   if (migradminstat_ih == 4) cout << "Fit failed!" << endl;
  //   double mix13_ih = fparm_deltacp_ih->BestFitValue(conv::kPhypMixing13);
  //   double mix23_ih = fparm_deltacp_ih->BestFitValue(conv::kPhypMixing23);
  //   double dmsq_ih = fparm_deltacp_ih->BestFitValue(conv::kPhypDm23sq);
  //   double deltacp_ih = fparm_deltacp_ih->BestFitValue(conv::kPhypCPV);
  //   double chisq_ih = fparm_deltacp_ih->GetBestFitStatistic();
   
  //   int bin_idx = h_chisq_deltacp_ih[iEvtToy]->FindBin(deltacp_ih);
  //   //double bincont = h_chisq_deltacp_ih[iEvtToy]->GetBinContent(bin_idx) / (nentries_deltacp_ih / nbins1_ih);
  //   double bincont = h_chisq_deltacp_ih[iEvtToy]->GetBinContent(bin_idx);
  //   if(bincont!=0){
  //     if(bincont>chisq_ih){
  // 	if(chisq_ih!=0){
  // 	  h_chisq_deltacp_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
  // 	  bestfit_points_ih_x[iEvtToy][bin_idx-1]  = mix13_ih;//sinsq2theta13;
  // 	  bestfit_points_ih_y[iEvtToy][bin_idx-1]  = deltacp_ih;
  // 	}
  //     }
  //   } 
  //   else{ 
  //     if(chisq_ih!=0){
  // 	h_chisq_deltacp_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
  // 	bestfit_points_ih_x[iEvtToy][bin_idx-1] = mix13_ih;//sinsq2theta13;
  // 	bestfit_points_ih_y[iEvtToy][bin_idx-1] = deltacp_ih;
  //     }
  //   }
  // }

  // TGraph ** gbestfit_ih_Toy = new TGraph*[ntoys_ih]; 

  // for(int itoy=0;itoy<ntoys_ih;itoy++){
  //   char graphname[200];
  //   sprintf(graphname,"gbestfit_ih_Toy%d",itoy);
  //   gbestfit_ih_Toy[itoy] = new TGraph( nbins1_ih, bestfit_points_ih_x[itoy] ,bestfit_points_ih_y[itoy] ); // dcp Vs mix13  
  //   gbestfit_ih_Toy[itoy]->SetName(graphname);
  //   //cout << graphname << endl;
  // }
  
  // cout << endl;

  // // Get the toy distribution for each dcp value
  // for(int ibin=0;ibin<nbins1_ih;ibin++){
  //   for(int itoy=0;itoy<ntoys_ih;itoy++){
  //     double bincont = h_chisq_deltacp_ih[itoy]->GetBinContent(ibin+1);
  //     h_chisq_deltacp_ih_bin[ibin]->Fill(bincont);
  //   }
  // }
    
  // // Compute the average sensitivity    
  // for(int itoy=0;itoy<ntoys_ih;itoy++){
  //   h_chisq_deltacp_Ave_ih->Add(h_chisq_deltacp_ih[itoy]);
  // }
  // h_chisq_deltacp_Ave_ih->Scale(1./ntoys_ih);

  // // Find the minimum chi2
  // double min_chisq_ih[ntoys_ih];
  // for(int itoy=0;itoy<ntoys_ih;itoy++){
  //   min_chisq_ih[itoy] = 999.;
  //   cout << "# Toy " << itoy << endl;    
  //   for(int i = 1; i <= nbins1_ih; i++){
  //     double bin = h_chisq_deltacp_ih[itoy]->GetBinContent(i);
  //     if (bin < min_chisq_ih[itoy]) min_chisq_ih[itoy] = bin;
  //   }
  //   cout << "For deltaCP, min chisq  = " << min_chisq_ih[itoy] << endl;
  //   cout << endl;
  // }

  // double min_chisq_Ave_ih=999.;
  // for(int i = 1; i <= nbins1_ih; i++){
  //   double bin = h_chisq_deltacp_Ave_ih->GetBinContent(i);
  //   if (bin < min_chisq_Ave_ih) min_chisq_Ave_ih = bin;
  // }

  // // Take RMS of dcp in each true bin
  // // Chi2 distribution is gaussian (many bins)

  // TCanvas **c_chisq_deltacp_ih_bin = new TCanvas*[nbins1_ih];
  // double *RMS_ih = new double[nbins1_ih];
  // double *RelErr_ih = new double[nbins1_ih];      

  // for(int ibin=0;ibin<nbins1_ih;ibin++){
  //   char canvname[200];
  //   sprintf(canvname,"c_chisq_deltacp_ih_bin%d",ibin);
  //   //c_chisq_deltacp_ih_bin = new TCanvas(canvname,canvname);
  //   //h_chisq_deltacp_ih_bin[ibin]->Draw();    
  //   RMS_ih[ibin] = h_chisq_deltacp_ih_bin[ibin]->GetRMS();
  //   RelErr_ih[ibin] = h_chisq_deltacp_ih_bin[ibin]->GetRMS() / h_chisq_deltacp_ih_bin[ibin]->GetMean();
  // }
  
  // // Sensi with error bars

  // TGraphErrors *g_chisq_deltacp_ih_Err = new TGraphErrors();
  // TGraph *g_chisq_deltacp_ih = new TGraph();
  // for(int ipt=0;ipt<nbins1_ih;ipt++){
  //   double contx = h_chisq_deltacp_Ave_ih->GetBinCenter(ipt+1);
  //   double conty = h_chisq_deltacp_Ave_ih->GetBinContent(ipt+1);
  //   if(doShiftMin) conty -= min_chisq_Ave_ih; // shift the chi2 wrt the minimum
  //   double errx = h_chisq_deltacp_Ave_ih->GetBinWidth(ipt+1)/2;
  //   double erry = RMS_ih[ipt];
  //   if(erry<=0.) continue;
  //   g_chisq_deltacp_ih->SetPoint(ipt,contx,conty);
  //   g_chisq_deltacp_ih_Err->SetPoint(ipt,contx,conty);
  //   g_chisq_deltacp_ih_Err->SetPointError(ipt,errx,erry);	
  //   //cout << contx << " +- " << errx << endl;
  // }
    
  // // Shift the chi2 wrt the minimum
  // if(doShiftMin){
  //   for(int itoy=0;itoy<ntoys_ih;itoy++){
  //     for(int i = 1; i <= nbins1_ih; i++){
  // 	double bin = h_chisq_deltacp_ih[itoy]->GetBinContent(i);
  // 	h_chisq_deltacp_ih[itoy]->SetBinContent(i,bin-min_chisq_ih[itoy]); // scanned global minimum
  // 	//h_chisq_deltacp_ih[itoy]->SetBinContent(i,bin-chisq_bf_ih); // fiited global minimum
  //     }
  //   }
  //   for(int i = 1; i <= nbins1_ih; i++){
  //     double bin = h_chisq_deltacp_Ave_ih->GetBinContent(i);
  //     h_chisq_deltacp_Ave_ih->SetBinContent(i,bin-min_chisq_Ave_ih); // scanned global minimum
  //   }
  // }

















  // Draw plots

  // Normal Hierarchy

   TCanvas *c_chisq_deltacp_nh = new TCanvas("c_chisq_deltacp_nh","c_chisq_deltacp_nh");
  double min = TMath::Min(0.,(double)g_chisq_deltacp_nh_Err->GetMinimum());
  min = TMath::Max(0.,min);
  double max = g_chisq_deltacp_nh_Err->GetMaximum() * 1.25;
  g_chisq_deltacp_nh_Err->SetMinimum(max);
  g_chisq_deltacp_nh_Err->SetMinimum(0);
  g_chisq_deltacp_nh_Err->SetFillColor(kRed);
  g_chisq_deltacp_nh_Err->SetFillStyle(3001);
  g_chisq_deltacp_nh_Err->DrawClone("a2");
  g_chisq_deltacp_nh->SetFillColor(kBlue);
  g_chisq_deltacp_nh->DrawClone("same");

  //TCanvas *cdeltacp_nh = new TCanvas("cdeltacp_nh","cdeltacp_nh");
  //h_chisq_deltacp_Ave_nh->SetLineColor(kBlue);
  //h_chisq_deltacp_Ave_nh->GetXaxis()->SetTitle("#delta_{CP} (NH)");
  //h_chisq_deltacp_Ave_nh->GetYaxis()->SetTitle("#Delta#chi^{2}");
  //h_chisq_deltacp_Ave_nh->Draw("histo");
  //h_chisq_deltacp_nh[199]->Draw("histo same");    
  
  return;



  


  for(int itoy=0;itoy<ntoys_nh;itoy++){
    //for(int itoy=1;itoy<5;itoy++){
    h_chisq_deltacp_nh[itoy]->Draw("c same");    
    h_chisq_deltacp_nh[itoy]->SetLineColor(itoy+1);
  }

  // TCanvas *cbestfit_nh = new TCanvas("cbestfit_nh","cbestfit_nh");
  // for(int itoy=0;itoy<ntoys_nh;itoy++){  
  //   gbestfit_nh_Toy[itoy]->GetXaxis()->SetTitle("#delta_{CP} (NH)");
  //   gbestfit_nh_Toy[itoy]->GetYaxis()->SetTitle("#Delta#chi^{2}");
  //   gbestfit_nh_Toy[itoy]->GetYaxis()->SetRangeUser(-3.2,3.2);
  //   gbestfit_nh_Toy[itoy]->SetLineWidth(2);
  //   gbestfit_nh_Toy[itoy]->SetLineColor(kBlue);
  //   if(itoy==0) gbestfit_nh_Toy[itoy]->Draw("ac");
  //   else gbestfit_nh_Toy[itoy]->Draw("c same");
  // }


  return;






  // Inverted Hierarchy

  TCanvas *c_chisq_deltacp_ih = new TCanvas("c_chisq_deltacp_ih","c_chisq_deltacp_ih");
  double min = TMath::Min(0.,(double)g_chisq_deltacp_ih_Err->GetMinimum());
  min = TMath::Max(0.,min);
  double max = g_chisq_deltacp_ih_Err->GetMaximum() * 1.25;
  g_chisq_deltacp_ih_Err->SetMinimum(max);
  g_chisq_deltacp_ih_Err->SetMinimum(0);
  g_chisq_deltacp_ih_Err->SetFillColor(kRed);
  g_chisq_deltacp_ih_Err->SetFillStyle(3001);
  g_chisq_deltacp_ih_Err->DrawClone("a2");
  g_chisq_deltacp_ih->SetFillColor(kBlue);
  g_chisq_deltacp_ih->DrawClone("same");

  // TCanvas *cdeltacp_ih = new TCanvas("cdeltacp_ih","cdeltacp_ih");
  // h_chisq_deltacp_Ave_ih->SetLineColor(kBlue);
  // h_chisq_deltacp_Ave_ih->GetXaxis()->SetTitle("#delta_{CP} (IH)");
  // h_chisq_deltacp_Ave_ih->GetYaxis()->SetTitle("#Delta#chi^{2}");
  // h_chisq_deltacp_Ave_ih->Draw("c");
  // for(int itoy=0;itoy<ntoys_ih;itoy++){
  //   h_chisq_deltacp_ih[itoy]->Draw("c same");    
  // }

  // TCanvas *cbestfit_ih = new TCanvas("cbestfit_ih","cbestfit_ih");
  // for(int itoy=0;itoy<ntoys_ih;itoy++){  
  //   gbestfit_ih_Toy[itoy]->GetXaxis()->SetTitle("#delta_{CP} (IH)");
  //   gbestfit_ih_Toy[itoy]->GetYaxis()->SetTitle("#Delta#chi^{2}");
  //   gbestfit_ih_Toy[itoy]->GetYaxis()->SetRangeUser(-3.2,3.2);
  //   gbestfit_ih_Toy[itoy]->SetLineWidth(2);
  //   gbestfit_ih_Toy[itoy]->SetLineColor(kBlue);
  //   if(itoy==0) gbestfit_ih_Toy[itoy]->Draw("ac");
  //   else gbestfit_ih_Toy[itoy]->Draw("c same");
  // }






  return;
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
