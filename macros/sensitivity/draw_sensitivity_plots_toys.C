//
// Sensitivity study w/ 1D or 2D plost
//

#define doShiftMin true
#define doNH true
#define doIH false

void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style);

void draw_sensitivity_plots_toys(
			    bool reactor = true, 
			    bool theta2theta = false,
			    const char * name_nh = "../inputs/toyexptsensi_PROVA_3toys_New",
			    const char * name_ih = "../inputs/toyexptsensi_PROVA"
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
    h_chisq_deltacp_nh_bin[ibin] = new TH1D(histoname,histoname,50,50,150);
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
    // char canvname[200];
    // sprintf(canvname,"c_chisq_deltacp_nh_bin%d",ibin);
    // c_chisq_deltacp_nh_bin = new TCanvas(canvname,canvname);
    // h_chisq_deltacp_nh_bin[ibin]->Draw();    
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
    if(erry<=0.) continue;
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

  // Draw plots

  TCanvas *c = new TCanvas("c","c");
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

  TCanvas *cdeltacp_nh = new TCanvas("cdeltacp_nh","cdeltacp_nh");
  h_chisq_deltacp_Ave_nh->SetLineColor(kBlue);
  h_chisq_deltacp_Ave_nh->GetXaxis()->SetTitle("#delta_{CP} (NH)");
  h_chisq_deltacp_Ave_nh->GetYaxis()->SetTitle("#Delta#chi^{2}");
  h_chisq_deltacp_Ave_nh->Draw("c");
  for(int itoy=0;itoy<ntoys_nh;itoy++){
    h_chisq_deltacp_nh[itoy]->Draw("c same");    
  }

  TCanvas *cbestfit_nh = new TCanvas("cbestfit_nh","cbestfit_nh");
  for(int itoy=0;itoy<ntoys_nh;itoy++){  
    gbestfit_nh_Toy[itoy]->GetXaxis()->SetTitle("#delta_{CP} (NH)");
    gbestfit_nh_Toy[itoy]->GetYaxis()->SetTitle("#Delta#chi^{2}");
    gbestfit_nh_Toy[itoy]->GetYaxis()->SetRangeUser(-3.2,3.2);
    gbestfit_nh_Toy[itoy]->SetLineWidth(2);
    gbestfit_nh_Toy[itoy]->SetLineColor(kBlue);
    if(itoy==0) gbestfit_nh_Toy[itoy]->Draw("ac");
    else gbestfit_nh_Toy[itoy]->Draw("c same");
  }









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
