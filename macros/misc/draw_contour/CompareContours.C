
// Produce VALOR contours with macros/misc/draw_contours.C
// and read the output files

//#define SwapXY // Swap x and y axis of contours // does this to the chisq2d, contour and the TMarker!

using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareContours()
{
  
  //////////////////////////// NEW VALOR contours /////////////////////
  
  // Global scan
  
  const int Nfiles = 2;
  const char *inname[Nfiles] = {


    //////////////////////////////////////////////
    //                                          //
    //      Run 1-7 -Asimov - Run 1-7 Flux      // 
    //                                          //
    //////////////////////////////////////////////

    // // sinsqth13 Vs dcp w/ reactors - Asimov A - Run 1-7
    // "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_A_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_A_reactors_ih_contour.root",

    // sinsqth13 Vs dcp w/ reactors - Asimov B - Run 1-7
    "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_B_reactors_nh_contour.root",
    "../TECH_NOTE_14_3_16/SENSI/RUN17/FluxRun17_Conference/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17_Asimov_B_reactors_ih_contour.root",
 

    //////////////////////////////
    //                          //
    //      Run 1-6 -Asimov     // 
    //                          //
    //////////////////////////////

    // sinsqth13 Vs dcp w/ reactors - Asimov A
    //"../TECH_NOTE_14_3_16/SENSI/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Asimov_A_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Asimov_A_reactors_ih_contour.root",

    // sinsqth13 Vs dcp w/ reactors - Asimov B
    //"../TECH_NOTE_14_3_16/SENSI/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Asimov_B_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Asimov_B_reactors_ih_contour.root",
    
    // sinsqth23 Vs dmsq w/ reactors - Asimov A
    //"../TECH_NOTE_14_3_16/SENSI/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Asimov_A_reactors_ih_contour.root",
    
    // // sinsqth23 Vs dmsq w/ reactors - Asimov B
    //"../TECH_NOTE_14_3_16/SENSI/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Asimov_B_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Asimov_B_reactors_ih_contour.root",



    //////////////////////////////
    //                          //
    //      Run 1-6 - Data      // 
    //                          //
    //////////////////////////////

    // // sinsqth23 Vs dmsq w/ reactors - Data
    // "../TECH_NOTE_14_3_16/DATA/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run16_Data_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run16_Data_reactors_ih_contour.root"

    // // sinsqth13 Vs dcp T2K only - Data
    //"../TECH_NOTE_14_3_16/DATA/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_data_t2konly_nh_contour.root",
    //"../TECH_NOTE_14_3_16/DATA/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_data_t2konly_ih_contour.root"

    // // sinsqth13 Vs dcp w/ reactors - Data
    // "../TECH_NOTE_14_3_16/DATA/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run16_Data_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run16_Data_reactors_ih_contour.root"

    // // dmsq Vs dcp w/ reactors - Data
    // "../TECH_NOTE_14_3_16/DATA/contours_dmsqVsdcp_react/fitmarg_joint_dmsqVsdcp_Run16_Data_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/contours_dmsqVsdcp_react/fitmarg_joint_dmsqVsdcp_Run16_Data_reactors_ih_contour.root"

    // // sinsqth23 Vs dcp w/ reactors - Data
    // "../TECH_NOTE_14_3_16/DATA/contours_sinsqth23Vsdcp_react/fitmarg_joint_sinsqth23Vsdcp_Run16_Data_reactors_nh_contour.root",
    // "../TECH_NOTE_14_3_16/DATA/contours_sinsqth23Vsdcp_react/fitmarg_joint_sinsqth23Vsdcp_Run16_Data_reactors_ih_contour.root"

    /////////////// NO FIT BIAS ////////////

    // sinsqth13 Vs dcp T2K reactors - Asimov A 
    //"../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth13Vsdcp_NoFitBias_Run16_Asimov_A_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth13Vsdcp_NoFitBias_Run16_Asimov_A_reactors_ih_contour.root"

    // sinsqth13 Vs dcp T2K reactors - Asimov B 
    //"../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth13Vsdcp_NoFitBias_Run16_Asimov_B_reactors_nh_contour.root"
    //"../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth13Vsdcp_NoFitBias_Run16_Asimov_B_reactors_ih_contour.root"

    // sinsqth23 Vs dmsq T2K reactors - Asimov A 
    //"../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth23Vsdmsq_NoFitBias_Run16_Asimov_A_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth23Vsdmsq_NoFitBias_Run16_Asimov_A_reactors_ih_contour.root"

    // sinsqth23 Vs dmsq T2K reactors - Asimov B
    //"../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth23Vsdmsq_NoFitBias_Run16_Asimov_B_reactors_nh_contour.root",
    //"../TECH_NOTE_14_3_16/SENSI/BANFF_NOFITBIAS_14_3_16/fitmarg_joint_sinsqth23Vsdmsq_NoFitBias_Run16_Asimov_B_reactors_ih_contour.root"
    
  };
  const int Ncont_68 = 2;
  const char *contour_68_name[Ncont_68] = { // different graph for each contour (+-)
    "g68_valor_0",
    "g68_valor_1"
  };
  const int Ncont_90 = 2;
  const char *contour_90_name[Ncont_90] = { // different graph for each contour (+-)
    "g90_valor_0",
    "g90_valor_1"
  };
  const char *latexname[Nfiles] = {
    "Normal Hierarchy",
    "Inverted Hierarchy",
    //"w/ Fit Bias",
    //"w/o Fit Bias",
  };
  const char *chisq_name = "chisq_map_2D";
  
  const char *bestfitpt_name = "best_fit_point_marker";
  
  bool draw_reactor_band = false;
  
  bool shiftglobal = false;


  ////////////////////////////////////////////////////////////////////////
  
  
  // Read contours
  
  TFile *filein[Nfiles];
  TGraph *g68_valor[Ncont_68][Nfiles] = {0};
  TGraph *g90_valor[Ncont_90][Nfiles] = {0};
  TH2D *hchisq_valor[Nfiles] = {0};
  TH2D *hchisq_glob_valor[Nfiles] = {0};
  TH2D *hdeltachisq_glob_valor_68[Nfiles] = {0};
  TH2D *hdeltachisq_glob_valor_90[Nfiles] = {0};
  TMarker *gbest_fit_global[Nfiles] = {0};

  for(int ifile=0;ifile<Nfiles;ifile++){

    filein[ifile] = new TFile(inname[ifile],"READ");
    
    for(int icont=0;icont<Ncont_68;icont++){
    
      g68_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_68_name[icont]);
#ifdef SwapXY
      if(g68_valor[icont][ifile] != NULL)
      SwapXY(g68_valor[icont][ifile]);
#endif
    }
    for(int icont=0;icont<Ncont_90;icont++){
      g90_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_90_name[icont]);
#ifdef SwapXY
      if(g90_valor[icont][ifile] != NULL)
      SwapXY(g90_valor[icont][ifile]);
#endif
    }
    hchisq_valor[ifile] = (TH2D*) filein[ifile]->Get(chisq_name);
    if(hchisq_valor[ifile] == NULL){
      std::cout << "Histogram not found! Exiting" << std::endl;
      exit(0);
    }
    
    hchisq_glob_valor[ifile] = (TH2D*) hchisq_valor[ifile]->Clone(Form("hchisq_glob_%d",ifile));

    hchisq_valor[ifile]->Reset(); // used to set the full range though tgraphs
#ifdef SwapXY
    hchisq_valor[ifile] = SwapXY(hchisq_valor[ifile]);
    hchisq_glob_valor[ifile] = SwapXY(hchisq_glob_valor[ifile]);
#endif
    gbest_fit_global[ifile] = (TMarker*)filein[ifile]->Get(bestfitpt_name);
#ifdef SwapXY
    SwapXY(gbest_fit_global[ifile]);
#endif
  }




  // Init histogram

  int NbinsX = hchisq_glob_valor[0]->GetNbinsX(); 
  double lowedgeX = hchisq_glob_valor[0]->GetXaxis()->GetBinLowEdge(1);
  double upedgeX = hchisq_glob_valor[0]->GetXaxis()->GetBinUpEdge(NbinsX);
  int NbinsY = hchisq_glob_valor[0]->GetNbinsY();
  double lowedgeY = hchisq_glob_valor[0]->GetYaxis()->GetBinLowEdge(1);
  double upedgeY = hchisq_glob_valor[0]->GetYaxis()->GetBinUpEdge(NbinsY);
  char *TitleX = hchisq_glob_valor[0]->GetXaxis()->GetTitle();
  char *TitleY = hchisq_glob_valor[0]->GetYaxis()->GetTitle();
  TH2D *hInit = new TH2D("hInit","hInit",NbinsX,lowedgeX,upedgeX,NbinsY,lowedgeY,upedgeY);
  hInit->GetXaxis()->SetTitle(TitleX);
  hInit->GetYaxis()->SetTitle(TitleY);


  // Get global best-fit contour

  double chisq_min_glob = 99999999.;
  double chisq_min[Nfiles] = 0.;
  for(int ifile=0;ifile<Nfiles;ifile++){
    chisq_min[ifile] = hchisq_glob_valor[ifile]->GetMinimum();
    if(chisq_min[ifile] < chisq_min_glob) chisq_min_glob = chisq_min[ifile];

    cout << "File " << inname[ifile] << " :" << endl;
    cout << "Min chisq = " << chisq_min[ifile] << endl;
    cout << endl;
  }
  cout << "Global min chisq = " << chisq_min_glob << endl;
  cout << endl;

  
  for(int ifile=0;ifile<Nfiles;ifile++){
    
    hdeltachisq_glob_valor_68[ifile] = (TH2D*) hchisq_glob_valor[ifile]->Clone(Form("hchisq_glob_68_%d",ifile));
    hdeltachisq_glob_valor_90[ifile] = (TH2D*) hchisq_glob_valor[ifile]->Clone(Form("hchisq_glob_90_%d",ifile));
    
    for(int i=0;i<hchisq_glob_valor[ifile]->GetNbinsX();i++){
      for(int j=0;j<hchisq_glob_valor[ifile]->GetNbinsY();j++){
	
	double cont = (double) hchisq_glob_valor[ifile]->GetBinContent(i+1,j+1);
	hdeltachisq_glob_valor_68[ifile]->SetBinContent(i+1,j+1,cont-chisq_min_glob);
	hdeltachisq_glob_valor_90[ifile]->SetBinContent(i+1,j+1,cont-chisq_min_glob);
	//hdeltachisq_glob_valor_68[ifile]->SetBinContent(i+1,j+1,cont-chisq_min[ifile]);
	//hdeltachisq_glob_valor_90[ifile]->SetBinContent(i+1,j+1,cont-chisq_min[ifile]);
	//cout << "cont = " << cont << endl;
      }
    }

    
    hdeltachisq_glob_valor_68[ifile]->SetContour(1);
    //hdeltachisq_glob_valor_68[ifile]->SetContourLevel(0,chisq_min+2.3);
    hdeltachisq_glob_valor_68[ifile]->SetContourLevel(0,2.3);
    hdeltachisq_glob_valor_68[ifile]->SetLineWidth(2);
    hdeltachisq_glob_valor_68[ifile]->SetLineColor(ifile+1);
    hdeltachisq_glob_valor_68[ifile]->SetLineStyle(2);
    
    hdeltachisq_glob_valor_90[ifile]->SetContour(1);
    //hdeltachisq_glob_valor_90[ifile]->SetContourLevel(0,chisq_min+4.61);
    hdeltachisq_glob_valor_90[ifile]->SetContourLevel(0,4.61);
    hdeltachisq_glob_valor_90[ifile]->SetLineWidth(2);
    hdeltachisq_glob_valor_90[ifile]->SetLineColor(ifile+1);
    hdeltachisq_glob_valor_90[ifile]->SetLineStyle(1);

  }
  

  



  // Global scan contours

  TCanvas *ccontours = new TCanvas("ccontours","ccontours");
  hchisq_valor[0]->Draw(); // to set the range
  
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int icont=0;icont<Ncont_90;icont++){
      if(g90_valor[icont][ifile]){
        g90_valor[icont][ifile]->SetLineWidth(2);
        g90_valor[icont][ifile]->SetLineColor(ifile+1);
        g90_valor[icont][ifile]->SetLineStyle(1);
        g90_valor[icont][ifile]->DrawClone("l same");
      }
    }
    for(int icont=0;icont<Ncont_68;icont++){
      if(g68_valor[icont][ifile]){
        g68_valor[icont][ifile]->SetLineWidth(2);
        g68_valor[icont][ifile]->SetLineColor(ifile+1);
        g68_valor[icont][ifile]->SetLineStyle(2);
        g68_valor[icont][ifile]->DrawClone("l same");
      }
    }
    gbest_fit_global[ifile]->Draw("same");
    gbest_fit_global[ifile]->SetMarkerColor(ifile+1);
  }
  TLegend *lcontours = new TLegend(0.5,0.25,0.8,0.65);
  lcontours->SetBorderSize(0);
  lcontours->SetFillColor(0);
  lcontours->SetLineColor(0);
  //lcontours->SetLineColorAlpha(kWhite,0.);
  lcontours->SetFillColor(kWhite);
  lcontours->SetFillStyle(1);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    char name[200];
    sprintf(name,"%s - 68%%CL",latexname[ifile]);
    lcontours->AddEntry(g68_valor[0][ifile],name,"l");
    sprintf(name,"%s - 90%%CL",latexname[ifile]);
    lcontours->AddEntry(g90_valor[0][ifile],name,"l");
  }
  
  // Draw reactors box
  
  //double X1 = kGlobFitSinsqTheta13 - kGlobFitSinsqTheta13*kGlobFitErrSinsqTheta13;
  //double X2 = kGlobFitSinsqTheta13 + kGlobFitSinsqTheta13*kGlobFitErrSinsqTheta13;
  //cout << kGlobFitSinsqTheta13 << " +- " << kGlobFitSinsqTheta13*kGlobFitErrSinsqTheta13 << endl;

  const double sinsqth13_nom = 0.0219;
  const double sinsqth13_err = 0.0012;
  double X1 = sinsqth13_nom + sinsqth13_err;
  double X2 = sinsqth13_nom - sinsqth13_err;
  
  double Y1 = hchisq_valor[0]->GetYaxis()->GetBinLowEdge(hchisq_valor[0]->GetYaxis()->GetFirst());
  //double Y2 = hInit->GetMaximum();
  double Y2 = hchisq_valor[0]->GetYaxis()->GetBinUpEdge(hchisq_valor[0]->GetYaxis()->GetLast());
 
  TBox *breactors = new TBox(X1,Y1,X2,Y2);
  breactors->SetFillColor(kYellow+1);
  breactors->SetLineColor(kYellow+1);
  //breactors->SetFillColorAlpha(kYellow+1,0.5);
  if(draw_reactor_band) breactors->Draw("");
  
  lcontours->Draw();
  
  //  TLatex *latcontours = new TLatex(0.075,2.,"Global scan");
  //  latcontours->SetTextSize(0.07);
  //  latcontours->Draw("same");
  





  // TCanvas *cDeltaChisq[Nfiles];
  // for(int ifile=0;ifile<Nfiles;ifile++){
  //   cDeltaChisq[ifile] = new TCanvas(Form("cDeltaChisq_%d",ifile),Form("cDeltaChisq_%d",ifile));
  //   if(ifile==0) hdeltachisq_glob_valor_68[ifile]->Draw("colz");
  //   else hdeltachisq_glob_valor_68[ifile]->Draw("colz same");
  //   hdeltachisq_glob_valor_90[ifile]->Draw("colz same");
  // }



  TCanvas *cContGlobMin = new TCanvas("cContGlobMin","cContGlobMin");
  hInit->Draw();
  for(int ifile=0;ifile<Nfiles;ifile++){
    hdeltachisq_glob_valor_68[ifile]->Draw("cont3 same");
    hdeltachisq_glob_valor_90[ifile]->Draw("cont3 same");

    gbest_fit_global[ifile]->Draw("same");
    gbest_fit_global[ifile]->SetMarkerColor(ifile+1);
  }
  if(draw_reactor_band) breactors->Draw("");
  lcontours->Draw();
  



  return;
}



void ConvertHistoFromSinsqthToSinsq2th(TGraph *graph){
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    double sinsq2th = util::FromSinsqthToSinsq2th(x);
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,sinsq2th,y);
  }
  return;
}

void SwapXY(TGraph *graph){
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,y,x);
  }
  return;
}

TH2D* SwapXY(TH2D *graph){
  
  int nbinsX = graph->GetNbinsX();
  int nbinsY = graph->GetNbinsY();
  double xmin = graph->GetXaxis()->GetXmin();
  double xmax = graph->GetXaxis()->GetXmax();
  double ymin = graph->GetYaxis()->GetXmin();
  double ymax = graph->GetYaxis()->GetXmax();
  
  TH2D* graph_out = new TH2D("h_temp","",nbinsY,ymin,ymax,nbinsX,xmin,xmax);
  graph_out->GetXaxis()->SetTitle(graph->GetYaxis()->GetTitle());
  graph_out->GetYaxis()->SetTitle(graph->GetXaxis()->GetTitle());
  
  for(int i=0;i<graph_out->GetNbinsX();i++){
    for(int j=0;j<graph_out->GetNbinsY();j++){  
      double bincont = graph->GetBinContent(j+1,i+1);
      graph_out->SetBinContent(i+1,j+1,bincont);
    }
  }
  
  return graph_out;
}

void SwapXY(TMarker *graph){
  double x = graph->GetX();
  double y = graph->GetY();
  graph->SetX(y);
  graph->SetY(x);
  
  return;
}
