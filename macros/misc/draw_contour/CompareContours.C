
// Produce VALOR contours with macros/misc/draw_contours.C 
// and read the output files

using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareContours() 
{
  //////////////////////////// OLD VALOR contours /////////////////////
  
  //const char *inname_rasterscan_old = "/Users/davidesgalaberna/Desktop/Run14_Study/T2K_run14_Lorena/T2K_data/data_fit_joint_final_nh_abs23_contours_with_2D.root";
  //const char *inname_rasterscan_old = "/Users/davidesgalaberna/Desktop/Run14_Study/T2K_run14_Lorena/T2K_data/data_joint_nh_abs23_contours_common_bf.root";
  const char *inname_rasterscan_old = "/Users/davidesgalaberna/Desktop/Run14_Study/T2K_run14_Lorena/T2K_data/data_joint_ih_abs31_contours_common_bf.root";
  
  //const char *inname_globalscan_old = "/Users/davidesgalaberna/Desktop/Run14_Study/T2K_run14_Lorena/T2K_data/data_joint_nh_abs23_contours_common_bf.root";

  //////////////////////////// NEW VALOR contours /////////////////////
  
  // Global scan

  const int Nfiles = 3;
  const char *inname[Nfiles] = {
    //"/Users/davidesgalaberna/Desktop/Run14_Study/CONTOURS/fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_NH_contour.root",
    //"/Users/davidesgalaberna/Desktop/Run14_Study/CONTOURS/fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_NH_FixMEC_contour.root",
    //"/Users/davidesgalaberna/Desktop/Run14_Study/CONTOURS/fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_NH_StatOnly_contour.root"
    "/Users/davidesgalaberna/Desktop/Run14_Study/CONTOURS/fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_IH_contour.root",
    "/Users/davidesgalaberna/Desktop/Run14_Study/CONTOURS/fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_IH_FixMEC_contour.root",
    "/Users/davidesgalaberna/Desktop/Run14_Study/CONTOURS/fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_IH_StatOnly_contour.root"
  };
  const int Ncont = 2;
  const char *contour_68_name[Ncont] = { // different graph for each contour (+-)
    "g68_valor_0",
    "g68_valor_1"
  };
  const char *contour_90_name[Ncont] = { // different graph for each contour (+-)
    "g90_valor_0",
    "g90_valor_1"
  };
  const char *latexname[Nfiles] = {
    "Stat+Syst",
    "Fix MEC",
    "Stat-only"
  };
  const char *chisq_name = "chisq_map_2D";

  const char *bestfitpt_name = "best_fit_point_marker";

  // Raster scan

  const int Ncl = 2;
  const char *contour_rasterscan_name[Ncl] = { // both 68% and 90% CL
    "h1_rasterscan", // 68% CL
    "h2_rasterscan"  // 90%CL
  };
  const char *latex_rasterscan_name[Ncl][Nfiles] = {
    {"Stat+Syst - 68CL", "Fix MEC - 68CL", "Stat-only - 68CL"},
    {"Stat+Syst - 90CL", "Fix MEC - 90CL", "Stat-only - 90CL"},
  };
  const char *bestfit_rasterscan_name = "best_fit_raster_scan";

  ////////////////////////////////////////////////////////////////////////

  // Read Old VALOR plots

  TFile *filein_rasterscan_old = new TFile(inname_rasterscan_old,"READ");
  TH2D* h_theta13_deltacp_rasterscan_old = (TH2D*)filein_rasterscan_old->Get("h_theta13_deltacp");
  TH2D* h_theta13_deltacp_rasterscan_68_old = (TH2D*)filein_rasterscan_old->Get("h_theta13_deltacp_68");
  TGraph *gbestfit_rasterscan_old = (TGraph*)filein_rasterscan_old->Get("gbestfit");
  //TMarker *bestfit_theta13_deltacp = (TMarker*)filein_rasterscan_old->Get("bestfit_theta13_deltacp");

  // // There's a bug in the contour generation: deltaCP = -0.9 is contained within
  // // 68% and 90% CL contour --> fill those bins by hand

  // // 90% CL
  // int binX_wrong_1 = h_theta13_deltacp_rasterscan_old->GetXaxis()->FindBin(0.02); 
  // int binX_wrong_2 = h_theta13_deltacp_rasterscan_old->GetXaxis()->FindBin(0.045); 
  // int binY_wrong = h_theta13_deltacp_rasterscan_old->GetYaxis()->FindBin(-0.9);
  // for(int iy=0;iy<h_theta13_deltacp_rasterscan_old->GetYaxis()->GetNbins();iy++){
  //   for(int ix=0;ix<h_theta13_deltacp_rasterscan_old->GetXaxis()->GetNbins();ix++){
  //     //double bincont = h_theta13_deltacp_rasterscan_old->GetBinContent(1,1);
  //     double bincont =160.;
  //     if( ((ix+1) < binX_wrong_1) && ((iy+1) == binY_wrong) ){	
  // 	h_theta13_deltacp_rasterscan_old->SetBinContent(ix+1,iy+1,bincont);
  //     }
  //     if( ((ix+1) > binX_wrong_2) && ((iy+1) == binY_wrong) ){	
  // 	h_theta13_deltacp_rasterscan_old->SetBinContent(ix+1,iy+1,bincont);
  //     }
  //   }
  // }
  // // 68% CL
  // int mybinX = h_theta13_deltacp_rasterscan_68_old->GetXaxis()->FindBin(0.045);
  // int mybinY = h_theta13_deltacp_rasterscan_68_old->GetYaxis()->FindBin(-1.);
  // int binX_68_wrong_1 = h_theta13_deltacp_rasterscan_68_old->GetXaxis()->FindBin(0.025); 
  // int binX_68_wrong_2 = h_theta13_deltacp_rasterscan_68_old->GetXaxis()->FindBin(0.04); 
  // int binY_68_wrong = h_theta13_deltacp_rasterscan_68_old->GetYaxis()->FindBin(-0.9);
  // for(int iy=0;iy<h_theta13_deltacp_rasterscan_68_old->GetYaxis()->GetNbins();iy++){
  //   for(int ix=0;ix<h_theta13_deltacp_rasterscan_68_old->GetXaxis()->GetNbins();ix++){
  //     double bincont = h_theta13_deltacp_rasterscan_68_old->GetBinContent(ix+1,iy);
  //     double bincont_left = 100.;
  //     double bincont_right = 200000.;
  //     if( ((ix+1) < binX_68_wrong_1) && ((iy+1) == binY_68_wrong) ){	
  // 	h_theta13_deltacp_rasterscan_68_old->SetBinContent(ix+1,iy+1,bincont_left);
  //     }
  //     if( ((ix+1) > binX_68_wrong_2) && ((iy+1) == binY_68_wrong) ){	
  // 	h_theta13_deltacp_rasterscan_68_old->SetBinContent(ix+1,iy+1,bincont_right);
  //     }
  //   }
  // }
  

  // Read contours 
  
  TFile *filein[Nfiles];
  TGraph *g68_valor[Ncont][Nfiles];
  TGraph *g90_valor[Ncont][Nfiles];
  TH2D *hchisq_valor[Nfiles];
  TH2D *hcont_rasterscan_valor[Ncont][Nfiles];
  TGraph *gbest_fit_raster_scan[Nfiles];
  TMarker *gbest_fit_global[Nfiles];

  for(int ifile=0;ifile<Nfiles;ifile++){
    filein[ifile] = new TFile(inname[ifile],"READ");
    for(int icont=0;icont<Ncont;icont++){
      g68_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_68_name[icont]); 
      g90_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_90_name[icont]); 
      hcont_rasterscan_valor[icont][ifile] = (TH2D*)filein[ifile]->Get(contour_rasterscan_name[icont]);
    }
    hchisq_valor[ifile] = (TH2D*)filein[ifile]->Get(chisq_name);
    hchisq_valor[ifile]->Reset(); // used to set the full range though tgraphs

    gbest_fit_raster_scan[ifile] = (TGraph*)filein[ifile]->Get(bestfit_rasterscan_name);
    gbest_fit_global[ifile] = (TMarker*)filein[ifile]->Get(bestfitpt_name);
  }
  

  // Raster scan contours

  TCanvas *ccontours_rasterscan = new TCanvas("ccontours_rasterscan","ccontours_rasterscan");
  hchisq_valor[0]->Draw(); // to set the range
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int icl=0;icl<Ncl;icl++){
      hcont_rasterscan_valor[icl][ifile]->SetLineWidth(2);
      hcont_rasterscan_valor[icl][ifile]->SetLineColor(ifile+1);
      if     (icl==0) hcont_rasterscan_valor[icl][ifile]->SetLineStyle(2);
      else if(icl==1) hcont_rasterscan_valor[icl][ifile]->SetLineStyle(1);
      hcont_rasterscan_valor[icl][ifile]->DrawCopy("cont3 same");
    }
    gbest_fit_raster_scan[ifile]->SetLineWidth(2);
    gbest_fit_raster_scan[ifile]->SetLineColor(ifile+1);
    gbest_fit_raster_scan[ifile]->SetMarkerColor(ifile+1);
    gbest_fit_raster_scan[ifile]->DrawClone("p same");
  }
  h_theta13_deltacp_rasterscan_68_old->Draw("cont3 same");
  h_theta13_deltacp_rasterscan_old->Draw("cont3 same");
  h_theta13_deltacp_rasterscan_68_old->SetLineColor(kBlue);
  h_theta13_deltacp_rasterscan_68_old->SetLineStyle(2);
  h_theta13_deltacp_rasterscan_old->SetLineColor(kBlue);
  gbestfit_rasterscan_old->Draw("l same");

  TLegend *lcontours_rasterscan = new TLegend(0.57,0.25,0.8,0.65);
  lcontours_rasterscan->SetFillColor(0);
  lcontours_rasterscan->SetLineColor(0);
  lcontours_rasterscan->SetTextFont(132);
  lcontours_rasterscan->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int icl=0;icl<Ncl;icl++){
      lcontours_rasterscan->AddEntry(hcont_rasterscan_valor[icl][ifile],latex_rasterscan_name[icl][ifile],"l");
    }
  }
  lcontours_rasterscan->Draw();
  TLatex *latcontours_rasterscan = new TLatex(0.07,2.,"Raster scan");
  latcontours_rasterscan->SetTextSize(0.07);
  latcontours_rasterscan->Draw("same");
  


  // Calculate difference between bestfit lines
  
  double *reldiff[Nfiles];
  double *absdiff[Nfiles];
  double reldiff_mean[Nfiles] = {0.};
  double absdiff_mean[Nfiles] = {0.};

  for(int ifile=0;ifile<Nfiles;ifile++){

    reldiff[ifile] = new double[gbest_fit_raster_scan[0]->GetN()];
    absdiff[ifile] = new double[gbest_fit_raster_scan[0]->GetN()];

    for(int ipt=0;ipt<gbest_fit_raster_scan[0]->GetN();ipt++){
      double x,y;
      gbest_fit_raster_scan[0]->GetPoint(ipt,x,y);
      double x_old,y_old;
      gbestfit_rasterscan_old->GetPoint(ipt,x_old,y_old);
      if(fabs(y-y_old)>1.E-4){
	cerr << endl;
	cerr << "The old and new points don't match in y" << endl;
	cerr << endl;
	exit(1);
      }
      
      reldiff[ifile][ipt] = (x_old - x) / x_old;
      absdiff[ifile][ipt] = (x_old - x);    

      reldiff_mean[ifile] += reldiff[ifile][ipt] / gbest_fit_raster_scan[0]->GetN();
      absdiff_mean[ifile] += absdiff[ifile][ipt] / gbest_fit_raster_scan[0]->GetN();
    }
  }

  
  
  // Global scan contours

  TCanvas *ccontours = new TCanvas("ccontours","ccontours");
  hchisq_valor[0]->Draw(); // to set the range
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int icont=0;icont<Ncont;icont++){
      g90_valor[icont][ifile]->SetLineWidth(2);
      g90_valor[icont][ifile]->SetLineColor(ifile+1);
      g90_valor[icont][ifile]->SetLineStyle(1);
      g90_valor[icont][ifile]->DrawClone("l same");
      g68_valor[icont][ifile]->SetLineWidth(2);
      g68_valor[icont][ifile]->SetLineColor(ifile+1);
      g68_valor[icont][ifile]->SetLineStyle(2);
      g68_valor[icont][ifile]->DrawClone("l same");
      gbest_fit_global[ifile]->Draw("same");
      gbest_fit_global[ifile]->SetMarkerColor(ifile+1);
    }
  }
  TLegend *lcontours = new TLegend(0.6,0.25,0.8,0.65);
  lcontours->SetFillColor(0);
  lcontours->SetLineColor(0);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    char name[200];
    sprintf(name,"%s - 68CL",latexname[ifile]);
    lcontours->AddEntry(g68_valor[0][ifile],name,"l");
    sprintf(name,"%s - 90CL",latexname[ifile]);
    lcontours->AddEntry(g90_valor[0][ifile],name,"l");
  }
  lcontours->Draw();
  
  TLatex *latcontours = new TLatex(0.075,2.,"Global scan");
  latcontours->SetTextSize(0.07);
  latcontours->Draw("same");
  


  // Print difference between best-fit lines

  cout << endl;
  cout << "Absolute difference (sinsqth13) Old/New best-fit raster scan:" << endl;
  cout << endl;
  for(int ifile=0;ifile<Nfiles;ifile++){
  cout << latexname[ifile] << endl;
    for(int ipt=0;ipt<gbest_fit_raster_scan[0]->GetN();ipt++){
      cout << "sinsqth13(Old)-sinsqth13(New) = " << absdiff[ifile][ipt] 
	   << " --> rel diff = " << reldiff[ifile][ipt] << endl;
    }
    cout << endl;
  }
  cout << endl;
  cout << "Average relative difference:" << endl;
  cout << endl;
  for(int ifile=0;ifile<Nfiles;ifile++){
    cout << latexname[ifile] << " --> " << reldiff_mean[ifile] << endl;
  }
  cout << endl;
  cout << "Average absolute difference:" << endl;
  cout << endl;
  for(int ifile=0;ifile<Nfiles;ifile++){
    cout << latexname[ifile] << " --> " << absdiff_mean[ifile] << endl;
  }
  cout << endl;

  
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
