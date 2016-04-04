
// Produce VALOR contours with macros/misc/draw_contours.C 
// and read the output files
//#define SwapXY // Swap x and y axis of contours
using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareContours() 
{

  //////////////////////////// NEW VALOR contours /////////////////////
  
  // Global scan

  const int Nfiles = 2;
  const char *inname[Nfiles] = {
    //"CurrentStat/fitmarg_joint_sinsqth13Vsdcp_react_nh_contour.root",
    //"SensiSummer2016/fitmarg_joint_sinsqth13Vsdcp_react_nh_Summer2016_contour.root"
    
    "cp13/toyfits_dcp_cnh-sinsqtheta[12_0.304,13_0.022,23_0.450],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_contour.root",
    "cp13/toyfits_dcp_cih-sinsqtheta[12_0.304,13_0.022,23_0.450],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_contour.root"
    
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
  };
  const char *chisq_name = "chisq_map_2D";

  const char *bestfitpt_name = "best_fit_point_marker";

  ////////////////////////////////////////////////////////////////////////


  // Read contours 
  
  TFile *filein[Nfiles];
  TGraph *g68_valor[Ncont_68][Nfiles];
  TGraph *g90_valor[Ncont_90][Nfiles];
  TH2D *hchisq_valor[Nfiles];
  TMarker *gbest_fit_global[Nfiles];

  for(int ifile=0;ifile<Nfiles;ifile++){
    filein[ifile] = new TFile(inname[ifile],"READ");
    for(int icont=0;icont<Ncont_68;icont++){
      g68_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_68_name[icont]);
#ifdef SwapXY
      SwapXY(g68_valor[icont][ifile]);
#endif
    }
    for(int icont=0;icont<Ncont_90;icont++){
      g90_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_90_name[icont]);
#ifdef SwapXY
      SwapXY(g90_valor[icont][ifile]);
#endif
    }
    hchisq_valor[ifile] = (TH2D*)filein[ifile]->Get(chisq_name);
    hchisq_valor[ifile]->Reset(); // used to set the full range though tgraphs
#ifdef SwapXY
  hchisq_valor[ifile] = SwapXY(hchisq_valor[ifile]);
#endif
    gbest_fit_global[ifile] = (TMarker*)filein[ifile]->Get(bestfitpt_name);
#ifdef SwapXY
    SwapXY(gbest_fit_global[ifile]);
#endif
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
  lcontours->SetFillColor(kWhite);
  lcontours->SetFillStyle(1);
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
  
//  TLatex *latcontours = new TLatex(0.075,2.,"Global scan");
//  latcontours->SetTextSize(0.07);
//  latcontours->Draw("same");  

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
  
  double ymax = graph->GetXaxis()->GetXmax();
  double xmax = graph->GetYaxis()->GetXmax();

  double ymin = graph->GetXaxis()->GetXmin();
  double xmin = graph->GetYaxis()->GetXmin();

  TH2D* graph_out = new TH2D("h_temp","",100,xmin,xmax,100,ymin,ymax);
  graph_out->GetXaxis()->SetTitle(graph->GetYaxis()->GetTitle());
  graph_out->GetYaxis()->SetTitle(graph->GetXaxis()->GetTitle());

  return graph_out;
}

void SwapXY(TMarker *graph){
  double x = graph->GetX();
  double y = graph->GetY();
  graph->SetX(y);
  graph->SetY(x);

  return;
}