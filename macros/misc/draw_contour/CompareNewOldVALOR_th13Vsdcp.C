

// Produce 1D chisq Vs dcp --> plot_contours_dcp.C (both NH and IH)
// Produce 2D chisq Vs dcp Vs sinsqth13 --> draw_contour.C 

//#define SwapXY

using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareNewOldVALOR_th13Vsdcp
( 
 ) 
{
  const bool draw_preliminary = true;
  const string str_preliminary = "T2K Run1-7c preliminary";

  const double x_low = 0.06;
  const double x_up = 0.2;
  const double y_low = -TMath::Pi();
  const double y_up = +TMath::Pi();

  const char *title_x = "sin^{2} 2#theta_{13}";
  const char *title_y = "#delta_{CP} (radians)";


  //////////////////////////// VALOR contours /////////////////////
  
  //const char *inname_nh = "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run17_Data_react_nh_contour.root";
  //const char *inname_ih = "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth13Vsdcp/fitmarg_joint_sinsqth13Vsdcp_Run17_Data_react_ih_contour.root";

  const char *inname_nh = "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17c_Data_react_nh_contour.root";
  const char *inname_ih = "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth13Vsdcp_react/fitmarg_joint_sinsqth13Vsdcp_Run17c_Data_react_ih_contour.root";

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
  
  const char *chisq_name = "chisq_map_2D";
  const char *bestfitpt_name = "best_fit_point_marker";
  
  const int colorT2K = kBlack;

  //////////////////////////// Run 1-4 contours /////////////////////
 
  const char *inname_run14_68CL_nh_A = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_NH_68CL_digitized_TN175v4_fig52_PartA.txt"; // IH
  const char *inname_run14_68CL_nh_B = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_NH_68CL_digitized_TN175v4_fig52_PartB.txt"; // IH
  const char *inname_run14_90CL_nh_A = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_NH_90CL_digitized_TN175v4_fig52_PartA.txt"; // IH
  const char *inname_run14_90CL_nh_B = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_NH_90CL_digitized_TN175v4_fig52_PartB.txt"; // IH
  const char *inname_run14_bf_nh   = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_NH_bf_digitized_TN175v4_fig52.txt"; // IH
  const char *inname_run14_68CL_ih_A = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_IH_68CL_digitized_TN175v4_fig54_PartA.txt"; // NH
  const char *inname_run14_68CL_ih_B = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_IH_68CL_digitized_TN175v4_fig54_PartB.txt"; // NH
  const char *inname_run14_90CL_ih_A = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_IH_90CL_digitized_TN175v4_fig54_PartA.txt"; // NH
  const char *inname_run14_90CL_ih_B = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_IH_90CL_digitized_TN175v4_fig54_PartB.txt"; // NH
  const char *inname_run14_bf_ih   = "../../../inputs/OfficialContours/T2K/Run14/digitized/sinsq2th13Vsdcp_react_data_fixed_IH_bf_digitized_TN175v4_fig52.txt"; // NH
  const int colorRun14 = kRed;

  ///////////////////////////////////////////////////////////////////

  // READ VALOR CONTOURS
  
  TFile *filein_nh;
  TFile *filein_ih;
  TGraph *g68_valor_nh[Ncont_68] = {0};
  TGraph *g68_valor_ih[Ncont_68] = {0};
  TGraph *g90_valor_nh[Ncont_90] = {0};
  TGraph *g90_valor_ih[Ncont_90] = {0};
  TH2D *hchisq_valor_nh = 0;
  TH2D *hchisq_valor_ih = 0;
  TMarker *gbest_fit_global_nh = 0;
  TMarker *gbest_fit_global_ih = 0;
   
  filein_nh = new TFile(inname_nh,"READ");
  for(int icont=0;icont<Ncont_68;icont++){
    g68_valor_nh[icont] = (TGraph*)filein_nh->Get(contour_68_name[icont]);
#ifdef SwapXY
    if(g68_valor_nh[icont] != NULL)
      g68_valor_nh[icont] = SwapXY(g68_valor_nh[icont]);
#endif
  }
  for(int icont=0;icont<Ncont_90;icont++){
    g90_valor_nh[icont] = (TGraph*)filein_nh->Get(contour_90_name[icont]);
#ifdef SwapXY
    if(g90_valor_nh[icont] != NULL)
      g90_valor_nh[icont] = SwapXY(g90_valor_nh[icont]);
#endif
  }
  hchisq_valor_nh = (TH2D*) filein_nh->Get(chisq_name);
  if(hchisq_valor_nh == NULL){
    std::cout << "Histogram not found! Exiting" << std::endl;
    exit(0);
  }  
  gbest_fit_global_nh = (TMarker*)filein_nh->Get(bestfitpt_name);
#ifdef SwapXY
  SwapXY(gbest_fit_global_nh);
#endif
 
  filein_ih = new TFile(inname_ih,"READ");
  for(int icont=0;icont<Ncont_68;icont++){
    g68_valor_ih[icont] = (TGraph*)filein_ih->Get(contour_68_name[icont]);
#ifdef SwapXY
    if(g68_valor_ih[icont] != NULL)
      g68_valor_ih[icont] = SwapXY(g68_valor_ih[icont]);
#endif
  }
  for(int icont=0;icont<Ncont_90;icont++){
    g90_valor_ih[icont] = (TGraph*)filein_ih->Get(contour_90_name[icont]);
#ifdef SwapXY
    if(g90_valor_ih[icont] != NULL)
      g90_valor_ih[icont] = SwapXY(g90_valor_ih[icont]);
#endif
  }
  gbest_fit_global_ih = (TMarker*)filein_ih->Get(bestfitpt_name);
#ifdef SwapXY
  SwapXY(gbest_fit_global_ih);
#endif


  // READ RUN 1-4 CONTOURS
  
  // TFile * infile_run14_nh = new TFile(inname_run14_nh,"READ");
  // TH2D *run14_NH_S14_68CL = (TH2D*) infile_run14_nh->Get("h_theta23_dmsq_68");
  // TH2D *run14_NH_S23M23_90CL = (TH2D*) infile_run14_nh->Get("h_theta23_dmsq");
  // TMarker *run14_gbest_fit_global_nh = (TMarker*) infile_run14_nh->Get("bestfit_theta23_dmsq");
  // TFile * infile_run14_ih = new TFile(inname_run14_ih,"READ");
  // TH2D *run14_IH_S23M23_68CL = (TH2D*) infile_run14_ih->Get("h_theta23_dmsq_68");
  // TH2D *run14_IH_S23M23_90CL = (TH2D*) infile_run14_ih->Get("h_theta23_dmsq");
  // TMarker *run14_gbest_fit_global_ih = (TMarker*) infile_run14_ih->Get("bestfit_theta23_dmsq");
  
  TGraph *run14_NH_S23M23_68CL_A = new TGraph(inname_run14_68CL_nh_A);
  TGraph *run14_NH_S23M23_68CL_B = new TGraph(inname_run14_68CL_nh_B);
  TGraph *run14_NH_S23M23_90CL_A = new TGraph(inname_run14_90CL_nh_A);
  TGraph *run14_NH_S23M23_90CL_B = new TGraph(inname_run14_90CL_nh_B);
  TGraph *run14_IH_S23M23_68CL_A = new TGraph(inname_run14_68CL_ih_A);
  TGraph *run14_IH_S23M23_68CL_B = new TGraph(inname_run14_68CL_ih_B);
  TGraph *run14_IH_S23M23_90CL_A = new TGraph(inname_run14_90CL_ih_A);
  TGraph *run14_IH_S23M23_90CL_B = new TGraph(inname_run14_90CL_ih_B);
  TGraph *run14_gbest_fit_global_nh = new TGraph(inname_run14_bf_nh); 
  TGraph *run14_gbest_fit_global_ih = new TGraph(inname_run14_bf_ih); 
  

  // Init histogram
  int NbinsX = hchisq_valor_nh->GetNbinsX();
  int NbinsY = hchisq_valor_nh->GetNbinsY();
  TH2D *hInit = new TH2D("hInit","hInit",NbinsX,x_low,x_up,NbinsY,y_low,y_up);
  hInit->GetXaxis()->SetTitle(title_x);
  hInit->GetYaxis()->SetTitle(title_y);
  hInit->Draw();

  double posX = hInit->GetXaxis()->GetBinCenter(2);
  double posY = hInit->GetYaxis()->GetBinCenter(2);
  TLatex *latPreliminary = new TLatex(posX,posY,str_preliminary.c_str());
  latPreliminary->SetTextSize(0.05);

  ///////////////////////////////////
  //                               //
  //   CONTOURS NORMAL HIERARCHY   // 
  //                               //
  ///////////////////////////////////

  TCanvas *ccontours_nh = new TCanvas("ccontours_nh","ccontours_nh");
  hInit->Draw();

  if(draw_preliminary) latPreliminary->Draw("same");

  //
  // Run 1-4
  //  
  run14_NH_S23M23_90CL_A->SetLineStyle(1);
  run14_NH_S23M23_90CL_A->SetLineColor(colorRun14);
  run14_NH_S23M23_90CL_A->SetLineWidth(2);
  run14_NH_S23M23_90CL_B->SetLineStyle(1);
  run14_NH_S23M23_90CL_B->SetLineColor(colorRun14);
  run14_NH_S23M23_90CL_B->SetLineWidth(2);
  run14_NH_S23M23_68CL_A->SetLineStyle(2);
  run14_NH_S23M23_68CL_A->SetLineColor(colorRun14);
  run14_NH_S23M23_68CL_A->SetLineWidth(2);
  run14_NH_S23M23_68CL_B->SetLineStyle(2);
  run14_NH_S23M23_68CL_B->SetLineColor(colorRun14);
  run14_NH_S23M23_68CL_B->SetLineWidth(2);
  run14_NH_S23M23_90CL_A ->Draw("l same"); 
  run14_NH_S23M23_90CL_B ->Draw("l same"); 
  run14_NH_S23M23_68CL_A ->Draw("l same"); 
  run14_NH_S23M23_68CL_B ->Draw("l same"); 
  run14_gbest_fit_global_nh->Draw("p same");
  run14_gbest_fit_global_nh->SetMarkerColor(colorRun14);
  run14_gbest_fit_global_nh->SetMarkerStyle(29);
  run14_gbest_fit_global_nh->SetMarkerSize(1.5);


  //
  // T2K VALOR
  //
  for(int icont=0;icont<Ncont_90;icont++){
    if(g90_valor_nh[icont]){
      //ScaleGraphAxis(g90_valor_nh[icont],1.,1000.);
      ConvertGraphFromSinsqthToSinsq2th(g90_valor_nh[icont]);
      g90_valor_nh[icont]->SetLineWidth(2);
      g90_valor_nh[icont]->SetLineColor(colorT2K);
      g90_valor_nh[icont]->SetLineStyle(1);
      g90_valor_nh[icont]->DrawClone("l same");
    }
  }
  for(int icont=0;icont<Ncont_68;icont++){
    if(g68_valor_nh[icont]){
      //ScaleGraphAxis(g68_valor_nh[icont],1.,1000.);
      ConvertGraphFromSinsqthToSinsq2th(g68_valor_nh[icont]);
      g68_valor_nh[icont]->SetLineWidth(2);
      g68_valor_nh[icont]->SetLineColor(colorT2K);
      g68_valor_nh[icont]->SetLineStyle(2);
      g68_valor_nh[icont]->DrawClone("l same");
    }
  }

  //ScaleMarkerAxis(gbest_fit_global_nh,1.,1000.);
  ConvertMarkerFromSinsqthToSinsq2th(gbest_fit_global_nh);
  gbest_fit_global_nh->Draw("same");
  gbest_fit_global_nh->SetMarkerColor(colorT2K);

  TLegend *lcontours_nh = new TLegend(0.5,0.5,0.8,0.65);
  lcontours_nh->SetBorderSize(0);
  lcontours_nh->SetFillColor(0);
  lcontours_nh->SetLineColor(0);
  lcontours_nh->SetFillColor(kWhite);
  lcontours_nh->SetFillStyle(1);
  lcontours_nh->SetTextFont(132);
  lcontours_nh->SetTextSize(0.05);
  char name[200];
  lcontours_nh->AddEntry(g68_valor_nh[0],"68.27%CL (-2ln L_{crit} = 2.3)","l");
  lcontours_nh->AddEntry(g90_valor_nh[0],"90%CL (-2ln L_{crit} = 4.61)","l");
  lcontours_nh->AddEntry(gbest_fit_global_nh,"Best-fit","p");
  lcontours_nh->Draw();

  TLegend *lcontours2_nh = new TLegend(0.5,0.5,0.8,0.65);
  lcontours2_nh->SetBorderSize(0);
  lcontours2_nh->SetFillColor(0);
  lcontours2_nh->SetLineColor(0);
  lcontours2_nh->SetFillColor(kWhite);
  lcontours2_nh->SetFillStyle(1);
  lcontours2_nh->SetTextFont(132);
  lcontours2_nh->SetTextSize(0.05);
  lcontours2_nh->AddEntry(g90_valor_nh[0],"T2K Run 1-7c","l");
  lcontours2_nh->AddEntry(run14_NH_S23M23_90CL_A,"T2K Run 1-4","l");
  lcontours2_nh->Draw();
  
  TLatex *latNormalHierarchy = new TLatex(posX,posY,"Normal Hierarchy");
  latNormalHierarchy->SetTextSize(0.06);
  latNormalHierarchy->Draw("same");


  ///////////////////////////////////
  //                               //
  //  CONTOURS INVERTED HIERARCHY  // 
  //                               //
  ///////////////////////////////////

  TCanvas *ccontours_ih = new TCanvas("ccontours_ih","ccontours_ih");
  hInit->Draw();
  if(draw_preliminary) latPreliminary->Draw("same");

  //
  // Run 1-4
  //  
  run14_IH_S23M23_90CL_A->SetLineStyle(1);
  run14_IH_S23M23_90CL_A->SetLineColor(colorRun14);
  run14_IH_S23M23_90CL_A->SetLineWidth(2);
  run14_IH_S23M23_90CL_B->SetLineStyle(1);
  run14_IH_S23M23_90CL_B->SetLineColor(colorRun14);
  run14_IH_S23M23_90CL_B->SetLineWidth(2);
  run14_IH_S23M23_68CL_A->SetLineStyle(2);
  run14_IH_S23M23_68CL_A->SetLineColor(colorRun14);
  run14_IH_S23M23_68CL_A->SetLineWidth(2);
  run14_IH_S23M23_68CL_B->SetLineStyle(2);
  run14_IH_S23M23_68CL_B->SetLineColor(colorRun14);
  run14_IH_S23M23_68CL_B->SetLineWidth(2);
  run14_IH_S23M23_90CL_A ->Draw("l same"); 
  run14_IH_S23M23_90CL_B ->Draw("l same"); 
  run14_IH_S23M23_68CL_A ->Draw("l same"); 
  run14_IH_S23M23_68CL_B ->Draw("l same"); 
  run14_gbest_fit_global_ih->Draw("p same");
  run14_gbest_fit_global_ih->SetMarkerColor(colorRun14);
  run14_gbest_fit_global_ih->SetMarkerStyle(29);
  run14_gbest_fit_global_ih->SetMarkerSize(1.5);


  //
  // T2K VALOR
  //
  for(int icont=0;icont<Ncont_90;icont++){
    if(g90_valor_ih[icont]){
      //ScaleGraphAxis(g90_valor_ih[icont],1.,1000.);
      ConvertGraphFromSinsqthToSinsq2th(g90_valor_ih[icont]);
      g90_valor_ih[icont]->SetLineWidth(2);
      g90_valor_ih[icont]->SetLineColor(colorT2K);
      g90_valor_ih[icont]->SetLineStyle(1);
      g90_valor_ih[icont]->DrawClone("l same");
    }
  }
  for(int icont=0;icont<Ncont_68;icont++){
    if(g68_valor_ih[icont]){
      //ScaleGraphAxis(g68_valor_ih[icont],1.,1000.);
      ConvertGraphFromSinsqthToSinsq2th(g68_valor_ih[icont]);
      g68_valor_ih[icont]->SetLineWidth(2);
      g68_valor_ih[icont]->SetLineColor(colorT2K);
      g68_valor_ih[icont]->SetLineStyle(2);
      g68_valor_ih[icont]->DrawClone("l same");
    }
  }

  //ScaleMarkerAxis(gbest_fit_global_ih,1.,1000.);
  ConvertMarkerFromSinsqthToSinsq2th(gbest_fit_global_ih);
  gbest_fit_global_ih->Draw("same");
  gbest_fit_global_ih->SetMarkerColor(colorT2K);

  TLegend *lcontours_ih = new TLegend(0.5,0.5,0.8,0.65);
  lcontours_ih->SetBorderSize(0);
  lcontours_ih->SetFillColor(0);
  lcontours_ih->SetLineColor(0);
  lcontours_ih->SetFillColor(kWhite);
  lcontours_ih->SetFillStyle(1);
  lcontours_ih->SetTextFont(132);
  lcontours_ih->SetTextSize(0.05);
  char name[200];
  lcontours_ih->AddEntry(g68_valor_ih[0],"68.27%CL (-2ln L_{crit} = 2.3)","l");
  lcontours_ih->AddEntry(g90_valor_ih[0],"90%CL (-2ln L_{crit} = 4.61)","l");
  lcontours_ih->AddEntry(gbest_fit_global_ih,"Best-fit","p");
  lcontours_ih->Draw();

  TLegend *lcontours2_ih = new TLegend(0.5,0.5,0.8,0.65);
  lcontours2_ih->SetBorderSize(0);
  lcontours2_ih->SetFillColor(0);
  lcontours2_ih->SetLineColor(0);
  lcontours2_ih->SetFillColor(kWhite);
  lcontours2_ih->SetFillStyle(1);
  lcontours2_ih->SetTextFont(132);
  lcontours2_ih->SetTextSize(0.05);
  lcontours2_ih->AddEntry(g90_valor_ih[0],"T2K Run 1-7c","l");
  lcontours2_ih->AddEntry(run14_IH_S23M23_90CL_A,"T2K Run 1-4","l");
  lcontours2_ih->Draw();
  
  posX = hInit->GetXaxis()->GetBinCenter(2);
  posY = hInit->GetYaxis()->GetBinCenter(2);
  TLatex *latInvertedHierarchy = new TLatex(posX,posY,"Inverted Hierarchy");
  latInvertedHierarchy->SetTextSize(0.06);
  latInvertedHierarchy->Draw("same");

  return;
}


void ScaleMarkerAxis(TMarker *marker,double scale_x,double scale_y){
  double x = marker->GetX();
  double y = marker->GetY();
  double x_new = x * scale_x;
  double y_new = y * scale_y;
  marker->SetX(x_new);
  marker->SetY(y_new);
  return;
}
void ScaleGraphAxis(TGraph *graph,double scale_x,double scale_y){
  double *x = graph->GetX();
  double *y = graph->GetY();
  for(int ipt=0;ipt<graph->GetN();ipt++){
    x[ipt] = x[ipt] * scale_x;
    y[ipt] = y[ipt] * scale_y;
    //cout << graph->GetN() << endl;
  }
  return;
}

void ConvertHistoFromSinsqthToSinsq2th(TGraph *graph){
  double *x = graph->GetX();
  double *y = graph->GetY();
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double sinsq2th = util::FromSinsqthToSinsq2th(x);
    x[ipt] = sinsq2th;
  }
  return;
}

void ConvertMarkerFromSinsqthToSinsq2th(TMarker *marker){
  double x = marker->GetX();
  double y = marker->GetY();
  double sinsq2th = util::FromSinsqthToSinsq2th(x);
  marker->SetX(sinsq2th);
    
  return;
}

TGraph *SwapXY(TGraph *graph){
  TGraph *graph_swap = new TGraph();
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);   
    graph_swap->SetPoint(ipt,y,x);
   }
  return graph_swap;
}

TH2D* SwapXY(TH2D *graph){
  
  int nbinsX = graph->GetNbinsX();
  int nbinsY = graph->GetNbinsY();
  double xmin = graph->GetXaxis()->GetXmin();
  double xmax = graph->GetXaxis()->GetXmax();
  double ymin = graph->GetYaxis()->GetXmin();
  double ymax = graph->GetYaxis()->GetXmax();
 
  //TH2D* graph_out = new TH2D("h_temp","",nbinsY,ymin,ymax,nbinsX,xmin,xmax);
 
  TArrayD *arrayx = graph->GetXaxis()->GetXbins();
  TArrayD *arrayy = graph->GetYaxis()->GetXbins();
  TH2D* graph_out = new TH2D("h_temp","",nbinsY,arrayy->GetArray(),nbinsX,arrayx->GetArray());

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
