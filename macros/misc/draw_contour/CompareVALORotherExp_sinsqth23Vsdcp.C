
// Produce VALOR contours with macros/misc/draw_contours.C
// and read the output files

//#define SwapXY // Swap x and y axis of contours // does this to the chisq2d, contour and the TMarker!

using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareVALORotherExp_sinsqth23Vsdcp()
{
  const bool draw_preliminary = true;
  TString str_preliminary = "T2K Run1-7c preliminary";
  
  const double x_low = 0.3;
  const double x_up = 0.7;
  const double y_low = -TMath::Pi();
  const double y_up = +TMath::Pi();
  const char *title_x = "sin^{2} #theta_{23}";
  const char *title_y = "#delta_{CP} (radians)";

  //////////////////////////// VALOR contours /////////////////////
  
  //const char *inname_nh = "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth23Vsdmsq/fitmarg_joint_sinsqth23Vsdmsq_Run17_Data_react_nh_contour.root";
  //const char *inname_ih = "../TECH_NOTE_14_3_16/DATA/RUN17b/contours_sinsqth23Vsdmsq/fitmarg_joint_sinsqth23Vsdmsq_Run17_Data_react_ih_contour.root";

  const char *inname_nh = "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth23Vsdcp_react/fitmarg_joint_sinsqth23Vsdcp_Run17c_Data_react_nh_contour.root";
  const char *inname_ih = "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth23Vsdcp_react/fitmarg_joint_sinsqth23Vsdcp_Run17c_Data_react_ih_contour.root";

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


  //////////////////////////// NOvA contours /////////////////////////

  //const char *inname_nova_nh = "../../../inputs/OfficialContours/OtherExp/NOvAFAResults.root";
  const char *inname_nova_nh = "../../../inputs/OfficialContours/OtherExp/nova_surfaces_neutrino2016.root";
  const int colorNOvA = kGreen+1;
  
  ////////////////////////////////////////////////////////////////////////


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
 

  
  // READ NOvA CONTOURS (only normal hierarchy)
  
  TFile * infile_nova_nh = new TFile(inname_nova_nh,"READ");
  
  TH2F *nova_ssth23_vs_delta_68CL_nh = (TH2F*) infile_nova_nh->Get("nue_ssth23_vs_delta_nh"); 
  nova_ssth23_vs_delta_68CL_nh->SetName("nova_ssth23_vs_delta_68CL_nh");
  nova_ssth23_vs_delta_68CL_nh->SetTitle("nova_ssth23_vs_delta_68CL_nh");  
  TH2F *nova_ssth23_vs_delta_90CL_nh = (TH2F*) nova_ssth23_vs_delta_68CL_nh->Clone("nova_ssth23_vs_delta_nh_90CL");
  nova_ssth23_vs_delta_90CL_nh->SetName("nova_ssth23_vs_delta_90CL_nh");
  nova_ssth23_vs_delta_90CL_nh->SetTitle("nova_ssth23_vs_delta_90CL_nh");

  double chisq_min = nova_ssth23_vs_delta_68CL_nh->GetMinimum();
  nova_ssth23_vs_delta_68CL_nh->SetContour(1);
  nova_ssth23_vs_delta_68CL_nh->SetContourLevel(0,chisq_min+2.3);
  nova_ssth23_vs_delta_90CL_nh->SetContour(1);
  nova_ssth23_vs_delta_90CL_nh->SetContourLevel(0,chisq_min+4.61);

  // Take graphs of the contours

  TMultiGraph *gnova_contour_68CL_nh = (TMultiGraph*) (get_contour("68CL_nova",nova_ssth23_vs_delta_68CL_nh))->Clone();
  TMultiGraph *gnova_contour_90CL_nh = (TMultiGraph*) (get_contour("90CL_nova",nova_ssth23_vs_delta_90CL_nh))->Clone();
  assert(gnova_contour_68CL_nh);
  assert(gnova_contour_90CL_nh);
  
  // Convert from dcp/pi to dcp
  ConvertHistoFromDcpOverPiToDcp(gnova_contour_68CL_nh);
  ConvertHistoFromDcpOverPiToDcp(gnova_contour_90CL_nh);
  
  // Swap dcp with sinsqth23
  SwapXY(gnova_contour_68CL_nh);
  SwapXY(gnova_contour_90CL_nh);
  
  TCanvas *cNova = new TCanvas("cNova","cNova"); 
  gnova_contour_68CL_nh->Draw("ac");
  gnova_contour_90CL_nh->Draw("c");
  

  cout << endl;
  cout << "TO DEBUG IT!!!" << endl;
  cout << " - SwapXY function for TMultiGraph" << endl;
  cout << " - Contour limits" << endl;
  cout << endl;

  return;




  //nova_ssth23_vs_delta_68CL_nh->Draw("colz same");
  //gnova_contour_68CL_nh->Draw("l same");
  //gnova_contour_90CL_nh->Draw("l same");




  
  // Init histogram
  int NbinsX = hchisq_valor_nh->GetNbinsX();
  int NbinsY = hchisq_valor_nh->GetNbinsY();
  TH2D *hInit = new TH2D("hInit","hInit",NbinsX,x_low,x_up,NbinsY,y_low,y_up);
  hInit->GetXaxis()->SetTitle(title_x);
  hInit->GetYaxis()->SetTitle(title_y);
  hInit->Draw();

  double posX = hInit->GetXaxis()->GetBinCenter(2);
  double posY = hInit->GetYaxis()->GetBinCenter(2);
  
  TLatex *latPreliminary = new TLatex(posX,posY,str_preliminary);
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
  // NOvA
  //
  nova_ssth23_vs_delta_90CL_nh->SetLineStyle(1);
  nova_ssth23_vs_delta_90CL_nh->SetLineColor(colorNOvA);
  nova_ssth23_vs_delta_90CL_nh->SetLineWidth(2);
  nova_ssth23_vs_delta_68CL_nh->SetLineStyle(2);
  nova_ssth23_vs_delta_68CL_nh->SetLineColor(colorNOvA);
  nova_ssth23_vs_delta_68CL_nh->SetLineWidth(2);
  //nova_ssth23_vs_delta_68CL_nh->Draw("l same");
  //nova_ssth23_vs_delta_90CL_nh->Draw("l same");
  nova_ssth23_vs_delta_68CL_nh->Draw("cont3 same");
  nova_ssth23_vs_delta_90CL_nh->Draw("cont3 same");
  
  //
  // T2K VALOR
  //
  for(int icont=0;icont<Ncont_90;icont++){
    if(g90_valor_nh[icont]){
      ScaleGraphAxis(g90_valor_nh[icont],1.,1000.);
      g90_valor_nh[icont]->SetLineWidth(2);
      g90_valor_nh[icont]->SetLineColor(colorT2K);
      g90_valor_nh[icont]->SetLineStyle(1);
      g90_valor_nh[icont]->DrawClone("l same");
    }
  }
  for(int icont=0;icont<Ncont_68;icont++){
    if(g68_valor_nh[icont]){
      ScaleGraphAxis(g68_valor_nh[icont],1.,1000.);
      g68_valor_nh[icont]->SetLineWidth(2);
      g68_valor_nh[icont]->SetLineColor(colorT2K);
      g68_valor_nh[icont]->SetLineStyle(2);
      g68_valor_nh[icont]->DrawClone("l same");
    }
  }
  ScaleMarkerAxis(gbest_fit_global_nh,1.,1000.);
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
  lcontours_nh->AddEntry(g68_valor_nh[0],"68%CL","l");
  lcontours_nh->AddEntry(g90_valor_nh[0],"90%CL","l");
  lcontours_nh->AddEntry(gbest_fit_global_nh,"T2K best-fit","p");
  lcontours_nh->Draw();
  
  double posX = hInit->GetXaxis()->GetBinCenter(2);
  double posY = hInit->GetYaxis()->GetBinCenter(2);
  TLatex *latT2K = new TLatex(posX,posY,"T2K");
  latT2K->SetTextSize(0.06);
  latT2K->SetTextColor(colorT2K);
  latT2K->Draw("same");
  TLatex *latNOvA = new TLatex(posX,posY,"NOvA (2016)");
  latNOvA->SetTextSize(0.06);
  latNOvA->SetTextColor(colorNOvA);
  latNOvA->Draw("same");

  TLatex *latNormalHierarchy = new TLatex(posX,posY,"Normal Hierarchy");
  latNormalHierarchy->SetTextSize(0.06);
  latNormalHierarchy->Draw("same");


  return;



  ///////////////////////////////////
  //                               //
  //  CONTOURS INVERTED HIERARCHY  // 
  //                               //
  ///////////////////////////////////

  TCanvas *ccontours_ih = new TCanvas("ccontours_ih","ccontours_ih");
  hInit->Draw();

  if(draw_preliminary) latPreliminary->Draw("same");


  // //
  // // NOvA
  // //
  // nova_ssth23_vs_delta_90CL_ih->SetLineStyle(1);
  // nova_ssth23_vs_delta_90CL_ih->SetLineColor(colorNOvA);
  // nova_ssth23_vs_delta_90CL_ih->SetLineWidth(2);
  // nova_ssth23_vs_delta_68CL_ih->SetLineStyle(2);
  // nova_ssth23_vs_delta_68CL_ih->SetLineColor(colorNOvA);
  // nova_ssth23_vs_delta_68CL_ih->SetLineWidth(2);
  // nova_ssth23_vs_delta_68CL_ih->Draw("l same");
  // nova_ssth23_vs_delta_90CL_ih->Draw("l same");
  
  
  //
  // T2K VALOR
  //
  for(int icont=0;icont<Ncont_90;icont++){
    if(g90_valor_ih[icont]){
      ScaleGraphAxis(g90_valor_ih[icont],1.,1000.);
      g90_valor_ih[icont]->SetLineWidth(2);
      g90_valor_ih[icont]->SetLineColor(colorT2K);
      g90_valor_ih[icont]->SetLineStyle(1);
      g90_valor_ih[icont]->DrawClone("l same");
    }
  }
  for(int icont=0;icont<Ncont_68;icont++){
    if(g68_valor_ih[icont]){
      ScaleGraphAxis(g68_valor_ih[icont],1.,1000.);
      g68_valor_ih[icont]->SetLineWidth(2);
      g68_valor_ih[icont]->SetLineColor(colorT2K);
      g68_valor_ih[icont]->SetLineStyle(2);
      g68_valor_ih[icont]->DrawClone("l same");
    }
  }
  ScaleMarkerAxis(gbest_fit_global_ih,1.,1000.);
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
  lcontours_ih->AddEntry(g68_valor_ih[0],"68%CL","l");
  lcontours_ih->AddEntry(g90_valor_ih[0],"90%CL","l");
  lcontours_ih->AddEntry(gbest_fit_global_ih,"T2K best-fit","p");
  lcontours_ih->Draw();
 
  latT2K->Draw("same");
  //latNOvA->Draw("same");
  
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
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    double x_new = x * scale_x;
    double y_new = y * scale_y;
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,x_new,y_new);
  }
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

void ConvertHistoFromDcpOverPiToDcp(TMultiGraph *mg){

  TList* grlist;
  grlist = mg->GetListOfGraphs();
  printf("This multi-graph contains %d graphs\n",grlist->GetSize());
  TIter next(grlist);
  TObject *obj;
  TGraph *graph;
  
  while ((obj = next())) {
    graph=(TGraph*)obj;
   
    for(int ipt=0;ipt<graph->GetN();ipt++){
      double x = 0.;
      double y = 0.;
      graph->GetPoint(ipt,x,y);
      double dcp = x * TMath::Pi() - TMath::Pi();
      graph->RemovePoint(ipt);
      graph->SetPoint(ipt,dcp,y);
    }
  }

  return;
}

void SwapXY(TMultiGraph *mg){

  TList* grlist;
  grlist = mg->GetListOfGraphs();
  printf("This multi-graph contains %d graphs\n",grlist->GetSize());
  TIter next(grlist);
  TObject *obj;
  TGraph *graph;
  
  while ((obj = next())) {
    graph=(TGraph*)obj;
    
    //
    // THIS PIECE OF CODE DOESN'T WORK
    // IT REMOVES HALF OF THE POINTS OF A TGRAPH
    //
    // for(int ipt=0;ipt<graph->GetN();ipt++){
    //   double x = 0.;
    //   double y = 0.;
    //   graph->GetPoint(ipt,x,y);
    //   graph->RemovePoint(ipt);
    //   graph->SetPoint(ipt,y,x);
    // }
    //
    
    //mg->RecursiveRemove(obj);
    //graph = SwapXY(graph);
    //mg->Add(graph);
  }
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


// void SwapXY(TGraph *graph){
//   for(int ipt=0;ipt<graph->GetN();ipt++){
//     double x = 0.;
//     double y = 0.;
//     graph->GetPoint(ipt,x,y);
//     graph->RemovePoint(ipt);
//     graph->SetPoint(ipt,y,x);
//   }
//   return;
// }

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


TMultiGraph * get_contour(char *tag,TH2F* h3=0){//,TMultiGraph *multig){
  
  TCanvas* c2 = new TCanvas("c2","c2",0,0,900,600);
  h3->Draw("CONT LIST");
  c2->Update();
  TObjArray *contours = 0;
  contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
 
  TList *list90 = (TList*)contours->At(0);
 
  TMultiGraph *multig = new TMultiGraph();
  char multigname[200];
  sprintf(multigname, "multig%s_%d");
  multig->SetTitle(multigname);
  multig->SetName(multigname);
  
  TGraph *g;
  int it=0;
  TIter next(list90);

  while ((g = (TGraph*)next()))
    {
      char gname[200];
      sprintf(gname, "g%s_%d",tag,it);
       
      g->SetLineWidth(2);
      g->SetLineStyle(1);
      g->SetLineColor(kRed);
      g->SetTitle(gname);
      g->SetName(gname);

      multig->Add(g);
      
      //cout << "it = " << it << endl;

      it++;
    }
  
  //  c2 = 0;
  //  contours = 0;
  //  list90 = 0;
  //  g = 0;
  //  delete c2;
  //  delete contours;
  //  delete list90;
  //  delete g;
  
  return multig;
}
