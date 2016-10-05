
// Produce VALOR contours with macros/misc/draw_contours.C
// and read the output files

//#define SwapXY // Swap x and y axis of contours // does this to the chisq2d, contour and the TMarker!

using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareVALORotherExp()
{
  const bool draw_preliminary = true;
  TString str_preliminary = "T2K Run1-7c preliminary";
  
  const bool draw_reactor_band = false;
  const double x_low = 0.3;
  const double x_up = 0.73;
  const double y_low = 1.7;
  const double y_up = 3.6;

  const char *title_x = "sin^{2} #theta_{23}";
  const char *title_y = "|#Delta m^{2}_{32}| (10^{-3} eV^{2}/c^{4})";

  //////////////////////////// VALOR contours /////////////////////
  
  //const char *inname_nh = "../TECH_NOTE_14_3_16/DATA/RUN17/contours_sinsqth23Vsdmsq/fitmarg_joint_sinsqth23Vsdmsq_Run17_Data_react_nh_contour.root";
  //const char *inname_ih = "../TECH_NOTE_14_3_16/DATA/RUN17/contours_sinsqth23Vsdmsq/fitmarg_joint_sinsqth23Vsdmsq_Run17_Data_react_ih_contour.root";

  const char *inname_nh = "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17c_Data_react_nh_contour.root";
  const char *inname_ih = "../TECH_NOTE_14_3_16/DATA/RUN17c/contours_sinsqth23Vsdmsq_react/fitmarg_joint_sinsqth23Vsdmsq_Run17c_Data_react_ih_contour.root";

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

  //////////////////////////// Super-K contours /////////////////////

  const char *inname_sk_nh = "../../../inputs/OfficialContours/OtherExp/sk.atmospheric.201505/sk1234.atmospheric.201505.nh.root";
  const char *inname_sk_ih = "../../../inputs/OfficialContours/OtherExp/sk.atmospheric.201505/sk1234.atmospheric.201505.ih.root";
  const int colorSK = kBlue;

  //////////////////////////// NOvA contours /////////////////////////

  //const char *inname_nova_nh = "../../../inputs/OfficialContours/OtherExp/NOvAFAResults.root";
  const char *inname_nova_nh = "../../../inputs/OfficialContours/OtherExp/nova_surfaces_neutrino2016.root";
  const int colorNOvA = kGreen+1;

  //////////////////////////// MINOS contours /////////////////////////
  
  const char *inname_minos = "../../../inputs/OfficialContours/OtherExp/MINOSPlus_ThreeFlavour_PreliminaryContours_June2014.root";
  const int colorMINOS = kRed;

  //////////////////////////// IceCube contours /////////////////////////

  const char *inname_icecube = "../../../inputs/OfficialContours/OtherExp/sinsqth23Vsdmsq_IceCube2015_ORIGINAL_arXiv:1410.7227v2_nh.txt"; // taken from IceCube collaboration
  //const int colorIceCube = kCyan;
  const int colorIceCube = kOrange;
  
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
      SwapXY(g68_valor_nh[icont]);
#endif
  }
  for(int icont=0;icont<Ncont_90;icont++){
    g90_valor_nh[icont] = (TGraph*)filein_nh->Get(contour_90_name[icont]);
#ifdef SwapXY
    if(g90_valor_nh[icont] != NULL)
      SwapXY(g90_valor_nh[icont]);
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
      SwapXY(g68_valor_ih[icont]);
#endif
  }
  for(int icont=0;icont<Ncont_90;icont++){
    g90_valor_ih[icont] = (TGraph*)filein_ih->Get(contour_90_name[icont]);
#ifdef SwapXY
    if(g90_valor_ih[icont] != NULL)
      SwapXY(g90_valor_ih[icont]);
#endif
  }
  gbest_fit_global_ih = (TMarker*)filein_ih->Get(bestfitpt_name);
#ifdef SwapXY
  SwapXY(gbest_fit_global_ih);
#endif
 


  // READ SK CONTOURS
  TFile * infile_sk_nh = new TFile(inname_sk_nh,"READ");
  TGraph * sk1234_NH_S23M23_99CL = (TGraph*) infile_sk_nh->Get("sk1234_NH_S23M23_99CL");
  TGraph * sk1234_NH_S23M23_95CL = (TGraph*) infile_sk_nh->Get("sk1234_NH_S23M23_95CL");
  TGraph * sk1234_NH_S23M23_90CL = (TGraph*) infile_sk_nh->Get("sk1234_NH_S23M23_90CL");
  TGraph * sk1234_NH_S23M23_68CL = (TGraph*) infile_sk_nh->Get("sk1234_NH_S23M23_68CL");
  TGraph * base_NH_S23M23        = (TGraph*) infile_sk_nh->Get("base_NH_S23M23");
  TFile * infile_sk_ih = new TFile(inname_sk_ih,"READ");
  TGraph * sk1234_IH_S23M23_99CL = (TGraph*) infile_sk_ih->Get("sk1234_IH_S23M23_99CL");
  TGraph * sk1234_IH_S23M23_95CL = (TGraph*) infile_sk_ih->Get("sk1234_IH_S23M23_95CL");
  TGraph * sk1234_IH_S23M23_90CL = (TGraph*) infile_sk_ih->Get("sk1234_IH_S23M23_90CL");
  TGraph * sk1234_IH_S23M23_68CL = (TGraph*) infile_sk_ih->Get("sk1234_IH_S23M23_68CL");
  TGraph * base_IH_S23M23        = (TGraph*) infile_sk_ih->Get("base_IH_S23M23");
  
  // READ NOvA CONTOURS (only normal hierarchy)
  TFile * infile_nova_nh = new TFile(inname_nova_nh,"READ");
  TH1F *nova_s23m23_68CL_nh = (TH1F*) infile_nova_nh->Get("numu_dmsq32_vs_ssth23_nh"); 
  nova_s23m23_90CL_nh = (TH1F*) nova_s23m23_68CL_nh->Clone("nova_s23m23_nh_90CL");
  TH1F *nova_s23m23_68CL_ih = (TH1F*) infile_nova_nh->Get("numu_dmsq32_vs_ssth23_ih"); 
  nova_s23m23_90CL_ih = (TH1F*) nova_s23m23_68CL_ih->Clone("nova_s23m23_ih_90CL");
  //TGraph *nova_s23m23_68CL_nh = (TGraph*) infile_nova_nh->Get("FC68percentCL");
  //TGraph *nova_s23m23_90CL_nh = (TGraph*) infile_nova_nh->Get("FC90percentCL");
  double chisq_min = nova_s23m23_68CL_nh->GetBinContent(nova_s23m23_68CL_nh->GetMinimumBin());
  nova_s23m23_68CL_nh->SetContour(1);
  nova_s23m23_68CL_nh->SetContourLevel(0,chisq_min+2.3);
  nova_s23m23_90CL_nh->SetContour(1);
  nova_s23m23_90CL_nh->SetContourLevel(0,chisq_min+4.61);
  
  // READ MINOS CONTOURS
  TFile * infile_minos = new TFile(inname_minos,"READ");
  TGraph *minos_s23m23_68CL_nh = (TGraph*) infile_minos->Get("minosPlus_nu2014_prelim_contour_sinsq23_dmsq32_normal_68CL");
  TGraph *minos_s23m23_90CL_nh = (TGraph*) infile_minos->Get("minosPlus_nu2014_prelim_contour_sinsq23_dmsq32_normal_90CL");
  TGraph *minos_s23m23_68CL_ih = (TGraph*) infile_minos->Get("minosPlus_nu2014_prelim_contour_sinsq23_dmsq32_inverted_68CL");
  TGraph *minos_s23m23_90CL_ih = (TGraph*) infile_minos->Get("minosPlus_nu2014_prelim_contour_sinsq23_dmsq32_inverted_90CL");

  // READ ICECUBE CONTOURS
  TGraph *icecube_s23m23_90CL_nh;
  icecube_s23m23_90CL_nh = new TGraph(inname_icecube);
  
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
  // Super-K
  //  
  ScaleGraphAxis(sk1234_NH_S23M23_90CL,1.,1000.);
  ScaleGraphAxis(sk1234_NH_S23M23_68CL,1.,1000.);
  sk1234_NH_S23M23_90CL->SetLineStyle(1);
  sk1234_NH_S23M23_90CL->SetLineColor(colorSK);
  sk1234_NH_S23M23_90CL->SetLineWidth(2);
  sk1234_NH_S23M23_68CL->SetLineStyle(2);
  sk1234_NH_S23M23_68CL->SetLineColor(colorSK);
  sk1234_NH_S23M23_68CL->SetLineWidth(2);
  sk1234_NH_S23M23_90CL ->Draw("C"); 
  //sk1234_NH_S23M23_68CL ->Draw("Csame"); 
  
  //
  // NOvA
  //
  // //ScaleGraphAxis(nova_s23m23_90CL,1.,1000.);
  // //ScaleGraphAxis(nova_s23m23_68CL,1.,1000.);
  nova_s23m23_90CL_nh->SetLineStyle(1);
  nova_s23m23_90CL_nh->SetLineColor(colorNOvA);
  nova_s23m23_90CL_nh->SetLineWidth(2);
  nova_s23m23_68CL_nh->SetLineStyle(2);
  nova_s23m23_68CL_nh->SetLineColor(colorNOvA);
  nova_s23m23_68CL_nh->SetLineWidth(2);
  //nova_s23m23_68CL_nh->Draw("l same");
  //nova_s23m23_90CL_nh->Draw("l same");
  nova_s23m23_68CL_nh->Draw("cont3 same");
  nova_s23m23_90CL_nh->Draw("cont3 same");
  
  //
  // MINOS
  //
  ScaleGraphAxis(minos_s23m23_90CL_nh,1.,1000.);
  ScaleGraphAxis(minos_s23m23_68CL_nh,1.,1000.);
  minos_s23m23_90CL_nh->SetLineStyle(1);
  minos_s23m23_90CL_nh->SetLineColor(colorMINOS);
  minos_s23m23_90CL_nh->SetLineWidth(2);
  minos_s23m23_68CL_nh->SetLineStyle(2);
  minos_s23m23_68CL_nh->SetLineColor(colorMINOS);
  minos_s23m23_68CL_nh->SetLineWidth(2);
  minos_s23m23_68CL_nh->Draw("l same");
  minos_s23m23_90CL_nh->Draw("l same");

  //
  // IceCube
  //
  ScaleGraphAxis(icecube_s23m23_90CL_nh,1.,1000.);
  icecube_s23m23_90CL_nh->SetLineStyle(1);
  icecube_s23m23_90CL_nh->SetLineColor(colorIceCube);
  icecube_s23m23_90CL_nh->SetLineWidth(2);
  icecube_s23m23_90CL_nh->Draw("l same");
  
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
  TLatex *latSK = new TLatex(posX,posY,"Super-K");
  latSK->SetTextSize(0.06);
  latSK->SetTextColor(colorSK);
  latSK->Draw("same");
  TLatex *latNOvA = new TLatex(posX,posY,"NOvA (2016)");
  latNOvA->SetTextSize(0.06);
  latNOvA->SetTextColor(colorNOvA);
  latNOvA->Draw("same");
  TLatex *latMINOS = new TLatex(posX,posY,"MINOS+");
  latMINOS->SetTextSize(0.06);
  latMINOS->SetTextColor(colorMINOS);
  latMINOS->Draw("same");
  TLatex *latIceCube = new TLatex(posX,posY,"IceCube");
  latIceCube->SetTextSize(0.06);
  latIceCube->SetTextColor(colorIceCube);
  latIceCube->Draw("same");

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
  // Super-K
  //  
  ScaleGraphAxis(sk1234_IH_S23M23_90CL,1.,1000.);
  ScaleGraphAxis(sk1234_IH_S23M23_68CL,1.,1000.);
  sk1234_IH_S23M23_90CL->SetLineStyle(1);
  sk1234_IH_S23M23_90CL->SetLineColor(colorSK);
  sk1234_IH_S23M23_90CL->SetLineWidth(2);
  sk1234_IH_S23M23_68CL->SetLineStyle(2);
  sk1234_IH_S23M23_68CL->SetLineColor(colorSK);
  sk1234_IH_S23M23_68CL->SetLineWidth(2);
  sk1234_IH_S23M23_90CL ->Draw("C"); 
  //sk1234_IH_S23M23_68CL ->Draw("Csame"); 

  // //
  // // NOvA
  // //
  // //ScaleGraphAxis(nova_s23m23_90CL,1.,1000.);
  // //ScaleGraphAxis(nova_s23m23_68CL,1.,1000.);
  // nova_s23m23_90CL_ih->SetLineStyle(1);
  // nova_s23m23_90CL_ih->SetLineColor(colorNOvA);
  // nova_s23m23_90CL_ih->SetLineWidth(2);
  // nova_s23m23_68CL_ih->SetLineStyle(2);
  // nova_s23m23_68CL_ih->SetLineColor(colorNOvA);
  // nova_s23m23_68CL_ih->SetLineWidth(2);
  // nova_s23m23_68CL_ih->Draw("l same");
  // nova_s23m23_90CL_ih->Draw("l same");
  
  //
  // MINOS
  //
  ScaleGraphAxis(minos_s23m23_90CL_ih,1.,-1000.);
  ScaleGraphAxis(minos_s23m23_68CL_ih,1.,-1000.);
  minos_s23m23_90CL_ih->SetLineStyle(1);
  minos_s23m23_90CL_ih->SetLineColor(colorMINOS);
  minos_s23m23_90CL_ih->SetLineWidth(2);
  minos_s23m23_68CL_ih->SetLineStyle(2);
  minos_s23m23_68CL_ih->SetLineColor(colorMINOS);
  minos_s23m23_68CL_ih->SetLineWidth(2);
  minos_s23m23_68CL_ih->Draw("l same");
  minos_s23m23_90CL_ih->Draw("l same");  
  
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
  latSK->Draw("same");
  //latNOvA->Draw("same");
  latMINOS->Draw("same");
  
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
