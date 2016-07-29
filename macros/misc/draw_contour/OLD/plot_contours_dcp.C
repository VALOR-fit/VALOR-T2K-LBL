
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;

void plot_contours_dcp
( 
 //char infilename_1[200] = "../marginalization/joint/fitmarg_deltacp_react_Run14_Asimov_A_NH_2Dtempl_OldPDG-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char infilename_1[200] = "../marginalization/joint/fitmarg_deltacp_react_Run14_Asimov_A_IH_2Dtempl-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char infilename_1[200] = "../marginalization/joint/fitmarg_deltaCP_t2kreact_Run16_Asimov_A_NH_2Dtempl-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char infilename_2[200] = "../marginalization/joint/fitmarg_deltaCP_t2kreact_Run16_Asimov_A_IH_2Dtempl-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",

 //char infilename_2[200] = "../../../../trunk/outputs/marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_OldJoint.root",
 //char infilename_2[200] = "../../../../trunk/outputs/marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_NewJoint.root",
 //char infilename_1[200] = "../../../../trunk/outputs/marginalization/joint_newtoys/fitmarg_joint_Run16_Asimov_A_FHCpRHC_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char infilename_2[200] = "../../../../trunk/outputs/marginalization/joint_newtoys/fitmarg_joint_Run16_Asimov_A_FHCpRHC_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",

 //char infilename_1[200] = "../marginalization/joint/fitmarg_deltaCP_Run14_t2konly_Asimov_A_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char infilename_2[200] = "../marginalization/joint/fitmarg_deltaCP_Run14_t2konly_Asimov_A_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char infilename_3[200] = "../marginalization/joint/fitmarg_deltaCP_Run16_t2konly_Asimov_A_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 //char infilename_4[200] = "../marginalization/joint/fitmarg_deltaCP_Run16_t2konly_Asimov_A_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",

 char infilename_1[200] = "../marginalization/joint/fitmarg_deltaCP_FullStat_t2kreact_Asimov_A_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 char infilename_2[200] = "../marginalization/joint/fitmarg_deltaCP_FullStat_t2kreact_Asimov_A_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 char infilename_3[200] = "../../../../trunk/outputs/marginalization/joint_newtoys/ASIMOV/fitmarg_deltaCP_FullStat_t2kreact_Asimov_A_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 char infilename_4[200] = "../../../../trunk/outputs/marginalization/joint_newtoys/ASIMOV/fitmarg_deltaCP_FullStat_t2kreact_Asimov_A_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",

 char tag_1[200] = "2D Templates - Full Stat - NH",
 char tag_2[200] = "2D Templates - Full Stat - IH",
 char tag_3[200] = "1D Templates - Full Stat - NH",
 char tag_4[200] = "1D Templates - Full Stat - IH",

 bool doGlobalMin  = true,

 const char tag[200] = "fitmarg_dcp_T2Kreact_DataRun14",
 bool SaveRootFile = false

)
{
  char name[200];

  TH1D *hchisq1D_1;
  TH1D *hchisq1D_2;
  TH1D *hchisq1D_3;
  TH1D *hchisq1D_4;
  TH1D *hdeltachisq1D_1;
  TH1D *hdeltachisq1D_2;
  TH1D *hdeltachisq1D_3;
  TH1D *hdeltachisq1D_4;
  TH1D *hdeltachisq1D_GlobMin_1;
  TH1D *hdeltachisq1D_GlobMin_2;
  TH1D *hdeltachisq1D_GlobMin_3;
  TH1D *hdeltachisq1D_GlobMin_4;
  
  // File 1
  TFile *filein_1 = new TFile(infilename_1,"READ");
  TTree *treein_1 = (TTree*)filein_1->Get("fitnt");  
  treein_1->SetBranchAddress("hchisq1D", &hchisq1D_1 );
  treein_1->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D_1 );
  treein_1->GetEntry(0);

  // File 2
  TFile *filein_2 = new TFile(infilename_2,"READ");
  TTree *treein_2 = (TTree*)filein_2->Get("fitnt");  
  treein_2->SetBranchAddress("hchisq1D", &hchisq1D_2 );
  treein_2->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D_2 );
  treein_2->GetEntry(0);

  // File 3
  TFile *filein_3 = new TFile(infilename_3,"READ");
  TTree *treein_3 = (TTree*)filein_3->Get("fitnt");  
  treein_3->SetBranchAddress("hchisq1D", &hchisq1D_3 );
  treein_3->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D_3 );
  treein_3->GetEntry(0);

  // File 4
  TFile *filein_4 = new TFile(infilename_4,"READ");
  TTree *treein_4 = (TTree*)filein_4->Get("fitnt");  
  treein_4->SetBranchAddress("hchisq1D", &hchisq1D_4 );
  treein_4->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D_4 );
  treein_4->GetEntry(0);

  
  // Shift the chisq wrt to the global minimum

  hdeltachisq1D_GlobMin_1 = (TH1D*)hchisq1D_1->Clone("hdeltachisq1D_GlobMin_1");
  hdeltachisq1D_GlobMin_1->Reset();
  hdeltachisq1D_GlobMin_2 = (TH1D*)hchisq1D_2->Clone("hdeltachisq1D_GlobMin_2");
  hdeltachisq1D_GlobMin_2->Reset();
  hdeltachisq1D_GlobMin_3 = (TH1D*)hchisq1D_3->Clone("hdeltachisq1D_GlobMin_3");
  hdeltachisq1D_GlobMin_3->Reset();
  hdeltachisq1D_GlobMin_4 = (TH1D*)hchisq1D_4->Clone("hdeltachisq1D_GlobMin_4");
  hdeltachisq1D_GlobMin_4->Reset();
  
  double ChisqMin_1 = hchisq1D_1->GetMinimum();
  double ChisqMin_2 = hchisq1D_2->GetMinimum();
  double ChisqMin_3 = hchisq1D_3->GetMinimum();
  double ChisqMin_4 = hchisq1D_4->GetMinimum();
  double ChisqMin_Glob = TMath::Min(ChisqMin_1,ChisqMin_2);
  ChisqMin_Glob = TMath::Min(ChisqMin_Glob,ChisqMin_3);
  ChisqMin_Glob = TMath::Min(ChisqMin_Glob,ChisqMin_4);
  double ChisqMin_Glob12 = TMath::Min(ChisqMin_1,ChisqMin_2);
  double ChisqMin_Glob34 = TMath::Min(ChisqMin_3,ChisqMin_4);
  cout << ChisqMin_1 << " " << ChisqMin_2 << " ";
  cout << ChisqMin_3 << " " << ChisqMin_4 << " ";
  cout << ChisqMin_Glob << " " << ChisqMin_Glob12 << " " << ChisqMin_Glob34 << endl;

  for(int i=0;i<hdeltachisq1D_GlobMin_1->GetNbinsX();i++){
    double cont = hchisq1D_1->GetBinContent(i+1); 
    if(doGlobalMin) hdeltachisq1D_GlobMin_1->SetBinContent(i+1,cont-ChisqMin_Glob12);
    else            hdeltachisq1D_GlobMin_1->SetBinContent(i+1,cont-ChisqMin_1);
  }
  for(int i=0;i<hdeltachisq1D_GlobMin_2->GetNbinsX();i++){
    double cont = hchisq1D_2->GetBinContent(i+1); 
    if(doGlobalMin) hdeltachisq1D_GlobMin_2->SetBinContent(i+1,cont-ChisqMin_Glob12);
    else            hdeltachisq1D_GlobMin_2->SetBinContent(i+1,cont-ChisqMin_2);
  }
  for(int i=0;i<hdeltachisq1D_GlobMin_3->GetNbinsX();i++){
    double cont = hchisq1D_3->GetBinContent(i+1); 
    if(doGlobalMin) hdeltachisq1D_GlobMin_3->SetBinContent(i+1,cont-ChisqMin_Glob34);
    else            hdeltachisq1D_GlobMin_3->SetBinContent(i+1,cont-ChisqMin_3);
  }
  for(int i=0;i<hdeltachisq1D_GlobMin_4->GetNbinsX();i++){
    double cont = hchisq1D_4->GetBinContent(i+1); 
    if(doGlobalMin) hdeltachisq1D_GlobMin_4->SetBinContent(i+1,cont-ChisqMin_Glob34);
    else            hdeltachisq1D_GlobMin_4->SetBinContent(i+1,cont-ChisqMin_4);
  }
  

  // Draw histograms

  TCanvas *cDeltaChisq = new TCanvas("cDeltaChisq","cDeltaChisq");
  hdeltachisq1D_GlobMin_1->SetLineColor(kRed);
  hdeltachisq1D_GlobMin_2->SetLineColor(kRed);
  hdeltachisq1D_GlobMin_3->SetLineColor(kBlack);
  hdeltachisq1D_GlobMin_4->SetLineColor(kBlack);
  hdeltachisq1D_GlobMin_1->SetLineWidth(2);
  hdeltachisq1D_GlobMin_2->SetLineWidth(2);
  hdeltachisq1D_GlobMin_3->SetLineWidth(2);
  hdeltachisq1D_GlobMin_4->SetLineWidth(2);
  hdeltachisq1D_GlobMin_1->SetLineStyle(1);
  hdeltachisq1D_GlobMin_2->SetLineStyle(2);
  hdeltachisq1D_GlobMin_3->SetLineStyle(1);
  hdeltachisq1D_GlobMin_4->SetLineStyle(2);
  hdeltachisq1D_GlobMin_1->GetYaxis()->SetRangeUser(0.,7.);
  hdeltachisq1D_GlobMin_1->Draw();
  hdeltachisq1D_GlobMin_2->Draw("same");
  hdeltachisq1D_GlobMin_3->Draw("same");
  hdeltachisq1D_GlobMin_4->Draw("same");

  TLegend *leg = new TLegend(0.3, 0.7, 0.5, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hdeltachisq1D_GlobMin_1,tag_1,"L");
  leg->AddEntry(hdeltachisq1D_GlobMin_2,tag_2,"L");
  leg->AddEntry(hdeltachisq1D_GlobMin_3,tag_3,"L");
  leg->AddEntry(hdeltachisq1D_GlobMin_4,tag_4,"L");
  leg->Draw();

  TLine *l1sigma = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),1.,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),1.);
  TLine *l90 = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),2.71,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),2.71);
  TLine *l2sigma = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),4.,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),4.);
  TLine *l3sigma = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),9.,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),9.);
  TLine *l4sigma = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),16.,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),16.);
  TLine *l5sigma = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),25.,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),25.);
  l1sigma->SetLineWidth(2);
  l90->SetLineWidth(2);
  l2sigma->SetLineWidth(2);
  l3sigma->SetLineWidth(2);
  l4sigma->SetLineWidth(2);
  l5sigma->SetLineWidth(2);
  l1sigma->SetLineColor(kBlue);
  l90->SetLineColor(kBlue);
  l2sigma->SetLineColor(kBlue);
  l3sigma->SetLineColor(kBlue);
  l4sigma->SetLineColor(kBlue);
  l5sigma->SetLineColor(kBlue);
  //l1sigma->Draw();
  //l90->Draw();
  l2sigma->Draw();
  l3sigma->Draw();
  l4sigma->Draw();
  l5sigma->Draw();
  
  TLatex *latLine1Sigma = new TLatex(-1.7,1.1,"1#sigma");
  TLatex *latLine90CL = new TLatex(-1.7,2.81,"90% CL");
  TLatex *latLine2Sigma = new TLatex(-1.7,4.1,"2#sigma");
  TLatex *latLine3Sigma = new TLatex(-1.7,9.1,"3#sigma");
  TLatex *latLine4Sigma = new TLatex(-1.7,16.1,"4#sigma");
  TLatex *latLine5Sigma = new TLatex(-1.7,25.1,"5#sigma");
  latLine1Sigma->SetTextColor(kBlue);
  latLine90CL->SetTextColor(kBlue);
  latLine2Sigma->SetTextColor(kBlue);
  latLine3Sigma->SetTextColor(kBlue);
  latLine4Sigma->SetTextColor(kBlue);
  latLine5Sigma->SetTextColor(kBlue);
  latLine1Sigma->SetTextSize(0.07);
  latLine90CL->SetTextSize(0.07);
  latLine2Sigma->SetTextSize(0.07);
  latLine3Sigma->SetTextSize(0.07);
  latLine4Sigma->SetTextSize(0.07);
  latLine5Sigma->SetTextSize(0.07);
  //latLine1Sigma->Draw("same");
  //latLine90CL->Draw("same");
  latLine2Sigma->Draw("same");
  latLine3Sigma->Draw("same");
  latLine4Sigma->Draw("same");
  latLine5Sigma->Draw("same");


  // Store in root file

  if(SaveRootFile){

    TFile* f2 = new TFile(Form("%s_contour.root",tag),"recreate");
    
    char name[200];
    sprintf(name,"hdeltachisq_%s",tag_1);
    hdeltachisq1D_GlobMin_1->Write(name);
    
    sprintf(name,"hdeltachisq_%s",tag_2);
    hdeltachisq1D_GlobMin_2->Write(name);

    sprintf(name,"hdeltachisq_%s",tag_3);
    hdeltachisq1D_GlobMin_3->Write(name);

    sprintf(name,"hdeltachisq_%s",tag_4);
    hdeltachisq1D_GlobMin_4->Write(name);

    f2->Close();
  }
  
  return;
}


