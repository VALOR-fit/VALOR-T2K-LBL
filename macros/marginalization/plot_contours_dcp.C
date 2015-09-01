
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
 // char infilename_1[200] = "../marginalization/joint/fitmarg_joint_dcp-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_1Dtempl.root",
 // char infilename_2[200] = "../marginalization/joint/fitmarg_joint_dcp-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_2Dtempl_ErecoVsTheta.root",
 // char infilename_3[200] = "../marginalization/joint/fitmarg_dcp_NH_T2Kreact_Run1to6_DataSetA_ErecoVsQ2-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 // char tag_1[200] = "E_{reco} (NH)",
 // char tag_2[200] = "E_{reco} Vs #theta (NH)",
 // char tag_3[200] = "E_{reco} Vs Q^{2} (NH)"

 char infilename_1[200] = "../marginalization/joint/fitmarg_joint_dcp_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_1Dtempl.root",
 char infilename_2[200] = "../marginalization/joint/fitmarg_joint_dcp_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_2Dtempl_ErecoVsTheta.root",
 char infilename_3[200] = "../marginalization/joint/fitmarg_dcp_IH_T2Kreact_Run1to6_DataSetA_ErecoVsQ2-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root",
 char tag_1[200] = "E_{reco} (NH)",
 char tag_2[200] = "E_{reco} Vs #theta (IH)",
 char tag_3[200] = "E_{reco} Vs Q^{2} (IH)"

)
{
  char name[200];

  TH1D *hchisq1D_1;
  TH1D *hchisq1D_2;
  TH1D *hchisq1D_3;
  TH1D *hdeltachisq1D_1;
  TH1D *hdeltachisq1D_2;
  TH1D *hdeltachisq1D_3;
  
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
 
  // Draw histograms

  TCanvas *cDeltaChisq = new TCanvas("cDeltaChisq","cDeltaChisq");
  hdeltachisq1D_1->SetLineColor(kRed);
  hdeltachisq1D_2->SetLineColor(kBlack);
  hdeltachisq1D_3->SetLineColor(kGreen);
  hdeltachisq1D_1->SetLineWidth(2);
  hdeltachisq1D_2->SetLineWidth(2);
  hdeltachisq1D_3->SetLineWidth(2);
  hdeltachisq1D_1->Draw();
  hdeltachisq1D_2->Draw("same");
  hdeltachisq1D_3->Draw("same");

  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hdeltachisq1D_1,tag_1,"L");
  leg->AddEntry(hdeltachisq1D_2,tag_2,"L");
  leg->AddEntry(hdeltachisq1D_3,tag_3,"L");
  leg->Draw();



  return;
}


