
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
 char infilename_1[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_OldJoint.root",
 char infilename_2[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_OldJoint.root",

 // char infilename_1[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_NewJoint.root",
 // char infilename_2[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_NewJoint.root",

 char tag_1[200] = "NH",
 char tag_2[200] = "IH",

 const char tag[200] = "fitmarg_dcp_T2Kreact_DataRun14",
 bool SaveRootFile = false

)
{
  char name[200];

  TH1D *hchisq1D_1;
  TH1D *hchisq1D_2;
  TH1D *hdeltachisq1D_1;
  TH1D *hdeltachisq1D_2;
  TH1D *hdeltachisq1D_GlobMin_1;
  TH1D *hdeltachisq1D_GlobMin_2;
  
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




  // Shift the chisq wrt to the global minimum

  hdeltachisq1D_GlobMin_1 = (TH1D*)hchisq1D_1->Clone("hdeltachisq1D_GlobMin_1");
  hdeltachisq1D_GlobMin_1->Reset();
  hdeltachisq1D_GlobMin_2 = (TH1D*)hchisq1D_2->Clone("hdeltachisq1D_GlobMin_2");
  hdeltachisq1D_GlobMin_2->Reset();
  
  double ChisqMin_1 = hchisq1D_1->GetMinimum();
  double ChisqMin_2 = hchisq1D_2->GetMinimum();
  double ChisqMin_Glob = TMath::Min(ChisqMin_1,ChisqMin_2);
  //cout << ChisqMin_1 << " " << ChisqMin_2 << " " << ChisqMin_Glob << endl;

  for(int i=0;i<hdeltachisq1D_GlobMin_1->GetNbinsX();i++){
    double cont = hchisq1D_1->GetBinContent(i+1); 
    hdeltachisq1D_GlobMin_1->SetBinContent(i+1,cont-ChisqMin_Glob);
  }
  for(int i=0;i<hdeltachisq1D_GlobMin_2->GetNbinsX();i++){
    double cont = hchisq1D_2->GetBinContent(i+1); 
    hdeltachisq1D_GlobMin_2->SetBinContent(i+1,cont-ChisqMin_Glob);
  }
  

  // Draw histograms

  TCanvas *cDeltaChisq = new TCanvas("cDeltaChisq","cDeltaChisq");
  hdeltachisq1D_GlobMin_1->SetLineColor(kRed);
  hdeltachisq1D_GlobMin_2->SetLineColor(kBlack);
  hdeltachisq1D_GlobMin_1->SetLineWidth(2);
  hdeltachisq1D_GlobMin_2->SetLineWidth(2);
  hdeltachisq1D_GlobMin_1->GetYaxis()->SetRangeUser(0.,7.);
  hdeltachisq1D_GlobMin_1->Draw();
  hdeltachisq1D_GlobMin_2->Draw("same");

  TLegend *leg = new TLegend(0.3, 0.7, 0.5, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hdeltachisq1D_GlobMin_1,tag_1,"L");
  leg->AddEntry(hdeltachisq1D_GlobMin_2,tag_2,"L");
  leg->Draw();


  TLine *l1sigma = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),1.,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),1.);
  TLine *l90 = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),2.71,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),2.71);
  TLine *l2sigma = new TLine(hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinLowEdge(1),4.,hdeltachisq1D_GlobMin_1->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin_1->GetNbinsX()),4.);
  l1sigma->SetLineWidth(2);
  l90->SetLineWidth(2);
  l2sigma->SetLineWidth(2);
  l1sigma->SetLineColor(kBlue);
  l90->SetLineColor(kBlue);
  l2sigma->SetLineColor(kBlue);
  l1sigma->Draw();
  l90->Draw();
  l2sigma->Draw();
  
  TLatex *latLine1Sigma = new TLatex(-1.7,1.1,"1#sigma");
  TLatex *latLine90CL = new TLatex(-1.7,2.81,"90% CL");
  TLatex *latLine2Sigma = new TLatex(-1.7,4.1,"2#sigma");
  latLine1Sigma->SetTextColor(kBlue);
  latLine90CL->SetTextColor(kBlue);
  latLine2Sigma->SetTextColor(kBlue);
  latLine1Sigma->SetTextSize(0.07);
  latLine90CL->SetTextSize(0.07);
  latLine2Sigma->SetTextSize(0.07);
  latLine1Sigma->Draw("same");
  latLine90CL->Draw("same");
  latLine2Sigma->Draw("same");


  // Store in root file

  if(SaveRootFile){

    TFile* f2 = new TFile(Form("%s_contour.root",tag),"recreate");
    
    char name[200];
    sprintf(name,"hdeltachisq_%s",tag_1);
    hdeltachisq1D_GlobMin_1->Write(name);
    
    sprintf(name,"hdeltachisq_%s",tag_2);
    hdeltachisq1D_GlobMin_2->Write(name);

    f2->Close();
  }
  
  return;
}


