
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

void compare_contours_dcp_with_Mach3
( 
 //char infilename_1[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_OldJoint.root",
 //char infilename_1[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_OldJoint.root",

 //char infilename_1[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_NewJoint.root",
 //char infilename_1[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Asimov_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_NewJoint.root",

 //char infilename_1[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Data_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_OldJoint.root",
 char infilename_1[200] = "../marginalization/joint_newtoys/fitmarg_dcp_T2Kreact_Run1to4_Data_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0_OldJoint.root",

 //char infilename_Mach3[200] = "../marginalization/Mach3/mach3_dcp_chisq/contours_1d_dcp.root",
 //char infilename_Mach3[200] = "../marginalization/Mach3/mach3_dcp_chisq/contours_1d_dcp_newRC.root",
 char infilename_Mach3[200] = "../marginalization/Mach3/mach3_dcp_chisq/compareold_1d_dcp.root", // from Canvas (Data)
 //char infilename_Mach3[200] = "../marginalization/Mach3/dcp_1d_compareFHCjoint.root", // from Canvas (Asimov), it should be new reactor constraint
 
 //char infilename_ptheta[200] = "../marginalization/ptheta/Run1-4POT_dcp_reactor_point1.root", // new reactor constraint (NH)
 char infilename_ptheta[200] = "../marginalization/ptheta/Run1-4POT_dcp_inv_reactor_point1.root", // new reactor constraint (IH)
 //char infilename_ptheta[200] = "../marginalization/ptheta/Run1-4POT_dcp_inv_oldreactor_point1.root", // old reactor constraint (NH, wrong label)
 //char infilename_ptheta[200] = "../marginalization/ptheta/Run1-4POT_dcp_oldreactor_point1.root", // old reactor constraint (IH, wrong label)
 
 char tag_1[200] = "NH",
)
{
  char name[200];

  TH1D *hchisq1D_1;
  TH1D *hdeltachisq1D_1;
  TH1D *hchisq1D_Mach3_1;
  TH1D *hchisq1D_ptheta_1;

  // File ptheta
  TFile *filein_ptheta = new TFile(infilename_ptheta,"READ");
  hchisq1D_ptheta_1 = (TH1D*)filein_ptheta->Get("SA_dcp"); 
  hchisq1D_ptheta_1->Draw();


  // File Mach3 (from Canvas)
  TFile *filein_Mach3 = new TFile(infilename_Mach3,"READ");
  TCanvas *c1 = (TCanvas*)filein_Mach3->Get("b"); 
  //hchisq1D_Mach3_1 = (TH1D*)c1->GetPrimitive("dcp_NH"); // data
  hchisq1D_Mach3_1 = (TH1D*)c1->GetPrimitive("dcp_IH"); // data
  //hchisq1D_Mach3_1 = (TH1D*)c1->GetPrimitive("h");
  c1->Clear();
  c1->Draw(); 
  hchisq1D_Mach3_1->Scale(1./hchisq1D_Mach3_1->Integral());
 
  // // File Mach3
  // TFile *filein_Mach3 = new TFile(infilename_Mach3,"READ");
  // hchisq1D_Mach3_1 = (TH1D*)filein_Mach3->Get("dcp_NH"); 
  // hchisq1D_Mach3_1->Scale(1./hchisq1D_Mach3_1->Integral());
  // hchisq1D_Mach3_1->Draw();
  // return;

  double min=+9999999.;
  // Turn likelihood into chisq
  for(int i=0;i<hchisq1D_Mach3_1->GetNbinsX();i++){
    double cont = hchisq1D_Mach3_1->GetBinContent(i+1);
    double chisq = 0.;
    if(cont>0.){
      chisq = -2*log(cont);
      hchisq1D_Mach3_1->SetBinContent(i+1,chisq);
      if(chisq<min) min = chisq;
    }
  }
  for(int i=0;i<hchisq1D_Mach3_1->GetNbinsX();i++){
    double cont = hchisq1D_Mach3_1->GetBinContent(i+1);
    hchisq1D_Mach3_1->SetBinContent(i+1,cont-min);
  }

  // File 1
  TFile *filein_1 = new TFile(infilename_1,"READ");
  TTree *treein_1 = (TTree*)filein_1->Get("fitnt");  
  treein_1->SetBranchAddress("hchisq1D", &hchisq1D_1 );
  treein_1->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D_1 );
  treein_1->GetEntry(0);

  // Draw histograms

  TCanvas *cDeltaChisq = new TCanvas("cDeltaChisq","cDeltaChisq");
  hdeltachisq1D_1->SetLineColor(kBlue);
  hdeltachisq1D_1->SetLineWidth(2);
  hdeltachisq1D_1->GetYaxis()->SetRangeUser(0.,7.);
  hdeltachisq1D_1->Draw();
  
  hchisq1D_Mach3_1->SetLineColor(kRed);
  hchisq1D_Mach3_1->SetLineWidth(2);
  hchisq1D_Mach3_1->Draw("same");
  
  hchisq1D_ptheta_1->SetLineColor(kBlack);
  hchisq1D_ptheta_1->SetLineWidth(2);
  hchisq1D_ptheta_1->Draw("same");

  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hdeltachisq1D_1,"VALOR","L");
  leg->AddEntry(hchisq1D_Mach3_1,"Mach3","L");
  leg->AddEntry(hchisq1D_ptheta_1,"p-#theta","L");
  leg->Draw();



  return;
}


