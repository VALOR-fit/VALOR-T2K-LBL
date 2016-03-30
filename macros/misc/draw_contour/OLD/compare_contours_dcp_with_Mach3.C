
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

 ///////////////// REACTORS ////////////////

 // // Asimov A - NH - reactors
 // char infilename_1[300] = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_nh_A_nexp_1_-1.601_-999.000-job0-_0_0.root", // VALOR
 // char infilename_Mach3[300] = "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov1_wRC_NH.root", // Mach3 
 // char infilename_ptheta[300] = "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_newreactor_Set1.root", // p-theta

 // // Asimov A - IH - reactors
 // char infilename_1[300] = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_ih_A_nexp_1_-1.601_-999.000-job0-_0_0.root", // VALOR
 // char infilename_Mach3[300] = "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov1_wRC_IH.root", // Mach3 
 // char infilename_ptheta[300] = "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_inv_newreactor_Set1.root", // p-theta

 // // Asimov B - NH - reactors
 // char infilename_1[300] = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_nh_B_nexp_1_0.000_-999.000-job0-_0_0.root", // VALOR
 // char infilename_Mach3[300] = "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov2_wRC_NH.root", // Mach3 
 // char infilename_ptheta[300] = "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_newreactor_Set2.root", // p-theta

 // Asimov B - IH - reactors
 char infilename_1[300] = "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_ih_B_nexp_1_0.000_-999.000-job0-_0_0.root", // VALOR
 char infilename_Mach3[300] = "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov2_wRC_IH.root", // Mach3 
 char infilename_ptheta[300] = "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_inv_newreactor_Set2.root", // p-theta


 ///////////////// T2K ONLY ////////////////

 // // Asimov A - NH 
 // char infilename_1[300] = "../TECH_NOTE_14_3_16/SENSI/files/valor/toyfits_dcp_anh-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // VALOR
 // char infilename_Mach3[300] = "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov1_woRC_NH.root", // Mach3 
 // char infilename_ptheta[300] = "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_Set1.root", // p-theta

 // // Asimov A - IH - t2k only
 // char infilename_1[300] = "../TECH_NOTE_14_3_16/SENSI/files/valor/toyfits_dcp_aih-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // VALOR
 // char infilename_Mach3[300] = "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov1_woRC_IH.root", // Mach3 
 // char infilename_ptheta[300] = "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_inv_Set1.root", // p-theta
 
 // // Asimov B - NH - t2k only
 // char infilename_1[300] = "../TECH_NOTE_14_3_16/SENSI/files/valor/toyfits_dcp_cnh-sinsqtheta[12_0.304,13_0.022,23_0.450],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // VALOR
 // char infilename_Mach3[300] = "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov2_woRC_NH.root", // Mach3 
 // char infilename_ptheta[300] = "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_Set2.root", // p-theta

 // // Asimov B - IH - t2k only
 // char infilename_1[300] = "../TECH_NOTE_14_3_16/SENSI/files/valor/toyfits_dcp_cih-sinsqtheta[12_0.304,13_0.022,23_0.450],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // VALOR
 // char infilename_Mach3[300] = "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov2_woRC_IH.root", // Mach3 
 // char infilename_ptheta[300] = "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_inv_Set2.root", // p-theta
 
 char tag_1[300] = "NH",
)
{
  char name[300];

  TH1D *hchisq1D_1;
  TH1D *hdeltachisq1D_1;
  TH1D *hchisq1D_Mach3_1;
  TH1D *hchisq1D_ptheta_1;

  // File ptheta
  TFile *filein_ptheta = new TFile(infilename_ptheta,"READ");
  hchisq1D_ptheta_1 = (TH1D*)filein_ptheta->Get("SA_dcp"); 
  hchisq1D_ptheta_1->Draw();

 
  // // File Mach3 (from Canvas)
  // TFile *filein_Mach3 = new TFile(infilename_Mach3,"READ");
  // TCanvas *c1 = (TCanvas*)filein_Mach3->Get("b"); 
  // //hchisq1D_Mach3_1 = (TH1D*)c1->GetPrimitive("dcp_NH"); // data
  // hchisq1D_Mach3_1 = (TH1D*)c1->GetPrimitive("dcp_IH"); // data
  // //hchisq1D_Mach3_1 = (TH1D*)c1->GetPrimitive("h");
  // c1->Clear();
  // c1->Draw(); 
  // hchisq1D_Mach3_1->Scale(1./hchisq1D_Mach3_1->Integral());
  // double min=+9999999.;
  // // Turn likelihood into chisq
  // for(int i=0;i<hchisq1D_Mach3_1->GetNbinsX();i++){
  //   double cont = hchisq1D_Mach3_1->GetBinContent(i+1);
  //   double chisq = 0.;
  //   if(cont>0.){
  //     chisq = -2*log(cont);
  //     hchisq1D_Mach3_1->SetBinContent(i+1,chisq);
  //     if(chisq<min) min = chisq;
  //   }
  // }
  // for(int i=0;i<hchisq1D_Mach3_1->GetNbinsX();i++){
  //   double cont = hchisq1D_Mach3_1->GetBinContent(i+1);
  //   hchisq1D_Mach3_1->SetBinContent(i+1,cont-min);
  // }

  // File Mach3

  TFile *filein_Mach3 = new TFile(infilename_Mach3,"READ");
  hchisq1D_Mach3_1 = (TH1D*)filein_Mach3->Get("h_dcp_chi2"); 
  //hchisq1D_Mach3_1->Scale(1./hchisq1D_Mach3_1->Integral());
  //hchisq1D_Mach3_1->Draw();
  hchisq1D_Mach3_1->SetFillColor(0);


  // File VALOR

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


