
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

void compare_contours_dcp_with_Mach3_hierarchy
( 
 const int Nfiles = 2,
)
{
  
  char tag_1[Nfiles][300] = {"NH", "IH"};

  ///////////////// REACTORS ////////////////

  // // Asimov A 
  // char infilename_1[Nfiles][300] = // VALOR
  //   {
  //     "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_nh_A_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_ih_A_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
  //   }; 
  // char infilename_Mach3[Nfiles][300] = // Mach3
  //   { 
  //     "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov1_wRC_NH.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov1_wRC_IH.root"  // IH 
  //   };
  // char infilename_ptheta[Nfiles][300] = // P-Theta
  //   { 
  //     "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_newreactor_Set1.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_inv_newreactor_Set1.root"  // IH
  //   };
  
  // // Asimov B 
  // char infilename_1[Nfiles][300] = // VALOR
  //   {
  //     "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_nh_B_nexp_1_0.000_-999.000-job0-_0_0.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/valor/fit_fitmarg_joint_dcp_react_ih_B_nexp_1_0.000_-999.000-job0-_0_0.root"  // IH
  //   }; 
  // char infilename_Mach3[Nfiles][300] = // Mach3
  //   { 
  //     "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov2_wRC_NH.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov2_wRC_IH.root"  // IH 
  //   };
  // char infilename_ptheta[Nfiles][300] = // P-Theta
  //   { 
  //     "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_inv_newreactor_Set2.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_newreactor_Set2.root"  // IH
  //   };

  ///////////////// T2K ONLY ////////////////

  // // Asimov A
  // char infilename_1[Nfiles][300] = // VALOR
  //   {
  //     "../TECH_NOTE_14_3_16/SENSI/files/valor/toyfits_dcp_anh-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/valor/toyfits_dcp_aih-sinsqtheta[12_0.304,13_0.022,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root"  // IH
  //   }; 
  // char infilename_Mach3[Nfiles][300] = // Mach3
  //   { 
  //     "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov1_woRC_NH.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov1_woRC_IH.root"  // IH 
  //   };
  // char infilename_ptheta[Nfiles][300] = // P-Theta
  //   { 
  //     "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_Set1.root", // NH
  //     "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_inv_Set1.root"  // IH
  //   };

  // Asimov B
  char infilename_1[Nfiles][300] = // VALOR
    {
      "../TECH_NOTE_14_3_16/SENSI/files/valor/toyfits_dcp_cnh-sinsqtheta[12_0.304,13_0.022,23_0.450],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root", // NH
      "../TECH_NOTE_14_3_16/SENSI/files/valor/toyfits_dcp_cih-sinsqtheta[12_0.304,13_0.022,23_0.450],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0.root"  // IH
    }; 
  char infilename_Mach3[Nfiles][300] = // Mach3
    { 
      "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov2_woRC_NH.root", // NH
      "../TECH_NOTE_14_3_16/SENSI/files/mach3/dcp1D_asimov2_woRC_IH.root"  // IH 
    };
  char infilename_ptheta[Nfiles][300] = // P-Theta
    { 
      "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_Set2.root", // NH
      "../TECH_NOTE_14_3_16/SENSI/files/ptheta/Run1-6POT_dcp_inv_Set2.root"  // IH
    };


  char name[300];

  TH1D *hchisq1D_1[Nfiles];
  TH1D *hdeltachisq1D_1[Nfiles];
  TH1D *hchisq1D_Mach3_1[Nfiles];
  TH1D *hchisq1D_ptheta_1[Nfiles];

  double min_valor  = 999999.;
  double min_ptheta = 999999.;
  double min_mach3  = 999999.;

  // File ptheta

  TFile *filein_ptheta[Nfiles]; 
  for(int ifile=0;ifile<Nfiles; ifile++){
    filein_ptheta[ifile] = new TFile(infilename_ptheta[ifile],"READ");
    hchisq1D_ptheta_1[ifile] = (TH1D*)filein_ptheta[ifile]->Get("SA_dcp"); 
    
    double min_histo = hchisq1D_ptheta_1[ifile]->GetBinContent(hchisq1D_ptheta_1[ifile]->GetMinimumBin());
    if(min_histo < min_ptheta) min_ptheta = min_histo;
    //cout << "min_histo = " << min_histo << endl;
  }
  
  for(int ifile=0;ifile<Nfiles; ifile++){
    for(int ibin=0;ibin<hchisq1D_ptheta_1[ifile]->GetNbinsX();ibin++){
      double bincont = hchisq1D_ptheta_1[ifile]->GetBinContent(ibin+1);
      hchisq1D_ptheta_1[ifile]->SetBinContent(ibin+1,bincont-min_ptheta);
    }
  }


  
    
  // File Mach3

  TFile *filein_Mach3[Nfiles];
  for(int ifile=0;ifile<Nfiles;ifile++){
    filein_Mach3[ifile] = new TFile(infilename_Mach3[ifile],"READ");
    hchisq1D_Mach3_1[ifile] = (TH1D*)filein_Mach3[ifile]->Get("h_dcp_chi2"); 
    hchisq1D_Mach3_1[ifile]->SetFillColor(0);

    double min_histo = hchisq1D_Mach3_1[ifile]->GetBinContent(hchisq1D_Mach3_1[ifile]->GetMinimumBin());
    if(min_histo < min_mach3) min_mach3 = min_histo;
    //cout << "min_histo = " << min_histo << endl;
  }
  
  for(int ifile=0;ifile<Nfiles; ifile++){
    for(int ibin=0;ibin<hchisq1D_Mach3_1[ifile]->GetNbinsX();ibin++){
      double bincont = hchisq1D_Mach3_1[ifile]->GetBinContent(ibin+1);
      hchisq1D_Mach3_1[ifile]->SetBinContent(ibin+1,bincont-min_mach3);
    }
  }


  // File VALOR

  TFile *filein_1[Nfiles];
  for(int ifile=0;ifile<Nfiles;ifile++){
    filein_1[ifile] = new TFile(infilename_1[ifile],"READ");
    TH1D *histo_chisq;
    TH1D *histo_deltachisq;
    TTree *treein_1 = (TTree*)filein_1[ifile]->Get("fitnt");  
    treein_1->SetBranchAddress("hchisq1D", &histo_chisq );
    //treein_1->SetBranchAddress("hdeltachisq1D", &histo_chisq );
    treein_1->GetEntry(0);

    sprintf(name,"hchisq1D_%d",ifile);
    hchisq1D_1[ifile] = (TH1D*)histo_chisq->Clone(name);
  
    double min_histo = hchisq1D_1[ifile]->GetBinContent(hchisq1D_1[ifile]->GetMinimumBin());
    if(min_histo < min_valor) min_valor = min_histo;
  }
  //min_valor = hchisq1D_1[0]->GetBinContent(hchisq1D_1[0]->GetMinimumBin());
  
  for(int ifile=0;ifile<Nfiles; ifile++){
    for(int ibin=0;ibin<hchisq1D_1[ifile]->GetNbinsX();ibin++){
      
      double bincont = hchisq1D_1[ifile]->GetBinContent(ibin+1);
      hchisq1D_1[ifile]->SetBinContent(ibin+1,bincont-min_valor);
      
      //double min = hchisq1D_1[ifile]->GetBinContent(hchisq1D_1[ifile]->GetMinimumBin());
      //cout << min << endl;
      // cout << "bincont = " << bincont 
      // 	   << " - min = " << min_valor 
      // 	   << " - shift = " <<  bincont - min_valor << endl;
      //hchisq1D_1[ifile]->SetBinContent(ibin+1,bincont-min_valor);
    }
  }

  
  
  
  // Draw histograms

  TCanvas *cDeltaChisq = new TCanvas("cDeltaChisq","cDeltaChisq");
  //hchisq1D_ptheta_1[0]->Draw();
  //hchisq1D_ptheta_1[1]->Draw("same");
  
  for(int ifile=0;ifile<Nfiles;ifile++){
 
    hchisq1D_1[ifile]->SetLineWidth(2);
    hchisq1D_1[ifile]->SetLineColor(kBlue);
    hchisq1D_1[ifile]->SetLineStyle(ifile+1);
    if(ifile==0) hchisq1D_1[ifile]->Draw();
    else         hchisq1D_1[ifile]->Draw("same");

    hchisq1D_Mach3_1[ifile]->SetLineWidth(2);
    hchisq1D_Mach3_1[ifile]->SetLineColor(kRed);
    hchisq1D_Mach3_1[ifile]->SetLineStyle(ifile+1);
    hchisq1D_Mach3_1[ifile]->Draw("same");

    hchisq1D_ptheta_1[ifile]->SetLineWidth(2);
    hchisq1D_ptheta_1[ifile]->SetLineColor(kBlack);
    hchisq1D_ptheta_1[ifile]->SetLineStyle(ifile+1);
    hchisq1D_ptheta_1[ifile]->Draw("same");
  }
  
  TLegend *leg = new TLegend(0.7, 0.4, 0.9, 0.7);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  for(int ifile=0;ifile<1;ifile++){
    leg->AddEntry(hchisq1D_1[ifile],"VALOR","L");
    leg->AddEntry(hchisq1D_Mach3_1[ifile],"Mach3","L");
    leg->AddEntry(hchisq1D_ptheta_1[ifile],"p-#theta","L");
  }
  leg->Draw();

  TLegend *leg_hierarchy = new TLegend(0.2, 0.7, 0.4, 0.9);
  leg_hierarchy->SetFillColor(0);
  leg_hierarchy->SetBorderSize(0);
  for(int ifile=0;ifile<2;ifile++){
    leg_hierarchy->AddEntry(hchisq1D_ptheta_1[ifile],tag_1[ifile],"L");
  }
  leg_hierarchy->Draw();


  return;
}


