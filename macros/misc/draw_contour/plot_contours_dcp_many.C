
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

void plot_contours_dcp_many
( 
 bool doGlobalMin  = true,

 //const char tag[200] = "fitmarg_dcp_T2Kreact_Comparison_Run16_Run17",
 //bool SaveRootFile = true
)
{

  // NIWG2015c Summer 2016 statistics
  const int Nfiles = 5;
  const int Nmh = 2;
  char *infilename[Nfiles][Nmh] = {
    // Current
    {
      "../marginalization/joint/NIWG2015c/CurrentStat/fit_fitmarg_joint_dcp_nh_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH 
      "../marginalization/joint/NIWG2015c/CurrentStat/fit_fitmarg_joint_dcp_ih_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
    }
    ,
    // Summer 2016
    {
      "../marginalization/joint/NIWG2015c/Summer2016/fit_fitmarg_joint_dcp_react_nh_Summer2016_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
      "../marginalization/joint/NIWG2015c/Summer2016/fit_fitmarg_joint_dcp_react_ih_Summer2016_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
    }
    ,
    // Summer 2016 100/0
    {
      "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_nh_Summer2016_100_0_nexp_1_-1.601_-999.000-job0-_0_0.root",      
      "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_ih_Summer2016_100_0_nexp_1_-1.601_-999.000-job0-_0_0.root",            
    }
    // Summer 2016 75/25
    {
      "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_nh_Summer2016_75_25_nexp_1_-1.601_-999.000-job0-_0_0.root",      
      "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_ih_Summer2016_75_25_nexp_1_-1.601_-999.000-job0-_0_0.root"
    }
    // Summer 2016 50/50
    {
      "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_nh_Summer2016_50_50_nexp_1_-1.601_-999.000-job0-_0_0.root",      
      "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_ih_Summer2016_50_50_nexp_1_-1.601_-999.000-job0-_0_0.root"
    }


  };



  
  char *latexname[Nfiles][Nmh] = {
    {"Current - NH", "Current - IH"}
    {"Summer Exp - NH", "Summer Exp - IH"}
    {"Summer 100/0 - NH", "Summer 100/0 - IH"}
    {"Summer 75/25 - NH", "Summer 75/25 - IH"}
    {"Summer 50/50 - NH", "Summer 50/50 - IH"}
  };
  
  char name[200];

  TFile *filein[Nfiles][Nmh];
  TTree *treein[Nfiles][Nmh];

  TH1F *hchisq1D[Nfiles][Nmh];
  TH1F *hdeltachisq1D[Nfiles][Nmh];
  TH1F *hdeltachisq1D_GlobMin[Nfiles][Nmh];
  
  TCanvas *cChisq; //[Nfiles][Nmh];

  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){      
    
      cout << endl;
      cout << "Get file " << infilename[ifile][imh] << endl;
      cout << endl;

      filein[ifile][imh] = new TFile(infilename[ifile][imh],"READ");
      if(!filein[ifile][imh]->IsOpen()){
	cerr << endl;
	cerr << "File " << infilename[ifile][imh] << " is not open!" << endl;
	cerr << endl;
	exit(1);
      }
      filein[ifile][imh]->cd();
      treein[ifile][imh] = (TTree*)filein[ifile][imh]->Get("fitnt");
      
      sprintf(name,"hchisq1D_file%d_mh%d",ifile,imh);
      hchisq1D[ifile][imh] = new TH1F(name,name,1,0,1);
      treein[ifile][imh]->SetBranchAddress("hchisq1D", &hchisq1D[ifile][imh] );
      
      sprintf(name,"hdeltachisq1D_file%d_mh%d",ifile,imh);
      hdeltachisq1D[ifile][imh] = new TH1F(name,name,1,0,1);
      treein[ifile][imh]->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D[ifile][imh] );
      
      // Get event
      treein[ifile][imh]->GetEntry(0);
      
      // Print chisq distributions
      cout << endl;
      cout << "Chisq file " << ifile << " - mh = " << imh << " :" << endl;
      for(int i=0;i<hchisq1D[ifile][imh]->GetNbinsX();i++){
	cout << hchisq1D[ifile][imh]->GetBinContent(i+1) << " ";
      }
      cout << endl;  
    
      sprintf(name,"hdeltachisq1D_GlobMin_file%d_mh%d",ifile,imh);
      hdeltachisq1D_GlobMin[ifile][imh] = (TH1F*)hchisq1D[ifile][imh]->Clone(name);
      hdeltachisq1D_GlobMin[ifile][imh]->Reset();
      
    } // end loop Nmh

    
    // Shift the chisq wrt to the global minimum
    
    double ChisqMin_0 = hchisq1D[ifile][0]->GetMinimum();
    double ChisqMin_1 = hchisq1D[ifile][1]->GetMinimum();
    double ChisqMin_Glob = ChisqMin_0;
    if( ChisqMin_1 < ChisqMin_Glob ) ChisqMin_Glob = ChisqMin_1;
    cout << endl;
    cout << "ChisqMin 0 = " << ChisqMin_0 << endl;
    cout << "ChisqMin 1 = " << ChisqMin_1 << endl;
    cout << "ChisqMin Glob = " << ChisqMin_Glob << endl;
    cout << endl;
    for(int imh=0;imh<Nmh;imh++){      
      double ChisqMin_Curr = hchisq1D[ifile][imh]->GetMinimum();
      for(int i=0;i<hdeltachisq1D_GlobMin[ifile][imh]->GetNbinsX();i++){
	double cont = hchisq1D[ifile][imh]->GetBinContent(i+1); 
	if(doGlobalMin) hdeltachisq1D_GlobMin[ifile][imh]->SetBinContent(i+1,cont-ChisqMin_Glob);
	else            hdeltachisq1D_GlobMin[ifile][imh]->SetBinContent(i+1,cont-ChisqMin_Curr);
      }
    }

  } // end loop Nfiles
  
  // Draw histograms
  
  // Chisq
 
  double Max_chisq = 0.;
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      if(hchisq1D[ifile][imh]->GetMaximum()>Max_chisq) 
	Max_chisq = hchisq1D[ifile][imh]->GetMaximum();
    }
  }
  
  sprintf(name,"cChisq1D");
  cChisq = new TCanvas(name,name);

  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      hchisq1D[ifile][imh]->SetMaximum(Max_chisq*1.2);
      hchisq1D[ifile][imh]->SetLineStyle(imh+1);
      hchisq1D[ifile][imh]->SetLineWidth(2);
      hchisq1D[ifile][imh]->SetLineColor(ifile+1);
      if((ifile+imh)==0) hchisq1D[ifile][imh]->Draw();
      else hchisq1D[ifile][imh]->Draw("same");
    }
  }

  TLegend *lchisq = new TLegend(0.6,0.25,0.8,0.65);
  lchisq->SetFillColor(0);
  lchisq->SetLineColor(0);
  lchisq->SetTextFont(132);
  lchisq->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      lchisq->AddEntry(hchisq1D[ifile][imh],latexname[ifile][imh],"l");
    }
  }
  lchisq->Draw();
  
  // Delta Chisq 

  double Max_deltachisq = 0.;
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      if(hdeltachisq1D[ifile][imh]->GetMaximum()>Max_deltachisq) 
	Max_deltachisq = hdeltachisq1D[ifile][imh]->GetMaximum();
    }
  }

  sprintf(name,"cDeltaChisq1D");
  cDeltaChisq = new TCanvas(name,name);

  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      hdeltachisq1D_GlobMin[ifile][imh]->SetMaximum(Max_deltachisq*1.2);
      hdeltachisq1D_GlobMin[ifile][imh]->SetLineStyle(imh+1);
      hdeltachisq1D_GlobMin[ifile][imh]->SetLineWidth(2);
      hdeltachisq1D_GlobMin[ifile][imh]->SetLineColor(ifile+1);
      if((ifile+imh)==0) hdeltachisq1D_GlobMin[ifile][imh]->Draw();
      else hdeltachisq1D_GlobMin[ifile][imh]->Draw("same");
    }
  }

  TLegend *ldeltachisq = new TLegend(0.6,0.25,0.8,0.65);
  ldeltachisq->SetFillColor(0);
  ldeltachisq->SetLineColor(0);
  ldeltachisq->SetTextFont(132);
  ldeltachisq->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      ldeltachisq->AddEntry(hdeltachisq1D_GlobMin[ifile][imh],latexname[ifile][imh],"l");
    }
  }
  ldeltachisq->Draw();
  
  TLine *l1sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),1.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),1.);
  TLine *l90 = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),2.71,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),2.71);
  TLine *l2sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),4.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),4.);
  TLine *l3sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),9.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),9.);
  TLine *l4sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),16.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),16.);
  TLine *l5sigma = new TLine(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1),25.,hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(hdeltachisq1D_GlobMin[0][0]->GetNbinsX()),25.);
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
  l1sigma->Draw();
  l90->Draw();
  l2sigma->Draw();
  //l3sigma->Draw();
  //l4sigma->Draw();
  //l5sigma->Draw();
  
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
  latLine1Sigma->Draw("same");
  latLine90CL->Draw("same");
  latLine2Sigma->Draw("same");
  //latLine3Sigma->Draw("same");
  //latLine4Sigma->Draw("same");
  //latLine5Sigma->Draw("same");



  // // Store in root file
  // if(SaveRootFile){
  //   TFile* f2 = new TFile(Form("%s_contour.root",tag),"recreate");    
  //   for(int ifile=0;ifile<Nfiles;ifile++){
  //     for(int imh=0;imh<Nmh;imh++){
  // 	//sprintf(name,"hdeltachisq_%i_%i",ifile);
  // 	hdeltachisq1D_GlobMin[ifile][imh]->Write();//(name);
  //     }
  //   }    
  //   f2->Close();
  // }
  
  return;
}


