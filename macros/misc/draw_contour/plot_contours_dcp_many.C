
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
 const bool doGlobalMin  = true,
 const bool same_mh_color = false, // use same color for NH and IH
 const bool same_mh_style = true, // use same style for NH and IH
 const bool draw_preliminary = true,

 const bool draw_1sigmaCL = false,
 const bool draw_90CL     = false,
 const bool draw_2sigmaCL = false,
 const bool draw_3sigmaCL = true,
 const bool draw_4sigmaCL = false,
 const bool draw_5sigmaCL = true,

 const bool draw_reactors = false
)
{

  // NIWG2015c Summer 2016 statistics
  const int Nfiles = 1;
  const int Nmh = 2;
  char *infilename[Nfiles][Nmh] = {

    // Run16 Data
    // {
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_dcp_react_nh_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_dcp_react_ih_nexp_1_-1.601_-999.000-job0-_0_0.root" // IH
    // },
    // {
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_sinsqth13_data_t2konly_nh_nexp_1_0.022_-999.000-job0-_0_0.root", 
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_sinsqth13_data_t2konly_ih_nexp_1_0.022_-999.000-job0-_0_0.root" 
    // },
    // {
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_sinsqth23_data_react_nh_nexp_1_0.528_-999.000-job0-_0_0.root", 
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_sinsqth23_data_react_ih_nexp_1_0.528_-999.000-job0-_0_0.root" 
    // },
    // {
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_dmsq_data_react_nh_nexp_1_0.003_-999.000-job0-_0_0.root",
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_dmsq_data_react_ih_nexp_1_0.003_-999.000-job0-_0_0.root"
    // },
    // {
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_dcp_t2konly_nh_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_dcp_t2konly_ih_nexp_1_-1.601_-999.000-job0-_0_0.root" // IH
    // },
   
    // // Run16 Data (FHC-only)
    // {
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_dcp_data_react_nh_fhconly_nexp_1_-1.601_-999.000-job0-_0_0.root",
    //   "../marginalization/joint/NIWG2015c/Data/fit_fitmarg_joint_dcp_data_react_ih_fhconly_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // },
    
    // // Current
    // {
    //   "../marginalization/joint/NIWG2015c/CurrentStat/fit_fitmarg_joint_dcp_nh_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH 
    //   "../marginalization/joint/NIWG2015c/CurrentStat/fit_fitmarg_joint_dcp_ih_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
    //   //"../marginalization/joint/NIWG2015c/CurrentStat/fit_fitmarg_joint_dcp_t2konly_nh_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH 
    //   //"../marginalization/joint/NIWG2015c/CurrentStat/fit_fitmarg_joint_dcp_t2konly_ih_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
    // }
    // ,
    // // Summer 2016
    // {
    //   "../marginalization/joint/NIWG2015c/Summer2016/fit_fitmarg_joint_dcp_react_nh_Summer2016_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
    //   "../marginalization/joint/NIWG2015c/Summer2016/fit_fitmarg_joint_dcp_react_ih_Summer2016_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
    //   //"../marginalization/joint/NIWG2015c/Summer2016/fit_fitmarg_joint_dcp_t2konly_nh_Summer2016_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
    //   //"../marginalization/joint/NIWG2015c/Summer2016/fit_fitmarg_joint_dcp_t2konly_ih_Summer2016_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
    // }
    // ,
    // // Summer 2016 100/0
    // {
    //   "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_nh_Summer2016_100_0_nexp_1_-1.601_-999.000-job0-_0_0.root",      
    //   "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_ih_Summer2016_100_0_nexp_1_-1.601_-999.000-job0-_0_0.root",            
    //   //"../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_t2konly_nh_Summer2016_100_0_nexp_1_-1.601_-999.000-job0-_0_0.root",      
    //   //"../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_t2konly_ih_Summer2016_100_0_nexp_1_-1.601_-999.000-job0-_0_0.root",            
    // }
    // // Summer 2016 75/25
    // {
    //   "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_nh_Summer2016_75_25_nexp_1_-1.601_-999.000-job0-_0_0.root",      
    //   "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_ih_Summer2016_75_25_nexp_1_-1.601_-999.000-job0-_0_0.root"
    //   //"../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_t2konly_nh_Summer2016_75_25_nexp_1_-1.601_-999.000-job0-_0_0.root",      
    //   //"../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_t2konly_ih_Summer2016_75_25_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // }
    // // Summer 2016 50/50
    // {
    //   "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_nh_Summer2016_50_50_nexp_1_-1.601_-999.000-job0-_0_0.root",      
    //   "../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_react_ih_Summer2016_50_50_nexp_1_-1.601_-999.000-job0-_0_0.root"
    //   //"../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_t2konly_nh_Summer2016_50_50_nexp_1_-1.601_-999.000-job0-_0_0.root",      
    //   //"../marginalization/joint/NIWG2015c/Summer2016_FHCRHCother/fit_fitmarg_joint_dcp_t2konly_ih_Summer2016_50_50_nexp_1_-1.601_-999.000-job0-_0_0.root"
    // }

    // High Statistics (125x10^20 POT for both Nu and NuBar) 
    {
      "../marginalization/joint/NIWG2015c/HighStat/fit_fitmarg_joint_dcp_react_nh_125nu125nubarPOT_nexp_1_-1.601_-999.000-job0-_0_0.root", // NH
      "../marginalization/joint/NIWG2015c/HighStat/fit_fitmarg_joint_dcp_react_ih_125nu125nubarPOT_nexp_1_-1.601_-999.000-job0-_0_0.root"  // IH
    }
  };
  
  char *latexname[Nfiles][Nmh] = {
    {"Normal Hierarchy","Inverted Hierarchy"},
    // {"NH","IH"},
    // {"NH all","IH all"},
    // {"NH FHC-only","IH FHC-only"}
    // {"Current - NH", "Current - IH"}
    // {"Summer Exp - NH", "Summer Exp - IH"}
    // {"Summer 100/0 - NH", "Summer 100/0 - IH"}
    // {"Summer 75/25 - NH", "Summer 75/25 - IH"}
    // {"Summer 50/50 - NH", "Summer 50/50 - IH"}
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
    
    int XminBin_0 = hchisq1D[ifile][0]->GetMinimumBin();
    double Xmin_0 = hchisq1D[ifile][0]->GetBinCenter(XminBin_0);
    int XminBin_1 = hchisq1D[ifile][1]->GetMinimumBin();
    double Xmin_1 = hchisq1D[ifile][1]->GetBinCenter(XminBin_1);

    cout << endl;
    cout << "Minimum for " << latexname[ifile][0] << " :" << endl;
    cout << " - X = " << Xmin_0 << endl;
    cout << " - ChisqMin " << latexname[ifile][0] << " = " << ChisqMin_0 << endl;
    cout << "Minimum for " << latexname[ifile][1] << " :" << endl;
    cout << " - X = " << Xmin_1 << endl;
    cout << " - ChisqMin " << latexname[ifile][1] << " = " << ChisqMin_1 << endl;
    cout << endl;
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
  

  // Init histogram

  int NbinsX = hdeltachisq1D_GlobMin[0][0]->GetNbinsX(); 
  int NbinsY = hdeltachisq1D_GlobMin[0][0]->GetNbinsY();
  double lowedgeX = hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(1);
  double upedgeX = hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinUpEdge(NbinsX);
  char *TitleX = hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetTitle();
  char *TitleY = hdeltachisq1D_GlobMin[0][0]->GetYaxis()->GetTitle();
  TH1D *hInit = new TH1D("hInit","hInit",NbinsX,lowedgeX,upedgeX);
  hInit->GetXaxis()->SetTitle(TitleX);
  hInit->GetYaxis()->SetTitle(TitleY);

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
  
  hInit->SetMaximum(Max_chisq*1.2);
  hInit->DrawCopy();

  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      hchisq1D[ifile][imh]->SetMaximum(Max_chisq*1.2);
      if(same_mh_style) hchisq1D[ifile][imh]->SetLineStyle(1);
      else hchisq1D[ifile][imh]->SetLineStyle(imh+1);
      hchisq1D[ifile][imh]->SetLineWidth(2);
      if(same_mh_color) hchisq1D[ifile][imh]->SetLineColor(ifile+1);
      else hchisq1D[ifile][imh]->SetLineColor(ifile+imh+1);
      if((ifile+imh)==0) hchisq1D[ifile][imh]->Draw("same");
      else hchisq1D[ifile][imh]->Draw("same");
    }
  }

  TLegend *lchisq = new TLegend(0.6,0.45,0.8,0.65);
  lchisq->SetFillColor(0);
  lchisq->SetLineColor(0);
  lchisq->SetTextFont(132);
  lchisq->SetTextSize(0.07);
  for(int ifile=0;ifile<Nfiles;ifile++){
    //for(int imh=0;imh<1;imh++){
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

  hInit->SetMaximum(Max_deltachisq*1.2);
  hInit->DrawCopy();

  // Draw reactors box
  const double sinsqth13_nom = 0.0219;
  const double sinsqth13_err = 0.0012;
  double X1 = sinsqth13_nom - sinsqth13_err;
  double X2 = sinsqth13_nom + sinsqth13_err;
  double Y1 = hInit->GetMinimum();
  //double Y2 = hInit->GetMaximum();
  double Y2 = 45.;
  TBox *breactors = new TBox(X1,Y1,X2,Y2);
  breactors->SetFillColor(kYellow+1);
  if(draw_reactors) breactors->Draw("");


  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int imh=0;imh<Nmh;imh++){
      hdeltachisq1D_GlobMin[ifile][imh]->SetMaximum(Max_deltachisq*1.2);
      if(same_mh_style) hdeltachisq1D_GlobMin[ifile][imh]->SetLineStyle(1);
      else hdeltachisq1D_GlobMin[ifile][imh]->SetLineStyle(imh+1);
      hdeltachisq1D_GlobMin[ifile][imh]->SetLineWidth(2);
      if(same_mh_color) hdeltachisq1D_GlobMin[ifile][imh]->SetLineColor(ifile+1);
      else hdeltachisq1D_GlobMin[ifile][imh]->SetLineColor(ifile+imh+1);
      if((ifile+imh)==0) hdeltachisq1D_GlobMin[ifile][imh]->Draw("same");
      else hdeltachisq1D_GlobMin[ifile][imh]->Draw("same");
    }
  }

  TLegend *ldeltachisq = new TLegend(0.6,0.45,0.8,0.65);
  ldeltachisq->SetFillColor(0);
  ldeltachisq->SetLineColor(0);
  ldeltachisq->SetTextFont(132);
  ldeltachisq->SetTextSize(0.07);
  for(int ifile=0;ifile<Nfiles;ifile++){
    //for(int imh=0;imh<1;imh++){
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
  if(draw_1sigmaCL) l1sigma->Draw();
  if(draw_90CL)     l90->Draw();
  if(draw_2sigmaCL) l2sigma->Draw();
  if(draw_3sigmaCL) l3sigma->Draw();
  if(draw_4sigmaCL) l4sigma->Draw();
  if(draw_5sigmaCL) l5sigma->Draw();
  
  TLatex *latLine1Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),1.1,"1#sigma");
  TLatex *latLine90CL = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),2.81,"90% CL");
  TLatex *latLine2Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),4.1,"2#sigma");
  TLatex *latLine3Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),9.1,"3#sigma");
  TLatex *latLine4Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),16.1,"4#sigma");
  TLatex *latLine5Sigma = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),25.1,"5#sigma");
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
  if(draw_1sigmaCL) latLine1Sigma->Draw("same");
  if(draw_90CL)     latLine90CL->Draw("same");
  if(draw_2sigmaCL) latLine2Sigma->Draw("same");
  if(draw_3sigmaCL) latLine3Sigma->Draw("same");
  if(draw_4sigmaCL) latLine4Sigma->Draw("same");
  if(draw_5sigmaCL) latLine5Sigma->Draw("same");


  TLatex *latPreliminary = new TLatex(hdeltachisq1D_GlobMin[0][0]->GetXaxis()->GetBinLowEdge(5),1.1,"PRELIMINARY");
  if(draw_preliminary) latPreliminary->Draw("same");


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


