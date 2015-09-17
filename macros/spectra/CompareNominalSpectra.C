
void CompareNominalSpectra
(
 //const char tag[200] = "1Re_OldFit_Asimov"
 //const char tag[200] = "1Re_Mach3_BestFit"
 const char tag[200] = "1Rmu_Mach3_Asimov",
 //const char tag[200] = "1Rmu_Mach3_BestFit",

 const double range_min = 0.,
 const double range_max = 2.
 )
{
  
  const int Nfiles = 2;
  const char *infilename[Nfiles] = {

    ///////////////////// 1Re //////////////////////
    
    // // 1Re - Run14 Mach3 Best-Fit
    // "/Users/davidesgalaberna/Desktop/Run14_Study/SPECTRA/new/sk_ereco_spectrum_1Re_prova-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2-1Re.root",
    // "/Users/davidesgalaberna/Desktop/Run14_Study/SPECTRA/old/sk_ereco_1Re_spectrum-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2.root"
    
    // // 1Re - Old fit Asimov
    // "/Users/davidesgalaberna/Desktop/Run14_Study/SPECTRA/new/sk_ereco_spectrum_1Re_prova-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2-1Re.root",
    // "/Users/davidesgalaberna/Desktop/Run14_Study/SPECTRA/old/sk_ereco_1Re_spectrum-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.root"

    ///////////////////// 1Rmu //////////////////////

    // 1Rmu - Run14 Mach3 Best-Fit
    "/Users/davidesgalaberna/Desktop/Run14_Study/SPECTRA/new/sk_ereco_spectrum_1Rmu_prova-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2-1Rmu.root",
    "/Users/davidesgalaberna/Desktop/Run14_Study/SPECTRA/old/sk_ereco_1Rmu_spectrum-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2.root"
    
    // // 1Rmu - Old fit Asimov
    // "/Users/davidesgalaberna/Desktop/Run14_Study/SPECTRA/new/sk_ereco_spectrum_1Rmu_prova-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2-1Rmu.root",
    // "/Users/davidesgalaberna/Desktop/Run14_Study/SPECTRA/old/sk_ereco_1Rmu_spectrum-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.root"

  };
  const char *tagfile[Nfiles] = {
    "New VALOR",
    "Old VALOR"
  };

  TFile *filein[Nfiles];
  TH1D *hErecoMode_Tot[Nfiles];
  TH1D *hErecoMode_NC[Nfiles];
  TH1D *hErecoMode_CCQE[Nfiles];
  TH1D *hErecoMode_CCnQE[Nfiles];
  TH1D *hErecoMode_NumuBeam[Nfiles]; 
  TH1D *hErecoMode_NumuBarBeam[Nfiles];
  TH1D *hErecoMode_NueBeam[Nfiles];
  TH1D *hErecoMode_NumuToNue[Nfiles];
  TH1D *hErecoMode_NumuBarToNueBar[Nfiles];
  TH1D *hErecoMode_NumuToNue[Nfiles];
  for(int ifile=0;ifile<Nfiles;ifile++){
    filein[ifile] = new TFile(infilename[ifile],"READ");
    hErecoMode_Tot[ifile]   = (TH1D*)filein[ifile]->Get("sk_1Retotal");
    if(!hErecoMode_Tot[ifile])  
      hErecoMode_Tot[ifile] = (TH1D*)filein[ifile]->Get("sk_1Rmutotal");
    hErecoMode_NC[ifile]   = (TH1D*)filein[ifile]->Get("hErecoMode_NC");
    hErecoMode_CCQE[ifile]  = (TH1D*)filein[ifile]->Get("hErecoMode_CCQE");
    hErecoMode_CCnQE[ifile] = (TH1D*)filein[ifile]->Get("hErecoMode_CCnQE");
    hErecoMode_NumuBeam[ifile] = (TH1D*)filein[ifile]->Get("hErecoMode_NumuBeam");
    hErecoMode_NumuBarBeam[ifile] = (TH1D*)filein[ifile]->Get("hErecoMode_NumuBarBeam");
    hErecoMode_NueBeam[ifile] = (TH1D*)filein[ifile]->Get("hErecoMode_NueBeam");
    hErecoMode_NumuToNue[ifile] = (TH1D*)filein[ifile]->Get("hErecoMode_NumuToNue");
    hErecoMode_NumuBarToNueBar[ifile] = (TH1D*)filein[ifile]->Get("hErecoMode_NumuBarToNueBar");
    hErecoMode_NumuToNue[ifile] = (TH1D*)filein[ifile]->Get("hErecoMode_NumuToNue");
  }
  
  // Draw histogram comparisons

  TCanvas *cTot = new TCanvas("cTot","cTot");
  DrawManyHisto(Nfiles,hErecoMode_Tot,tagfile,"Total",range_min,range_max);
  TCanvas *cNC = new TCanvas("cNC","cNC");
  DrawManyHisto(Nfiles,hErecoMode_NC,tagfile,"NC",range_min,range_max);
  TCanvas *cCCQE = new TCanvas("cCCQE","cCCQE");
  DrawManyHisto(Nfiles,hErecoMode_CCQE,tagfile,"CCQE",range_min,range_max);
  TCanvas *cCCnQE = new TCanvas("cCCnQE","cCCnQE");
  DrawManyHisto(Nfiles,hErecoMode_CCnQE,tagfile,"CCnQE",range_min,range_max);

  TCanvas *cNumuBeam = new TCanvas("cNumuBeam","cNumuBeam");
  DrawManyHisto(Nfiles,hErecoMode_NumuBeam,tagfile,"CC #nu_{#mu}",range_min,range_max);
  TCanvas *cNumuBarBeam = new TCanvas("cNumuBarBeam","cNumuBarBeam");
  DrawManyHisto(Nfiles,hErecoMode_NumuBarBeam,tagfile,"CC #bar{#nu}_{#mu}",range_min,range_max);
  TCanvas *cNueBeam = new TCanvas("cNueBeam","cNueBeam");
  DrawManyHisto(Nfiles,hErecoMode_NueBeam,tagfile,"CC #nu_{e}",range_min,range_max);
  TCanvas *cNumuToNue = new TCanvas("cNumuToNue","cNumuToNue");
  DrawManyHisto(Nfiles,hErecoMode_NumuToNue,tagfile,"CC #nu_{#mu} #rightarrow #nu_{e}",range_min,range_max);
  TCanvas *cNumuBarToNueBar = new TCanvas("cNumuBarToNueBar","cNumuBarToNueBar");
  DrawManyHisto(Nfiles,hErecoMode_NumuBarToNueBar,tagfile,"CC #bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}",range_min,range_max);
  TCanvas *cNumuToNue = new TCanvas("cNumuToNue","cNumuToNue");
  DrawManyHisto(Nfiles,hErecoMode_NumuToNue,tagfile,"CC #nu_{#mu} #rightarrow #nu_{e}",range_min,range_max);
 
  // Check
  for(int ifile=0;ifile<Nfiles;ifile++){
    double tot = hErecoMode_Tot[ifile]->Integral();
    double tot2 = 
      hErecoMode_NC[ifile]->Integral() +
      hErecoMode_CCQE[ifile]->Integral() +
      hErecoMode_CCnQE[ifile]->Integral();
    double tot3 = 
      hErecoMode_NC[ifile]->Integral() +
      hErecoMode_NumuBeam[ifile]->Integral() + 
      hErecoMode_NumuBarBeam[ifile]->Integral() + 
      hErecoMode_NueBeam[ifile]->Integral() + 
      hErecoMode_NumuBarToNueBar[ifile]->Integral() +
      hErecoMode_NumuToNue[ifile]->Integral();

    if( (tot-tot2)>1.E-5 ){
      cout << "file " << ifile << endl;
      cout << "tot = " << tot << " - tot2 = " << tot2 << endl;
      exit(1);
    }
    if( (tot-tot3)>1.E-5 ){
      cout << "file " << ifile << endl;
      cout << "tot = " << tot << " - tot3 = " << tot3 << endl;
      exit(1);
    }
  }

  // Print out # of events
  cout << endl;
  for(int ifile=0;ifile<Nfiles;ifile++){
    cout << tagfile[ifile] << ":" << endl
	 << " Tot # events = " << hErecoMode_Tot[ifile]->Integral() << endl
	 << " NC # events = " << hErecoMode_NC[ifile]->Integral() << endl
	 << " CCQE # events = " << hErecoMode_CCQE[ifile]->Integral() << endl
	 << " CCnQE # events = " << hErecoMode_CCnQE[ifile]->Integral() << endl
	 << " CC NumuBeam # events = " << hErecoMode_NumuBeam[ifile]->Integral() << endl 
	 << " CC NumuBarBeam # events = " << hErecoMode_NumuBarBeam[ifile]->Integral() << endl
	 << " CC Nue+NueBar Beam # events = " << hErecoMode_NueBeam[ifile]->Integral() << endl
	 << " CC NumuBarToNueBar # events = " << hErecoMode_NumuBarToNueBar[ifile]->Integral() << endl
	 << " CC NumuToNue # events = " << hErecoMode_NumuToNue[ifile]->Integral() 
	 << endl <<  endl;
  }

  // Print out relative difference wrt file ifile
  cout << endl;
  int ifile_ref = 0;
  for(int ifile=0;ifile<Nfiles;ifile++){
    if(ifile==ifile_ref) continue;
    cout << "Relative difference: " << tagfile[ifile_ref] << " / " << tagfile[ifile] << endl
      	 << " All events = " 
	 << (hErecoMode_Tot[ifile_ref]->Integral() - hErecoMode_Tot[ifile]->Integral()) / hErecoMode_Tot[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << " NC events = " 
	 << (hErecoMode_NC[ifile_ref]->Integral() - hErecoMode_NC[ifile]->Integral()) / hErecoMode_NC[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << " CCQE events = " 
	 << (hErecoMode_CCQE[ifile_ref]->Integral() - hErecoMode_CCQE[ifile]->Integral()) / hErecoMode_CCQE[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << " CCnQE events = "
	 << (hErecoMode_CCnQE[ifile_ref]->Integral() - hErecoMode_CCnQE[ifile]->Integral()) / hErecoMode_CCnQE[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << " CC NumuBeam events = " 
	 << (hErecoMode_NumuBeam[ifile_ref]->Integral() - hErecoMode_NumuBeam[ifile]->Integral()) / hErecoMode_NumuBeam[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << " CC NumuBarBeam events = "
	 << (hErecoMode_NumuBarBeam[ifile_ref]->Integral() - hErecoMode_NumuBarBeam[ifile]->Integral()) / hErecoMode_NumuBarBeam[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << " CC Nue+NueBar Beam events = "
	 << (hErecoMode_NueBeam[ifile_ref]->Integral() - hErecoMode_NueBeam[ifile]->Integral()) / hErecoMode_NueBeam[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << " CC NumuBarToNueBar events = " 
	 << (hErecoMode_NumuBarToNueBar[ifile_ref]->Integral() - hErecoMode_NumuBarToNueBar[ifile]->Integral()) / hErecoMode_NumuBarToNueBar[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << " CC NumuToNue events = " 
	 << (hErecoMode_NumuToNue[ifile_ref]->Integral() - hErecoMode_NumuToNue[ifile]->Integral()) / hErecoMode_NumuToNue[ifile_ref]->Integral() * 100.
	 << " %" << endl
	 << endl <<  endl;
  }

  // Save pdf
  char name[200];
  sprintf(name,"%s_Tot.pdf",tag);             cTot->Print(name);
  sprintf(name,"%s_NC.pdf",tag);              cNC->Print(name);
  sprintf(name,"%s_CCQE.pdf",tag);            cCCQE->Print(name);
  sprintf(name,"%s_CCnQE.pdf",tag);           cCCnQE->Print(name);
  sprintf(name,"%s_NumuBeam.pdf",tag);        cNumuBeam->Print(name);
  sprintf(name,"%s_NumuBarBeam.pdf",tag);     cNumuBarBeam->Print(name);
  sprintf(name,"%s_NueBeam.pdf",tag);         cNueBeam->Print(name);
  sprintf(name,"%s_NumuToNue.pdf",tag);       cNumuToNue->Print(name);
  sprintf(name,"%s_NumuBarToNueBar.pdf",tag); cNumuBarToNueBar->Print(name);
  sprintf(name,"%s_NumuToNue.pdf",tag);       cNumuToNue->Print(name);

  return;
}

void DrawManyHisto(const int Nfiles,TH1D **hErecoMode,const char **tagfile,const char *type,const double rmin,const double rmax){

  double max = 0;
  hErecoMode[0]->GetXaxis()->SetRangeUser(rmin,rmax);
  hErecoMode[0]->SetFillColor(0);
  hErecoMode[0]->SetLineColor(1);
  if(hErecoMode[0]->GetMaximum() >  max) 
    max = hErecoMode[0]->GetMaximum();
  for(int ifile=1;ifile<Nfiles;ifile++){
    hErecoMode[ifile]->SetFillColor(0);
    hErecoMode[ifile]->SetLineColor(ifile+1);
    if(hErecoMode[ifile]->GetMaximum() >  max) 
      max = hErecoMode[ifile]->GetMaximum();
  }

  hErecoMode[0]->SetMaximum(max*1.2);
  hErecoMode[0]->DrawCopy();
  for(int ifile=1;ifile<Nfiles;ifile++){
    hErecoMode[ifile]->SetMaximum(max*1.2);
    hErecoMode[ifile]->DrawCopy("same");
  }

  TLegend *ltot = new TLegend(0.65,0.45,0.8,0.65);
  ltot->SetFillColor(0);
  ltot->SetLineColor(0);
  ltot->SetTextFont(132);
  ltot->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    char name[200];
    sprintf(name,"%s",tagfile[ifile]);
    ltot->AddEntry(hErecoMode[ifile],name,"l");
  }
  ltot->Draw();
  double x = hErecoMode[0]->GetBinCenter(3);
  double y = hErecoMode[0]->GetMaximum() * 0.8;
  TLatex *lat = new TLatex(x,y,type);
  lat->SetTextSize(0.07);
  lat->Draw("same");
}
