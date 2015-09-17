
void CompareMargSpectra()
{
  // VALOR histograms
  const int Nfiles = 2; //3;
  const char *infilename[Nfiles] = {
    "../marginalization/joint_newtoys/Run14_Study/spectrum_bestfit_sinsqth13Vsdcp_T2Konly_Run14_DataNH.root",
    //"../marginalization/joint_newtoys/Run14_Study/spectrum_bestfit_sinsqth13Vsdcp_T2Konly_Run14_DataNH_FixMEC.root",
    //"../marginalization/joint_newtoys/Run14_Study/spectrum_bestfit_sinsqth13Vsdcp_T2Konly_Run14_DataNH_StatOnly.root"

    "../marginalization/joint_newtoys/Run14_Study/spectrum_bestfit_sinsqth13Vsdcp_T2Konly_Run14_DataIH.root"
    //"../marginalization/joint_newtoys/Run14_Study/spectrum_bestfit_sinsqth13Vsdcp_T2Konly_Run14_DataIH_FixMEC.root",
    //"../marginalization/joint_newtoys/Run14_Study/spectrum_bestfit_sinsqth13Vsdcp_T2Konly_Run14_DataIH_StatOnly.root"
  };  
  const char *HistoName[Nfiles] = {
    "h_bestfit_marg",
    "h_bestfit_marg"
    //"h_bestfit_marg"
  };
  const char *AnalysisMode[Nfiles] = {
    "Stat+Syst - T2K only - NH",
    //"Fix MEC - T2K only - NH",
    //"Stat-only - T2K only - NH"

    "Stat+Syst - T2K only - IH"
    //"Fix MEC - T2K only - IH",
    //"Stat-only - T2K only - IH"
  };

  // Mach3 histograms
  const int Nfiles_mach3 = 2;
  const char *infilename_mach3[Nfiles_mach3] = {
    "/Users/davidesgalaberna/Desktop/Run14_Study/kirsty_files/new_4th13_bestfit_rates.root",
    "/Users/davidesgalaberna/Desktop/Run14_Study/kirsty_files/jointfit_woRC_4xth13step_erecspectra.root"
  };
  const char *HistoName_mach3[Nfiles_mach3] = {
    "nue_nominal_hist", 
    "nue_hist_1d"
  };
  const char *AnalysisMode_mach3[Nfiles_mach3] = {
    "Mach3 - T2K only - profiled",
    "Mach3 - T2K only - marginalized"
  };

  //

  TH1D *h_spectrum_data;
  TH1D *h_spectrum_bestfit[Nfiles];
  TFile *fin[Nfiles];
  TFile *fin_mach3[Nfiles_mach3];
  TH1D *h_spectrum_bestfit_mach3[Nfiles_mach3];
  
  // Get best-fit spectra
  for(int ifile=0;ifile<Nfiles;ifile++){
    fin[ifile] = new TFile(infilename[ifile],"READ");
    h_spectrum_bestfit[ifile] = (TH1D*)fin[ifile]->Get(HistoName[ifile]);
  }
  
  // Get Mach3 spectra
  for(int ifile=0;ifile<Nfiles_mach3;ifile++){
    fin_mach3[ifile] = new TFile(infilename_mach3[ifile],"READ");
    h_spectrum_bestfit_mach3[ifile] = (TH1D*)fin_mach3[ifile]->Get(HistoName_mach3[ifile]);
  }
  
  // Get data
  h_spectrum_data = (TH1D*)fin[0]->Get("hErecoToy");
  
  TCanvas *c = new TCanvas("c","c");
  h_spectrum_data->Draw("pe");
  for(int ifile=0;ifile<Nfiles;ifile++){
    h_spectrum_bestfit[ifile]->Draw("same");
    h_spectrum_bestfit[ifile]->SetLineColor(ifile+1);
    h_spectrum_bestfit[ifile]->SetLineWidth(2);
  }
  
  h_spectrum_bestfit_mach3[0]->Draw("histo same");
  h_spectrum_bestfit_mach3[1]->Draw("histo same");
  h_spectrum_bestfit_mach3[0]->SetLineColor(kCyan);
  h_spectrum_bestfit_mach3[1]->SetLineColor(kViolet);

  TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);
  for(int ifile=0;ifile<Nfiles;ifile++){
    leg->AddEntry(h_spectrum_bestfit[ifile],AnalysisMode[ifile],"l");
  }
  for(int ifile=0;ifile<Nfiles_mach3;ifile++){
    leg->AddEntry(h_spectrum_bestfit_mach3[ifile],AnalysisMode_mach3[ifile],"l");
  }
  leg->Draw();

  return;
}
