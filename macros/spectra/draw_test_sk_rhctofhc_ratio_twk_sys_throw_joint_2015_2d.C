
void draw_test_sk_rhctofhc_ratio_twk_sys_throw_joint_2015_2d
(
 char *tag = "postfit",
 )
{

  string name_skdetfsi[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_SKDet_FSI.root",
    "SKDet+FSI+SI"
  };
  string name_flux[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_BANFF_flux.root",
    "Flux"
  };
  string name_xsec[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_BANFF_xsec.root",
    "XSec Tot"
  };

  string name_banfftot[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_BANFF.root",
    "Flux+XSec"
  };
  
  string name_skdetfsi_banfftot[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_skdetfsi_banff.root",
    "Flux+XSec+SDet+FSI+SI"
  }
  string name_osc[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_Osc_param.root",
    "Oscillations"
  }

  string name_all[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_all.root",
    "All"
  }

  string name_xsec_ncother[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_BANFF_xsec.root",
    "XSec NC other (uncorr)"
  };

  string name_xsec_nc1gamma[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_BANFF_xsec.root",
    "XSec NC $1\\gamma$ (uncorr)"
  };

  string name_xsec_nuenumu[2] = {
    "Ratio_FHC_RHC/sk_rhctofhc_ratio-twk_sys_throw_sinsqtheta23_-99.000,dcp_-1.601,dmsq_Normal23_-99.000000eV2_BANFF_xsec.root",
    "XSec $\\nu_e$ / $\\nu_{\\mu}$ (uncorr)"
  };




  // Check empty and full files
  //bool 
  


  TH1D *hskdetfsi_elike; 
  TH1D *hskdetfsi_mulike; 
  TH1D *hflux_elike;
  TH1D *hflux_mulike;
  TH1D *hxsec_ncother_elike;
  TH1D *hxsec_ncother_mulike;
  TH1D *hxsec_nc1gamma_elike;
  TH1D *hxsec_nc1gamma_mulike;
  TH1D *hxsec_nuenumu_elike;
  TH1D *hxsec_nuenumu_mulike;
  TH1D *hxsec_elike;
  TH1D *hxsec_mulike;
  TH1D *hbanfftot_elike;
  TH1D *hbanfftot_mulike;
  TH1D *hskdetfsi_banfftot_elike;
  TH1D *hskdetfsi_banfftot_mulike;
  TH1D *hosc_elike;
  TH1D *hosc_mulike;
  TH1D *hall_elike;
  TH1D *hall_mulike;

  char *title_histo_elike = "Re_ratio";
  char *title_histo_mulike = "Rmu_ratio";

  TFile *file_skdetfsi = new TFile(name_skdetfsi[0].c_str(),"READ");
  if(!file_skdetfsi->IsOpen()){ cerr << "File " << name_skdetfsi[0] << " is not open" << endl; exit(1);}
  hskdetfsi_elike = (TH1D*) file_skdetfsi->Get(title_histo_elike);
  hskdetfsi_mulike = (TH1D*) file_skdetfsi->Get(title_histo_mulike);

  TFile *file_flux = new TFile(name_flux[0].c_str(),"READ");
  if(!file_flux->IsOpen()){ cerr << "File " << name_flux[0] << " is not open" << endl; exit(1);}
  hflux_elike = (TH1D*) file_flux->Get(title_histo_elike);
  hflux_mulike = (TH1D*) file_flux->Get(title_histo_mulike);

  TFile *file_xsec_ncother = new TFile(name_xsec_ncother[0].c_str(),"READ");
  if(!file_xsec_ncother->IsOpen()){ cerr << "File " << name_xsec_ncother[0] << " is not open" << endl; exit(1);}
  hxsec_ncother_elike = (TH1D*) file_xsec_ncother->Get(title_histo_elike);
  hxsec_ncother_mulike = (TH1D*) file_xsec_ncother->Get(title_histo_mulike);

  TFile *file_xsec_nc1gamma = new TFile(name_xsec_nc1gamma[0].c_str(),"READ");
  if(!file_xsec_nc1gamma->IsOpen()){ cerr << "File " << name_xsec_nc1gamma[0] << " is not open" << endl; exit(1);}
  hxsec_nc1gamma_elike = (TH1D*) file_xsec_nc1gamma->Get(title_histo_elike);
  hxsec_nc1gamma_mulike = (TH1D*) file_xsec_nc1gamma->Get(title_histo_mulike);

  TFile *file_xsec_nuenumu = new TFile(name_xsec_nuenumu[0].c_str(),"READ");
  if(!file_xsec_nuenumu->IsOpen()){ cerr << "File " << name_xsec_nuenumu[0] << " is not open" << endl; exit(1);}
  hxsec_nuenumu_elike = (TH1D*) file_xsec_nuenumu->Get(title_histo_elike);
  hxsec_nuenumu_mulike = (TH1D*) file_xsec_nuenumu->Get(title_histo_mulike);
  
  TFile *file_xsec = new TFile(name_xsec[0].c_str(),"READ");
  if(!file_xsec->IsOpen()){ cerr << "File " << name_xsec[0] << " is not open" << endl; exit(1);}
  hxsec_elike = (TH1D*) file_xsec->Get(title_histo_elike);
  hxsec_mulike = (TH1D*) file_xsec->Get(title_histo_mulike);

  TFile *file_banfftot = new TFile(name_banfftot[0].c_str(),"READ");
  if(!file_banfftot->IsOpen()){ cerr << "File " << name_banfftot[0] << " is not open" << endl; exit(1);}
  hbanfftot_elike = (TH1D*) file_banfftot->Get(title_histo_elike);
  hbanfftot_mulike = (TH1D*) file_banfftot->Get(title_histo_mulike);

  TFile *file_skdetfsi_banfftot = new TFile(name_skdetfsi_banfftot[0].c_str(),"READ");
  if(!file_skdetfsi_banfftot->IsOpen()){ cerr << "File " << name_skdetfsi_banfftot[0] << " is not open" << endl; exit(1);}
  hskdetfsi_banfftot_elike = (TH1D*) file_skdetfsi_banfftot->Get(title_histo_elike);
  hskdetfsi_banfftot_mulike = (TH1D*) file_skdetfsi_banfftot->Get(title_histo_mulike);

  TFile *file_osc = new TFile(name_osc[0].c_str(),"READ");
  if(!file_osc->IsOpen()){ cerr << "File " << name_osc[0] << " is not open" << endl; exit(1);}
  hosc_elike = (TH1D*) file_osc->Get(title_histo_elike);
  hosc_mulike = (TH1D*) file_osc->Get(title_histo_mulike);

  TFile *file_all = new TFile(name_all[0].c_str(),"READ");
  if(!file_all->IsOpen()){ cerr << "File " << name_all[0] << " is not open" << endl; exit(1);}
  hall_elike = (TH1D*) file_all->Get(title_histo_elike);
  hall_mulike = (TH1D*) file_all->Get(title_histo_mulike);
  

  TCanvas *c_ratio_elike = new TCanvas("c_ratio_elike","c_ratio_elike");
  hall_elike->Draw("");
  hskdetfsi_elike->Draw("same");
  hflux_elike->Draw("same");
  //hxsec_ncother_elike->Draw("same");
  //hxsec_nc1gamma_elike->Draw("same");
  //hxsec_nuenumu_elike->Draw("same");
  hxsec_elike->Draw("same");
  hbanfftot_elike->Draw("same");
  hskdetfsi_banfftot_elike->Draw("same");
  //hosc_elike->Draw("same");

  hall_elike->GetXaxis()->SetTitle("Ratio of # of events FHC / RHC");
  hall_elike->GetYaxis()->SetTitle("Entries");

  hall_elike->SetLineColor(kBlack);
  hskdetfsi_elike->SetLineColor(kRed);
  hflux_elike->SetLineColor(kBlue);
  hxsec_ncother_elike->SetLineColor(kGreen);
  hxsec_nc1gamma_elike->SetLineColor(kGreen);
  hxsec_nuenumu_elike->SetLineColor(kGreen);
  hxsec_elike->SetLineColor(kGreen);
  hbanfftot_elike->SetLineColor(kCyan);
  hskdetfsi_banfftot_elike->SetLineColor(kViolet);
  hosc_elike->SetLineColor(kGreen+2);

  hall_elike->SetLineWidth(2);
  hskdetfsi_elike->SetLineWidth(2);
  hflux_elike->SetLineWidth(2);
  hxsec_ncother_elike->SetLineWidth(2);
  hxsec_nc1gamma_elike->SetLineWidth(2);
  hxsec_nuenumu_elike->SetLineWidth(2);
  hxsec_elike->SetLineWidth(2);
  hbanfftot_elike->SetLineWidth(2);
  hskdetfsi_banfftot_elike->SetLineWidth(2);
  hosc_elike->SetLineWidth(2);

  TLegend *l_ratio_elike = new TLegend(0.7,0.6,1.0,0.9);
  l_ratio_elike->SetFillColor(0);
  l_ratio_elike->SetLineColor(0);
  l_ratio_elike->SetTextFont(132);
  l_ratio_elike->SetTextSize(0.06);
  l_ratio_elike->AddEntry(hall_elike,"All","l");
  l_ratio_elike->AddEntry(hskdetfsi_elike,"SKDet+FSI+SI","l");
  l_ratio_elike->AddEntry(hflux_elike,"Flux","l");
  l_ratio_elike->AddEntry(hxsec_ncother_elike,"XSec NC other","l");
  l_ratio_elike->AddEntry(hxsec_nc1gamma_elike,"XSec NC 1#gamma","l");
  l_ratio_elike->AddEntry(hxsec_nuenumu_elike,"XSec #nu_{e}/#nu_{#mu}","l");
  l_ratio_elike->AddEntry(hxsec_elike,"XSec","l");
  l_ratio_elike->AddEntry(hbanfftot_elike,"Flux+XSec","l");
  l_ratio_elike->AddEntry(hskdetfsi_banfftot_elike,"Flux+XSec+SKDet+FSI+SI","l");
  //l_ratio_elike->AddEntry(hosc_elike,"Oscillations","l");
  l_ratio_elike->Draw();
  

  TCanvas *c_ratio_mulike = new TCanvas("c_ratio_mulike","c_ratio_mulike");
  hall_mulike->Draw("");
  hskdetfsi_mulike->Draw("same");
  hflux_mulike->Draw("same");
  hxsec_ncother_mulike->Draw("same");
  hxsec_nc1gamma_mulike->Draw("same");
  hxsec_nuenumu_mulike->Draw("same");
  hxsec_mulike->Draw("same");
  hbanfftot_mulike->Draw("same");
  hskdetfsi_banfftot_mulike->Draw("same");
  //hosc_mulike->Draw("same");

  hall_mulike->GetXaxis()->SetTitle("Ratio of # of events FHC / RHC");
  hall_mulike->GetYaxis()->SetTitle("Entries");

  //hall_mulike->SetMaximum(2100);

  hall_mulike->SetLineColor(kBlack);
  hskdetfsi_mulike->SetLineColor(kRed);
  hflux_mulike->SetLineColor(kBlue);
  hxsec_ncother_mulike->SetLineColor(kGreen);
  hxsec_nc1gamma_mulike->SetLineColor(kGreen);
  hxsec_nuenumu_mulike->SetLineColor(kGreen);
  hxsec_mulike->SetLineColor(kGreen);
  hbanfftot_mulike->SetLineColor(kCyan);
  hskdetfsi_banfftot_mulike->SetLineColor(kViolet);
  hosc_mulike->SetLineColor(kGreen+2);

  hall_mulike->SetLineWidth(2);
  hskdetfsi_mulike->SetLineWidth(2);
  hflux_mulike->SetLineWidth(2);
  hxsec_ncother_mulike->SetLineWidth(2);
  hxsec_nc1gamma_mulike->SetLineWidth(2);
  hxsec_nuenumu_mulike->SetLineWidth(2);
  hxsec_mulike->SetLineWidth(2);
  hbanfftot_mulike->SetLineWidth(2);
  hskdetfsi_banfftot_mulike->SetLineWidth(2);
  hosc_mulike->SetLineWidth(2);

  TLegend *l_ratio_mulike = new TLegend(0.7,0.6,1.0,0.9);
  l_ratio_mulike->SetFillColor(0);
  l_ratio_mulike->SetLineColor(0);
  l_ratio_mulike->SetTextFont(132);
  l_ratio_mulike->SetTextSize(0.06);
  l_ratio_mulike->AddEntry(hall_mulike,"All","l");
  l_ratio_mulike->AddEntry(hskdetfsi_mulike,"SKDet+FSI+SI","l");
  l_ratio_mulike->AddEntry(hflux_mulike,"Flux","l");
  l_ratio_mulike->AddEntry(hxsec_ncother_mulike,"XSec NC other","l");
  l_ratio_mulike->AddEntry(hxsec_nc1gamma_mulike,"XSec NC 1#gamma","l");
  l_ratio_mulike->AddEntry(hxsec_nuenumu_mulike,"XSec #nu_{e}/#nu_{#mu}","l");
  l_ratio_mulike->AddEntry(hxsec_mulike,"XSec","l");
  l_ratio_mulike->AddEntry(hbanfftot_mulike,"Flux+XSec","l");
  l_ratio_mulike->AddEntry(hskdetfsi_banfftot_mulike,"Flux+XSec+SKDet+FSI+SI","l");
  //l_ratio_mulike->AddEntry(hosc_mulike,"Oscillations","l");
  l_ratio_mulike->Draw();


  // Get variables

  double Ratio_Mean_skdetfsi_mulike = hskdetfsi_mulike->GetMean();
  double Ratio_Mean_skdetfsi_elike = hskdetfsi_elike->GetMean();
  double Ratio_RMS_skdetfsi_mulike = hskdetfsi_mulike->GetRMS();
  double Ratio_RMS_skdetfsi_elike = hskdetfsi_elike->GetRMS();
  double Ratio_Err_skdetfsi_mulike = Ratio_RMS_skdetfsi_mulike / Ratio_Mean_skdetfsi_mulike;
  double Ratio_Err_skdetfsi_elike = Ratio_RMS_skdetfsi_elike / Ratio_Mean_skdetfsi_elike;

  double Ratio_Mean_flux_mulike = hflux_mulike->GetMean();
  double Ratio_Mean_flux_elike = hflux_elike->GetMean();
  double Ratio_RMS_flux_mulike = hflux_mulike->GetRMS();
  double Ratio_RMS_flux_elike = hflux_elike->GetRMS();
  double Ratio_Err_flux_mulike = Ratio_RMS_flux_mulike / Ratio_Mean_flux_mulike;
  double Ratio_Err_flux_elike = Ratio_RMS_flux_elike / Ratio_Mean_flux_elike;

  double Ratio_Mean_xsec_ncother_mulike = hxsec_ncother_mulike->GetMean();
  double Ratio_Mean_xsec_ncother_elike = hxsec_ncother_elike->GetMean();
  double Ratio_RMS_xsec_ncother_mulike = hxsec_ncother_mulike->GetRMS();
  double Ratio_RMS_xsec_ncother_elike = hxsec_ncother_elike->GetRMS();
  double Ratio_Err_xsec_ncother_mulike = Ratio_RMS_xsec_ncother_mulike / Ratio_Mean_xsec_ncother_mulike;
  double Ratio_Err_xsec_ncother_elike = Ratio_RMS_xsec_ncother_elike / Ratio_Mean_xsec_ncother_elike;

  double Ratio_Mean_xsec_nc1gamma_mulike = hxsec_nc1gamma_mulike->GetMean();
  double Ratio_Mean_xsec_nc1gamma_elike = hxsec_nc1gamma_elike->GetMean();
  double Ratio_RMS_xsec_nc1gamma_mulike = hxsec_nc1gamma_mulike->GetRMS();
  double Ratio_RMS_xsec_nc1gamma_elike = hxsec_nc1gamma_elike->GetRMS();
  double Ratio_Err_xsec_nc1gamma_mulike = Ratio_RMS_xsec_nc1gamma_mulike / Ratio_Mean_xsec_nc1gamma_mulike;
  double Ratio_Err_xsec_nc1gamma_elike = Ratio_RMS_xsec_nc1gamma_elike / Ratio_Mean_xsec_nc1gamma_elike;

  double Ratio_Mean_xsec_nuenumu_mulike = hxsec_nuenumu_mulike->GetMean();
  double Ratio_Mean_xsec_nuenumu_elike = hxsec_nuenumu_elike->GetMean();
  double Ratio_RMS_xsec_nuenumu_mulike = hxsec_nuenumu_mulike->GetRMS();
  double Ratio_RMS_xsec_nuenumu_elike = hxsec_nuenumu_elike->GetRMS();
  double Ratio_Err_xsec_nuenumu_mulike = Ratio_RMS_xsec_nuenumu_mulike / Ratio_Mean_xsec_nuenumu_mulike;
  double Ratio_Err_xsec_nuenumu_elike = Ratio_RMS_xsec_nuenumu_elike / Ratio_Mean_xsec_nuenumu_elike;

  double Ratio_Mean_xsec_mulike = hxsec_mulike->GetMean();
  double Ratio_Mean_xsec_elike = hxsec_elike->GetMean();
  double Ratio_RMS_xsec_mulike = hxsec_mulike->GetRMS();
  double Ratio_RMS_xsec_elike = hxsec_elike->GetRMS();
  double Ratio_Err_xsec_mulike = Ratio_RMS_xsec_mulike / Ratio_Mean_xsec_mulike;
  double Ratio_Err_xsec_elike = Ratio_RMS_xsec_elike / Ratio_Mean_xsec_elike;

  double Ratio_Mean_banfftot_mulike = hbanfftot_mulike->GetMean();
  double Ratio_Mean_banfftot_elike = hbanfftot_elike->GetMean();
  double Ratio_RMS_banfftot_mulike = hbanfftot_mulike->GetRMS();
  double Ratio_RMS_banfftot_elike = hbanfftot_elike->GetRMS();
  double Ratio_Err_banfftot_mulike = Ratio_RMS_banfftot_mulike / Ratio_Mean_banfftot_mulike;
  double Ratio_Err_banfftot_elike = Ratio_RMS_banfftot_elike / Ratio_Mean_banfftot_elike;

  double Ratio_Mean_skdetfsi_banfftot_mulike = hskdetfsi_banfftot_mulike->GetMean();
  double Ratio_Mean_skdetfsi_banfftot_elike = hskdetfsi_banfftot_elike->GetMean();
  double Ratio_RMS_skdetfsi_banfftot_mulike = hskdetfsi_banfftot_mulike->GetRMS();
  double Ratio_RMS_skdetfsi_banfftot_elike = hskdetfsi_banfftot_elike->GetRMS();
  double Ratio_Err_skdetfsi_banfftot_mulike = Ratio_RMS_skdetfsi_banfftot_mulike / Ratio_Mean_skdetfsi_banfftot_mulike;
  double Ratio_Err_skdetfsi_banfftot_elike = Ratio_RMS_skdetfsi_banfftot_elike / Ratio_Mean_skdetfsi_banfftot_elike;

  double Ratio_Mean_osc_mulike = hosc_mulike->GetMean();
  double Ratio_Mean_osc_elike = hosc_elike->GetMean();
  double Ratio_RMS_osc_mulike = hosc_mulike->GetRMS();
  double Ratio_RMS_osc_elike = hosc_elike->GetRMS();
  double Ratio_Err_osc_mulike = Ratio_RMS_osc_mulike / Ratio_Mean_osc_mulike;
  double Ratio_Err_osc_elike = Ratio_RMS_osc_elike / Ratio_Mean_osc_elike;

  double Ratio_Mean_all_mulike = hall_mulike->GetMean();
  double Ratio_Mean_all_elike = hall_elike->GetMean();
  double Ratio_RMS_all_mulike = hall_mulike->GetRMS();
  double Ratio_RMS_all_elike = hall_elike->GetRMS();
  double Ratio_Err_all_mulike = Ratio_RMS_all_mulike / Ratio_Mean_all_mulike;
  double Ratio_Err_all_elike = Ratio_RMS_all_elike / Ratio_Mean_all_elike;


  // Write the table

  ofstream * outfilep = 0;
  char name[200];
  sprintf(name,"sk_rhctofhc_ratio-twk_sys_throw_%s.txt",tag);
  outfilep = new ofstream(name,ios_base::trunc);
  ofstream & outfile = *outfilep;
  
  stringstream streamer;
  
  streamer << std::fixed;
  streamer << std::setprecision(3);
  
  streamer << endl;
  streamer << "\\begin{table}"   << endl;
  streamer << "\\begin{tabular}{| c || c | c || c | c |}" << endl;
  streamer << "\\hline" << endl;
  streamer << "Source of uncertainty & "
	   << "$\\mu$-like Mean & "
	   << "$\\mu$-like RMS/Mean &"
	   << "e-like Mean &"
	   << "e-like RMS/Mean "
	   << "\\\\"
	   << endl;

  streamer << "\\hline" << endl;
   
  streamer << name_skdetfsi[1] 
	   << " & " << Ratio_Mean_skdetfsi_mulike
	   << " & " << Ratio_Err_skdetfsi_mulike
	   << " & " << Ratio_Mean_skdetfsi_elike
	   << " & " << Ratio_Err_skdetfsi_elike
	   << " \\\\ \\hline"
	   << endl;
  
  streamer << "\\hline" << endl;

  streamer << name_flux[1] 
	   << " & " << Ratio_Mean_flux_mulike
	   << " & " << Ratio_Err_flux_mulike
	   << " & " << Ratio_Mean_flux_elike
	   << " & " << Ratio_Err_flux_elike
	   << " \\\\ \\hline"
	   << endl;

  streamer << name_xsec_ncother[1] 
	   << " & " << Ratio_Mean_xsec_ncother_mulike
	   << " & " << Ratio_Err_xsec_ncother_mulike
	   << " & " << Ratio_Mean_xsec_ncother_elike
	   << " & " << Ratio_Err_xsec_ncother_elike
	   << " \\\\ \\hline"
	   << endl;

  streamer << name_xsec_nc1gamma[1] 
	   << " & " << Ratio_Mean_xsec_nc1gamma_mulike
	   << " & " << Ratio_Err_xsec_nc1gamma_mulike
	   << " & " << Ratio_Mean_xsec_nc1gamma_elike
	   << " & " << Ratio_Err_xsec_nc1gamma_elike
	   << " \\\\ \\hline"
	   << endl;

  streamer << name_xsec_nuenumu[1] 
	   << " & " << Ratio_Mean_xsec_nuenumu_mulike
	   << " & " << Ratio_Err_xsec_nuenumu_mulike
	   << " & " << Ratio_Mean_xsec_nuenumu_elike
	   << " & " << Ratio_Err_xsec_nuenumu_elike
	   << " \\\\ \\hline"
	   << endl;

  streamer << name_xsec[1] 
	   << " & " << Ratio_Mean_xsec_mulike
	   << " & " << Ratio_Err_xsec_mulike
	   << " & " << Ratio_Mean_xsec_elike
	   << " & " << Ratio_Err_xsec_elike
	   << " \\\\ \\hline"
	   << endl;

  streamer << name_banfftot[1] 
	   << " & " << Ratio_Mean_banfftot_mulike
	   << " & " << Ratio_Err_banfftot_mulike
	   << " & " << Ratio_Mean_banfftot_elike
	   << " & " << Ratio_Err_banfftot_elike
	   << " \\\\ \\hline"
	   << endl;

  streamer << " \\hline" << endl;

  streamer << name_skdetfsi_banfftot[1] 
	   << " & " << Ratio_Mean_skdetfsi_banfftot_mulike
	   << " & " << Ratio_Err_skdetfsi_banfftot_mulike
	   << " & " << Ratio_Mean_skdetfsi_banfftot_elike
	   << " & " << Ratio_Err_skdetfsi_banfftot_elike
	   << " \\\\ \\hline"
	   << endl;

  streamer << " \\hline" << endl;

  streamer << name_osc[1] 
	   << " & " << Ratio_Mean_osc_mulike
	   << " & " << Ratio_Err_osc_mulike
	   << " & " << Ratio_Mean_osc_elike
	   << " & " << Ratio_Err_osc_elike
	   << " \\\\ \\hline"
	   << endl;

  streamer << " \\hline" << endl;
  
  streamer << name_all[1] 
	   << " & " << Ratio_Mean_all_mulike
	   << " & " << Ratio_Err_all_mulike
	   << " & " << Ratio_Mean_all_elike
	   << " & " << Ratio_Err_all_elike
	   << " \\\\ \\hline"
	   << endl;
  
  streamer << "\\end{tabular}" << endl;
  streamer << "\\end{table}" << endl;

  cout << streamer.str() << endl;


  // get size of file
  long size;
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);
  // write to outfile
  outfile.write (streamer.str().c_str(),size);
  outfile.close();

  
}
