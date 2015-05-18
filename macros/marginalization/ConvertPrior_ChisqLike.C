
//
// Take the T2K+reactor joint fit chisq distribution as input and make the likelihood distribution 
// to be used as prior probability distribution for nuebar osc analysis 
//

void ConvertPrior_ChisqLike
(
 const char * tag    = "Priors_deltaCP_run14_IH",   
 
 const char * filein_prior_chisq_nh = "../../inputs/priors/joint_T2Kreact_run14/valor/data_fit_joint_final_reactor_nh_dcp_abs23.root",
 const char * filein_prior_chisq_ih = "../../inputs/priors/joint_T2Kreact_run14/valor/data_fit_joint_final_reactor_ih_dcp_abs31.root",
 
 const char * namein_prior_chisq = "h_chisq_deltacp"
 )
{
  // Normal Hierarchy
 
  // Get chisq prior
  TFile *fPrior_Chisqdcp_nh = new TFile(filein_prior_chisq_nh,"READ");
  TH1D *hPrior_Chisqdcp_nh = (TH1D*) fPrior_Chisqdcp_nh->Get(namein_prior_chisq);

  double binwidth = 2*TMath::Pi() / (hPrior_Chisqdcp_nh->GetNbinsX()-1); // Change the X axis: from -1/+1 to -3.14/+3.14
  
  TH1D *hPrior_Chisqdcp_newaxis_nh = new TH1D("hPrior_Chisqdcp_newaxis_nh","hPrior_Chisqdcp_newaxis_nh",hPrior_Chisqdcp_nh->GetNbinsX(),-TMath::Pi()-binwidth/2.,+TMath::Pi()+binwidth/2.);

  cout << endl;
  cout << "binwidth = " << binwidth << endl;
  cout << "pi = " << TMath::Pi() << endl;
  cout << "last bin center = " << hPrior_Chisqdcp_newaxis_nh->GetBinCenter(1) << endl;
  cout << "first bin center = " << hPrior_Chisqdcp_newaxis_nh->GetBinCenter( hPrior_Chisqdcp_newaxis_nh->GetNbinsX() ) << endl;
  cout << endl;

  for(int ichisq=0;ichisq<hPrior_Chisqdcp_nh->GetNbinsX();ichisq++){
    double bincont = hPrior_Chisqdcp_nh->GetBinContent(ichisq+1);
    hPrior_Chisqdcp_newaxis_nh->SetBinContent(ichisq+1,bincont);
  }
  
  TH1D *hPrior_Likedcp_nh = hPrior_Chisqdcp_newaxis_nh->Clone("hPrior_Ldcp_nh");
  
  for(int ichisq=0;ichisq<hPrior_Chisqdcp_newaxis_nh->GetNbinsX();ichisq++){
    double chisq = hPrior_Chisqdcp_newaxis_nh->GetBinContent(ichisq+1);
    double like = exp(-chisq/2.);
    hPrior_Likedcp_nh->SetBinContent(ichisq+1,like);
  }


  // Inverted Hierarchy
 
  // Get chisq prior
  TFile *fPrior_Chisqdcp_ih = new TFile(filein_prior_chisq_ih,"READ");
  TH1D *hPrior_Chisqdcp_ih = (TH1D*) fPrior_Chisqdcp_ih->Get(namein_prior_chisq);
 
  double binwidth = 2*TMath::Pi() / (hPrior_Chisqdcp_ih->GetNbinsX()-1); // Change the X axis: from -1/+1 to -3.14/+3.14
  
  TH1D *hPrior_Chisqdcp_newaxis_ih = new TH1D("hPrior_Chisqdcp_newaxis_ih","hPrior_Chisqdcp_newaxis_ih",hPrior_Chisqdcp_ih->GetNbinsX(),-TMath::Pi()-binwidth/2.,+TMath::Pi()+binwidth/2.);

  cout << endl;
  cout << "binwidth = " << binwidth << endl;
  cout << "pi = " << TMath::Pi() << endl;
  cout << "last bin center = " << hPrior_Chisqdcp_newaxis_ih->GetBinCenter(1) << endl;
  cout << "first bin center = " << hPrior_Chisqdcp_newaxis_ih->GetBinCenter( hPrior_Chisqdcp_newaxis_ih->GetNbinsX() ) << endl;
  cout << endl;

  for(int ichisq=0;ichisq<hPrior_Chisqdcp_ih->GetNbinsX();ichisq++){
    double bincont = hPrior_Chisqdcp_ih->GetBinContent(ichisq+1);
    hPrior_Chisqdcp_newaxis_ih->SetBinContent(ichisq+1,bincont);
  }
  
  TH1D *hPrior_Likedcp_ih = hPrior_Chisqdcp_newaxis_ih->Clone("hPrior_Ldcp_ih");
  
  for(int ichisq=0;ichisq<hPrior_Chisqdcp_newaxis_ih->GetNbinsX();ichisq++){
    double chisq = hPrior_Chisqdcp_newaxis_ih->GetBinContent(ichisq+1);
    double like = exp(-chisq/2.);
    hPrior_Likedcp_ih->SetBinContent(ichisq+1,like);
  }

  
  // Get Mass Hierarchy probabilities
  
  double bin_minchisq_nh = hPrior_Chisqdcp_nh->GetMinimumBin(); 
  double bin_minchisq_ih = hPrior_Chisqdcp_ih->GetMinimumBin(); 
  //cout << bin_minchisq_nh << " " << bin_minchisq_ih << endl;

  double chisq_nh = hPrior_Chisqdcp_nh->GetBinContent(bin_minchisq_nh);
  double chisq_ih = hPrior_Chisqdcp_ih->GetBinContent(bin_minchisq_ih);
  //cout << chisq_nh << " " << chisq_ih << endl;
  
  double deltachisq_mh = 0;
  if(chisq_nh > chisq_ih) deltachisq_mh = chisq_nh - chisq_ih;
  else deltachisq_mh = chisq_ih - chisq_nh;
  
  double deltalike_mh = exp(-deltachisq_mh/2.);

  cout << endl;
  cout << "Prob (NH) = " << deltalike_mh << endl;
  cout << "Prob (IH) = " << 1 - deltalike_mh << endl;
  cout << endl;

  return;

  
  // Draw priors

  TCanvas *cPrior_Chisqdcp = new TCanvas("cPrior_Chisqdcp","cPrior_Chisqdcp");
  double max = TMath::Max(hPrior_Chisqdcp_nh->GetMaximum(),hPrior_Chisqdcp_ih->GetMaximum());
  hPrior_Chisqdcp_nh->SetMaximum(max*1.2);
  hPrior_Chisqdcp_ih->SetMaximum(max*1.2);
  hPrior_Chisqdcp_nh->SetLineColor(kBlack);
  hPrior_Chisqdcp_ih->SetLineColor(kRed);
  hPrior_Chisqdcp_nh->Draw();
  hPrior_Chisqdcp_ih->Draw("same");

  TCanvas *cPrior_Likedcp = new TCanvas("cPrior_Likedcp","cPrior_Likedcp");
  max = TMath::Max(hPrior_Likedcp_nh->GetMaximum(),hPrior_Likedcp_ih->GetMaximum());
  hPrior_Likedcp_nh->SetMaximum(max*1.2);
  hPrior_Likedcp_ih->SetMaximum(max*1.2);
  hPrior_Likedcp_nh->SetLineColor(kBlack);
  hPrior_Likedcp_ih->SetLineColor(kRed);
  hPrior_Likedcp_ih->Draw();
  hPrior_Likedcp_nh->Draw("same");


 

  char name[200];
  sprintf(name,"%s.root",tag);

  TFile *fout_Likedcp = new TFile(name,"RECREATE");
  hPrior_Likedcp_nh->Write();
  hPrior_Likedcp_ih->Write();
  hPrior_Chisqdcp_nh->Write();
  hPrior_Chisqdcp_ih->Write();
  hPrior_Chisqdcp_newaxis_nh->Write();
  hPrior_Chisqdcp_newaxis_ih->Write();
  
  return;

}
