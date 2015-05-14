
// Take the chisq distribution as input and make the likelihood distribution 
// to be used as prior probability distribution for osc analyses

void ConvertPrior_ChisqLike
(
 const char * tag    = "Priors_deltaCP_run14_IH",   
 
 const char * filein_prior_chisq = "../../inputs/priors/joint_T2Kreact_run14/data_fit_joint_final_reactor_ih_dcp_abs31.root",
 //const char * filein_prior_chisq = "../../inputs/priors/joint_T2Kreact_run14/data_fit_joint_final_reactor_nh_dcp_abs23.root",
 
 const char * namein_prior_chisq = "h_chisq_deltacp"
 )
{
  // Get chisq prior
  TFile *fPrior_Chisqdcp = new TFile(filein_prior_chisq,"READ");
  TH1D *hPrior_Chisqdcp = (TH1D*) fPrior_Chisqdcp->Get(namein_prior_chisq);

  // Change the X axis: from -1/+1 to -3.14/+3.14
  double binwidth = 2*TMath::Pi() / (hPrior_Chisqdcp->GetNbinsX()-1);
  
  TH1D *hPrior_Chisqdcp_newaxis = new TH1D("hPrior_Chisqdcp_newaxis","hPrior_Chisqdcp_newaxis",hPrior_Chisqdcp->GetNbinsX(),-TMath::Pi()-binwidth/2.,+TMath::Pi()+binwidth/2.);

  cout << endl;
  cout << "binwidth = " << binwidth << endl;
  cout << "pi = " << TMath::Pi() << endl;
  cout << "last bin center = " << hPrior_Chisqdcp_newaxis->GetBinCenter(1) << endl;
  cout << "first bin center = " << hPrior_Chisqdcp_newaxis->GetBinCenter( hPrior_Chisqdcp_newaxis->GetNbinsX() ) << endl;
  cout << endl;

  for(int ichisq=0;ichisq<hPrior_Chisqdcp->GetNbinsX();ichisq++){
    double bincont = hPrior_Chisqdcp->GetBinContent(ichisq+1);
    hPrior_Chisqdcp_newaxis->SetBinContent(ichisq+1,bincont);
  }
  

  TH1D *hPrior_Likedcp = hPrior_Chisqdcp_newaxis->Clone("hPrior_Ldcp");
  
  for(int ichisq=0;ichisq<hPrior_Chisqdcp_newaxis->GetNbinsX();ichisq++){
    double chisq = hPrior_Chisqdcp_newaxis->GetBinContent(ichisq+1);
    double like = exp(-chisq/2.);
    hPrior_Likedcp->SetBinContent(ichisq+1,like);
  }
  //hPrior_Chisqdcp->Draw();
  hPrior_Likedcp->Draw();
 
  char name[200];
  sprintf(name,"%s.root",tag);

  TFile *fout_Likedcp = new TFile(name,"RECREATE");
  hPrior_Likedcp->Write();
  hPrior_Chisqdcp->Write();
  hPrior_Chisqdcp_newaxis->Write();
  
  return;

}
