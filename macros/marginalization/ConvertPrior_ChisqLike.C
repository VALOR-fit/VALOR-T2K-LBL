
// Take the chisq distribution as input and make the likelihood distribution 
// to be used as prior probability distribution for osc analyses

void ConvertPrior_ChisqLike
(
 const char * tag    = "Priors_deltaCP_run14",   
 const char * filein_prior_chisq = "../../inputs/priors/deltaCP_T2Kreact_run14/data_fit_joint_final_reactor_nh_dcp_abs23.root",
 const char * namein_prior_chisq = "h_chisq_deltacp"
 )
{
  // Get chisq prior
  TFile *fPrior_Chisqdcp = new TFile(filein_prior_chisq,"READ");
  TH1D *hPrior_Chisqdcp = (TH1D*) fPrior_Chisqdcp->Get(namein_prior_chisq);
  TH1D *hPrior_Likedcp = hPrior_Chisqdcp->Clone("hPrior_Ldcp");
  
  for(int ichisq=0;ichisq<hPrior_Chisqdcp->GetNbinsX();ichisq++){
    double chisq = hPrior_Chisqdcp->GetBinContent(ichisq+1);
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
  
  return;

}
