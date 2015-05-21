
//
// Take the T2K+reactor joint fit chisq distribution as input and make the likelihood distribution 
// to be used as prior probability distribution for nuebar osc analysis 
//

void ConvertPrior_ChisqLike_Valor
(
 const char * tag    = "Priors_T2Kreact_run14_Valor",   
 
 const char * filein_prior_chisq_nh    = "../../inputs/priors/joint_T2Kreact_run14/valor/data_final_reactor_abs23_nh_contours_common_bf.root",
 const char * filein_prior_chisq_ih    = "../../inputs/priors/joint_T2Kreact_run14/valor/data_final_reactor_abs31_ih_contours_common_bf.root",
 const char * namein_prior_chisq_dcp   = "h_deltacp",
 const char * namein_prior_chisq_23    = "h_theta23_dmsq"
 )
{
  TH1D *hPrior_Likedcp_nh;
  TH1D *hPrior_Chisqdcp_nh;
  TH1D *hPrior_Chisqdcp_newaxis_nh;

  TH1D *hPrior_Likedcp_ih;
  TH1D *hPrior_Chisqdcp_ih;
  TH1D *hPrior_Chisqdcp_newaxis_ih;

  TH2D *hPrior_Like23_nh;
  TH2D *hPrior_Chisq23_nh;
  TH2D *hPrior_Chisq23_newaxis_nh;

  TH2D *hPrior_Like23_ih;
  TH2D *hPrior_Chisq23_ih;
  TH2D *hPrior_Chisq23_newaxis_ih;

  TH1D *hPrior_Like_mh;
    
  // Normal Hierarchy
 
  // Get chisq prior
  TFile *fPrior_Chisqdcp_nh = new TFile(filein_prior_chisq_nh,"READ");
  hPrior_Chisqdcp_nh = (TH1D*) fPrior_Chisqdcp_nh->Get(namein_prior_chisq_dcp);
  hPrior_Chisq23_nh = (TH2D*) fPrior_Chisqdcp_nh->Get(namein_prior_chisq_23);

  // deltacp
  double binwidth = 2*TMath::Pi() / (hPrior_Chisqdcp_nh->GetNbinsX()-1); // Change the X axis: from -1/+1 to -3.14/+3.14
  hPrior_Chisqdcp_newaxis_nh = new TH1D("hPrior_Chisqdcp_newaxis_nh","hPrior_Chisqdcp_newaxis_nh",hPrior_Chisqdcp_nh->GetNbinsX(),-TMath::Pi()-binwidth/2.,+TMath::Pi()+binwidth/2.);
  // cout << endl;
  // cout << "binwidth = " << binwidth << endl;
  // cout << "pi = " << TMath::Pi() << endl;
  // cout << "last bin center = " << hPrior_Chisqdcp_newaxis_nh->GetBinCenter(1) << endl;
  // cout << "first bin center = " << hPrior_Chisqdcp_newaxis_nh->GetBinCenter( hPrior_Chisqdcp_newaxis_nh->GetNbinsX() ) << endl;
  // cout << endl;

  for(int ichisq=0;ichisq<hPrior_Chisqdcp_nh->GetNbinsX();ichisq++){
    double bincont = hPrior_Chisqdcp_nh->GetBinContent(ichisq+1);
    hPrior_Chisqdcp_newaxis_nh->SetBinContent(ichisq+1,bincont);
  } 
  hPrior_Likedcp_nh = (TH1D*) hPrior_Chisqdcp_newaxis_nh->Clone("hPrior_Ldcp_nh");
  hPrior_Likedcp_nh->SetTitle("Likelihood prior for dcp NH");

  for(int ichisq=0;ichisq<hPrior_Chisqdcp_newaxis_nh->GetNbinsX();ichisq++){
    double chisq = hPrior_Chisqdcp_newaxis_nh->GetBinContent(ichisq+1);
    double like = exp(-chisq/2.);
    hPrior_Likedcp_nh->SetBinContent(ichisq+1,like);
  }

  // 23
  TAxis *axis23_nh_x = (TAxis*) hPrior_Chisq23_nh->GetXaxis();
  TAxis *axis23_nh_y = (TAxis*) hPrior_Chisq23_nh->GetYaxis();

  int nbins_x = hPrior_Chisq23_nh->GetNbinsX();
  double lowedge_x = axis23_nh_x->GetBinLowEdge(1);
  double upedge_x = axis23_nh_x->GetBinUpEdge(nbins_x);

  int nbins_y = hPrior_Chisq23_nh->GetNbinsY();
  double lowedge_y = axis23_nh_y->GetBinLowEdge(1) * 1.E-3;
  double upedge_y = axis23_nh_y->GetBinUpEdge(nbins_y) * 1.E-3;

  hPrior_Chisq23_newaxis_nh = new TH2D("hPrior_Chisq23_newaxis_nh","hPrior_Chisq23_newaxis_nh",nbins_x,lowedge_x,upedge_x,nbins_y,lowedge_y,upedge_y);
  // cout << endl;
  // cout << "# bins x = " << nbins_x << endl;
  // cout << "Low edge x = " << lowedge_x << endl;
  // cout << "Up edge x = " << upedge_x << endl;
  // cout << "# bins y = " << nbins_y << endl;
  // cout << "Low edge y = " << lowedge_y << endl;
  // cout << "Up edge y = " << upedge_y << endl;
  // cout << endl;

  for(int ibinX=0;ibinX<nbins_x;ibinX++){
    for(int ibinY=0;ibinY<nbins_y;ibinY++){
      double bincont = hPrior_Chisq23_nh->GetBinContent(ibinX+1,ibinY+1);
      double mincont = hPrior_Chisq23_nh->GetMinimum();
      bincont -= mincont;
      
      hPrior_Chisq23_newaxis_nh->SetBinContent(ibinX+1,ibinY+1,bincont);

      // if(bincont<0){
      // 	cout << bincont << " ";
      // 	cout << "Bin: "
      // 	     << hPrior_Chisq23_nh->GetXaxis()->GetBinCenter(ibinX+1) << " "
      // 	     << hPrior_Chisq23_nh->GetYaxis()->GetBinCenter(ibinY+1) << " "
      // 	     << endl;
      // }
      
    } 
  }
  hPrior_Like23_nh = (TH2D*) hPrior_Chisq23_newaxis_nh->Clone("hPrior_L23_nh");
  hPrior_Like23_nh->SetTitle("Likelihood prior for 23 NH");

  for(int ichisqX=0;ichisqX<hPrior_Chisq23_newaxis_nh->GetNbinsX();ichisqX++){
    for(int ichisqY=0;ichisqY<hPrior_Chisq23_newaxis_nh->GetNbinsY();ichisqY++){
      double chisq = hPrior_Chisq23_newaxis_nh->GetBinContent(ichisqX+1,ichisqY+1);
      double like = exp(-chisq/2.);
      //cout << ichisqX << " - " << ichisqY << " --> chisq = " << chisq << ", like = " << like << endl;
      hPrior_Like23_nh->SetBinContent(ichisqX+1,ichisqY+1,like);
    }
  }

  // Inverted Hierarchy
 
  // Get chisq prior
  TFile *fPrior_Chisqdcp_ih = new TFile(filein_prior_chisq_ih,"READ");
  hPrior_Chisqdcp_ih = (TH1D*) fPrior_Chisqdcp_ih->Get(namein_prior_chisq_dcp);
  hPrior_Chisq23_ih = (TH2D*) fPrior_Chisqdcp_ih->Get(namein_prior_chisq_23);

  double binwidth = 2*TMath::Pi() / (hPrior_Chisqdcp_ih->GetNbinsX()-1); // Change the X axis: from -1/+1 to -3.14/+3.14
  
  hPrior_Chisqdcp_newaxis_ih = new TH1D("hPrior_Chisqdcp_newaxis_ih","hPrior_Chisqdcp_newaxis_ih",hPrior_Chisqdcp_ih->GetNbinsX(),-TMath::Pi()-binwidth/2.,+TMath::Pi()+binwidth/2.);
  // cout << endl;
  // cout << "binwidth = " << binwidth << endl;
  // cout << "pi = " << TMath::Pi() << endl;
  // cout << "last bin center = " << hPrior_Chisqdcp_newaxis_ih->GetBinCenter(1) << endl;
  // cout << "first bin center = " << hPrior_Chisqdcp_newaxis_ih->GetBinCenter( hPrior_Chisqdcp_newaxis_ih->GetNbinsX() ) << endl;
  // cout << endl;

  for(int ichisq=0;ichisq<hPrior_Chisqdcp_ih->GetNbinsX();ichisq++){
    double bincont = hPrior_Chisqdcp_ih->GetBinContent(ichisq+1);
    hPrior_Chisqdcp_newaxis_ih->SetBinContent(ichisq+1,bincont);
  }
  
  hPrior_Likedcp_ih = (TH1D*) hPrior_Chisqdcp_newaxis_ih->Clone("hPrior_Ldcp_ih");
  hPrior_Likedcp_ih->SetTitle("Likelihood prior for dcp IH");

  for(int ichisq=0;ichisq<hPrior_Chisqdcp_newaxis_ih->GetNbinsX();ichisq++){
    double chisq = hPrior_Chisqdcp_newaxis_ih->GetBinContent(ichisq+1);
    double like = exp(-chisq/2.);
    hPrior_Likedcp_ih->SetBinContent(ichisq+1,like);
  }


  // 23
  TAxis *axis23_ih_x = (TAxis*) hPrior_Chisq23_ih->GetXaxis();
  TAxis *axis23_ih_y = (TAxis*) hPrior_Chisq23_ih->GetYaxis();

  int nbins_x = hPrior_Chisq23_ih->GetNbinsX();
  double lowedge_x = axis23_ih_x->GetBinLowEdge(1);
  double upedge_x = axis23_ih_x->GetBinUpEdge(nbins_x);

  int nbins_y = hPrior_Chisq23_ih->GetNbinsY();
  double lowedge_y = axis23_ih_y->GetBinLowEdge(1) * 1.E-3;
  double upedge_y = axis23_ih_y->GetBinUpEdge(nbins_y) * 1.E-3;

  hPrior_Chisq23_newaxis_ih = new TH2D("hPrior_Chisq23_newaxis_ih","hPrior_Chisq23_newaxis_ih",nbins_x,lowedge_x,upedge_x,nbins_y,lowedge_y,upedge_y);
  // cout << endl;
  // cout << "# bins x = " << nbins_x << endl;
  // cout << "Low edge x = " << lowedge_x << endl;
  // cout << "Up edge x = " << upedge_x << endl;
  // cout << "# bins y = " << nbins_y << endl;
  // cout << "Low edge y = " << lowedge_y << endl;
  // cout << "Up edge y = " << upedge_y << endl;
  // cout << endl;

  for(int ibinX=0;ibinX<nbins_x;ibinX++){
    for(int ibinY=0;ibinY<nbins_y;ibinY++){
      double bincont = hPrior_Chisq23_ih->GetBinContent(ibinX+1,ibinY+1);
      double mincont = hPrior_Chisq23_ih->GetMinimum();
      bincont -= mincont;

      hPrior_Chisq23_newaxis_ih->SetBinContent(ibinX+1,ibinY+1,bincont);

      // if(bincont<0){
      // 	cout << bincont << " ";
      // 	cout << "Bin: "
      // 	     << hPrior_Chisq23_ih->GetXaxis()->GetBinCenter(ibinX+1) << " "
      // 	     << hPrior_Chisq23_ih->GetYaxis()->GetBinCenter(ibinY+1) << " "
      // 	     << endl;
      // }
    } 
  }
  hPrior_Like23_ih = (TH2D*) hPrior_Chisq23_newaxis_ih->Clone("hPrior_L23_ih");
  hPrior_Like23_ih->SetTitle("Likelihood prior for 23 IH");

  for(int ichisqX=0;ichisqX<hPrior_Chisq23_newaxis_ih->GetNbinsX();ichisqX++){
    for(int ichisqY=0;ichisqY<hPrior_Chisq23_newaxis_ih->GetNbinsY();ichisqY++){
      double chisq = hPrior_Chisq23_newaxis_ih->GetBinContent(ichisqX+1,ichisqY+1);
      double like = exp(-chisq/2.);
      hPrior_Like23_ih->SetBinContent(ichisqX+1,ichisqY+1,like);
    }
  }

  
  // Get Mass Hierarchy probabilities from deltaCP distribution
  
  hPrior_Like_mh = new TH1D("hPrior_Like_mh","Likelihood prior MH: bin 1 --> NH, bin 2 --> IH",2,0,2);

  int bin_minchisq_nh = hPrior_Chisqdcp_nh->GetMinimumBin(); 
  int bin_minchisq_ih = hPrior_Chisqdcp_ih->GetMinimumBin(); 
  //cout << bin_minchisq_nh << " " << bin_minchisq_ih << endl;

  double chisq_nh = hPrior_Chisqdcp_nh->GetBinContent(bin_minchisq_nh);
  double chisq_ih = hPrior_Chisqdcp_ih->GetBinContent(bin_minchisq_ih);
  //cout << chisq_nh << " " << chisq_ih << endl;
  
  double deltachisq_mh = 0;
  if(chisq_nh > chisq_ih) deltachisq_mh = chisq_nh - chisq_ih;
  else deltachisq_mh = chisq_ih - chisq_nh;
  
  double deltalike_mh = exp(-deltachisq_mh/2.);
  double likeprior_nh = deltalike_mh;
  double likeprior_ih = 1 - deltalike_mh;

  hPrior_Like_mh->SetBinContent(1,likeprior_nh);
  hPrior_Like_mh->SetBinContent(2,likeprior_ih);

  cout << endl;
  cout << "Get Mass Hierarchy probabilities from deltaCP distribution:" << endl;
  cout << "DeltaChisq (MH) = " << deltachisq_mh << endl; 
  cout << "DeltaLike (MH) = " << deltalike_mh << endl; 
  cout << "Prob (NH) = " << likeprior_nh << endl;
  cout << "Prob (IH) = " << likeprior_ih << endl;
  cout << endl;

  // Get Mass Hierarchy probabilities from 23 distribution
  
  bin_minchisq_nh = hPrior_Chisq23_nh->GetMinimumBin(); 
  bin_minchisq_ih = hPrior_Chisq23_ih->GetMinimumBin(); 
  //cout << bin_minchisq_nh << " " << bin_minchisq_ih << endl;

  chisq_nh = hPrior_Chisq23_nh->GetBinContent(bin_minchisq_nh);
  chisq_ih = hPrior_Chisq23_ih->GetBinContent(bin_minchisq_ih);
  //cout << chisq_nh << " " << chisq_ih << endl;
  
  deltachisq_mh = 0;
  if(chisq_nh > chisq_ih) deltachisq_mh = chisq_nh - chisq_ih;
  else deltachisq_mh = chisq_ih - chisq_nh;
  
  deltalike_mh = exp(-deltachisq_mh/2.);

  cout << endl;
  cout << "Get Mass Hierarchy probabilities from 23 distribution:" << endl;
  cout << "DeltaChisq (MH) = " << deltachisq_mh << endl; 
  cout << "DeltaLike (MH) = " << deltalike_mh << endl; 
  cout << "Prob (NH) = " << deltalike_mh << endl;
  cout << "Prob (IH) = " << 1 - deltalike_mh << endl;
  cout << endl;
  

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

  TCanvas *cPrior_Chisq23 = new TCanvas("cPrior_Chisq23","cPrior_Chisq23");
  hPrior_Chisq23_nh->SetLineColor(kBlack);
  hPrior_Chisq23_ih->SetLineColor(kRed);
  hPrior_Chisq23_nh->Draw("lego");
  hPrior_Chisq23_ih->Draw("lego same");

  TCanvas *cPrior_Like23 = new TCanvas("cPrior_Like23","cPrior_Like23");
  hPrior_Like23_nh->SetLineColor(kBlack);
  hPrior_Like23_ih->SetLineColor(kRed);
  hPrior_Like23_nh->Draw("lego");
  hPrior_Like23_ih->Draw("lego same");

  TCanvas *cPrior_Like_mh = new TCanvas("cPrior_Like_mh","cPrior_Like_mh");
  hPrior_Like_mh->Draw();

  // Store histogram in output root file

  char name[200];
  sprintf(name,"%s.root",tag);

  TFile *fout_Likedcp = new TFile(name,"RECREATE");

  hPrior_Likedcp_nh->Write();
  hPrior_Likedcp_ih->Write();
  hPrior_Chisqdcp_nh->Write();
  hPrior_Chisqdcp_ih->Write();
  hPrior_Chisqdcp_newaxis_nh->Write();
  hPrior_Chisqdcp_newaxis_ih->Write();

  hPrior_Like23_nh->Write();
  hPrior_Like23_ih->Write();
  hPrior_Chisq23_nh->Write();
  hPrior_Chisq23_ih->Write();
  hPrior_Chisq23_newaxis_nh->Write();
  hPrior_Chisq23_newaxis_ih->Write();
  
  hPrior_Like_mh->Write();
  
  return;
}
