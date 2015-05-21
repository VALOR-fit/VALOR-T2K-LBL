
//
// Take the T2K+reactor joint fit chisq distribution as input and make the likelihood distribution 
// to be used as prior probability distribution for nuebar osc analysis 
//

using namespace std;
using namespace valor;
using namespace valor::conv;

void ConvertPrior_ChisqLike_Mach3
(
 const char * tag    = "Priors_T2Kreact_run14_Mach3",   
 
 const char * filein_prior_chisq_dcp = "../../inputs/priors/joint_T2Kreact_run14/mach3/deltacp_T2Kreact_run14_Mach3.root",
 const char * filein_prior_chisq_23 = "../../inputs/priors/joint_T2Kreact_run14/mach3/parameters23_T2Kreact_run14_Mach3.root",

 const double prior_nh = 0.684
 )
{

  // deltacp prior (Normal Hierarchy)
  
  TH1D *hPrior_dcp_nh;
  TH1D *hPrior_dcp_ih;
  
  TFile *fPrior_dcp_nh = new TFile(filein_prior_chisq_dcp,"READ");
  hPrior_dcp_nh = (TH1D*) fPrior_dcp_nh->Get("dcp_hist");
  hPrior_dcp_nh->SetName("hPrior_Ldcp_nh");
  hPrior_dcp_nh->SetTitle("hPrior_Ldcp_nh");


  // deltacp prior (Inverted Hierarchy)

  TFile *fPrior_dcp_ih = new TFile(filein_prior_chisq_dcp,"READ");
  hPrior_dcp_ih = (TH1D*) fPrior_dcp_ih->Get("dcp_hist");
  hPrior_dcp_ih->SetName("hPrior_Ldcp_ih");
  hPrior_dcp_ih->SetTitle("hPrior_Ldcp_ih");

  
  // 23 prior (Normal Hierarchy)
  
  TH2D *hPrior_23_nh;

  TFile *fPrior_23_nh = new TFile(filein_prior_chisq_23,"READ");
  hPrior_23_nh = (TH2D*) fPrior_23_nh->Get("par23_hist_normal");
  hPrior_23_nh->SetName("hPrior_L23_nh");
  hPrior_23_nh->SetTitle("hPrior_L23_nh");

  
  // 23 prior (Inverted Hierarchy)

  TH2D *hPrior_23_ih;

  TFile *fPrior_23_ih = new TFile(filein_prior_chisq_23,"READ");
  TH2D *hPrior_23_ih_old = (TH2D*) fPrior_23_ih->Get("par23_hist_inverted");
  hPrior_23_ih_old->SetName("hPrior_L23_ih");
  hPrior_23_ih_old->SetTitle("hPrior_L23_ih");

  TAxis *axis_mix23 = (TAxis*) hPrior_23_ih_old->GetXaxis();
  TAxis *axis_dmsq  = (TAxis*) hPrior_23_ih_old->GetYaxis();
  
  TAxis *axis_mix23 = (TAxis*) hPrior_23_ih_old->GetXaxis();
  int nbins_mix23 = hPrior_23_ih_old->GetNbinsX();
  double lowedge_mix23 = axis_mix23->GetBinLowEdge(1);
  double upedge_mix23 = axis_mix23->GetBinUpEdge(nbins_mix23);
  
  // dmsq for inverted hierarchy gets negative values and uses dmsq23 --> to convert!!!
  int nbins_dmsq = hPrior_23_ih_old->GetNbinsY();
  double upedge_dmsq = fabs( axis_dmsq->GetBinLowEdge(1) );
  double lowedge_dmsq = fabs( axis_dmsq->GetBinUpEdge(nbins_dmsq) );

  // Mach3 uses dmsq23 for inverted hierarchy as well --> dmsq31 = dmsq23+dmsq12
  upedge_dmsq += conv::kT2KReacBestFitDmsq12;   
  lowedge_dmsq += conv::kT2KReacBestFitDmsq12;

  hPrior_23_ih = new TH2D("hPrior_L23_ih","hPrior_L23_ih",nbins_mix23,lowedge_mix23,upedge_mix23,nbins_dmsq,lowedge_dmsq,upedge_dmsq);

  // fill the histogram inverting the order for dmsq
  
  for(int ibinX = 0; ibinX<hPrior_23_ih_old->GetNbinsX(); ibinX++){
    for(int ibinY = 0; ibinY<hPrior_23_ih_old->GetNbinsY(); ibinY++){
      
      double bincont = hPrior_23_ih_old->GetBinContent(ibinX+1,ibinY+1);
      
      int newbinX = ibinX+1; // mixing23
      int newbinY = hPrior_23_ih_old->GetNbinsY() - ibinY; // dmsq
      //cout << ibinY+1 << " - " << newbinY << endl;

      hPrior_23_ih->SetBinContent(newbinX,newbinY,bincont);
    }
  }
  
  // Prior for Mass Hierarchy

  TH1D *hPrior_Like_mh = new TH1D("hPrior_Like_mh","hPrior_Like_mh",2,0,2);
  hPrior_Like_mh->SetBinContent(1,prior_nh);
  hPrior_Like_mh->SetBinContent(2,1-prior_nh);
  
  cout << endl;
  cout << "Prob (NH) = " << prior_nh << endl;
  cout << "Prob (IH) = " << 1 - prior_nh << endl;
  cout << endl;


  // Draw priors

  TCanvas *cPrior_dcp = new TCanvas("cPrior_dcp","cPrior_dcp");
  double max = TMath::Max(hPrior_dcp_nh->GetMaximum(),hPrior_dcp_ih->GetMaximum());
  hPrior_dcp_nh->SetMaximum(max*1.2);
  hPrior_dcp_ih->SetMaximum(max*1.2);
  hPrior_dcp_nh->SetLineColor(kBlack);
  hPrior_dcp_ih->SetLineColor(kRed);
  hPrior_dcp_nh->Draw();
  hPrior_dcp_ih->Draw("same");

  TCanvas *cPrior_23_nh = new TCanvas("cPrior_23_nh","cPrior_23_nh");
  hPrior_23_nh->Draw("colz");

  TCanvas *cPrior_23_ih_old = new TCanvas("cPrior_23_ih_old","cPrior_23_ih_old");
  hPrior_23_ih_old->Draw("colz");

  TCanvas *cPrior_23_ih = new TCanvas("cPrior_23_ih","cPrior_23_ih");
  hPrior_23_ih->Draw("colz");

  TCanvas *cPrior_Like_mh = new TCanvas("cPrior_Like_mh","cPrior_Like_mh");
  hPrior_Like_mh->Draw();

  // Store histogram in output root file

  char name[200];
  sprintf(name,"%s.root",tag);

  TFile *fout_Likedcp = new TFile(name,"RECREATE");
  hPrior_dcp_nh->Write();
  hPrior_dcp_ih->Write();
  hPrior_23_nh->Write();
  hPrior_23_ih->Write();
  hPrior_Like_mh->Write();
  
  return;
}
