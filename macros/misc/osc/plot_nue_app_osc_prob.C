//
// ------------------------------------------------------------------------------------------------
// Plot numu survival probability and nue appearance probability as function of angle theta23
// for a fixed energy and rest of parameters, for 3F in vacuum and with matter effects. 
//
// Lorena Escudero
// ------------------------------------------------------------------------------------------------
//

double NuMuNuMuOscProbApprox(double theta, double mixing23, double dm23sq, double energy){
  if(theta){
    double aux = TMath::ASin(TMath::Sqrt(mixing23));
    aux *= 2;
    double sinsq2theta23 = TMath::Power(sin(aux),2);
  }
  else{
    double sinsq2theta23 = mixing23;
  }
  double L = 295.;
  double sinsq = TMath::Power(sin((1.27*dm23sq*L)/(energy)),2);
  double prob = 1.0;
  prob -= sinsq2theta23 * sinsq;
  //  cout << "prob approx = " << prob << endl;
  return prob;
}

double NuMuNuEOscProbApprox(double theta, double mixing23, double mixing13, double dm23sq, double dm12sq, double energy){
  if(theta) {
    double sinsqtheta23 = mixing23;
    double aux = TMath::ASin(TMath::Sqrt(mixing13));
    aux *= 2;
    double sinsq2theta13 = TMath::Power(sin(aux),2);
  }
  else {
    double sinsq2theta13 = mixing13;
    double aux = TMath::ASin(TMath::Sqrt(mixing23));
    aux /= 2;
    double sinsqtheta23 = TMath::Power(sin(aux),2);
  }
  double L = 295.;
  double dm13sq = dm12sq + dm23sq;
  double sinsq = TMath::Power(sin((1.27*dm13sq*L)/(energy)),2);
  double prob = sinsq2theta13*sinsqtheta23*sinsq;
  return prob;
}
void plot_nue_app_osc_prob(
			   double       input_mixing12     = 0.306,                                        
			   double       input_mixing13     = 0.025658,
			   double       input_mixing23     = 0.5,                                        
			   double       input_dcp          = 0.0,//0.78539816,
			   double       input_dm12sq       = 7.5E-5,                                       
			   double       input_dm23sq       = 2.4E-3,   
			   bool         normal = true,
			   int          n_bins = 10000,
			   double       E_min = 0.0,
			   double       E_max = 2.0
)
{
 

  TH1D * numu_nue_probt           = new TH1D("numu_nue_probt", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);
    
  TH1D * numu_nue_prob0           = new TH1D("numu_nue_prob0", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_prob1           = new TH1D("numu_nue_prob1", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_prob2           = new TH1D("numu_nue_prob2", " P_{#nu_{#mu}->#nu_{e}} ",
					     n_bins,E_min,E_max);

  TH1D * numu_nue_prob3           = new TH1D("numu_nue_prob3", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_prob4           = new TH1D("numu_nue_prob4", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_prob5           = new TH1D("numu_nue_prob5", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_ratio_prob0           = new TH1D("numu_nue_ratio_prob0", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_ratio_prob1           = new TH1D("numu_nue_ratio_prob1", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_ratio_prob2           = new TH1D("numu_nue_ratio_prob2", " P_{#nu_{#mu}->#nu_{e}} ",
					     n_bins,E_min,E_max);

  TH1D * numu_nue_ratio_prob3           = new TH1D("numu_nue_ratio_prob3", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_ratio_prob4           = new TH1D("numu_nue_ratio_prob4", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  TH1D * numu_nue_ratio_prob5           = new TH1D("numu_nue_ratio_prob5", " P_{#nu_{#mu}->#nu_{e}} ",
					      n_bins,E_min,E_max);

  
 double E_step = (E_max-E_min)/n_bins;


  double c12 = sqrt(1-input_mixing12);
  double s12 = sqrt(input_mixing12);

  double c13 = sqrt(1-input_mixing13);
  double s13 = sqrt(input_mixing13);

  double c23 = sqrt(1-input_mixing23);
  double s23 = sqrt(input_mixing23);
  
  if(!normal) input_dm23sq = -input_dm23sq;
  double input_dm13sq = input_dm23sq + input_dm12sq;

  double density = 2.6;
  double baseline = 295;

  for(int ir = 1; ir <= n_bins; ir++) {
      
    double Energy = ir*E_step;
    double phi12 = (1.27*baseline*input_dm12sq)/(Energy);
    double phi13 = (1.27*baseline*input_dm13sq)/(Energy);
    double phi23 = (1.27*baseline*input_dm23sq)/(Energy);
    double a = 7.56E-5 *density*Energy;

    double P_leading = 4*c13*c13*s13*s13*s23*s23*sin(phi13)*sin(phi13); //- 8.*c13*c13*s12*s13*s23*(s12*s13*s23)*cos(phi23)*sin(phi13)*sin(phi12);
    double P_matter_1 = (P_leading * 2*a*(1-2*s13*s13))/input_dm13sq;
    double P_matter_2 = (-8*c13*c13*s13*s13*s23*s23*(1-2*s13*s13)*cos(phi23)*sin(phi13)*1.27*a*baseline)/(Energy);
    double P_cpc = 8.*c13*c13*s12*s13*s23*(c12*c23*cos(input_dcp)-s12*s13*s23)*cos(phi23)*sin(phi13)*sin(phi12);
      //double P_cpc = 8.*c13*c13*s12*s13*s23*(c12*c23*cos(input_dcp))*cos(phi23)*sin(phi13)*sin(phi12);
    //-8.*s12*s12*c13*c13**c12*c23*s12*s23*s13*cos(input_dcp)*sin(phi12)*sin(phi12);
    double P_cpv = -8.*c13*c13*c12*c23*s12*s13*s23*sin(input_dcp)*sin(phi23)*sin(phi13)*sin(phi12);
    double P_solar = 4.*s12*s12*c13*c13*(c12*c12*c23*c23+s12*s12*s23*s23*s13*s13-2.*c12*c23*s12*s23*s13*cos(input_dcp))*sin(phi12)*sin(phi12);
    //double P_solar = 4.*s12*s12*c13*c13*(c12*c12*c23*c23+s12*s12*s23*s23*s13*s13)*sin(phi12)*sin(phi12);

    double P_matter = P_matter_1 + P_matter_2;
    double P_total = P_leading + P_matter +P_cpc + P_cpv + P_solar;
    numu_nue_probt->Fill(Energy,P_total);

    numu_nue_prob0->Fill(Energy,P_leading);
    numu_nue_ratio_prob0->Fill(Energy,P_leading/P_total);

    numu_nue_prob1->Fill(Energy,P_matter);
    numu_nue_ratio_prob1->Fill(Energy,P_matter/P_leading);

    numu_nue_prob2->Fill(Energy,P_cpc);
    numu_nue_ratio_prob2->Fill(Energy,P_cpc/P_leading);

    numu_nue_prob3->Fill(Energy,P_cpv);
    numu_nue_ratio_prob3->Fill(Energy,P_cpv/P_leading);

    numu_nue_prob4->Fill(Energy,P_solar);
    numu_nue_ratio_prob4->Fill(Energy,P_solar/P_leading);

  }
  TCanvas * c = new TCanvas("c","",20,20,700,1000);
  c->Divide(1,2);
  c->cd(1);
  TLine *line = new TLine(E_min,0.0,E_max,0.0);
  numu_nue_probt->SetLineColor(kBlack);
  numu_nue_probt->GetYaxis()->SetTitle("P_{#nu_{#mu}->#nu_{e}}");
  numu_nue_probt->GetYaxis()->SetRangeUser(-0.06,0.06);
  numu_nue_probt->GetXaxis()->SetRangeUser(E_min,E_max);
  numu_nue_probt->Draw("c");

  numu_nue_prob0->SetLineColor(kBlue);
  numu_nue_prob0->Draw("csame");

  numu_nue_prob1->SetLineColor(kMagenta);
  numu_nue_prob1->Draw("csame");

  numu_nue_prob2->SetLineColor(kGreen);
  numu_nue_prob2->Draw("csame");

  numu_nue_prob3->SetLineColor(kRed);
  numu_nue_prob3->Draw("csame");

  numu_nue_prob4->SetLineColor(kOrange);
  numu_nue_prob4->Draw("csame");

  line->Draw("same");
  TLegend * legend = new TLegend(0.7, 0.2, 0.9, 0.47);
  legend->AddEntry(numu_nue_probt,"Total","L");
  legend->AddEntry(numu_nue_prob0,"Leading term","L");
  legend->AddEntry(numu_nue_prob1,"Matter term","L");
  legend->AddEntry(numu_nue_prob2,"CPC term","L");
  legend->AddEntry(numu_nue_prob3,"CPV term","L");
  legend->AddEntry(numu_nue_prob4,"Solar term","L");

  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  legend->Draw("same");

  
  c->cd(2);
  numu_nue_ratio_prob1->SetLineColor(kMagenta);
  numu_nue_ratio_prob1->GetYaxis()->SetTitle("Ratio to leading term");
  numu_nue_ratio_prob1->GetYaxis()->SetRangeUser(-1.,1.);
  numu_nue_ratio_prob1->Draw("c");
  
  numu_nue_ratio_prob2->SetLineColor(kGreen);
  numu_nue_ratio_prob2->Draw("csame");

  numu_nue_ratio_prob3->SetLineColor(kRed);
  numu_nue_ratio_prob3->Draw("csame");
  
  numu_nue_ratio_prob4->SetLineColor(kOrange);
  numu_nue_ratio_prob4->Draw("csame");
  
  line->Draw("same");
  legend->Draw("same");
  c->SaveAs("Numu_to_nue_osc_prob_terms.eps");

}
