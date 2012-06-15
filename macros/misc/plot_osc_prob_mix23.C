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
void plot_osc_prob_mix23(
			 bool theta = true   /*  This parameter selects sinsq(theta) if true or sinsq(2theta) if false */ 
)
{
  using namespace std;
  using namespace t2kana;

  const char * hypothesis         = "3f_osc";                                       /* options: "2f_numu_disap", "3f_osc" */

  if (theta)   Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  else         Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;

  if (theta){
    double       input_mixing12     = 0.32;                                         /* value from nue appearance TN52 */
    double       input_mixing13     = 0.025;                                           /* */
    double       input_dcp          = 0.00;                                           /* */
    double       input_dm12sq       = 7.6E-5;                                         /* value from nue appearance TN52 */
    double       input_dm23sq       = 2.4E-3;   
    
    double mix23min = 0.25;
    double mix23max = 0.75;
  } 
  else {
    double       input_mixing12     = 0.8704;                                         /* value from nue appearance TN52 */
    double       input_mixing13     = 0.1;                                           /* */
    double       input_dcp          = 0.00;                                           /* */
    double       input_dm12sq       = 7.6E-5;                                         /* value from nue appearance TN52 */
    double       input_dm23sq       = 2.4E-3;   
  
    double mix23min = 0.8;
    double mix23max = 1.0;
  }

  double       input_energy       = 1.0;
  const int n = 1000;
  double dmix23   = (mix23max-mix23min)/(n-1);


  PhysHypothesis ph(hypothesis);
  ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
  ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
  ph.SetCurrValue(conv::kPhypDeltaCP,        input_dcp);
  ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
  ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  cout << "Set hypothesis" << endl;

  TH1D * numu_numu_prob           = new TH1D("numu_numu_prob", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_numu_prob_vac       = new TH1D("numu_numu_prob_vac", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_numu_prob_vac2       = new TH1D("numu_numu_prob_vac2", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_numu_prob_vac3       = new TH1D("numu_numu_prob_vac3", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_numu_prob_vac4       = new TH1D("numu_numu_prob_vac4", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_numu_prob_vac5       = new TH1D("numu_numu_prob_vac5", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * approx_numu_numu_prob    = new TH1D("approx_numu_numu_prob", " P_{#nu_{#mu}->#nu_{#mu}} ",
                                             n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_numu_term12         = new TH1D("numu_numu_term12", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_numu_term13         = new TH1D("numu_numu_term13", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_numu_term23         = new TH1D("numu_numu_term23", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);



  TH1D * numu_nue_prob            = new TH1D("numu_nue_prob", " P_{#nu_{#mu}->#nu_{e}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * numu_nue_prob_vac        = new TH1D("numu_nue_prob_vac", " P_{#nu_{#mu}->#nu_{e}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);

  TH1D * approx_numu_nue_prob     = new TH1D("approx_numu_nue_prob", " P_{#nu_{#mu}->#nu_{e}} ",
					     n, mix23min - 0.5*dmix23, mix23max + 0.5*dmix23);



  if (theta) {
    numu_numu_prob     ->GetXaxis() -> SetTitle ("sin^{2}#theta_{23}");
    numu_numu_prob_vac ->GetXaxis() -> SetTitle ("sin^{2}#theta_{23}");
    numu_nue_prob      ->GetXaxis() -> SetTitle ("sin^{2}#theta_{23}");
    numu_nue_prob_vac  ->GetXaxis() -> SetTitle ("sin^{2}#theta_{23}");
    numu_numu_term13   ->GetXaxis() -> SetTitle ("sin^{2}#theta_{23}");
  }
  else {
    numu_numu_prob     ->GetXaxis() -> SetTitle ("sin^{2}2#theta_{23}");
    numu_numu_prob_vac ->GetXaxis() -> SetTitle ("sin^{2}2#theta_{23}");
    numu_nue_prob      ->GetXaxis() -> SetTitle ("sin^{2}2#theta_{23}");
    numu_nue_prob_vac  ->GetXaxis() -> SetTitle ("sin^{2}2#theta_{23}");
    numu_numu_term13   ->GetXaxis() -> SetTitle ("sin^{2}2#theta_{23}");
  }

  NuOsc3Flavour fNuOsc3Flavour;

  for(int i=1; i<=n; i++) {
    double mix23 = mix23min + (i-1)*dmix23;
    ph.SetCurrValue(conv::kPhypMixing23, mix23);    

    if (theta) fNuOsc3Flavour.SetOscParametersSinSqTAll(input_mixing12,input_mixing13,mix23,input_dm12sq,input_dm23sq,input_dcp);    
    else fNuOsc3Flavour.SetOscParameters(input_mixing12,input_mixing13,mix23,input_dm12sq,input_dm23sq,input_dcp);

    double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(14,    14,    input_energy);
    double p_numu_numu_vac   = fNuOsc3Flavour.CalcVacuumOscProb(14,    14,    input_energy);
    numu_numu_prob    ->SetBinContent(i,p_numu_numu);
    numu_numu_prob_vac->SetBinContent(i,p_numu_numu_vac);
    //    cout << "numu_numu_prob=" << p_numu_numu << endl;

    double term12 = fNuOsc3Flavour.GetTerm12();  
    double term13 = fNuOsc3Flavour.GetTerm13();                                                                                         
    double term23 = fNuOsc3Flavour.GetTerm23();  

    numu_numu_term12 ->SetBinContent(i,1.-4*term12);
    numu_numu_term13 ->SetBinContent(i,1.-4*term13);
    numu_numu_term23 ->SetBinContent(i,1.-4*term23);

    double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(14,    12,     input_energy);    
    double p_numu_nue_vac    = fNuOsc3Flavour.CalcVacuumOscProb(14,    12,     input_energy);    
    numu_nue_prob    ->SetBinContent(i,p_numu_nue);
    numu_nue_prob_vac->SetBinContent(i,p_numu_nue_vac);
  
    double approx_numu_numu = NuMuNuMuOscProbApprox(theta, mix23, input_dm23sq, input_energy);
    approx_numu_numu_prob -> SetBinContent(i,approx_numu_numu);

    double approx_numu_nue = NuMuNuEOscProbApprox(theta, mix23, input_mixing13, input_dm23sq, input_dm12sq, input_energy);
    approx_numu_nue_prob -> SetBinContent(i,approx_numu_nue);
  }

  gStyle->SetOptStat(0);
  gStyle->SetFillColor(kWhite);

  TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  TH1F * h1 = (TH1F*) c1->DrawFrame(mix23min,0.,mix23max,1.);
  c1->SetFillColor(kWhite);
  h1->SetTitle("Oscillation probabilities");
  h1->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  h1->GetYaxis()->SetTitle("P_{#nu_{#mu}->#nu_{#mu}}");
  numu_numu_prob->Draw("c");
  numu_numu_prob_vac-> SetLineColor(kRed);
  numu_numu_prob_vac->Draw("csame");
  approx_numu_numu_prob->SetLineStyle(2);
  approx_numu_numu_prob->Draw("csame");

  TLegend * legend1 = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);
  legend1->AddEntry(numu_numu_prob, "P_{#nu_{#mu}->#nu_{#mu}} 3F matter" ,"L");
  legend1->AddEntry(numu_numu_prob_vac, "P_{#nu_{#mu}->#nu_{#mu}} 3F vacuum" ,"L");
  legend1->AddEntry(approx_numu_numu_prob, "P_{#nu_{#mu}->#nu_{#mu}} approx." ,"L");
  legend1->Draw("same");
  c1->Update();

  TCanvas * c2 = new TCanvas("c2","",20,20,500,500);
  TH1F * h2 = (TH1F*) c2->DrawFrame(mix23min,0.,mix23max,1.);
  c2->SetFillColor(kWhite);
  h2->SetTitle("Oscillation probabilities");
  h2->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  h2->GetYaxis()->SetTitle("P_{#nu_{#mu}->#nu_{e}}");
  numu_nue_prob->Draw("c");
  numu_nue_prob_vac->SetLineColor(kBlue);
  numu_nue_prob_vac->Draw("csame");
  approx_numu_nue_prob->SetLineStyle(2);
  approx_numu_nue_prob->Draw("csame");
  TLegend * legend2 = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);
  legend2->AddEntry(numu_nue_prob, "P_{#nu_{#mu}->#nu_{e}} 3F matter" ,"L");
  legend2->AddEntry(numu_nue_prob_vac, "P_{#nu_{#mu}->#nu_{e}} 3F vacuum" ,"L");
  legend2->AddEntry(approx_numu_nue_prob, "P_{#nu_{#mu}->#nu_{e}} approx." ,"L");
  legend2->Draw("same");
  c2->Update();

  TCanvas * c3 = new TCanvas("c3","",20,20,500,500);
  TH1F * h3 = (TH1F*) c3->DrawFrame(mix23min,0.,mix23max,1.);
  c3->SetFillColor(kWhite);
  h3->SetTitle("Oscillation probabilities");
  h3->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  h3->GetYaxis()->SetTitle("P_{#nu_{#mu}->#nu_{e}}");
  numu_numu_term13 -> SetLineColor(kBlue);
  numu_numu_term13 -> Draw("c");
  numu_numu_term23 -> SetLineColor(kRed);
  numu_numu_term23 -> Draw("csame");
  numu_numu_term12 -> Draw("csame");
  TLegend * legend3 = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend3->SetFillColor(0);
  legend3->SetBorderSize(0);
  legend3->AddEntry(numu_numu_term12, "1 - 4 (|U_{#mu1}|^{2}|U_{#mu2}|^{2}sin^{2}(1.27#Deltam^{2}_{12}L/E))" ,"L");
  legend3->AddEntry(numu_numu_term13, "1 - 4 (|U_{#mu1}|^{2}|U_{#mu3}|^{2}sin^{2}(1.27#Deltam^{2}_{13}L/E))" ,"L");
  legend3->AddEntry(numu_numu_term23, "1 - 4 (|U_{#mu2}|^{2}|U_{#mu3}|^{2}sin^{2}(1.27#Deltam^{2}_{23}L/E))" ,"L");
  legend3->Draw("same");
}
