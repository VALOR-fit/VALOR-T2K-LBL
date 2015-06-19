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
void plot_osc_prob_deltacp(
			   double       input_mixing12     = 0.306,                                        
			   double       input_mixing13     = 0.0251,                                        
			   double       input_mixing23     = 0.5,                                        
			   double       input_dcp          = 0.0,                                          
			   double       input_dm12sq       = 7.5E-5,                                       
			   double       input_dm23sq       = 2.4E-3,   
			   bool normal = true  
)
{
  using namespace std;
  using namespace valor;

  TStopwatch w;
  w.Start();

  const char * hypothesis         = "3f_osc";                                      

  TGaxis::SetMaxDigits(3);

  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  //Opt::Instance()->gOscDmsqScheme        = conv::kNuDmsqFogliLisi;
  if (normal) Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;   
  else      Opt::Instance()->gMassHierarchyAssumed = conv::kNmhInverted;


  PhysHypothesis ph(hypothesis);
  ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
  ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
  ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
  ph.SetCurrValue(conv::kPhypCPV    ,        input_dcp);
  ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
  ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);

 
  TH1D * numu_numu_prob1           = new TH1D("numu_numu_prob1", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob2           = new TH1D("numu_numu_prob2", " P_{#nu_{#mu}->#nu_{#mu}} ",
					     util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob3           = new TH1D("numu_numu_prob3", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob4           = new TH1D("numu_numu_prob4", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob5           = new TH1D("numu_numu_prob5", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob_vac1       = new TH1D("numu_numu_prob_vac1", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob_vac2       = new TH1D("numu_numu_prob_vac2", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob_vac3       = new TH1D("numu_numu_prob_vac3", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob_vac4       = new TH1D("numu_numu_prob_vac4", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_numu_prob_vac5       = new TH1D("numu_numu_prob_vac5", " P_{#nu_{#mu}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));


  TH1D * numubar_numubar_prob1           = new TH1D("numubar_numubar_prob1", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob2           = new TH1D("numubar_numubar_prob2", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob3           = new TH1D("numubar_numubar_prob3", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob4           = new TH1D("numubar_numubar_prob4", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob5           = new TH1D("numubar_numubar_prob5", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob_vac1       = new TH1D("numubar_numubar_prob_vac1", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob_vac2       = new TH1D("numubar_numubar_prob_vac2", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob_vac3       = new TH1D("numubar_numubar_prob_vac3", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob_vac4       = new TH1D("numubar_numubar_prob_vac4", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_numubar_prob_vac5       = new TH1D("numubar_numubar_prob_vac5", " P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));


  TH1D * nue_nue_prob1           = new TH1D("nue_nue_prob1", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob2           = new TH1D("nue_nue_prob2", " P_{#nu_{e}->#nu_{e}} ",
					     util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob3           = new TH1D("nue_nue_prob3", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob4           = new TH1D("nue_nue_prob4", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob5           = new TH1D("nue_nue_prob5", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob_vac1       = new TH1D("nue_nue_prob_vac1", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob_vac2       = new TH1D("nue_nue_prob_vac2", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob_vac3       = new TH1D("nue_nue_prob_vac3", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob_vac4       = new TH1D("nue_nue_prob_vac4", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_nue_prob_vac5       = new TH1D("nue_nue_prob_vac5", " P_{#nu_{e}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));


  TH1D * nuebar_nuebar_prob1           = new TH1D("nuebar_nuebar_prob1", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob2           = new TH1D("nuebar_nuebar_prob2", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob3           = new TH1D("nuebar_nuebar_prob3", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob4           = new TH1D("nuebar_nuebar_prob4", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob5           = new TH1D("nuebar_nuebar_prob5", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob_vac1       = new TH1D("nuebar_nuebar_prob_vac1", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob_vac2       = new TH1D("nuebar_nuebar_prob_vac2", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob_vac3       = new TH1D("nuebar_nuebar_prob_vac3", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob_vac4       = new TH1D("nuebar_nuebar_prob_vac4", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_nuebar_prob_vac5       = new TH1D("nuebar_nuebar_prob_vac5", " P_{#bar{#nu_{e}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));



  TH1D * numu_nue_prob1           = new TH1D("numu_nue_prob1", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob2           = new TH1D("numu_nue_prob2", " P_{#nu_{#mu}->#nu_{e}} ",
					     util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob3           = new TH1D("numu_nue_prob3", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob4           = new TH1D("numu_nue_prob4", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob5           = new TH1D("numu_nue_prob5", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob_vac1       = new TH1D("numu_nue_prob_vac1", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob_vac2       = new TH1D("numu_nue_prob_vac2", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob_vac3       = new TH1D("numu_nue_prob_vac3", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob_vac4       = new TH1D("numu_nue_prob_vac4", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numu_nue_prob_vac5       = new TH1D("numu_nue_prob_vac5", " P_{#nu_{#mu}->#nu_{e}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));


  TH1D * numubar_nuebar_prob1           = new TH1D("numubar_nuebar_prob1", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob2           = new TH1D("numubar_nuebar_prob2", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob3           = new TH1D("numubar_nuebar_prob3", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob4           = new TH1D("numubar_nuebar_prob4", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob5           = new TH1D("numubar_nuebar_prob5", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob_vac1       = new TH1D("numubar_nuebar_prob_vac1", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob_vac2       = new TH1D("numubar_nuebar_prob_vac2", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob_vac3       = new TH1D("numubar_nuebar_prob_vac3", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob_vac4       = new TH1D("numubar_nuebar_prob_vac4", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * numubar_nuebar_prob_vac5       = new TH1D("numubar_nuebar_prob_vac5", " P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));


 TH1D * nue_numu_prob1           = new TH1D("nue_numu_prob1", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob2           = new TH1D("nue_numu_prob2", " P_{#nu_{e}->#nu_{#mu}} ",
					     util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob3           = new TH1D("nue_numu_prob3", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob4           = new TH1D("nue_numu_prob4", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob5           = new TH1D("nue_numu_prob5", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob_vac1       = new TH1D("nue_numu_prob_vac1", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob_vac2       = new TH1D("nue_numu_prob_vac2", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob_vac3       = new TH1D("nue_numu_prob_vac3", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob_vac4       = new TH1D("nue_numu_prob_vac4", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nue_numu_prob_vac5       = new TH1D("nue_numu_prob_vac5", " P_{#nu_{e}->#nu_{#mu}} ",
					      util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));


  TH1D * nuebar_numubar_prob1           = new TH1D("nuebar_numubar_prob1", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob2           = new TH1D("nuebar_numubar_prob2", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob3           = new TH1D("nuebar_numubar_prob3", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob4           = new TH1D("nuebar_numubar_prob4", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob5           = new TH1D("nuebar_numubar_prob5", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob_vac1       = new TH1D("nuebar_numubar_prob_vac1", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob_vac2       = new TH1D("nuebar_numubar_prob_vac2", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob_vac3       = new TH1D("nuebar_numubar_prob_vac3", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob_vac4       = new TH1D("nuebar_numubar_prob_vac4", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * nuebar_numubar_prob_vac5       = new TH1D("nuebar_numubar_prob_vac5", " P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}} ",
						    util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));


  NuOsc3Flavour fNuOsc3Flavour;

  //first value of deltacp
  double dcp_1 = -conv::kpi;
  ph.SetCurrValue(conv::kPhypCPV    , dcp_1);

  //loop over Mu-like Ereco bin edges
  for(int ir = 1; ir <= util::NumErecoBins(conv::kSKsFCFV1RingMuLike); ir++) {
   
    double Energy = 0.5 * (util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir-1) + util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir));
    //    cout << " Energy  = " << Energy << endl;
    
    fNuOsc3Flavour.SetOscParametersSinSqTAll(input_mixing12,input_mixing13,input_mixing23,input_dm12sq,input_dm23sq,dcp_1);
   
    double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(14,    14,    Energy);
    double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(-14,  -14,    Energy);
    double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(12,    12,    Energy);
    double p_nuebar_nuebar   = fNuOsc3Flavour.CalcMatterOscProb(-12,  -12,    Energy);
    double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(14,    12,    Energy);
    double p_numubar_nuebar  = fNuOsc3Flavour.CalcMatterOscProb(-14,  -12,    Energy);
    double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(12,    14,    Energy);
    double p_nuebar_numubar  = fNuOsc3Flavour.CalcMatterOscProb(-12,  -14,    Energy);
    
    numu_numu_prob1->Fill(Energy,p_numu_numu);
    numubar_numubar_prob1->Fill(Energy,p_numubar_numubar);
    nue_nue_prob1->Fill(Energy,p_nue_nue);
    nuebar_nuebar_prob1->Fill(Energy,p_nuebar_nuebar);
    numu_nue_prob1->Fill(Energy,p_numu_nue);
    numubar_nuebar_prob1->Fill(Energy,p_numubar_nuebar);
    nue_numu_prob1->Fill(Energy,p_nue_numu);
    nuebar_numubar_prob1->Fill(Energy,p_nuebar_numubar);

  }


  //first value of deltacp
  double dcp_2 = -conv::kpi/2.;
  ph.SetCurrValue(conv::kPhypCPV    , dcp_2);

  //loop over Mu-like Ereco bin edges
  for(int ir = 1; ir <= util::NumErecoBins(conv::kSKsFCFV1RingMuLike); ir++) {
   
    double Energy = 0.5 * (util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir-1) + util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir));
    //    cout << " Energy  = " << Energy << endl;
    
    fNuOsc3Flavour.SetOscParametersSinSqTAll(input_mixing12,input_mixing13,input_mixing23,input_dm12sq,input_dm23sq,dcp_2);
   
    double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(14,    14,    Energy);
    double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(-14,  -14,    Energy);
    double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(12,    12,    Energy);
    double p_nuebar_nuebar   = fNuOsc3Flavour.CalcMatterOscProb(-12,  -12,    Energy);
    double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(14,    12,    Energy);
    double p_numubar_nuebar  = fNuOsc3Flavour.CalcMatterOscProb(-14,  -12,    Energy);
    double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(12,    14,    Energy);
    double p_nuebar_numubar  = fNuOsc3Flavour.CalcMatterOscProb(-12,  -14,    Energy);
    
    numu_numu_prob2->Fill(Energy,p_numu_numu);
    numubar_numubar_prob2->Fill(Energy,p_numubar_numubar);
    nue_nue_prob2->Fill(Energy,p_nue_nue);
    nuebar_nuebar_prob2->Fill(Energy,p_nuebar_nuebar);
    numu_nue_prob2->Fill(Energy,p_numu_nue);
    numubar_nuebar_prob2->Fill(Energy,p_numubar_nuebar);
    nue_numu_prob2->Fill(Energy,p_nue_numu);
    nuebar_numubar_prob2->Fill(Energy,p_nuebar_numubar);

  }

  //first value of deltacp
  double dcp_3 = 0.0;
  ph.SetCurrValue(conv::kPhypCPV    , dcp_3);

  //loop over Mu-like Ereco bin edges
  for(int ir = 1; ir <= util::NumErecoBins(conv::kSKsFCFV1RingMuLike); ir++) {
   
    double Energy = 0.5 * (util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir-1) + util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir));
    //    cout << " Energy  = " << Energy << endl;
    
    fNuOsc3Flavour.SetOscParametersSinSqTAll(input_mixing12,input_mixing13,input_mixing23,input_dm12sq,input_dm23sq,dcp_3);
   
    double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(14,    14,    Energy);
    double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(-14,  -14,    Energy);
    double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(12,    12,    Energy);
    double p_nuebar_nuebar   = fNuOsc3Flavour.CalcMatterOscProb(-12,  -12,    Energy);
    double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(14,    12,    Energy);
    double p_numubar_nuebar  = fNuOsc3Flavour.CalcMatterOscProb(-14,  -12,    Energy);
    double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(12,    14,    Energy);
    double p_nuebar_numubar  = fNuOsc3Flavour.CalcMatterOscProb(-12,  -14,    Energy);
    
    numu_numu_prob3->Fill(Energy,p_numu_numu);
    numubar_numubar_prob3->Fill(Energy,p_numubar_numubar);
    nue_nue_prob3->Fill(Energy,p_nue_nue);
    nuebar_nuebar_prob3->Fill(Energy,p_nuebar_nuebar);
    numu_nue_prob3->Fill(Energy,p_numu_nue);
    numubar_nuebar_prob3->Fill(Energy,p_numubar_nuebar);
    nue_numu_prob3->Fill(Energy,p_nue_numu);
    nuebar_numubar_prob3->Fill(Energy,p_nuebar_numubar);

  }


  //first value of deltacp
  double dcp_4 = conv::kpi/2.;
  ph.SetCurrValue(conv::kPhypCPV    , dcp_4);

  //loop over Mu-like Ereco bin edges
  for(int ir = 1; ir <= util::NumErecoBins(conv::kSKsFCFV1RingMuLike); ir++) {
   
    double Energy = 0.5 * (util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir-1) + util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir));
    //    cout << " Energy  = " << Energy << endl;
    
    fNuOsc3Flavour.SetOscParametersSinSqTAll(input_mixing12,input_mixing13,input_mixing23,input_dm12sq,input_dm23sq,dcp_4);
   
    double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(14,    14,    Energy);
    double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(-14,  -14,    Energy);
    double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(12,    12,    Energy);
    double p_nuebar_nuebar   = fNuOsc3Flavour.CalcMatterOscProb(-12,  -12,    Energy);
    double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(14,    12,    Energy);
    double p_numubar_nuebar  = fNuOsc3Flavour.CalcMatterOscProb(-14,  -12,    Energy);
    double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(12,    14,    Energy);
    double p_nuebar_numubar  = fNuOsc3Flavour.CalcMatterOscProb(-12,  -14,    Energy);
    
    numu_numu_prob4->Fill(Energy,p_numu_numu);
    numubar_numubar_prob4->Fill(Energy,p_numubar_numubar);
    nue_nue_prob4->Fill(Energy,p_nue_nue);
    nuebar_nuebar_prob4->Fill(Energy,p_nuebar_nuebar);
    numu_nue_prob4->Fill(Energy,p_numu_nue);
    numubar_nuebar_prob4->Fill(Energy,p_numubar_nuebar);
    nue_numu_prob4->Fill(Energy,p_nue_numu);
    nuebar_numubar_prob4->Fill(Energy,p_nuebar_numubar);

  }


    //first value of deltacp
  double dcp_5 = conv::kpi;
  ph.SetCurrValue(conv::kPhypCPV    , dcp_5);

  //loop over Mu-like Ereco bin edges
  for(int ir = 1; ir <= util::NumErecoBins(conv::kSKsFCFV1RingMuLike); ir++) {
   
    double Energy = 0.5 * (util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir-1) + util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,ir));
    //    cout << " Energy  = " << Energy << endl;
    
    fNuOsc3Flavour.SetOscParametersSinSqTAll(input_mixing12,input_mixing13,input_mixing23,input_dm12sq,input_dm23sq,dcp_5);
   
    double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(14,    14,    Energy);
    double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(-14,  -14,    Energy);
    double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(12,    12,    Energy);
    double p_nuebar_nuebar   = fNuOsc3Flavour.CalcMatterOscProb(-12,  -12,    Energy);
    double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(14,    12,    Energy);
    double p_numubar_nuebar  = fNuOsc3Flavour.CalcMatterOscProb(-14,  -12,    Energy);
    double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(12,    14,    Energy);
    double p_nuebar_numubar  = fNuOsc3Flavour.CalcMatterOscProb(-12,  -14,    Energy);
    
    numu_numu_prob5->Fill(Energy,p_numu_numu);
    numubar_numubar_prob5->Fill(Energy,p_numubar_numubar);
    nue_nue_prob5->Fill(Energy,p_nue_nue);
    nuebar_nuebar_prob5->Fill(Energy,p_nuebar_nuebar);
    numu_nue_prob5->Fill(Energy,p_numu_nue);
    numubar_nuebar_prob5->Fill(Energy,p_numubar_nuebar);
    nue_numu_prob5->Fill(Energy,p_nue_numu);
    nuebar_numubar_prob5->Fill(Energy,p_nuebar_numubar);

  }

  w.Stop();
  cout << "Took " << w.CpuTime() << " to calculate the " << 5 * util::NumErecoBins(conv::kSKsFCFV1RingMuLike) << " oscillation probabilities" << endl;
 
  gStyle->SetOptStat(0);
  //  gStyle->SetFillColor(kWhite);

  TCanvas * c_numu_numu = new TCanvas("c_numu_numu","",20,20,500,500);
  TH1F * h2 = (TH1F*) c_numu_numu->DrawFrame(util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,0),-0.5,util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,util::NumErecoBins(conv::kSKsFCFV1RingMuLike)+1),1.5);
  //  h2->Draw("c");
  c_numu_numu->SetFillColor(kWhite);
  numu_numu_prob1->SetTitle("P_{#nu_{#mu}->#nu_{#mu}}");
  numu_numu_prob1->GetYaxis()->SetTitle("P_{#nu_{#mu}->#nu_{#mu}}");
  numu_numu_prob1->GetYaxis()->SetTitleOffset(1.0);
  numu_numu_prob1->GetYaxis()->SetRangeUser(-0.2,1.2);
  numu_numu_prob1->GetXaxis()->SetRangeUser(0.0,2.0);
  numu_numu_prob1->GetXaxis()->SetTitle("Energy (GeV)");
  numu_numu_prob1->SetLineColor(kRed);
  numu_numu_prob1->SetLineWidth(2);
  numu_numu_prob1->Draw("c");
  numu_numu_prob2->SetLineColor(kOrange);
  numu_numu_prob2->SetLineWidth(2);
  numu_numu_prob2->Draw("csame");
  numu_numu_prob3->SetLineColor(kBlack);
  numu_numu_prob3->SetLineWidth(2);
  numu_numu_prob3->Draw("csame");
  numu_numu_prob4->SetLineColor(kGreen);
  numu_numu_prob4->SetLineWidth(2);
  numu_numu_prob4->SetLineStyle(kDashed);
  numu_numu_prob4->Draw("csame");
  numu_numu_prob5->SetLineColor(kBlue);
  numu_numu_prob5->SetLineWidth(2);
  numu_numu_prob5->SetLineStyle(kDashed);
  numu_numu_prob5->Draw("csame");

  TLegend * legend = new TLegend(0.7, 0.7, 0.95, 0.95);
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  legend->AddEntry(numu_numu_prob1,"#delta_{CP} = -#pi","L");
  legend->AddEntry(numu_numu_prob2,"#delta_{CP} = -#pi/2","L");
  legend->AddEntry(numu_numu_prob3,"#delta_{CP} = 0","L");
  legend->AddEntry(numu_numu_prob4,"#delta_{CP} = #pi/2","L");
  legend->AddEntry(numu_numu_prob5,"#delta_{CP} = #pi","L");
  legend->Draw("csame");

  c_numu_numu->SaveAs("Osc_prob_dcp_numu_numu.eps");


  TCanvas * c_numubar_numubar = new TCanvas("c_numubar_numubar","",20,20,500,500);
  TH1F * h_numubar_numubar = (TH1F*) c_numubar_numubar->DrawFrame(util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,0),-0.1,util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,util::NumErecoBins(conv::kSKsFCFV1RingMuLike)+1),1.5);
  c_numubar_numubar->SetFillColor(kWhite);
  numubar_numubar_prob1->SetTitle("P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}}");
  numubar_numubar_prob1->GetYaxis()->SetTitle("P_{#bar{#nu_{#mu}}->#bar{#nu_{#mu}}}");
  numubar_numubar_prob1->GetYaxis()->SetTitleOffset(1.0);
  numubar_numubar_prob1->GetYaxis()->SetRangeUser(-0.2,1.2);
  numubar_numubar_prob1->GetXaxis()->SetRangeUser(0.0,2.0);
  numubar_numubar_prob1->GetXaxis()->SetTitle("Energy (GeV)");
  numubar_numubar_prob1->SetLineColor(kRed);
  numubar_numubar_prob1->SetLineWidth(2);
  numubar_numubar_prob1->Draw("c");
  numubar_numubar_prob2->SetLineColor(kOrange);
  numubar_numubar_prob2->SetLineWidth(2);
  numubar_numubar_prob2->Draw("csame");
  numubar_numubar_prob3->SetLineColor(kBlack);
  numubar_numubar_prob3->SetLineWidth(2);
  numubar_numubar_prob3->Draw("csame");
  numubar_numubar_prob4->SetLineColor(kGreen);
  numubar_numubar_prob4->SetLineWidth(2);
  numubar_numubar_prob4->SetLineStyle(kDashed);
  numubar_numubar_prob4->Draw("csame");
  numubar_numubar_prob5->SetLineColor(kBlue);
  numubar_numubar_prob5->SetLineWidth(2);
  numubar_numubar_prob5->SetLineStyle(kDashed);
  numubar_numubar_prob5->Draw("csame");

  legend->Draw("csame");

  c_numubar_numubar->SaveAs("Osc_prob_dcp_numubar_numubar.eps");


  TCanvas * c_nue_nue = new TCanvas("c_nue_nue","",20,20,500,500);
  TH1F * h_nue_nue = (TH1F*) c_nue_nue->DrawFrame(util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,0),-0.1,util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,util::NumErecoBins(conv::kSKsFCFV1RingMuLike)+1),1.5);
  c_nue_nue->SetFillColor(kWhite);
  nue_nue_prob1->SetTitle("P_{#nu_{e}->#nu_{e}}");
  nue_nue_prob1->GetYaxis()->SetTitle("P_{#nu_{e}->#nu_{e}}");
  nue_nue_prob1->GetYaxis()->SetTitleOffset(1.0);
  nue_nue_prob1->GetYaxis()->SetRangeUser(-0.2,1.2);
  nue_nue_prob1->GetXaxis()->SetRangeUser(0.0,2.0);
  nue_nue_prob1->GetXaxis()->SetTitle("Energy (GeV)");
  nue_nue_prob1->SetLineColor(kRed);
  nue_nue_prob1->SetLineWidth(2);
  nue_nue_prob1->Draw("c");
  nue_nue_prob2->SetLineColor(kOrange);
  nue_nue_prob2->SetLineWidth(2);
  nue_nue_prob2->Draw("csame");
  nue_nue_prob3->SetLineColor(kBlack);
  nue_nue_prob3->SetLineWidth(2);
  nue_nue_prob3->Draw("csame");
  nue_nue_prob4->SetLineColor(kGreen);
  nue_nue_prob4->SetLineWidth(2);
  nue_nue_prob4->SetLineStyle(kDashed);
  nue_nue_prob4->Draw("csame");
  nue_nue_prob5->SetLineColor(kBlue);
  nue_nue_prob5->SetLineWidth(2);
  nue_nue_prob5->SetLineStyle(kDashed);
  nue_nue_prob5->Draw("csame");

  legend->Draw("csame");

  c_nue_nue->SaveAs("Osc_prob_dcp_nue_nue.eps");



  TCanvas * c_nuebar_nuebar = new TCanvas("c_nuebar_nuebar","",20,20,500,500);
  TH1F * h_nuebar_nuebar = (TH1F*) c_nuebar_nuebar->DrawFrame(util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,0),-0.1,util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,util::NumErecoBins(conv::kSKsFCFV1RingMuLike)+1),1.5);
  c_nuebar_nuebar->SetFillColor(kWhite);
  nuebar_nuebar_prob1->SetTitle("P_{#bar{#nu_{e}}->#bar{#nu_{e}}}");
  nuebar_nuebar_prob1->GetYaxis()->SetTitle("P_{#bar{#nu_{e}}->#bar{#nu_{e}}}");
  nuebar_nuebar_prob1->GetYaxis()->SetTitleOffset(1.0);
  nuebar_nuebar_prob1->GetYaxis()->SetRangeUser(-0.2,1.2);
  nuebar_nuebar_prob1->GetXaxis()->SetRangeUser(0.0,2.0);
  nuebar_nuebar_prob1->GetXaxis()->SetTitle("Energy (GeV)");
  nuebar_nuebar_prob1->SetLineColor(kRed);
  nuebar_nuebar_prob1->SetLineWidth(2);
  nuebar_nuebar_prob1->Draw("c");
  nuebar_nuebar_prob2->SetLineColor(kOrange);
  nuebar_nuebar_prob2->SetLineWidth(2);
  nuebar_nuebar_prob2->Draw("csame");
  nuebar_nuebar_prob3->SetLineColor(kBlack);
  nuebar_nuebar_prob3->SetLineWidth(2);
  nuebar_nuebar_prob3->Draw("csame");
  nuebar_nuebar_prob4->SetLineColor(kGreen);
  nuebar_nuebar_prob4->SetLineWidth(2);
  nuebar_nuebar_prob4->SetLineStyle(kDashed);
  nuebar_nuebar_prob4->Draw("csame");
  nuebar_nuebar_prob5->SetLineColor(kBlue);
  nuebar_nuebar_prob5->SetLineWidth(2);
  nuebar_nuebar_prob5->SetLineStyle(kDashed);
  nuebar_nuebar_prob5->Draw("csame");

  legend->Draw("csame");

  c_nuebar_nuebar->SaveAs("Osc_prob_dcp_nuebar_nuebar.eps");



  TCanvas * c_numu_nue = new TCanvas("c_numu_nue","",20,20,500,500);
  TH1F * h_numu_nue = (TH1F*) c_numu_nue->DrawFrame(util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,0),-0.1,util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,util::NumErecoBins(conv::kSKsFCFV1RingMuLike)+1),1.5);
  c_numu_nue->SetFillColor(kWhite);
  numu_nue_prob1->SetTitle("P_{#nu_{#mu}->#nu_{e}}");
  numu_nue_prob1->GetYaxis()->SetTitle("P_{#nu_{#mu}->#nu_{e}}");
  numu_nue_prob1->GetYaxis()->SetTitleOffset(1.0);
  numu_nue_prob1->GetYaxis()->SetRangeUser(-0.2,1.2);
  numu_nue_prob1->GetXaxis()->SetRangeUser(0.0,2.0);
  numu_nue_prob1->GetXaxis()->SetTitle("Energy (GeV)");
  numu_nue_prob1->SetLineColor(kRed);
  numu_nue_prob1->SetLineWidth(2);
  numu_nue_prob1->Draw("c");
  numu_nue_prob2->SetLineColor(kOrange);
  numu_nue_prob2->SetLineWidth(2);
  numu_nue_prob2->Draw("csame");
  numu_nue_prob3->SetLineColor(kBlack);
  numu_nue_prob3->SetLineWidth(2);
  numu_nue_prob3->Draw("csame");
  numu_nue_prob4->SetLineColor(kGreen);
  numu_nue_prob4->SetLineWidth(2);
  numu_nue_prob4->SetLineStyle(kDashed);
  numu_nue_prob4->Draw("csame");
  numu_nue_prob5->SetLineColor(kBlue);
  numu_nue_prob5->SetLineWidth(2);
  numu_nue_prob5->SetLineStyle(kDashed);
  numu_nue_prob5->Draw("csame");

  legend->Draw("csame");

  c_numu_nue->SaveAs("Osc_prob_dcp_numu_nue.eps");


  TCanvas * c_numubar_nuebar = new TCanvas("c_numubar_nuebar","",20,20,500,500);
  TH1F * h_numubar_nuebar = (TH1F*) c_numubar_nuebar->DrawFrame(util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,0),-0.1,util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,util::NumErecoBins(conv::kSKsFCFV1RingMuLike)+1),1.5);
  c_numubar_nuebar->SetFillColor(kWhite);
  numubar_nuebar_prob1->SetTitle("P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}}");
  numubar_nuebar_prob1->GetYaxis()->SetTitle("P_{#bar{#nu_{#mu}}->#bar{#nu_{e}}}");
  numubar_nuebar_prob1->GetYaxis()->SetTitleOffset(1.0);
  numubar_nuebar_prob1->GetYaxis()->SetRangeUser(-0.2,1.2);
  numubar_nuebar_prob1->GetXaxis()->SetRangeUser(0.0,2.0);
  numubar_nuebar_prob1->GetXaxis()->SetTitle("Energy (GeV)");
  numubar_nuebar_prob1->SetLineColor(kRed);
  numubar_nuebar_prob1->SetLineWidth(2);
  numubar_nuebar_prob1->Draw("c");
  numubar_nuebar_prob2->SetLineColor(kOrange);
  numubar_nuebar_prob2->SetLineWidth(2);
  numubar_nuebar_prob2->Draw("csame");
  numubar_nuebar_prob3->SetLineColor(kBlack);
  numubar_nuebar_prob3->SetLineWidth(2);
  numubar_nuebar_prob3->Draw("csame");
  numubar_nuebar_prob4->SetLineColor(kGreen);
  numubar_nuebar_prob4->SetLineWidth(2);
  numubar_nuebar_prob4->SetLineStyle(kDashed);
  numubar_nuebar_prob4->Draw("csame");
  numubar_nuebar_prob5->SetLineColor(kBlue);
  numubar_nuebar_prob5->SetLineWidth(2);
  numubar_nuebar_prob5->SetLineStyle(kDashed);
  numubar_nuebar_prob5->Draw("csame");

  legend->Draw("csame");

  c_numubar_nuebar->SaveAs("Osc_prob_dcp_numubar_nuebar.eps");


   TCanvas * c_nue_numu = new TCanvas("c_nue_numu","",20,20,500,500);
  TH1F * h_nue_numu = (TH1F*) c_nue_numu->DrawFrame(util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,0),-0.1,util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,util::NumErecoBins(conv::kSKsFCFV1RingMuLike)+1),1.5);
  c_nue_numu->SetFillColor(kWhite);
  nue_numu_prob1->SetTitle("P_{#nu_{e}->#nu_{#mu}}");
  nue_numu_prob1->GetYaxis()->SetTitle("P_{#nu_{e}->#nu_{#mu}}");
  nue_numu_prob1->GetYaxis()->SetTitleOffset(1.0);
  nue_numu_prob1->GetYaxis()->SetRangeUser(-0.2,1.2);
  nue_numu_prob1->GetXaxis()->SetRangeUser(0.0,2.0);
  nue_numu_prob1->GetXaxis()->SetTitle("Energy (GeV)");
  nue_numu_prob1->SetLineColor(kRed);
  nue_numu_prob1->SetLineWidth(2);
  nue_numu_prob1->Draw("c");
  nue_numu_prob2->SetLineColor(kOrange);
  nue_numu_prob2->SetLineWidth(2);
  nue_numu_prob2->Draw("csame");
  nue_numu_prob3->SetLineColor(kBlack);
  nue_numu_prob3->SetLineWidth(2);
  nue_numu_prob3->Draw("csame");
  nue_numu_prob4->SetLineColor(kGreen);
  nue_numu_prob4->SetLineWidth(2);
  nue_numu_prob4->SetLineStyle(kDashed);
  nue_numu_prob4->Draw("csame");
  nue_numu_prob5->SetLineColor(kBlue);
  nue_numu_prob5->SetLineWidth(2);
  nue_numu_prob5->SetLineStyle(kDashed);
  nue_numu_prob5->Draw("csame");

  legend->Draw("csame");

  c_nue_numu->SaveAs("Osc_prob_dcp_nue_numu.eps");


  TCanvas * c_nuebar_numubar = new TCanvas("c_nuebar_numubar","",20,20,500,500);
  TH1F * h_nuebar_numubar = (TH1F*) c_nuebar_numubar->DrawFrame(util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,0),-0.1,util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike,util::NumErecoBins(conv::kSKsFCFV1RingMuLike)+1),1.5);
  c_nuebar_numubar->SetFillColor(kWhite);
  nuebar_numubar_prob1->SetTitle("P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}}");
  nuebar_numubar_prob1->GetYaxis()->SetTitle("P_{#bar{#nu_{e}}->#bar{#nu_{#mu}}}");
  nuebar_numubar_prob1->GetYaxis()->SetTitleOffset(1.0);
  nuebar_numubar_prob1->GetYaxis()->SetRangeUser(-0.2,1.2);
  nuebar_numubar_prob1->GetXaxis()->SetRangeUser(0.0,2.0);
  nuebar_numubar_prob1->GetXaxis()->SetTitle("Energy (GeV)");
  nuebar_numubar_prob1->SetLineColor(kRed);
  nuebar_numubar_prob1->SetLineWidth(2);
  nuebar_numubar_prob1->Draw("c");
  nuebar_numubar_prob2->SetLineColor(kOrange);
  nuebar_numubar_prob2->SetLineWidth(2);
  nuebar_numubar_prob2->Draw("csame");
  nuebar_numubar_prob3->SetLineColor(kBlack);
  nuebar_numubar_prob3->SetLineWidth(2);
  nuebar_numubar_prob3->Draw("csame");
  nuebar_numubar_prob4->SetLineColor(kGreen);
  nuebar_numubar_prob4->SetLineWidth(2);
  nuebar_numubar_prob4->SetLineStyle(kDashed);
  nuebar_numubar_prob4->Draw("csame");
  nuebar_numubar_prob5->SetLineColor(kBlue);
  nuebar_numubar_prob5->SetLineWidth(2);
  nuebar_numubar_prob5->SetLineStyle(kDashed);
  nuebar_numubar_prob5->Draw("csame");

  legend->Draw("csame");

  c_nuebar_numubar->SaveAs("Osc_prob_dcp_nuebar_numubar.eps");

  const char * filename_out = "file_osc_prob_deltacp.root";
  TFile * f_out = new TFile(filename_out, "recreate");
  numu_numu_prob1->Write();
  numu_numu_prob2->Write();
  numu_numu_prob3->Write();
  numu_numu_prob4->Write();
  numu_numu_prob5->Write();

  numubar_numubar_prob1->Write();
  numubar_numubar_prob2->Write();
  numubar_numubar_prob3->Write();
  numubar_numubar_prob4->Write();
  numubar_numubar_prob5->Write();

  nue_nue_prob1->Write();
  nue_nue_prob2->Write();
  nue_nue_prob3->Write();
  nue_nue_prob4->Write();
  nue_nue_prob5->Write();

  nuebar_nuebar_prob1->Write();
  nuebar_nuebar_prob2->Write();
  nuebar_nuebar_prob3->Write();
  nuebar_nuebar_prob4->Write();
  nuebar_nuebar_prob5->Write();

  numu_nue_prob1->Write();
  numu_nue_prob2->Write();
  numu_nue_prob3->Write();
  numu_nue_prob4->Write();
  numu_nue_prob5->Write();

  numubar_nuebar_prob1->Write();
  numubar_nuebar_prob2->Write();
  numubar_nuebar_prob3->Write();
  numubar_nuebar_prob4->Write();
  numubar_nuebar_prob5->Write();

  nue_numu_prob1->Write();
  nue_numu_prob2->Write();
  nue_numu_prob3->Write();
  nue_numu_prob4->Write();
  nue_numu_prob5->Write();

  nuebar_numubar_prob1->Write();
  nuebar_numubar_prob2->Write();
  nuebar_numubar_prob3->Write();
  nuebar_numubar_prob4->Write();
  nuebar_numubar_prob5->Write();

}
