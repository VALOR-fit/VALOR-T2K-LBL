//
// ------------------------------------------------------------------------------------------------
// Plot numu -> numu survival probability (computed in 2-flavour framework) as a function of
// neutrino energy for the given mixing23 and dm23sq parameters.
// Superimpose the range of 3-flavour framework predictions obtained for the same mixing23 and 
// dm23sq params and for random variations of mixing12, mixing13, dm23sq and delta_cp.
//
// C.Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void plot_numu_surv_prob_2f_vs_3frange(double mixing23=1.0, double dm23sq=2.32E-3)
{
  using namespace std;
  using namespace t2kana;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Constants
  const int    numu                 = conv::kPdgNuMu;
  const int    nrand                = 4000;
  const int    ne                   = 500;
  const double emin                 = 0.250;
  const double emax                 = 10.000;
  const double de                   = (emax-emin)/(ne-1);

  NuOsc3Flavour calc_osc3f;

  double e[ne];
  double numu_surv_prob_2f[ne];
  double numu_surv_prob_3fv_upper[ne];
  double numu_surv_prob_3fv_lower[ne];
  double numu_surv_prob_3fm_upper[ne];
  double numu_surv_prob_3fm_lower[ne];

  for(int i = 0; i < ne; i++) {

    e[i] = emin + i*de;
    //cout << "E = " << e[i] << " GeV" << endl;

    numu_surv_prob_2f[i] = util::NuMuSurvProbOsc2Flavour(mixing23, dm23sq, e[i]);

    numu_surv_prob_3fv_upper[i] = -999;
    numu_surv_prob_3fv_lower[i] = +999;
    numu_surv_prob_3fm_upper[i] = -999;
    numu_surv_prob_3fm_lower[i] = +999;

    for(int irand = 0; irand < nrand; irand++) {

       double dm12sq   = util::RandomizeDm12sq();
       double mixing12 = util::RandomizeMixing12();
       double mixing13 = util::RandomizeMixing13();
       double delta_cp = util::RandomizeDeltaCP();

       calc_osc3f.SetOscParameters(mixing12,mixing13,mixing23,dm12sq,dm23sq,delta_cp);

       double p3fv = calc_osc3f.CalcVacuumOscProb(numu,numu,e[i]);
       if(p3fv < 0. || p3fv > 1. || TMath::IsNaN(p3fv)) {
	 cout << "WARNING! P_{numu->numu} ("
              << "E = " << e[i] << " GeV"
              << ", mixing12 = "  << mixing12 
              << ", mixing13 = "  << mixing13 
              << ", mixing23 = "  << mixing23 
              << ", Dm12^2 = "    << dm12sq << " eV^2"
              << ", |Dm23^2| = "  << dm23sq << " eV^2"
              << ", delta_cp = "  << delta_cp << " rad) = " << p3fv << endl;
         continue;
       }
       numu_surv_prob_3fv_upper[i] = TMath::Max(numu_surv_prob_3fv_upper[i], p3fv);
       numu_surv_prob_3fv_lower[i] = TMath::Min(numu_surv_prob_3fv_lower[i], p3fv);

       double p3fm = calc_osc3f.CalcMatterOscProb(numu,numu,e[i]);
       if(p3fm < 0. || p3fm > 1. || TMath::IsNaN(p3fm)) {
	 cout << "WARNING! P_{numu->numu} ("
              << "E = " << e[i] << " GeV"
              << ", mixing12 = "  << mixing12 
              << ", mixing13 = "  << mixing13 
              << ", mixing23 = "  << mixing23 
              << ", Dm12^2 = "    << dm12sq << " eV^2"
              << ", |Dm23^2| = "  << dm23sq << " eV^2"
              << ", delta_cp = "  << delta_cp << " rad) = " << p3fm << endl;
         continue;
       }
       numu_surv_prob_3fm_upper[i] = TMath::Max(numu_surv_prob_3fm_upper[i], p3fm);
       numu_surv_prob_3fm_lower[i] = TMath::Min(numu_surv_prob_3fm_lower[i], p3fm);
    }//irand
  }//i

  TGraph * gr_numu_surv_prob_2f         = new TGraph(ne, e, numu_surv_prob_2f);
  TGraph * gr_numu_surv_prob_3fv_upper  = new TGraph(ne, e, numu_surv_prob_3fv_upper);
  TGraph * gr_numu_surv_prob_3fv_lower  = new TGraph(ne, e, numu_surv_prob_3fv_lower);
  TGraph * gr_numu_surv_prob_3fm_upper  = new TGraph(ne, e, numu_surv_prob_3fm_upper);
  TGraph * gr_numu_surv_prob_3fm_lower  = new TGraph(ne, e, numu_surv_prob_3fm_lower);

  TCanvas * c1 = new TCanvas("c","",10,10,800,600);
  TLegend * legend1 = new TLegend(0.3, 0.2, 0.85, 0.5);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);
  TH1F * h1 = (TH1F*) c1->DrawFrame(emin,0.,emax,1.);
  h1->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  h1->GetYaxis()->SetTitle("P_{#nu_{#mu} #rightarrow #nu_{#mu}}");
  gr_numu_surv_prob_3fv_upper->SetLineStyle(kDashed);
  gr_numu_surv_prob_3fv_lower->SetLineStyle(kDashed);
  gr_numu_surv_prob_3fm_upper->SetLineColor(kRed);
  gr_numu_surv_prob_3fm_lower->SetLineColor(kRed);
  gr_numu_surv_prob_3fm_upper->SetLineStyle(kDotted);
  gr_numu_surv_prob_3fm_lower->SetLineStyle(kDotted);
  legend1->SetHeader(Form("|#Deltam_{23}^{2}|=%.3e eV^{2}, sin^{2}2#theta_{23}=%.3e",dm23sq,mixing23));
  legend1->AddEntry(gr_numu_surv_prob_2f,        "2-flavour", "L");
  legend1->AddEntry(gr_numu_surv_prob_3fv_upper, "3-flavour vaccum (range for #theta_{12}, #theta_{13}, #delta_{CP}, #Deltam_{12}^{2} variation)", "L");
  legend1->AddEntry(gr_numu_surv_prob_3fm_upper, "3-flavour matter (range for #theta_{12}, #theta_{13}, #delta_{CP}, #Deltam_{12}^{2} variation)", "L");
  legend1->Draw();
  gr_numu_surv_prob_2f->Draw("l");
//  gr_numu_surv_prob_3fv_upper->Draw("l");
//  gr_numu_surv_prob_3fv_lower->Draw("l");
  gr_numu_surv_prob_3fm_upper->Draw("l");
  gr_numu_surv_prob_3fm_lower->Draw("l");
}
