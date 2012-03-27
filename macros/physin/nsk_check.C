//
// ------------------------------------------------------------------------------------------------
// Generates histograms of nominal SK spectrum with/without BANFF tuning
// with either Numu disappearance, three flavour oscillatin or no oscillation
//
// Adapted from test_nsk.C by Steve Dennis
// ------------------------------------------------------------------------------------------------

#include <iomanip>
using namespace std;
using namespace t2kana;
using namespace t2kana::conv;
  
void nsk_check
(
  bool withbanff  = true,
  PhysHypothesis_t phyp = kHypNuMuDisapp2F,
  double mixing13 = 0.1,
  double mixing23 = 1,
  double mixing12 = 0.86,
  double dm23sq = 2.32E-3,
  double dm12sq = 7.58E-5,
  double deltacp = 0,
  t2kana::conv::NuMassHierarchy_t hier = t2kana::conv::kNmhNormal, //or t2kana::conv::kNMhInverted
  const char * config_file = "analysis_configuration.xml"
)
{
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  if(hier == t2kana::conv::kNmhInverted)
    dm23sq = -dm23sq;

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  if (!withbanff) {
    skexp.GetBANFF()->RestoreNominalMC();
  }
  
  // generate predictions
  if (phyp == kHypNuMuDisapp2F || phyp == kHypOsc3F) {
    // fit params
    PhysHypothesis ph(phyp);

    // generate predictions
    ph.SetCurrValue(t2kana::conv::kPhypMixing23,      mixing23);
    ph.SetCurrValue(t2kana::conv::kPhypMixing13,      mixing13);
    ph.SetCurrValue(t2kana::conv::kPhypMixing12,      mixing12);
    ph.SetCurrValue(t2kana::conv::kPhypDm23sq,        dm23sq);
    ph.SetCurrValue(t2kana::conv::kPhypDm12sq,        dm12sq);
    ph.SetCurrValue(t2kana::conv::kPhypDeltaCP,       deltacp);
  
    skexp.Generate(ph);
  }
  else { /* Adapted from SuperkExpectation */
    // update SK MC templates
    skexp.GetSuperkMC()->UpdateTemplates();

    // update Ereco spectrum prediction
    skexp.GetSuperkMC()->UpdateErecoSpectra();

    // copy histograms
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
      const int nbins = util::NumErecoBins(s);
      for(int ibin = 1; ibin <= nbins; ibin++) {
         double ns = skexp.GetSuperkMC()->CurrErecoSpectrum(s)->GetBinContent(ibin);
         skexp.CurrPredictedErecoSpectrum(s)->SetBinContent(ibin, ns);
         for(int im = 0; im < knm; im++) {
             SuperkMode_t m = kSKReactionModes[is][im];
             if(m == kSKmUndefined) continue;
             double nsm = skexp.GetSuperkMC()->CurrErecoSpectrum(s,m)->GetBinContent(ibin);
             skexp.CurrPredictedErecoSpectrum(s,m)->SetBinContent(ibin, nsm);       
         }//m	
      }//ibin
    }//s
  }

  // integrate Ereco distribution
  double nsk = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
  double nskmucc = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE)->Integral();
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi)->Integral();
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth_mu)->Integral();
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE)->Integral();
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCC1pi)->Integral();
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCoth_mu)->Integral();
  double nsknc = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuNC)->Integral();
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarNC)->Integral();
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuENC)->Integral();
  double nskebeamcc = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCQE)->Integral();
  nskebeamcc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECC1pi)->Integral();
  nskebeamcc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCoth_mu)->Integral();
  double nskesigcc = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCQE)->Integral();
  nskesigcc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECC1pi)->Integral();
  nskesigcc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCoth)->Integral();


  //output in form of table 16, TN52
  //~ cout << "Hierarchy is " << t2kana::util::AsString(hier) << " dm23sq is " << dm23sq << endl;
  //~ cout << "     " << nskmucc+nskebeamcc+nsknc << " beam total" << endl;
  //~ cout << "     " << nskmucc << " numu+numubar CC" << endl;
  //~ cout << "     " << nskebeamcc << " nue beam CC" << endl;
  //~ cout << "     " << nsknc << " all NC" << endl;
  //~ cout << "     " << nskesigcc << " nue signal CC" << endl;
  //~ cout << "Number of predicted events at SK is: " << nsk << "     = " << nskmucc+nsknc+nskebeamcc+nskesigcc << endl;
  
  ostringstream outfile_name;
  
  string phyp_short;
  if (phyp == kHypOsc3F) {
    phyp_short = "3A";
  }
  else if (phyp == kHypNuMuDisapp2F) {
    phyp_short = "2A";
  }
  else {
    phyp_short = "N";
  }
  
  outfile_name << "N-" << phyp_short << "-" << (withbanff ? "W" : "N");
  TFile outfile( (outfile_name.str() + ".root").c_str() , "RECREATE");
  
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Write("Total");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE)->Write("NuMuCCQE");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi)->Write("NuMuCC1pi");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth_mu)->Write("NuMuCCoth_mu");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE)->Write("NuMuBarCCQE");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCC1pi)->Write("NuMuBarCC1pi");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCoth_mu)->Write("NuMuBarCCoth_mu");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuNC)->Write("NuMuNC");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarNC)->Write("NuMuBarNC");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuENC)->Write("NuENC");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCQE)->Write("NuECCQE");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECC1pi)->Write("NuECC1pi");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCoth_mu)->Write("NuECCoth_mu");
  //~ skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCQE)->Write("SigNuECCQE");
  //~ skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECC1pi)->Write("SigNuECC1pi");
  //~ skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCoth)->Write("SigNuECCoth");
  
  TH1D * tot_hist = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike);
  
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  int highn = tot_hist->GetNbinsX();
  double emax  = tot_hist->GetBinCenter(highn)+ tot_hist->GetBinWidth(highn)/2.;
  
  TH1F * hframe = (TH1F*)c->DrawFrame(0,0, emax, tot_hist->GetMaximum());
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Number of events");
  c->SetLogy();
  
  tot_hist->Draw();
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth_mu)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCC1pi)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCoth_mu)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuNC)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarNC)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuENC)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCQE)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECC1pi)->Draw("same");
  skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCoth_mu)->Draw("same");
  
  c->Update();
  c->Write("All");
  
  ofstream returnstr( (outfile_name.str() + ".txt").c_str() );
  returnstr << "EReco_min\tEreco_max\tn_events\n";
  
  //~ returnstr.precision(5);
  returnstr.setf(ios::fixed,ios::floatfield);

  for (int ii = 1 ; ii <= tot_hist->GetNbinsX() ; ii++) {
    
    double bin_width = tot_hist->GetBinWidth(ii);
    double bin_center = tot_hist->GetBinCenter(ii);
    double bin_min = bin_center - bin_width/2.;
    double bin_max = bin_center + bin_width/2.;;
    
    double n_events = tot_hist->GetBinContent(ii);
    
    returnstr << setprecision(2) << bin_min << "         \t" << bin_max
              << "         \t"   << setprecision(5) << n_events << endl;
              
    //~ returnstr.unsetf(ios::floatfield);
  }
  
  
}
