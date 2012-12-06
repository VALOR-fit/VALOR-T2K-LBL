//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of 
// cross section tweaking dials.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_sk_ereco_spectrum_twk_xsec
(
  const char * hypothesis         = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12     = -99,                            /* */
  double       input_mixing13     = -99,                            /* */
  double       input_mixing23     = -99,                            /* */
  double       input_deltacp      = -99,                            /* */
  double       input_dm12sq       = -99,                            /* */
  double       input_dm23sq       = -99,                            /* */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }
  
  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Physics hypothesis: \n" << ph << endl;
  cout << *Opt::Instance() << endl;


  // Tweaks
  const int ntwkmodes = 13;
  string twklabel[ntwkmodes] = {
    "FermiMom",
    "SpecFunc",
    "Wshape",
    "ShapeCCoth",
    "PiLessDelta",
    "ENuShape",
    "NormCCQEE1",
    "NormCCQEE2",
    "NormCC1piE1",
    "NormCCcoh",
    "NormNCcoh",
    "NormNCoth",
    "NormNuEToNuMu"
  };
  const int ntwkdials = 5;
  double twkdial[ntwkdials] = { -1.0, -0.5, 0.0, +0.5, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2,       1,       1      };
  int hcolor [ntwkdials] = { kBlue,   kBlue,   kBlack,  kRed,    kRed   };
  int hstyle [ntwkdials] = { kDotted, kDashed, kSolid,  kDashed, kDotted};
  
  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);
  
  // Calculate spectra
  TH1D * hEreco   [ntwkmodes][ntwkdials];
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
       // reset previously set tweaks
       skexp.ResetSystematicVariations(); 
       // tweak spectrum generation model
       if      (imode==0)  skexp.GetNuXSec()->TweakFermiMom      ( twkdial[idial] );
       else if (imode==1)  skexp.GetNuXSec()->TweakSpecFunc      ( twkdial[idial] );
       else if (imode==2)  skexp.GetNuXSec()->TweakWshape        ( twkdial[idial] );
       else if (imode==3)  skexp.GetNuXSec()->TweakShapeCCoth    ( twkdial[idial] );
       else if (imode==4)  skexp.GetNuXSec()->TweakPiLessDelta   ( twkdial[idial] );
       else if (imode==5)  skexp.GetNuXSec()->TweakENuShape      ( twkdial[idial] );
       else if (imode==6)  skexp.GetNuXSec()->TweakNormCCQEE1    ( twkdial[idial] );
       else if (imode==7)  skexp.GetNuXSec()->TweakNormCCQEE2    ( twkdial[idial] );
       else if (imode==8)  skexp.GetNuXSec()->TweakNormCC1piE1   ( twkdial[idial] );
       else if (imode==9)  skexp.GetNuXSec()->TweakNormCCcoh     ( twkdial[idial] );
       else if (imode==10) skexp.GetNuXSec()->TweakNormNCcoh     ( twkdial[idial] );
       else if (imode==11) skexp.GetNuXSec()->TweakNormNCoth     ( twkdial[idial] );
       else if (imode==12) skexp.GetNuXSec()->TweakNormNuEToNuMu ( twkdial[idial] );
       // generate predictions  
       skexp.Generate(ph);
       // store spectrum
       hEreco[imode][idial] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
       hEreco[imode][idial]->SetDirectory(0);
    }//idial
  }//imode
  // get fractional err
  TH1D * hErecoErr[ntwkmodes][ntwkdials];
  double maxerrpos = -1.;
  double maxerrneg = +1.;
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
       hErecoErr[imode][idial] = new TH1D("","",conv::kNErecoBins1Rmu,conv::kErecoBinEdges1Rmu);
       hErecoErr[imode][idial]->SetDirectory(0);
       for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
         double ncurr  = hEreco[imode][idial]->GetBinContent(ibin);
         double ndeflt = hEreco[imode][deflt]->GetBinContent(ibin);
         if(ndeflt>0) {
           double err = (ncurr-ndeflt)/ndeflt;
           maxerrpos = TMath::Max(maxerrpos, err);
           maxerrneg = TMath::Min(maxerrneg, err);          
           hErecoErr[imode][idial]->SetBinContent(ibin, (ncurr-ndeflt)/ndeflt);
         }
       }//ibin
    }//idial
  }//imode

  //
  // Draw
  //
  for(int imode = 0; imode < ntwkmodes; imode++) {
    for(int idial = 0; idial < ntwkdials; idial++) {
       hEreco   [imode][idial]->SetLineWidth(hwidth[idial]);
       hEreco   [imode][idial]->SetLineColor(hcolor[idial]);
       hEreco   [imode][idial]->SetLineStyle(hstyle[idial]);
       hErecoErr[imode][idial]->SetLineWidth(hwidth[idial]);
       hErecoErr[imode][idial]->SetLineColor(hcolor[idial]);
       hErecoErr[imode][idial]->SetLineStyle(hstyle[idial]);
    }
  }
  TCanvas * c[ntwkmodes];
  gStyle->SetOptStat(0);
  for(int imode = 0; imode < ntwkmodes; imode++) {
    c[imode] = new TCanvas(Form("c%d",imode),"",20,20,500,500);
    c[imode]->Divide(1,2);
    c[imode]->cd(1);
    TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetHeader(twklabel[imode].c_str());
    cout << " *** " << twklabel[imode] << endl;
    for(int idial = 0; idial < ntwkdials; idial++) {
       const char * opt   = (idial==0) ? "" : "same";
       const char * label = Form("tweak dial = %.2f", twkdial[idial]);
       if(idial==0) {
        hEreco[imode][idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
        hEreco[imode][idial]->GetYaxis()->SetTitle("Number of events");
       }
       hEreco[imode][idial]->Draw(opt);
       legend->AddEntry(hEreco[imode][idial], label, "L");
       cout << " * " << label << " : NumEvents = " << hEreco[imode][idial]->Integral() << endl;
    }//idial
    legend->Draw();
    c[imode]->cd(2);
    TH1F * hframe = (TH1F*) c[imode]->GetPad(2)->DrawFrame(
       0., 1.2*maxerrneg, hErecoErr[imode][deflt]->GetXaxis()->GetXmax(), 1.2*maxerrpos);
    hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
    hframe->GetYaxis()->SetTitle("Fractional variation");
    hframe->Draw();
    for(int idial = 0; idial < ntwkdials; idial++) {
      if(idial==deflt) continue;
      hErecoErr[imode][idial]->Draw("same");
    }
    c[imode]->Update();
  }//imode

  // Write-out in file
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_xsec-%s.root",oscparams.c_str()),"recreate");
  for(int imode = 0; imode < ntwkmodes; imode++) {
     for(int idial = 0; idial < ntwkdials; idial++) {
       const char * syst = "";
       syst = twklabel[imode].c_str();
       const char * dial = "";
       if( TMath::Abs(twkdial[idial]-1.0) < conv::kASmallNum ) dial="plus_one_sigma";
       if( TMath::Abs(twkdial[idial]-0.5) < conv::kASmallNum ) dial="plus_half_sigma";
       if( TMath::Abs(twkdial[idial]    ) < conv::kASmallNum ) dial="nominal";
       if( TMath::Abs(twkdial[idial]+1.0) < conv::kASmallNum ) dial="minus_one_sigma";
       if( TMath::Abs(twkdial[idial]+0.5) < conv::kASmallNum ) dial="minus_half_sigma";
       const char * title = Form("sk_ereco_spectrum_twk_%s_%s",syst,dial);
       hEreco[imode][idial]->Write(title);
     }
  }
  f.Close();

}
