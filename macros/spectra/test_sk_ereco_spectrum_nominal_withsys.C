//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include "Conventions.h"
#include "Utilities.h"
#include "SuperkExpectation.h"
#include "Opt.h"
#include "BANFF.h"
#include "FitParamList.h"

#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH1.h"
#include "TSystem.h"
#include "TString.h"

#ifndef __CINT__
void test_sk_ereco_spectrum_nominal_withsyst
(
  const char * sk_sample            ,
  //const char * hypothesis           ,
  double       input_mixing12       ,
  double       input_mixing13       ,
  double       input_mixing23       ,
  double       input_deltacp        ,
  double       input_dm12sq         ,
  double       input_dm23sq         ,
  bool         disable_banff_tuning ,
  const char * config_file          ,
  bool         print_spectrum       ,
  //sk
  double       input_sk_escale ,
  double       input_sk_ccqeoth ,
  double       input_sk_ccqerc0 ,
  double       input_sk_ccqerc1 ,
  double       input_sk_ccqerc2 ,
  double       input_sk_ccnqe ,
  double       input_sk_nc ,
  double       input_sk_nue ,
  //banff flux
  double       input_banff_00 ,
  double       input_banff_01 ,
  double       input_banff_02 ,
  double       input_banff_03 ,
  double       input_banff_04 ,
  double       input_banff_05 ,
  double       input_banff_06 ,
  double       input_banff_07 ,
  double       input_banff_08 ,
  double       input_banff_09 ,
  double       input_banff_10 ,
  double       input_banff_11 ,
  double       input_banff_12 ,
  double       input_banff_13 ,
  double       input_banff_14 ,
  double       input_banff_15 ,
  //banff xsec
  double       input_banff_16 ,
  double       input_banff_17 ,
  double       input_banff_18 ,
  double       input_banff_19 ,
  double       input_banff_20 ,
  double       input_banff_21 ,
  double       input_banff_22 ,
  //xsec
  double       input_xsec_pf ,
  double       input_xsec_sf ,
  double       input_xsec_wshape ,
  double       input_xsec_ccoth ,
  double       input_xsec_piless ,
  double       input_xsec_binde ,
  double       input_xsec_cccoh ,
  double       input_xsec_nc1pic ,
  double       input_xsec_ncoth ,
  double       input_xsec_nuetonumu ,
  double       input_xsec_antinutonu ,
  //fsi
  double       input_fsi_0 ,
  double       input_fsi_1 ,
  double       input_fsi_2 ,
  double       input_fsi_3 ,
  double       input_fsi_4 ,
  double       input_fsi_5
 );
#endif

void test_sk_ereco_spectrum_nominal_withsyst
(
  const char * sk_sample            = "1Rmu",                         /* options: "1Rmu", "1Re" */
  //const char * hypothesis           = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12       = -99,                            /* */
  double       input_mixing13       = -99,                            /* */
  double       input_mixing23       = -99,                            /* */
  double       input_deltacp        = -99,                            /* */
  double       input_dm12sq         = -99,                            /* */
  double       input_dm23sq         = -99,                            /* */
  bool         disable_banff_tuning = false,                          /* */
  const char * config_file          = "analysis_configuration.xml",   /* */
  bool         print_spectrum       = false,                          /* */
  //sk
  double       input_sk_escale = 0,
  double       input_sk_ccqeoth = 0,
  double       input_sk_ccqerc0 = 0,
  double       input_sk_ccqerc1 = 0,
  double       input_sk_ccqerc2 = 0,
  double       input_sk_ccnqe = 0,
  double       input_sk_nc = 0,
  double       input_sk_nue = 0,
  //banff flux
  double       input_banff_00 = 0,
  double       input_banff_01 = 0,
  double       input_banff_02 = 0,
  double       input_banff_03 = 0,
  double       input_banff_04 = 0,
  double       input_banff_05 = 0,
  double       input_banff_06 = 0,
  double       input_banff_07 = 0,
  double       input_banff_08 = 0,
  double       input_banff_09 = 0,
  double       input_banff_10 = 0,
  double       input_banff_11 = 0,
  double       input_banff_12 = 0,
  double       input_banff_13 = 0,
  double       input_banff_14 = 0,
  double       input_banff_15 = 0,
  //banff xsec
  double       input_banff_16 = 0,
  double       input_banff_17 = 0,
  double       input_banff_18 = 0,
  double       input_banff_19 = 0,
  double       input_banff_20 = 0,
  double       input_banff_21 = 0,
  double       input_banff_22 = 0,
  //xsec
  double       input_xsec_pf = 0,
  double       input_xsec_sf = 0,
  double       input_xsec_wshape = 0,
  double       input_xsec_ccoth = 0,
  double       input_xsec_piless = 0,
  double       input_xsec_binde = 0,
  double       input_xsec_cccoh = 0,
  double       input_xsec_nc1pic = 0,
  double       input_xsec_ncoth = 0,
  double       input_xsec_nuetonumu = 0,
  double       input_xsec_antinutonu = 0,
  //fsi
  double       input_fsi_0 = 0,
  double       input_fsi_1 = 0,
  double       input_fsi_2 = 0,
  double       input_fsi_3 = 0,
  double       input_fsi_4 = 0,
  double       input_fsi_5 = 0
)
{
  using namespace std;
  using namespace t2kana;

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
  FitParamList * fparm = FitParamList::NamedParameterList("nuosc3f_numudisapp_fix4+nuisance");
  if(fparm->PhysicsHypothesis() == conv::kHypNuMuDisapp2F) {
    fparm->Update(conv::kPhypMixing23,  input_mixing23);
    fparm->Update(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(fparm->PhysicsHypothesis() == conv::kHypOsc3F) {
    fparm->Update(conv::kPhypMixing12,  input_mixing12);
    fparm->Update(conv::kPhypMixing13,  input_mixing13);
    fparm->Update(conv::kPhypMixing23,  input_mixing23);
    fparm->Update(conv::kPhypDeltaCP,   input_deltacp);
    fparm->Update(conv::kPhypDm12sq,    input_dm12sq);
    fparm->Update(conv::kPhypDm23sq,    input_dm23sq);
  }
  //set systematics
  fparm->Update(conv::kNuisXSecSpecFunc, 0);//spec func starts at 0.5!
  //sk
  fparm->Update(conv::kNuisSuperkEnergyScale,input_sk_escale);
  fparm->Update(conv::kNuisSuperkEff1RmuCCQEoth,input_sk_ccqeoth);
  fparm->Update(conv::kNuisSuperkEff1RmuCCQErc0u,input_sk_ccqerc0);
  fparm->Update(conv::kNuisSuperkEff1RmuCCQErc1u,input_sk_ccqerc1);
  fparm->Update(conv::kNuisSuperkEff1RmuCCQErc2u,input_sk_ccqerc2);
  fparm->Update(conv::kNuisSuperkEff1RmuCCnonQE,input_sk_ccnqe);
  fparm->Update(conv::kNuisSuperkEff1RmuNC,input_sk_nc);
  fparm->Update(conv::kNuisSuperkEff1RmuNueCC,input_sk_nue);
  //banff flux
  fparm->Update(conv::kNuisBANFF00,input_banff_00);
  fparm->Update(conv::kNuisBANFF01,input_banff_01);
  fparm->Update(conv::kNuisBANFF02,input_banff_02);
  fparm->Update(conv::kNuisBANFF03,input_banff_03);
  fparm->Update(conv::kNuisBANFF04,input_banff_04);
  fparm->Update(conv::kNuisBANFF05,input_banff_05);
  fparm->Update(conv::kNuisBANFF06,input_banff_06);
  fparm->Update(conv::kNuisBANFF07,input_banff_07);
  fparm->Update(conv::kNuisBANFF08,input_banff_08);
  fparm->Update(conv::kNuisBANFF09,input_banff_09);
  fparm->Update(conv::kNuisBANFF10,input_banff_10);
  fparm->Update(conv::kNuisBANFF11,input_banff_11);
  fparm->Update(conv::kNuisBANFF12,input_banff_12);
  fparm->Update(conv::kNuisBANFF13,input_banff_13);
  fparm->Update(conv::kNuisBANFF14,input_banff_14);
  fparm->Update(conv::kNuisBANFF15,input_banff_15);
  //banff xsec
  fparm->Update(conv::kNuisBANFFMAQE,input_banff_16);
  fparm->Update(conv::kNuisBANFFMARES,input_banff_17);
  fparm->Update(conv::kNuisBANFFCCQEE0,input_banff_18);
  fparm->Update(conv::kNuisBANFFCCQEE1,input_banff_19);
  fparm->Update(conv::kNuisBANFFCCQEE2,input_banff_20);
  fparm->Update(conv::kNuisBANFFCC1piE0,input_banff_21);
  fparm->Update(conv::kNuisBANFFCC1piE1,input_banff_22);
  //xsec
  fparm->Update(conv::kNuisXSecFermiMom,input_xsec_pf);
  fparm->Update(conv::kNuisXSecSpecFunc,input_xsec_sf);
  fparm->Update(conv::kNuisXSecWshape,input_xsec_wshape);
  fparm->Update(conv::kNuisXSecShapeCCoth,input_xsec_ccoth);
  fparm->Update(conv::kNuisXSecPiLessDelta,input_xsec_piless);
  fparm->Update(conv::kNuisXSecBindE,input_xsec_binde);
  fparm->Update(conv::kNuisXSecNormCCcoh,input_xsec_cccoh);
  fparm->Update(conv::kNuisXSecNormNC1piC,input_xsec_nc1pic);
  fparm->Update(conv::kNuisXSecNormNCoth,input_xsec_ncoth);
  fparm->Update(conv::kNuisXSecNormNuEToNuMu,input_xsec_nuetonumu);
  fparm->Update(conv::kNuisXSecNormAntiNuToNu,input_xsec_antinutonu);
  //fsi
  fparm->Update(conv::kNuisFSI000,input_fsi_0);
  fparm->Update(conv::kNuisFSI001,input_fsi_1);
  fparm->Update(conv::kNuisFSI002,input_fsi_2);
  fparm->Update(conv::kNuisFSI003,input_fsi_3);
  fparm->Update(conv::kNuisFSI004,input_fsi_4);
  fparm->Update(conv::kNuisFSI005,input_fsi_5);

  string oscparams = util::AsString(fparm->PhysicsHypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Plotting spectra for sample: " << util::AsString(s) << endl;
  cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
  cout << "- Fit param list: " << fparm << endl;
  cout << *Opt::Instance() << endl;

  // Histogram booking
  const int nbins = util::NumErecoBins(s);
  const double * ereco_bin_edges = util::ErecoBinEdges(s);
  TH1D * hErecoTot = new TH1D("","",nbins,ereco_bin_edges);
  TH1D * hErecoMode[conv::knm]; //< contribution from each mode  
  const int nbins_coarse = util::NumVariableXsEtrueBins(s,conv::kSKmNuMuCCQE);
  const double * ereco_bin_edges_coarse = util::VariableXsEtrueBinEdges(s,conv::kSKmNuMuCCQE);
  TH1 * hErecoTot_coarse;
  TH1 * hErecoMode_coarse[conv::knm]; //< contribution from each mode

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // Disable BANFF tuning if requested
  if(disable_banff_tuning) {
    cout << "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..." << endl;
    skexp.GetBANFF()->RestoreNominalMC();
  }




  // Generate prediction
  cout << "- Generating predictions..." << endl;
  skexp.Generate(*fparm);

  // Get energy spectra
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) << endl;
    TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
    assert(h);
    hErecoMode[im] = new TH1D(*h);
    hErecoMode_coarse[im] = hErecoMode[im]->Rebin(nbins_coarse,Form("%s_coarse",hErecoMode[im]->GetName()),ereco_bin_edges_coarse);
  }
  hErecoTot = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  hErecoTot_coarse = hErecoTot->Rebin(nbins_coarse,Form("%s_coarse",hErecoTot->GetName()),ereco_bin_edges_coarse);

  //
  // Draw
  //

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TLegend * legend = new TLegend(0.6, 0.4, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
  if(input_mixing23 == 1.0) hErecoTot->GetYaxis()->SetRangeUser(0,5);
  hErecoTot->Draw();
  hErecoTot->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot->GetYaxis()->SetTitle("Number of events");
  THStack * stack = new THStack("stack","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack->Add(hErecoMode[im]);
  }
  stack->Draw("same");
  legend->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend->Draw();
  c->Update();
  c->SaveAs(Form("sk_ereco_spectrum-%f-%f%s.eps",input_mixing23,input_dm23sq,(disable_banff_tuning ? "_nobanff" : "")));
  hErecoTot->GetXaxis()->UnZoom();

  //draw coarse binned histograms
  TCanvas * c_coarse = new TCanvas("c_coarse","",20,20,500,500);
  TLegend * legend_coarse = new TLegend(0.6, 0.4, 0.9, 0.9);
  legend_coarse->SetName("legend_coarse");
  legend_coarse->SetFillColor(0);
  legend_coarse->SetBorderSize(0);
  hErecoTot_coarse->GetXaxis()->SetRangeUser(0,7-conv::kASmallNum);
  hErecoTot_coarse->Draw();
  hErecoTot_coarse->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot_coarse->GetYaxis()->SetTitle("Number of events");
  THStack * stack_coarse = new THStack("stack_coarse","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack_coarse->Add(hErecoMode_coarse[im]);
  }
  stack_coarse->Draw("same");
  legend_coarse->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend_coarse->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend_coarse->Draw();
  c_coarse->Update();
  c_coarse->SaveAs(Form("sk_ereco_spectrum-%f-%f%s_coarse.eps",input_mixing23,input_dm23sq,(disable_banff_tuning ? "_nobanff" : "")));


  //can save this to file &/or print it
  stringstream streamer;

  streamer << fixed;
  double ntot = hErecoTot->Integral();
  double ntot_coarse = hErecoTot->Integral();
  if(TMath::Abs(ntot - ntot_coarse) > conv::kASmallNum) {
    cerr << "Integrals are different between differently binned histograms "
	 << "50MeV: " << ntot
	 << "BANFF binning: " << ntot_coarse << endl;
    exit(1);
  }

  streamer << "Total number of events: " 
       << setw(10) << setfill(' ') << setprecision(6) << ntot << endl;
  const int nslices=6;
  double slice_max[nslices] = { 0.50, 0.75, 1.00, 1.25, 1.50, 2.00 };
  for(int islice=0; islice<nslices; islice++) {
    double esmax = slice_max[islice];
    streamer << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << esmax << " GeV): " 
         << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(esmax)) << endl;
  }
  for(int im = 0; im < conv::knm; im++) {
     conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
     if(m == conv::kSKmUndefined) continue;
     double nmode = hErecoMode[im]->Integral();
     streamer << "Contrib from mode: " 
          << setw(40) << setfill(' ') << util::AsStringTex(m) << " : " 
          << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
          << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  }

  //
  // Save histograms in a ROOT file
  //

  TFile f(Form("sk_ereco_spectrum-%s%s.root",oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : "")) ,"recreate");
  hErecoTot -> Write(util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str());
  hErecoTot_coarse->Write(Form("%s_coarse",util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str()));
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    hErecoMode[im] -> Write(util::AsStringCompact(s,m,"sk_","").c_str());
    hErecoMode_coarse[im] -> Write(Form("sk_%s_coarse",util::AsStringCompact(s,m).c_str()));
  }
  stack->Write();
  legend->Write();
  stack_coarse->Write();
  legend_coarse->Write();
  f.Close();

  //
  // Print-out numbers requested during the analysis review
  //

  if(print_spectrum) {  
    int nprec=5;
    streamer << fixed;

    streamer << endl;
    for(int ibin = 1; ibin <= hErecoTot_coarse->GetXaxis()->GetNbins(); ibin++) {
      streamer << setw(6) << setfill(' ') << ibin
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinLowEdge(ibin)
	       << " - "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinUpEdge(ibin)
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot_coarse->GetBinContent(ibin)
	       << endl;
    }//ibin

    streamer << endl;
    for(int ibin = 1; ibin <= hErecoTot->GetXaxis()->GetNbins(); ibin++) {
      streamer << setw(6) << setfill(' ') << ibin
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinLowEdge(ibin)
	       << " - "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinUpEdge(ibin)
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot->GetBinContent(ibin)
	       << endl;
    }//ibin
  }

  //actually print
  cout << streamer.str() << endl;
  
  // get size of streamer
  long size;
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);

  //write to file
  TString outname = Form("sk_ereco_spectrum-%s%s.txt",oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : ""));
  ofstream outfile (outname);
  outfile.write (streamer.str().c_str(),size);
  outfile.close();

}
