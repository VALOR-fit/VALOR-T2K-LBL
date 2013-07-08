//
// ------------------------------------------------------------------------------------------------
// Use the minimum of a chisq map to seed the fitter, to find the true minimum of the data
// Output is the values and pulls of the fitted parameters, and a deltachisq map
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void fit_xml_seeded_from_chisq_map
(
  const char * xmlfile,                                            /* input XML file with toy MC experiments */
  const int    dataset_id,                                         /* which dataset to use*/
  const char * chisqfile,                                          /* root file where the chisq map is */
  const int    nfits,                                              /* number of fits to do (each with different seeds of the systematic parameters */
  const char * config_file         = "analysis_configuration.xml", /* analysis configuration file */
  const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisanceselection",    /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option          = "minuit",                     /* fitter option (see run_fitter.C) */
  const char * chisq_surface_name  = "",
  const char * mixing_mode         = "2theta",
  const char * hierarchy_mode      = "normal"
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);

  // Global options
  if(mixing_mode == "theta")
    Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  else
    Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  if(hierarchy_mode == "inverted")
    Opt::Instance()->gMassHierarchyAssumed = conv::kNmhInverted;
  else
    Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiment " << dataset_id << " from: " << xmlfile2 << endl;
  cout << "   - Reading chisq map from: " << chisqfile << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit param list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit param list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Open the file to read the chisq map from it
  TFile inputfile(chisqfile,"READ");
  TH2D *hchisq = 0;
  inputfile.GetObject(Form("chisqwithsyst_%d",dataset_id),hchisq);
  if(!hchisq) {
    inputfile.GetObject(chisq_surface_name,hchisq);
  }
  if(!hchisq)
    return;
  TH2D & hdchisq = *hchisq;
  
  // Find the minimum, and convert into mixing23 & dm23sq. This will give the seeds for the osc parameters
  int globalbin = hchisq->GetMinimumBin();
  int xbin, ybin, zbin;
  hchisq->GetBinXYZ(globalbin, xbin, ybin, zbin);
  double input_mixing23 = hchisq->GetXaxis()->GetBinCenter(xbin);
  double input_dm23sq   = hchisq->GetYaxis()->GetBinCenter(ybin);

  // Seed the oscillation parameters
  fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
  fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);

  // random numbers
  TRandom3 rand(0);

  // Get the BANFF errors & bestfits
  double errs_banff_maqe   = skexp->GetBANFF()->OneSigmaError(conv::kNuisBANFFMAQE);
  double errs_banff_ccqeE0 = skexp->GetBANFF()->OneSigmaError(conv::kNuisBANFFCCQEE0);
  double bf_banff_maqe    = skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFMAQE);
  double bf_banff_ccqeE0  = skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCCQEE0);

  // Variables to save bestfits & errors in
  double bestfit_chisq = 9999999;
  double bestfit_chisq5 = 9999999;

  double bestfit_mixing12 = -1, err_mixing12 = -1;
  double bestfit_mixing13 = -1, err_mixing13 = -1;
  double bestfit_mixing23 = -1, err_mixing23 = -1;
  double bestfit_deltacp = -1, err_deltacp = -1;
  double bestfit_dm12sq = -1, err_dm12sq = -1;
  double bestfit_dm23sq = -1, err_dm23sq = -1;
  double bestfit_banff_00 = -1, err_banff_00 = -1;
  double bestfit_banff_01 = -1, err_banff_01 = -1;
  double bestfit_banff_02 = -1, err_banff_02 = -1;
  double bestfit_banff_03 = -1, err_banff_03 = -1;
  double bestfit_banff_04 = -1, err_banff_04 = -1;
  double bestfit_banff_05 = -1, err_banff_05 = -1;
  double bestfit_banff_06 = -1, err_banff_06 = -1;
  double bestfit_banff_07 = -1, err_banff_07 = -1;
  double bestfit_banff_08 = -1, err_banff_08 = -1;
  double bestfit_banff_09 = -1, err_banff_09 = -1;
  double bestfit_banff_10 = -1, err_banff_10 = -1;
  double bestfit_banff_11 = -1, err_banff_11 = -1;
  double bestfit_banff_12 = -1, err_banff_12 = -1;
  double bestfit_banff_13 = -1, err_banff_13 = -1;
  double bestfit_banff_14 = -1, err_banff_14 = -1;
  double bestfit_banff_15 = -1, err_banff_15 = -1;
  double bestfit_banff_maqe = -1, err_banff_maqe = -1;
  double bestfit_banff_mares = -1, err_banff_mares = -1;
  double bestfit_banff_ccqeE0 = -1, err_banff_ccqeE0 = -1;
  double bestfit_banff_ccqeE1 = -1, err_banff_ccqeE1 = -1;
  double bestfit_banff_ccqeE2 = -1, err_banff_ccqeE2 = -1;
  double bestfit_banff_cc1piE0 = -1, err_banff_cc1piE0 = -1;
  double bestfit_banff_cc1piE1 = -1, err_banff_cc1piE1 = -1;
  double bestfit_xsec_pf = -1, err_xsec_pf = -1;
  double bestfit_xsec_sf = -1, err_xsec_sf = -1;
  double bestfit_xsec_binde = -1, err_xsec_binde = -1;
  double bestfit_xsec_wshape = -1, err_xsec_wshape = -1;
  double bestfit_xsec_ccothshape = -1, err_xsec_ccothshape = -1;
  double bestfit_xsec_cc1pishape = -1, err_xsec_cc1pishape = -1;
  double bestfit_xsec_pilessdelta = -1, err_xsec_pilessdelta = -1;
  double bestfit_xsec_normcccoh = -1, err_xsec_normcccoh = -1;
  double bestfit_xsec_normnc1piC = -1, err_xsec_normnc1piC = -1;
  double bestfit_xsec_normncoth = -1, err_xsec_normncoth = -1;
  double bestfit_xsec_normnuetonumu = -1, err_xsec_normnuetonumu = -1;
  double bestfit_xsec_normantinutonu = -1, err_xsec_normantinutonu = -1;
  double bestfit_skdet_energyscale = -1, err_skdet_energyscale = -1;
  double bestfit_skdet_ccqerc0 = -1, err_skdet_ccqerc0 = -1;
  double bestfit_skdet_ccqerc1 = -1, err_skdet_ccqerc1 = -1;
  double bestfit_skdet_ccqerc2 = -1, err_skdet_ccqerc2 = -1;
  double bestfit_skdet_ccqeoth = -1, err_skdet_ccqeoth = -1;
  double bestfit_skdet_ccnqe = -1, err_skdet_ccnqe = -1;
  double bestfit_skdet_nc = -1, err_skdet_nc = -1;
  double bestfit_skdet_nuecc = -1, err_skdet_nuecc = -1;
  double bestfit_fsi_000 = -1, err_fsi_000 = -1;
  double bestfit_fsi_001 = -1, err_fsi_001 = -1;
  double bestfit_fsi_002 = -1, err_fsi_002 = -1;
  double bestfit_fsi_003 = -1, err_fsi_003 = -1;
  double bestfit_fsi_004 = -1, err_fsi_004 = -1;
  double bestfit_fsi_005 = -1, err_fsi_005 = -1;


  double bestfit_mixing23 = -1;
  double bestfit_dm23sq   = -1;
  double bestfit_sf       = -1;
  double bestfit_piless   = -1;
  double bestfit_binde    = -1;
  double bestfit_ccnqeeff = -1;
  double bestfit_nceff    = -1;
  //double bestfit_escale   = -1;
  double bestfit_maqe     = -1;
  double bestfit_ccqeE0   = -1;
  double err_mixing23 = -1;
  double err_dm23sq   = -1;
  double err_sf       = -1;
  double err_piless   = -1;
  double err_binde    = -1;
  double err_ccnqeeff = -1;
  double err_nceff    = -1;
  //double err_escale   = -1;
  double err_maqe     = -1;
  double err_ccqeE0   = -1;

  // Get superk spectrum observation (toy experiment)
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(xmlfile2.Data(),dataset_id);

  //Get spectra
  const int      n_ereco_bins    = util::NumErecoBins (conv::kSKsFCFV1RingMuLike);
  const double * ereco_bin_edges = util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike);
  TH1D * bestfit_spectrum = new TH1D("","", n_ereco_bins, ereco_bin_edges);
  TH1D observation_spectrum = (*(skobs->ErecoSpectrum(conv::kSKsFCFV1RingMuLike)));

  for(int ifit = 0; ifit < nfits; ifit++) {

    // Randomise seeds
    if(fit_param_list == "nuosc3f_numudisapp_fix4+nuisanceselection") {
      fparm->SetInitValue(conv::kNuisXSecSpecFunc,         rand.Uniform(0,1));
      fparm->SetInitValue(conv::kNuisXSecPiLessDelta,      rand.Gaus(0,1));
      fparm->SetInitValue(conv::kNuisXSecBindE,            rand.Gaus(0,1));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCnonQE, rand.Gaus(0,1));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuNC,      rand.Gaus(0,1));
      fparm->SetInitValue(conv::kNuisSuperkEnergyScale,    rand.Gaus(0,1));
      fparm->SetInitValue(conv::kNuisBANFFMAQE,            rand.Gaus(bf_banff_maqe,   errs_banff_maqe));
      fparm->SetInitValue(conv::kNuisBANFFCCQEE0,          rand.Gaus(bf_banff_ccqeE0, errs_banff_ccqeE0));
    }
    else if(fit_param_list == "nuosc3f_numudisapp_fix4+nuisance") {
      skexp->TweakAllSystematics();
      fparm->SetInitValue(conv::kNuisBANFF00, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF00));
      fparm->SetInitValue(conv::kNuisBANFF01, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF01));
      fparm->SetInitValue(conv::kNuisBANFF02, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF02));
      fparm->SetInitValue(conv::kNuisBANFF03, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF03));
      fparm->SetInitValue(conv::kNuisBANFF04, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF04));
      fparm->SetInitValue(conv::kNuisBANFF05, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF05));
      fparm->SetInitValue(conv::kNuisBANFF06, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF06));
      fparm->SetInitValue(conv::kNuisBANFF07, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF07));
      fparm->SetInitValue(conv::kNuisBANFF08, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF08));
      fparm->SetInitValue(conv::kNuisBANFF09, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF09));
      fparm->SetInitValue(conv::kNuisBANFF10, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF10));
      fparm->SetInitValue(conv::kNuisBANFF11, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF11));
      fparm->SetInitValue(conv::kNuisBANFF12, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF12));
      fparm->SetInitValue(conv::kNuisBANFF13, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF13));
      fparm->SetInitValue(conv::kNuisBANFF14, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF14));
      fparm->SetInitValue(conv::kNuisBANFF15, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF15));
      fparm->SetInitValue(conv::kNuisBANFFMAQE, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFMAQE));
      fparm->SetInitValue(conv::kNuisBANFFMARES, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFMARES));
      fparm->SetInitValue(conv::kNuisBANFFCCQEE0, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE0));
      fparm->SetInitValue(conv::kNuisBANFFCCQEE1, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE1));
      fparm->SetInitValue(conv::kNuisBANFFCCQEE2, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE2));
      fparm->SetInitValue(conv::kNuisBANFFCC1piE0, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCC1piE0));
      fparm->SetInitValue(conv::kNuisBANFFCC1piE1, skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCC1piE1));
      fparm->SetInitValue(conv::kNuisXSecFermiMom, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecFermiMom));
      fparm->SetInitValue(conv::kNuisXSecSpecFunc, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecSpecFunc));
      fparm->SetInitValue(conv::kNuisXSecBindE, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecBindE));
      fparm->SetInitValue(conv::kNuisXSecWshape, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecWshape));
      fparm->SetInitValue(conv::kNuisXSecShapeCCoth, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecShapeCCoth));
      fparm->SetInitValue(conv::kNuisXSecPiLessDelta, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecPiLessDelta));
      fparm->SetInitValue(conv::kNuisXSecNormCCcoh, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormCCcoh));
      fparm->SetInitValue(conv::kNuisXSecNormNC1piC, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormNC1piC));
      fparm->SetInitValue(conv::kNuisXSecNormNCoth, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormNCoth));
      fparm->SetInitValue(conv::kNuisXSecNormNuEToNuMu, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormNuEToNuMu));
      fparm->SetInitValue(conv::kNuisXSecNormAntiNuToNu, skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormAntiNuToNu));
      fparm->SetInitValue(conv::kNuisSuperkEnergyScale, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEnergyScale));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCQErc0u, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc0u));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCQErc1u, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc1u));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCQErc2u, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc2u));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCQEoth, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQEoth));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCnonQE, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCnonQE));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuNC, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuNC));
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuNueCC, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuNueCC));
      fparm->SetInitValue(conv::kNuisFSI000, skexp->GetFSI()->CurrentValue(conv::kNuisFSI000));
      fparm->SetInitValue(conv::kNuisFSI001, skexp->GetFSI()->CurrentValue(conv::kNuisFSI001));
      fparm->SetInitValue(conv::kNuisFSI002, skexp->GetFSI()->CurrentValue(conv::kNuisFSI002));
      fparm->SetInitValue(conv::kNuisFSI003, skexp->GetFSI()->CurrentValue(conv::kNuisFSI003));
      fparm->SetInitValue(conv::kNuisFSI004, skexp->GetFSI()->CurrentValue(conv::kNuisFSI004));
      fparm->SetInitValue(conv::kNuisFSI005, skexp->GetFSI()->CurrentValue(conv::kNuisFSI005));
      skexp->ResetSystematicVariations();
    }
    // Set the osc params to be within 3 bins of the chisq surfaces' best fit
    fparm->SetInitValue(conv::kPhypMixing23, rand.Uniform(input_mixing23 - 0.015,  TMath::Min(input_mixing23 + 0.015, 1.0)));
    fparm->SetInitValue(conv::kPhypDm23sq,   rand.Uniform(input_dm23sq   - 0.0003, input_dm23sq - 0.0003));
      
    // Re-seed fitter (reset current values to specified init values)
    fparm->ResetCurrToInitValues();

    // Get fitter
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

    // Do the fit
    fitter->DoTheFit();
    int fitstatus = fitter->FitStatusCode();

    // Get best-fit params (if its a better minimum)
    if(fit_inputs->BestFitChisquare() < bestfit_chisq) {
      bestfit_chisq = fit_inputs->BestFitChisquare();

      fit_inputs->SetChisqDefinition(true);
      bestfit_chisq5 = fit_inputs->BestFitChisquare();
      fit_inputs->SetChisqDefinition(false);

      bestfit_mixing12 = fparm->BestFitValue(conv::kPhypMixing12);
      bestfit_mixing13 = fparm->BestFitValue(conv::kPhypMixing13);
      bestfit_mixing23 = fparm->BestFitValue(conv::kPhypMixing23);
      bestfit_deltacp  = fparm->BestFitValue(conv::kPhypDeltaCP);
      bestfit_dm12sq   = fparm->BestFitValue(conv::kPhypDm12sq);
      bestfit_dm23sq   = fparm->BestFitValue(conv::kPhypDm23sq);

      err_mixing12 = fparm->MigradError(conv::kPhypMixing12);
      err_mixing13 = fparm->MigradError(conv::kPhypMixing13);
      err_mixing23 = fparm->MigradError(conv::kPhypMixing23);
      err_deltacp  = fparm->MigradError(conv::kPhypDeltaCP);
      err_dm12sq   = fparm->MigradError(conv::kPhypDm12sq);
      err_dm23sq   = fparm->MigradError(conv::kPhypDm23sq);

      if(fparm->Exists(conv::kNuisBANFF00)) bestfit_banff_00 = fparm->BestFitValue(conv::kNuisBANFF00);
      if(fparm->Exists(conv::kNuisBANFF01)) bestfit_banff_01 = fparm->BestFitValue(conv::kNuisBANFF01);
      if(fparm->Exists(conv::kNuisBANFF02)) bestfit_banff_02 = fparm->BestFitValue(conv::kNuisBANFF02);
      if(fparm->Exists(conv::kNuisBANFF03)) bestfit_banff_03 = fparm->BestFitValue(conv::kNuisBANFF03);
      if(fparm->Exists(conv::kNuisBANFF04)) bestfit_banff_04 = fparm->BestFitValue(conv::kNuisBANFF04);
      if(fparm->Exists(conv::kNuisBANFF05)) bestfit_banff_05 = fparm->BestFitValue(conv::kNuisBANFF05);
      if(fparm->Exists(conv::kNuisBANFF06)) bestfit_banff_06 = fparm->BestFitValue(conv::kNuisBANFF06);
      if(fparm->Exists(conv::kNuisBANFF07)) bestfit_banff_07 = fparm->BestFitValue(conv::kNuisBANFF07);
      if(fparm->Exists(conv::kNuisBANFF08)) bestfit_banff_08 = fparm->BestFitValue(conv::kNuisBANFF08);
      if(fparm->Exists(conv::kNuisBANFF09)) bestfit_banff_09 = fparm->BestFitValue(conv::kNuisBANFF09);
      if(fparm->Exists(conv::kNuisBANFF10)) bestfit_banff_10 = fparm->BestFitValue(conv::kNuisBANFF10);
      if(fparm->Exists(conv::kNuisBANFF11)) bestfit_banff_11 = fparm->BestFitValue(conv::kNuisBANFF11);
      if(fparm->Exists(conv::kNuisBANFF12)) bestfit_banff_12 = fparm->BestFitValue(conv::kNuisBANFF12);
      if(fparm->Exists(conv::kNuisBANFF13)) bestfit_banff_13 = fparm->BestFitValue(conv::kNuisBANFF13);
      if(fparm->Exists(conv::kNuisBANFF14)) bestfit_banff_14 = fparm->BestFitValue(conv::kNuisBANFF14);
      if(fparm->Exists(conv::kNuisBANFF15)) bestfit_banff_15 = fparm->BestFitValue(conv::kNuisBANFF15);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) bestfit_banff_maqe = fparm->BestFitValue(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) bestfit_banff_mares = fparm->BestFitValue(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) bestfit_banff_ccqeE0 = fparm->BestFitValue(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) bestfit_banff_ccqeE1 = fparm->BestFitValue(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) bestfit_banff_ccqeE2 = fparm->BestFitValue(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) bestfit_banff_cc1piE0 = fparm->BestFitValue(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) bestfit_banff_cc1piE1 = fparm->BestFitValue(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) bestfit_xsec_pf = fparm->BestFitValue(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) bestfit_xsec_sf = fparm->BestFitValue(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecBindE)) bestfit_xsec_binde = fparm->BestFitValue(conv::kNuisXSecBindE);
      if(fparm->Exists(conv::kNuisXSecWshape)) bestfit_xsec_wshape = fparm->BestFitValue(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) bestfit_xsec_ccothshape = fparm->BestFitValue(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) bestfit_xsec_pilessdelta = fparm->BestFitValue(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) bestfit_xsec_normcccoh = fparm->BestFitValue(conv::kNuisXSecNormCCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) bestfit_xsec_normnc1piC = fparm->BestFitValue(conv::kNuisXSecNormNC1piC);
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) bestfit_xsec_normncoth = fparm->BestFitValue(conv::kNuisXSecNormNCoth);
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) bestfit_xsec_normnuetonumu = fparm->BestFitValue(conv::kNuisXSecNormNuEToNuMu);
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) bestfit_xsec_normantinutonu = fparm->BestFitValue(conv::kNuisXSecNormAntiNuToNu);
      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) bestfit_skdet_energyscale = fparm->BestFitValue(conv::kNuisSuperkEnergyScale);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) bestfit_skdet_ccqerc0 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc0u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) bestfit_skdet_ccqerc1 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc1u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) bestfit_skdet_ccqerc2 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc2u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) bestfit_skdet_ccqeoth = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQEoth);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) bestfit_skdet_ccnqe = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) bestfit_skdet_nc = fparm->BestFitValue(conv::kNuisSuperkEff1RmuNC);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) bestfit_skdet_nuecc = fparm->BestFitValue(conv::kNuisSuperkEff1RmuNueCC);
      if(fparm->Exists(conv::kNuisFSI000)) bestfit_fsi_000 = fparm->BestFitValue(conv::kNuisFSI000);
      if(fparm->Exists(conv::kNuisFSI001)) bestfit_fsi_001 = fparm->BestFitValue(conv::kNuisFSI001);
      if(fparm->Exists(conv::kNuisFSI002)) bestfit_fsi_002 = fparm->BestFitValue(conv::kNuisFSI002);
      if(fparm->Exists(conv::kNuisFSI003)) bestfit_fsi_003 = fparm->BestFitValue(conv::kNuisFSI003);
      if(fparm->Exists(conv::kNuisFSI004)) bestfit_fsi_004 = fparm->BestFitValue(conv::kNuisFSI004);
      if(fparm->Exists(conv::kNuisFSI005)) bestfit_fsi_005 = fparm->BestFitValue(conv::kNuisFSI005);

      if(fparm->Exists(conv::kNuisBANFF00)) err_banff_00 = fparm->MigradError(conv::kNuisBANFF00);
      if(fparm->Exists(conv::kNuisBANFF01)) err_banff_01 = fparm->MigradError(conv::kNuisBANFF01);
      if(fparm->Exists(conv::kNuisBANFF02)) err_banff_02 = fparm->MigradError(conv::kNuisBANFF02);
      if(fparm->Exists(conv::kNuisBANFF03)) err_banff_03 = fparm->MigradError(conv::kNuisBANFF03);
      if(fparm->Exists(conv::kNuisBANFF04)) err_banff_04 = fparm->MigradError(conv::kNuisBANFF04);
      if(fparm->Exists(conv::kNuisBANFF05)) err_banff_05 = fparm->MigradError(conv::kNuisBANFF05);
      if(fparm->Exists(conv::kNuisBANFF06)) err_banff_06 = fparm->MigradError(conv::kNuisBANFF06);
      if(fparm->Exists(conv::kNuisBANFF07)) err_banff_07 = fparm->MigradError(conv::kNuisBANFF07);
      if(fparm->Exists(conv::kNuisBANFF08)) err_banff_08 = fparm->MigradError(conv::kNuisBANFF08);
      if(fparm->Exists(conv::kNuisBANFF09)) err_banff_09 = fparm->MigradError(conv::kNuisBANFF09);
      if(fparm->Exists(conv::kNuisBANFF10)) err_banff_10 = fparm->MigradError(conv::kNuisBANFF10);
      if(fparm->Exists(conv::kNuisBANFF11)) err_banff_11 = fparm->MigradError(conv::kNuisBANFF11);
      if(fparm->Exists(conv::kNuisBANFF12)) err_banff_12 = fparm->MigradError(conv::kNuisBANFF12);
      if(fparm->Exists(conv::kNuisBANFF13)) err_banff_13 = fparm->MigradError(conv::kNuisBANFF13);
      if(fparm->Exists(conv::kNuisBANFF14)) err_banff_14 = fparm->MigradError(conv::kNuisBANFF14);
      if(fparm->Exists(conv::kNuisBANFF15)) err_banff_15 = fparm->MigradError(conv::kNuisBANFF15);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) err_banff_maqe = fparm->MigradError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) err_banff_mares = fparm->MigradError(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) err_banff_ccqeE0 = fparm->MigradError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) err_banff_ccqeE1 = fparm->MigradError(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) err_banff_ccqeE2 = fparm->MigradError(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) err_banff_cc1piE0 = fparm->MigradError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) err_banff_cc1piE1 = fparm->MigradError(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) err_xsec_pf = fparm->MigradError(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) err_xsec_sf = fparm->MigradError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecBindE)) err_xsec_binde = fparm->MigradError(conv::kNuisXSecBindE);
      if(fparm->Exists(conv::kNuisXSecWshape)) err_xsec_wshape = fparm->MigradError(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) err_xsec_ccothshape = fparm->MigradError(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) err_xsec_pilessdelta = fparm->MigradError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) err_xsec_normcccoh = fparm->MigradError(conv::kNuisXSecNormCCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) err_xsec_normnc1piC = fparm->MigradError(conv::kNuisXSecNormNC1piC);
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) err_xsec_normncoth = fparm->MigradError(conv::kNuisXSecNormNCoth);
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) err_xsec_normnuetonumu = fparm->MigradError(conv::kNuisXSecNormNuEToNuMu);
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) err_xsec_normantinutonu = fparm->MigradError(conv::kNuisXSecNormAntiNuToNu);
      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) err_skdet_energyscale = fparm->MigradError(conv::kNuisSuperkEnergyScale);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) err_skdet_ccqerc0 = fparm->MigradError(conv::kNuisSuperkEff1RmuCCQErc0u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) err_skdet_ccqerc1 = fparm->MigradError(conv::kNuisSuperkEff1RmuCCQErc1u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) err_skdet_ccqerc2 = fparm->MigradError(conv::kNuisSuperkEff1RmuCCQErc2u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) err_skdet_ccqeoth = fparm->MigradError(conv::kNuisSuperkEff1RmuCCQEoth);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) err_skdet_ccnqe = fparm->MigradError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) err_skdet_nc = fparm->MigradError(conv::kNuisSuperkEff1RmuNC);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) err_skdet_nuecc = fparm->MigradError(conv::kNuisSuperkEff1RmuNueCC);
      if(fparm->Exists(conv::kNuisFSI000)) err_fsi_000 = fparm->MigradError(conv::kNuisFSI000);
      if(fparm->Exists(conv::kNuisFSI001)) err_fsi_001 = fparm->MigradError(conv::kNuisFSI001);
      if(fparm->Exists(conv::kNuisFSI002)) err_fsi_002 = fparm->MigradError(conv::kNuisFSI002);
      if(fparm->Exists(conv::kNuisFSI003)) err_fsi_003 = fparm->MigradError(conv::kNuisFSI003);
      if(fparm->Exists(conv::kNuisFSI004)) err_fsi_004 = fparm->MigradError(conv::kNuisFSI004);
      if(fparm->Exists(conv::kNuisFSI005)) err_fsi_005 = fparm->MigradError(conv::kNuisFSI005);

      // Print-out best-fit params
      fparm->ShowOnlyBestFitValuesInNextPrint();
      cout << "- Best-fit oscillation params (below) yield chisq = " << bestfit_chisq 
	   << *fparm;
      //Get the bestfit spectrum
      for(int ibin = 1; ibin <= n_ereco_bins; ibin++)
	bestfit_spectrum->SetBinContent(ibin, skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingMuLike)[ibin-1]);

      stringstream streamer;
      streamer << "\\chi^2_{bf} \t" << bestfit_chisq << endl;
      streamer << "\\chi^2_{gof} \t" << bestfit_chisq5 << endl << endl;

      streamer << util::AsStringPaperRooTex(conv::kPhypMixing12) << "\t" << bestfit_mixing12 << "\t" << err_mixing12 << endl;
      streamer << util::AsStringPaperRooTex(conv::kPhypMixing13) << "\t" << bestfit_mixing13 << "\t" << err_mixing13 << endl;
      streamer << util::AsStringPaperRooTex(conv::kPhypMixing23) << "\t" << bestfit_mixing23 << "\t" << err_mixing23 << endl;
      streamer << util::AsStringPaperRooTex(conv::kPhypDeltaCP) << "\t" << bestfit_deltacp  << "\t" << err_deltacp  << endl;
      streamer << util::AsStringPaperRooTex(conv::kPhypDm12sq) << "\t" << bestfit_dm12sq   << "\t" << err_dm12sq   << endl;
      streamer << util::AsStringPaperRooTex(conv::kPhypDm23sq) << "\t" << bestfit_dm23sq   << "\t" << err_dm23sq   << endl;

      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) streamer << util::AsStringPaperRooTex(conv::kNuisSuperkEnergyScale) << "\t" << bestfit_skdet_energyscale << "\t" << err_skdet_energyscale << endl;
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) streamer << util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc0u) << "\t" << bestfit_skdet_ccqerc0 << "\t" << err_skdet_ccqerc0 << endl;
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) streamer << util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc1u) << "\t" << bestfit_skdet_ccqerc1 << "\t" << err_skdet_ccqerc1 << endl;
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) streamer << util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc2u) << "\t" << bestfit_skdet_ccqerc2 << "\t" << err_skdet_ccqerc2 << endl;
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) streamer << util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQEoth) << "\t" << bestfit_skdet_ccqeoth << "\t" << err_skdet_ccqeoth << endl;
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) streamer << util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCnonQE) << "\t" << bestfit_skdet_ccnqe << "\t" << err_skdet_ccnqe << endl;
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) streamer << util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuNC) << "\t" << bestfit_skdet_nc << "\t" << err_skdet_nc << endl;
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) streamer << util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuNueCC) << "\t" << bestfit_skdet_nuecc << "\t" << err_skdet_nuecc << endl;
      if(fparm->Exists(conv::kNuisBANFF00)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF00) << "\t" << bestfit_banff_00 << "\t" << err_banff_00 << endl;
      if(fparm->Exists(conv::kNuisBANFF01)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF01) << "\t" << bestfit_banff_01 << "\t" << err_banff_01 << endl;
      if(fparm->Exists(conv::kNuisBANFF02)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF02) << "\t" << bestfit_banff_02 << "\t" << err_banff_02 << endl;
      if(fparm->Exists(conv::kNuisBANFF03)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF03) << "\t" << bestfit_banff_03 << "\t" << err_banff_03 << endl;
      if(fparm->Exists(conv::kNuisBANFF04)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF04) << "\t" << bestfit_banff_04 << "\t" << err_banff_04 << endl;
      if(fparm->Exists(conv::kNuisBANFF05)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF05) << "\t" << bestfit_banff_05 << "\t" << err_banff_05 << endl;
      if(fparm->Exists(conv::kNuisBANFF06)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF06) << "\t" << bestfit_banff_06 << "\t" << err_banff_06 << endl;
      if(fparm->Exists(conv::kNuisBANFF07)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF07) << "\t" << bestfit_banff_07 << "\t" << err_banff_07 << endl;
      if(fparm->Exists(conv::kNuisBANFF08)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF08) << "\t" << bestfit_banff_08 << "\t" << err_banff_08 << endl;
      if(fparm->Exists(conv::kNuisBANFF09)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF09) << "\t" << bestfit_banff_09 << "\t" << err_banff_09 << endl;
      if(fparm->Exists(conv::kNuisBANFF10)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF10) << "\t" << bestfit_banff_10 << "\t" << err_banff_10 << endl;
      if(fparm->Exists(conv::kNuisBANFF11)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF11) << "\t" << bestfit_banff_11 << "\t" << err_banff_11 << endl;
      if(fparm->Exists(conv::kNuisBANFF12)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF12) << "\t" << bestfit_banff_12 << "\t" << err_banff_12 << endl;
      if(fparm->Exists(conv::kNuisBANFF13)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF13) << "\t" << bestfit_banff_13 << "\t" << err_banff_13 << endl;
      if(fparm->Exists(conv::kNuisBANFF14)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF14) << "\t" << bestfit_banff_14 << "\t" << err_banff_14 << endl;
      if(fparm->Exists(conv::kNuisBANFF15)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFF15) << "\t" << bestfit_banff_15 << "\t" << err_banff_15 << endl;
      if(fparm->Exists(conv::kNuisBANFFMAQE)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFFMAQE) << "\t" << bestfit_banff_maqe << "\t" << err_banff_maqe << endl;
      if(fparm->Exists(conv::kNuisBANFFMARES)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFFMARES) << "\t" << bestfit_banff_mares << "\t" << err_banff_mares << endl;
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE0) << "\t" << bestfit_banff_ccqeE0 << "\t" << err_banff_ccqeE0 << endl;
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE1) << "\t" << bestfit_banff_ccqeE1 << "\t" << err_banff_ccqeE1 << endl;
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE2) << "\t" << bestfit_banff_ccqeE2 << "\t" << err_banff_ccqeE2 << endl;
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFFCC1piE0) << "\t" << bestfit_banff_cc1piE0 << "\t" << err_banff_cc1piE0 << endl;
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) streamer << util::AsStringPaperRooTex(conv::kNuisBANFFCC1piE1) << "\t" << bestfit_banff_cc1piE1 << "\t" << err_banff_cc1piE1 << endl;
      if(fparm->Exists(conv::kNuisXSecFermiMom)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecFermiMom) << "\t" << bestfit_xsec_pf << "\t" << err_xsec_pf << endl;
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecSpecFunc) << "\t" << bestfit_xsec_sf << "\t" << err_xsec_sf << endl;
      if(fparm->Exists(conv::kNuisXSecBindE)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecBindE) << "\t" << bestfit_xsec_binde << "\t" << err_xsec_binde << endl;
      if(fparm->Exists(conv::kNuisXSecWshape)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecWshape) << "\t" << bestfit_xsec_wshape << "\t" << err_xsec_wshape << endl;
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecShapeCCoth) << "\t" << bestfit_xsec_ccothshape << "\t" << err_xsec_ccothshape << endl;
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecPiLessDelta) << "\t" << bestfit_xsec_pilessdelta << "\t" << err_xsec_pilessdelta << endl;
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecNormCCcoh) << "\t" << bestfit_xsec_normcccoh << "\t" << err_xsec_normcccoh << endl;
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecNormNC1piC) << "\t" << bestfit_xsec_normnc1piC << "\t" << err_xsec_normnc1piC << endl;
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecNormNCoth) << "\t" << bestfit_xsec_normncoth << "\t" << err_xsec_normncoth << endl;
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecNormNuEToNuMu) << "\t" << bestfit_xsec_normnuetonumu << "\t" << err_xsec_normnuetonumu << endl;
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) streamer << util::AsStringPaperRooTex(conv::kNuisXSecNormAntiNuToNu) << "\t" << bestfit_xsec_normantinutonu << "\t" << err_xsec_normantinutonu << endl;
      if(fparm->Exists(conv::kNuisFSI000)) streamer << util::AsStringPaperRooTex(conv::kNuisFSI000) << "\t" << bestfit_fsi_000 << "\t" << err_fsi_000 << endl;
      if(fparm->Exists(conv::kNuisFSI001)) streamer << util::AsStringPaperRooTex(conv::kNuisFSI001) << "\t" << bestfit_fsi_001 << "\t" << err_fsi_001 << endl;
      if(fparm->Exists(conv::kNuisFSI002)) streamer << util::AsStringPaperRooTex(conv::kNuisFSI002) << "\t" << bestfit_fsi_002 << "\t" << err_fsi_002 << endl;
      if(fparm->Exists(conv::kNuisFSI003)) streamer << util::AsStringPaperRooTex(conv::kNuisFSI003) << "\t" << bestfit_fsi_003 << "\t" << err_fsi_003 << endl;
      if(fparm->Exists(conv::kNuisFSI004)) streamer << util::AsStringPaperRooTex(conv::kNuisFSI004) << "\t" << bestfit_fsi_004 << "\t" << err_fsi_004 << endl;
      if(fparm->Exists(conv::kNuisFSI005)) streamer << util::AsStringPaperRooTex(conv::kNuisFSI005) << "\t" << bestfit_fsi_005 << "\t" << err_fsi_005 << endl;

      // write the table
      long size;
      // get size of file
      streamer.seekg(0,ifstream::end);
      size=streamer.tellg();
      streamer.seekg(0);
      ofstream outfile (Form("%s.data",chisqfile),ios_base::out);
      // write to outfile
      outfile.write (streamer.str().c_str(),size);
      outfile.close();
    }
    delete fit_inputs;
    delete fitter;
  }//ifit
  delete skobs;

  //Output file
  //inputfile.Close();
  TFile outputfile(Form("fitresults_%d.root",dataset_id),"RECREATE");

  // Now we have the best fits, lets plot them
  for(int ix = 1; ix <= hdchisq.GetNbinsX(); ix++) {
    for(int iy = 1; iy <= hdchisq.GetNbinsY(); iy++) {
      hdchisq.SetBinContent(ix, iy, hdchisq.GetBinContent(ix, iy) - bestfit_chisq);
    }//iy
  }//ix

  gStyle->SetOptStat(0);

  TCanvas *c0 = new TCanvas();
  hdchisq.Draw("colz");
  c0->SaveAs(Form("dchisq_%d.pdf",dataset_id));
  double constchisq[2] = {2.30, 4.61};
  hdchisq.SetContour(2, constchisq);
  hdchisq.SetDirectory(0);
  hdchisq.Write("hdchisq");
  TCanvas *c1 = new TCanvas();
  hdchisq.SetTitle("");
  hdchisq.Draw("cont1");
  c1->SaveAs(Form("contour_%d.pdf",dataset_id));

  // Also plot the pulls & fitted values of the systematics
  TH1D * hsystpull = new TH1D("hsystpull","hsystpull",fparm->Size()-6,0,fparm->Size()-6);
  TH1D * hsystfit = new TH1D("hsystfit","hsystfit",fparm->Size()-6,0,fparm->Size()-6);
  if(fit_param_list == "nuosc3f_numudisapp_fix4+nuisanceselection") {
    hsystpull->GetXaxis()->SetBinLabel(1, "SF");
    hsystpull->GetXaxis()->SetBinLabel(2, "Piless-Delta");
    hsystpull->GetXaxis()->SetBinLabel(3, "Binding-E");
    hsystpull->GetXaxis()->SetBinLabel(4, "CCnQE-eff");
    hsystpull->GetXaxis()->SetBinLabel(5, "NC-eff");
    hsystpull->GetXaxis()->SetBinLabel(6, "MAQE");
    hsystpull->GetXaxis()->SetBinLabel(7, "CCQE-E0");
    hsystpull->GetXaxis()->SetBinLabel(8,"E-scale");
    hsystpull->SetBinContent(1, bestfit_xsec_sf          / err_xsec_sf);
    hsystpull->SetBinContent(2, bestfit_xsec_pilessdelta      / err_xsec_pilessdelta);
    hsystpull->SetBinContent(3, bestfit_xsec_binde       / err_xsec_binde);
    hsystpull->SetBinContent(4, bestfit_skdet_ccnqe    / err_skdet_ccnqe);
    hsystpull->SetBinContent(5, bestfit_skdet_nc       / err_skdet_nc);
    hsystpull->SetBinContent(6, (bestfit_banff_maqe   - bf_banff_maqe  )  / err_banff_maqe);
    hsystpull->SetBinContent(7, (bestfit_banff_ccqeE0 - bf_banff_ccqeE0)  / err_banff_ccqeE0);
    hsystpull->SetBinContent(8 ,bestfit_skdet_energyscale      / err_skdet_energyscale);
    hsystfit->GetXaxis()->SetBinLabel(1, "SF");
    hsystfit->GetXaxis()->SetBinLabel(2, "Piless-Delta");
    hsystfit->GetXaxis()->SetBinLabel(3, "Binding-E");
    hsystfit->GetXaxis()->SetBinLabel(4, "CCnQE-eff");
    hsystfit->GetXaxis()->SetBinLabel(5, "NC-eff");
    hsystfit->GetXaxis()->SetBinLabel(6, "MAQE");
    hsystfit->GetXaxis()->SetBinLabel(7, "CCQE-E0");
    hsystfit->GetXaxis()->SetBinLabel(8,"E-scale");
    hsystfit->SetBinContent(1, bestfit_xsec_sf          / err_xsec_sf);
    hsystfit->SetBinContent(2, bestfit_xsec_pilessdelta      / err_xsec_pilessdelta);
    hsystfit->SetBinContent(3, bestfit_xsec_binde       / err_xsec_binde);
    hsystfit->SetBinContent(4, bestfit_skdet_ccnqe    / err_skdet_ccnqe);
    hsystfit->SetBinContent(5, bestfit_skdet_nc       / err_skdet_nc);
    hsystfit->SetBinContent(6, (bestfit_banff_maqe   - bf_banff_maqe  )  / err_banff_maqe);
    hsystfit->SetBinContent(7, (bestfit_banff_ccqeE0 - bf_banff_ccqeE0)  / err_banff_ccqeE0);
    hsystfit->SetBinContent(8 ,bestfit_skdet_energyscale      / err_skdet_energyscale);
  }
  else if(fit_param_list == "nuosc3f_numudisapp_fix4+nuisance") {
    int ip = 1;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEnergyScale).c_str());
    hsystpull->SetBinContent(ip, (bestfit_skdet_energyscale - 0) / err_skdet_energyscale);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQEoth).c_str());
    hsystpull->SetBinContent(ip, (bestfit_skdet_ccqeoth - 0) / err_skdet_ccqeoth);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc0u).c_str());
    hsystpull->SetBinContent(ip, (bestfit_skdet_ccqerc0 - 0) / err_skdet_ccqerc0);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc1u).c_str());
    hsystpull->SetBinContent(ip, (bestfit_skdet_ccqerc1 - 0) / err_skdet_ccqerc1);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc2u).c_str());
    hsystpull->SetBinContent(ip, (bestfit_skdet_ccqerc2 - 0) / err_skdet_ccqerc2);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCnonQE).c_str());
    hsystpull->SetBinContent(ip, (bestfit_skdet_ccnqe - 0) / err_skdet_ccnqe);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuNC).c_str());
    hsystpull->SetBinContent(ip, (bestfit_skdet_nc - 0) / err_skdet_nc);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuNueCC).c_str());
    hsystpull->SetBinContent(ip, (bestfit_skdet_nuecc - 0) / err_skdet_nuecc);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF00).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_00 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF00)) / err_banff_00);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF01).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_01 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF01)) / err_banff_01);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF02).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_02 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF02)) / err_banff_02);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF03).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_03 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF03)) / err_banff_03);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF04).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_04 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF04)) / err_banff_04);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF05).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_05 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF05)) / err_banff_05);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF06).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_06 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF06)) / err_banff_06);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF07).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_07 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF07)) / err_banff_07);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF08).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_08 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF08)) / err_banff_08);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF09).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_09 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF09)) / err_banff_09);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF10).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_10 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF10)) / err_banff_10);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF11).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_11 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF11)) / err_banff_11);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF12).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_12 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF12)) / err_banff_12);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF13).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_13 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF13)) / err_banff_13);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF14).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_14 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF14)) / err_banff_14);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF15).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_15 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFF15)) / err_banff_15);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFMAQE).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_maqe - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFMAQE)) / err_banff_maqe);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFMARES).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_mares - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFMARES)) / err_banff_mares);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE0).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_ccqeE0 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCCQEE0)) / err_banff_ccqeE0);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE1).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_ccqeE1 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCCQEE1)) / err_banff_ccqeE1);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE2).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_ccqeE2 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCCQEE2)) / err_banff_ccqeE2);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCC1piE0).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_cc1piE0 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCC1piE0)) / err_banff_cc1piE0);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCC1piE1).c_str());
    hsystpull->SetBinContent(ip, (bestfit_banff_cc1piE1 - skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCC1piE1)) / err_banff_cc1piE1);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecFermiMom).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_pf - 0) / err_xsec_pf);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecSpecFunc).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_sf - 0) / err_xsec_sf);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecBindE).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_binde - 0) / err_xsec_binde);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecWshape).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_wshape - 0) / err_xsec_wshape);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecShapeCCoth).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_ccothshape - 0) / err_xsec_ccothshape);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecPiLessDelta).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_pilessdelta - 0) / err_xsec_pilessdelta);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormCCcoh).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_normcccoh - 0) / err_xsec_normcccoh);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormNC1piC).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_normnc1piC - 0) / err_xsec_normnc1piC);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormNCoth).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_normncoth - 0) / err_xsec_normncoth);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormNuEToNuMu).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_normnuetonumu - 0) / err_xsec_normnuetonumu);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormAntiNuToNu).c_str());
    hsystpull->SetBinContent(ip, (bestfit_xsec_normantinutonu - 0) / err_xsec_normantinutonu);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI000).c_str());
    hsystpull->SetBinContent(ip, (bestfit_fsi_000 - 1) / err_fsi_000);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI001).c_str());
    hsystpull->SetBinContent(ip, (bestfit_fsi_001 - 1) / err_fsi_001);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI002).c_str());
    hsystpull->SetBinContent(ip, (bestfit_fsi_002 - 1) / err_fsi_002);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI003).c_str());
    hsystpull->SetBinContent(ip, (bestfit_fsi_003 - 1) / err_fsi_003);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI004).c_str());
    hsystpull->SetBinContent(ip, (bestfit_fsi_004 - 1) / err_fsi_004);
    ip++;
    hsystpull->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI005).c_str());
    hsystpull->SetBinContent(ip, (bestfit_fsi_005 - 1) / err_fsi_005);
    ip++;
    //best fits
    ip = 1;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEnergyScale).c_str());
    hsystfit->SetBinContent(ip, (bestfit_skdet_energyscale));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQEoth).c_str());
    hsystfit->SetBinContent(ip, (bestfit_skdet_ccqeoth));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc0u).c_str());
    hsystfit->SetBinContent(ip, (bestfit_skdet_ccqerc0));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc1u).c_str());
    hsystfit->SetBinContent(ip, (bestfit_skdet_ccqerc1));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCQErc2u).c_str());
    hsystfit->SetBinContent(ip, (bestfit_skdet_ccqerc2));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuCCnonQE).c_str());
    hsystfit->SetBinContent(ip, (bestfit_skdet_ccnqe));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuNC).c_str());
    hsystfit->SetBinContent(ip, (bestfit_skdet_nc));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisSuperkEff1RmuNueCC).c_str());
    hsystfit->SetBinContent(ip, (bestfit_skdet_nuecc));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF00).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_00));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF01).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_01));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF02).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_02));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF03).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_03));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF04).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_04));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF05).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_05));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF06).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_06));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF07).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_07));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF08).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_08));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF09).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_09));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF10).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_10));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF11).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_11));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF12).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_12));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF13).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_13));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF14).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_14));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFF15).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_15));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFMAQE).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_maqe));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFMARES).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_mares));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE0).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_ccqeE0));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE1).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_ccqeE1));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCCQEE2).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_ccqeE2));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCC1piE0).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_cc1piE0));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisBANFFCC1piE1).c_str());
    hsystfit->SetBinContent(ip, (bestfit_banff_cc1piE1));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecFermiMom).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_pf));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecSpecFunc).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_sf));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecBindE).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_binde));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecWshape).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_wshape));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecShapeCCoth).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_ccothshape));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecPiLessDelta).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_pilessdelta));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormCCcoh).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_normcccoh));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormNC1piC).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_normnc1piC));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormNCoth).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_normncoth));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormNuEToNuMu).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_normnuetonumu));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisXSecNormAntiNuToNu).c_str());
    hsystfit->SetBinContent(ip, (bestfit_xsec_normantinutonu));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI000).c_str());
    hsystfit->SetBinContent(ip, (bestfit_fsi_000));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI001).c_str());
    hsystfit->SetBinContent(ip, (bestfit_fsi_001));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI002).c_str());
    hsystfit->SetBinContent(ip, (bestfit_fsi_002));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI003).c_str());
    hsystfit->SetBinContent(ip, (bestfit_fsi_003));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI004).c_str());
    hsystfit->SetBinContent(ip, (bestfit_fsi_004));
    ip++;
    hsystfit->GetXaxis()->SetBinLabel(ip, util::AsStringPaperRooTex(conv::kNuisFSI005).c_str());
    hsystfit->SetBinContent(ip, (bestfit_fsi_005));
    ip++;
  }
  TCanvas *c2 = new TCanvas();
  hsystpull->Draw();
  c2->SaveAs(Form("pulls_%d.pdf",dataset_id));
  hsystpull->Write();
  hsystfit->Draw();
  c2->SaveAs(Form("fits_%d.pdf",dataset_id));
  hsystfit->Write();
  
  //And plot the best fit spectra with the observation
  TCanvas *c3 = new TCanvas();
  bestfit_spectrum->GetXaxis()->SetRangeUser(0,9.99);
  observation_spectrum.GetXaxis()->SetRangeUser(0,9.99);
  observation_spectrum.Draw("E1");
  bestfit_spectrum->Draw("same");
  c3->SaveAs(Form("spectrum_%d.pdf",dataset_id));
  bestfit_spectrum->Write("bf_spectrum");
  observation_spectrum.Write("obs_spectrum");

  //And not forgetting the value & error of the mixing angles
  TVectorD v(6);
  v[0] = bestfit_mixing23;
  v[1] = err_mixing23;
  v[2] = bestfit_dm23sq;
  v[3] = err_dm23sq;
  v[4] = bestfit_chisq;
  v[5] = bestfit_chisq5;
  cout << endl << "Bestfit of oscillation parameters yield chisq: " << bestfit_chisq << endl
       << "\tMixing23: " << bestfit_mixing23 << " +- " << err_mixing23 << endl
       << "\tDm32sq:   "   << bestfit_dm23sq << " +- "   << err_dm23sq << endl
       << "\t5 bin chisq: " << bestfit_chisq5 << endl;
  v.Write("oscparams");

  //inputfile.Close();
  outputfile.Close();
}
