//
// ------------------------------------------------------------------------------------------------
// Print chi^2 for the given dataset using the prediction with the input oscillation params
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
//#include "Conventions.h"

void print_dataset_chisq
(
  double       input_mixing12      = -99,                         /* */
  double       input_mixing13      = -99,                         /* */
  double       input_mixing23      = -99,                         /* */
  double       input_deltacp       = -99,                         /* */
  double       input_dm12sq        = -99,                         /* */
  double       input_dm23sq        = -99,                         /* */
  const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisance",    /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  int          dataset_id          = 3,                           /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file           = "inputs/data.xml",           /* */
  const char * config_file         = "analysis_configuration.xml",/* analysis configuration file */
  const int    nfits               = 10                           /* number of fits to do (if systematics included in FitParamList, throwing syst to get smaller chisq */
)
{
  using namespace std;
  using namespace t2kana;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  cout << *Opt::Instance() << endl;

  // Get superk spectrum predictor & configure
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Set osc params
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
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "Specified physics hypothesis: " << util::AsString(hyp) << endl;
  if(hyp == conv::kHypOsc3F) {
    fparm->SetInitValue(conv::kPhypMixing12, input_mixing12);
    fparm->SetInitValue(conv::kPhypMixing13, input_mixing13);
    fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
    fparm->SetInitValue(conv::kPhypDeltaCP,  input_deltacp);
    fparm->SetInitValue(conv::kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
  }
  else
  if(hyp == conv::kHypNuMuDisapp2F) {
    fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
    fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
  }
  else {
    return;
  }
  //fix the osc params
  fparm->FixParameter(conv::kPhypMixing23);
  fparm->FixParameter(conv::kPhypDm23sq);

  // Generate prediction
  skexp->Generate(PhysHypothesis(*fparm));
//TH1D * hExp = new TH1D(*skexp->CurrPredictedErecoSpectrum());

  // Get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);
//TH1D * hObs = new TH1D(*skobs->ErecoSpectrum());

  double chisq = 9999999;

  if(fparm->GetSystMode() > 0) {

    TRandom3 rand(0);

    double errs_banff_maqe   = skexp->GetBANFF()->OneSigmaError(conv::kNuisBANFFMAQE);
    double errs_banff_ccqeE0 = skexp->GetBANFF()->OneSigmaError(conv::kNuisBANFFCCQEE0);
    double bf_banff_maqe    = skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFMAQE);
    double bf_banff_ccqeE0  = skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCCQEE0);

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

      // Get fit inputs object (chisq is calculated there)
      FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

      //need to fit for systematics
      FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,"minuit");
      fitter->DoTheFit();
      if(fit_inputs->BestFitChisquare() < chisq) {
	chisq = fit_inputs->BestFitChisquare();

	// report
	fparm->ShowOnlyCurrValuesInNextPrint();
	cout << *fparm;
      }
      delete fit_inputs;
      delete fitter;
    }//ifit
  }//fit for syst?
  else {
    // Get fit inputs object (chisq is calculated there)
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

    chisq = fit_inputs->Chisquare();

    // report
    fparm->ShowOnlyCurrValuesInNextPrint();
    cout << *fparm;
  }

  cout << endl << "chisq = " << chisq << endl;

}
