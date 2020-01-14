//____________________________________________________________________________
/*!

\app       valor_dune-tweak_pdf_sigma_var

\brief     Vary the MC spectra by setting systematics to +/- 1 and 3 sigma

\options     --output              : Where to save tweaked histograms to. [def: sigma_var_[PARAMS].root]
             --physics-hypothesis  : Physics hypothesis [def: "3f_osc"]
             --fit-param-opt       : Options used when setting up the FitParamList
                                     [def: "NONE"]
             --config              : Analysis configuration file [required]
             --mesg-thresholds     : Message thresholds [def: "verbosity.xml"]
             --params              : Which parameters to tweak [def: allsyst]
             --kinevar             : The kine var to store 1D PDFs for

             --input-mixing12,etc  : Option to give input values to osc params (in .arg)

             --load-prior-from-file : Option for file to give constrained prior matrix [ use default prior ]
             --prior-fit-number     : Option to specify the entry in the tree for above file [ 0 ]
             
             [--nthreads ]         : Number of CPU threads to use [ std::thread::hardware_concurrency ]
             

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
           Steve Dennis <s.r.dennis \at liverpool.ac.uk>
           University of Liverpool
           Francis Bench <f.bench \at liverpool.ac.uk>
           University of Liverpool

\created   June 2019

\copyright (c) 2010-2019, The VALOR Fitting Group, All Rights Reserved

*/
//____________________________________________________________________________

// root
#include "TFile.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"

// fmwk
#include "valor/Conventions/Opt.h"
#include "valor/Conventions/Names.h"
#include "valor/Conventions/ParameterPool.h"
#include "valor/Conventions/SystType.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/FitCore/CorrelatedSystematics.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/FitDataABC.h"
#include "valor/FitCore/FitInputs.h"
#include "valor/FitCore/FitterABC.h"
#include "valor/FitCore/UtlFit.h"
#include "valor/FitCore/CompoundFitElements.h"
#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/Utilities/OscParamValues.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/VThreadPool.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Utilities/UtlParam.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Ntuple/StdFitNtp.h"

using namespace std;
using namespace valor;
using namespace valor::conv;


bool t2k = false;


void SpecialCaseTweak(PhysParamId_t p, FitParamList * fp, double& tweak_p1, double& tweak_m1, double& tweak_p3, double& tweak_m3);
void GeneratePDFNominal(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv,
                        TFile * outfile, bool test_force_untuned);
void GeneratePDFTweaked(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv,
                        unsigned int ip, PhysParamId_t p, double tweak_dial, TFile * outfile,
                        std::string tag, bool test_force_untuned);

int main(int argc, char ** argv)
{

  valor::util::PrintBanner(Env::Instance()->Banner(),500);
 
  CmdLnArgList arg;
  arg.Read(argc,argv,"$VALOR_EXPT/config/apps/tweak_pdf_sigma_var.arg:$VALOR_EXPT/config/apps/seeds.arg:$VALOR_EXPT/config/apps/oscparams.arg");
  LOG("valor_dune", pNOTICE) << arg;
  
  // Read the command line
  int n_threads           = arg.GetInt("nthreads");
  string hypname          = arg.GetString("physics-hypothesis");
  string fplopt           = arg.GetString("fit-param-opt");
  string output_rootfile  = arg.GetString("output");
  string params_to_tweak  = arg.GetString("params");
  string kinevar_str      = arg.GetString("kinevar");

  // If true, then force all systematic parameters to their untuned values. This is for testing weights as is intended to be used with e.g. the DEBUG_SINGLE_BIN preprocessor directive in GSSFM
  bool test_force_untuned = arg.GetBool("test-force-untuned");

  KineVar_t kv = util::KineVarFromString(util::RemoveAllWhitespace(kinevar_str));
  if(kv == kvUndefined) {
    cerr << "ERROR: input kine var string [" << util::RemoveAllWhitespace(kinevar_str)
         << "] was not understood. Exiting." << std::endl;
    util::Die(-1);
  }

  string read_prior_from_file = arg.GetString("load-prior-from-file");
  int prior_fit_number = arg.GetInt("prior-fit-number");
  
  std::string mesg_thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  Messenger::Instance()->Configure(mesg_thresholds);
  
  // Set up the thread pool
  if(n_threads>0) {//without n-cores specified, use as many as the host has
    util::VThreadPool::Instance()->SetNThreads(n_threads);
  }
  
  // Set the random seeds
  Rnd::Instance()->SetSeeds(arg);
  
  // And OA parameters
  Opt::Instance()->gMassHierarchyAssumed = ( (arg.GetBool("inverted"))     ? conv::kNmhInverted      : conv::kNmhNormal);

  conv::NuOscDmsqScheme_t dmsq_scheme = util::DmsqSchemeFromString(arg.GetString("dmsq-scheme"));
  if(dmsq_scheme == kNuDmsqUndefined) {
    LOG("valor_dune", pFATAL) << "Could not parse the input Dmsq scheme: [" << arg.GetString("dmsq-scheme") << "]";
    util::Die(-1);
  }
  Opt::Instance()->gOscDmsqScheme = dmsq_scheme;

  // TODO: Update this with the new input param search
  double mixing12 = arg.GetDouble("input-mixing12");
  double mixing13 = arg.GetDouble("input-mixing13");
  double mixing23 = arg.GetDouble("input-mixing23");
  double dm12sq   = arg.GetDouble("input-dm12sq");
  double dm23sq   = arg.GetDouble("input-dm23sq");
  double deltacp  = arg.GetDouble("input-deltaCP");
  
  // Config
  std::string config_path = util::CompleteConfigPath(arg.GetString("config"));

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config_path << endl;
  cout << "- Physics hypothesis : " << hypname << endl;
  
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config_path.c_str());
  FitModelABC * model = elements->FitModel();
  CompoundFitModel * cm = dynamic_cast<CompoundFitModel*>(model);
  assert(cm);

  // Get the parameters to vary
  std::vector<PhysParamId_t> param_set = util::GetParameterSet(params_to_tweak, true);    // Only get exact matches to any input param names (otherwise problems can be caused)
  assert(param_set.size());

  // Make the output file
  TFile * outfile = NULL;
  output_rootfile = util::ReplaceAll(output_rootfile,"[PARAMS]",params_to_tweak) + "_" + util::AsString(kv) + ".root";
  outfile = new TFile(output_rootfile.c_str(),"RECREATE");
  LOG("valor_dune",pNOTICE)<<"Writing to file: "<<output_rootfile;

  const CorrelatedSystematics * syst = cm->GetSystematics();
  FitParamList * fp = util::GetFitParameterList(hypname.c_str(), fplopt.c_str());

  if (fp->Exists(conv::kPhypMixing12)) { fp->Update(conv::kPhypMixing12,mixing12); }
  if (fp->Exists(conv::kPhypMixing13)) { fp->Update(conv::kPhypMixing13,mixing13); }
  if (fp->Exists(conv::kPhypMixing23)) { fp->Update(conv::kPhypMixing23,mixing23); }
  if (fp->Exists(conv::kPhypDm12sq  )) { fp->Update(conv::kPhypDm12sq  ,dm12sq  ); }
  if (fp->Exists(conv::kPhypDm23sq  )) { fp->Update(conv::kPhypDm23sq  ,dm23sq  ); }
  if (fp->Exists(conv::kPhypCPV     )) { fp->Update(conv::kPhypCPV     ,deltacp ); }

  // Generate pdfs
  GeneratePDFNominal(model, cm, fp, kv, outfile, test_force_untuned);
  for (unsigned int ip = 0 ; ip < param_set.size() ; ip++) {
    PhysParamId_t p = param_set[ip];
    if(!util::IsSystematic(p)) { continue; }
    if(ParameterPool::Instance()->GetParam(p).GetType() == conv::kSystDisabled) { continue; }
    double bestfit = fp->InitValue(p);

    if(!syst->Contains(p)) {
      LOG("valor_dune", pFATAL) << "Could not get the 1 sigma error for param \'"
                                << util::NameRegistry::Instance()->ParamAsString(p) << "\'";
      util::Die(-1);
    }
    double sigma_err = syst->OneSigmaError(p);

    double tweak_p1 = bestfit + (1 * sigma_err);
    double tweak_m1 = bestfit - (1 * sigma_err);
    double tweak_p3 = bestfit + (3 * sigma_err);
    double tweak_m3 = bestfit - (3 * sigma_err);

    // Apply any special cases here (e.g. to force the tweaks to a certain value)
    SpecialCaseTweak(p, fp, tweak_p1, tweak_m1, tweak_p3, tweak_m3);

    LOG("valor_dune", pNOTICE)
      << "Tweak systematic parameter: "  << util::NameRegistry::Instance()->ParamAsString(p)
      << "\n  Nominal value = "          << bestfit
      << "\n  sigma error = "            << sigma_err
      << "\n  nominal +1 sigma tweak = " << tweak_p1
      << "\n  nominal -1 sigma tweak = " << tweak_m1
      << "\n  nominal +3 sigma tweak = " << tweak_p3
      << "\n  nominal -3 sigma tweak = " << tweak_m3;

    GeneratePDFTweaked(model, cm, fp, kv, ip, p, tweak_p1, outfile, "p1", test_force_untuned);
    GeneratePDFTweaked(model, cm, fp, kv, ip, p, tweak_m1, outfile, "m1", test_force_untuned);
    GeneratePDFTweaked(model, cm, fp, kv, ip, p, tweak_p3, outfile, "p3", test_force_untuned);
    GeneratePDFTweaked(model, cm, fp, kv, ip, p, tweak_m3, outfile, "m3", test_force_untuned);
  } // loop over params

  outfile->Close();
  delete outfile;

  LOG("valor_dune",pNOTICE)<<"Finished writing to file: "<<output_rootfile;
  std::cout << "Done!" << std::endl;
}
//________________________________________________________________________________________
void SpecialCaseTweak(conv::PhysParamId_t p, FitParamList * fp, double& tweak_p1, double& tweak_m1, double& tweak_p3, double& tweak_m3)
{
  if(t2k) {
    if(p == util::NameRegistry::Instance()->PhysParamFromString("2p2h_shape_O")) {
      double bestfit = fp->InitValue(p);
      double one_sigma_err   = 0.5;
      double three_sigma_err = 1.0;
      tweak_p1 = bestfit + one_sigma_err;
      tweak_m1 = bestfit - one_sigma_err;
      tweak_p3 = bestfit + three_sigma_err;
      tweak_m3 = bestfit - three_sigma_err;
    }
  }
}
//________________________________________________________________________________________
void GeneratePDFNominal(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv, TFile * outfile, bool test_force_untuned)
{
  cout << "\n\nGenerating nominal PDF:" << std::endl;
  outfile->cd();
  if(test_force_untuned) {
    model->ForceSystematicsToUntuned();
  }
  else {
    model->ResetSystematics("all");
  }
  fp->FromPhysParamList(model->Tweaks(),"current-silent");

  // Generate predictions
  cout << "with fit param list:" << std::endl;
  fp->Print(cout);
  model->Generate(*fp,"+SB");

  // Array to store results
  TObjArray * array = new TObjArray();
  array->SetOwner(true);

  // Write stuff out
  Option_t * write_opt = "pdfcs_rf-pdfts_rf+SB";
  if (util::IsTrue(kv)) { write_opt = "pdfcs_rf+SB"; }    // The true PDFs don't have totals in GSSFM, only the mode breakdown // TODO: util::binning::Make1DHistogram currently doesn't work with truth...
  model->Write(array,"nom",write_opt);

  // Pick results wanted
  for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
    GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
    Detector_t d = this_model->GetDetector();
    BeamConf_t b = this_model->GetBeam();
    Sample_t   s = this_model->GetSample();

    // Get the total PDF
    std::string model_name = util::nm::MdlName(d, b, s, util::KineSpaceFromKineVar(kv));
    if (util::IsReco(kv)) {
      TString name_tot     = TString::Format("nom_pdfts_sb_rf_%s_%s", model_name.c_str(), util::AsString(kv).c_str());
      TH1D * curr_pdf_tot = (TH1D*) array->FindObject(name_tot.Data());
      if (!curr_pdf_tot) {
        cerr<<"Could not find histogram "<<name_tot<<endl;
        exit(-1);
      }
      else {
        TString save_name = TString::Format("env_sb_%s_%s", model_name.c_str(), "nom");
        curr_pdf_tot->Write(save_name.Data());
      }
    }


    // Get the PDF broken down by mode contributions
    for(conv::ReacMode_t m : this_model->GetReactionModes()) {
      TString name_contrib = TString::Format("nom_pdfcs_sb_rf_%s_%s_%s", model_name.c_str(), util::AsString(kv).c_str(), util::AsString(m).c_str());
      TH1D * curr_pdf_contrib = (TH1D*) array->FindObject(name_contrib.Data());
      if (!curr_pdf_contrib) {
        cerr<<"Could not find histogram "<<name_contrib<<endl;
        exit(-1);
      }
      TString save_name = TString::Format("env_sb_%s_%s_%s", model_name.c_str(), util::AsString(m).c_str(), "nom");
      curr_pdf_contrib->Write(save_name.Data());
    } //imode
  } //imodel
  array->Clear("C");
  delete array;
}
//________________________________________________________________________________________
void GeneratePDFTweaked(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv,
                        unsigned int ip, PhysParamId_t p, double tweak_dial, TFile * outfile, std::string tag, bool test_force_untuned)
{
  outfile->cd();
  std::string param_name = util::NameRegistry::Instance()->ParamAsString(p);
  cout << "\n\nGenerating " << tag << " PDF with tweak dial = " << tweak_dial << " for param " << param_name << std::endl;
  if(test_force_untuned) {
    model->ForceSystematicsToUntuned();
  }
  else {
    model->ResetSystematics("all");
  }

  model->TweakSystematic(p, tweak_dial, "");
  fp->FromPhysParamList(model->Tweaks(),"current-silent");

  // Generate predictions
  model->Generate(*fp,"+SB");

  // Array to store results
  TObjArray * array = new TObjArray();
  array->SetOwner(true);

  // Write stuff out
  Option_t * write_opt = "pdfcs_rf-pdfts_rf+SB";
  if (util::IsTrue(kv)) { write_opt = "pdfcs_rf+SB"; }    // The true PDFs don't have totals in GSSFM, only the mode breakdown
  model->Write(array,"curr",write_opt);

  // Pick results wanted
  for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
    GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
    Detector_t d = this_model->GetDetector();
    BeamConf_t b = this_model->GetBeam();
    Sample_t   s = this_model->GetSample();

    // Get the total PDF
    std::string model_name = util::nm::MdlName(d, b, s, util::KineSpaceFromKineVar(kv));
    if (util::IsReco(kv)) {
      TString name = TString::Format("curr_pdfts_sb_rf_%s_%s", model_name.c_str(), util::AsString(kv).c_str());
      TH1D * curr_pdf = (TH1D*) array->FindObject(name.Data());
      if (!curr_pdf) {
        cerr<<"Could not find histogram "<<name<<endl;
        exit(-1);
      }
      TString save_name = TString::Format("env_sb_%s_%s_%s", model_name.c_str(), param_name.c_str(), tag.c_str());
      curr_pdf->Write(save_name.Data());
    }

    for (conv::ReacMode_t m : this_model->GetReactionModes()) {
      TString name_contrib = TString::Format("curr_pdfcs_sb_rf_%s_%s_%s", model_name.c_str(), util::AsString(kv).c_str(), util::AsString(m).c_str());
      TH1D * curr_pdf_contrib = (TH1D*) array->FindObject(name_contrib.Data());
      if (!curr_pdf_contrib) {
        cerr<<"Could not find histogram "<<name_contrib<<endl;
        exit(-1);
      }
      TString save_name = TString::Format("env_sb_%s_%s_%s_%s", model_name.c_str(), util::AsString(m).c_str(),
                        param_name.c_str(), tag.c_str());
      curr_pdf_contrib->Write(save_name.Data());
    } //mode
  } //imodel
  array->Clear("C");
  delete array;
}
//________________________________________________________________________________________
