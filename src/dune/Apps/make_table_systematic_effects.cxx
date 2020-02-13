//____________________________________________________________________________
/*!

\app       valor_dune-make_table_systematic_effects

\brief     Vary the MC spectra by throwing sets of systematics and output a table containing the mean an RMS event
           rate in all samples.

\options     --output              : Where to save tweaked histograms to. [def: tweaks_[PARAMS].root]
             --ntoys               : Number of toy-MC experiments [def: 200]
             --physics-hypothesis  : Physics hypothesis [def: "3f_osc"]
             --fit-param-opt       : Options used when setting up the FitParamList
                                     [def: "NONE"]
             --config              : Analysis configuration file [required]
             --mesg-thresholds     : Message thresholds [def: "verbosity.xml"]
             --enu-reco-max        : Maximum Ereco to show on the plot [def: 10]
             --params              : Which parameters to tweak [def: allsyst]

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

\created   July 2019

\copyright (c) 2010-2019, The VALOR Fitting Group, All Rights Reserved

*/
//____________________________________________________________________________

// stl
#include <sstream>
#include <fstream>
#include <iomanip>

// root
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
#include "TPRegexp.h"

// fmwk
#include "valor/Conventions/Opt.h"
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
#include "valor/FitCore/UtlOscillation.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/VThreadPool.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Utilities/SingleSampleKineBinning.h"
#include "valor/Utilities/OscParamValues.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Ntuple/StdFitNtp.h"
#include "valor/KineBinning/VALORBinning.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

const double enu_reco_max = 10;

int main(int argc, char ** argv)
{

  valor::util::PrintBanner(Env::Instance()->Banner(),500);
 
  CmdLnArgList arg;
  arg.Read(argc,argv,"$VALOR_EXPT/config/apps/make_table_systematic_effects.arg:$VALOR_EXPT/config/apps/seeds.arg:$VALOR_EXPT/config/apps/oscparams.arg");
  LOG("valor_dune", pNOTICE) << arg;
  
  // Read the command line
  int ntoys               = arg.GetInt("ntoys");
  int n_threads           = arg.GetInt("nthreads");
  string hypname          = arg.GetString("physics-hypothesis");
  string fplopt           = arg.GetString("fit-param-opt");
  string output_rootfile  = arg.GetString("output");
  string params_to_tweak  = arg.GetString("params");
  bool rethrow_unphysical = arg.GetBool("rethrow-unphysical");
  
  std::string mesg_thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  Messenger::Instance()->Configure(mesg_thresholds);

  LOG("valor_dune",pINFO) << "Will rethrow unphysical toys";
  
  // Set up the thread pool
  if(n_threads>0) {//without n-cores specified, use as many as the host has
    util::VThreadPool::Instance()->SetNThreads(n_threads);
  }
  
  // Set the random seeds
  Rnd::Instance()->SetSeeds(arg);

  //change if something different to default (i.e. normal, sinsqtheta, first octant and dmsq12abs23 (NH) / dmsq12abs31 (IH))
  Opt::Instance()->gMassHierarchyAssumed = ( (arg.GetBool("inverted"))     ? conv::kNmhInverted      : conv::kNmhNormal);
  Opt::Instance()->gOscMixingParamChoice = ( (arg.GetBool("sinsq2theta"))  ? conv::kNuMixSinsq2Theta : conv::kNuMixSinsqTheta);
  Opt::Instance()->gOscOctantChoice      = ( (arg.GetBool("second-octant"))? conv::kNuOctantTwo      : conv::kNuOctantOne);

  conv::NuOscDmsqScheme_t dmsq_scheme = util::DmsqSchemeFromString(arg.GetString("dmsq-scheme"));
  if(dmsq_scheme == kNuDmsqUndefined) {
    LOG("valor_dune", pFATAL) << "Could not parse the input Dmsq scheme: [" << arg.GetString("dmsq-scheme") << "]";
    util::Die(-1);
  }
  Opt::Instance()->gOscDmsqScheme = dmsq_scheme;

  // Config
  std::string config_path = util::CompleteConfigPath(arg.GetString("config"));

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config_path << endl;
  cout << "- Physics hypothesis used : " << hypname << endl;

  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config_path.c_str());
  FitModelABC * model = elements->FitModel();
  CompoundFitModel * cm = dynamic_cast<CompoundFitModel*>(model);
  assert(cm);
  const util::NameRegistry * nr = util::NameRegistry::Instance();

  // Array to store the number of events in each sample for each toy
  double events [cm->GetNSubModels()][ntoys];

  std::vector<PhysParamId_t> param_set = util::GetParameterSet(params_to_tweak);
  assert(param_set.size());

  std::map<PhysParamId_t,bool> syst_param_map;
  std::vector<PhysParamId_t> osc_params;
  for (unsigned int ip = 0 ; ip < param_set.size() ; ip++) {
    PhysParamId_t p = param_set[ip];
    if(!util::IsSystematic(p)) { osc_params.push_back(p); }
    else                    { syst_param_map[p] = true; }
  }

  FitParamList * fp = util::GetFitParameterList(hypname.c_str(), fplopt.c_str());
  util::oscillation::SetOscillationParameterPriors(arg, *fp);

  // Generate pdfs
  for(int itoy = 0; itoy < ntoys; itoy++) {
    // Randomise the systematics in the model, and copy them to the FPL
    model->ResetSystematics("all");
    model->RandomizeSystematics(syst_param_map);
    fp->FromPhysParamList(model->Tweaks(),"current-silent");

    // Set input parameters
    RgKeyList input_args = arg.FindKeys("input-");
    TPRegexp input_match("input-(\\w*)");
    for (RgKey argname : input_args) { // RgKey is just a typename for std::string right now
      TObjArray * matches = input_match.MatchS(argname.c_str());

      if (matches->GetEntries() != 2) {
        LOG("valor_dune",pFATAL) << "Cannot parse argument "<<argname<<" to get parameter name";
        util::Die(-1);
      }
      TObjString * tos_param_name = dynamic_cast<TObjString*>(matches->At(1));
      std::string param_name(tos_param_name->String().Data());
      delete matches;

      conv::PhysParamId_t param = nr->PhysParamFromString(param_name);
      if (param == conv::kPhypUndefined) {
        LOG("valor_dune",pFATAL) << "Cannot find parameter to fix from argument string "<<argname;
        util::Die(-1);
      }

      double param_val = arg.GetDouble(argname);
      if ( std::fabs(param_val + 999. ) < conv::kASmallNum ) {
        param_val = Opt::Instance()->GetGlobalBestFit(param);
      }
      if(fp->Exists(param)) {
        fp->Update      (param, param_val);
        fp->SetInitValue(param, param_val);
      }
    }
    
    for(PhysParamId_t p : osc_params) {
      // TODO: Can this be removed now that we have osc params in the config files?
      std::cerr << "Osc param randomisation is not yet supported. Exiting" << std::endl;
      util::Die(-1);

      double rnd_value = util::oscillation::GetRandomizedOscParam(p,*fp);
      fp->Update(p, rnd_value);
      if(p == conv::kPhypMixing14    || p == conv::kPhypMixing15    || p == conv::kPhypMixing24    ||
         p == conv::kPhypMixing25    || p == conv::kPhypMixing34    || p == conv::kPhypMixing35    ||
         p == conv::kPhypDm41sq      || p == conv::kPhypDm51sq      || p == conv::kPhypDelta12     ||
         p == conv::kPhypDelta13     || p == conv::kPhypDelta15     || p == conv::kPhypDelta34     ||
         p == conv::kPhypDelta35     || p == conv::kPhypNSIepsMuTau || p == conv::kPhypMixingEE    ||
         p == conv::kPhypMixingMuMu  || p == conv::kPhypMixingMuE) {
        LOG("valor_dune", pWARN) << "Randomising sterile params may not work as intended";
      }
    }

    fp->FromPhysParamList(model->Tweaks(),"current-silent");

    // Check that the parameter values are physical. Rethrow if they aren't
    if (!util::AllWithinLimits(*fp) && rethrow_unphysical) {
      itoy--;
      continue;
    }

    if(itoy % 100 == 0) { cout << "** Toy-MC: #" << itoy << endl; }
    
    // Generate predictions
    model->Generate(*fp,"+SB");

    // Array to store results
    TObjArray * array = new TObjArray();
    array->SetOwner(true);

    // Write stuff out
    model->Write(array,"curr","pdft_rf +SB");
    
    // Pick results wanted 
    for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
      GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
      BeamConf_t b = this_model->GetBeam();
      Sample_t   s = this_model->GetSample();
      const SingleSampleKineBinning & binning = *this_model->GetBinning();
      const VALORBinning & reco_smp = *(binning.RecoSmpBinning());
      const std::vector<KineVar_t>  & reco_vars = reco_smp.KineVars();
      conv::KineSpace_t ks = util::KineSpaceFromKineVars(reco_vars);
      TString name = TString::Format("curr_pdft_sb_rf_%s_%s", this_model->Name().c_str(), util::AsString(ks).c_str());
      TH1 * curr_pdf = (TH1*) array->FindObject(name.Data());
      if (!curr_pdf) {
        cerr<<"Could not find histogram "<<name<<endl;
        exit(-1);
      }
      events[imodel][itoy] = curr_pdf->Integral();
    } //imodel
    array->Clear("C");
    delete array;  
  }//itoy

  // Setup a streamer which we'll use to write the Tex table to the output file
  stringstream streamer;
  streamer << fixed << setprecision(2);

  streamer << "\\begin{table} \n\\centering" << endl;
  streamer << "\t\\begin{tabular}{| c | c | c | c |} \\hline" << endl;
  streamer << "\t\tSample & Mean & RMS & 100\\times RMS/mean \\\\ \\hline" << endl;

  for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
    GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
    BeamConf_t b = this_model->GetBeam();
    Sample_t   s = this_model->GetSample();
    std::string submodel_name = this_model->Name();

    // Copy the events from the big array to one that's just for this submodel
    double submodel_events [ntoys];
    memcpy(submodel_events, events[imodel], sizeof(submodel_events));

    // Write the mean, RMS and their ratio to the Tex file
    double mean = TMath::Mean(ntoys, submodel_events);
    double rms  = TMath::RMS (ntoys, submodel_events);
    double frac = 100 * rms / mean;
    streamer << "\t\t" << submodel_name << " & " << mean << " & " << rms << " & " << frac << " \\\\" << endl;
  }//imodel

  streamer << "\t\t\\hline\n\t\\end{tabular}" << endl;
  streamer << "\t\\caption{ Randomisation statistics from " << ntoys
           << " pseudo-data sets, with the following set of parameters varied: " << params_to_tweak << " }" << endl;
  streamer << "\\end{table}" << endl;

  output_rootfile = util::ReplaceAll(output_rootfile,"[PARAMS]",params_to_tweak);
  LOG("valor_dune",pNOTICE)<<"Writing to file: "<<output_rootfile;

  streamer.seekg(0, ifstream::end);
  long size = streamer.tellg();
  streamer.seekg(0);
  ofstream outfile (output_rootfile);
  outfile.write (streamer.str().c_str(), size);
  outfile.close();
  
  std::cout <<"Done!"<<std::endl;
}
