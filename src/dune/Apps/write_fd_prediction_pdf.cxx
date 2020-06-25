//____________________________________________________________________________                                                       
/*!                                                                                                                                  
                                                                                                                                      
\app       valor_dune-write_fd_prediction_pdf  

\brief     Write out generic single sample PDFs and MC templates and make fd predictions from nd fits 

\options
           --physics-hypothesis  : physics hypothesis (see UtlParam.cxx)
                                     [def: "3f_osc"]
           --fit-param-opt       : Options used when setting up the FitParamList
                                     [def: "NONE"]
           --config              : Analysis configuration file [required]
           --write-opt           : Options defining what to write-out [def:pdft_vf]
           --output-name         : Output file name [def:test]
           --mesg-thresholds     : Message thresholds [def: "verbosity.xml"] 
           --input               : Input values of the oscillation parameters
           --fix                 : Oscillation parameters to be fixed 
           --pot-override        : Fix all samples to given POT (use 1E21 when making template files)
           --force-dst           : Ignore the XML requesting you use prebuilt templates and read the DST anyway (useful for making the prebuilt templates)
           --generate-fd-pred    : Generate the far detector prediction? If so, this requires the next parameter too [options: generate-fd-pred || no-fd-pred (def)]
           --fit-output-file     : File containing the fit output - this is where the FitParamList comes from to Generate the FD prediction
           --fit-output-entry    : Desired entry in the tree [def: 0]

\author    VALOR DUNE people, aka Costas, Steve, Lorena
\created   March 27, 2014

\copyright (c) 2010-2017, The VALOR Fitting Group, All Rights Reserved

*/
//____________________________________________________________________________ 

// root
#include "TFile.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include "TPRegexp.h"

// fmwk
#include "valor/Conventions/Opt.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/FitDataABC.h"
#include "valor/FitCore/FitInputs.h"
#include "valor/FitCore/FitterABC.h"
#include "valor/FitCore/UtlFit.h"
#include "valor/FitCore/CompoundFitElements.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "valor/FitElements/GenericSingleSampleFitData.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Utilities/OscParamValues.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Ntuple/StdFitNtp.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

int main(int argc, char ** argv)
{

  valor::util::PrintBanner(Env::Instance()->Banner(),500);

  CmdLnArgList arg;
  arg.Read(argc,argv,"$VALOR_EXPT/config/apps/write_fd_prediction_pdf.arg:$VALOR_EXPT/config/apps/oscparams.arg:$VALOR_EXPT/config/apps/seeds.arg");
  LOG("write_fd_prediction_pdf", pNOTICE) << arg;

  string hypname          = arg.GetString("physics-hypothesis");
  string fplopt           = arg.GetString("fit-param-opt");
  string write_opt        = arg.GetString("write-opt");
  string mesg_thresholds  = arg.GetString("mesg-thresholds");
  string output_name      = arg.GetString("output-name");
  string config_file      = arg.GetString("config");
  string fit_output       = arg.GetString("fit-output-file");
  int    output_entry     = arg.GetInt("fit-output-entry");
  
  // If pot_override is set to a positive number, every sample will be initialised with the same POT - this number
  double pot_override     = arg.GetDouble("pot-override");
  if (pot_override > 0 ) {
    Opt::Instance()->MutableGlobalRegistry().Set("pot_override",pot_override);
  }
  
  TString readopt;
  // Check if we want to force the fit to read trees rather than templates, ignoring the config file
  bool force_dst = arg.GetBool("force-dst");
  if (force_dst) {
    readopt += " forcedst";
  }
 
  Messenger::Instance()->Configure(util::CompleteConfigPath(mesg_thresholds));

  bool generate_fd = arg.GetBool("generate-fd-pred");
  if(generate_fd){
    if(fit_output == "NONE"){
      LOG("write_fd_prediction_pdf",pFATAL) << "FD prediction requested but no fit result file given";
      util::Die(-1);
    }
  }

  // Set the random seeds                                                                                                          
  Rnd::Instance()->SetSeeds(arg);

  LOG("write_fd_prediction_pdf", pNOTICE) << *Rnd::Instance();

  // Create file to store outputs
  TString outfilename = TString::Format("%s.root", output_name.c_str());

  // Get ND model & configure                
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(util::CompleteConfigPath(config_file),readopt.Data());
  FitModelABC * model = elements->FitModel();
  const util::NameRegistry * nr = util::NameRegistry::Instance();

  // Get the specified fit parameter list (specifies which parameters are floated/fixed,                                              
  // parameter limits, initial values etc).
  FitParamList *fp = 0;
  
  // If we have not been asked to generate the FD prediction from fit result
  if(!generate_fd) fp = util::GetFitParameterList(hypname.c_str(), fplopt.c_str());
  // If we have, get the FPL from the fit result output file
  else {
    // Get FitParamList from the fit result file
    TFile f(fit_output.c_str(), "read");
    TTree * fitnt = (TTree*)f.Get("valor_dune_nd_ntp");
    fitnt->SetBranchAddress("best_fit_params", &fp);
    fitnt->GetEntry(output_entry);
    LOG("write_fd_prediction_pdf",pNOTICE) << "FitParamList read from fit result";

  }
  if(fp == 0){ // Make sure that worked
    LOG("write_fd_prediction_pdf",pFATAL) << "FitParamList not filled";
    util::Die(-1);
  }

  //change if something different to default (i.e. normal, sinsqtheta, first octant and dmsq12abs23)                                  
  Opt::Instance()->gMassHierarchyAssumed = ( (arg.GetBool("inverted"))     ? conv::kNmhInverted      : conv::kNmhNormal);
  Opt::Instance()->gOscMixingParamChoice = ( (arg.GetBool("sinsq2theta"))  ? conv::kNuMixSinsq2Theta : conv::kNuMixSinsqTheta);
  Opt::Instance()->gOscOctantChoice      = ( (arg.GetBool("second-octant"))? conv::kNuOctantTwo      : conv::kNuOctantOne);
  Opt::Instance()->gOscDmsqScheme        = ( (arg.GetBool("foglilisi"))    ? conv::kNuDmsqFogliLisi  : conv::kNuDmsq12andAbs23);

  // Update values of oscillation parameters 
  RgKeyList input_args = arg.FindKeys("input-");
  TPRegexp input_match("input-(\\w*)");
  for (RgKey argname : input_args) { // RgKey is just a typename for std::string right now
    TObjArray * matches = input_match.MatchS(argname.c_str());

    if (matches->GetEntries() != 2) {
      LOG("write_fd_prediction_pdf",pFATAL) << "Cannot parse argument "<<argname<<" to get parameter name";
      util::Die(-1);
    }
    TObjString * tos_param_name = dynamic_cast<TObjString*>(matches->At(1));
    std::string param_name(tos_param_name->String().Data());
    delete matches;

    conv::PhysParamId_t param = nr->PhysParamFromString(param_name);
    if (param == conv::kPhypUndefined) {
      LOG("write_fd_prediction_pdf",pFATAL) << "Cannot find parameter to fix from argument string "<<argname;
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

  // Generate predictions                                                                                                         
  Option_t * generate_opt = "+SB"; // also generate non-default (for fitting) predictions in extra binning schemes                
  model->Generate(*fp,generate_opt);

  CompoundFitModel * cmodel = dynamic_cast<CompoundFitModel*>(model);
  for (unsigned int i = 0 ; i < cmodel->GetNSubModels() ; i++) {
    GenericSingleSampleFitModel * gmodel = dynamic_cast<GenericSingleSampleFitModel*>(cmodel->GetSubModel(i));
    gmodel->UpdateObservableBinningFromSimple();
  }

  // Write stuff out                              
  model->Write(outfilename,"valor",write_opt.c_str());

  std::cout<<"Done. Have a good day!"<<std::endl << "GO VALOR!!!" << std::endl;
  
  fp->Print(cout);

}
