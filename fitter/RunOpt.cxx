//____________________________________________________________________________
/*
  \copyright (c) 2010-2014, The VALOR group.

  >> ChangeLog:

*/
//____________________________________________________________________________

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

#include <string>

#include "RunOpt.h"
#include "Conventions/Opt.h"
#include "Conventions/Constants.h"
#include "Conventions/PhysParamId.h"
#include "Utilities/Env.h"

#include "UtlFitParamListFactory.h"

#include "TMath.h"
#include "TString.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

using namespace valor;
using namespace valor::conv;

//____________________________________________________________________________
namespace valor {
  ostream & operator << (ostream & stream, const RunOpt & opt)
  {
    opt.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
valor::RunOpt::RunOpt()
{
  Init();
}
//____________________________________________________________________________
valor::RunOpt::RunOpt(int argc, char ** argv, int parseopt)
{
  Init();
  ReadFromCommandLine(argc, argv, parseopt);
  Print(cout);
}
//____________________________________________________________________________  
valor::RunOpt::~RunOpt()
{
  delete fParser;
}
//____________________________________________________________________________
/* uncomment this if we want to make the class a singleton again (also remove the public constructor/destructor)
RunOpt * valor::RunOpt::fInstance = 0;
//____________________________________________________________________________
valor::RunOpt::RunOpt()
{
  fInstance = 0;

  this->Init();
}
//____________________________________________________________________________
valor::RunOpt::~RunOpt()
{
  fInstance = 0;
}
//____________________________________________________________________________
RunOpt * valor::RunOpt::Instance()
{
  if(fInstance == 0) {
    static valor::RunOpt::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new RunOpt;
  }
  return fInstance;
}
*/
//____________________________________________________________________________
void valor::RunOpt::Init(void)
{
  fParser     = 0;
  fParseType  = kParseUndefined;
  fHierStr[0] = "nh";
  fHierStr[1] = "ih";
  //common options
  fConfigFile      = util::FormFilePath("", "config");
  fFileTag         = "tag";
  fRun             = 0;
  //toy options
  fFirstExperiment     = 0;
  fLastExperiment      = 0;
  fNExperiments        = 1;
  fStatFluctuations    = true;
  fSystTweakMode       = 1;
  fRethrowIfUnphysical = true;
  fVerboseSystThrows   = false; // Print out the output systematic tweaks
  //global options
  fGlobalDmsqScheme     = "auto";
  fGlobalMassHierarchy  = "normal";
  fGlobalMixParamScheme = "sinsqtheta";
  fGlobalMixParamOctant = "one";
  fGlobalCPVPhaseScheme = "delta";
  //3f osc options
  InitOscParams(-99);
  fHypothesis = "3f_osc";

  
  // Priors

  fPriorsFile = ""; //use this for default OFF
  //fPriorsFile = util::FormFilePath("priors/joint_T2Kreact_run14/Priors_T2Kreact_run14_Mach3_22Maj15.root", "inputs"); //use this for default ON

  fMargMH       = false;
  fPriorProbNH  = 0.5;

  // NEW
  fMargMixing12 = false;
  fMargMixing13 = false;
  fMargMixing23 = false;
  fMargDm12sq   = false;
  fMargDm23sq   = false;
  fMargDeltaCP  = false;

  fPriorMixing12 = "";
  fPriorMixing13 = "";
  fPriorMixing23 = "";
  fPriorDm12sq   = "";
  fPriorDm23sq   = "";
  fPriorDeltaCP  = "";
  
  fMixModeMixing12 = "";
  fMixModeMixing23 = "";
  fMixModeMixing13 = "";

  for(int i=0;i<2;i++){
    fPriorParMixing12[i] = conv::kBadResult;
    fPriorParMixing13[i] = conv::kBadResult;
    fPriorParMixing23[i] = conv::kBadResult;
    fPriorParDm12sq[i]   = conv::kBadResult;
    fPriorParDm23sq[i]   = conv::kBadResult;
    fPriorParDeltaCP[i]  = conv::kBadResult;
  }
  fSynchronizeIHtoNH = true;
  
  // fit marg options
  fXmlfile = "";
  fRootfile = "";
  fUsing_xml_file = false;
  fDoToys = true;
  fFit_param_list = "nuosc3f+nuisance";
  fInfilename     = "submit_prova_TemplVArray";
  fNfiles         = 9;
  fNtoysMarg      = 100;
  fFreeparname0   = "deltacp";
  fFreeparmin0    = -3.14159265358979312e+00;
  fFreeparmax0    = +3.14159265358979312e+00;
  fNptScan0       = 9;
  fFreeparname1   = "";
  fFreeparmin1    = 0.;
  fFreeparmax1    = 0.;
  fNptScan1       = 1;
  fDoMargMH       = false;
  fHierarchy_marg = "normal";
  fTypeError      = "";
  fDrawFullChisq  = true;
  fTeststat       = "rate+shape";
}
//____________________________________________________________________________
void valor::RunOpt::InitOscParams(double value)
{
  for(int i = 0; i < 2; i++) {
    fMix12 [i] = value;
    fMix13 [i] = value;
    fMix23 [i] = value;
    fCPV   [i] = value;
    fDm12sq[i] = value;
    fDm23sq[i] = value;
  }//i
  fBetaBar     = value;
  
  //randomizing the prior makes no sense, so only set it if it's not -999
  //if(! (TMath::Abs(value + 999) < kASmallNum))
  //fPriorProbNH = value;
}
//____________________________________________________________________________
void valor::RunOpt::ReadFromCommandLine(int argc, char ** argv, int parseopt)
{
  //save which options we're using
  fParseType = parseopt;
  if(fParseType & kParseUndefined) {
    cerr << "kParseUndefined given to valor::RunOpt::ReadFromCommandLine()" << endl;
    exit(-1);
  }

  //get the parser
  fParser = new CmdLnArgParser(argc,argv);

  if(fParser->OptionExists("help") || fParser->OptionExists('h')) {
    Usage();
    exit(1);
  }
  
  //read the options of each group
  if(parseopt & kParseCommon)
    ReadCommonOptions();
  if(parseopt & kParseToys)
    ReadToyOptions();
  if(parseopt & kParseGlobal)
    ReadGlobalOptions();
  if(parseopt & kParse3FOscParams)
    Read3FOscParamsOptions();

  // NEW FIT MARG
  if(parseopt & kParseFitMarg)
    ReadFitMargOptions();
  //
}
//____________________________________________________________________________
void valor::RunOpt::ReadCommonOptions()
{
  if(fParser->OptionExists("file-tag"))
    fFileTag = fParser->ArgAsString("file-tag");
  
  if(fParser->OptionExists("run"))
    fRun = fParser->ArgAsInt("run");

  if(fParser->OptionExists("config-file"))
    fConfigFile = util::FormFilePath(fParser->ArgAsString("config-file").c_str(), "config");
}
//____________________________________________________________________________
void valor::RunOpt::ReadToyOptions()
{
  if(fParser->OptionExists("first-experiment")) {
    if(!fParser->OptionExists("last-experiment")) {
      cerr << "last-experiment is required if first-experiment is given" << endl;
      exit(1);
    }
    if(fParser->OptionExists("n-experiments")) {
      cerr << "n-experiments is not allowed if first-experiment is given" << endl;
      exit(1);
    }    
    fFirstExperiment = fParser->ArgAsInt("first-experiment");
  }
  
  if(fParser->OptionExists("last-experiment")) {
    if(!fParser->OptionExists("first-experiment")) {
      cerr << "first-experiment is required if last-experiment is given" << endl;
      exit(1);
    }
    fLastExperiment = fParser->ArgAsInt("last-experiment");
  }
  
  if(fParser->OptionExists("n-experiments")) {
    fNExperiments = fParser->ArgAsInt("n-experiments");
    //set the 'other' number of experiments option
    fFirstExperiment = 0;
    fLastExperiment  = fNExperiments - 1;
  }
  else
      fNExperiments = fLastExperiment - fFirstExperiment + 1;

  if( fParser->OptionExists("enable-stat-fluct") &&
      fParser->OptionExists("disable-stat-fluct")
      ){
    cerr << endl;
    cerr << "Statistical fluctuations are both enabled and disabled" << endl;
    cerr << endl;
    exit(1);
  }
  else if(fParser->OptionExists("enable-stat-fluct"))
    fStatFluctuations = true;
  else if(fParser->OptionExists("disable-stat-fluct"))
    fStatFluctuations = false;
  else{
    cerr << endl;
    cerr << "You must enable or disable the statistical fluctuations" << endl;
    cerr << "Use --enable-stat-fluct or --disable-stat-fluct" << endl;
    cerr << endl;
    exit(1);
  }

  if(fParser->OptionExists("syst-tweak-mode"))
    fSystTweakMode = fParser->ArgAsInt("syst-tweak-mode");
  
  if( fParser->OptionExists("enable-rethrow-unphysical-toys") &&
      fParser->OptionExists("disable-rethrow-unphysical-toys")
      ){
    cerr << endl;
    cerr << "Rethrowing of physical and unphysical toys is both enabled and disabled" << endl;
    cerr << endl;
    exit(1);
  }
  else if(fParser->OptionExists("enable-rethrow-unphysical-toys"))
    fRethrowIfUnphysical = true;
  else if(fParser->OptionExists("disable-rethrow-unphysical-toys"))
    fRethrowIfUnphysical = false;
  else{
    cerr << endl;
    cerr << "You must enable or disable the rethrowing of unphysical toys" << endl;
    cerr << "Use --enable-rethrow-unphysical-toys or --disable-rethrow-unphysical-toys" << endl;
    cerr << endl;
    exit(1);
  }

  if(fParser->OptionExists("enable-verbose-syst-throws"))
    fVerboseSystThrows = true;
  else if(fParser->OptionExists("disable-verbose-syst-throws"))
    fVerboseSystThrows = false;
}
//____________________________________________________________________________
void valor::RunOpt::ReadGlobalOptions()
{
  if(fParser->OptionExists("dmsq-scheme"))
    fGlobalDmsqScheme = fParser->ArgAsString("dmsq-scheme");

  if(fParser->OptionExists("mix-param-scheme"))
    fGlobalMixParamScheme = fParser->ArgAsString("mix-param-scheme");

  if(fParser->OptionExists("cpv-phase-scheme"))
    fGlobalCPVPhaseScheme = fParser->ArgAsString("cpv-phase-scheme");

  if(fParser->OptionExists("mass-hierarchy"))
    fGlobalMassHierarchy = fParser->ArgAsString("mass-hierarchy");

  if(fParser->OptionExists("mix-param-octant"))
    fGlobalMixParamOctant = fParser->ArgAsString("mix-param-octant");

  cout << "ReadGlobalOptions" << endl;

  Opt::Instance()->SetValues(fGlobalDmsqScheme.c_str(),     fGlobalMassHierarchy.c_str(),
			   fGlobalMixParamScheme.c_str(), fGlobalMixParamOctant.c_str(),
			   fGlobalCPVPhaseScheme.c_str());
}
//____________________________________________________________________________
void valor::RunOpt::Read3FOscParamsOptions()
{

  if(fParser->OptionExists("hypothesis"))
    fHypothesis = fParser->ArgAsString("hypothesis");

  //
  // Keep the following order:
  // - randomize-all-oscparams
  // - default-all-oscparams
  // - set each single osc params
  //

  if(fParser->OptionExists("randomize-all-oscparams"))
    InitOscParams(-999);

  if(fParser->OptionExists("default-all-oscparams"))
    InitOscParams(-99);
   
  for(int i = 0; i < 2; i++) {
    if(fParser->OptionExists(Form("mixing12-%s", fHierStr[i].c_str())))
      fMix12[i] = fParser->ArgAsDouble(Form("mixing12-%s", fHierStr[i].c_str()));
    
    if(fParser->OptionExists(Form("mixing13-%s", fHierStr[i].c_str())))
      fMix13[i] = fParser->ArgAsDouble(Form("mixing13-%s", fHierStr[i].c_str()));
    
    if(fParser->OptionExists(Form("mixing23-%s", fHierStr[i].c_str())))
      fMix23[i] = fParser->ArgAsDouble(Form("mixing23-%s", fHierStr[i].c_str()));
    
    if(fParser->OptionExists(Form("cpv-%s", fHierStr[i].c_str())))            
      fCPV[i] = fParser->ArgAsDouble(Form("cpv-%s", fHierStr[i].c_str()));
    
    if(fParser->OptionExists(Form("dm12sq-%s", fHierStr[i].c_str())))
      fDm12sq[i] = fParser->ArgAsDouble(Form("dm12sq-%s", fHierStr[i].c_str()));
    
    if(fParser->OptionExists(Form("dm23sq-%s", fHierStr[i].c_str())))
      fDm23sq[i] = fParser->ArgAsDouble(Form("dm23sq-%s", fHierStr[i].c_str()));
  }//i


  if( fParser->OptionExists("enable-synch-osc-ih-to-nh") && 
      fParser->OptionExists("disable-synch-osc-ih-to-nh") 
      ){
    cerr << endl;
    cerr << "You both enabled and disabled the synchronization of IH to NH osc params" << endl;
    cerr << endl;
    exit(1);
  }
  if(fParser->OptionExists("enable-synch-osc-ih-to-nh"))
    fSynchronizeIHtoNH = true;
  else if( fParser->OptionExists("disable-synch-osc-ih-to-nh") )
    fSynchronizeIHtoNH = false;
  else{
    cerr << endl;
    cerr << "You didn't enabled or disabled the synchronization of IH to NH osc params" << endl;
    cerr << "Use --enable-synch-osc-ih-to-nh or --disable-synch-osc-ih-to-nh" << endl;
    cerr << endl;
    exit(1);
  }

  if(fParser->OptionExists("beta-bar"))
    fBetaBar = fParser->ArgAsDouble("beta-bar");

  if(fParser->OptionExists("normal-mass-hierarchy-prior"))
    fPriorProbNH = fParser->ArgAsDouble("normal-mass-hierarchy-prior");
  
  if( fParser->OptionExists("enable-marg-mass-hierarchy") &&
      fParser->OptionExists("disable-marg-mass-hierarchy")
      ){
    cerr << endl;
    cerr << "You both disabled and enabled the marginalization of the mass hierarchy" << endl;
    cerr << endl;
    exit(1);
  }
  else if(fParser->OptionExists("enable-marg-mass-hierarchy"))
    fMargMH = true;
  else if(fParser->OptionExists("disable-marg-mass-hierarchy"))
    fMargMH = false;
  else{
    cerr << endl;
    cerr << "You must enable or disable the marginalization of mass hierarchy" << endl;
    cerr << "Use --enable-marg-mass-hierarchy or --disable-marg-mass-hierarchy" << endl;
    cerr << endl;
    exit(1);
  }

  // Set the priors

  // NEW 
  
  if( fParser->OptionExists("prior-t2konly") && 
      fParser->OptionExists("prior-reactors")
      ){
    cerr << endl;
    cerr << "Both T2K-only and T2K+reactors priors are called" << endl;
    cerr << "You can call only one of them" << endl;
    cerr << endl;
    exit(1);
  }
  else if(fParser->OptionExists("prior-reactors")){ 

    fMargMixing12        = true;
    fPriorMixing12       = "gaus";   
    fPriorParMixing12[0] = 0.846;
    fPriorParMixing12[1] = 0.021;
    fMixModeMixing12     = "sinsq2th";

    fMargMixing13        = true;
    fPriorMixing13       = "gaus";   
    fPriorParMixing13[0] = 0.093;
    fPriorParMixing13[1] = 0.008;
    fMixModeMixing13     = "sinsq2th";

    fMargMixing23        = true;
    fPriorMixing23       = "flat";   
    fPriorParMixing23[0] = 0.3;
    fPriorParMixing23[1] = 0.7;
    fMixModeMixing23     = "sinsqth";

    fMargDm23sq          = true;
    fPriorDm23sq         = "flat";   
    fPriorParDm23sq[0]   = 2E-3;
    fPriorParDm23sq[1]   = 3E-3;

    fMargDm12sq          = true;
    fPriorDm12sq         = "gaus";   
    fPriorParDm12sq[0]   = 7.53E-5;
    fPriorParDm12sq[1]   = 0.18E-5;

    fMargDeltaCP         = true;
    fPriorDeltaCP        = "flat";   
    fPriorParDeltaCP[0]  = -TMath::Pi();
    fPriorParDeltaCP[1]  = +TMath::Pi();
  }
  else if(fParser->OptionExists("prior-t2konly")){

    fMargMixing12        = true;
    fPriorMixing12       = "gaus";   
    fPriorParMixing12[0] = 0.846;
    fPriorParMixing12[1] = 0.021;
    fMixModeMixing12     = "sinsq2th";

    fMargMixing13        = true;
    fPriorMixing13       = "gaus";   
    fPriorParMixing13[0] = 0.095;
    fPriorParMixing13[1] = 0.010;
    fMixModeMixing13     = "sinsq2th";

    fMargMixing23        = true;
    fPriorMixing23       = "flat";   
    fPriorParMixing23[0] = 0.3;
    fPriorParMixing23[1] = 0.7;
    fMixModeMixing23     = "sinsqth";

    fMargDm23sq          = true;
    fPriorDm23sq         = "flat";   
    fPriorParDm23sq[0]   = 2E-3;
    fPriorParDm23sq[1]   = 3E-3;

    fMargDm12sq          = true;
    fPriorDm12sq         = "gaus";   
    fPriorParDm12sq[0]   = 7.53E-5;
    fPriorParDm12sq[1]   = 0.18E-5;

    fMargDeltaCP         = true;
    fPriorDeltaCP        = "flat";   
    fPriorParDeltaCP[0]  = -TMath::Pi();
    fPriorParDeltaCP[1]  = +TMath::Pi();
  }
  else{
    if(!fParser->OptionExists("default-all-oscparams")){
      cerr << endl;
      cerr << "Priors are not defined" << endl;
      cerr << endl;
      exit(1);
    }
  }
  

  if(fParser->OptionExists("prior-mixing12")){
    fMargMixing12 = true;
    vector<string> input = fParser->ArgAsStringTokens("prior-mixing12",",");
    fPriorMixing12       = input[0];   
    fPriorParMixing12[0] = (double)atof(input[1].c_str());
    fPriorParMixing12[1] = (double)atof(input[2].c_str());
    fMixModeMixing12     = input[3];
  }
  if(fParser->OptionExists("prior-mixing13")){
    fMargMixing13 = true;
    vector<string> input = fParser->ArgAsStringTokens("prior-mixing13",",");
    fPriorMixing13       = input[0];   
    fPriorParMixing13[0] = (double)atof(input[1].c_str());
    fPriorParMixing13[1] = (double)atof(input[2].c_str());
    fMixModeMixing13     = input[3];
  }
  if(fParser->OptionExists("prior-mixing23")){
    fMargMixing23 = true;
    vector<string> input = fParser->ArgAsStringTokens("prior-mixing23",",");
    fPriorMixing23       = input[0];   
    fPriorParMixing23[0] = (double)atof(input[1].c_str());
    fPriorParMixing23[1] = (double)atof(input[2].c_str());
    fMixModeMixing23     = input[3];
  }
  if(fParser->OptionExists("prior-dmsq23")){
    fMargDm23sq = true;
    vector<string> input = fParser->ArgAsStringTokens("prior-dmsq23",",");
    fPriorDm23sq       = input[0];   
    fPriorParDm23sq[0] = (double)atof(input[1].c_str());
    fPriorParDm23sq[1] = (double)atof(input[2].c_str());
  }
  if(fParser->OptionExists("prior-dmsq12")){
    fMargDm12sq = true;
    vector<string> input = fParser->ArgAsStringTokens("prior-dmsq12",",");
    fPriorDm12sq       = input[0];   
    fPriorParDm12sq[0] = (double)atof(input[1].c_str());
    fPriorParDm12sq[1] = (double)atof(input[2].c_str());
  }
  if(fParser->OptionExists("prior-deltacp")){
    fMargDeltaCP = true;
    vector<string> input = fParser->ArgAsStringTokens("prior-deltacp",",");
    fPriorDeltaCP       = input[0];   
    fPriorParDeltaCP[0] = (double)atof(input[1].c_str());
    fPriorParDeltaCP[1] = (double)atof(input[2].c_str());
  }

  // Synchronize IH to NH osc params
  if(fSynchronizeIHtoNH){
    fMix12[1]  = fMix12[0];
    fMix13[1]  = fMix13[0];
    fMix23[1]  = fMix23[0];
    fCPV[1]    = fCPV[0];
    fDm12sq[1] = fDm12sq[0];
    fDm23sq[1] = fDm23sq[0];
  }

  // Set Osc param values
  util::OscParams::Instance()->SetValues(fMix12[0], fMix13[0], fMix23[0], fCPV[0], fDm12sq[0], fDm23sq[0],
					 fBetaBar,
					 fMix12[1], fMix13[1], fMix23[1], fCPV[1], fDm12sq[1], fDm23sq[1],
					 fPriorProbNH, fMargMH
					 );

  // Set priors of osc params to defaults, or set them to be randomised
  if(fMargMixing23) util::OscParams::Instance()->SetPrior(conv::kPhypMixing23, fPriorMixing23,  fPriorParMixing23[0], fPriorParMixing23[1], fMixModeMixing23);      
  if(fMargDm23sq)   util::OscParams::Instance()->SetPrior(conv::kPhypDm23sq,   fPriorDm23sq,    fPriorParDm23sq[0],   fPriorParDm23sq[1],   ""              );  
  if(fMargMixing13) util::OscParams::Instance()->SetPrior(conv::kPhypMixing13, fPriorMixing13,  fPriorParMixing13[0], fPriorParMixing13[1], fMixModeMixing13); // Convert to sinsqth
  if(fMargDeltaCP)  util::OscParams::Instance()->SetPrior(conv::kPhypCPV,      fPriorDeltaCP,  fPriorParDeltaCP[0],  fPriorParDeltaCP[1],  ""              ); 
  if(fMargMixing12) util::OscParams::Instance()->SetPrior(conv::kPhypMixing12, fPriorMixing12, fPriorParMixing12[0], fPriorParMixing12[1], fMixModeMixing12); // Convert to sinsqth
  if(fMargDm12sq)   util::OscParams::Instance()->SetPrior(conv::kPhypDm12sq,   fPriorDm12sq,   fPriorParDm23sq[0],   fPriorParDm23sq[1],   ""              ); 

  //if(fParser->OptionExists("disable-priors")) {
  //fPriorsFile = "";
  //}
  //else if(fParser->OptionExists("priors-file")) {
  //if(fParser->ArgAsString("priors-file").compare("off"))
  //fPriorsFile = "";
  //else
  //fPriorsFile = util::FormFilePath(fParser->ArgAsString("priors-file").c_str(), "inputs");
  //}
}


//____________________________________________________________________________
void valor::RunOpt::ReadFitMargOptions()
{
  if(fParser->OptionExists("do-toys")){
    fDoToys = true;
  }
  else{
    fDoToys = false;

    if(fParser->OptionExists("use-xml-file")){
      fXmlfile = fParser->ArgAsString("use-xml-file");
      fUsing_xml_file = true;
    }
    else if(fParser->OptionExists("use-root-file")){
      fRootfile = fParser->ArgAsString("use-root-file");
      fUsing_xml_file = false;
    }
    else{
      cerr << endl;
      cerr << "You not doing toys on the fly but you didn't pass neither root nor xml file" << endl;
      cerr << endl;
      exit(1);
    }
  }
  
  if(fParser->OptionExists("fit-param-list"))
    fFit_param_list = fParser->ArgAsString("fit-param-list");

  if(fParser->OptionExists("marg-file"))
    fInfilename = fParser->ArgAsString("marg-file");
  
  if(fParser->OptionExists("nfiles-marg"))
    fNfiles = fParser->ArgAsInt("nfiles-marg");

  if(fParser->OptionExists("ntoys-marg"))
    fNtoysMarg = fParser->ArgAsInt("ntoys-marg");
  
  if(fParser->OptionExists("freepar0-name"))
    fFreeparname0 = fParser->ArgAsString("freepar0-name");
  
  if(fParser->OptionExists("freepar1-name"))
    fFreeparname1 = fParser->ArgAsString("freepar1-name");
  
  if(fParser->OptionExists("freepar0-min"))
    fFreeparmin0 = fParser->ArgAsDouble("freepar0-min");

  if(fParser->OptionExists("freepar1-min"))
    fFreeparmin1 = fParser->ArgAsDouble("freepar1-min");

  if(fParser->OptionExists("freepar0-max"))
    fFreeparmax0 = fParser->ArgAsDouble("freepar0-max");

  if(fParser->OptionExists("freepar1-max"))
    fFreeparmax1 = fParser->ArgAsDouble("freepar1-max");


  // to check

  if(fParser->OptionExists("freepar0-npt"))
    fNptScan0 = fParser->ArgAsInt("freepar0-npt");
 
  if(fParser->OptionExists("freepar1-npt"))
    fNptScan1 = fParser->ArgAsInt("freepar1-npt");
 
  if(fParser->OptionExists("marginalize-mh"))
    fDoMargMH = true;
  else
    fDoMargMH = false;
  
  if(fParser->OptionExists("mass-hierarchy-marg"))
    fHierarchy_marg = fParser->ArgAsString("mass-hierarchy-marg");
  
  if(fParser->OptionExists("type-error"))
    fTypeError = fParser->ArgAsString("type-error");
  
  if(fParser->OptionExists("full-chisq"))
    fDrawFullChisq = true;
  else
    fDrawFullChisq = false;

  if(fParser->OptionExists("rate-shape"))
    fTeststat = "rate+shape";
  else if(fParser->OptionExists("rate-only"))
    fTeststat = "rate-only";
  else{
    cerr << endl;
    cerr << "The test statistic for marginalization is not defined" << endl;
    cerr << endl;
    exit(1);
  }
}
//____________________________________________________________________________
void valor::RunOpt::Usage(const int w)
{
  cout << std::left;
  cout << endl;
  std::stringstream note_stream;

  cout << " " << setw(w) << "--help -h" << "Print this help and exit" << endl
       << endl;

  if(fParseType & kParseCommon) {
    cout << "Common running options:" << endl
	 << " " << setw(w) << "--config-file VALUE"      << "analysis configuration file to use (in $T2KOSCANA/config/)" << endl
	 << " " << setw(w) << "--file-tag VALUE"         << "tag for output file(s)" << endl
	 << " " << setw(w) << "--run VALUE"              << "job number for output file(s)" << endl
	 << endl;
  }//kParseCommon
  if(fParseType & kParseToys) {
    cout << "Toy running options:" << endl
	 << " " << setw(w) << "--first-experiment VALUE"         << "first experiment to read/write" << endl
	 << " " << setw(w) << "--last-experiment VALUE"          << "last experiment to read/write (option required if --first-experiment given)" << endl
	 << " " << setw(w) << "--n-experiments VALUE"            << "number of experiments to read/write. Use this or first/last above" << endl
	 << " " << setw(w) << "--enable-stat-fluct"              << "use statistical fluctuations" << endl
	 << " " << setw(w) << "--syst-tweak-mode VALUE"          << "which systematics to throw" << endl
	 << " " << setw(w) << "--enable-rethrow-unphysical-toys" << "rethrow unphysical toys" << endl
	 << " " << setw(w) << "--enable-verbose-syst-throws"     << "print the SKExp everytime you tweak systematics" << endl
	 << endl;    
  }//kParseToys
  if(fParseType & kParseGlobal) {
    cout << "Global running options:" << endl
	 << " " << setw(w) << "--dmsq-scheme VALUE"        << "auto abs23 abs31 fogli" << endl
	 << " " << setw(w) << "--mix-param-scheme VALUE"   << "sinsqtheta sinsq2theta" << endl
	 << " " << setw(w) << "--cpv-phase-scheme VALUE"   << "delta sindelta" << endl
	 << " " << setw(w) << "--mass-hierarchy VALUE"     << "normal inverted" << endl
	 << " " << setw(w) << "--mix-param-octant VALUE"   << "one two" << endl
	 << endl;
  }//kParseGlobal
  if(fParseType & kParse3FOscParams) {
    cout << "3 flavour running options:" << endl
	 << " " << setw(w) << "--hypothesis VALUE"         << "3f_osc 3f_oscbar 3f_osc_normbar 3f_oscbar_normbar 2f_numu_disap" << endl
	 << " " << setw(w) << "--randomize-all-oscparams"  << "Randomize all oscillation parameters (i.e. set to -999)" << endl
	 << " " << setw(w) << "--default-all-oscparams"    << "Set all oscillation parameters to defaults (i.e. set to -99)" << endl;
    for(int i = 0; i < 2; i++) {
      cout << " " << setw(w) << Form("--mixing12-%s VALUE", fHierStr[i].c_str()) << endl
	   << " " << setw(w) << Form("--mixing13-%s VALUE", fHierStr[i].c_str()) << endl
	   << " " << setw(w) << Form("--mixing23-%s VALUE", fHierStr[i].c_str()) << endl
	   << " " << setw(w) << Form("--cpv-%s VALUE",      fHierStr[i].c_str()) << endl
	   << " " << setw(w) << Form("--dm12sq-%s VALUE",   fHierStr[i].c_str()) << endl
	   << " " << setw(w) << Form("--dm23sq-%s VALUE",   fHierStr[i].c_str()) << endl;
    }//i
    cout << " " << setw(w) << "--beta-bar VALUE" << endl
	 << " " << setw(w) << "--normal-mass-hierarchy-prior VALUE" << "The prior probability of normal hierarchy" << endl
	 << " " << setw(w) << "--enable-marg-mass-hierarchy"        << "Marginalise the mass hierarchy" << endl
	 << " " << setw(w) << "--priors-file VALUE"                 << "The file containing 2D priors (in $T2KOSCANA/inputs/). Use 'off' for no priors" << endl
	 << " " << setw(w) << "--disable-priors"                    << "Turn off the use of priors. Overrides --priors-file" << endl
	 << endl;
    note_stream << "Note: options are read in the order above, so you can randomize all then set e.g. beta-bar to something specific" << endl;
  }//kParse3FOscParams

  cout << note_stream.str();
  cout << "Note: all 'enable' options can be turned off by substituting 'disable'" << endl
       << "Note: 'enable' takes precedence" << endl
       << endl;
}
//____________________________________________________________________________
void valor::RunOpt::Print(ostream & stream, const int w) const
{
  cout << std::left;
  cout << endl;
  
  if(fParseType & kParseCommon) {
    stream << "Common running options:" << endl
	   << " - " << setw(w) << "Config file:"     << fConfigFile << endl
	   << " - " << setw(w) << "Output file tag:" << fFileTag << endl
	   << " - " << setw(w) << "Run number:"      << fRun << endl
	   << endl;
  }//kParseCommon
  if(fParseType & kParseToys) {
    stream << "Toy running options:" << endl
	   << " - " << setw(w) << "First experiment:"         << fFirstExperiment << endl
	   << " - " << setw(w) << "Last experiment:"          << fLastExperiment << endl
	   << " - " << setw(w) << "N experiments:"            << fNExperiments << endl
	   << " - " << setw(w) << "Statistical fluctuations:" << ((fStatFluctuations) ? "YES " : "NO") << endl
	   << " - " << setw(w) << "Randomize systematics:"    << util::SysThrowIdToString(fSystTweakMode) << endl
	   << " - " << setw(w) << "Rethrow unphysical:"       << ((fRethrowIfUnphysical) ? "YES " : "NO") << endl
	   << " - " << setw(w) << "Verbose skexp::Tweak():"   << ((fVerboseSystThrows) ? "YES " : "NO") << endl
	   << endl;
  }//kParseToys
  if(fParseType & kParseGlobal) {
    stream << "Global running options:" << endl
	   << " - " << setw(w) << "Dmsq scheme:"              << fGlobalDmsqScheme << endl
	   << " - " << setw(w) << "Mixing parameters scheme:" << fGlobalMixParamScheme << endl
	   << " - " << setw(w) << "CPV phase scheme:"         << fGlobalCPVPhaseScheme << endl
	   << " - " << setw(w) << "Mass hierarchy:"           << fGlobalMassHierarchy << endl
	   << " - " << setw(w) << "Mixing parameter octant:"  << fGlobalMixParamOctant << endl
	   << endl;
  }//kParseGlobal
  if(fParseType & kParse3FOscParams) {
    stream << "3f oscillation running options:" << endl
	   << " - " << setw(w) << "Hypothesis:" << fHypothesis << endl;
    for(int i = 0; i < 2; i++) {
      stream << " - " << setw(w) << Form("%s_%s", util::RawParamAsString(kPhypMixing12).c_str(), fHierStr[i].c_str()) << fMix12[i] << endl
	     << " - " << setw(w) << Form("%s_%s", util::RawParamAsString(kPhypMixing13).c_str(), fHierStr[i].c_str()) << fMix13[i] << endl
	     << " - " << setw(w) << Form("%s_%s", util::RawParamAsString(kPhypMixing23).c_str(), fHierStr[i].c_str()) << fMix23[i] << endl
	     << " - " << setw(w) << Form("%s_%s", util::RawParamAsString(kPhypCPV).c_str(),      fHierStr[i].c_str()) << fCPV[i] << endl
	     << " - " << setw(w) << Form("%s_%s", util::RawParamAsString(kPhypDm12sq).c_str(),   fHierStr[i].c_str()) << fDm12sq[i] << endl
	     << " - " << setw(w) << Form("%s_%s", util::RawParamAsString(kPhypDm23sq).c_str(),   fHierStr[i].c_str()) << fDm23sq[i] << endl;
    }//i
    stream << " - " << setw(w) << "Beta bar:"                    << fBetaBar << endl
	   << " - " << setw(w) << "Normal hierarchy prior prob:" << fPriorProbNH << endl
	   << " - " << setw(w) << "Marginalise mass hierarchy:"  << ((fMargMH) ? "YES " : "NO") << endl
	   << " - " << setw(w) << "Using priors:"                << (fPriorsFile.empty() ? "NO" : Form("YES: %s", fPriorsFile.c_str())) << endl
	   << endl;
  }//kParse3FOscParams
}
//___________________________________________________________________________


