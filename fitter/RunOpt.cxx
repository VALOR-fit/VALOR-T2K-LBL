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

#include "RunOpt.h"
#include "Conventions/Opt.h"
#include "Conventions/Constants.h"
#include "Conventions/PhysParamId.h"
#include "Utilities/Env.h"

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
  fConfigFile      = util::FormFilePath("analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml", "config");
  fFileTag         = "tag";
  fRun             = 0;
  //toy options
  fFirstExperiment     = 0;
  fLastExperiment      = 0;
  fNExperiments        = 1;
  fStatFluctuations    = true;
  fSystTweakMode       = 1;
  fRethrowIfUnphysical = true;
  fVerboseSystThrows   = false;
  //global options
  fGlobalDmsqScheme     = "auto";
  fGlobalMassHierarchy  = "normal";
  fGlobalMixParamScheme = "sinsqtheta";
  fGlobalMixParamOctant = "one";
  fGlobalCPVPhaseScheme = "delta";
  //3f osc options
  InitOscParams(-99);
  fHypothesis = "3f_osc";
  fMargMH     = true;
  fPriorsFile = util::FormFilePath("priors/joint_T2Kreact_run14/Priors_T2Kreact_run14_Mach3_22Maj15.root", "inputs"); //use this for default ON
  //fPriorsFile = ""; //use this for default OFF
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
  if(! (TMath::Abs(value + 999) < kASmallNum))
    fPriorProbNH = value;
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

  if(fParser->OptionExists("enable-stat-fluct"))
    fStatFluctuations = true;
  else if(fParser->OptionExists("disable-stat-fluct"))
    fStatFluctuations = false;

  if(fParser->OptionExists("syst-tweak-mode"))
    fSystTweakMode = fParser->ArgAsInt("syst-tweak-mode");

  if(fParser->OptionExists("enable-rethrow-unphysical-toys"))
    fRethrowIfUnphysical = true;
  else if(fParser->OptionExists("disable-rethrow-unphysical-toys"))
    fRethrowIfUnphysical = false;

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

  Opt::Instance()->SetValues(fGlobalDmsqScheme.c_str(),     fGlobalMassHierarchy.c_str(),
			   fGlobalMixParamScheme.c_str(), fGlobalMixParamOctant.c_str(),
			   fGlobalCPVPhaseScheme.c_str());
}
//____________________________________________________________________________
void valor::RunOpt::Read3FOscParamsOptions()
{
  if(fParser->OptionExists("hypothesis"))
    fHypothesis = fParser->ArgAsString("hypothesis");

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

  if(fParser->OptionExists("beta-bar"))
    fBetaBar = fParser->ArgAsDouble("beta-bar");

  if(fParser->OptionExists("normal-mass-hierarchy-prior"))
    fPriorProbNH = fParser->ArgAsDouble("normal-mass-hierarchy-prior");
  
  if(fParser->OptionExists("enable-marg-mass-hierarchy"))
    fMargMH = true;
  else if(fParser->OptionExists("disable-marg-mass-hierarchy"))
    fMargMH = false;

  if(fParser->OptionExists("disable-priors")) {
    fPriorsFile = "";
  }
  else if(fParser->OptionExists("priors-file")) {
    if(fParser->ArgAsString("priors-file").compare("off"))
      fPriorsFile = "";
    else
      fPriorsFile = util::FormFilePath(fParser->ArgAsString("priors-file").c_str(), "inputs");
  }
  util::OscParams::Instance()->SetValues(fMix12[0], fMix13[0], fMix23[0], fCPV[0], fDm12sq[0], fDm23sq[0],
					 fBetaBar,
					 fMix12[1], fMix13[1], fMix23[1], fCPV[1], fDm12sq[1], fDm23sq[1],
					 fPriorProbNH, fMargMH
					 );
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


