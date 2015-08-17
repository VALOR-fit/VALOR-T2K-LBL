//____________________________________________________________________________
/*!

\class    valor::RunOpt

\brief    Some common VALOR macro options
          Based on $GENIE/src/Utils/RunOpt

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  14 Aug 2015

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _RUN_OPT_H_
#define _RUN_OPT_H_

#include "Parser/CmdLnArgParser.h"
#include "UtlForMacros.h"

#include <iostream>
#include <string>

using std::ostream;
using std::string;
using std::cerr;
using std::endl;

namespace valor {

// Type of which options to parse
typedef enum EParseType
{   
  kParseUndefined   = (1 << 0),
  kParseCommon      = (1 << 1),   // config file, run number, tag
  kParseToys        = (1 << 2),   // number of toys. throw stats? how to tweak systematics
  kParseGlobal      = (1 << 3),   // options to set in Opt
  kParse3FOscParams = (1 << 4)    // 3 flavour osc param options for both hierarchies (including MH marginalisation & priors)
} 
ParseType_t;
  
class RunOpt
{
public:
  RunOpt();
  RunOpt(int argc, char ** argv, int parseopt);
  virtual ~RunOpt();
  
  /*uncomment this if we want to make the class a singleton again (also remove the public constructor/destructor)
  static RunOpt * Instance(void);
  */
  
  // Read options from the command line. Call from all VALOR macros
  virtual void ReadFromCommandLine(int argc, char ** argv, int parseopt);

  // Get options set
  //common options
  string ConfigFile      () const { if(fParseType & kParseCommon) return fConfigFile; cerr << "kParseCommon not set" << endl; exit(-1); }
  string FileTag         () const { if(fParseType & kParseCommon) return fFileTag;    cerr << "kParseCommon not set" << endl; exit(-1); }
  int    Run             () const { if(fParseType & kParseCommon) return fRun;        cerr << "kParseCommon not set" << endl; exit(-1); }
  //toy options
  int    FirstExperiment     () const { if(fParseType & kParseToys) return fFirstExperiment;     cerr << "kParseToys not set" << endl; exit(-1); }
  int    LastExperiment      () const { if(fParseType & kParseToys) return fLastExperiment;      cerr << "kParseToys not set" << endl; exit(-1); }
  int    NExperiments        () const { if(fParseType & kParseToys) return fNExperiments;        cerr << "kParseToys not set" << endl; exit(-1); }
  bool   StatFluctuations    () const { if(fParseType & kParseToys) return fStatFluctuations;    cerr << "kParseToys not set" << endl; exit(-1); }
  int    SystTweakMode       () const { if(fParseType & kParseToys) return fSystTweakMode;       cerr << "kParseToys not set" << endl; exit(-1); }
  bool   RethrowIfUnphysical () const { if(fParseType & kParseToys) return fRethrowIfUnphysical; cerr << "kParseToys not set" << endl; exit(-1); }
  bool   VerboseSystThrows   () const { if(fParseType & kParseToys) return fVerboseSystThrows;   cerr << "kParseToys not set" << endl; exit(-1); }
  //global options
  string GlobalDmsqScheme    () const { if(fParseType & kParseGlobal) return fGlobalDmsqScheme;     cerr << "kParseGlobal not set" << endl; exit(-1); }
  string GlobalMassHierarchy () const { if(fParseType & kParseGlobal) return fGlobalMassHierarchy;  cerr << "kParseGlobal not set" << endl; exit(-1); }
  string GlobalMixParamScheme() const { if(fParseType & kParseGlobal) return fGlobalMixParamScheme; cerr << "kParseGlobal not set" << endl; exit(-1); }
  string GlobalMixParamOctant() const { if(fParseType & kParseGlobal) return fGlobalMixParamOctant; cerr << "kParseGlobal not set" << endl; exit(-1); }
  string GlobalCPVPhaseScheme() const { if(fParseType & kParseGlobal) return fGlobalCPVPhaseScheme; cerr << "kParseGlobal not set" << endl; exit(-1); }
  // 3f osc options
  string Hypothesis   () const { if(fParseType & kParse3FOscParams) return fHypothesis; cerr << "kParse3FOscParams not set" << endl; exit(-1); }
    
  // Print 
  virtual void   Print (ostream & stream, const int w = 30) const;
  friend ostream & operator << (ostream & stream, const RunOpt & opt);

 protected:
  // usage
  virtual void Usage(const int w = 40);

  void Init (void);
  void InitOscParams(double value = -99);

  // These functions get the parameters from the parser and set the
  //  relevant private members
  void ReadCommonOptions();
  void ReadToyOptions();
  void ReadGlobalOptions();
  void Read3FOscParamsOptions();

  // keep the parser
  CmdLnArgParser * fParser;
  // keep track of what mode we were read with
  int fParseType;
  // stores nh and ih
  string fHierStr[2];

  //Keep track of the options we've read
  // common options
  string fConfigFile;
  string fFileTag;
  int    fRun;
  // toy options
  int  fFirstExperiment;
  int  fLastExperiment;
  int  fNExperiments;
  bool fStatFluctuations;
  int  fSystTweakMode;
  bool fRethrowIfUnphysical;
  bool fVerboseSystThrows;
  // global options
  string fGlobalDmsqScheme;
  string fGlobalMassHierarchy;
  string fGlobalMixParamScheme;
  string fGlobalMixParamOctant;
  string fGlobalCPVPhaseScheme;
  // 3f osc options
  string fHypothesis;
  double fMix12 [2];
  double fMix13 [2];
  double fMix23 [2];
  double fCPV   [2];
  double fDm12sq[2];
  double fDm23sq[2];
  double fBetaBar;
  double fPriorProbNH;
  bool   fMargMH;
  string fPriorsFile;

  /* uncomment this if we want to make the class a singleton again (also remove the public constructor/destructor)
  // Self
  static RunOpt * fInstance;

  // Singleton class: constructors are private
  RunOpt();
  RunOpt(const RunOpt & opt);
  virtual ~RunOpt();

  // Clean-up
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (RunOpt::fInstance != 0) {
            delete RunOpt::fInstance;
            RunOpt::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
  */
};

ostream & operator << (ostream & stream, const RunOpt & opt);
}      // valor namespace
#endif // _RUN_OPT_H_
