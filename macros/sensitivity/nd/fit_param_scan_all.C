//
// ------------------------------------------------------------------------------------------------
// The macro just runs fit_param_scan.C for all model parameters and merges the output files.
// See fit_param_scan.C for details.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//
	
void fit_param_scan_all
(
  const char * tag            = "scan",                                /* tag for output file */
  double       pmin           = 0.5,                                   /* min parameter value */
  double       pmax           = 1.5,                                   /* max parameter value */
  int          npv            = 11,                                    /* number of points between minimum and maximum */
  int          ntoymc         = 1,                                     /* number of toy-MC experiments */
  bool         stat_fluct     = false,                                 /* include statistical fluctuation in toy MC expts? */
  const char * param_values   = "nominal",                             /* decide how to set the values of params used for toy-MC generation (see UtlParam.cxx) */
  const char * param_list     = "standard",                            /* model parameter list (see UtlParam.cxx) */
  const char * config_file    = "dune_ndfit-2015v1-FGTv1.xml"    /* analysis configuration file */
)

{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  //
  // run fit_param_scan.C for each parameter
  //

  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(util::CompleteConfigPath(config_file),"");
  const std::vector<Parameter> & params = ParameterPool::Instance()->GetParams();

  for(int ip = 0; ip < param.size(); ip++) {
     PhysParamId_t p = params.At(ip);
      if(p == kPhypUndefined) continue;
      if(!util::IsSystematic(p)) continue;

     ostringstream command;
     command << ".x macros/fit/nd/fit_param_scan.C(\"" << tag << "\",\"" << util::NameRegistry::Instance()->ParamAsString(p) << "\","
             << pmin << "," << pmax << "," << npv << "," << ntoymc << "," << stat_fluct << ",\""
             << param_values << "\",\"" << param_list << "\",\"" << config_file << "\")";
     gROOT->ProcessLine(command.str().c_str());
  }


  // merge output files
  gSystem->Exec(Form("hadd %s.root %s-*",tag,tag));
}
