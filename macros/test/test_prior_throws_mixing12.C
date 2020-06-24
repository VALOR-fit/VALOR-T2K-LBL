
// For testing the XML priors

// Francis Bench

// These functions are direct copies of those in util::oscillation,
// but for some reason valor_dune can't find this namespace...
void SetPriorsAll(std::string file_priors_full, FitParamList& fpl);
bool SetSinglePrior(TXMLEngine* xmleng, XMLNodePointer_t top_node, FitParamList& fpl);
double GetRandomizedOscParam(conv::PhysParamId_t p, FitParamList& fpl);

void test_prior_throws_mixing12
()
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using namespace valor;
  using namespace valor::conv;

  const int nthrows = 10;

  //Messenger::Instance()->Configure(util::CompleteConfigPath("verbosity.xml"));

  Opt::Instance()->gMassHierarchyAssumed = kNmhNormal;
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsqTheta;
  Opt::Instance()->gOscOctantChoice      = kNuOctantTwo;
  Opt::Instance()->gOscDmsqScheme        = kNuDmsqAuto;

  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(util::CompleteConfigPath("analysis/t2k_analysis_configuration_oa2018_run1-9d_v1_postfit.xml"));
  FitModelABC* model = elements->FitModel();

  FitParamList* fpl = util::GetFitParameterList("3f_osc", "");

  std::string file_priors_full = util::CompleteConfigPath("priors/t2k_run_1-9d_reactor_fc.xml");
  SetPriorsAll(file_priors_full, *fpl);

  Rnd::Instance()->SetSeeds(10, 20, 30, 40, 50); // toymcstat, toymcsyst, toymcosc, fitter, other

  for(int ithrow = 0; ithrow < nthrows; ++ithrow) {
    cout << "Throw #" << ithrow << ": sinsq(theta12) = " << GetRandomizedOscParam(kPhypMixing12, *fpl) << endl;
  }
}




void SetPriorsAll(std::string file_priors_full, FitParamList& fpl)
{
  // Parse the input XML priors file
  TXMLEngine* xmleng = new TXMLEngine;
  XMLDocPointer_t xmldoc = xmleng->ParseFile(file_priors_full.c_str());
  if (xmldoc==0) {
    cerr << "Can not parse input XML file [" << file_priors_full << "]" << endl;
    exit(1);
  }
  XMLNodePointer_t root_node = xmleng->DocGetRootElement(xmldoc);
  const char * list_name_cstr = xmleng->GetAttr(root_node,"name");
  cout << "Reading priors from element list named" << list_name_cstr << endl;

  XMLNodePointer_t curr_node = xmleng->GetChild(root_node);
  int npriors_set = 0;
  while(curr_node) {
    string curr_node_name = xmleng->GetNodeName(curr_node);
    if(!strcmp(curr_node_name.c_str(),"prior")) {
      bool success = SetSinglePrior(xmleng,curr_node,fpl);
      if (!success) {
        cerr << "Failed to read prior list. Exiting." << endl;
        exit(1);
      }
      npriors_set++;
    }
    curr_node = xmleng->GetNext(curr_node);
  }
}


bool SetSinglePrior(TXMLEngine* xmleng, XMLNodePointer_t top_node, FitParamList& fpl)
{
  // First read the information that's common to each prior type
  string param_names    = "";  // Comma-separated list of simple parameter names
  string prior_type     = "";  // 'flat', 'gaus', or 'user'
  string opt            = "";
  XMLNodePointer_t curr_node = xmleng->GetChild(top_node);
  while(curr_node) {
    string curr_node_name = xmleng->GetNodeName(curr_node);
    if(xmleng->GetNodeContent(curr_node)) {
      string content = util::TrimSpaces(xmleng->GetNodeContent(curr_node));
      if(curr_node_name == "params") { param_names = content; }
      if(curr_node_name == "type"  ) { prior_type = content; }
      if(curr_node_name == "scheme") { opt        = content; }
    }
    curr_node = xmleng->GetNext(curr_node);
  }

  // Now read the information that's specific to each prior type
  string prior_attribute1 = "";    // If type = flat, this is the lower end of the range. If type = gaus, this is the mean.  If type = user, this is the NH prior file
  string prior_attribute2 = "";    // If type = flat, this is the upper end of the range. If type = gaus, this is the sigma. If type = user, this is the IH prior file
  string prior_attribute3 = "";    // If type = user, this is the prior hist name. Otherwise this stays as an empty string
  curr_node = xmleng->GetChild(top_node);   // Put the current node back to the top of the list
  while(curr_node) {
    string curr_node_name = xmleng->GetNodeName(curr_node);
    if(xmleng->GetNodeContent(curr_node)) {
      string content = util::TrimSpaces(xmleng->GetNodeContent(curr_node));
      if(prior_type == "flat") {
        if(curr_node_name == "low" ) { prior_attribute1 = content; }
        if(curr_node_name == "high") { prior_attribute2 = content; }
      }
      else if(prior_type == "gaus") {
        if(curr_node_name == "mean" ) { prior_attribute1 = content; }
        if(curr_node_name == "sigma") { prior_attribute2 = content; }
      }
      else if(prior_type == "user") {
        if(curr_node_name == "file_nh"  ) { prior_attribute1 = content; }
        if(curr_node_name == "file_ih"  ) { prior_attribute2 = content; }
        if(curr_node_name == "hist_name") { prior_attribute3 = content; }
      }
      else {
        cerr << "Could not understand the prior type [" << prior_type << "]";
        exit(1);
      }
    }
    curr_node = xmleng->GetNext(curr_node);
  }

  // Check that we read the correct number of attributes for each prior type
  if(prior_type == "flat" && (prior_attribute1 == "" || prior_attribute2 == "" || prior_attribute3 != "")) {
    cerr << "Something went wrong when setting the flat prior named [" << param_names << "]";
    exit(1);
  }
  else if(prior_type == "gaus" && (prior_attribute1 == "" || prior_attribute2 == "" || prior_attribute3 != "")) {
    cerr << "Something went wrong when setting the Gaussian prior named [" << param_names << "]";
    exit(1);
  }
  else if(prior_type == "user" && (prior_attribute1 == "" || prior_attribute2 == "" || prior_attribute3 == "")) {
    cerr << "Something went wrong when setting the user prior named [" << param_names << "]";
    exit(1);
  }

  // Get the parameter(s) that we'll apply this prior to
  const util::NameRegistry * nr = util::NameRegistry::Instance();
  std::vector<std::string> tokens = util::Tokens(param_names, ",");
  std::vector<conv::PhysParamId_t> params;
  for(std::string param_name : tokens) {
    param_name = util::RemoveAllWhitespace(param_name);
    conv::PhysParamId_t p = nr->PhysParamFromString(param_name);
    if(p == conv::kPhypUndefined) {
      cerr << "Could not parse parameter [" << param_name << "]";
      exit(1);
    }
    params.push_back(p);
  }

  // Get the mixing/dmsq/phase scheme(s) for this prior (1 scheme per parameter)
  std::vector<std::string> scheme_tokens = util::Tokens(opt, ",");
  std::vector<std::string> schemes;
  for(std::string scheme_name : scheme_tokens) {
    scheme_name = util::RemoveAllWhitespace(scheme_name);
    schemes.push_back(scheme_name);
  }
  if(schemes.size() != params.size()) {
    cerr << "The prior named [" << param_names
                                  << "] has a mismatch in the number of schemes and parameters defined";
    exit(1);
  }

  // Check that the number of parameters is as expected
  if((prior_type == "flat" || prior_type == "gaus") && params.size() != 1) {
    cerr << "The number of parameters specified for flat or Gaussian prior named [" << param_names << "] was != 1";
    exit(1);
  }
  else if(!(params.size() == 1 || params.size() == 2)){
    cerr << "The number of parameters specified for user prior named [" << param_names << "] was not 1 or 2";
    exit(1);
  }

  // Set the prior
  if(prior_type == "flat") {
    double low  = atof(prior_attribute1.c_str());
    double high = atof(prior_attribute2.c_str());
    fpl.SetPriorMode   (params[0], "flat"    );
    fpl.SetMinValueFlat(params[0], low       );
    fpl.SetMaxValueFlat(params[0], high      );
    fpl.SetPriorScheme (params[0], schemes[0]);
    cout << "Setting flat prior for param " << param_names
                                   << ": [" << low << ", " << high << "], using scheme: "
                                   << schemes[0];
  }
  else if(prior_type == "gaus") {
    double mean  = atof(prior_attribute1.c_str());
    double sigma = atof(prior_attribute2.c_str());
    fpl.SetPriorMode  (params[0], "gaus"    );
    fpl.SetMeanGaus   (params[0], mean      );
    fpl.SetSigmaGaus  (params[0], sigma     );
    fpl.SetPriorScheme(params[0], schemes[0]);
    cout << "Setting Gaussian prior for param " << param_names
                                   << ": mean = " << mean << ", sigma = " << sigma
                                   << ", using scheme: " << schemes[0];
  }
  else {
    // TODO: This will not work if we throw the mass hierarchy!
    std::string prior_hist_file_name = "";
    if(Opt::Instance()->gMassHierarchyAssumed == util::HierarchyFromString("normal")) {
      prior_hist_file_name = prior_attribute1;
    }
    else {
      prior_hist_file_name = prior_attribute2;
    }
    TString base_dir = Env::Instance()->TopLevelDataDir().c_str();
    TFile* prior_root_file = new TFile((base_dir + "/" + prior_hist_file_name.c_str()).Data(), "read");
    if(!prior_root_file->IsOpen()) {
      cerr << "Could not find prior hist file [" << prior_hist_file_name << "]";
      exit(1);
    }

    if(params.size() == 1) {
      TH1D* prior_hist = (TH1D*) prior_root_file->Get(prior_attribute3.c_str());
      if(!prior_hist) {
        cerr << "Could not find 1D prior hist [" << prior_attribute3 << "] in file ["
                                      << prior_hist_file_name << "]";
        exit(1);
      }
      fpl.SetPriorMode  (params[0], "user");
      fpl.SetPriorHisto (params[0], prior_hist);
      fpl.SetPriorScheme(params[0], schemes[0]);
      cout << "Setting 1D hist prior for param " << param_names
                                     << ": file = " << prior_hist_file_name
                                     << ", using scheme: " << schemes[0];
    }
    else {
      // Set the axis titles to the parameter names they correspond to, so we know which is which when we come to throw from the hist
      TH2D* prior_hist = (TH2D*) prior_root_file->Get(prior_attribute3.c_str());
      if(!prior_hist) {
        cerr << "Could not find 1D prior hist [" << prior_attribute3 << "] in file ["
                                      << prior_hist_file_name << "]";
        exit(1);
      }
      prior_hist->GetXaxis()->SetTitle(nr->ParamAsStringVerySimple(params[0]).c_str());
      prior_hist->GetYaxis()->SetTitle(nr->ParamAsStringVerySimple(params[1]).c_str());
      fpl.SetPriorMode   (params[0], "user"    );
      fpl.SetPriorHisto2D(params[0], prior_hist);
      fpl.SetPriorScheme (params[0], schemes[0]);
      fpl.SetPriorMode   (params[1], "user"    );
      fpl.SetPriorHisto2D(params[1], prior_hist);
      fpl.SetPriorScheme (params[1], schemes[1]);
      cout << "Setting 2D hist prior for params " << param_names
                                     << ": file = " << prior_hist_file_name
                                     << "; Associating x-axis with " << nr->ParamAsStringVerySimple(params[0])
                                     << " (scheme = " << schemes[0] << ")"
                                     << "; Associating y-axis with " << nr->ParamAsStringVerySimple(params[1])
                                     << " (scheme = " << schemes[1] << ")";
    }
    prior_root_file->Close();
    delete prior_root_file;
  }

  return true;
}

double GetRandomizedOscParam(conv::PhysParamId_t p, FitParamList& fpl)
{
  // Here we throw the parameter according to the scheme set for its prior.
  // After throwing, the value is converted to whichever scheme we're using internally.
  string mode = fpl.GetPriorModeName(p);
  string prior_scheme_name  = fpl.GetPriorSchemeName(p);
  string target_scheme_name = ""; // Leave this as an empty string to allow util::ConvertOscParamToUsedScheme to set the target scheme according to the values set in Opt::Instance

  double value = -999;

  if(mode == "flat")
  {
    // Get the limits imposed by the prior
    double min = fpl.GetMinValueFlat(p);
    double max = fpl.GetMaxValueFlat(p);

    value = min - 1.0;
    do
    {
      value = Rnd::Instance()->RndToyMCosc().Uniform(min, max);
    }
    while(value < min || value > max);
  }
  else if(mode == "gaus")
  {
    double mean = fpl.GetMeanGaus(p);
    double sigma = fpl.GetSigmaGaus(p);
    // Get the absolute physical limits, these are hard-coded assuming some conventions,
    // so are converted to the conventions used for the throw in e.g. OscParamMinValue
    double min = util::OscParamMinValue(p, prior_scheme_name);
    double max = util::OscParamMaxValue(p, prior_scheme_name);

    value = min - 1.0;
    do
    {
      value = Rnd::Instance()->RndToyMCosc().Gaus(mean, sigma);
    }
    while(value < min || value > max);
  }
  else if(fpl.GetPriorHisto(p)) {
    double min = util::OscParamMinValue(p, prior_scheme_name);
    double max = util::OscParamMaxValue(p, prior_scheme_name);
    value = min - 1.0;
    do
    {
      value = fpl.GetPriorHisto(p)->GetRandom();
    }
    while(value < min || value > max);
  }
  else {
    cerr << "Invalid prior type, can not throw oscillation parameter, "
                                  << util::NameRegistry::Instance()->ParamAsStringVerySimple(p);
    exit(1);
  }

  return util::ConvertOscParamToUsedScheme(p, value, prior_scheme_name, target_scheme_name);
}