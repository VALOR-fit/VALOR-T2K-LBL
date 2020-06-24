//
// ------------------------------------------------------------------------------------------------
// Take two trees of toy throw values, one SDK, one legacy. For the selected parameter, replace the legacy throws
// with the SDK throws and save to an output 'hybrid' tree. This tree can then be used to make SDK toys.
//
// Francis Bench
// Oct. 2019
// ------------------------------------------------------------------------------------------------
//

//___________________________________________________________________________
void make_sdk_legacy_hybrid_toy_throw_tree(
  const char * sdk_tree_file    = "migrated_run1-9d_fc_normal_cpv_reactor_0_1_40k_combined_toy_throws.root",
  const char * legacy_tree_file = "fit_fc_deltacp_normal_nexp_40000_-3.142_-999.000-job0-_0_39999_toy_throws.root",
  const char * out_tag          = "hybrid_ebo_run1-9d_fc_normal_cpv_reactor_0_1_40k_toy_throws.root",
  const char * replace_param    = "ebo",
  const int    n_toys           = 40000,
  const char * param_map        = "maps/t2xk_sdk_param_name_map.arg",
  const char * config           = "analysis/t2k_analysis_configuration_oa2018_run1-9d_v1_postfit.xml",
  const char * hypname          = "3f_osc",
  const char * fplopt           = "NONE"
)
{
    CompoundFitElements * elements = new CompoundFitElements();
    elements->Read(util::CompleteConfigPath(config), "nomc nosplines");
    FitModelABC* model = elements->FitModel();
    FitParamList* fpl = util::GetFitParameterList(hypname, fplopt);

    // Set up an array to hold all of the parameter values that we'll read/write
    double sdk_params[conv::knppmax];
    double leg_params[conv::knppmax];
    for(int ip = 0; ip < conv::knppmax; ip++){
        sdk_params[ip] = -999;
        leg_params[ip] = -999;
    }

    // Open the parameter map file and build a map of the SDK param names to the legacy param names
    std::map<std::string, std::string> par_name_map;
    std::ifstream map(util::CompleteConfigPath(param_map));
    std::string arg_file_line = "";
    while (std::getline(map, arg_file_line))
    {
        // Find the read parameter name that corresponds to this parameter
        vector<std::string> tokens = util::Tokens(arg_file_line, ",");
        std::string sdk_name = util::GetStringBetweenDelimiters(util::RemoveAllWhitespace(tokens[0]), "'", "'");
        std::string leg_name = util::GetStringBetweenDelimiters(util::RemoveAllWhitespace(tokens[1]), "'", "'");
        par_name_map[sdk_name] = leg_name;
    }

    // Read in the SDK tree
    TFile * sdk_file = new TFile(sdk_tree_file, "READ");
    TTree * tree_sdk_pars = (TTree*) sdk_file->Get("marg_parameter_values");
    tree_sdk_pars->SetDirectory(0);

    // Read in the legacy tree
    TFile * leg_file = new TFile(legacy_tree_file, "READ");
    TTree * tree_leg_pars = (TTree*) leg_file->Get("marg_parameter_values");
    tree_leg_pars->SetDirectory(0);

    // Set up the output hybrid tree
    TFile * out_file = new TFile(out_tag, "RECREATE");
    TTree * tree_out_pars = new TTree("marg_parameter_values", "Input parameters used for marginalisation");
    tree_out_pars->SetDirectory(0);

    // Set up the branches to read/write
    for (int ip = 0 ; ip < conv::knppmax ; ip++) {
        conv::PhysParamId_t p = (conv::PhysParamId_t)ip;
        if (p == conv::kPhypUndefined) { continue; }
        if (!fpl->Exists(p)) { continue; }
        std::string sdk_param_name = util::NameRegistry::Instance()->ParamAsStringVerySimple(p);

        if(par_name_map.find(sdk_param_name) == par_name_map.end()) {
          std::cout << "ERROR: map does not contain key [" << sdk_param_name << "]" << std::endl;
          util::Die(-1);
        }
        std::string leg_param_name = par_name_map.at(sdk_param_name);

        // Branches to read
        int branch_okay_sdk = tree_sdk_pars->SetBranchAddress(sdk_param_name.c_str(), &(sdk_params[p]));
        int branch_okay_leg = tree_leg_pars->SetBranchAddress(leg_param_name.c_str(), &(leg_params[p]));
        if(branch_okay_leg != 0) {
            std::cout << "ERROR: Could not read legacy input tree branch for param " << leg_param_name << std::endl;
            util::Die(-1);
        }
        if(branch_okay_sdk != 0) {
            std::cout << "ERROR: Could not read SDK input tree branch for param " << sdk_param_name << std::endl;
            util::Die(-1);
        }

        // Branch to write
        std::string leg_branch_info = leg_param_name + "/D";
        tree_out_pars->Branch(leg_param_name.c_str(), &(leg_params[ip]), leg_branch_info.c_str());
    }

    // Read the toy from the trees and copy over the relevant info from the SDK tree to the output legacy tree
    for(int itoy = 0; itoy < n_toys; ++itoy) {
        tree_sdk_pars->GetEntry(itoy);
        tree_leg_pars->GetEntry(itoy);
        for (int ip = 0 ; ip < conv::knppmax ; ip++) {
            conv::PhysParamId_t p = (conv::PhysParamId_t)ip;
            if (p == conv::kPhypUndefined) { continue; }
            if (!fpl->Exists(p)) { continue; }

            if (sdk_params[p] + 999 < conv::kASmallNum ||
                leg_params[p] + 999 < conv::kASmallNum) {
                std::cout << "The thrown value for param ["
                          << util::NameRegistry::Instance()->ParamAsStringVerySimple(p)
                          << "] could not be found in the SDK and/or legacy tree." << std::endl;
                util::Die(-1);
            }

            // Replace legacy param throws with SDK param throws
            std::string param_sdk_name = util::NameRegistry::Instance()->ParamAsStringVerySimple(p);
            if( param_sdk_name == (std::string)replace_param ||
                ((std::string)replace_param == "dm23sq-mixing23" && (param_sdk_name == "dm23sq" || param_sdk_name == "mixing23")) ) {

                double param_value = sdk_params[p];

                // Put in special cases here, where the param value needs changing too
                if (p == util::NameRegistry::Instance()->PhysParamFromString("sk_e_scale")) {
                    param_value -= 1;
                }

                if(itoy % 1000 == 0) {
                  std::cout << "Reading toy #" << itoy << ": replacing param [" << param_sdk_name
                            << "] in the legacy tree. Leg value = " << leg_params[p] << ", replaced with SDK value = "
                            << param_value << std::endl;
                }

                leg_params[p] = param_value;
            } // replacing this param
        } // ip
        tree_out_pars->Fill();
    } // itoy

    out_file->cd();
    tree_out_pars->Write();
    out_file->Close();
    delete out_file;

    sdk_file->Close();
    leg_file->Close();
    delete sdk_file;
    delete leg_file;

    std::cout << "Done!" << std::endl;
}