#!/usr/bin/env python

# Add parameters to the BANFF matrix
#
# Andy Chappell
# Francis Bench

import ROOT
import argparse
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

#_________________________________________________________________________________
def extend_banff(N, banff_matrix, banff_nominal, banff_names):
    size = len(banff_nominal)
    # Copy existing BANFF matrix to new matrix with space for new variables
    new_matrix = ROOT.TMatrixD(N, N)
    for x in range(size):
        for y in range(size):
            new_matrix[y][x] = banff_matrix[y][x]

    # Copy existing nominal values to new vector with space for new variables
    new_nominal = ROOT.TVectorD(N)
    for i in range(size):
        new_nominal[i] = banff_nominal[i]

    # Copy existing parameter names to new TObjArray with space for new variables
    new_names = ROOT.TObjArray(N)
    for name in banff_names:
        new_names.Add(ROOT.TObjString(name))

    return new_matrix, new_nominal, new_names
#_________________________________________________________________________________
def copy_eb_to_banff(banff_matrix, banff_nominal, banff_names):
    # NOTE: untuned values are added in reorder_banff_matrix.py, as the file supplied
    #       by the BANFF group does not contain untuned values, and this macro needs
    #       to output a file in the same format as the BANFF file.
    params_to_add = ["EB_dial_O_nu", "EB_dial_O_nubar"]
    nominal_vals = {
        "EB_dial_O_nu"    : 4.0,
        "EB_dial_O_nubar" : 0.0
    }
    untuned_vals = {
        "EB_dial_O_nu"    : 0.0,
        "EB_dial_O_nubar" : 0.0
    }
    error_vals = {
        "EB_dial_O_nu"    : 6.0,
        "EB_dial_O_nubar" : 6.0
    }
    correlation_vals = {
        "EB_dial_O_nu"    : { "EB_dial_O_nubar" : 0.70 },
        "EB_dial_O_nubar" : { "EB_dial_O_nu"    : 0.70 }
    }

    return copy_params_to_banff(banff_matrix, banff_nominal, banff_names, params_to_add,
                                nominal_vals, error_vals, correlation_vals)

#_________________________________________________________________________________
def copy_intermediate_run10_params_to_banff(banff_matrix, banff_nominal, banff_names):
    # NOTE: untuned values are added in reorder_banff_matrix.py, as the file supplied
    #       by the BANFF group does not contain untuned values, and this macro needs
    #       to output a file in the same format as the BANFF file.
    params_to_add = ["FEFQE", "FEFQEH", "FEFINEL", "FEFABS", "FEFCX", "CC_norm_nu", "CC_norm_nubar",
                     "CC_DIS_MultPi_Norm_Nu", "CC_DIS_MultPi_Norm_Nubar", "CC_Misc", "ISO_BKG_LowPPi",
                     "2p2h_Edep_lowEnu", "2p2h_Edep_highEnu", "2p2h_Edep_lowEnubar", "2p2h_Edep_highEnubar",
                     "Q2_norm_0", "Q2_norm_1", "Q2_norm_2", "Q2_norm_3", "Q2_norm_4",
                     "Q2_norm_5", "Q2_norm_6", "Q2_norm_7",
                     "EB_dial_O_nu", "EB_dial_O_nubar", "CC_BY_DIS", "CC_BY_MPi", "CC_AGKY_Mult"]
    nominal_vals = {
        "FEFQE"                    : 1.069,
        "FEFQEH"                   : 1.824,
        "FEFINEL"                  : 1.002,
        "FEFABS"                   : 1.404,
        "FEFCX"                    : 0.697,
        "CC_norm_nu"               : 1.0,
        "CC_norm_nubar"            : 1.0,
        "CC_DIS_MultPi_Norm_Nu"    : 1.0,
        "CC_DIS_MultPi_Norm_Nubar" : 1.0,
        "CC_Misc"                  : 1.0,
        "ISO_BKG_LowPPi"           : 0.96,
        "2p2h_Edep_lowEnu"         : 1.0,
        "2p2h_Edep_highEnu"        : 1.0,
        "2p2h_Edep_lowEnubar"      : 1.0,
        "2p2h_Edep_highEnubar"     : 1.0,
        "Q2_norm_0"                : 0.495,
        "Q2_norm_1"                : 0.695,
        "Q2_norm_2"                : 0.78,
        "Q2_norm_3"                : 0.89,
        "Q2_norm_4"                : 0.93,
        "Q2_norm_5"                : 1.0,
        "Q2_norm_6"                : 1.0,
        "Q2_norm_7"                : 1.0,
        "EB_dial_O_nu"             : 4.0,
        "EB_dial_O_nubar"          : 0.0,
        "CC_BY_DIS"                : 0.0,
        "CC_BY_MPi"                : 0.0,
        "CC_AGKY_Mult"             : 0.0
    }
    untuned_vals = {
        "FEFQE"                    : 1.069,
        "FEFQEH"                   : 1.824,
        "FEFINEL"                  : 1.002,
        "FEFABS"                   : 1.404,
        "FEFCX"                    : 0.697,
        "CC_norm_nu"               : 1.0,
        "CC_norm_nubar"            : 1.0,
        "CC_DIS_MultPi_Norm_Nu"    : 1.0,
        "CC_DIS_MultPi_Norm_Nubar" : 1.0,
        "CC_Misc"                  : 1.0,
        "ISO_BKG_LowPPi"           : 1.3,
        "2p2h_Edep_lowEnu"         : 1.0,
        "2p2h_Edep_highEnu"        : 1.0,
        "2p2h_Edep_lowEnubar"      : 1.0,
        "2p2h_Edep_highEnubar"     : 1.0,
        "Q2_norm_0"                : 1.0,
        "Q2_norm_1"                : 1.0,
        "Q2_norm_2"                : 1.0,
        "Q2_norm_3"                : 1.0,
        "Q2_norm_4"                : 1.0,
        "Q2_norm_5"                : 1.0,
        "Q2_norm_6"                : 1.0,
        "Q2_norm_7"                : 1.0,
        "EB_dial_O_nu"             : 0.0,
        "EB_dial_O_nubar"          : 0.0,
        "CC_BY_DIS"                : 0.0,
        "CC_BY_MPi"                : 0.0,
        "CC_AGKY_Mult"             : 0.0
    }
    error_vals = {
        "FEFQE"                    : 0.313,
        "FEFQEH"                   : 0.859,
        "FEFINEL"                  : 1.101,
        "FEFABS"                   : 0.432,
        "FEFCX"                    : 0.305,
        "CC_norm_nu"               : 0.02,
        "CC_norm_nubar"            : 0.01,
        "CC_DIS_MultPi_Norm_Nu"    : 0.035,
        "CC_DIS_MultPi_Norm_Nubar" : 0.065,
        "CC_Misc"                  : 1.0,
        "ISO_BKG_LowPPi"           : 1.3,
        "2p2h_Edep_lowEnu"         : 1.0,
        "2p2h_Edep_highEnu"        : 1.0,
        "2p2h_Edep_lowEnubar"      : 1.0,
        "2p2h_Edep_highEnubar"     : 1.0,
        "Q2_norm_0"                : 0.165,
        "Q2_norm_1"                : 0.145,
        "Q2_norm_2"                : 0.13,
        "Q2_norm_3"                : 0.15,
        "Q2_norm_4"                : 0.16,
        "Q2_norm_5"                : 0.11,
        "Q2_norm_6"                : 0.18,
        "Q2_norm_7"                : 0.40,
        "EB_dial_O_nu"             : 6.0,
        "EB_dial_O_nubar"          : 6.0,
        "CC_BY_DIS"                : 1.0,
        "CC_BY_MPi"                : 1.0,
        "CC_AGKY_Mult"             : 1.0
    }
    correlation_vals = {
        "FEFQE"                    : { "FEFQEH" : -0.03, "FEFINEL" :  0.03, "FEFABS"  : -0.27, "FEFCX"  : -0.20 },
        "FEFQEH"                   : { "FEFQE"  : -0.03, "FEFINEL" : -0.53, "FEFABS"  :  0.03, "FEFCX"  :  0.01 },
        "FEFINEL"                  : { "FEFQE"  :  0.03, "FEFQEH"  : -0.53, "FEFABS"  : -0.03, "FEFCX"  : -0.01 },
        "FEFABS"                   : { "FEFQE"  : -0.27, "FEFQEH"  :  0.03, "FEFINEL" : -0.03, "FEFCX"  : -0.31 },
        "FEFCX"                    : { "FEFQE"  : -0.20, "FEFQEH"  :  0.01, "FEFINEL" : -0.01, "FEFABS" : -0.31 },
        "CC_norm_nu"               : { "CC_norm_nubar"            : -0.999999 },
        "CC_norm_nubar"            : { "CC_norm_nu"               : -0.999999 },
        "CC_DIS_MultPi_Norm_Nu"    : { "CC_DIS_MultPi_Norm_Nubar" : 0.999999  },
        "CC_DIS_MultPi_Norm_Nubar" : { "CC_DIS_MultPi_Norm_Nu"    : 0.999999  },
        "EB_dial_O_nu"             : { "EB_dial_O_nubar"          : 0.70      },
        "EB_dial_O_nubar"          : { "EB_dial_O_nu"             : 0.70      }
    }

    return copy_params_to_banff(banff_matrix, banff_nominal, banff_names, params_to_add,
                                nominal_vals, error_vals, correlation_vals)
#_________________________________________________________________________________
def copy_params_to_banff(banff_matrix, banff_nominal, banff_names, params_to_add, nominal_vals, error_vals, correlation_vals):
    # Get extended matrices/vectors with the already-existing params as they were + space for the new params
    N = len(banff_nominal) + len(params_to_add)
    new_matrix, new_nominal, new_names = extend_banff(N, banff_matrix, banff_nominal, banff_names)

    # Add the new param names
    for name in params_to_add:
        new_names.Add(ROOT.TObjString(name))

    # Add the new param nominal values
    for i, name in enumerate(params_to_add):
        new_nominal[len(banff_nominal) + i] = nominal_vals[name]

    # Get a single list of 1 sigma errors
    new_errors = []
    for i in range(0,len(banff_nominal)):
        new_errors.append(ROOT.TMath.Sqrt(banff_matrix[i][i]))
    for name in params_to_add:
        new_errors.append(error_vals[name])

    # Calculate the new covariance values from the correlations and errors using
    # cov_{ij} = corr_{ij} * err_{i} * err_{j}
    for i, name_i in enumerate(new_names):
        # Convert the param names from TObjStrings to strings (otherwise string comparisons will fail)
        name_i = name_i.GetName()
        for j, name_j in enumerate(new_names):
            # Loop over all params (new + old) because new params might have correlations with old params.
            name_j = name_j.GetName()
            if name_i not in params_to_add and name_j not in params_to_add:
                # If both params are old, then we don't need to do anything
                continue

            corr_ij = 0.0
            if i == j:
                corr_ij = 1.0
            elif name_i in correlation_vals:
                if name_j in correlation_vals[name_i]:
                    corr_ij = correlation_vals[name_i][name_j]

            cov_ij = corr_ij * new_errors[i] * new_errors[j]
            new_matrix[i][j] = cov_ij

    return new_matrix, new_nominal, new_names
#_________________________________________________________________________________
def add_param_to_banff(banff_filename, outputname, var):
    banff_file = ROOT.TFile(banff_filename, "READ")

    prefit_matrix   = banff_file.Get("prefit_cov")
    prefit_nominal  = banff_file.Get("prefit_params")
    postfit_matrix  = banff_file.Get("postfit_cov")
    postfit_nominal = banff_file.Get("postfit_params")
    banff_names     = banff_file.Get("param_list")

    new_prefit_matrix   = None
    new_prefit_nominal  = None
    new_postfit_matrix  = None
    new_postfit_nominal = None
    new_names           = None

    if var == "eb":
        new_prefit_matrix , new_prefit_nominal , new_names = copy_eb_to_banff(prefit_matrix , prefit_nominal , banff_names)
        new_postfit_matrix, new_postfit_nominal, _         = copy_eb_to_banff(postfit_matrix, postfit_nominal, banff_names)
    elif var == "intermediate_run10":
        new_prefit_matrix , new_prefit_nominal , new_names = copy_intermediate_run10_params_to_banff(prefit_matrix , prefit_nominal , banff_names)
        new_postfit_matrix, new_postfit_nominal, _         = copy_intermediate_run10_params_to_banff(postfit_matrix, postfit_nominal, banff_names)
    else:
        print "Currently only adding Eb or FSI params to the BANFF matrix are supported. Exiting."
        sys.exit()

    print "Saving extended BANFF matrix to file", outputname
    new_banff_file = ROOT.TFile(outputname, "RECREATE")
    new_names          .Write("param_list", 1)
    new_prefit_nominal .Write("prefit_params")
    new_prefit_matrix  .Write("prefit_cov")
    new_postfit_nominal.Write("postfit_params")
    new_postfit_matrix .Write("postfit_cov")

    # Copy objects that don't change
    for mode in [ "sk_numode", "sk_anumode" ]:
        for selection in [ "numu", "numub", "nue", "nueb" ]:
            name = mode + "_" + selection + "_bins"
            banff_file.cd()
            obj = banff_file.Get(name)
            new_banff_file.cd()
            obj.Write(name)

    new_banff_file.Close()
    new_banff_file.Close()

    print "Done!"
#_________________________________________________________________________________
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert the new FSI parameters into the unprocessed BANFF matrix")
    parser.add_argument("-b", "--banff", type=str, help="The filename for the original BANFF matrix", required=True)
    parser.add_argument("-v", "--var", type=str, help="The name of the variable to add", required=True)
    parser.add_argument("--output", type=str, help="Name of the output file", required=True)
    args = parser.parse_args()
    if args.var == "eb" or args.var == "intermediate_run10":
        add_param_to_banff(args.banff, args.output, args.var)
    else:
        print "Currently only adding Eb, or intermediate run 10 params to the BANFF matrix are supported. Exiting."
        sys.exit()
