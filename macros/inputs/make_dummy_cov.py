#!/usr/bin/env python

# Make a dummy covariance matrix from arbitrary hardcoded info
#
# Example usage:
# python $VALOR_EXPT/src/macros/inputs/make_dummy_cov.py --output dummy_covariance_matrix.root
#
# Francis Bench

import ROOT
import argparse
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

#_________________________________________________________________________________
def make_cov(outputname):
    # All of these values have been chosen completely arbitrarily
    params_to_add = ["flux_numu00", "flux_numu01", "MAQE", "e_scale"]
    nominal_vals = {
        "flux_numu00" : 1.20,
        "flux_numu01" : 1.10,
        "MAQE"        : 0.90,
        "e_scale"     : 1.00
    }
    untunued_vals = {
        "flux_numu00" : 1.00,
        "flux_numu01" : 1.00,
        "MAQE"        : 1.00,
        "e_scale"     : 1.00
    }
    error_vals = {
        "flux_numu00" : 0.10,
        "flux_numu01" : 0.10,
        "MAQE"        : 0.10,
        "e_scale"     : 0.01
    }
    correlation_vals = {
        "flux_numu00" : { "flux_numu01" : 0.3},
        "flux_numu01" : { "flux_numu00" : 0.3},
    }

    N = len(params_to_add)
    new_matrix  = ROOT.TMatrixD(N, N)
    new_nominal = ROOT.TVectorD(N)
    new_untuned = ROOT.TVectorD(N)
    new_names   = ROOT.TObjArray(N)

    # Add the new param names
    for name in params_to_add:
        new_names.Add(ROOT.TObjString(name))

    # Add the new param nominal values
    for i, name in enumerate(params_to_add):
        new_nominal[i] = nominal_vals[name]

    # Add the new param untuned values
    for i, name in enumerate(params_to_add):
        new_untuned[i] = untunued_vals[name]

    # Get a single list of 1 sigma errors
    new_errors = []
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

    print "Saving covariance matrix to file", outputname
    new_file = ROOT.TFile(outputname, "RECREATE")
    new_names  .Write("param_list", 1)
    new_nominal.Write("params")
    new_untuned.Write("untuned_central_values")
    new_matrix .Write("cov")
    new_file.Close()
    print "Done!"

#_________________________________________________________________________________
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make a covariance matrix from hardcoded info")
    parser.add_argument("--output", type=str, help="Name of the output file", required=True)
    args = parser.parse_args()
    make_cov(args.output)
